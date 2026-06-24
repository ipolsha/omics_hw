import pandas as pd
import numpy as np
import pyBigWig
import bioframe
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

print("Загрузка данных...")
gtf_file = 'chm13v2.0_main_protein_coding_chrNames.gtf.gz'
peak_file = 'MoPh7_H3K27Ac_peaks.narrowPeak'
bw_file = 'MoPh7.rnaseq.STAR.bw'
cov_file = 'data/bismark/MoPh7_1_bismark_bt2_pe.bismark.cov.gz'

# читаем GTF (это GFF3...) 
genes = pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                    usecols=[0, 2, 3, 4, 6, 8], 
                    names=['chrom', 'feature', 'start', 'end', 'strand', 'attr'])

# извлекаем Parent (для CDS/exon) ИЛИ ID (для mRNA/transcript)
# GFF3: строки CDS имеют Parent=rna-NM_..., строки mRNA имеют ID=rna-NM_...
genes['rna_id'] = genes['attr'].str.extract(r'Parent=(rna-[^;,]+)')
# если Parent нет, берём ID (для строк типа mRNA)
genes['rna_id'].fillna(genes['attr'].str.extract(r'ID=(rna-[^;,]+)')[0], inplace=True)

# только строки с валидным rna_id
genes = genes.dropna(subset=['rna_id'])

print(f"  Найдено записей с rna_id: {len(genes)}")

# группируем по rna_id для получения границ транскрипта
transcripts = genes.groupby(['rna_id', 'chrom', 'strand']).agg(
    start=('start', 'min'),
    end=('end', 'max')
).reset_index()

print(f"  Собрано уникальных транскриптов: {len(transcripts)}")

# TSS: для '+' — это start, для '−' — это end
transcripts['tss'] = np.where(transcripts['strand'] == '+', 
                               transcripts['start'], 
                               transcripts['end'])

# промоторы: 1000 bp перед TSS
# для '+': [TSS-1000, TSS]
# для '−': [TSS, TSS+1000]
promoters = pd.DataFrame({
    'chrom': transcripts['chrom'],
    'start': np.where(transcripts['strand'] == '+', 
                      transcripts['tss'] - 1000, 
                      transcripts['tss']),
    'end': np.where(transcripts['strand'] == '+', 
                    transcripts['tss'], 
                    transcripts['tss'] + 1000)
})

# убираем промоторы с отрицательными координатами
promoters = promoters[promoters['start'] >= 0].drop_duplicates().reset_index(drop=True)
print(f"  Сформировано промоторов: {len(promoters)}")

# унификация хромосом (убираем 'chr' если есть)
def unify_chroms(df):
    df['chrom'] = df['chrom'].astype(str).str.replace('^chr', '', regex=True)
    return df

promoters = unify_chroms(promoters)

# пики H3K27Ac
peaks = pd.read_csv(peak_file, sep='\t', header=None, 
                    usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
peaks = unify_chroms(peaks)
print(f"  Найдено пиков H3K27Ac: {len(peaks)}")

# метилирование
meth = pd.read_csv(cov_file, sep='\t', header=None,
                   names=['chrom', 'start', 'end', 'meth_pct', 'count_meth', 'count_unmeth'])
meth = unify_chroms(meth)
meth['coverage'] = meth['count_meth'] + meth['count_unmeth']
meth = meth[meth['coverage'] >= 5]
print(f"  Сайтов метилирования (покрытие >=5): {len(meth)}")

print("\nПересечение данных...")

# H3K27Ac overlap
prom_peaks = bioframe.overlap(promoters, peaks, how='left')
# ищу колонки динамически
pcols = [c for c in prom_peaks.columns if c.startswith('chrom') or c.startswith('start') or c.startswith('end')]
p_chrom = [c for c in prom_peaks.columns if 'chrom' in c][0]
p_start = [c for c in prom_peaks.columns if 'start' in c][0]
p_end = [c for c in prom_peaks.columns if 'end' in c][0]
peak_start = [c for c in prom_peaks.columns if 'start' in c][-1]

prom_peaks['has_h3k27ac'] = prom_peaks[peak_start].notna()
res_df = prom_peaks.groupby([p_chrom, p_start, p_end])['has_h3k27ac'].any().reset_index()
res_df.columns = ['chrom', 'start', 'end', 'has_h3k27ac']

# метилирование overlap
meth_prom = bioframe.overlap(promoters, meth, how='inner')
print(f"  Пересечений промоторов с метилированием: {len(meth_prom)}")

if len(meth_prom) == 0:
    print("!!!!Нет пересечений метилирования с промоторами!!!!")
else:
    m_chrom = [c for c in meth_prom.columns if 'chrom' in c][0]
    m_start = [c for c in meth_prom.columns if 'start' in c][0]
    m_end = [c for c in meth_prom.columns if 'end' in c][0]
    m_pct = [c for c in meth_prom.columns if 'meth_pct' in c][0]
    
    prom_meth = meth_prom.groupby([m_chrom, m_start, m_end])[m_pct].mean().reset_index()
    prom_meth.columns = ['chrom', 'start', 'end', 'meth_pct']
    
    final_df = res_df.merge(prom_meth, on=['chrom', 'start', 'end'], how='left')

print("\nСчитывание RNA-seq из BigWig...")
bw = pyBigWig.open(bw_file)
bw_chroms = set(bw.chroms().keys())

def get_rnaseq(row):
    c = row['chrom']
    # в bigwig обычно 'chr1', 'chr2' и т.д.
    chrom_to_try = c if c in bw_chroms else f"chr{c}"
    if chrom_to_try not in bw_chroms:
        return 0.0
    try:
        val = bw.stats(chrom_to_try, int(row['start']), int(row['end']), 
                       type='mean', exact=True)
        return val[0] if val and val[0] is not None else 0.0
    except:
        return 0.0

final_df['rnaseq_signal'] = final_df.apply(get_rnaseq, axis=1)
bw.close()

# finally
final_df = final_df.dropna(subset=['meth_pct']).copy()
final_df['log_rnaseq'] = np.log1p(final_df['rnaseq_signal'])

h3k_count = final_df['has_h3k27ac'].sum()
total = len(final_df)
print(f"\nИТОГО: {total} промоторов проанализировано")
print(f"   С H3K27Ac: {h3k_count} ({100*h3k_count/total:.1f}%)")

# graphs
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sns.boxplot(data=final_df, x='has_h3k27ac', y='meth_pct', ax=axes[0], palette=['#e74c3c', '#2ecc71'])
axes[0].set_title('Метилирование ДНК в промоторах')
axes[0].set_xticklabels(['Без H3K27Ac', 'С H3K27Ac'])
axes[0].set_ylabel('% метилирования')

sns.boxplot(data=final_df, x='has_h3k27ac', y='log_rnaseq', ax=axes[1], palette=['#e74c3c', '#2ecc71'])
axes[1].set_title('Экспрессия генов (log RNA-seq)')
axes[1].set_xticklabels(['Без H3K27Ac', 'С H3K27Ac'])
axes[1].set_ylabel('log₁(RNA-seq сигнал)')

plt.tight_layout()
plt.savefig('integration_results.png', dpi=300)
print("Графики в integration_results.png")

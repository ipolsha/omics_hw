import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import bioframe
import warnings
warnings.filterwarnings('ignore')

print("Чтение данных метилирования Bismark...")
cov_file = "data/bismark/MoPh7_1_bismark_bt2_pe.bismark.cov.gz"
meth = pd.read_csv(cov_file, sep='\t', header=None,
                   names=['chrom', 'start', 'end', 'meth_pct', 'count_meth', 'count_unmeth'])

# фильтр шум, только сайты с покрытием >= 5
meth['coverage'] = meth['count_meth'] + meth['count_unmeth']
meth = meth[meth['coverage'] >= 5].copy()

print("Скачивание аннотации генов GTF (T2T/hs1)...")
url = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz"
# чтение GTF, пропуская заголовочные строки с #
genes = pd.read_csv(url, sep='\t', comment='#', header=None,
                    usecols=[0, 2, 3, 4, 6],
                    names=['chrom', 'feature', 'start', 'end', 'strand'])

# только транскрипты на основных хромосомах
genes = genes[genes['feature'] == 'transcript'].copy()
valid_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
genes = genes[genes['chrom'].isin(valid_chroms)]

# определение промоторов из day5
genes['tss'] = genes.apply(lambda row: row['start'] if row['strand'] == '+' else row['end'], axis=1)
promoters = pd.DataFrame({
    'chrom': genes['chrom'],
    'start': genes.apply(lambda row: row['tss'] - 1000 if row['strand'] == '+' else row['tss'], axis=1),
    'end': genes.apply(lambda row: row['tss'] if row['strand'] == '+' else row['tss'] + 1000, axis=1)
})

print("Пересечение координат, ищем метилирование в промоторах...")
# метилирование внутри промоторов
meth_prom = bioframe.overlap(meth, promoters, how='inner')
meth_col = 'meth_pct_1' if 'meth_pct_1' in meth_prom.columns else 'meth_pct'

#  средние значения
avg_bg = meth['meth_pct'].mean()
avg_prom = meth_prom[meth_col].mean()
print(f"\nСреднее метилирование промоторов: {avg_prom:.2f}%")
print(f"Среднее метилирование по всему геному: {avg_bg:.2f}%")

print("Отрисовка графика...")
# случайная выборка для скорости
bg_sample = meth.sample(n=min(50000, len(meth)))['meth_pct'].values
prom_sample = meth_prom[meth_col].dropna().sample(n=min(50000, len(meth_prom))).values

plot_df = pd.DataFrame({
    'Метилирование (%)': list(prom_sample) + list(bg_sample),
    'Регион': ['Промоторы (Активные зоны)'] * len(prom_sample) + ['Весь геном (Фон)'] * len(bg_sample)
})

plt.figure(figsize=(9, 6))

# Добавлен параметр cut=0, чтобы скрипки не вылезали за реальные пределы данных (0-100%)
sns.violinplot(x='Регион', y='Метилирование (%)', data=plot_df, palette=['#2ecc71', '#95a5a6'], inner='quartile', cut=0)

plt.title('Распределение уровня метилирования\n(Промоторы vs Весь геном)', fontsize=14, fontweight='bold')
plt.ylabel('% метилированных аллелей', fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Принудительно расширяем ось Y вверх, чтобы освободить место
plt.ylim(-5, 120)

# Настройки для полупрозрачного белого фона под текстом
bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85)

# Добавляем текст на высоту 108 с фоновой подложкой
plt.text(0, 108, f"Среднее: {avg_prom:.1f}%", ha='center', fontsize=12, fontweight='bold', color='#27ae60', bbox=bbox_props)
plt.text(1, 108, f"Среднее: {avg_bg:.1f}%", ha='center', fontsize=12, fontweight='bold', color='#7f8c8d', bbox=bbox_props)

plt.tight_layout()
plt.savefig('methylation_presentation_plot.png', dpi=300)
print("\nГрафик сохранен в methylation_presentation_plot.png")

#!/bin/bash

# остановка при ошибке
set -e

mkdir -p data/raw data/trimmed data/juicer
mkdir -p results/fastqc_raw results/cutadapt results/hic

samples=("MoPh11 S86" "MoPh14 S87" "MoPh15 S88")

# цикл для каждого образца
for info in "${samples[@]}"; do
    # делит строку нна SAMPLE и S_NUM
    read -r SAMPLE S_NUM <<< "$info"
    
    echo "Начало обработки $SAMPLE"

    echo "Загрузка сырых ридов..."
    URL_R1="https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20${SAMPLE}_${S_NUM}_L001_R1_001.fastq.gz"
    URL_R2="https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20${SAMPLE}_${S_NUM}_L001_R2_001.fastq.gz"

    wget --no-check-certificate -O "data/raw/${SAMPLE}_R1.fastq.gz" "$URL_R1"
    wget --no-check-certificate -O "data/raw/${SAMPLE}_R2.fastq.gz" "$URL_R2"

    echo "FastQC..."
    fastqc "data/raw/${SAMPLE}_R1.fastq.gz" "data/raw/${SAMPLE}_R2.fastq.gz" -o results/fastqc_raw

    echo "Обрезка адаптеров..."
    cutadapt \
        -q 20 \
        -m 70 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -o "data/trimmed/${SAMPLE}_R1.trimmed.fastq.gz" \
        -p "data/trimmed/${SAMPLE}_R2.trimmed.fastq.gz" \
        "data/raw/${SAMPLE}_R1.fastq.gz" \
        "data/raw/${SAMPLE}_R2.fastq.gz" \
        > "results/cutadapt/${SAMPLE}.cutadapt.log" 2>&1


    echo "Подготовка директорий для Juicer..."
    mkdir -p "data/juicer/${SAMPLE}/fastq"
    ln -sf "$(pwd)/data/trimmed/${SAMPLE}_R1.trimmed.fastq.gz" "data/juicer/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz"
    ln -sf "$(pwd)/data/trimmed/${SAMPLE}_R2.trimmed.fastq.gz" "data/juicer/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz"

    echo "Juicer..."
    bash tools/juicer/scripts/juicer.sh \
        -D "$(pwd)/tools/juicer" \
        -d "$(pwd)/data/juicer/${SAMPLE}" \
        -g T2T_human \
        -z "$(pwd)/data/reference/T2T_human.fna" \
        -p "$(pwd)/data/reference/chrom.sizes" \
        -y "$(pwd)/data/reference/restriction_sites_DpnII.txt" \
        -s DpnII \
        -t 4

    echo "Копирование .hic карты в results/hic/..."
    cp "data/juicer/${SAMPLE}/aligned/inter_30.hic" "results/hic/${SAMPLE}.inter_30.hic"

    echo "Обработка $SAMPLE завершена."
done

echo "Пайплайн завершил работу. Карты в results/hic/"

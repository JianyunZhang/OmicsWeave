#!/bin/bash

# Bulk RNA-Seq Analysis Pipeline
# Generated automatically

# --- CONFIGURATION ---
# Set the path to where your raw data is located (e.g., mounted network drive or local copy)
RAW_DATA_DIR="./01.raw_data"
OUTPUT_DIR="./analysis_results"
GENOME_INDEX="/path/to/genome/star_index_mm39"  # UPDATE THIS PATH
GTF_FILE="/path/to/genome/gencode.vM27.annotation.gtf"  # UPDATE THIS PATH
THREADS=8

# Create output directories
mkdir -p $OUTPUT_DIR/01.fastqc
mkdir -p $OUTPUT_DIR/02.trimmed
mkdir -p $OUTPUT_DIR/03.alignment
mkdir -p $OUTPUT_DIR/04.counts

# Check for required tools
for tool in fastqc fastp STAR featureCounts; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is not installed or not in PATH."
        exit 1
    fi
done

echo "Starting Analysis..."

#Processing Sample: 52_brain
echo "Processing 52_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/52_brain/52_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/52_brain/52_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/52_brain/52_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/52_brain/52_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/52_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/52_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/52_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/52_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/52_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/52_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/52_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 52_gast_muscle
echo "Processing 52_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/52_gast_muscle/52_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/52_gast_muscle/52_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/52_gast_muscle/52_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/52_gast_muscle/52_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/52_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/52_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/52_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/52_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/52_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/52_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/52_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 52_heart
echo "Processing 52_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/52_heart/52_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/52_heart/52_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/52_heart/52_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/52_heart/52_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/52_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/52_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/52_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/52_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/52_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/52_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/52_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 52_kidney
echo "Processing 52_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/52_kidney/52_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/52_kidney/52_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/52_kidney/52_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/52_kidney/52_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/52_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/52_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/52_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/52_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/52_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/52_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/52_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 52_liver
echo "Processing 52_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/52_liver/52_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/52_liver/52_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/52_liver/52_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/52_liver/52_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/52_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/52_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/52_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/52_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/52_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/52_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/52_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 53_brain
echo "Processing 53_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/53_brain/53_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/53_brain/53_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/53_brain/53_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/53_brain/53_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/53_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/53_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/53_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/53_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/53_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/53_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/53_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 53_gast_muscle
echo "Processing 53_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/53_gast_muscle/53_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/53_gast_muscle/53_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/53_gast_muscle/53_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/53_gast_muscle/53_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/53_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/53_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/53_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/53_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/53_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/53_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/53_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 53_heart
echo "Processing 53_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/53_heart/53_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/53_heart/53_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/53_heart/53_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/53_heart/53_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/53_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/53_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/53_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/53_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/53_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/53_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/53_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 53_kidney
echo "Processing 53_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/53_kidney/53_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/53_kidney/53_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/53_kidney/53_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/53_kidney/53_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/53_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/53_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/53_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/53_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/53_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/53_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/53_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 53_liver
echo "Processing 53_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/53_liver/53_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/53_liver/53_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/53_liver/53_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/53_liver/53_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/53_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/53_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/53_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/53_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/53_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/53_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/53_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 54_brain
echo "Processing 54_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/54_brain/54_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/54_brain/54_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/54_brain/54_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/54_brain/54_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/54_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/54_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/54_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/54_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/54_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/54_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/54_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 54_gast_muscle
echo "Processing 54_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/54_gast_muscle/54_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/54_gast_muscle/54_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/54_gast_muscle/54_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/54_gast_muscle/54_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/54_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/54_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/54_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/54_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/54_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/54_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/54_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 54_heart
echo "Processing 54_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/54_heart/54_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/54_heart/54_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/54_heart/54_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/54_heart/54_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/54_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/54_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/54_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/54_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/54_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/54_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/54_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 54_kidney
echo "Processing 54_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/54_kidney/54_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/54_kidney/54_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/54_kidney/54_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/54_kidney/54_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/54_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/54_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/54_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/54_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/54_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/54_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/54_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 54_liver
echo "Processing 54_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/54_liver/54_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/54_liver/54_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/54_liver/54_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/54_liver/54_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/54_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/54_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/54_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/54_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/54_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/54_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/54_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 54_sq_adipose
echo "Processing 54_sq_adipose..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/54_sq_adipose/54_sq_adipose_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/54_sq_adipose/54_sq_adipose_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/54_sq_adipose/54_sq_adipose_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/54_sq_adipose/54_sq_adipose_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/54_sq_adipose_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/54_sq_adipose_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/54_sq_adipose_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/54_sq_adipose_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/54_sq_adipose_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/54_sq_adipose_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/54_sq_adipose_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 63_brain
echo "Processing 63_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/63_brain/63_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/63_brain/63_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/63_brain/63_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/63_brain/63_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/63_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/63_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/63_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/63_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/63_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/63_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/63_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 63_gast_muscle
echo "Processing 63_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/63_gast_muscle/63_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/63_gast_muscle/63_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/63_gast_muscle/63_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/63_gast_muscle/63_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/63_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/63_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/63_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/63_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/63_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/63_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/63_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 63_heart
echo "Processing 63_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/63_heart/63_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/63_heart/63_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/63_heart/63_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/63_heart/63_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/63_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/63_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/63_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/63_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/63_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/63_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/63_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 63_kidney
echo "Processing 63_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/63_kidney/63_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/63_kidney/63_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/63_kidney/63_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/63_kidney/63_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/63_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/63_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/63_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/63_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/63_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/63_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/63_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 63_liver
echo "Processing 63_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/63_liver/63_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/63_liver/63_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/63_liver/63_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/63_liver/63_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/63_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/63_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/63_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/63_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/63_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/63_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/63_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 64_brain
echo "Processing 64_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/64_brain/64_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/64_brain/64_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/64_brain/64_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/64_brain/64_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/64_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/64_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/64_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/64_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/64_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/64_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/64_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 64_gast_muscle
echo "Processing 64_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/64_gast_muscle/64_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/64_gast_muscle/64_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/64_gast_muscle/64_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/64_gast_muscle/64_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/64_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/64_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/64_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/64_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/64_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/64_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/64_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 64_heart
echo "Processing 64_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/64_heart/64_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/64_heart/64_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/64_heart/64_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/64_heart/64_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/64_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/64_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/64_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/64_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/64_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/64_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/64_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 64_kidney
echo "Processing 64_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/64_kidney/64_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/64_kidney/64_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/64_kidney/64_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/64_kidney/64_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/64_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/64_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/64_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/64_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/64_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/64_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/64_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 64_liver
echo "Processing 64_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/64_liver/64_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/64_liver/64_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/64_liver/64_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/64_liver/64_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/64_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/64_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/64_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/64_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/64_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/64_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/64_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 64_sq_adipose
echo "Processing 64_sq_adipose..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/64_sq_adipose/64_sq_adipose_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/64_sq_adipose/64_sq_adipose_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/64_sq_adipose/64_sq_adipose_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/64_sq_adipose/64_sq_adipose_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/64_sq_adipose_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/64_sq_adipose_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/64_sq_adipose_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/64_sq_adipose_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/64_sq_adipose_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/64_sq_adipose_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/64_sq_adipose_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 65_brain
echo "Processing 65_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/65_brain/65_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/65_brain/65_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/65_brain/65_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/65_brain/65_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/65_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/65_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/65_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/65_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/65_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/65_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/65_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 65_gast_muscle
echo "Processing 65_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/65_gast_muscle/65_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/65_gast_muscle/65_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/65_gast_muscle/65_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/65_gast_muscle/65_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/65_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/65_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/65_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/65_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/65_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/65_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/65_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 65_heart
echo "Processing 65_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/65_heart/65_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/65_heart/65_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/65_heart/65_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/65_heart/65_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/65_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/65_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/65_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/65_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/65_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/65_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/65_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 65_kidney
echo "Processing 65_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/65_kidney/65_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/65_kidney/65_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/65_kidney/65_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/65_kidney/65_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/65_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/65_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/65_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/65_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/65_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/65_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/65_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 65_liver
echo "Processing 65_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/65_liver/65_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/65_liver/65_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/65_liver/65_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/65_liver/65_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/65_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/65_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/65_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/65_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/65_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/65_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/65_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 65_sq_adipose
echo "Processing 65_sq_adipose..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/65_sq_adipose/65_sq_adipose_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/65_sq_adipose/65_sq_adipose_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/65_sq_adipose/65_sq_adipose_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/65_sq_adipose/65_sq_adipose_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/65_sq_adipose_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/65_sq_adipose_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/65_sq_adipose_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/65_sq_adipose_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/65_sq_adipose_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/65_sq_adipose_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/65_sq_adipose_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 73_brain
echo "Processing 73_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/73_brain/73_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/73_brain/73_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/73_brain/73_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/73_brain/73_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/73_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/73_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/73_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/73_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/73_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/73_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/73_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 73_gast_muscle
echo "Processing 73_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/73_gast_muscle/73_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/73_gast_muscle/73_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/73_gast_muscle/73_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/73_gast_muscle/73_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/73_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/73_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/73_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/73_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/73_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/73_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/73_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 73_heart
echo "Processing 73_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/73_heart/73_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/73_heart/73_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/73_heart/73_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/73_heart/73_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/73_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/73_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/73_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/73_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/73_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/73_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/73_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 73_kidney
echo "Processing 73_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/73_kidney/73_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/73_kidney/73_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/73_kidney/73_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/73_kidney/73_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/73_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/73_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/73_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/73_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/73_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/73_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/73_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 73_liver
echo "Processing 73_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/73_liver/73_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/73_liver/73_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/73_liver/73_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/73_liver/73_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/73_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/73_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/73_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/73_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/73_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/73_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/73_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 73_sq_adipose
echo "Processing 73_sq_adipose..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/73_sq_adipose/73_sq_adipose_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/73_sq_adipose/73_sq_adipose_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/73_sq_adipose/73_sq_adipose_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/73_sq_adipose/73_sq_adipose_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/73_sq_adipose_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/73_sq_adipose_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/73_sq_adipose_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/73_sq_adipose_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/73_sq_adipose_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/73_sq_adipose_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/73_sq_adipose_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 74_brain
echo "Processing 74_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/74_brain/74_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/74_brain/74_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/74_brain/74_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/74_brain/74_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/74_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/74_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/74_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/74_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/74_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/74_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/74_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 74_gast_muscle
echo "Processing 74_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/74_gast_muscle/74_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/74_gast_muscle/74_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/74_gast_muscle/74_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/74_gast_muscle/74_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/74_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/74_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/74_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/74_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/74_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/74_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/74_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 74_heart
echo "Processing 74_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/74_heart/74_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/74_heart/74_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/74_heart/74_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/74_heart/74_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/74_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/74_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/74_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/74_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/74_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/74_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/74_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 74_kidney
echo "Processing 74_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/74_kidney/74_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/74_kidney/74_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/74_kidney/74_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/74_kidney/74_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/74_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/74_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/74_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/74_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/74_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/74_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/74_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 74_liver
echo "Processing 74_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/74_liver/74_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/74_liver/74_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/74_liver/74_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/74_liver/74_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/74_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/74_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/74_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/74_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/74_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/74_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/74_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 74_sq_adipose
echo "Processing 74_sq_adipose..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/74_sq_adipose/74_sq_adipose_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/74_sq_adipose/74_sq_adipose_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/74_sq_adipose/74_sq_adipose_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/74_sq_adipose/74_sq_adipose_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/74_sq_adipose_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/74_sq_adipose_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/74_sq_adipose_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/74_sq_adipose_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/74_sq_adipose_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/74_sq_adipose_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/74_sq_adipose_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 75_brain
echo "Processing 75_brain..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/75_brain/75_brain_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/75_brain/75_brain_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/75_brain/75_brain_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/75_brain/75_brain_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/75_brain_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/75_brain_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/75_brain_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/75_brain_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/75_brain_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/75_brain_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/75_brain_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 75_gast_muscle
echo "Processing 75_gast_muscle..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/75_gast_muscle/75_gast_muscle_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/75_gast_muscle/75_gast_muscle_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/75_gast_muscle/75_gast_muscle_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/75_gast_muscle/75_gast_muscle_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/75_gast_muscle_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/75_gast_muscle_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/75_gast_muscle_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/75_gast_muscle_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/75_gast_muscle_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/75_gast_muscle_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/75_gast_muscle_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 75_heart
echo "Processing 75_heart..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/75_heart/75_heart_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/75_heart/75_heart_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/75_heart/75_heart_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/75_heart/75_heart_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/75_heart_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/75_heart_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/75_heart_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/75_heart_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/75_heart_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/75_heart_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/75_heart_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 75_kidney
echo "Processing 75_kidney..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/75_kidney/75_kidney_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/75_kidney/75_kidney_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/75_kidney/75_kidney_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/75_kidney/75_kidney_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/75_kidney_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/75_kidney_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/75_kidney_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/75_kidney_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/75_kidney_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/75_kidney_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/75_kidney_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 75_liver
echo "Processing 75_liver..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/75_liver/75_liver_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/75_liver/75_liver_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/75_liver/75_liver_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/75_liver/75_liver_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/75_liver_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/75_liver_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/75_liver_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/75_liver_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/75_liver_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/75_liver_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/75_liver_ \
     --outSAMtype BAM SortedByCoordinate

#Processing Sample: 75_sq_adipose
echo "Processing 75_sq_adipose..."
echo "  - Running FastQC..."
fastqc -t $THREADS -o $OUTPUT_DIR/01.fastqc ${RAW_DATA_DIR}/75_sq_adipose/75_sq_adipose_S1_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/75_sq_adipose/75_sq_adipose_S1_L001_R2_001.fastq.gz
echo "  - Running fastp (trimming)..."
fastp -i ${RAW_DATA_DIR}/75_sq_adipose/75_sq_adipose_S1_L001_R1_001.fastq.gz -I ${RAW_DATA_DIR}/75_sq_adipose/75_sq_adipose_S1_L001_R2_001.fastq.gz \
      -o $OUTPUT_DIR/02.trimmed/75_sq_adipose_R1.clean.fq.gz \
      -O $OUTPUT_DIR/02.trimmed/75_sq_adipose_R2.clean.fq.gz \
      -h $OUTPUT_DIR/02.trimmed/75_sq_adipose_fastp.html \
      -j $OUTPUT_DIR/02.trimmed/75_sq_adipose_fastp.json \
      --thread $THREADS
echo "  - Running STAR (alignment)..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $OUTPUT_DIR/02.trimmed/75_sq_adipose_R1.clean.fq.gz $OUTPUT_DIR/02.trimmed/75_sq_adipose_R2.clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/03.alignment/75_sq_adipose_ \
     --outSAMtype BAM SortedByCoordinate

# 4. featureCounts (Quantification)
echo "Running featureCounts on all BAM files..."
featureCounts -T $THREADS -p -t exon -g gene_id \
              -a $GTF_FILE \
              -o $OUTPUT_DIR/04.counts/counts.txt \
              $OUTPUT_DIR/03.alignment/*_Aligned.sortedByCoord.out.bam

echo "Analysis Pipeline Completed!"

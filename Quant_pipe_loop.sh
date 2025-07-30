#!/bin/bash

# Exit on error
set -e

# Check input arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_directory> <genome_directory> <output_directory>"
    exit 1
fi

# Assign input arguments to variables
INPUT_DIR=$1
GENOME_DIR=$2
OUTPUT_DIR=$3

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all FASTQ files in the input directory
for INPUT_FASTQ in "$INPUT_DIR"/*.fq.gz; do
    # Check if the file is a regular file and not a directory
    if [ -f "$INPUT_FASTQ" ]; then
        # Extract basename without extension
        BASENAME=$(basename "$INPUT_FASTQ" | cut -d. -f1)

        # Step 1: Quality Control with FastQC
        echo "Running FastQC for $BASENAME..."
        fastqc "$INPUT_FASTQ" -o "$OUTPUT_DIR"

        # Step 2: Trimming with Cutadapt
        R1_trimmed="$OUTPUT_DIR/${BASENAME}_trimmed.fastq"
        echo "Trimming with Cutadapt for $BASENAME..."
        cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 "$INPUT_FASTQ" | \
        cutadapt -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
        cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o "$R1_trimmed" -

        # Step 3: Alignment and Gene Counting with STAR
        echo "Aligning and counting with STAR for $BASENAME..."
        STAR --genomeDir "$GENOME_DIR" --readFilesIn "$R1_trimmed" --outFileNamePrefix "$OUTPUT_DIR/${BASENAME}_" --quantMode GeneCounts --alignIntronMax 10000 --outSAMmultNmax 20 --outSAMtype BAM SortedByCoordinate
    fi
done

echo "Pipeline completed successfully!"

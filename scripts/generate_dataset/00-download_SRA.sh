#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --array=4-12
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/gvangeest/Courses/20220228_ISCTR/log/download_output_%j.txt
#SBATCH --error=/data/users/gvangeest/Courses/20220228_ISCTR/log/download_error_%j.txt
#SBATCH --job-name=download
#SBATCH --partition=pall


# activate conda environment
module add UHTS/Analysis/sratoolkit/2.10.7

# get sample information from the table
INFO=$(head -n "$SLURM_ARRAY_TASK_ID" sample_info.csv | tail -n 1)

# split relevant parts
ID=$(echo $INFO | cut -d "," -f 1)
SAMPLE=$(echo $INFO | cut -d "," -f 2)
TYPE=$(echo $INFO | cut -d "," -f 3)
LINK=$(echo $INFO | cut -d "," -f 5)

# move to data folder
mkdir -p data/reads/
cd data/reads/

# create directory structure
OUTDIR="${SAMPLE}/"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# download file
# wget -O "${ID}.sra" "$LINK"

# convert to fastq
fastq-dump -O . --gzip --split-files "${ID}.sra"
# rename files for cellranger
# 1 = I1
# 2 = R1
# 3 = R2
mv "${ID}_1.fastq.gz" "${SAMPLE}_S1_L001_I1_001.fastq.gz"
mv "${ID}_2.fastq.gz" "${SAMPLE}_S1_L001_R1_001.fastq.gz"
mv "${ID}_3.fastq.gz" "${SAMPLE}_S1_L001_R2_001.fastq.gz"

# remove original file
rm "${ID}.sra"

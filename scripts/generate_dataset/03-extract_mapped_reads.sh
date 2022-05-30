#!/usr/bin/env bash

#SBATCH --time=05:00:00
#SBATCH --array=2-3
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=12
#SBATCH --output=/data/users/gvangeest/Courses/20220228_ISCTR/log/extract_output_%j.txt
#SBATCH --error=/data/users/gvangeest/Courses/20220228_ISCTR/log/extract_error_%j.txt
#SBATCH --job-name=extract
#SBATCH --partition=pall

module add UHTS/Analysis/seqtk/1.2
module add UHTS/Analysis/samtools/1.10

# get sample information from the table
INFO=$(head -n "$SLURM_ARRAY_TASK_ID" sample_info.csv | tail -n 1)

# split relevant parts
ID=$(echo $INFO | cut -d "," -f 1)
SAMPLE=$(echo $INFO | cut -d "," -f 2)
TYPE=$(echo $INFO | cut -d "," -f 3)

# extract aligned read names (note that only R2 is aligned to genome)
samtools view -F 4 -h -@ $(($SLURM_CPUS_PER_TASK - 1)) \
results/cellranger/${SAMPLE}/outs/possorted_genome_bam.bam 21 22 \
| samtools sort -@ $(($SLURM_CPUS_PER_TASK - 1)) -u \
| samtools fastq -@ $(($SLURM_CPUS_PER_TASK - 1)) \
| sed -n '1~4p' | sed 's/^@//' > results/${SAMPLE}_aligned_read_names.txt

# extract reads from original files
mkdir -p data/reads/subsampled

seqtk subseq \
data/reads/${SAMPLE}/${SAMPLE}_S1_L001_R1_001.fastq.gz \
results/${SAMPLE}_aligned_read_names.txt \
| seqtk sample -s1 - 1000000 \
| gzip \
> data/reads/subsampled/${SAMPLE}_S1_L001_R1_001.fastq.gz

seqtk subseq \
data/reads/${SAMPLE}/${SAMPLE}_S1_L001_R2_001.fastq.gz \
results/${SAMPLE}_aligned_read_names.txt \
| seqtk sample -s1 - 1000000 \
| gzip \
> data/reads/subsampled/${SAMPLE}_S1_L001_R2_001.fastq.gz

seqtk subseq \
data/reads/${SAMPLE}/${SAMPLE}_S1_L001_I1_001.fastq.gz \
results/${SAMPLE}_aligned_read_names.txt \
| seqtk sample -s1 - 1000000 \
| gzip \
> data/reads/subsampled/${SAMPLE}_S1_L001_I1_001.fastq.gz

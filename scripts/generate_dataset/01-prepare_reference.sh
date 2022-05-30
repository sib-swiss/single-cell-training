#!/usr/bin/env bash

#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --output=/data/users/gvangeest/Courses/20220228_ISCTR/log/ref_output_%j.txt
#SBATCH --error=/data/users/gvangeest/Courses/20220228_ISCTR/log/ref_error_%j.txt
#SBATCH --job-name=reference
#SBATCH --partition=pall

mkdir -p data/reference/
cd data/reference/


#### Download reference genome ####

# download two chromosomes only (to reduce size of the data)
wget --no-check-certificate -O chr21.fa.gz https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
wget --no-check-certificate -O chr22.fa.gz https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

# concatenate
cat chr21.fa.gz chr22.fa.gz > genome.fa.gz 
rm chr21.fa.gz chr22.fa.gz 

# unzip
gunzip genome.fa.gz


#### Download and prepare gene annotation ####

# download from ENSEMBL
wget --no-check-certificate -O annotation.gtf.gz https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz

# decompress (for cellranger)
gunzip annotation.gtf.gz

grep -E "^21|^22|^#" annotation.gtf > annotation.21.22.gtf

rm annotation.gtf

#### cellranger indexing ####
module add UHTS/SingleCell/cellranger/6.0.1

cellranger mkref \
  --nthreads="${SLURM_CPUS_PER_TASK}" \
  --genome=cellranger_index \
  --fasta=genome.fa \
  --genes=annotation.21.22.gtf

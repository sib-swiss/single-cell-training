#!/usr/bin/env bash

#SBATCH --time=3-00:00:00
#SBATCH --array=2-3
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=24
#SBATCH --output=/data/users/gvangeest/Courses/20220228_ISCTR/log/cellranger_output_%j.txt
#SBATCH --error=/data/users/gvangeest/Courses/20220228_ISCTR/log/cellranger_error_%j.txt
#SBATCH --job-name=cellranger
#SBATCH --partition=pall

echo "Start: $(date)"

PROJDIR=/data/users/gvangeest/Courses/20220228_ISCTR

cd $PROJDIR/make_data

# get sample information from the table
INFO=$(head -n "$SLURM_ARRAY_TASK_ID" sample_info.csv | tail -n 1)

# split relevant parts
SAMPLE=$(echo $INFO | cut -d "," -f 2)
TYPE=$(echo $INFO | cut -d "," -f 3)
ID=$(echo $INFO | cut -d "," -f 1)

echo "Processing $SAMPLE"

# rename SRR ids to sample
cd $PROJDIR/make_data/data/reads/"${SAMPLE}"
rename "${ID}" "${SAMPLE}" *.fastq.gz

# reference genome
cd $PROJDIR/make_data/
REF="data/reference/cellranger_index"

# run cellranger count pipeline
module add UHTS/SingleCell/cellranger/6.0.1

cellranger count \
  --id="${SAMPLE}" \
  --sample="${SAMPLE}" \
  --transcriptome="${REF}" \
  --fastqs="data/reads/${SAMPLE}/" \
  --localcores="${SLURM_CPUS_PER_TASK}" \
  --mempercore="${SLURM_MEM_PER_CPU}"

# move to output directory
mkdir -p results/cellranger/
mv "${SAMPLE}" results/cellranger/

echo "End: $(date)"

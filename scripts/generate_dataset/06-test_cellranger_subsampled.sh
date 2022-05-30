#!/usr/bin/env bash

#SBATCH --time=3-00:00:00
#SBATCH --array=0-1
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/users/gvangeest/Courses/20220228_ISCTR/log/test_cellranger_output_%j.txt
#SBATCH --error=/data/users/gvangeest/Courses/20220228_ISCTR/log/test_cellranger_error_%j.txt
#SBATCH --job-name=cellranger
#SBATCH --partition=pall

echo "Start: $(date)"

PROJDIR=/data/users/gvangeest/Courses/20220228_ISCTR

SAMPLES=(ETV6-RUNX1_1 ETV6-RUNX1_2)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "Processing $SAMPLE"

# reference genome
REF="$PROJDIR/make_data/data/reference/cellranger_index"

OUTDIR="$PROJDIR/make_data/results/cellranger_subsampled"
mkdir -p $OUTDIR
cd $OUTDIR

# run cellranger count pipeline
module add UHTS/SingleCell/cellranger/6.0.1

cellranger count \
--id="${SAMPLE}" \
--sample="${SAMPLE}" \
--transcriptome="${REF}" \
--fastqs="$PROJDIR/make_data/data/reads/subsampled/" \
--localcores="${SLURM_CPUS_PER_TASK}" \
--mempercore="${SLURM_MEM_PER_CPU}"

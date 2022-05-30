#!/usr/bin/env bash


#SBATCH --time=01:00:00
#SBATCH --array=2-12
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/users/gvangeest/Courses/20220228_ISCTR/log/move_counts_output_%j.txt
#SBATCH --error=/data/users/gvangeest/Courses/20220228_ISCTR/log/move_counts_error_%j.txt
#SBATCH --job-name=move_counts
#SBATCH --partition=pall

# get sample information from the table
INFO=$(head -n "$SLURM_ARRAY_TASK_ID" sample_info.csv | tail -n 1)

# split relevant parts
ID=$(echo $INFO | cut -d "," -f 1)
SAMPLE=$(echo $INFO | cut -d "," -f 2)
TYPE=$(echo $INFO | cut -d "," -f 3)

# count matrices downloaded from GSE132509 (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132509&format=file)
cd count_matrices
# tar xvf GSE132509_RAW.tar 

mkdir -p "${SAMPLE}"/outs/filtered_feature_bc_matrix

for FILE in GSM*_"${SAMPLE}".*.gz
do 
    EXT=`echo $FILE | cut -f 2-4 -d "."`
    mv $FILE "${SAMPLE}"/outs/filtered_feature_bc_matrix/"${EXT}"
    
    mv "${SAMPLE}"/outs/filtered_feature_bc_matrix/genes.tsv.gz \
    "${SAMPLE}"/outs/filtered_feature_bc_matrix/features.tsv.gz
done


#!/usr/bin/env bash
# run this commands in the admin container
# e.g. docker run ubuntu -it --rm -v data:/data ubuntu

sudo su
cd /home/rstudio/data
wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1677264067&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NzcyNjQwNjd9fX1dfQ__&Signature=L45--zcMjqMTYgW4KS7ClV5CGlZjBc24lZ9ft-KC0xlH3NtHCTfLAnhT0gT30ptsJG52kE2cwBQkXO~THLW7RJP8eUyELXY25PnUdWXw5IUnd9EK5THWbooA8Wdb77SVq3Sd3cCffuRCaIU2YsY~ej--Ba8uU8y8aW8w9OnuUEgaPSbuxH9mMArp5Jc7nDMfmIyZn~1mRFhNNTC5ucVvX5JLcYKukIN71pXlVnNoEt9AJ6UcCGUXteCts6NYeAdSVFL1blVb6d~uOpjk7TwzehJFfURO7dDtYHFaOP0F6tSqcVnXoYiMH8DHDIe2FK9nYjplFz8Wds20Ozf3irvQ6Q__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -xzvf cellranger-7.1.0.tar.gz

# now users need to add the cellranger directory to their path with
# export PATH=$PATH:/data/cellranger-7.1.0
# for users data is mounted to the root directory

# also download the reference data
wget https://single-cell-transcriptomics.s3.eu-central-1.amazonaws.com/cellranger_index.tar.gz
tar -xvf cellranger_index.tar.gz
rm cellranger_index.tar.gz

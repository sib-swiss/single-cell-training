FROM rocker/rstudio:4.1.2
#FROM rocker/rstudio:4.0.5

RUN sudo apt-get update

RUN sudo apt-get install -y \
libxt6 \
libgl1-mesa-glx \
libcurl4-openssl-dev \
libglpk-dev \
libxml2-dev

COPY install_packages.R ./

RUN Rscript install_packages.R

# add cellranger to PATH
# check installation Seurat

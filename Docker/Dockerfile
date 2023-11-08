FROM rocker/rstudio:4.3.2

# Script for to set up multiple users
# COPY batch_user_creation.sh /rocker_scripts/experimental/batch_user_creation.sh

# CMD [ "sh", "-c", "/rocker_scripts/experimental/batch_user_creation.sh && /init" ]

RUN sudo apt-get update

RUN sudo apt-get install -y \
libxt6 \
libgl1-mesa-glx \
libcurl4-openssl-dev \
libglpk-dev \
libxml2-dev \
libproj-dev \
libudunits2-dev \
libgdal-dev \
tree \
libcairo2-dev \
libxt-dev

COPY install_packages.R ./

RUN Rscript install_packages.R

# add cellranger to PATH
# check installation Seurat

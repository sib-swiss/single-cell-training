FROM rocker/rstudio:4.0.5

RUN sudo apt-get update

RUN sudo apt-get install -y libxml2-dev

RUN sudo apt-get install -y libglpk-dev

RUN sudo apt-get install -y libcurl4-openssl-dev

RUN sudo apt-get install -y libxtst6

COPY install_packages.R ./

RUN Rscript install_packages.R

docker run \
-p 8787:8787 \
--name rstudio_admin \
-e PASSWORD=test123 \
-v $PWD:/home/rstudio \
geertvangeest/single-cell-rstudio:2022.2

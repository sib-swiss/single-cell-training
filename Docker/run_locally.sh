docker run \
--rm \
-p 8787:8787 \
-e PASSWORD=test \
-v $PWD:/home/rstudio \
geertvangeest/single-cell-rstudio:2021.6

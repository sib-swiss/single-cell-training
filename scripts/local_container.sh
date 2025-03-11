docker run \
-p 8787:8787 \
--rm \
--name rstudio_admin \
-e PASSWORD=test123 \
-v $PWD:/home/rstudio \
geertvangeest/ngs-longreads-vscode:latest

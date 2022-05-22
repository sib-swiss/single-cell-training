## To deploy

https://leangaurav.medium.com/simplest-https-setup-nginx-reverse-proxy-letsencrypt-ssl-certificate-aws-cloud-docker-4b74569b3c61

Make sure the IP is correctly associated to the domain (`sib-training-test.ml`)

First host the app over http

On the server, generate the certificates:

```sh
docker compose -f docker-compose-le.yaml up --build
```

Then run: 


```sh
docker compose up --build -d nginx
```


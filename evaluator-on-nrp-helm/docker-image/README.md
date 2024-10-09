# Creating a Docker image


Here is an example of creating an image and pushing it to GitHub. This image is used for the Evaluator server pod.

```
IMAGE_NAME='tmp-evaluator-ambertools-distributed-server-v11'
docker build --platform linux/amd64 --tag $IMAGE_NAME .
docker tag $IMAGE_NAME "ghcr.io/lilyminium/scrapbook-projects:${IMAGE_NAME}"
docker push "ghcr.io/lilyminium/scrapbook-projects:${IMAGE_NAME}"
```

The Dockerfile and run-server.sh scripts copy directly from the Dask docker image, which as of this time of writing is BSD-3 licensed.
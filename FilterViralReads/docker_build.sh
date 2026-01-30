#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define the image tag here for easy updating
IMAGE_TAG="fleharty/viral-bam-filter:v2"

echo "Building image: $IMAGE_TAG for linux/amd64..."
docker build --platform linux/amd64 -t "$IMAGE_TAG" .

echo "Pushing image to Docker Hub..."
docker push "$IMAGE_TAG"

echo "Success! Image pushed: $IMAGE_TAG"



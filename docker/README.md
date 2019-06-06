# Dockerfile for PetIBM-0.4

The present directory contains the `Dockerfile` used to create the Docker image `petibm:0.4.1-GPU-OpenMPI-ubuntu`.
The Docker image is available on DockerHub; to pull the image:

```shell
docker pull barbagroup/petibm:0.4.1-GPU-OpenMPI-ubuntu
```

In this image, PetIBM (0.4.1) was installed along with its dependencies:

* CUDA Toolkit 8.0
* OpenMPI 1.10.2
* PETSc 3.10.2
* AmgX 2.0 (commit [732338c](https://github.com/NVIDIA/AMGX/tree/732338c32e30ad87f9b71244346346f66fc3f735))
* AmgXWrapper 1.4

Docker containers created out of this image can run PetIBM applications on either CPU or GPU devices.
Note that if you want to run a PetIBM application using GPUs, you have to create the container using the utility `nvidia-docker` (instead of `docker`):

```shell
nvidia-docker run -it barbagroup/petibm:0.4.1-GPU-OpenMPI-ubuntu /bin/bash
```

If you wish to re-build the image locally:

```shell
docker build --tag=mypetibm:latest --file=Dockerfile-0.4.1-GPU-OpenMPI-ubuntu .
```
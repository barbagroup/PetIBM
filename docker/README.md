# Dockerfile for PetIBM-0.4.2

The present directory contains the `Dockerfile` used to create the Docker image `petibm:0.4.2-GPU-OpenMPI-xenial`.
The Docker image, based on Ubuntu 16.04 (Xenial) is available on DockerHub; to pull the image:

```shell
docker pull barbagroup/petibm:0.4.2-GPU-OpenMPI-xenial
```

In this image, PetIBM (0.4.2) was installed along with its dependencies:

* CUDA Toolkit 10.1 (requires CUDA Driver Version >= 418.39)
* OpenMPI 3.0.1
* PETSc 3.11.3
* AmgX 2.0 (commit [a46b311][732338c](https://github.com/NVIDIA/AMGX/tree/a46b3112bc563592b8d794ba95e57350d282d584))
* AmgXWrapper 1.4

Docker containers based on this image can run PetIBM applications on CPU and GPU.
Note that if you want to run a PetIBM application using GPUs, you have to create the container using the utility [`nvidia-docker v2`](https://github.com/NVIDIA/nvidia-docker) (instead of `docker`):

```shell
nvidia-docker run -it barbagroup/petibm:0.4.2-GPU-OpenMPI-xenial /bin/bash
```

If you want to re-build the image locally:

```shell
docker build --tag=mypetibm:latest --file=Dockerfile-0.4.2-GPU-OpenMPI-xenial .
```

The Docker image `barbagroup/petibm:0.4.2-GPU-OpenMPI-xenial` does not contain the source files and builds.
This was intended to make the Docker image as light as possible.
If you wish to have the source files and builds, you can use the bigger Docker image `barbagroup/petibm:0.4.2-GPU-OpenMPI-xenial-devel`.

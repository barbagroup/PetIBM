## Dependencies
----------------

* [OpenMPI 1.8](https://www.open-mpi.org/software/ompi/v1.8/)
* [CUDA 6.5](https://developer.nvidia.com/cuda-toolkit-65)
* [PETSc](https://www.mcs.anl.gov/petsc/) (versions above 3.5 has been tested)
* [AmgX](https://developer.nvidia.com/amgx)
  (you can download from [here](https://developer.nvidia.com/rdp/assets/amgx-trial-download))
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) (optional)

Currently, NVIDIA AmgX only provides compiled binary library. The binary library
seems to be compiled with OpenMPI 1.8 and CUDA 6.5 before being released by NVIDIA, 
so the CUDA version is required to be exactly 6.5. However, for OnpeMPI, we have
tried version 1.10, and it worked. The only one thing of using OpenMPI 1.10 is 
that you have to create a fake OpenMPI 1.8 library which is a soft link to 1.10.
For example, using

```bash
$ ln -s ${PATH_TO_libmpi.so} ${PATH_TO_ANYWHERE_YOU_LIKE}/libmpi.so.1
```

Doxygen is required only when you want to build API manual.

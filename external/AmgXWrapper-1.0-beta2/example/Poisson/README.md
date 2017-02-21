To build:

```
$ cmake . -DCUDA_TOOLKIT_ROOT_DIR=<installation of CUDA 6.5> -DAMGX_DIR=<installation of AmgX> -DPETSC_DIR=<installation of PETSc>

$ make
```

Exectuable file will be located in the folder `bin`. To run, for example, a 100x100x100 3D Poisson benchmark using AmgX solver and two GPUs on one node:

```
$ mpiexec -map-by ppr:2:node -n 2 bin/Poisson -caseName test -mode GPU -Nx 100 -Ny 100 -Nz 100 -cfgFileName configs/AmgX_SolverOptions_classical.info
```

Or if PETSc KSP solver and 16 CPUs are preferred:

```
$ mpiexec -map-by ppr:16:node -n 16 bin/Poisson -caseName test -mode PETSc -Nx 100 -Ny 100 -Nz 100 -cfgFileName configs/PETSc_SolverOptions_hypre.info -optFileName performance.txt
```

For 2D Poisson benchmarks, do not specify `Nz` or set `Nz` as zero. For example:

```
$ mpiexec -map-by ppr:2:node -n 2 bin/Poisson -caseName test -mode GPU -Nx 100 -Ny 100 -cfgFileName configs/AmgX_SolverOptions_classical.info
```

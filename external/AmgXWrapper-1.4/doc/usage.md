## Usage
---------

Using AmgX is simple (and this is what we intended to do!).

**Reminder**: before running your application, make sure your executable binary
can find necessary library (OpenMPI 1.8, CUDA 6.5, AmgX) and the license file 
for AmgX during runtime. You can use environmental variable `LM_LICENSE_FILE`
to indicate the path of the license.

### Step 1

In your PETSc 
application, after `PetscInitialize(...)`, declare an instance of `AmgXSolver`.

```c++
AmgXSolver      solver;
```

### Step 2

Next, at any location you like, initialize this instance.

```c++
solver.initialize(comm, mode, configFile);
```

1. `comm` is your preferred MPI communicator. Typically, this is `MPI_COMM_WORLD`
or `PETSC_COMM_WORLD`.

2. `mode` is a `std::string` indicating AmgX solver mode. `mode` string is 
composed of four letters. The first one letter (lowercase) indicates if AmgX
solver will run on GPU (`d`) or CPU (`h`). The second letter (uppercase) indicates
the precision of matrix entries: float (`F`) or double (`D`). The third one
(uppercase) indicates the precision of the vectors. The last one indicates the
integer type of the indices of matrix or vector entries (currently only 32bit
integer is supported, so only `I` is available). For example, `mode = dDDI` means
the AmgX solver will run on GPUs, and using double precision floating numbers 
for matrix and vector entries.

3. `configFile` is a `std::string` indicating the path of AmgX configuration file. 
Please see AmgX Reference Manual for writing configuration files. For example, 
[here](../example/poisson/configs/AmgX_SolverOptions_Classical.info) 
is a typical solver and preconditioner settings we used in PetIBM simulations.

The return values of all AmgXWrapper functions are `PetscErrorCode`, so you can
call these functions and check the returns in PETSc style:

```c++
ierr = solver.initialize(comm, mode, configFile); CHKERRQ(ierr);
```

You can also combine the first and this step:

```c++
AmgXSolver      solver(comm, mode, configFile);
```

### Step 3

After you finish creating the coefficient matrix using PETSc, upload the 
matrix to GPUs and bind to the `solver` with

```c++
ierr = solver.setA(A); CHKERRQ(ierr);
```

`A` is the coefficient matrix with data type PETSc `Mat`.

No matter how many GPUs and how many CPU cores or nodes you are using, `setA` 
will handle data gathering/scattering for you.

### Step 4

After you have the right-hand-side vector, you can solve the system now:

```c++
ierr = solver.solve(lhs, rhs); CHKERRQ(ierr);
```

`lhs` is a PETSc `Vec` type vector holding the initial guess of solution. Also,
`lhs` will become the resulting solution in the end of solving. `rhs` is a
PETSc `Vec` type vector for the right hand side. Matrix `A`, `lhs`, and `rhs` 
must be created with the same MPI communicator.

For time-marching problems, the same linear system may be solved again and again
with different right-hand-side values. In this case, simply call `solver.solve`
again with updated right-hand-side vector.

### Step 5 (optional)

If you are interested in the number of iterations used in a solve, you can use

```c++
int     iters;
ierr = solver.getIters(iters); CHKERRQ(ierr);
```

The `iters` is the number of iterations in the last solve.

`solver.getResidual(...)` can be used to obtain the residual at any iteration
in the last solve. For example, to get the residual at the second iteration,
use

```c++
double      res;
ierr = solver.getResidual(2, res); CHKERRQ(ierr);
```

Combining `getIters` and `getResidual`, you can then obtain the final residual.

### Step 6

Finalization can be done manually:

```c++
ierr = solver.finalize(); CHKERRQ(ierr);
```

The destructor of the class `AmgXSolver` also does the same thing. This benefits
when using pointer. For example,

```c++
AmgXSolver  *solver = new AmgXSolver(comm, mode, configFile);

solver->setA(A);

solver->solve(lhs, rhs);

delete solver;
```

**Note:**  
If you care about memory leak, you must finalize all `AmgXSolver` instances
before calling `PetscFinalize()`. This is because there are some PETSc data
in `AmgXSolver` instances.

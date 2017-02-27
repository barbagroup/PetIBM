# AmgXWrapper

This wrapper simplifies the usage of AmgX, especially when using AmgX together with other libraries, like PETSc. Currently we only support its use with PETSc, but in the future we may extend it to work with other libraries.

The usage is simple. Just follow the procedure: ***initialize -> set matrix A -> solve -> finalize***. 

For example,

```c++
int main(int argc, char **argv)
{
    // initialize matrix A, RHS, etc using PETSc
    ...

    // create an instance of the solver wrapper
    AmgXSolver    solver;
    // initialize the instance
    solver.initialize(...);    
    // set matrix A. Currently it only accept PETSc AIJ matrix
    solver.setA(A);    
    // solve. x and rhs are PETSc vectors. unkns will be the final result in the end
    solver.solve(unks, rhs);    
    // get number of iterations
    int         iters = solver. getIters();    
    // get residual at specific iteration
    double      res = solver.getResidual(iters);    
    // finalization
    solver.finalize();
 
    // other codes
    ....

    return 0;
}
```


Note, this wrapper is specifically designed for our CFD solver -- **[PetIBM](https://github.com/barbagroup/PetIBM)**, so it may lack some features.  We are trying to make it more general.

## Update: 

### 01/08/2016

Now the wrapper has better performance on the cases that MPI ranks are more than GPUs.

It has poor performance when we have multiple MPI ranks calling AmgX solvers on one GPU. 
So, in order to keep only one MPI rank on each one GPU, this wrapper will now only call AmgX solvers on part of MPI ranks when there are more MPI ranks than GPUs. Data on other ranks will be transferred to those launching AmgX solvers.

For example, 
if there are totally 18 MPI ranks and 6 GPUs (while 12 MPI ranks and 2 GPUs on node 1, and 6 MPI ranks and 4 GPUs on node 2), 
only the rank 0, 6, 12, 14, 16, 17 will call AmgX solvers.
Data on rank 1-5 will be transfer to rank 0; data on rank 7-11 will go to rank 6; 13 to 12; and 15 to 14.
This causes some penalties on computing time because of data transfer between MPI ranks.
However, the overall performance is much better than calling AmgX solvers on ALL ranks directly.

The criterion to redistribute data is based on an assumption that original data are equally distributed on each rank, i.e., each rank holds similar number of rows of a parallel matrix. It's not necessary to have row indices continuous on successive ranks. 

## Future work:
* Add support of updating entries of a matrix A which is already on CUDA device
* Support other matrix structures other than AIJ
* Add mechanisms for debugging and error handling.

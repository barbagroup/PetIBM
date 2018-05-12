## Test
--------

To test the code, you can use the [Poisson example](../example/poisson). This
example solves a Poisson equation which we already know the exact solution.
Running this example, you will see the error norms of the numerical solution
against the exact solution. By doing grid convergence tests, you can also check
the order of convergence. If you want to solve other Poisson equations, you are
free to modify the hard-coded equation in the source code of functions
`generateRHS` and `generateExt`.

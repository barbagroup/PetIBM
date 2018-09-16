/* A simple coding example for using PetIBM API to build a Navier-Stokes solver.
 *
 * The numerical method used is based on the following paper:
 *
 * Perot, J. B. (1993). An analysis of the fractional step method. Journal of
 * Computational Physics, 108(1), 51-58.
 *
 * The discretization method is central finite difference and staggered grids.
 *
 * \file examples/api_examples/liddrivencavity2d/main.cpp
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

// PetIBM relies on YAML to simulation configurations, such as mesh settings,
// time step size, boundary conditions, etc. So we first include YAML header.
#include <yaml-cpp/yaml.h>

// End-users of a solver can provide simulation configurations through a YAML
// text file, so we need some parsing functions to parse the content of that
// text file, and convert then to a `YAML::Node` in our code. A `YAML::Node` is
// something similar to a `std::map` in C++ or a dictionary in Python. The
// required parsing functions are in the header `parser.h` and under the
// namespace `parser`.
#include <petibm/parser.h>

// The header `mesh.h` defines classes and functions regarding mesh objects
// under the namespace `mesh`. And it also defines a data type called `Mesh` to
// unify the member function interface for all different kinds of meshes. This
// data type, `Mesh`, is under the namespace `type`. In other words, with
// `type::Mesh`, high-level CFD code developers don't have to worry about what's
// the underlying mesh. A factory function, `mesh::createMesh`, can create an
// instance of `type::Mesh` with correct underlying mesh types according to the
// simulation configurations provided by end-users.
#include <petibm/mesh.h>

// The header `boundary.h` defines classes and functions regarding boundary
// conditions and ghost points. A data type called `Boundary` is defined under
// the namespace `type`. An instance of `type::Boundary` can help developers
// to call its member functions in a collective way under parallel environments.
#include <petibm/boundary.h>

// This header defines classes and functions regarding objects holding
// solutions. A data type called `Solution` is defined under `type` namespace to
// unify the interface to member functions of different kinds of solution
// structure.
#include <petibm/solution.h>

// This header defines functions to create different kinds of operators (i.e.,
// matrices). In PetIBM, we try to use linear algebra to express calculations.
// So even the non-linear operators, such as the convection term, are designed
// as matrix-free operators (matrix-free matrices).
#include <petibm/operators.h>

// `linsolver.h` defines classes and functions for linear solvers. It also
// defines a top-level data type `LinSolver` under the namespace `type` to unify
// the interface of different kinds of underlying linear solvers. The factory
// function `createLinSolver` will create the correct underlying linear solver
// for a `type::LinSolver` instance according to the simulation configurations.
#include <petibm/linsolver.h>

int main(int argc, char **argv)
{
    // We will use this variable to hold returned error code for all functions.
    PetscErrorCode ierr;

    // Initilaize MPI and PETSc
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    /*
     * Step 0: pre-define some parameters.
     *
     * Though all parameters should be able to be specified by the YAML text
     * file provided by end-users and hence are available in the `YAML::Node`
     * object later, we hard-coded them here to simplify this example code.
     */

    // Time marching parameters
    PetscInt start = 0;   // starting step of time marching
    PetscInt end = 2000;  // end step of time marching
    PetscReal dt = 0.01;  // time-step size
    PetscReal t = 0.0;    // beginning time value

    // flow parameter
    PetscReal nu = 0.01;  // flow viscosity

    /*
     * Step 1: read simulation configurations through a YAML text file.
     */

    // A YAML::Node to hold simulation configurations provided by end-users
    YAML::Node config;

    // Convert the content of simulation configurations from a text file to a
    // YAML::Node. The function `getSettings` will read and parse a text file
    // called "config.yaml" in the working directory, which by default is "./".
    // The working directory can be changed to other location through a command-
    // line argument `-directory` when end-users launch the solver. But in this
    // example, we will not use that for simplicity.
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    /*
     * Step 2: create a mesh object according to the configurations in the
     *         `config` object.
     */

    // Declare the mesh object
    petibm::type::Mesh mesh;

    // Create an instance of `type::Mesh` with correct underlying mesh type
    // according to the simulation parameters.
    ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh);
    CHKERRQ(ierr);

    // Write out the coordinates of grid points to a HDF file (grid.h5) under
    // the current directory.
    ierr = mesh->write("./grid"); CHKERRQ(ierr);

    /*
     * Step 3: create an object to hold solutions according the mesh.
     */

    // Declare the solution object
    petibm::type::Solution soln;

    ierr = petibm::solution::createSolution(mesh, soln); CHKERRQ(ierr);

    // Set the solution to be the initial conditions at this point. The initial
    // condition should exist in the simulation configurations.
    ierr = soln->applyIC(config); CHKERRQ(ierr);

    /*
     * Step 4: create a boundary object to hold all boundary conditions and
     *         information of ghost points according to the simulation
     *         configurations.
     */

    // Declare the boundary object
    petibm::type::Boundary bd;

    // Create the boundary object
    ierr = petibm::boundary::createBoundary(mesh, config, bd); CHKERRQ(ierr);

    // Apply initial conditions to the ghost points. This implies the current
    // state of the Solution object should be the initial conditions.
    ierr = bd->setGhostICs(soln); CHKERRQ(ierr);

    /*
     * Step 5: create operators (matrices).
     */

    // Declare operators (matrices)
    Mat L;     // Laplacian operator for velocity points
    Mat Lbc;   // Laplacian correction operator due to boundary conditions
    Mat G;     // Gradient (gradient of pressure to velocity points)
    Mat D;     // Divergence (divergence of velocity to pressure points)
    Mat Dbc;   // Divergence correction operator due to boundary conditions
    Mat N;     // Convection operator (non-linear term; matrix-free matrix)
    Mat A;     // Matrix resulting from implicit temporal treatment
    Mat BN;    // An approximation of inverse A matrix
    Mat BNG;   // Projection operator
    Mat DBNG;  // Poisson matrix for pressure points

    // create divergence operator
    ierr = petibm::operators::createDivergence(mesh, bd, D, Dbc, PETSC_FALSE);
    CHKERRQ(ierr);

    // create gradient operator
    ierr = petibm::operators::createGradient(mesh, G, PETSC_FALSE);
    CHKERRQ(ierr);

    // create Laplacian operator
    ierr = petibm::operators::createLaplacian(mesh, bd, L, Lbc); CHKERRQ(ierr);

    // create convection operator
    ierr = petibm::operators::createConvection(mesh, bd, N);

    // create the combined operator A; we use PETSc's function here
    ierr = MatDuplicate(L, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
    ierr = MatScale(A, -nu); CHKERRQ(ierr);
    ierr = MatShift(A, 1.0 / dt); CHKERRQ(ierr);

    // create the 1st-order approximated inverse A
    ierr = petibm::operators::createBnHead(L, dt, nu, 1, BN); CHKERRQ(ierr);

    // create BNG (= BN * G)
    ierr = MatMatMult(BN, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNG);
    CHKERRQ(ierr);

    // create DBNG (= D * BNG)
    ierr = MatMatMult(D, BNG, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DBNG);
    CHKERRQ(ierr);

    // DBNG is a singular matrix in this case, so we utilize PETSc's NullSpace
    {
        MatNullSpace nsp;  // a null space object from PETSc
        ierr =
            MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, nullptr, &nsp);
        CHKERRQ(ierr);
        ierr = MatSetNullSpace(DBNG, nsp); CHKERRQ(ierr);
        ierr = MatSetNearNullSpace(DBNG, nsp); CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(&nsp); CHKERRQ(ierr);
    }

    /*
     * Step 6: create vectors for holding values other than velocity and
     * pressure
     */

    // Declare vectors
    Vec dP;    // pressure increment
    Vec bc1;   // boundary correction values for velocity (resulting from Lbc *
               // velocity)
    Vec rhs1;  // right-hand side for velocity system
    Vec rhs2;  // right-hand side for Poisson system

    // create dp; it has the same structure as the pressure, so we just copy the
    // structure
    ierr = VecDuplicate(soln->pGlobal, &dP); CHKERRQ(ierr);

    // create bc1; it has the same structure as the velocity, so we just copy
    // the structure
    ierr = VecDuplicate(soln->UGlobal, &bc1); CHKERRQ(ierr);

    // create rhs1; it has the same structure as the velocity, so we just copy
    // the structure
    ierr = VecDuplicate(soln->UGlobal, &rhs1); CHKERRQ(ierr);

    // create rhs2; it has the same structure as the pressure, so we just copy
    // the structure
    ierr = VecDuplicate(soln->pGlobal, &rhs2); CHKERRQ(ierr);

    /*
     * Step 7: create linear solver objects. End-users can specify the type of
     *         linear solver to use, e.g. PETSc (CPU) or AmgX (GPU) in the
     *         YAML text file they provided.
     */

    // Declare solver objects
    petibm::type::LinSolver vSolver;  // linear solver for velocity system
    petibm::type::LinSolver pSolver;  // linear solver for Poisson system

    // create vSolver
    ierr = petibm::linsolver::createLinSolver("velocity", config, vSolver);
    CHKERRQ(ierr);

    // bind the coefficient matrix to vSolver
    ierr = vSolver->setMatrix(A); CHKERRQ(ierr);

    // create pSolver
    ierr = petibm::linsolver::createLinSolver("poisson", config, pSolver);
    CHKERRQ(ierr);

    // bind the coefficient matrix to pSolver
    ierr = pSolver->setMatrix(DBNG); CHKERRQ(ierr);

    /*
     * Step 8: time marching.
     *
     * To simplify the code, we use simple 1st explicit Euler method for
     * convection term and 1st implicit Euler for diffusion term. For more
     * complicated temporal integration, developers can use a data type
     * `type::TimeIntegration` from PetIBM. Here we omit this because we try
     * to keep this example simple.
     */
    for (int iter = start; iter <= end; iter++)
    {
        /*
         * Step 8-1: prepare the linear system for the velocity
         *
         * rhs1 = \frac{u^n}{\Delta t} - Gradient(p^n) - Convection(u^n) +
         *      (correction from the implicit diffusion term).
         */
        {
            // Initialize rhs1 as the pressure gradient on velocity points
            ierr = MatMult(G, soln->pGlobal, rhs1); CHKERRQ(ierr);

            // Add convection term
            ierr = MatMultAdd(N, soln->UGlobal, rhs1, rhs1);

            // Add explicit temporal derivation term, \frac{u^n}{\Delta t}
            ierr = VecAXPBY(rhs1, 1.0 / dt, -1.0, soln->UGlobal);
            CHKERRQ(ierr);

            // Add boundary condition correction of Laplacian at time n+1
            {
                // Update the equation between ghost points and solution points
                // based on the current state of solution
                ierr = bd->updateEqs(soln, dt); CHKERRQ(ierr);

                // Use the current state of solution to calculate the correction
                // values at the next time step
                ierr = MatMult(Lbc, soln->UGlobal, bc1); CHKERRQ(ierr);

                // Add the correction term to rhs1. We explicitly write down
                // 1.0 here is to remind that this correction is coming from
                // the implicit treatment of diffusion term in implicit 1st
                // Euler method. So 1.0 is the implicit coefficient of the 1st
                // implicit Euler method.
                ierr = VecAXPY(rhs1, 1.0 * nu, bc1); CHKERRQ(ierr);
            }
        }

        /*
         * Step 8-2: solve the velocity system, A \times u^{*) = rhs1
         */
        ierr = vSolver->solve(soln->UGlobal, rhs1); CHKERRQ(ierr);

        /*
         * Step 8-2: prepare the linear system for the Poisson system
         */
        {
            // Initialize rhs2 with divergence of velocity on pressure points
            ierr = MatMult(D, soln->UGlobal, rhs2); CHKERRQ(ierr);

            // Add the correction of divergence due to boundary conditions.
            ierr = MatMultAdd(Dbc, soln->UGlobal, rhs2, rhs2); CHKERRQ(ierr);
        }

        /*
         * Step 8-3: solve the Poisson system, DBNG \times dP = rhs2
         */
        ierr = pSolver->solve(dP, rhs2); CHKERRQ(ierr);

        /*
         * Step 8-4: projection step
         */
        {
            // Re-use rhs1 temporarily to hold BNG \times dP
            ierr = MatMult(BNG, dP, rhs1); CHKERRQ(ierr);

            // Correct the velocity to make it divergence free
            ierr = VecAXPY(soln->UGlobal, -1.0, rhs1); CHKERRQ(ierr);

            // Add pressure increment to the pressure
            ierr = VecAXPY(soln->pGlobal, 1.0, dP); CHKERRQ(ierr);
        }

        /*
         * Step 8-5: update the values of ghost points based on the equation
         *           between ghost points and solution points using the current
         *           state of solution.
         */
        ierr = bd->updateGhostValues(soln); CHKERRQ(ierr);
    }

    /*
     * Step 9: output solution to a HDF5 file (./solution/0002000.h5).
     */
    ierr = soln->write("./solution/0002000"); CHKERRQ(ierr);

    /*
     * Step 10: clean up objects to avoid memory leak.
     */
    ierr = pSolver->destroy(); CHKERRQ(ierr);
    ierr = vSolver->destroy(); CHKERRQ(ierr);

    ierr = VecDestroy(&rhs2); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs1); CHKERRQ(ierr);
    ierr = VecDestroy(&bc1); CHKERRQ(ierr);
    ierr = VecDestroy(&dP); CHKERRQ(ierr);
    ierr = MatDestroy(&DBNG); CHKERRQ(ierr);
    ierr = MatDestroy(&BNG); CHKERRQ(ierr);
    ierr = MatDestroy(&BN); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = MatDestroy(&N); CHKERRQ(ierr);
    ierr = MatDestroy(&Dbc); CHKERRQ(ierr);
    ierr = MatDestroy(&D); CHKERRQ(ierr);
    ierr = MatDestroy(&G); CHKERRQ(ierr);
    ierr = MatDestroy(&Lbc); CHKERRQ(ierr);
    ierr = MatDestroy(&L); CHKERRQ(ierr);

    ierr = bd->destroy(); CHKERRQ(ierr);
    ierr = soln->destroy(); CHKERRQ(ierr);
    ierr = mesh->destroy(); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}

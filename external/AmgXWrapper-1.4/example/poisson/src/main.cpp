/**
 * \file main.cpp
 * \brief An example and benchmark of AmgX and PETSc with Poisson system.
 *
 * This example solves a Poisson equation in 2D and 3D. We use standard central
 * difference with uniform grid within each direction. 

 * The Poisson equation we solve here is
 *      \nabla^2 u(x, y) = -8\pi^2 \cos{2\pi x} \cos{2\pi y}
 * for 2D. And
 *      \nabla^2 u(x, y, z) = -12\pi^2 \cos{2\pi x} \cos{2\pi y} \cos{2\pi z}
 * for 3D.
 *
 * The exact solutions are
 *      u(x, y) = \cos{2\pi x} \cos{2\pi y}
 * for 2D. And
 *      u(x, y, z) = \cos{2\pi x} \cos{2\pi y} \cos{2\pi z}
 * for 3D.
 *
 * The domain size is Lx = Ly = Lz = 1.
 * 
 * The boundary condition is all-Neumann BC, except that we pin a point as a
 * reference point (i.e., apply Dirichlet BC to that point) to avoid singular
 * matrix. We choose the point represented by the first row in matrix A as our 
 * reference point.
 *
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2015-02-01
 */


// STL
# include <cstring>

// PETSc
# include <petscsys.h>
# include <petscmat.h>
# include <petscvec.h>
# include <petscksp.h>

// AmgXWrapper
# include <AmgXSolver.hpp>

// headers
# include "StructArgs.hpp"
# include "io.hpp"
# include "factories.hpp"
# include "fixSingularMat.hpp"
# include "createKSP.hpp"
# include "solve.hpp"


int main(int argc, char **argv)
{
    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    StructArgs          args;   // a structure containing CMD arguments

    PetscScalar         Lx = 1.0, // Lx
                        Ly = 1.0, // Ly
                        Lz = 1.0; // Lz

    PetscScalar         dx,     // dx, calculated using Lx=1.0
                        dy,     // dy, calculated using Ly=1.0
                        dz;     // dy, calculated using Ly=1.0

    DM                  grid;   // DM object

    Vec                 x,      // x-coordinates
                        y,      // y-coordinates
                        z;      // z-coordinates

    Vec                 lhs,      // unknowns
                        rhs,    // RHS
                        u_exact,// exact solution
                        err;    // errors

    Mat                 A;      // coefficient matrix

    KSP                 ksp;    // PETSc KSP solver instance

    AmgXSolver          amgx;   // AmgX wrapper instance

    PetscMPIInt         size,   // MPI size
                        myRank; // rank of current process

    PetscClassId        solvingID,
                        warmUpID;

    PetscLogEvent       solvingEvent,
                        warmUpEvent;





    // initialize PETSc and MPI
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);


    // obtain the rank and size of MPI
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &myRank); CHKERRQ(ierr);


    // get parameters from command-line arguments
    ierr = args.getArgs(); CHKERRQ(ierr);


    // pring case information
    ierr = printHeader(args); CHKERRQ(ierr);

    // create DMDA object
    if (args.Nz > 0)
    {
        ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                DMDA_STENCIL_STAR, 
                args.Nx, args.Ny, args.Nz,
                PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                1, 1, nullptr, nullptr, nullptr, &grid); CHKERRQ(ierr);

        ierr = DMSetUp(grid); CHKERRQ(ierr);
        
        ierr = DMDASetUniformCoordinates(grid, 0., Lx, 0., Ly, 0., Lz); CHKERRQ(ierr);
    }
    else if (args.Nz == 0)
    {
        ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, 
                DMDA_STENCIL_STAR, 
                args.Nx, args.Ny,
                PETSC_DECIDE, PETSC_DECIDE,
                1, 1, nullptr, nullptr, &grid); CHKERRQ(ierr);

        ierr = DMSetUp(grid); CHKERRQ(ierr);

        ierr = DMDASetUniformCoordinates(grid, 0., Lx, 0., Ly, 0., 0.); CHKERRQ(ierr);
    }
    else
    {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "The number of cells in z direction can not be negetive. Your "
                "input is %d\n", args.Nz);
    }
    ierr = DMSetMatType(grid, MATAIJ); CHKERRQ(ierr);

            

    // create vectors (x, y, p, b, u) and matrix A
    {
        ierr = DMCreateGlobalVector(grid, &x); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(grid, &y); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(grid, &lhs); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(grid, &rhs); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(grid, &u_exact); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(grid, &err); CHKERRQ(ierr);
        ierr = DMCreateMatrix(grid, &A); CHKERRQ(ierr);

        ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);    CHKERRQ(ierr);
        ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);  CHKERRQ(ierr);
    }



    // setup the system
    if (args.Nz == 0) // 2D
    {
        // set values of dx, dy, dx, x, y, and z
        ierr = generateGrid(grid, args.Nx, args.Ny, Lx, Ly, dx, dy, x, y); CHKERRQ(ierr);

        // set values of RHS -- the vector rhs
        ierr = generateRHS(grid, x, y, rhs); CHKERRQ(ierr);

        // generate exact solution
        ierr = generateExt(grid, x, y, u_exact); CHKERRQ(ierr);

        // initialize and set up the coefficient matrix A
        ierr = generateA(grid, dx, dy, A); CHKERRQ(ierr);
    }
    else if (args.Nz > 0) // 3D
    {
        // create a PETSc Vec for z direction
        ierr = DMCreateGlobalVector(grid, &z); CHKERRQ(ierr);

        // set values of dx, dy, dx, x, y, and z
        ierr = generateGrid(grid, args.Nx, args.Ny, args.Nz,
                Lx, Ly, Lz, dx, dy, dz, x, y, z); CHKERRQ(ierr);

        // set values of RHS -- the vector rhs
        ierr = generateRHS(grid, x, y, z, rhs); CHKERRQ(ierr);

        // generate exact solution
        ierr = generateExt(grid, x, y, z, u_exact); CHKERRQ(ierr);

        // initialize and set up the coefficient matrix A
        ierr = generateA(grid, dx, dy, dz, A); CHKERRQ(ierr);
    }



    // pin a point with Dirichlet BC to resolve sinular mat due to all-Neumann BC
    ierr = fixSingularMat(A, rhs, u_exact); CHKERRQ(ierr);



    // register a PETSc event for warm-up and solving
    ierr = PetscClassIdRegister("SolvingClass", &solvingID); CHKERRQ(ierr);
    ierr = PetscClassIdRegister("WarmUpClass", &warmUpID); CHKERRQ(ierr);
    ierr = PetscLogEventRegister("Solving", solvingID, &solvingEvent); CHKERRQ(ierr);
    ierr = PetscLogEventRegister("WarmUp", warmUpID, &warmUpEvent); CHKERRQ(ierr);



    // create a solver and solve based whether it is AmgX or PETSc
    if (std::strcmp(args.mode, "PETSc") == 0) // PETSc mode
    {
        ierr = createKSP(ksp, A, grid, args.cfgFileName); CHKERRQ(ierr);

        ierr = solve(ksp, A, lhs, rhs, u_exact, err, 
                args, warmUpEvent, solvingEvent); CHKERRQ(ierr);

        // destroy KSP
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    }
    else // AmgX mode
    { 
        // AmgX GPU mode
        if (std::strcmp(args.mode, "AmgX_GPU") == 0)
        {
            ierr = amgx.initialize(PETSC_COMM_WORLD, "dDDI", args.cfgFileName);
            CHKERRQ(ierr);
        }
        // AmgX CPU mode (not supported yet)
        //else if (std::strcmp(args.mode, "AmgX_CPU") == 0)
        //    amgx.initialize(PETSC_COMM_WORLD, "hDDI", args.cfgFileName);
        else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE,
            "Invalid mode: %s\n", args.mode); CHKERRQ(ierr);


        ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = amgx.setA(A); CHKERRQ(ierr);

        ierr = solve(amgx, A, lhs, rhs, u_exact, err, 
                args, warmUpEvent, solvingEvent); CHKERRQ(ierr);

        // destroy solver
        ierr = amgx.finalize(); CHKERRQ(ierr);
    }


    // output a file for petsc performance
    if (args.optFileBool == PETSC_TRUE)
    {
        PetscViewer         viewer; // PETSc viewer

        std::strcat(args.optFileName ,".log");

        ierr = PetscViewerASCIIOpen(
                PETSC_COMM_WORLD, args.optFileName, &viewer); CHKERRQ(ierr);
        ierr = PetscLogView(viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    }
    

    // destroy vectors, matrix, dmda
    {
        ierr = VecDestroy(&x); CHKERRQ(ierr);
        ierr = VecDestroy(&y); CHKERRQ(ierr);
        ierr = VecDestroy(&lhs); CHKERRQ(ierr);
        ierr = VecDestroy(&rhs); CHKERRQ(ierr);
        ierr = VecDestroy(&u_exact); CHKERRQ(ierr);
        ierr = VecDestroy(&err); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        if (args.Nz > 0) {ierr = VecDestroy(&z); CHKERRQ(ierr); }

        ierr = DMDestroy(&grid); CHKERRQ(ierr);
    }


    // print a footer
    ierr = printFooter(args); CHKERRQ(ierr);

    // finalize PETSc
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}


/**
 * @file Poisson.cpp
 * @brief An example and benchmark of AmgX and PETSc
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version beta
 * @date 2015-02-01
 */
# include "headers.hpp"

static std::string help = "Test PETSc and AmgX solvers.";

int main(int argc, char **argv)
{
    StructArgs          args;   // a structure containing CMD arguments

    PetscScalar           Lx = 1.0, // Lx
                        Ly = 1.0, // Ly
                        Lz = 1.0; // Lz

    PetscScalar           dx,     // dx, calculated using Lx=1.0
                        dy,     // dy, calculated using Ly=1.0
                        dz;     // dy, calculated using Ly=1.0

    DM                  grid;   // DM object

    Vec                 x,      // x-coordinates
                        y,      // y-coordinates
                        z;      // z-coordinates

    Vec                 u,      // unknowns
                        rhs,    // RHS
                        bc,     // boundary conditions
                        u_exact,// exact solution
                        err;    // errors

    Mat                 A;      // coefficient matrix

    KSP                 ksp;    // PETSc KSP solver instance

    AmgXSolver          amgx;   // AmgX wrapper instance

    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    PetscMPIInt         size,   // MPI size
                        myRank; // rank of current process

    PetscClassId        solvingID,
                        warmUpID;

    PetscLogEvent       solvingEvent,
                        warmUpEvent;





    // initialize PETSc and MPI
    ierr = PetscInitialize(&argc, &argv, nullptr, help.c_str());             CHK;
    ierr = PetscLogDefaultBegin();                                           CHK;


    // obtain the rank and size of MPI
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);                           CHK;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &myRank);                         CHK;


    // get parameters from command-line arguments
    ierr = args.getArgs();                                                   CHK;


    // pring case information
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;
        for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "=");
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;
        for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;

        if (myRank == 0) args.print();
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;

        for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;
        for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "=");
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;

    }

    // create DMDA object
    if (args.Nz > 0)
    {
        ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                DMDA_STENCIL_STAR, 
                args.Nx, args.Ny, args.Nz,
                PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                1, 1, nullptr, nullptr, nullptr, &grid);                     CHK;

        ierr = DMDASetUniformCoordinates(grid, 0., 1., 0., 1., 0., 1.);      CHK;
    }
    else if (args.Nz == 0)
    {
        ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, 
                DMDA_STENCIL_STAR, 
                args.Nx, args.Ny,
                PETSC_DECIDE, PETSC_DECIDE,
                1, 1, nullptr, nullptr, &grid);                              CHK;

        ierr = DMDASetUniformCoordinates(grid, 0., 1., 0., 1., 0., 0.);      CHK;
    }

    ierr = PetscObjectSetName((PetscObject) grid, "DMDA grid");              CHK;
    ierr = DMSetMatType(grid, MATAIJ);                                       CHK;

            

    // create vectors (x, y, p, b, u) and matrix A
    {
        ierr = DMCreateGlobalVector(grid, &x);                               CHK;
        ierr = DMCreateGlobalVector(grid, &y);                               CHK;
        ierr = DMCreateGlobalVector(grid, &u);                               CHK;
        ierr = DMCreateGlobalVector(grid, &rhs);                             CHK;
        ierr = DMCreateGlobalVector(grid, &bc);                              CHK;
        ierr = DMCreateGlobalVector(grid, &u_exact);                         CHK;
        ierr = DMCreateGlobalVector(grid, &err);                             CHK;
        ierr = DMCreateMatrix(grid, &A);                                     CHK;

        ierr = PetscObjectSetName((PetscObject) x, "x coordinates");         CHK;
        ierr = PetscObjectSetName((PetscObject) y, "y coordinates");         CHK;
        ierr = PetscObjectSetName((PetscObject) u, "vec for unknowns");      CHK;
        ierr = PetscObjectSetName((PetscObject) rhs, "RHS");                 CHK;
        ierr = PetscObjectSetName((PetscObject) bc, "boundary conditions");  CHK;
        ierr = PetscObjectSetName((PetscObject) u_exact, "exact solution");  CHK;
        ierr = PetscObjectSetName((PetscObject) err, "errors");              CHK;
        ierr = PetscObjectSetName((PetscObject) A, "matrix A");              CHK;

        ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);    CHK;
        ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);  CHK;
    }




    if (args.Nz == 0)
    {
        // set values of dx, dy, dx, x, y, and z
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        ierr = generateGrid(grid, args.Nx, args.Ny, 
                Lx, Ly, dx, dy, x, y);                                       CHK;

        // set values of RHS -- the vector rhs
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        ierr = generateRHS(grid, x, y, rhs);                                 CHK;

        // generate exact solution
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        ierr = generateExt(grid, x, y, u_exact);                             CHK;

        // initialize and set up the coefficient matrix A
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        ierr = generateA(grid, dx, dy, A);                                   CHK;
    }
    else if (args.Nz > 0)
    {

        ierr = DMCreateGlobalVector(grid, &z);                               CHK;
        ierr = PetscObjectSetName((PetscObject) z, "z coordinates");         CHK;

        // set values of dx, dy, dx, x, y, and z
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        ierr = generateGrid(grid, args.Nx, args.Ny, args.Nz,
                Lx, Ly, Lz, dx, dy, dz, x, y, z);                            CHK;

        // set values of RHS -- the vector rhs
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        ierr = generateRHS(grid, x, y, z, rhs);                              CHK;

        // generate exact solution
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        ierr = generateExt(grid, x, y, z, u_exact);                          CHK;

        // initialize and set up the coefficient matrix A
        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        ierr = generateA(grid, dx, dy, dz, A);                               CHK;
    }


    // handle the issue for all-Neumann BC matrix
    ierr = MPI_Barrier(PETSC_COMM_WORLD);                                    CHK;
    ierr = applyNeumannBC(A, rhs, u_exact);                                  CHK;


    // register a PETSc event for warm-up and solving
    ierr = PetscClassIdRegister("SolvingClass", &solvingID);                 CHK;
    ierr = PetscClassIdRegister("WarmUpClass", &warmUpID);                   CHK;
    ierr = PetscLogEventRegister("Solving", solvingID, &solvingEvent);       CHK;
    ierr = PetscLogEventRegister("WarmUp", warmUpID, &warmUpEvent);          CHK;

    // create a solver and solve based whether it is AmgX or PETSc
    if (std::strcmp(args.mode, "PETSc") == 0) // PETSc mode
    {
        ierr = createKSP(ksp, A, grid, args.cfgFileName);                    CHK;

        ierr = solve(ksp, A, u, rhs, u_exact, err, 
                args, warmUpEvent, solvingEvent);                            CHK;

        // destroy KSP
        ierr = KSPDestroy(&ksp);                                             CHK;

    }
    else // AmgX mode
    {
        if (std::strcmp(args.mode, "AmgX") == 0) // AmgX GPU mode
            amgx.initialize(PETSC_COMM_WORLD, "dDDI", args.cfgFileName);
        else // AmgX CPU mode (not yet implemented in the wrapper) and other mode
        {   
            std::cerr << "Invalid mode." << std::endl;
            exit(EXIT_FAILURE);
        }


        ierr = MPI_Barrier(PETSC_COMM_WORLD);                                CHK;
        amgx.setA(A);

        ierr = solve(amgx, A, u, rhs, u_exact, err, 
                args, warmUpEvent, solvingEvent);                            CHK;

        // destroy solver
        ierr = amgx.finalize();                                              CHK;

    }


    // output a file for petsc performance
    if (args.optFileBool == PETSC_TRUE)
    {
        PetscViewer         viewer; // PETSc viewer

        std::strcat(args.optFileName ,".log");

        ierr = PetscViewerASCIIOpen(
                PETSC_COMM_WORLD, args.optFileName, &viewer);                CHK;
        ierr = PetscLogView(viewer);                                         CHK;
        ierr = PetscViewerDestroy(&viewer);                                  CHK;
    }
    
    // output VTK file for post-processing
    if (args.VTKFileBool == PETSC_TRUE)
    {
        PetscViewer         viewerVTK;

        std::strcat(args.VTKFileName ,".vts");

        ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewerVTK);              CHK;
        ierr = PetscViewerSetType(viewerVTK, PETSCVIEWERVTK);                CHK;
        ierr = PetscViewerSetFormat(viewerVTK, PETSC_VIEWER_VTK_VTS);        CHK;
        ierr = PetscViewerFileSetMode(viewerVTK, FILE_MODE_WRITE);           CHK;
        ierr = PetscViewerFileSetName(viewerVTK, args.VTKFileName);          CHK;

        ierr = VecView(u_exact, viewerVTK);                                  CHK;

        ierr = PetscViewerFileSetMode(viewerVTK, FILE_MODE_APPEND);          CHK;

        ierr = VecView(u, viewerVTK);                                        CHK;
        ierr = VecView(err, viewerVTK);                                      CHK;

        ierr = PetscViewerDestroy(&viewerVTK);                               CHK;
    }
    

    // destroy vectors, matrix, dmda
    {
        ierr = VecDestroy(&x);                                               CHK;
        ierr = VecDestroy(&y);                                               CHK;
        ierr = VecDestroy(&u);                                               CHK;
        ierr = VecDestroy(&rhs);                                             CHK;
        ierr = VecDestroy(&bc);                                              CHK;
        ierr = VecDestroy(&u_exact);                                         CHK;
        ierr = VecDestroy(&err);                                             CHK;

        ierr = MatDestroy(&A);                                               CHK;

        if (args.Nz > 0) 
        {
            ierr = VecDestroy(&z);                                           CHK;
        }

        ierr = DMDestroy(&grid);                                             CHK;
    }


    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;
        for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "=");
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;
        for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;
        ierr = PetscPrintf(PETSC_COMM_WORLD, 
                "End of %s\n", args.caseName);                               CHK;
        for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;
        for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "=");
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                          CHK;
    }

    // finalize PETSc
    ierr = PetscFinalize();                                                  CHK;

    return 0;
}


/**
 * \file init.cpp
 * \brief definition of some member functions of the class AmgXSolver.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2016-01-08
 */


// CUDA
# include <cuda_runtime.h>

// AmgXSolver
# include "AmgXSolver.hpp"


// definition of AmgXSolver::AmgXSolver
AmgXSolver::AmgXSolver(const MPI_Comm &comm,
        const std::string &modeStr, const std::string &cfgFile)
{
    initialize(comm, modeStr, cfgFile);
}


// definition of overload AmgXSolver::~AmgXSolver
AmgXSolver::~AmgXSolver()
{
    if (isInitialized) finalize();
}


// definition of AmgXSolver::initialize
PetscErrorCode AmgXSolver::initialize(const MPI_Comm &comm,
        const std::string &modeStr, const std::string &cfgFile)
{
    PetscErrorCode      ierr;

    PetscFunctionBeginUser;

    // if this instance has already been initialized, skip
    if (isInitialized) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
            "This AmgXSolver instance has been initialized on this process.");

    // increase the number of AmgXSolver instances
    count += 1;

    // get the name of this node
    int     len;
    char    name[MPI_MAX_PROCESSOR_NAME];
    ierr = MPI_Get_processor_name(name, &len); CHK;
    nodeName = name;

    // get the mode of AmgX solver
    ierr = setMode(modeStr); CHK;

    // initialize communicators and corresponding information
    ierr = initMPIcomms(comm); CHK;

    // only processes in gpuWorld are required to initialize AmgX
    if (gpuProc == 0)
    {
        ierr = initAmgX(cfgFile); CHK;
    }

    // a bool indicating if this instance is initialized
    isInitialized = true;

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::initMPIcomms
PetscErrorCode AmgXSolver::initMPIcomms(const MPI_Comm &comm)
{
    PetscErrorCode      ierr;

    PetscFunctionBeginUser;

    // duplicate the global communicator
    ierr = MPI_Comm_dup(comm, &globalCpuWorld); CHK;
    ierr = MPI_Comm_set_name(globalCpuWorld, "globalCpuWorld"); CHK;

    // get size and rank for global communicator
    ierr = MPI_Comm_size(globalCpuWorld, &globalSize); CHK;
    ierr = MPI_Comm_rank(globalCpuWorld, &myGlobalRank); CHK;


    // Get the communicator for processors on the same node (local world)
    ierr = MPI_Comm_split_type(globalCpuWorld, 
            MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &localCpuWorld); CHK;
    ierr = MPI_Comm_set_name(localCpuWorld, "localCpuWorld"); CHK;

    // get size and rank for local communicator
    ierr = MPI_Comm_size(localCpuWorld, &localSize); CHK;
    ierr = MPI_Comm_rank(localCpuWorld, &myLocalRank); CHK;


    // set up the variable nDevs
    ierr = setDeviceCount(); CHK;


    // set up corresponding ID of the device used by each local process
    ierr = setDeviceIDs(); CHK;
    ierr = MPI_Barrier(globalCpuWorld); CHK;


    // split the global world into a world involved in AmgX and a null world
    ierr = MPI_Comm_split(globalCpuWorld, gpuProc, 0, &gpuWorld); CHK;

    // get size and rank for the communicator corresponding to gpuWorld
    if (gpuWorld != MPI_COMM_NULL)
    {
        ierr = MPI_Comm_set_name(gpuWorld, "gpuWorld"); CHK;
        ierr = MPI_Comm_size(gpuWorld, &gpuWorldSize); CHK;
        ierr = MPI_Comm_rank(gpuWorld, &myGpuWorldRank); CHK;
    }
    else // for those can not communicate with GPU devices
    {
        gpuWorldSize = MPI_UNDEFINED;
        myGpuWorldRank = MPI_UNDEFINED;
    }
    

    // split local world into worlds corresponding to each CUDA device
    ierr = MPI_Comm_split(localCpuWorld, devID, 0, &devWorld); CHK;
    ierr = MPI_Comm_set_name(devWorld, "devWorld"); CHK;

    // get size and rank for the communicator corresponding to myWorld
    ierr = MPI_Comm_size(devWorld, &devWorldSize); CHK;
    ierr = MPI_Comm_rank(devWorld, &myDevWorldRank); CHK;

    ierr = MPI_Barrier(globalCpuWorld); CHK;

    return 0;
}


// definition of AmgXSolver::setDeviceCount
PetscErrorCode AmgXSolver::setDeviceCount()
{
    PetscFunctionBeginUser;

    // get the number of devices that AmgX solvers can use
    switch (mode)
    {
        case AMGX_mode_dDDI: // for GPU cases, nDevs is the # of local GPUs
        case AMGX_mode_dDFI: // for GPU cases, nDevs is the # of local GPUs
        case AMGX_mode_dFFI: // for GPU cases, nDevs is the # of local GPUs
            // get the number of total cuda devices
            CHECK(cudaGetDeviceCount(&nDevs));

            // Check whether there is at least one CUDA device on this node
            if (nDevs == 0) SETERRQ1(MPI_COMM_WORLD, PETSC_ERR_SUP_SYS,
                    "There is no CUDA device on the node %s !\n", nodeName.c_str());
            break;
        case AMGX_mode_hDDI: // for CPU cases, nDevs is the # of local processes
        case AMGX_mode_hDFI: // for CPU cases, nDevs is the # of local processes
        case AMGX_mode_hFFI: // for CPU cases, nDevs is the # of local processes
        default:
            nDevs = localSize;
            break;
    }

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::setDeviceIDs
PetscErrorCode AmgXSolver::setDeviceIDs()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // set the ID of device that each local process will use
    if (nDevs == localSize) // # of the devices and local precosses are the same
    {
        devID = myLocalRank;
        gpuProc = 0;
    }
    else if (nDevs > localSize) // there are more devices than processes
    {
        ierr = PetscPrintf(localCpuWorld, "CUDA devices on the node %s "
                "are more than the MPI processes launched. Only %d CUDA "
                "devices will be used.\n", nodeName.c_str(), localSize); CHK;

        devID = myLocalRank;
        gpuProc = 0;
    }
    else // there more processes than devices
    {
        int     nBasic = localSize / nDevs,
                nRemain = localSize % nDevs;

        if (myLocalRank < (nBasic+1)*nRemain)
        {
            devID = myLocalRank / (nBasic + 1);
            if (myLocalRank % (nBasic + 1) == 0)  gpuProc = 0;
        }
        else
        {
            devID = (myLocalRank - (nBasic+1)*nRemain) / nBasic + nRemain;
            if ((myLocalRank - (nBasic+1)*nRemain) % nBasic == 0) gpuProc = 0;
        }
    }

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::initAmgX
PetscErrorCode AmgXSolver::initAmgX(const std::string &cfgFile)
{
    PetscFunctionBeginUser;

    // only the first instance (AmgX solver) is in charge of initializing AmgX
    if (count == 1)
    {
        // initialize AmgX
        AMGX_SAFE_CALL(AMGX_initialize());

        // intialize AmgX plugings
        AMGX_SAFE_CALL(AMGX_initialize_plugins());

        // only the master process can output something on the screen
        AMGX_SAFE_CALL(AMGX_register_print_callback(
                    [](const char *msg, int length)->void
                    {PetscPrintf(PETSC_COMM_WORLD, "%s", msg);})); 

        // let AmgX to handle errors returned
        AMGX_SAFE_CALL(AMGX_install_signal_handler());
    }

    // create an AmgX configure object
    AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, cfgFile.c_str()));

    // let AmgX handle returned error codes internally
    AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));

    // create an AmgX resource object, only the first instance is in charge
    if (count == 1) AMGX_resources_create(&rsrc, cfg, &gpuWorld, 1, &devID);

    // create AmgX vector object for unknowns and RHS
    AMGX_vector_create(&AmgXP, rsrc, mode);
    AMGX_vector_create(&AmgXRHS, rsrc, mode);

    // create AmgX matrix object for unknowns and RHS
    AMGX_matrix_create(&AmgXA, rsrc, mode);

    // create an AmgX solver object
    AMGX_solver_create(&solver, rsrc, mode, cfg);

    // obtain the default number of rings based on current configuration
    AMGX_config_get_default_number_of_rings(cfg, &ring);

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::finalize
PetscErrorCode AmgXSolver::finalize()
{
    PetscErrorCode      ierr;

    PetscFunctionBeginUser;

    // skip if this instance has not been initialized
    if (! isInitialized)
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "This AmgXWrapper has not been initialized. "
                "Please initialize it before finalization.\n"); CHK;

        PetscFunctionReturn(0);
    }

    // only processes using GPU are required to destroy AmgX content
    if (gpuProc == 0)
    {
        // destroy solver instance
        AMGX_solver_destroy(solver);

        // destroy matrix instance
        AMGX_matrix_destroy(AmgXA);

        // destroy RHS and unknown vectors
        AMGX_vector_destroy(AmgXP);
        AMGX_vector_destroy(AmgXRHS);

        // only the last instance need to destroy resource and finalizing AmgX
        if (count == 1)
        {
            AMGX_resources_destroy(rsrc);
            AMGX_SAFE_CALL(AMGX_config_destroy(cfg));

            AMGX_SAFE_CALL(AMGX_finalize_plugins());
            AMGX_SAFE_CALL(AMGX_finalize());
        }
        else
        {
            AMGX_config_destroy(cfg);
        }

        // destroy gpuWorld
        ierr = MPI_Comm_free(&gpuWorld); CHK;
    }

    // destroy PETSc objects
    ierr = VecScatterDestroy(&scatterLhs); CHK;
    ierr = VecScatterDestroy(&scatterRhs); CHK;
    ierr = VecDestroy(&redistLhs); CHK;
    ierr = VecDestroy(&redistRhs); CHK;

    // re-set necessary variables in case users want to reuse 
    // the variable of this instance for a new instance
    gpuProc = MPI_UNDEFINED;
    ierr = MPI_Comm_free(&globalCpuWorld); CHK;
    ierr = MPI_Comm_free(&localCpuWorld); CHK;
    ierr = MPI_Comm_free(&devWorld); CHK;

    // decrease the number of instances
    count -= 1;

    // change status
    isInitialized = false;

    PetscFunctionReturn(0);
}

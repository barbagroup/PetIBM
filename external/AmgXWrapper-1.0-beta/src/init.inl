/**
 * @file init.inl
 * @brief Definition of some member functions of the class AmgXSolver
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2016-01-08
 */
# include "AmgXSolver.hpp"


/**
 * @brief Initialization of AmgXSolver instance
 *
 * @param comm MPI communicator for all processes
 * @param _mode The mode this solver will run in. Please refer to AmgX manual.
 * @param cfg_file A file containing the configurations of this solver
 *
 * @return Currently meaningless. May be error codes in the future.
 */
int AmgXSolver::initialize(MPI_Comm comm, 
        const std::string &_mode, const std::string &cfg_file)
{
    // get the number of total cuda devices
    CHECK(cudaGetDeviceCount(&nDevs));

    // Check whether there is at least one CUDA device on this node
    if (nDevs == 0) 
    {
        std::cerr << "There are no CUDA devices on the node " 
                  << nodeName << " !!" << std::endl;

        exit(EXIT_FAILURE);
    }

    // initialize other communicators
    initMPIcomms(comm);

    if (gpuProc == 0) initAmgX(_mode, cfg_file);

    return 0;
}


/**
 * @brief Initialize different MPI communicators
 *
 * This private initializes all communicators needed:
 *
 * 1. globalCpuWorld is the communicator for all cpu processors, which is 
 *    usually the same with MPI_COMM_WORLD
 * 2. localCpuWorld is the communicaotr for cpu processors located in the same 
 *    node, which can use CUDA devices on that node
 * 3. gpuWorld is the communicator for the cpu processors that involved in AmgX
 *    function calls, i.e., those who will call AmgX functions
 * 4. devWorld is the communicator for the cpu processors that share one CUDA
 *    device. The first processor in this communicator will also be in gpuWorld.
 *    And all other processors in this communicator will transfer their data to 
 *    that fisrt processor.
 * 
 * XXXXXXXSize and XXXXXXXRank represent the size and rank obtained from these
 * communicators.
 *
 * @param comm A communicator for all processes
 *
 * @return Currently meaningless. May be error codes in the future.
 */
int AmgXSolver::initMPIcomms(MPI_Comm &comm)
{
    // Get the name of this node
    MPI_Get_processor_name(nodeName, &nodeNameLen);

    // Copy communicator for all processors
    globalCpuWorld = comm;

    // get size and rank for global communicator
    MPI_Comm_size(globalCpuWorld, &globalSize);
    MPI_Comm_rank(globalCpuWorld, &myGlobalRank);

    // Get the communicator for processors on the same node (local world)
    MPI_Comm_split_type(globalCpuWorld, 
            MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &localCpuWorld);

    // get size and rank for local communicator
    MPI_Comm_size(localCpuWorld, &localSize);
    MPI_Comm_rank(localCpuWorld, &myLocalRank);


    if (nDevs == localSize)
    {
        devID = myLocalRank;
        gpuProc = 0;
    }
    else if (nDevs > localSize)
    {
        if (myLocalRank == 0)
            std::cout << "CUDA devices on the node " << nodeName 
                      << " are more than the MPI processes launched. " 
                      << "Only " << localSize << " CUDA devices will be used."
                      << std::endl;

        devID = myLocalRank;
        gpuProc = 0;
    }
    else
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

    MPI_Barrier(globalCpuWorld);

    // split the global world into a world involved in AmgX and a null world
    MPI_Comm_split(globalCpuWorld, gpuProc, 0, &gpuWorld);

    // get size and rank for the communicator corresponding to gpuWorld
    if (gpuWorld != MPI_COMM_NULL)
    {
        MPI_Comm_size(gpuWorld, &gpuWorldSize);
        MPI_Comm_rank(gpuWorld, &myGpuWorldRank);
    }
    else
    {
        gpuWorldSize = MPI_UNDEFINED;
        myGpuWorldRank = MPI_UNDEFINED;
    }
    

    // split local world into worlds corresponding to each CUDA device
    MPI_Comm_split(localCpuWorld, devID, 0, &devWorld);

    // get size and rank for the communicator corresponding to myWorld
    MPI_Comm_size(devWorld, &devWorldSize);
    MPI_Comm_rank(devWorld, &myDevWorldRank);

    MPI_Barrier(globalCpuWorld);

    return 0;
}



/**
 * @brief Initialize the AmgX library
 *
 * This function initializes the current instance (solver). Based on the count, 
 * only the instance initialized first is in charge of initializing AmgX and the 
 * resource instance.
 *
 * @param _mode The mode this solver will run in. Please refer to AmgX manual.
 * @param _cfg A file containing the configurations of this solver
 *
 * @return Currently meaningless. May be error codes in the future.
 */
int AmgXSolver::initAmgX(const std::string &_mode, const std::string &_cfg)
{
    count += 1;

    // only the first instance (AmgX solver) is in charge of initializing AmgX
    if (count == 1)
    {
        // initialize AmgX
        AMGX_SAFE_CALL(AMGX_initialize());

        // intialize AmgX plugings
        AMGX_SAFE_CALL(AMGX_initialize_plugins());

        // use user-defined output mechanism. only the master process can output
        // something on the screen
        if (myGpuWorldRank == 0) 
        { 
            AMGX_SAFE_CALL(
                AMGX_register_print_callback(&(AmgXSolver::print_callback))); 
        }
        else 
        { 
            AMGX_SAFE_CALL(
                AMGX_register_print_callback(&(AmgXSolver::print_none))); 
        }

        // let AmgX to handle errors returned
        AMGX_SAFE_CALL(AMGX_install_signal_handler());
    }

    // create an AmgX configure object
    AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, _cfg.c_str()));

    // let AmgX handle returned error codes internally
    AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));

    // create an AmgX resource object, only the first instance is in charge
    if (count == 1) AMGX_resources_create(&rsrc, cfg, &gpuWorld, 1, &devID);

    // set mode
    setMode(_mode);

    // create AmgX vector object for unknowns and RHS
    AMGX_vector_create(&AmgXP, rsrc, mode);
    AMGX_vector_create(&AmgXRHS, rsrc, mode);

    // create AmgX matrix object for unknowns and RHS
    AMGX_matrix_create(&AmgXA, rsrc, mode);

    // create an AmgX solver object
    AMGX_solver_create(&solver, rsrc, mode, cfg);

    // obtain the default number of rings based on current configuration
    AMGX_config_get_default_number_of_rings(cfg, &ring);

    isInitialized = true;

    return 0;
}


/**
 * @brief Finalizing the instance.
 *
 * This function destroys AmgX data. The instance last destroyed also needs to 
 * destroy shared resource instance and finalizing AmgX.
 *
 * @return Currently meaningless. May be error codes in the future.
 */
int AmgXSolver::finalize()
{

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

        // change status
        isInitialized = false;

        count -= 1;
    }

    return 0;
}


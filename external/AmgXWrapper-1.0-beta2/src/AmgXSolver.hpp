/**
 * @file AmgXSolver.hpp
 * @brief Declaration of the class AmgXSolver
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2015-09-01
 */

# pragma once

# include <iostream>
# include <string>
# include <cstring>
# include <functional>
# include <vector>
# include <cstdlib>

# include <cuda_runtime.h>
# include <amgx_c.h>

# include <petscmat.h>
# include <petscvec.h>
# include <petscis.h>

# include "check.hpp"

/**
 * @brief A wrapper class for an interface between PETSc and AmgX
 *
 * This class is a wrapper of AmgX library. PETSc user only need to pass 
 * PETSc matrix and vectors into the AmgXSolver instance to solve their problem.
 *
 */
class AmgXSolver
{
    public:

        /// default constructor
        AmgXSolver() = default;

        /// initialization of instance
        int initialize(MPI_Comm comm, 
                const std::string &_mode, const std::string &cfg_file);

        /// finalization
        int finalize();

        /// convert PETSc matrix into AmgX matrix and pass it to solver
        int setA(const Mat &A);

        /// solve the problem, soultion vector will be updated in the end
        int solve(Vec &p, Vec &b);

        /// Get the number of iterations of last solve phase
        int getIters();

        /// Get the residual at a specific iteration in last solve phase
        double getResidual(const int &iter);

        /// get the memory usage on device
        int getMemUsage();


    private:



        int                     nDevs,      /*< # of cuda devices*/
                                devID;

        int                     gpuProc = MPI_UNDEFINED;

        MPI_Comm                globalCpuWorld,
                                localCpuWorld,
                                gpuWorld,
                                devWorld;

        int                     globalSize,
                                localSize,
                                gpuWorldSize,
                                devWorldSize;

        int                     myGlobalRank,
                                myLocalRank,
                                myGpuWorldRank,
                                myDevWorldRank;

        int                     nodeNameLen;
        char                    nodeName[MPI_MAX_PROCESSOR_NAME];



        static int              count;      /*!< only one instance allowed*/
        int                     ring;       /*< a parameter used by AmgX*/

        AMGX_Mode               mode;               /*< AmgX mode*/
        AMGX_config_handle      cfg = nullptr;      /*< AmgX config object*/
        AMGX_matrix_handle      AmgXA = nullptr;    /*< AmgX coeff mat*/
        AMGX_vector_handle      AmgXP = nullptr,    /*< AmgX unknowns vec*/
                                AmgXRHS = nullptr;  /*< AmgX RHS vec*/
        AMGX_solver_handle      solver = nullptr;   /*< AmgX solver object*/
        static AMGX_resources_handle   rsrc;        /*< AmgX resource object*/



        bool                    isInitialized = false;  /*< as its name*/



        /// set up the mode of AmgX solver
        int setMode(const std::string &_mode);

        /// a printing function using stdout
        static void print_callback(const char *msg, int length);

        /// a printing function that prints nothing, used by AmgX
        static void print_none(const char *msg, int length);


        int initMPIcomms(MPI_Comm &comm);
        int initAmgX(const std::string &_mode, const std::string &_cfg);


        VecScatter      redistScatter = nullptr;
        Vec             redistLhs = nullptr,
                        redistRhs = nullptr;

        int getLocalMatRawData(Mat &localA, PetscInt &localN,
                std::vector<PetscInt> &row, std::vector<Petsc64bitInt> &col,
                std::vector<PetscScalar> &data);

        int getDevIS(const Mat &A, IS &devIS);
        int getLocalA(const Mat &A, const IS &devIS, Mat &localA);
        int redistMat(const Mat &A, const IS &devIS, Mat &newA);
        int getPartVec(
                const IS &devIS, const PetscInt &N, std::vector<PetscInt> &partVec);
        int destroyLocalA(const Mat &A, Mat &localA);
        int getVecScatter(const Mat &A1, const Mat &A2, const IS &devIS);

        int solve_real(Vec &p, Vec &b);

};

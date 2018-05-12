/**
 * \file AmgXSolver.hpp
 * \brief definition of class AmgXSolver.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2015-09-01
 */


# pragma once

// STL
# include <string>
# include <vector>

// AmgX
# include <amgx_c.h>

// PETSc
# include <petscmat.h>
# include <petscvec.h>


// for checking CUDA function call
# define CHECK(call)                                                        \
{                                                                           \
    const cudaError_t       error = call;                                   \
    if (error != cudaSuccess)                                               \
    {                                                                       \
        SETERRQ4(PETSC_COMM_WORLD, PETSC_ERR_SIG,                           \
            "Error: %s:%d, code:%d, reason: %s\n",                          \
            __FILE__, __LINE__, error, cudaGetErrorString(error));          \
    }                                                                       \
}


// abbreviation for check PETSc returned error code
# define CHK CHKERRQ(ierr)


/**
 * @brief A wrapper class for coupling PETSc and AmgX.
 *
 * This class is a wrapper of AmgX library for PETSc. PETSc users only need to
 * pass a PETSc matrix and vectors into the AmgXSolver instance to solve their
 * problems.
 *
 */
class AmgXSolver
{
    public:

        /** \brief default constructor. */
        AmgXSolver() = default;


        /**
         * \brief construct a AmgXSolver instance.
         *
         * \param comm [in] MPI communicator.
         * \param modeStr [in] a string, target mode of AmgX (e.g., dDDI).
         * \param cfgFile [in] a string indicate the path to AmgX configuration file.
         */
        AmgXSolver(const MPI_Comm &comm, 
                const std::string &modeStr, const std::string &cfgFile);


        /** \brief destructor. */
        ~AmgXSolver();


        /**
         * \brief initialize a AmgXSolver instance.
         *
         * \param comm [in] MPI communicator.
         * \param modeStr [in] a string, target mode of AmgX (e.g., dDDI).
         * \param cfgFile [in] a string indicate the path to AmgX configuration file.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode initialize(const MPI_Comm &comm,
                const std::string &modeStr, const std::string &cfgFile);


        /**
         * \brief Finalizing the instance.
         *
         * This function destroys AmgX data. The instance last destroyed is also 
         * in charge of destroying shared resource object and finalizing AmgX.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode finalize();


        /**
         * \brief set up the matrix used by AmgX.
         *
         * This function will automatically convert PETSc matrix to AmgX matrix.
         * If the AmgX solver is set to be the GPU one, we also redestribute the
         * matrix in this function and upload it to GPUs.
         *
         * Note: currently we can only handle AIJ format.
         *
         * \param A [in] a PETSc Mat.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode setA(const Mat &A);


        /**
         * \brief solve the linear system.
         *
         * `p` vector will be used as initial guess and will be updated to the 
         * solution in the end of solving.
         *
         * For cases that use more MPI processes than the number of GPUs, this 
         * function will do data gathering before solving and data scattering 
         * after the solving.
         *
         * \param p [in, out] a PETSc Vec object representing unknowns.
         * \param b [in] a PETSc Vec representing right hand side.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode solve(Vec &p, Vec &b);


        /**
         * \brief get the number of iterations of the last solving.
         *
         * \param iter [out] returned number of iterations.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode getIters(int &iter);


        /**
         * \brief get the residual at a specific iteration during the last solving.
         *
         * \param iter [in] the target iteration.
         * \param res [out] the returned residual.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode getResidual(const int &iter, double &res);


    private:

        /** \brief the current count of AmgXSolver instances. */
        static int              count;

        /** \brief a flag indicating if this instance has been initialized. */
        bool                    isInitialized = false;

        /** \brief the name of the node that this process belongs to. */
        std::string             nodeName;




        /** \brief number of local devices used by AmgX.*/
        PetscMPIInt             nDevs;

        /** \brief the ID of corresponding device used by this process. */
        PetscMPIInt             devID;

        /** \brief a flag indicating if this process will talk with GPU. */
        PetscMPIInt             gpuProc = MPI_UNDEFINED;

        /** \brief a communicator for global world. */
        MPI_Comm                globalCpuWorld;

        /** \brief a communicator for local world (i.e., in-node). */
        MPI_Comm                localCpuWorld;

        /** \brief a communicator for processes that can talk to GPUs. */
        MPI_Comm                gpuWorld;

        /** \brief a communicator for processes using the same devices. */
        MPI_Comm                devWorld;

        /** \brief size of the `globalCpuWorld`. */
        PetscMPIInt             globalSize;

        /** \brief size of the `localCpuWorld`. */
        PetscMPIInt             localSize;

        /** \brief size of the `gpuWorld`. */
        PetscMPIInt             gpuWorldSize;

        /** \brief size of the `devWorld`. */
        PetscMPIInt             devWorldSize;

        /** \brief rank of this process in the `globalCpuWorld`. */
        PetscMPIInt             myGlobalRank;

        /** \brief rank of this process in the `localCpuWorld`. */
        PetscMPIInt             myLocalRank;

        /** \brief rank of this process in the `gpuWorld`. */
        PetscMPIInt             myGpuWorldRank;

        /** \brief rank of this process in the `devWorld`. */
        PetscMPIInt             myDevWorldRank;




        /** \brief a parameter used by AmgX. */
        int                     ring;

        /** \brief AmgX solver mode. */
        AMGX_Mode               mode;

        /** \brief AmgX config object. */
        AMGX_config_handle      cfg = nullptr;

        /** \brief AmgX matrix object. */
        AMGX_matrix_handle      AmgXA = nullptr;

        /** \brief AmgX vector object representing unknowns. */
        AMGX_vector_handle      AmgXP = nullptr;

        /** \brief AmgX vector object representing RHS. */
        AMGX_vector_handle      AmgXRHS = nullptr;

        /** \brief AmgX solver object. */
        AMGX_solver_handle      solver = nullptr;

        /** \brief AmgX resource object. */
        static AMGX_resources_handle   rsrc;




        /** \brief a VecScatter for gathering/scattering between original PETSc 
         *         Vec and temporary redistributed PETSc Vec.*/
        VecScatter              scatterLhs = nullptr;

        /** \brief a VecScatter for gathering/scattering between original PETSc 
         *         Vec and temporary redistributed PETSc Vec.*/
        VecScatter              scatterRhs = nullptr;

        /** \brief a temporary PETSc Vec holding redistributed unknowns. */
        Vec                     redistLhs = nullptr;

        /** \brief a temporary PETSc Vec holding redistributed RHS. */
        Vec                     redistRhs = nullptr;




        /**
         * \brief set the AmgX solver mode from user-provided string.
         *
         * Available modes are: dDDI, dDFI, dFFI, hDDI, hDFI, hFFI.
         *
         * \param modeStr [in] a `std::string`.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode setMode(const std::string &modeStr);


        /**
         * \brief get the number of devices on thie node.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode setDeviceCount();


        /**
         * \brief set the ID of the corresponfing device used by this process.
         *
         * \return PetscErrorCode. 
         */
        PetscErrorCode setDeviceIDs();


        /**
         * \brief initialize different MPI communicators.
         *
         * The `comm` provided by users will be duplicated and saved to the 
         * member `globalCpuWorld`.
         *
         * \param comm [in] global communicator.
         *
         * \return PetscErrorCode. 
         */
        PetscErrorCode initMPIcomms(const MPI_Comm &comm);


        /**
         * \brief perform necessary initialization of AmgX.
         *
         * This function initializes AmgX for current instance. Based on the 
         * `count`, only the instance initialized first is in charge of 
         * initializing AmgX and the resource instance.
         *
         * \param cfgFile [in] the path to AmgX solver configuration file.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode initAmgX(const std::string &cfgFile);


        /**
         * \brief get IS for the row indices that processes in gpuWorld will held.
         *
         * \param A [in] PETSc matrix.
         * \param devIS [out] PETSc IS.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode getDevIS(const Mat &A, IS &devIS);


        /**
         * \brief get local sequential PETSc Mat of redistributed matrix.
         *
         * \param A [in] original PETSc Mat.
         * \param devIS [in] PETSc IS representing redistributed row indices.
         * \param localA [out] local sequential redistributed matrix.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode getLocalA(const Mat &A, const IS &devIS, Mat &localA);


        /**
         * \brief redistribute matrix.
         *
         * \param A [in] original PETSc Mat object.
         * \param devIS [in] PETSc IS representing redistributed rows.
         * \param newA [out] redistributed matrix.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode redistMat(const Mat &A, const IS &devIS, Mat &newA);


        /**
         * \brief get scatterLhs and scatterRhs.
         *
         * \param A1 [in] original PETSc Mat object.
         * \param A2 [in] redistributed PETSc Mat object.
         * \param devIS [in] PETSc IS representing redistributed row indices.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode getVecScatter(const Mat &A1, const Mat &A2, const IS &devIS);


        /**
         * \brief get data of compressed row layout of local sparse matrix.
         *
         * \param localA [in] sequential local redistributed PETSc Mat.
         * \param localN [out] number of local rows.
         * \param row [out] row vector in compressed row layout.
         * \param col [out] col vector in compressed row layout.
         * \param data [out] data vector in compressed row layout.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode getLocalMatRawData(const Mat &localA, 
                PetscInt &localN, std::vector<PetscInt> &row, 
                std::vector<PetscInt64> &col, std::vector<PetscScalar> &data);


        /**
         * \brief destroy the sequential local redistributed matrix.
         *
         * \param A [in] the original PETSc Mat.
         * \param localA [in, out] the local matrix. It will return null pointer.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode destroyLocalA(const Mat &A, Mat &localA);


        /**
         * \brief get a partition vector required by AmgX.
         *
         * \param devIS [in] PETSc IS representing redistributed row indices.
         * \param N [in] total number of rows in global matrix.
         * \param partVec [out] partition vector.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode getPartVec(const IS &devIS, 
                const PetscInt &N, std::vector<PetscInt> &partVec);


        /**
         * \brief function that actually solve the system.
         *
         * This function won't check if the process is involved in AmgX solver.
         * So if calling this function with processes not in the `gpuWorld`, 
         * something bad will happen. This function hence won't do data 
         * gathering/scattering, either.
         *
         * \param p [in, out] PETSc Vec for unknowns.
         * \param b [in] PETSc Vec for RHS.
         *
         * \return PetscErrorCode.
         */
        PetscErrorCode solve_real(Vec &p, Vec &b);
};

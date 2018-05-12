/**
 * \file setA.cpp
 * \brief definition of member functions regarding to setting A in AmgXSolver.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2016-01-08
 */


// STD
# include <cstring>

// AmgXSolver
# include "AmgXSolver.hpp"


// definition of AmgXSolver::setA
PetscErrorCode AmgXSolver::setA(const Mat &A)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    Mat                 localA;

    IS                  devIS;

    PetscInt            nGlobalRows,
                        nLocalRows;

    std::vector<PetscInt>       row;
    std::vector<PetscInt64>     col;
    std::vector<PetscScalar>    data;
    std::vector<PetscInt>       partVec;


    // get number of rows in global matrix
    ierr = MatGetSize(A, &nGlobalRows, nullptr); CHK;

    // get the row indices of redistributed matrix owned by processes in gpuWorld
    ierr = getDevIS(A, devIS); CHK;

    // get sequential local portion of redistributed matrix
    ierr = getLocalA(A, devIS, localA); CHK;

    // get compressed row layout of the local Mat
    ierr = getLocalMatRawData(localA, nLocalRows, row, col, data); CHK;

    // destroy local matrix
    ierr = destroyLocalA(A, localA); CHK;

    // get a partition vector required by AmgX
    ierr = getPartVec(devIS, nGlobalRows, partVec); CHK;


    // upload matrix A to AmgX
    if (gpuWorld != MPI_COMM_NULL)
    {
        ierr = MPI_Barrier(gpuWorld); CHK;

        AMGX_matrix_upload_all_global(
                AmgXA, nGlobalRows, nLocalRows, row[nLocalRows], 
                1, 1, row.data(), col.data(), data.data(), 
                nullptr, ring, ring, partVec.data());

        // bind the matrix A to the solver
        ierr = MPI_Barrier(gpuWorld); CHK;
        AMGX_solver_setup(solver, AmgXA);

        // connect (bind) vectors to the matrix
        AMGX_vector_bind(AmgXP, AmgXA);
        AMGX_vector_bind(AmgXRHS, AmgXA);
    }
    ierr = MPI_Barrier(globalCpuWorld); CHK;

    // destroy temporary PETSc objects
    ierr = ISDestroy(&devIS); CHK;

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::getDevIS
PetscErrorCode AmgXSolver::getDevIS(const Mat &A, IS &devIS)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    IS                  tempIS;

    // get index sets of A locally owned by each process
    // note that devIS is now a serial IS on each process
    ierr = MatGetOwnershipIS(A, &devIS, nullptr); CHK;

    // concatenate index sets that belong to the same devWorld
    // note that now devIS is a parallel IS of communicator devWorld
    ierr = ISOnComm(devIS, devWorld, PETSC_USE_POINTER, &tempIS); CHK;
    ierr = ISDestroy(&devIS); CHK;

    // all gather in order to have all indices belong to a devWorld on the 
    // leading rank of that devWorld. THIS IS NOT EFFICIENT!!
    // note that now devIS is again a serial IS on each process
    ierr = ISAllGather(tempIS, &devIS); CHK;
    ierr = ISDestroy(&tempIS); CHK;

    // empty devIS on ranks other than the leading ranks in each devWorld 
    if (myDevWorldRank != 0) 
        ierr = ISGeneralSetIndices(devIS, 0, nullptr, PETSC_COPY_VALUES); CHK;

    // devIS is not guaranteed to be sorted. We sort it here.
    ierr = ISSort(devIS); CHK;

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::getLocalA
PetscErrorCode AmgXSolver::getLocalA(const Mat &A, const IS &devIS, Mat &localA)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    MatType             type;

    // get the Mat type
    ierr = MatGetType(A, &type); CHK;

    // check whether the Mat type is supported
    if (std::strcmp(type, MATSEQAIJ) == 0) // sequential AIJ
    {
        // make localA point to the same memory space as A does
        localA = A;
    }
    else if (std::strcmp(type, MATMPIAIJ) == 0)
    {
        Mat                 tempA;

        // redistribute matrix and also get corresponding scatters.
        ierr = redistMat(A, devIS, tempA); CHK;

        // get local matrix from redistributed matrix
        ierr = MatMPIAIJGetLocalMat(tempA, MAT_INITIAL_MATRIX, &localA); CHK;

        // destroy redistributed matrix
        if (tempA == A)
        {
            tempA = nullptr;
        }
        else
        {
            ierr = MatDestroy(&tempA); CHK;
        }
    }
    else
    {
        SETERRQ1(globalCpuWorld, PETSC_ERR_ARG_WRONG,
                "Mat type %s is not supported!\n", type);
    }

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::redistMat
PetscErrorCode AmgXSolver::redistMat(const Mat &A, const IS &devIS, Mat &newA)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    if (gpuWorldSize == globalSize) // no redistributation required
    {
        newA = A;
    }
    else
    {
        IS      is;

        // re-set the communicator of devIS to globalCpuWorld
        ierr = ISOnComm(devIS, globalCpuWorld, PETSC_USE_POINTER, &is); CHK;

        // redistribute the matrix A to newA
        ierr = MatGetSubMatrix(A, is, is, MAT_INITIAL_MATRIX, &newA); CHK;

        // get VecScatters between original data layout and the new one
        ierr = getVecScatter(A, newA, is); CHK;

        // destroy the temporary IS
        ierr = ISDestroy(&is); CHK;
    }

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::getVecScatter
PetscErrorCode AmgXSolver::getVecScatter(
        const Mat &A1, const Mat &A2, const IS &devIS)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    Vec                 tempLhs;
    Vec                 tempRhs;

    ierr = MatCreateVecs(A1, &tempLhs, &tempRhs); CHK;
    ierr = MatCreateVecs(A2, &redistLhs, &redistRhs); CHK;

    ierr = VecScatterCreate(tempLhs, devIS, redistLhs, devIS, &scatterLhs); CHK;
    ierr = VecScatterCreate(tempRhs, devIS, redistRhs, devIS, &scatterRhs); CHK;
    
    ierr = VecDestroy(&tempRhs); CHK;
    ierr = VecDestroy(&tempLhs); CHK;

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::getLocalMatRawData
PetscErrorCode AmgXSolver::getLocalMatRawData(const Mat &localA,
        PetscInt &localN, std::vector<PetscInt> &row,
        std::vector<PetscInt64> &col, std::vector<PetscScalar> &data)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    const PetscInt      *rawCol, 
                        *rawRow;

    PetscScalar         *rawData;

    PetscInt            rawN;

    PetscBool           done;

    // get row and column indices in compressed row format
    ierr = MatGetRowIJ(localA, 0, PETSC_FALSE, PETSC_FALSE,
            &rawN, &rawRow, &rawCol, &done); CHK;

    // rawN will be returned by MatRestoreRowIJ, so we have to copy it
    localN = rawN;

    // check if the function worked
    if (! done)
        SETERRQ(globalCpuWorld, PETSC_ERR_SIG, "MatGetRowIJ did not work!");

    // get data
    ierr = MatSeqAIJGetArray(localA, &rawData); CHK;

    // copy values to STL vector. Note: there is an implicit conversion from 
    // PetscInt to PetscInt64 for the column vector
    col.assign(rawCol, rawCol+rawRow[localN]);
    row.assign(rawRow, rawRow+localN+1);
    data.assign(rawData, rawData+rawRow[localN]);


    // return ownership of memory space to PETSc
    ierr = MatRestoreRowIJ(localA, 0, PETSC_FALSE, PETSC_FALSE,
            &rawN, &rawRow, &rawCol, &done); CHK;

    // check if the function worked
    if (! done)
        SETERRQ(globalCpuWorld, PETSC_ERR_SIG, "MatRestoreRowIJ did not work!");

    // return ownership of memory space to PETSc
    ierr = MatSeqAIJRestoreArray(localA, &rawData); CHK;

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::destroyLocalA
PetscErrorCode AmgXSolver::destroyLocalA(const Mat &A, Mat &localA)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    MatType             type;

    // Get the Mat type
    ierr = MatGetType(A, &type); CHK;

    // when A is sequential, we can not destroy the memory space
    if (std::strcmp(type, MATSEQAIJ) == 0)
    {
        localA = nullptr;
    }
    // for parallel case, localA can be safely destroyed
    else if (std::strcmp(type, MATMPIAIJ) == 0)
    {
        ierr = MatDestroy(&localA); CHK;
    }

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::getPartVec
PetscErrorCode AmgXSolver::getPartVec(
        const IS &devIS, const PetscInt &N, std::vector<PetscInt> &partVec)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    VecScatter          scatter;
    Vec                 tempMPI,
                        tempSEQ;
    
    PetscInt            n;

    PetscScalar         *tempPartVec; 

    ierr = ISGetLocalSize(devIS, &n); CHK;

    if (gpuWorld != MPI_COMM_NULL)
    {
        ierr = VecCreateMPI(gpuWorld, n, N, &tempMPI); CHK;
    
        IS      is;
        ierr = ISOnComm(devIS, gpuWorld, PETSC_USE_POINTER, &is); CHK;
        ierr = VecISSet(tempMPI, is, (PetscScalar) myGpuWorldRank); CHK;
        ierr = ISDestroy(&is); CHK;

        ierr = VecScatterCreateToAll(tempMPI, &scatter, &tempSEQ); CHK;
        ierr = VecScatterBegin(scatter, 
                tempMPI, tempSEQ, INSERT_VALUES, SCATTER_FORWARD); CHK;
        ierr = VecScatterEnd(scatter, 
                tempMPI, tempSEQ, INSERT_VALUES, SCATTER_FORWARD); CHK;
        ierr = VecScatterDestroy(&scatter); CHK;
        ierr = VecDestroy(&tempMPI); CHK;

        ierr = VecGetArray(tempSEQ, &tempPartVec); CHK;

        partVec.assign(tempPartVec, tempPartVec+N);

        ierr = VecRestoreArray(tempSEQ, &tempPartVec); CHK;

        ierr = VecDestroy(&tempSEQ); CHK;
    }
    ierr = MPI_Barrier(globalCpuWorld); CHK;

    PetscFunctionReturn(0);
}

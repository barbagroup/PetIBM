/**
 * @file setA.inl
 * @brief Definition of some member functions of the class AmgXSolver
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2016-01-08
 */
# include "AmgXSolver.hpp"


/**
 * @brief A function convert PETSc Mat into AmgX matrix and bind it to solver
 *
 * This function will first extract the raw data from PETSc Mat and convert the
 * column index into 64bit integers. It also create a partition vector that is 
 * required by AmgX. Then, it upload the raw data to AmgX. Finally, it binds
 * the AmgX matrix to the AmgX solver.
 *
 * Be cautious! It lacks mechanism to check whether the PETSc Mat is AIJ format
 * and whether the PETSc Mat is using the same MPI communicator as the 
 * AmgXSolver instance.
 *
 * @param A A PETSc Mat. The coefficient matrix of a system of linear equations.
 * The matrix must be AIJ format and using the same MPI communicator as AmgX.
 *
 * @return Currently meaningless. May be error codes in the future.
 */
int AmgXSolver::setA(const Mat &A)
{
    PetscErrorCode      ierr;


    Mat                 localA;

    IS                  devIS;

    PetscInt            nGlobalRows,
                        nLocalRows;

    std::vector<PetscInt>       row;
    std::vector<Petsc64bitInt>  col;
    std::vector<PetscScalar>    data;
    std::vector<PetscInt>       partVec;


    // Get number of rows in global matrix
    ierr = MatGetSize(A, &nGlobalRows, nullptr);                             CHK;

    ierr = getDevIS(A, devIS);                                               CHK;
    ierr = getLocalA(A, devIS, localA);                                      CHK;
    ierr = getLocalMatRawData(localA, nLocalRows, row, col, data);           CHK;
    ierr = destroyLocalA(A, localA);                                         CHK;

    ierr = getPartVec(devIS, nGlobalRows, partVec);                          CHK;


    // upload matrix A to AmgX
    if (gpuWorld != MPI_COMM_NULL)
    {
        MPI_Barrier(gpuWorld);
        AMGX_matrix_upload_all_global(
                AmgXA, nGlobalRows, nLocalRows, row[nLocalRows], 
                1, 1, row.data(), col.data(), data.data(), 
                nullptr, ring, ring, partVec.data());

        // bind the matrix A to the solver
        MPI_Barrier(gpuWorld);
        AMGX_solver_setup(solver, AmgXA);

        // connect (bind) vectors to the matrix
        AMGX_vector_bind(AmgXP, AmgXA);
        AMGX_vector_bind(AmgXRHS, AmgXA);
    }
    MPI_Barrier(globalCpuWorld);

    return 0;
}


int AmgXSolver::getLocalA(const Mat &A, const IS &devIS, Mat &localA)
{
    PetscErrorCode      ierr;
    MatType             type;

    Mat                 tempA;

    // Get the Mat type
    ierr = MatGetType(A, &type);                                             CHK;

    // Check whether the Mat type is supported
    if (std::strcmp(type, MATSEQAIJ) == 0)
    {
        // make localA point to the same memory as A does
        localA = A;
    }
    else if (std::strcmp(type, MATMPIAIJ) == 0)
    {
        redistMat(A, devIS, tempA);

        ierr = MatMPIAIJGetLocalMat(
                tempA, MAT_INITIAL_MATRIX, &localA);                         CHK;

        if (tempA == A) 
        {
            tempA = nullptr;
        }
        else
        {
            ierr = MatDestroy(&tempA);                                       CHK;
        }
    }
    else
    {
        std::cerr << "Mat type " << type 
                  << " is not supported!" << std::endl;
        exit(0);
    }

    return 0;
}


int AmgXSolver::getDevIS(const Mat &A, IS &devIS)
{
    PetscErrorCode      ierr;

    ierr = MatGetOwnershipIS(A, &devIS, nullptr);                            CHK;

    ierr = ISOnComm(devIS, 
            devWorld, PETSC_USE_POINTER, &devIS);                            CHK;

    ierr = ISAllGather(devIS, &devIS);                                       CHK;

    if (myDevWorldRank != 0) 
        ierr = ISGeneralSetIndices(devIS, 0, nullptr, PETSC_COPY_VALUES);    CHK;

    ierr = ISSort(devIS);                                                    CHK;

    return 0;
}


int AmgXSolver::redistMat(const Mat &A, const IS &devIS, Mat &newA)
{
    PetscErrorCode      ierr;

    if (gpuWorldSize == globalSize)
    {
        newA = A;
    }
    else
    {
        ierr = MatGetSubMatrix(
            A, devIS, nullptr, MAT_INITIAL_MATRIX, &newA);                   CHK;

        ierr = getVecScatter(A, newA, devIS);                                CHK;
    }

    return 0;
}


int AmgXSolver::getPartVec(
        const IS &devIS, const PetscInt &N, std::vector<PetscInt> &partVec)
{
    PetscErrorCode      ierr;

    VecScatter          scatter;
    Vec                 tempMPI,
                        tempSEQ;
    
    PetscInt            n;

    PetscScalar         *tempPartVec; 

    ierr = ISGetLocalSize(devIS, &n);                                        CHK;

    if (gpuWorld != MPI_COMM_NULL)
    {
        ierr = VecCreateMPI(gpuWorld, n, N, &tempMPI);                       CHK;
    
        ierr = VecISSet(tempMPI, devIS, (PetscScalar) myGpuWorldRank);       CHK;

        ierr = VecScatterCreateToAll(tempMPI, &scatter, &tempSEQ);           CHK;
        ierr = VecScatterBegin(scatter, 
                tempMPI, tempSEQ, INSERT_VALUES, SCATTER_FORWARD);           CHK;
        ierr = VecScatterEnd(scatter, 
                tempMPI, tempSEQ, INSERT_VALUES, SCATTER_FORWARD);           CHK;

        ierr = VecScatterDestroy(&scatter);                                  CHK;
        ierr = VecDestroy(&tempMPI);                                         CHK;

        ierr = VecGetArray(tempSEQ, &tempPartVec);                           CHK;

        partVec.assign(tempPartVec, tempPartVec+N);

        ierr = VecRestoreArray(tempSEQ, &tempPartVec);                       CHK;

        ierr = VecDestroy(&tempSEQ);                                         CHK;
    }
    MPI_Barrier(globalCpuWorld);

    return 0;
}


int AmgXSolver::destroyLocalA(const Mat &A, Mat &localA)
{
    PetscErrorCode      ierr;
    MatType             type;

    // Get the Mat type
    ierr = MatGetType(A, &type);                                             CHK;

    // Check whether the Mat type is supported
    if (std::strcmp(type, MATSEQAIJ) == 0)
    {
        localA = nullptr;
    }
    else if (std::strcmp(type, MATMPIAIJ) == 0)
    {
        ierr = MatDestroy(&localA);                                          CHK;
    }

    return 0;
}


int AmgXSolver::getLocalMatRawData(Mat &localA, PetscInt &localN,
        std::vector<PetscInt> &row, std::vector<Petsc64bitInt> &col,
        std::vector<PetscScalar> &data)
{
    PetscErrorCode      ierr;
    PetscInt            tempN;
    const PetscInt      *rawCol, 
                        *rawRow;
    PetscScalar         *rawData;
    PetscBool           done;

    ierr = MatGetRowIJ(localA, 0, PETSC_FALSE, PETSC_FALSE,
            &tempN, &rawRow, &rawCol, &done);                                CHK;

    if (! done)
    {
        std::cerr << "MatGetRowIJ did not work!" << std::endl;
        exit(0);
    }

    ierr = MatSeqAIJGetArray(localA, &rawData);                              CHK;

    localN = tempN;

    col.assign(rawCol, rawCol+rawRow[localN]);
    row.assign(rawRow, rawRow+localN+1);
    data.assign(rawData, rawData+rawRow[localN]);


    ierr = MatRestoreRowIJ(localA, 0, PETSC_FALSE, PETSC_FALSE,
            &tempN, &rawRow, &rawCol, &done);                                CHK;

    // check whether MatRestoreRowIJ worked
    if (! done)
    {
        std::cerr << "MatRestoreRowIJ did not work!" << std::endl;
        exit(0);
    }

    ierr = MatSeqAIJRestoreArray(localA, &rawData);                          CHK;

    return 0;
}


int AmgXSolver::getVecScatter(const Mat &A1, const Mat &A2, const IS &devIS)
{
    PetscErrorCode      ierr;

    Vec                 tempV;

    ierr = MatCreateVecs(A1, nullptr, &tempV);                               CHK;
    ierr = MatCreateVecs(A2, nullptr, &redistRhs);                           CHK;
    ierr = MatCreateVecs(A2, nullptr, &redistLhs);                           CHK;

    ierr = VecScatterCreate(tempV, devIS, redistLhs, devIS, &redistScatter); CHK;

    ierr = VecDestroy(&tempV);                                               CHK;

    return 0;
}

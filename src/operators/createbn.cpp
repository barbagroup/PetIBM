/**
 * \file createbn.cpp
 * \brief Definition of functions for creating approximated inverse A.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// STL
# include <cmath>

// here goes PETSc headers
# include <petscmat.h>


namespace petibm
{
namespace operators
{

// implementation of createBnHead
PetscErrorCode createBnHead(const Mat &Op, const PetscReal &dt, 
        const PetscReal &coeff, const PetscInt &n, Mat &BnHead)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    if (n < 1) SETERRQ(
            PETSC_COMM_WORLD, 56, "The order of Bn can not be smaller than 1.");


    // a lazy and slow way to create BnHead. But we don't have to extract 
    // parallel distribution information.
    ierr = MatDuplicate(Op, MAT_DO_NOT_COPY_VALUES, &BnHead); CHKERRQ(ierr);
    ierr = MatSetOption(BnHead, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(BnHead, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

    // assume the BnHead is a diagonal matrix, which is only true for 1st order
    // so if the order is not 1, this will be inefficient
    ierr = MatSeqAIJSetPreallocation(BnHead, 1, nullptr); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(BnHead, 1, nullptr, 0, nullptr); CHKERRQ(ierr);

    // the first term; only diagonal values exist, and the value is dt
    // this one works only if the original diagonal is zero. Be careful.
    ierr = MatAssemblyBegin(BnHead, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(BnHead, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatShift(BnHead, dt); CHKERRQ(ierr);

    // if the order is not greater than 1, exit this function. No need to go 
    // through the following codes and save temporary memory usage.
    if (n == 1) PetscFunctionReturn(0);

    // the following code calculates each term, and add back to BnHead
    Mat     rightMat, tempResult;

    // calculate the matrix for each term if order > 1, and add back to BnHead
    for(PetscInt term=2; term<=n; ++term)
    {
        // the most right matrix of seccessive mat-mat-multiplication
        ierr = MatDuplicate(Op, MAT_COPY_VALUES, &rightMat); CHKERRQ(ierr);

        // successive mat-mat-multiplications
        for(PetscInt c=2; c<term; ++c)
        {
            // do tempResult = Op * rightMat once
            ierr = MatMatMult(Op, rightMat, MAT_INITIAL_MATRIX, 
                    PETSC_DEFAULT, &tempResult); CHKERRQ(ierr);

            // release the memory space where rightMat is pointing to
            ierr = MatDestroy(&rightMat); CHKERRQ(ierr);

            // make the pointer (Mat is actually a pointer) rightMat point 
            // to the same memory space pointed by tempResult
            rightMat = tempResult;
        }

        // the coefficient of this term
        PetscReal   a = std::pow(dt, term) * std::pow(coeff, term-1);

        // BnHead = BnHead + a * rightMat
        ierr = MatAXPY(BnHead, a, rightMat, DIFFERENT_NONZERO_PATTERN); 
        CHKERRQ(ierr);

        // assemble matrix
        ierr = MatAssemblyBegin(BnHead, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(BnHead, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

        // release the memory space pointed by rightMat (and also by tempResult)
        ierr = MatDestroy(&rightMat); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


// implementation of createBn
PetscErrorCode createBn(const Mat &Op, const Mat &R, const Mat &MHead,
        const PetscReal &dt, const PetscReal &coeff, const PetscInt &n, 
        Mat &Bn)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // get BnHead first
    ierr = createBnHead(Op, dt, coeff, n, Bn); CHKERRQ(ierr);

    Vec                 RDiag, MHeadDiagInv;

    // get the diagonal values of MHead^(-1)
    ierr = MatCreateVecs(MHead, nullptr, &MHeadDiagInv); CHKERRQ(ierr);
    ierr = MatGetDiagonal(MHead, MHeadDiagInv); CHKERRQ(ierr); 
    ierr = VecReciprocal(MHeadDiagInv); CHKERRQ(ierr);

    // get the diagonal values of R
    ierr = MatCreateVecs(R, nullptr, &RDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(R, RDiag); CHKERRQ(ierr); 

    // get Bn using Bn = R * BnHead * MHead^(-1)
    ierr = MatDiagonalScale(Bn, RDiag, MHeadDiagInv); CHKERRQ(ierr);

    // release memory space
    ierr = VecDestroy(&RDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&MHeadDiagInv); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// implementation of createBn
PetscErrorCode createBn(const Mat &Op, const Mat &M, const PetscReal &dt, 
        const PetscReal &coeff, const PetscInt &n, Mat &Bn)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    Vec                 MDiagInv;

    if (n < 1) SETERRQ(
            PETSC_COMM_WORLD, 56, "The order of Bn can not be smaller than 1.");


    // a lazy and slow way to create Bn. But we don't have to extract 
    // parallel distribution information.
    ierr = MatDuplicate(Op, MAT_DO_NOT_COPY_VALUES, &Bn); CHKERRQ(ierr);
    ierr = MatSetOption(Bn, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(Bn, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

    // assume the Bn is a diagonal matrix, which is only true for 1st order
    // so if the order is not 1, this will be inefficient
    ierr = MatSeqAIJSetPreallocation(Bn, 1, nullptr); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(Bn, 1, nullptr, 0, nullptr); CHKERRQ(ierr);

    // the first term; only diagonal values exist, and the value is dt
    // this one works only if the original diagonal is zero. Be careful.
    ierr = MatCreateVecs(M, nullptr, &MDiagInv); CHKERRQ(ierr);
    ierr = MatGetDiagonal(M, MDiagInv); CHKERRQ(ierr);
    ierr = VecReciprocal(MDiagInv); CHKERRQ(ierr);
    ierr = MatDiagonalSet(Bn, MDiagInv, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatScale(Bn, dt); CHKERRQ(ierr);

    // if the order is not greater than 1, exit this function. No need to go 
    // through the following codes and save temporary memory usage.
    if (n == 1) PetscFunctionReturn(0);


    // the following code calculates each term, and add back to Bn
    Mat     leftMat, rightMat, tempResult;

    // the left matrix we need for mat-mat-multiplication
    ierr = MatDuplicate(Op, MAT_COPY_VALUES, &leftMat); CHKERRQ(ierr);
    ierr = MatDiagonalScale(leftMat, MDiagInv, nullptr); CHKERRQ(ierr);

    for(PetscInt term=2; term<=n; ++term)
    {
        // the most right matrix of seccessive mat-mat-multiplication
        ierr = MatDuplicate(leftMat, MAT_COPY_VALUES, &rightMat); CHKERRQ(ierr);

        // successive mat-mat-multiplications
        for(PetscInt c=2; c<term; ++c)
        {
            // do tempResult = Op * rightMat once
            ierr = MatMatMult(leftMat, rightMat, MAT_INITIAL_MATRIX, 
                    PETSC_DEFAULT, &tempResult); CHKERRQ(ierr);

            // release the memory space where rightMat is pointing to
            ierr = MatDestroy(&rightMat); CHKERRQ(ierr);

            // make the pointer (Mat is actually a pointer) rightMat point 
            // to the same memory space pointed by tempResult
            rightMat = tempResult;
        }

        // do a final diagonal scaling on the right
        ierr = MatDiagonalScale(rightMat, nullptr, MDiagInv); CHKERRQ(ierr);

        // the coefficient of this term
        PetscReal   a = std::pow(dt, term) * std::pow(coeff, term-1);

        // Bn = Bn + a * rightMat
        ierr = MatAXPY(Bn, a, rightMat, DIFFERENT_NONZERO_PATTERN); 
        CHKERRQ(ierr);

        // assemble matrix
        ierr = MatAssemblyBegin(Bn, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Bn, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

        // release the memory space pointed by rightMat (and also by tempResult)
        ierr = MatDestroy(&rightMat); CHKERRQ(ierr);
    }

    // release the memory space
    ierr = MatDestroy(&leftMat); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

} // end of namespace operators
} // end of namespace petibm

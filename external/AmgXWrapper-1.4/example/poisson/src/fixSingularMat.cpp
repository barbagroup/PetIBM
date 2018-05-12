/**
 * @file fixSingularMat.cpp
 * @brief Definition of function fixSingularMat(...)
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @date 2016-02-13
 */


// header
# include "fixSingularMat.hpp"


// definition of fixSingularMat
PetscErrorCode fixSingularMat(Mat &A, Vec &RHS, const Vec &exact)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    PetscInt            row[1] = {0};

    PetscInt            n;

    PetscReal           scale;

    Vec                 diag;

    ierr = MatCreateVecs(A, nullptr, &diag); CHKERRQ(ierr);

    ierr = MatGetDiagonal(A, diag); CHKERRQ(ierr);

    ierr = VecGetSize(diag, &n); CHKERRQ(ierr);

    ierr = VecSum(diag, &scale); CHKERRQ(ierr);

    scale /= double(n);

    ierr = MatZeroRowsColumns(A, 1, row, scale, exact, RHS); CHKERRQ(ierr);

    ierr = VecDestroy(&diag); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

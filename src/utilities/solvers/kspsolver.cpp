/*! Implementation of the methods of the class `KSPSolver`.
 * \file kspsolver.cpp
 */

# include <cstring>

#include "kspsolver.h"


/*!
 * \brief Creates the KSP solver.
 */
PetscErrorCode KSPSolver::create(
        const Mat &A, const DM &daPack, const PetscBool &split)
{
  PetscErrorCode ierr;

  ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, 
                                options.c_str(), PETSC_FALSE); CHKERRQ(ierr);

  // create KSP for intermediate fluxes system
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(ksp, prefix.c_str()); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  // ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
  ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERRQ(ierr);


  DMType    type;
  ierr = DMGetType(daPack, &type); CHKERRQ(ierr);

  if (split)
  {
      if (std::strcmp(type, DMCOMPOSITE) == 0)
      {
          PetscInt  numDM;
          ierr = DMCompositeGetNumberDM(daPack, &numDM); CHKERRQ(ierr);

          ierr = setSplitPC(numDM, daPack); CHKERRQ(ierr);
      }
      else
      {
          ierr = PetscPrintf(PETSC_COMM_WORLD, "%s\n", type); CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD, "%s\n", DMCOMPOSITE); CHKERRQ(ierr);
          SETERRQ(PETSC_COMM_WORLD, 56, "Preconditioner is set to be split, "
                  "but the DM passed in is not a DMComposite!");
      }
  }

  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);
  ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  return 0;
} // create


PetscErrorCode KSPSolver::setSplitPC(const PetscInt &numDM, const DM &daPack)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    IS      *is;
    ierr = DMCompositeGetGlobalISs(daPack, &is); CHKERRQ(ierr);

    PC      pc;
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCFIELDSPLIT); CHKERRQ(ierr);
    ierr = PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE); CHKERRQ(ierr);

    for(int i=0; i<numDM; ++i)
    {
        ierr = PCFieldSplitSetIS(pc, nullptr, is[i]); CHKERRQ(ierr);
    }


    KSP     *childKsp;
    ierr = PCFieldSplitGetSubKSP(pc, nullptr, &childKsp); CHKERRQ(ierr);

    DM      *dms = new DM[numDM];
    ierr = DMCompositeGetEntriesArray(daPack, dms); CHKERRQ(ierr);

    for(int i=0; i<numDM; ++i)
    {
        ierr = KSPSetType(childKsp[i], KSPBCGS); CHKERRQ(ierr);
        ierr = KSPSetDM(childKsp[i], dms[i]); CHKERRQ(ierr);
        ierr = KSPSetDMActive(childKsp[i], PETSC_FALSE); CHKERRQ(ierr);
        ierr = KSPSetReusePreconditioner(childKsp[i], PETSC_TRUE); CHKERRQ(ierr);

        PC  pcField;
        ierr = KSPGetPC(childKsp[i], &pcField); CHKERRQ(ierr);
        ierr = PCSetType(pcField, PCBJACOBI); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


/*!
 * \brief Solves the system.
 */
PetscErrorCode KSPSolver::solve(Vec &x, Vec &b)
{
  PetscErrorCode ierr;

  ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
  if (reason < 0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, 
            "\nERROR: %s solver diverged due to reason: %d\n", 
            prefix.substr(0, prefix.size()-1).c_str(), reason); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);
    exit(1);
  }

  return 0;
} // solve


/*!
 * \brief Get the number of iterations performed.
 */
PetscErrorCode KSPSolver::getIters(PetscInt &iters)
{
  PetscErrorCode ierr;

  ierr = KSPGetIterationNumber(ksp, &iters); CHKERRQ(ierr);

  return 0;
} // getIters

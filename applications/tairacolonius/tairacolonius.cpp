/**
 * \file navierstokes.cpp
 * \brief Implementation of the class \c TairaColoniusSolver.
 */

// PETSc
# include <petscviewerhdf5.h>

// solver header
# include "tairacolonius.h"

using namespace petibm;


TairaColoniusSolver::TairaColoniusSolver(
        const petibm::type::Mesh &inMesh, const petibm::type::Boundary &inBC,
        const petibm::type::BodyPack &inBodies,const YAML::Node &node)
{
    initialize(inMesh, inBC, inBodies, node);
} // TairaColoniusSolver


PetscErrorCode TairaColoniusSolver::initialize(
        const petibm::type::Mesh &inMesh, const petibm::type::Boundary &inBC,
        const petibm::type::BodyPack &inBodies,const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    bodies = inBodies;
    ierr = PetscLogStageRegister("integrateForces", &stageIntegrateForces); CHKERRQ(ierr);

    ierr = NavierStokesSolver::initialize(inMesh, inBC, node); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // initialize


// create operators, i.e., PETSc Mats
PetscErrorCode TairaColoniusSolver::createOperators()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    
    Mat     R, MHat, BN, GH[2], DE[2]; // temporary operators
    Vec     RDiag, MHatDiag; // a temporary vectors
    IS      is[2];
    
    // create basic operators
    ierr = operators::createDivergence(
            mesh, bc, DE[0], DCorrection, PETSC_FALSE); CHKERRQ(ierr);
    ierr = operators::createGradient(mesh, GH[0], PETSC_FALSE); CHKERRQ(ierr);
    ierr = operators::createLaplacian(mesh, bc, L, LCorrection); CHKERRQ(ierr);
    ierr = operators::createConvection(mesh, bc, N); CHKERRQ(ierr);
    
    // create combined operator: A
    ierr = MatDuplicate(L, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
    ierr = MatScale(A, -diffCoeffs->implicitCoeff * nu); CHKERRQ(ierr);
    ierr = MatShift(A, 1.0/dt); CHKERRQ(ierr);

    // create diagonal matrix R and get a Vec holding the diagonal
    ierr = petibm::operators::createR(mesh, R); CHKERRQ(ierr);
    ierr = MatCreateVecs(R, nullptr, &RDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(R, RDiag); CHKERRQ(ierr);
    
    // create diagonal matrix MHat and get a Vec holding the diagonal
    ierr = petibm::operators::createMHead(mesh, MHat); CHKERRQ(ierr);
    ierr = MatCreateVecs(MHat, nullptr, &MHatDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(MHat, MHatDiag); CHKERRQ(ierr);
    
    // create a Delta operator and its transpose (equal to H)
    ierr = petibm::operators::createDelta(mesh, bc, bodies, DE[1]); CHKERRQ(ierr);
    ierr = MatTranspose(DE[1], MAT_INITIAL_MATRIX, &GH[1]); CHKERRQ(ierr);
    
    // create operator E
    ierr = MatDiagonalScale(DE[1], nullptr, RDiag); CHKERRQ(ierr);
    ierr = MatDiagonalScale(DE[1], nullptr, MHatDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&RDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&MHatDiag); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&MHat); CHKERRQ(ierr);
    
    // create -H (H operator in GH combination has a minus sign)
    ierr = MatScale(GH[1], -1.0); CHKERRQ(ierr);

    // get combined operator G; R is used temporarily
    ierr = MatCreateNest(mesh->comm, 1, nullptr, 2, nullptr, GH, &R); CHKERRQ(ierr);
    ierr = MatConvert(R, MATAIJ, MAT_INITIAL_MATRIX, &G); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&GH[0]); CHKERRQ(ierr);
    ierr = MatDestroy(&GH[1]); CHKERRQ(ierr);
    
    // get combined operator D; R is used temporarily; also get ISs
    ierr = MatCreateNest(mesh->comm, 2, nullptr, 1, nullptr, DE, &R); CHKERRQ(ierr);
    ierr = MatConvert(R, MATAIJ, MAT_INITIAL_MATRIX, &D); CHKERRQ(ierr);
    ierr = MatNestGetISs(R, is, nullptr); CHKERRQ(ierr);
    ierr = ISDuplicate(is[0], &isDE[0]); CHKERRQ(ierr);
    ierr = ISDuplicate(is[1], &isDE[1]); CHKERRQ(ierr);
    ierr = ISCopy(is[0], isDE[0]); CHKERRQ(ierr);
    ierr = ISCopy(is[1], isDE[1]); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&DE[0]); CHKERRQ(ierr);
    ierr = MatDestroy(&DE[1]); CHKERRQ(ierr);
    
    // create combined operator: BNG
    ierr = operators::createBnHead(
            L, dt, diffCoeffs->implicitCoeff * nu, 1, BN); CHKERRQ(ierr);
    ierr = MatMatMult(
            BN, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNG); CHKERRQ(ierr);
    
    // create combined operator: DBNG
    ierr = MatMatMult(
            D, BNG, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DBNG); CHKERRQ(ierr);
    
    // set null space to DBNG
    ierr = setNullSpace(); CHKERRQ(ierr);
    
    // destroy temporary operator
    ierr = MatDestroy(&BN); CHKERRQ(ierr);
    
    // register auto-destroy objects
    ierr = PetscObjectRegisterDestroy((PetscObject) DBNG); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) BNG); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) A); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) N); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) LCorrection); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) L); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) G); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) DCorrection); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) D); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleOperators


// create vectors
PetscErrorCode TairaColoniusSolver::createVectors()
{
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr;
    Vec temp;
    const PetscReal *data;
    
    // P; combination of pGlobal and forces
    ierr = MatCreateVecs(G, &P, nullptr); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) P); CHKERRQ(ierr);
    
    // swap pGlobal and P, so we can reuse functions from Navier-Stokes solver
    temp = solution->pGlobal;
    solution->pGlobal = P;
    P = temp;
    temp = PETSC_NULL;
    
    // reset the underlying data of P to the pressure portion in pGlobal
    ierr = VecGetSubVector(solution->pGlobal, isDE[0], &temp); CHKERRQ(ierr);
    ierr = VecGetArrayRead(temp, &data); CHKERRQ(ierr);
    ierr = VecReplaceArray(P, data); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(temp, &data); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(solution->pGlobal, isDE[0], &temp); CHKERRQ(ierr);
    
    // other vectors follow the same routine in N-S solver
    ierr = NavierStokesSolver::createVectors(); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


PetscErrorCode TairaColoniusSolver::setNullSpace()
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    std::string type;
    
    ierr = pSolver->getType(type); CHKERRQ(ierr);
    
    if (type == "PETSc KSP")
    {
        Vec n, phiPortion;
        MatNullSpace nsp;
        ierr = MatCreateVecs(DBNG, &n, nullptr); CHKERRQ(ierr);
        ierr = VecSet(n, 0.0); CHKERRQ(ierr);
        ierr = VecGetSubVector(n, isDE[0], &phiPortion); CHKERRQ(ierr);
        ierr = VecSet(phiPortion, 1.0 / std::sqrt(mesh->pN)); CHKERRQ(ierr);
        ierr = VecRestoreSubVector(n, isDE[0], &phiPortion); CHKERRQ(ierr);
        ierr = MatNullSpaceCreate(
                mesh->comm, PETSC_FALSE, 1, &n, &nsp); CHKERRQ(ierr);
        ierr = MatSetNullSpace(DBNG, nsp); CHKERRQ(ierr);
        ierr = MatSetNearNullSpace(DBNG, nsp); CHKERRQ(ierr);
        ierr = VecDestroy(&n); CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(&nsp); CHKERRQ(ierr);
        isRefP = PETSC_FALSE;
    }
    else if (type == "NVIDIA AmgX")
    {
        PetscInt row[1] = {0};
        ierr = MatZeroRowsColumns(
                DBNG, 1, row, 1.0, nullptr, nullptr); CHKERRQ(ierr);
        isRefP = PETSC_TRUE;
    }
    else
    {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "Could not recognize the type of linear solver: %s\n",
                type.c_str());
    }
    
    PetscFunctionReturn(0);
}


// prepare tight-hand-side vector for Poisson system
PetscErrorCode TairaColoniusSolver::assembleRHSPoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSPoisson); CHKERRQ(ierr);
    
    // get Du* (which is equal to (D_{interior} + D_{boundary})u*
    ierr = MatMult(D, solution->UGlobal, rhs2); CHKERRQ(ierr);
    
    // note, bc2 is a subset of rhs2
    Vec     bc2;
    ierr = VecGetSubVector(rhs2, isDE[0], &bc2); CHKERRQ(ierr);
    ierr = MatMultAdd(DCorrection, solution->UGlobal, bc2, bc2); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(rhs2, isDE[0], &bc2); CHKERRQ(ierr);
    
    // apply reference pressure 
    if (isRefP)
    {
        ierr = VecSetValue(rhs2, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(rhs2); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(rhs2); CHKERRQ(ierr);
    }
    
    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleRHSPoisson


// output solutions to the user provided file
PetscErrorCode TairaColoniusSolver::write(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    
    Vec temp = PETSC_NULL;
    
    // let solution->pGlobal point to P, so that we can use solution->write
    temp = solution->pGlobal;
    solution->pGlobal = P;

    ierr = NavierStokesSolver::write(filePath); CHKERRQ(ierr);
    
    // restore pointers
    solution->pGlobal = temp;
    temp = PETSC_NULL;

    PetscFunctionReturn(0);
} // write


// write averaged forces
PetscErrorCode TairaColoniusSolver::writeIntegratedForces(
            const PetscReal &t, const std::string &filePath)
{
    PetscErrorCode ierr;
    PetscViewer         viewer;
    petibm::type::RealVec2D     fAvg;
    static PetscFileMode mode = FILE_MODE_WRITE;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageIntegrateForces); CHKERRQ(ierr);
    
    // get sub section f and calculate averaged forces
    Vec f;
    ierr = VecGetSubVector(solution->pGlobal, isDE[1], &f); CHKERRQ(ierr);
    ierr = bodies->calculateAvgForces(f, fAvg); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(solution->pGlobal, isDE[1], &f); CHKERRQ(ierr);
    
    ierr = PetscLogStagePop(); CHKERRQ(ierr);
    
    // create ASCII viewer
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, mode); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    
    // write current time
    ierr = PetscViewerASCIIPrintf(viewer, "%10.8e", t); CHKERRQ(ierr);
    
    // write forces body by body
    for(unsigned int i=0; i<bodies->nBodies; ++i)
    {
        for(unsigned int d=0; d<mesh->dim; ++d)
        {
            ierr = PetscViewerASCIIPrintf(
                    viewer, "\t%10.8e", fAvg[i][d]); CHKERRQ(ierr);
        }
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);
    
    // destroy
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    // next we'll just append data
    mode = FILE_MODE_APPEND;

    PetscFunctionReturn(0);
} // writeIntegratedForces


// manual finalization
PetscErrorCode TairaColoniusSolver::finalize()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    bodies.~shared_ptr();
    ierr = ISDestroy(&isDE[0]); CHKERRQ(ierr);
    ierr = ISDestroy(&isDE[1]); CHKERRQ(ierr);
    ierr = NavierStokesSolver::finalize(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // finalize

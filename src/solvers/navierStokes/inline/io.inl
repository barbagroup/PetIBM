/***************************************************************************//**
 * \file io.inl
 * \author Olivier Mesnard (mesnardo@gwu), Anush Krishnan (anush@bu.edu)
 * \brief Implementation of I/O methods of the class NavierStokesSolver.
 */

#include "types.h"

#include <petscviewerhdf5.h>


/**
 * \brief Read fluxes from file.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readFluxes()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n[time-step %d] Reading fluxes from file... ",
                     timeStep); CHKERRQ(ierr);

  // get solution directory: 7 characters long, time-step preprend by leading zeros
  std::stringstream ss;
  ss << parameters->directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str();

  // get access to the individual vectors of the composite vector
  // depending on the dimension of the problem
  Vec qxGlobal, qyGlobal, qzGlobal;
  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }

  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string fileExtension;
  if (parameters->fileFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->fileFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }
  
  // read fluxes in x-direction
  std::string filePath = solutionDirectory + "/qx." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) qxGlobal, "qx"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(qxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  // read fluxes in y-direction
  filePath = solutionDirectory + "/qy." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) qyGlobal, "qy"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(qyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  if (dim == 3)
  {
    // read fluxes in z-direction
    filePath = solutionDirectory + "/qz." + fileExtension;
    ierr = PetscObjectSetName((PetscObject) qzGlobal, "qz"); CHKERRQ(ierr);
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    ierr = VecLoad(qzGlobal, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }

  if (dim == 2)
  {
    ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);
  
  return 0;
} // readFluxes


/**
 * \brief Reads the pressure field from saved numerical solution file.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readLambda()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n[time-step %d] Reading pressure from file... ",
                     timeStep); CHKERRQ(ierr);

  // get solution directory: 7 characters long, time-step preprend by leading zeros
  std::stringstream ss;
  ss << parameters->directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str();

  // get access to the pressure vector from the composite vector
  Vec phi;
  ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string fileExtension;
  if (parameters->fileFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->fileFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  // read pressure field
  std::string filePath = solutionDirectory + "/phi." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) phi, "phi"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // readLambda


/**
 * \brief Writes the numerical solution into respective files.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeData()
{
  PetscErrorCode ierr;

  ierr = writeIterationCounts(); CHKERRQ(ierr);

  if (timeStep%parameters->nsave == 0)
  {
    ierr = writeFluxes(); CHKERRQ(ierr);
    ierr = writeLambda(); CHKERRQ(ierr);
  }

  return 0;
} // writeData


/**
 * \brief Writes grid stations of the different field variables in HDF5 files.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeGrids()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n[time-step %d] Writing grids into files... ",
                     timeStep); CHKERRQ(ierr);

  // create subfolder `grids`
  std::string subfolder = parameters->directory + "/grids";
  mkdir(subfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  std::string filePath = subfolder + "/grid_qx.h5";
  ierr = mesh->write(filePath, STAGGERED_MODE_X); CHKERRQ(ierr);
  filePath = subfolder + "/grid_qy.h5";
  ierr = mesh->write(filePath, STAGGERED_MODE_Y); CHKERRQ(ierr);
  if (dim == 3)
  {
    filePath = subfolder + "/grid_qz.h5";
    ierr = mesh->write(filePath, STAGGERED_MODE_Z); CHKERRQ(ierr);
  }
  filePath = subfolder + "/grid_phi.h5";
  ierr = mesh->write(filePath, CELL_CENTERED); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // writeGrids


/**
 * \brief Writes fluxes into files located in the time-step directory.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeFluxes()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n[time-step %d] Writing fluxes into file... ",
                     timeStep); CHKERRQ(ierr);

  // create the solution directory
  std::stringstream ss;
  ss << parameters->directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str();
  mkdir(solutionDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // get access to the individual vectors of the composite vector
  // depending on the dimension of the problem
  Vec qxGlobal, qyGlobal, qzGlobal;
  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }

  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string fileExtension;
  if (parameters->fileFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->fileFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  // write fluxes in x-direction
  std::string filePath = solutionDirectory + "/qx." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) qxGlobal, "qx"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(qxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  // write fluxes in y-direction
  filePath = solutionDirectory + "/qy." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) qyGlobal, "qy"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(qyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  if (dim == 3)
  {
    // write fluxes in z-direction
    filePath = solutionDirectory + "/qz." + fileExtension;
    ierr = PetscObjectSetName((PetscObject) qzGlobal, "qz"); CHKERRQ(ierr);
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    ierr = VecView(qzGlobal, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }

  if (dim == 2)
  {
    ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // writeFluxes


/**
 * \brief Writes the pressure field into file located in solution directory.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeLambda()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n[time-step %d] Writing pressure into file... ",
                     timeStep); CHKERRQ(ierr);

  // create the solution directory
  std::stringstream ss;
  ss << parameters->directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str();
  mkdir(solutionDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // get access to the pressure vector from the composite vector
  Vec phi;
  ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string fileExtension;
  if (parameters->fileFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->fileFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  // write pressure field
  std::string filePath = solutionDirectory + "/phi." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) phi, "phi"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // writeLambda


/**
 * \brief Writes the iteration count for each KSP solver into a file.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeIterationCounts()
{
  PetscErrorCode ierr;

  PetscInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  if (rank == 0)
  {
    PetscInt countVelocitySolver,
             countPoissonSolver;
    std::string filePath = parameters->directory + "/iterationCounts.txt";
    if (timeStep == 1)
    {
      iterationCountsFile.open(filePath.c_str());
    }
    else
    {
      iterationCountsFile.open(filePath.c_str(), std::ios::out | std::ios::app);
    }
    ierr = KSPGetIterationNumber(ksp1, &countVelocitySolver); CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp2, &countPoissonSolver); CHKERRQ(ierr);
    iterationCountsFile << timeStep << '\t' \
                        << countVelocitySolver << '\t' \
                        << countPoissonSolver << std::endl;
    iterationCountsFile.close();
  }

  return 0;
} // writeIterationCounts


/**
 * \brief Code-development helper: outputs vectors to files.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::helperOutputVectors()
{
  PetscErrorCode ierr;
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d] Code-development: saving vectors to files... ", timeStep); CHKERRQ(ierr);

  // create the output directory
  std::string outputDirectory = parameters->directory + "/outputs";
  mkdir(outputDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  PetscViewer viewer;
  // bc1
  std::string filePath = outputDirectory + "/bc1.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(bc1, viewer); CHKERRQ(ierr);
  // H
  filePath = outputDirectory + "/H.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(H, viewer); CHKERRQ(ierr);
  // rn
  filePath = outputDirectory + "/rn.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(rn, viewer); CHKERRQ(ierr);
  // rhs1
  filePath = outputDirectory + "/rhs1.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(rhs1, viewer); CHKERRQ(ierr);
  // q
  filePath = outputDirectory + "/q.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(q, viewer); CHKERRQ(ierr);
  // qx, qy, qz
  Vec qxGlobal, qyGlobal, qzGlobal;
  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);  
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }
  filePath = outputDirectory + "/qx.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(qxGlobal, viewer); CHKERRQ(ierr);
  filePath = outputDirectory + "/qy.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(qyGlobal, viewer); CHKERRQ(ierr);
  if (dim == 3)
  {
    filePath = outputDirectory + "/qz.output";
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
    ierr = VecView(qzGlobal, viewer); CHKERRQ(ierr);
  }
  // r2
  filePath = outputDirectory + "/r2.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(r2, viewer); CHKERRQ(ierr);
  // rhs2
  filePath = outputDirectory + "/rhs2.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(rhs2, viewer); CHKERRQ(ierr);
  // lambda
  filePath = outputDirectory + "/lambda.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = VecView(lambda, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // helperOutputVectors


/**
 * \brief Code-development helper: outputs matrices to files.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::helperOutputMatrices()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d] Code-development: saving matrices to files... ", timeStep); CHKERRQ(ierr);

  // create the output directory
  std::string outputDirectory = parameters->directory + "/outputs";
  mkdir(outputDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  PetscViewer viewer;
  // A
  std::string filePath = outputDirectory + "/A.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = MatView(A, viewer); CHKERRQ(ierr);
  // QT
  filePath = outputDirectory + "/QT.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = MatView(QT, viewer); CHKERRQ(ierr);
  // BNQ
  filePath = outputDirectory + "/BNQ.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = MatView(BNQ, viewer); CHKERRQ(ierr);
  // QTBNQ
  filePath = outputDirectory + "/QTBNQ.output";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = MatView(QTBNQ, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // helperOutputMatrices

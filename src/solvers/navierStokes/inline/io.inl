/***************************************************************************//**
 * \file io.inl
 * \author Olivier Mesnard (mesnardo@gwu), Anush Krishnan (anush@bu.edu)
 * \brief Implementation of I/O methods of the class NavierStokesSolver.
 */

#include "types.h"

#include <petscviewerhdf5.h>


/**
 * \brief Reads the numerical solution from files.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readData()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n[time-step %d] Reading numerical solution from files... ",
                     timeStep); CHKERRQ(ierr);

  // get solution directory
  std::stringstream ss;
  ss << parameters->directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str();

  if (parameters->outputFlux)
  {
    ierr = readFluxes(solutionDirectory); CHKERRQ(ierr);
  }
  else if (parameters->outputVelocity)
  {
    ierr = readVelocities(solutionDirectory); CHKERRQ(ierr);
  }
  ierr = readLambda(solutionDirectory); CHKERRQ(ierr);
  if (parameters->startStep > 0)
  {
    ierr = readConvectiveTerms(solutionDirectory); CHKERRQ(ierr);
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // readData


/**
 * \brief Read flux fields from files.
 *
 * \param directory Directory where to read the flux fields.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readFluxes(std::string directory)
{
  PetscErrorCode ierr;
  Vec qxGlobal, qyGlobal, qzGlobal;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // get type of viewer depending on output format prescribed
  if (parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }
  
  // read fluxes in x-direction
  filePath = directory + "/qx." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) qxGlobal, "qx"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(qxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  // read fluxes in y-direction
  filePath = directory + "/qy." + fileExtension;
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
    filePath = directory + "/qz." + fileExtension;
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
  
  PetscFunctionReturn(0);
} // readFluxes


/**
 * \brief Read velocity fields from files.
 *
 * \param directory Directory where to read the velocity fields.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readVelocities(std::string directory)
{
  PetscErrorCode ierr;
  Vec u, uxGlobal, uyGlobal, uzGlobal;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // define the type of viewer and the file extension
  if (parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  ierr = VecDuplicate(q, &u); CHKERRQ(ierr);
  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, u, &uxGlobal, &uyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, u, &uxGlobal, &uyGlobal, &uzGlobal); CHKERRQ(ierr);
  }
  
  // read x-component of the velocity
  filePath = directory + "/ux." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) uxGlobal, "ux"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(uxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  // read y-component of the velocity
  filePath = directory + "/uy." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) uyGlobal, "uy"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(uyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  if (dim == 3)
  {
    // read z-component of the velocity
    filePath = directory + "/uz." + fileExtension;
    ierr = PetscObjectSetName((PetscObject) uzGlobal, "uz"); CHKERRQ(ierr);
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    ierr = VecLoad(uzGlobal, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }

  if (dim == 2)
  {
    ierr = DMCompositeRestoreAccess(qPack, u, &uxGlobal, &uyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeRestoreAccess(qPack, u, &uxGlobal, &uyGlobal, &uzGlobal); CHKERRQ(ierr);
  }

  // convert velocity into flux
  ierr = VecPointwiseDivide(q, u, RInv); CHKERRQ(ierr);

  ierr = VecDestroy(&u); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
} // readVelocities


/**
 * \brief Reads convective terms.
 *
 * \param directory Directory where to read the convective terms
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readConvectiveTerms(std::string directory)
{
  PetscErrorCode ierr;
  Vec HxGlobal, HyGlobal, HzGlobal;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // get type of viewer depending on output format prescribed
  if (parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, H, &HxGlobal, &HyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, H, &HxGlobal, &HyGlobal, &HzGlobal); CHKERRQ(ierr);
  }
  
  // read convective terms in x-direction
  filePath = directory + "/Hx." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) HxGlobal, "Hx"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(HxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  // read convective terms in y-direction
  filePath = directory + "/Hy." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) HyGlobal, "Hy"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(HyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  if (dim == 3)
  {
    // read convective terms in z-direction
    filePath = directory + "/Hz." + fileExtension;
    ierr = PetscObjectSetName((PetscObject) HzGlobal, "Hz"); CHKERRQ(ierr);
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    ierr = VecLoad(HzGlobal, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }

  if (dim == 2)
  {
    ierr = DMCompositeRestoreAccess(qPack, H, &HxGlobal, &HyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeRestoreAccess(qPack, H, &HxGlobal, &HyGlobal, &HzGlobal); CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
} // readConvectiveTerms


/**
 * \brief Reads the pressure field from file.
 *
 * \param directory Directory where to read the pressure field.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readLambda(std::string directory)
{
  PetscErrorCode ierr;
  Vec phi;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // define the type of viewer and the file extension
  if (parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  // read pressure field
  filePath = directory + "/phi." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) phi, "phi"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // readLambda


/**
 * \brief Writes the numerical solution into respective files.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeData()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = writeIterationCounts(); CHKERRQ(ierr);

  if (timeStep%parameters->nsave == 0 || timeStep%parameters->nrestart == 0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "\n[time-step %d] Writing numerical solution into files... ",
                       timeStep); CHKERRQ(ierr);

    // create solution directory
    std::stringstream ss;
    ss << parameters->directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
    std::string solutionDirectory = ss.str();
    mkdir(solutionDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (parameters->outputFlux)
    {
      ierr = writeFluxes(solutionDirectory); CHKERRQ(ierr);
    }
    if (parameters->outputVelocity)
    {
      ierr = writeVelocities(solutionDirectory); CHKERRQ(ierr);
    }
    ierr = writeLambda(solutionDirectory); CHKERRQ(ierr);
    if (timeStep%parameters->nrestart == 0)
    {
      ierr = writeConvectiveTerms(solutionDirectory); CHKERRQ(ierr);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // writeData


#ifdef PETSC_HAVE_HDF5
/**
 * \brief Writes grid stations of the different field variables in HDF5 files.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeGrids()
{
  PetscErrorCode ierr;
  std::string gridsDirectory, filePath;

  PetscFunctionBeginUser;

  // create subfolder `grids`
  gridsDirectory = parameters->directory + "/grids";
  mkdir(gridsDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  filePath = gridsDirectory + "/staggered-x.h5";
  ierr = mesh->write(filePath, STAGGERED_MODE_X); CHKERRQ(ierr);
  filePath = gridsDirectory + "/staggered-y.h5";
  ierr = mesh->write(filePath, STAGGERED_MODE_Y); CHKERRQ(ierr);
  if (dim == 3)
  {
    filePath = gridsDirectory + "/staggered-z.h5";
    ierr = mesh->write(filePath, STAGGERED_MODE_Z); CHKERRQ(ierr);
  }
  filePath = gridsDirectory + "/cell-centered.h5";
  ierr = mesh->write(filePath, CELL_CENTERED); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // writeGrids
#endif


/**
 * \brief Writes flux fields into files.
 *
 * \param directory Directory where to write the flux fields.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeFluxes(std::string directory)
{
  PetscErrorCode ierr;
  Vec qxGlobal, qyGlobal, qzGlobal;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // define the type of viewer and the file extension
  if (parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }

  // write fluxes in x-direction
  filePath = directory + "/qx." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) qxGlobal, "qx"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(qxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  // write fluxes in y-direction
  filePath = directory + "/qy." + fileExtension;
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
    filePath = directory + "/qz." + fileExtension;
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

  PetscFunctionReturn(0);
} // writeFluxes


/**
 * \brief Writes velocity fields into files.
 *
 * \param directory Directory where to write the velocity fields.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeVelocities(std::string directory)
{
  PetscErrorCode ierr;
  Vec u, uxGlobal, uyGlobal, uzGlobal;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // define the type of viewer and the file extension
  if (parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  // convert flux into velocity
  ierr = VecDuplicate(q, &u);
  ierr = VecPointwiseMult(u, q, RInv); CHKERRQ(ierr);
  
  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, u, &uxGlobal, &uyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, u, &uxGlobal, &uyGlobal, &uzGlobal); CHKERRQ(ierr);
  }

  // write x-component of the velocity
  filePath = directory + "/ux." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) uxGlobal, "ux"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(uxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  // write y-component of the velocity
  filePath = directory + "/uy." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) uyGlobal, "uy"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(uyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  if (dim == 3)
  {
    // write z-component of the velocity
    filePath = directory + "/uz." + fileExtension;
    ierr = PetscObjectSetName((PetscObject) uzGlobal, "uz"); CHKERRQ(ierr);
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    ierr = VecView(uzGlobal, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }

  if (dim == 2)
  {
    ierr = DMCompositeRestoreAccess(qPack, u, &uxGlobal, &uyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeRestoreAccess(qPack, u, &uxGlobal, &uyGlobal, &uzGlobal); CHKERRQ(ierr);
  }

  ierr = VecDestroy(&u); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // writeVelocities


/**
 * \brief Writes the convective terms into files.
 *
 * \param directory Directory where to write the convective terms.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeConvectiveTerms(std::string directory)
{
  PetscErrorCode ierr;
  Vec HxGlobal, HyGlobal, HzGlobal;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // define the type of viewer and the file extension
  if (parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  if (dim == 2)
  {
    ierr = DMCompositeGetAccess(qPack, H, &HxGlobal, &HyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeGetAccess(qPack, H, &HxGlobal, &HyGlobal, &HzGlobal); CHKERRQ(ierr);
  }

  // convective terms in x-direction
  filePath = directory + "/Hx." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) HxGlobal, "Hx"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(HxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  // convective terms in y-direction
  filePath = directory + "/Hy." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) HyGlobal, "Hy"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(HyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  if (dim == 3)
  {
    // convective terms in z-direction
    filePath = directory + "/Hz." + fileExtension;
    ierr = PetscObjectSetName((PetscObject) HzGlobal, "Hz"); CHKERRQ(ierr);
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    ierr = VecView(HzGlobal, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }

  if (dim == 2)
  {
    ierr = DMCompositeRestoreAccess(qPack, H, &HxGlobal, &HyGlobal); CHKERRQ(ierr);
  }
  else if (dim == 3)
  {
    ierr = DMCompositeRestoreAccess(qPack, H, &HxGlobal, &HyGlobal, &HzGlobal); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // writeConvectiveTerms


/**
 * \brief Writes the pressure field into file.
 *
 * \param directory Directory where to write the pressure field.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeLambda(std::string directory)
{
  PetscErrorCode ierr;
  Vec phi;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // define the type of viewer and the file extension
  if (parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  // write pressure field
  filePath = directory + "/phi." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) phi, "phi"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // writeLambda


/**
 * \brief Writes the iteration count for each KSP solver into a file.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeIterationCounts()
{
  PetscErrorCode ierr;

  PetscMPIInt rank;
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
    ierr = velocity->getIters(countVelocitySolver); CHKERRQ(ierr);
    ierr = poisson->getIters(countPoissonSolver); CHKERRQ(ierr);
    iterationCountsFile << timeStep << '\t' \
                        << countVelocitySolver << '\t' \
                        << countPoissonSolver << std::endl;
    iterationCountsFile.close();
  }

  return 0;
} // writeIterationCounts

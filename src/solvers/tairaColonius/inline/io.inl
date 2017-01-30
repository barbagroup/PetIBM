/***************************************************************************//**
 * \file io.inl
 * \author Olivier Mesnard (mesnardo@gwu), Anush Krishnan (anush@bu.edu)
 * \brief Implementation of I/O methods of the class `TairaColoniusSolver`.
 */


/**
 * \brief Reads the pressure field and body forces from files.
 *
 * \param directory Directory where to read the solutions.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::readLambda(std::string directory)
{
  PetscErrorCode ierr;
  Vec phi, fTilde;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  // define the type of viewer and the file extension
  if (NavierStokesSolver<dim>::parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (NavierStokesSolver<dim>::parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }
  
  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  // read pressure field
  filePath = directory + "/phi." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) phi, "phi"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecLoad(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  // read body forces if restarting the simulation
  if (NavierStokesSolver<dim>::timeStep > 0)
  {
    filePath = directory + "/fTilde." + fileExtension;
    ierr = PetscObjectSetName((PetscObject) fTilde, "fTilde"); CHKERRQ(ierr);
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    ierr = VecLoad(fTilde, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // readLambda


/**
 * \brief Writes the numerical solution into files.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeData()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = calculateForcesTC(); CHKERRQ(ierr);
  ierr = writeForces(); CHKERRQ(ierr);
  ierr = NavierStokesSolver<dim>::writeData(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // writeData


/**
 * \brief Writes the pressure field and the body forces into files.
 *
 * \param directory Directory where to write the solutions.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeLambda(std::string directory)
{
  PetscErrorCode ierr;
  Vec phi, fTilde;
  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string filePath, fileExtension;

  PetscFunctionBeginUser;

  if (NavierStokesSolver<dim>::parameters->outputFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (NavierStokesSolver<dim>::parameters->outputFormat == "binary")
  {
    viewerType = PETSCVIEWERBINARY;
    fileExtension = "dat";
  }

  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  // write pressure field
  filePath = directory + "/phi." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) phi, "phi"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  // write body forces
  filePath = directory + "/fTilde." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) fTilde, "fTilde"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(fTilde, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // writeLambda


/*!
 * \brief Writes force in each direction acting on the body.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeForces()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  if (rank == 0)
  {
    std::string filePath = NavierStokesSolver<dim>::parameters->directory + "/forces.txt";
    if (NavierStokesSolver<dim>::timeStep == 1)
    {
      forcesFile.open(filePath.c_str());
    }
    else
    {
      forcesFile.open(filePath.c_str(), std::ios::out | std::ios::app);
    }
    forcesFile << NavierStokesSolver<dim>::timeStep*NavierStokesSolver<dim>::parameters->dt;
    for (PetscInt bIdx=0; bIdx<numBodies; bIdx++)
    {
      for (PetscInt d=0; d<dim; d++)
      {
        forcesFile << '\t' << bodies[bIdx].forces[d];
      }
    }
    forcesFile << std::endl;
    forcesFile.close();
  }

  PetscFunctionReturn(0);
} // writeForces

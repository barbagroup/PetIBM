/***************************************************************************//**
 * \file io.inl
 * \author Olivier Mesnard (mesnardo@gwu), Anush Krishnan (anush@bu.edu)
 * \brief Implementation of I/O methods of the class `TairaColoniusSolver`.
 */


/**
 * \brief Reads the pressure field and body forces from saved numerical solution file.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::readLambda()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n[time-step %d] Reading pressure and body forces from file... ",
                     NavierStokesSolver<dim>::timeStep); CHKERRQ(ierr);

  // get solution directory: 7 characters long, time-step preprend by leading zeros
  std::stringstream ss;
  ss << NavierStokesSolver<dim>::parameters->directory << "/" << std::setfill('0') << std::setw(7) << NavierStokesSolver<dim>::timeStep;
  std::string solutionDirectory = ss.str();

  // get access to the pressure vector from the composite vector
  Vec phi, fTilde;
  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string fileExtension;
  if (NavierStokesSolver<dim>::parameters->fileFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (NavierStokesSolver<dim>::parameters->fileFormat == "binary")
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

  // read body forces iff restarting the simulation
  if (NavierStokesSolver<dim>::timeStep > 0)
  {
    filePath = solutionDirectory + "/fTilde." + fileExtension;
    ierr = PetscObjectSetName((PetscObject) fTilde, "fTilde"); CHKERRQ(ierr);
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    ierr = VecLoad(fTilde, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // readLambda


/**
 * \brief Writes the numerical solution into respective files.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeData()
{
  PetscErrorCode ierr;

  ierr = calculateForce(); CHKERRQ(ierr);
  ierr = writeForces(); CHKERRQ(ierr);
  ierr = NavierStokesSolver<dim>::writeData(); CHKERRQ(ierr);

  return 0;
} // writeData


/**
 * \brief Writes the pressure field and the body forces into files 
 *        located in solution directory.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeLambda()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n[time-step %d] Writing pressure and body forces into file... ",
                     NavierStokesSolver<dim>::timeStep); CHKERRQ(ierr);

  // create the solution directory
  std::stringstream ss;
  ss << NavierStokesSolver<dim>::parameters->directory << "/" << std::setfill('0') << std::setw(7) << NavierStokesSolver<dim>::timeStep;
  std::string solutionDirectory = ss.str();
  mkdir(solutionDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // get access to the pressure vector from the composite vector
  Vec phi, fTilde;
  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  PetscViewer viewer;
  PetscViewerType viewerType;
  std::string fileExtension;
  if (NavierStokesSolver<dim>::parameters->fileFormat == "hdf5")
  {
    viewerType = PETSCVIEWERHDF5;
    fileExtension = "h5";
  }
  else if (NavierStokesSolver<dim>::parameters->fileFormat == "binary")
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

  // write body forces
  filePath = solutionDirectory + "/fTilde." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) fTilde, "fTilde"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(fTilde, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // writeLambda


/**
 * \brief Writes force in each direction acting on the body.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeForces()
{
  PetscErrorCode ierr;
  PetscInt rank;
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
    for (PetscInt i=0; i<dim; i++)
    {
      forcesFile << '\t' << force[i];
    }
    forcesFile << std::endl;
    forcesFile.close();
  }

  return 0;
} // writeForces
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
  PetscViewer viewer;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d] Reading pressure and body forces from file... ", NavierStokesSolver<dim>::timeStep); CHKERRQ(ierr);

  // get solution directory: 7 characters long, time-step preprend by leading zeros
  std::stringstream ss;
  ss << NavierStokesSolver<dim>::parameters->directory << "/" << std::setfill('0') << std::setw(7) << NavierStokesSolver<dim>::timeStep;
  std::string solutionDirectory = ss.str();

  // get access to the pressure vector from the composite vector
  Vec phi, fTilde;
  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  // read pressure field
  std::string filePath = solutionDirectory + "/phi.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
  ierr = VecLoad(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);


  // read body forces iff restarting the simulation
  if (NavierStokesSolver<dim>::timeStep > 0)
  {
    filePath = solutionDirectory + "/fTilde.dat";
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
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
  PetscViewer viewer;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d] Writing pressure and body forces into file... ", NavierStokesSolver<dim>::timeStep); CHKERRQ(ierr);

  // create the solution directory
  std::stringstream ss;
  ss << NavierStokesSolver<dim>::parameters->directory << "/" << std::setfill('0') << std::setw(7) << NavierStokesSolver<dim>::timeStep;
  std::string solutionDirectory = ss.str();
  mkdir(solutionDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // get access to the pressure vector from the composite vector
  Vec phi, fTilde;
  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

  // write pressure field
  std::string filePath = solutionDirectory + "/phi.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  ierr = VecView(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  // write body forces
  filePath = solutionDirectory + "/fTilde.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
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
/***************************************************************************//**
 * \file io.inl
 * \author Olivier Mesnard (mesnardo@gwu), Anush Krishnan (anush@bu.edu)
 * \brief Implementation of I/O methods of the class NavierStokesSolver.
 */


/**
 * \brief Prints information about the simulation.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::printInfo()
{
  PetscErrorCode ierr;
  ierr = mesh->printInfo(); CHKERRQ(ierr);
  ierr = flow->printInfo(); CHKERRQ(ierr);
  ierr = parameters->printInfo(); CHKERRQ(ierr);
  return 0;
} // printInfo


/**
 * \brief Read fluxes from file.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readFluxes()
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d] Reading fluxes... ", timeStep); CHKERRQ(ierr);

  // get solution directory: 7 characters long, time-step preprend by leading zeros
  std::stringstream ss;
  ss << directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
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
  
  // read fluxes in x-direction
  std::string filePath = solutionDirectory + "/qx.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
  ierr = VecLoad(qxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  // read fluxes in y-direction
  filePath = solutionDirectory + "/qy.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
  ierr = VecLoad(qyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  if (dim == 3)
  {
    // read fluxes in z-direction
    filePath = solutionDirectory + "/qz.dat";
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
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
  PetscViewer viewer;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d] Reading pressure... ", timeStep); CHKERRQ(ierr);

  // get solution directory: 7 characters long, time-step preprend by leading zeros
  std::stringstream ss;
  ss << directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str();

  // get access to the pressure vector from the composite vector
  Vec phi;
  ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  // read pressure field
  std::string filePath = solutionDirectory + "/phi.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
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
 * \brief Writes fluxes into files located in the time-step directory.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeFluxes()
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d] Writing fluxes... ", timeStep); CHKERRQ(ierr);

  // create the solution directory
  std::stringstream ss;
  ss << directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
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

  // write fluxes in x-direction
  std::string filePath = solutionDirectory + "/qx.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  ierr = VecView(qxGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  // write fluxes in y-direction
  filePath = solutionDirectory + "/qy.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  ierr = VecView(qyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  if (dim == 3)
  {
    // write fluxes in z-direction
    filePath = solutionDirectory + "/qz.dat";
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
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
  PetscViewer viewer;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d] Writing pressure... ", timeStep); CHKERRQ(ierr);

  // create the solution directory
  std::stringstream ss;
  ss << directory << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str();
  mkdir(solutionDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // get access to the pressure vector from the composite vector
  Vec phi;
  ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  // write pressure field
  std::string filePath = solutionDirectory + "/phi.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  ierr = VecView(phi, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // writeLambda


/**
 * \brief Writes the grid into file.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeGrid()
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nWriting grid... "); CHKERRQ(ierr);

  PetscInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  if (rank == 0)
  {
    std::ofstream streamFile(directory + "/grid.txt");
    if (dim == 2)
    {
      streamFile << mesh->nx << '\t' << mesh->ny << '\n';
    }
    else if (dim == 3)
    {
      streamFile << mesh->nx << '\t' << mesh->ny << '\t' << mesh->nz << '\n';
    }
    for (std::vector<PetscReal>::const_iterator i=mesh->x.begin(); i!=mesh->x.end(); ++i)
      streamFile << *i << '\n';
    for (std::vector<PetscReal>::const_iterator i=mesh->y.begin(); i!=mesh->y.end(); ++i)
      streamFile << *i << '\n';
    if (dim == 3)
    {  
      for (std::vector<PetscReal>::const_iterator i=mesh->z.begin(); i!=mesh->z.end(); ++i)
        streamFile << *i << '\n';
    }
    streamFile.close();
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n", timeStep); CHKERRQ(ierr);

  return 0;
} // writeGrid


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
    std::string filePath = directory + "/iterationCounts.txt";
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
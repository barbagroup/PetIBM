/*! Implementation of I/O methods of the class `LiEtAlSolver`.
 * \file io.inl
 */


/*!
 * \brief Writes the numerical solution into respective files.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::writeData()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = calculateForces(); CHKERRQ(ierr);
  // ierr = calculateForces2(); CHKERRQ(ierr);
  ierr = writeForces(); CHKERRQ(ierr);

  ierr = writeIterationCounts(); CHKERRQ(ierr);

  if (NavierStokesSolver<dim>::timeStep%NavierStokesSolver<dim>::parameters->nsave == 0 ||
      NavierStokesSolver<dim>::timeStep%NavierStokesSolver<dim>::parameters->nrestart == 0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "\n[time-step %d] Writing numerical solution into files... ",
                       NavierStokesSolver<dim>::timeStep); CHKERRQ(ierr);

    // create solution directory
    std::stringstream ss;
    ss << NavierStokesSolver<dim>::parameters->directory << "/" << std::setfill('0') << std::setw(7) << NavierStokesSolver<dim>::timeStep;
    std::string solutionDirectory = ss.str();
    mkdir(solutionDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (NavierStokesSolver<dim>::parameters->outputFlux)
    {
      ierr = NavierStokesSolver<dim>::writeFluxes(solutionDirectory); CHKERRQ(ierr);
    }
    if (NavierStokesSolver<dim>::parameters->outputVelocity)
    {
      ierr = NavierStokesSolver<dim>::writeVelocities(solutionDirectory); CHKERRQ(ierr);
    }
    ierr = NavierStokesSolver<dim>::writeLambda(solutionDirectory); CHKERRQ(ierr);
    if (NavierStokesSolver<dim>::timeStep%NavierStokesSolver<dim>::parameters->nrestart == 0)
    {
      ierr = NavierStokesSolver<dim>::writeConvectiveTerms(solutionDirectory); CHKERRQ(ierr);
    }
    ierr = writeLagrangianForces(solutionDirectory); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // writeData



/*!
 * \brief Writes the Lagrangian forces into a file.
 *
 * \param directory Directory where to write the solutions.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::writeLagrangianForces(std::string directory)
{
  PetscErrorCode ierr;
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

  // write body forces
  filePath = directory + "/fTilde." + fileExtension;
  ierr = PetscObjectSetName((PetscObject) fTilde, "fTilde"); CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr); 
  ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
  ierr = VecView(fTilde, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // writeLambda


/*!
 * \brief Writes the iteration count for each KSP solver into a file.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::writeIterationCounts()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  if (rank == 0)
  {
    PetscInt countVelocitySolver,
             countPoissonSolver,
             countForceSolver;
    std::string filePath = NavierStokesSolver<dim>::parameters->directory + "/iterationCounts.txt";
    if (NavierStokesSolver<dim>::timeStep == 1)
    {
      NavierStokesSolver<dim>::iterationCountsFile.open(filePath.c_str());
    }
    else
    {
      NavierStokesSolver<dim>::iterationCountsFile.open(filePath.c_str(), std::ios::out | std::ios::app);
    }
    ierr = NavierStokesSolver<dim>::velocity->getIters(countVelocitySolver); CHKERRQ(ierr);
    ierr = NavierStokesSolver<dim>::poisson->getIters(countPoissonSolver); CHKERRQ(ierr);
    ierr = forces->getIters(countForceSolver); CHKERRQ(ierr);
    NavierStokesSolver<dim>::iterationCountsFile << NavierStokesSolver<dim>::timeStep << '\t' \
                                                 << countVelocitySolver << '\t' \
                                                 << countPoissonSolver << '\t' \
                                                 << countForceSolver << std::endl;
    NavierStokesSolver<dim>::iterationCountsFile.close();
  }

  PetscFunctionReturn(0);
} // writeIterationCounts


/*!
 * \brief Writes force in each direction acting on the body.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::writeForces()
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
    // for (PetscInt d=0; d<dim; d++)
    // {
    //   forcesFile << '\t' << bodyForces[d];
    // }
    forcesFile << std::endl;
    forcesFile.close();
  }

  PetscFunctionReturn(0);
} // writeForces

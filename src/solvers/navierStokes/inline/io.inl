/***************************************************************************//**
 * \file io.inl
 * \author Olivier Mesnard (mesnardo@gwu), Anush Krishnan (anush@bu.edu)
 * \brief Implementation of I/O methods of the class NavierStokesSolver.
 */


/**
 * \brief Prints the simulation parameters
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::printSimulationInfo()
{
  PetscErrorCode ierr;
  PetscInt       rank, maxits;
  PC             pc;
  PCType         pcType;
  KSPType        kspType;
  PetscReal      rtol, abstol;
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Flow\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "solver: "); CHKERRQ(ierr);
  switch (simParams->solverType)
  {
    case NAVIER_STOKES :
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Navier-Stokes\n"); CHKERRQ(ierr); 
      break;
    case TAIRA_COLONIUS: 
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Taira & Colonius (2007)\n"); CHKERRQ(ierr);
      break;
    default:
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Unrecognized solver!\n"); CHKERRQ(ierr);
      break;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "viscosity: %g\n", flowDesc->nu); CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Mesh\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  if (dim == 3)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "size: %d x %d x %d\n", mesh->nx, mesh->ny, mesh->nz); CHKERRQ(ierr);
  }
  else
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "size: %d x %d\n", mesh->nx, mesh->ny); CHKERRQ(ierr);
  }
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Time-stepping\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "convection: "); CHKERRQ(ierr);
  switch (simParams->convectionScheme)
  {
    case EULER_EXPLICIT   :
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Explicit Euler\n"); CHKERRQ(ierr);
      break;
    case ADAMS_BASHFORTH_2:
      ierr = PetscPrintf(PETSC_COMM_WORLD, "2nd-order Adams-Bashforth\n"); CHKERRQ(ierr);
      break;
    default:
      break;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "diffusion : "); CHKERRQ(ierr);
  switch (simParams->diffusionScheme)
  {
    case EULER_EXPLICIT:
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Explicit Euler\n"); CHKERRQ(ierr);
      break;
    case EULER_IMPLICIT:
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Implicit Euler\n"); CHKERRQ(ierr);
      break;
    case CRANK_NICOLSON:
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Crank-Nicolson\n"); CHKERRQ(ierr);
      break;
    default:
      break;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "time-increment      : %g\n", simParams->dt); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "starting time-step  : %d\n", simParams->startStep); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "number of time-steps: %d\n", simParams->nt); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "saving-interval     : %d\n", simParams->nsave); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for intermediate velocity\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = KSPGetType(ksp1, &kspType); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "solver: %s\n", kspType); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp1, &pc); CHKERRQ(ierr);
  ierr = PCGetType(pc, &pcType); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "preconditioner: %s\n", pcType); CHKERRQ(ierr);
  ierr = KSPGetTolerances(ksp1, &rtol, &abstol, NULL, &maxits); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %g\n", rtol); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %g\n", abstol); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "maximum iterations: %d\n", maxits); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for pressure-force\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = KSPGetType(ksp2, &kspType); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "solver: %s\n", kspType); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp2, &pc); CHKERRQ(ierr);
  ierr = PCGetType(pc, &pcType); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "preconditioner: %s\n", pcType); CHKERRQ(ierr);
  ierr = KSPGetTolerances(ksp2, &rtol, &abstol, NULL, &maxits); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %g\n", rtol); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %g\n", abstol); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "maximum iterations: %d\n", maxits); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);

  return 0;
} // printSimulationInfo


/**
 * \brief Read fluxes from file.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readFluxes()
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nReading fluxes at time step %d...\n", timeStep); CHKERRQ(ierr);

  // get solution directory: 7 characters long, time-step preprend by leading zeros
  std::stringstream ss;
  ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str()

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

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nReading pressure at time step %d...\n", timeStep); CHKERRQ(ierr);

  // get solution directory: 7 characters long, time-step preprend by leading zeros
  std::stringstream ss;
  ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
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

  if (timeStep%simParams->nsave == 0)
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

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nWriting fluxes at time step %d...\n", timeStep); CHKERRQ(ierr);

  // create the solution directory
  std::stringstream ss;
  ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str()
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
  std::string filePath = solutionDirectory + "/qy.dat";
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filePath.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  ierr = VecView(qyGlobal, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  if (dim == 3)
  {
    // write fluxes in z-direction
    std::string filePath = solutionDirectory + "/qz.dat";
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

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nWriting pressure at time step %d...\n", timeStep); CHKERRQ(ierr);

  // create the solution directory
  std::stringstream ss;
  ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
  std::string solutionDirectory = ss.str()
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

  return 0;
} // writeLambda


/**
 * \brief Writes the grid into file.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeGrid()
{
  PetscErrorCode ierr;
  PetscInt rank;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nWriting grid...\n"); CHKERRQ(ierr);

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  if (rank == 0)
  {
    std::ofstream stream = caseFolder + "/grid.txt";
    if (dim == 2)
    {
      stream << mesh->nx << '\t' << mesh->ny << '\n';
    }
    else if (dim == 3)
    {
      stream << mesh->nx << '\t' << mesh->ny << '\t' << mesh->nz << '\n';
    }
    for (std::vector<PetscReal>::const_iterator i=mesh->x.begin(); i!=mesh->x.end(); ++i)
      stream << *i << '\n';
    for (std::vector<PetscReal>::const_iterator i=mesh->y.begin(); i!=mesh->y.end(); ++i)
      stream << *i << '\n';
    if (dim == 3)
      for (std::vector<PetscReal>::const_iterator i=mesh->z.begin(); i!=mesh->z.end(); ++i)
        stream << *i << '\n';
    stream.close();
  }

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
    std::string filePath = caseFolder + "iterationCounts.txt";
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
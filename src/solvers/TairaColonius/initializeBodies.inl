/***************************************************************************//**
 * \file initializeBodies.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `initializeBodies` of \c TairaColoniusSolver.
 */


/**
 * \brief Initializes the immersed boundaries by parsing the input file `bodies.yaml`.
 *
 * Two-dimensional simulations.
 */
template <>
PetscErrorCode TairaColoniusSolver<2>::initializeBodies()
{
  PetscErrorCode ierr;
  PetscInt       rank;
  PetscInt       totalPoints;
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  if (rank == 0)
  {
    std::string bodiesFile = caseFolder + "/bodies.yaml";
    YAML::Node nodes = YAML::LoadFile(bodiesFile);
    const YAML::Node &node = nodes[0];

    std::string type = node["type"].as<std::string>();

    if (type == "circle")
    {
      PetscReal xc = node["circleOptions"][0].as<PetscReal>();
      PetscReal yc = node["circleOptions"][1].as<PetscReal>();
      PetscReal R = node["circleOptions"][2].as<PetscReal>();
      PetscInt numPoints = node["circleOptions"][3].as<PetscInt>();

      x.reserve(numPoints);
      y.reserve(numPoints);
      for (PetscInt i=0; i<numPoints; i++)
      {
        x.push_back(xc + R*cos(2.0*PETSC_PI*i/numPoints));
        y.push_back(yc + R*sin(2.0*PETSC_PI*i/numPoints));
      }
    }
    else if (type == "points")
    {
      PetscInt numPoints;
      PetscReal xCoord, yCoord;
      std::string pointsFile = caseFolder + "/" + node["pointsFile"].as<std::string>();
      std::cout << "Initiliazing body: reading coordinates from: " << pointsFile << std::endl;
      std::ifstream infile(pointsFile.c_str());
      infile >> numPoints;
      x.reserve(numPoints);
      y.reserve(numPoints);
      for (PetscInt i=0; i<numPoints; i++)
      {
        infile >> xCoord >> yCoord;
        x.push_back(xCoord);
        y.push_back(yCoord);
      }
      infile.close();
    }
    else if (type == "lineSegment")
    {
      PetscReal startX = node["segmentOptions"][0].as<PetscReal>();
      PetscReal startY = node["segmentOptions"][1].as<PetscReal>();
      PetscReal endX = node["segmentOptions"][2].as<PetscReal>();
      PetscReal endY = node["segmentOptions"][3].as<PetscReal>();
      PetscInt numPoints = node["segmentOptions"][4].as<PetscInt>();
      // initialize line segment
      x.reserve(numPoints);
      y.reserve(numPoints);
      PetscReal xi;
      for (PetscInt i=0; i<numPoints; i++)
      {
        xi = (i+0.5)/numPoints;
        x.push_back((1-xi)*startX + xi*endX);
        y.push_back((1-xi)*startY + xi*endY);
      }
    }
    else
    {
      std::cout << "\nERROR: Unknown type of body" << std::endl;
      exit(0);
    }
    totalPoints = x.size();
  }

  // broadcast total number of body points to all processes
  ierr = MPI_Bcast(&totalPoints, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

  ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

  // allocate memory
  x.resize(totalPoints);
  y.resize(totalPoints);

  // broadcast vectors to all processes
  ierr = MPI_Bcast(&x.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
  ierr = MPI_Bcast(&y.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

  return 0;
}

/**
 * \brief Initializes the immersed boundaries by parsing the input file `bodies.yaml`.
 *
 * Three-dimensional simulations.
 */
template <>
PetscErrorCode TairaColoniusSolver<3>::initializeBodies()
{
  PetscErrorCode ierr;
  PetscInt       rank;
  PetscInt       totalPoints;
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  if (rank == 0)
  {
    std::string bodiesFile = caseFolder + "/bodies.yaml";
    YAML::Node nodes = YAML::LoadFile(bodiesFile);
    const YAML::Node &node = nodes[0];

    std::string type = node["type"].as<std::string>();

    if (type == "quad")
    {
      PetscInt nXi = node["quadOptions"][0].as<PetscInt>(),
               nEta = node["quadOptions"][0].as<PetscInt>();
      PetscReal corners[4][3];

      x.reserve(nXi*nEta);
      y.reserve(nXi*nEta);
      z.reserve(nXi*nEta);

      for (size_t d=0; d<3; d++)
      {
        corners[0][d] = node["bottomLeft"][d].as<PetscReal>();
        corners[1][d] = node["bottomRight"][d].as<PetscReal>();
        corners[2][d] = node["topRight"][d].as<PetscReal>();
        corners[3][d] = node["topLeft"][d].as<PetscReal>();
      }
      
      PetscReal xi, eta;
      for (PetscInt j=0; j<nEta; j++)
      {
        eta = (j+0.5)/nEta;
        for (PetscInt i=0; i<nXi; i++)
        {
          xi = (i+0.5)/nXi;
          x.push_back((1-xi)*(1-eta)*corners[0][0] 
                      + xi*(1-eta)*corners[1][0] 
                      + xi*eta*corners[2][0] 
                      + (1-xi)*eta*corners[3][0]);
          y.push_back((1-xi)*(1-eta)*corners[0][1] 
                      + xi*(1-eta)*corners[1][1] 
                      + xi*eta*corners[2][1] 
                      + (1-xi)*eta*corners[3][1]);
          z.push_back((1-xi)*(1-eta)*corners[0][2] 
                      + xi*(1-eta)*corners[1][2] 
                      + xi*eta*corners[2][2] 
                      + (1-xi)*eta*corners[3][2]);
        }
      }
    }
    else if (type == "points")
    {
      std::string pointsFile = caseFolder + "/" + node["pointsFile"].as<std::string>();
      std::cout << "Initiliazing body: reading coordinates from: " << pointsFile << std::endl;
      std::ifstream infile(pointsFile.c_str());
      PetscInt numPoints;
      PetscReal xCoord, yCoord, zCoord;
      infile >> numPoints;
      x.reserve(numPoints);
      y.reserve(numPoints);
      z.reserve(numPoints);
      for (PetscInt i=0; i<numPoints; i++)
      {
        infile >> xCoord >> yCoord >> zCoord;
        x.push_back(xCoord);
        y.push_back(yCoord);
        z.push_back(zCoord);
      }
      infile.close();
    }
    else
    {
      std::cout << "\nERROR: Unknown type of body" << std::endl;
      exit(0);
    }
    totalPoints = x.size();
  }

  // broadcast total number of body points to all processes
  ierr = MPI_Bcast(&totalPoints, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

  ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

  // allocate memory
  x.resize(totalPoints);
  y.resize(totalPoints);
  z.resize(totalPoints);

  // broadcast vectors to all processes
  ierr = MPI_Bcast(&x.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
  ierr = MPI_Bcast(&y.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
  ierr = MPI_Bcast(&z.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

  return 0;
}
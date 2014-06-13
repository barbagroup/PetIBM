#include "yaml-cpp/yaml.h"
#include <fstream>

template <>
PetscErrorCode TairaColoniusSolver<2>::initialiseBodies()
{
	PetscErrorCode ierr;
	PetscInt       rank;
	PetscInt       totalPoints;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if(rank==0)
	{
		std::string bodiesFile = caseFolder + "/bodies.yaml";
		std::ifstream f(bodiesFile.c_str());
		YAML::Parser  parser(f);
		YAML::Node    doc;
		parser.GetNextDocument(doc);
		const YAML::Node &node = doc[0];

		std::string type;
		node["type"] >> type;

		if (type == "circle")
		{
			PetscReal cx, cy, R;
			PetscInt  numPoints;
			node["circleOptions"][0] >> cx;
			node["circleOptions"][1] >> cy;
			node["circleOptions"][2] >> R;
			node["circleOptions"][3] >> numPoints;

			x.reserve(numPoints);
			y.reserve(numPoints);
			for(PetscInt i=0; i<numPoints; i++)
			{
				x.push_back(cx + R*cos(i*2*PETSC_PI/numPoints));
				y.push_back(cy + R*sin(i*2*PETSC_PI/numPoints));
			}
		}
		else if (type == "points")
		{
			std::string fname;
			PetscInt    numPoints;
			PetscReal   xCoord, yCoord;
			node["pointsFile"] >> fname;
			fname = caseFolder + "/" + fname;
			std::cout << "Reading body data from file: " << fname << std::endl;

			std::ifstream file(fname.c_str());
			file >> numPoints;
			x.reserve(numPoints);
			y.reserve(numPoints);
			for(PetscInt i=0; i<numPoints; i++)
			{
				file >> xCoord >> yCoord;
				x.push_back(xCoord);
				y.push_back(yCoord);
			}
			file.close();
		}
		else if (type == "lineSegment")
		{
			PetscReal startX, startY, endX, endY;
			PetscInt  numPoints;
			node["segmentOptions"][0] >> startX;
			node["segmentOptions"][1] >> endX;
			node["segmentOptions"][2] >> startY;
			node["segmentOptions"][3] >> endY;
			node["segmentOptions"][4] >> numPoints;
			// initialise line segment
			x.reserve(numPoints);
			y.reserve(numPoints);
			for(PetscInt i=0; i<numPoints; i++)
			{
				PetscReal xi = (i+0.5)/numPoints;
				x.push_back((1-xi)*startX + xi*endX);
				y.push_back((1-xi)*startY + xi*endY);
			}
		}
		else
		{
			std::cout << "[E]: unknown Body type\n";
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

	ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of body points: %d\n", x.size()); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode TairaColoniusSolver<3>::initialiseBodies()
{
	PetscErrorCode ierr;
	PetscInt       rank;
	PetscInt       totalPoints;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if(rank==0)
	{
		std::string bodiesFile = caseFolder + "/bodies.yaml";
		std::ifstream f(bodiesFile.c_str());
		YAML::Parser  parser(f);
		YAML::Node    doc;
		parser.GetNextDocument(doc);
		const YAML::Node &node = doc[0];

		std::string type;
		node["type"] >> type;

		if (type == "quad")
		{
			PetscInt  nXi, nEta;
			PetscReal corners[4][3];
			node["quadOptions"][0] >> nXi;
			node["quadOptions"][1] >> nEta;

			x.reserve(nXi*nEta);
			y.reserve(nXi*nEta);
			z.reserve(nXi*nEta);

			for(size_t d=0; d<3; d++)
			{
				node["bottomLeft"][d] >> corners[0][d];
				node["bottomRight"][d] >> corners[1][d];
				node["topRight"][d] >> corners[2][d];
				node["topLeft"][d] >> corners[3][d];
			}
			
			PetscReal xi, eta;

			for(PetscInt j=0; j<nEta; j++)
			{
				eta = (j+0.5)/nEta;
				for(PetscInt i=0; i<nXi; i++)
				{
					xi = (i+0.5)/nXi;
					x.push_back((1-xi)*(1-eta)*corners[0][0] + xi*(1-eta)*corners[1][0] + xi*eta*corners[2][0] + (1-xi)*eta*corners[3][0]);
					y.push_back((1-xi)*(1-eta)*corners[0][1] + xi*(1-eta)*corners[1][1] + xi*eta*corners[2][1] + (1-xi)*eta*corners[3][1]);
					z.push_back((1-xi)*(1-eta)*corners[0][2] + xi*(1-eta)*corners[1][2] + xi*eta*corners[2][2] + (1-xi)*eta*corners[3][2]);
				}
			}
		}
		else if (type == "points")
		{
			std::string fname;
			PetscInt    numPoints;
			PetscReal   xCoord, yCoord, zCoord;
			node["pointsFile"] >> fname;
			fname = caseFolder + "/" + fname;
			std::cout << "Reading body data from file: " << fname << std::endl;
			// initialise points
			std::ifstream file(fname.c_str());
			file >> numPoints;
			x.reserve(numPoints);
			y.reserve(numPoints);
			z.reserve(numPoints);
			for(PetscInt i=0; i<numPoints; i++)
			{
				file >> xCoord >> yCoord >> zCoord;
				x.push_back(xCoord);
				y.push_back(yCoord);
				z.push_back(zCoord);
			}
			file.close();
		}
		else
		{
			std::cout << "[E]: unknown Body type\n";
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

	ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of body points: %d\n", x.size()); CHKERRQ(ierr);

	return 0;
}

#include "yaml-cpp/yaml.h"
#include <fstream>
template <>
void TairaColoniusSolver<2>::initialiseBodies()
{
	PetscInt rank, numProcs;
	PetscInt totalPoints;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

	numBoundaryPointsOnProcess.resize(numProcs);
	boundaryPointIndices.resize(numProcs);
	startGlobalIndices.resize(numProcs);

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
			// initialise circle
			x.reserve(numPoints);
			y.reserve(numPoints);
			for(PetscInt i=0; i<numPoints; i++)
			{
				x.push_back(cx + R*cos(i*2*PETSC_PI/numPoints));
				y.push_back(cy + R*sin(i*2*PETSC_PI/numPoints));
			}
		}
	/*	else if (type == "points")
		{
			string fname;
			node["pointsFile"] >> fname;
			fname = "bodyFiles/" + fname;
			std::cout << fname << std::endl;
			// initialise points
			std::ifstream file(fname.c_str());
			file >> Body.numPoints;
			Body.X.resize(Body.numPoints);
			Body.Y.resize(Body.numPoints);
			for(int i=0; i<Body.numPoints; i++)
			{
				file >> Body.X[i] >> Body.Y[i];
			}
			file.close();
		}
		else if (type == "lineSegment")
		{
			real startX, startY, endX, endY;
			int numPoints;
			node["segmentOptions"][0] >> startX;
			node["segmentOptions"][1] >> endX;
			node["segmentOptions"][2] >> startY;
			node["segmentOptions"][3] >> endY;
			node["segmentOptions"][4] >> numPoints;
			Body.numPoints = numPoints;
			// initialise line segment
		}*/
		else
			std::cout << "[E]: unknown Body type\n";
		totalPoints = x.size();
	}

	// broadcast total number of body points to all processes
	MPI_Bcast(&totalPoints, 1, MPIU_INT, 0, PETSC_COMM_WORLD);

	MPI_Barrier(PETSC_COMM_WORLD);

	// allocate memory
	x.resize(totalPoints);
	y.resize(totalPoints);

	// broadcast vectors to all processes
	MPI_Bcast(&x.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&y.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD);

	bodyGlobalIndices.resize(x.size());

	PetscPrintf(PETSC_COMM_WORLD, "Number of body points: %d\n", x.size());
}

template <>
void TairaColoniusSolver<3>::initialiseBodies()
{
	PetscInt rank, numProcs;
	PetscInt totalPoints;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

	numBoundaryPointsOnProcess.resize(numProcs);
	boundaryPointIndices.resize(numProcs);
	startGlobalIndices.resize(numProcs);

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
	/*	else if (type == "points")
		{
			string fname;
			node["pointsFile"] >> fname;
			fname = "bodyFiles/" + fname;
			std::cout << fname << std::endl;
			// initialise points
			std::ifstream file(fname.c_str());
			file >> Body.numPoints;
			Body.X.resize(Body.numPoints);
			Body.Y.resize(Body.numPoints);
			for(int i=0; i<Body.numPoints; i++)
			{
				file >> Body.X[i] >> Body.Y[i];
			}
			file.close();
		}
		else if (type == "lineSegment")
		{
			real startX, startY, endX, endY;
			int numPoints;
			node["segmentOptions"][0] >> startX;
			node["segmentOptions"][1] >> endX;
			node["segmentOptions"][2] >> startY;
			node["segmentOptions"][3] >> endY;
			node["segmentOptions"][4] >> numPoints;
			Body.numPoints = numPoints;
			// initialise line segment
		}*/
		else
			std::cout << "[E]: unknown Body type\n";
		totalPoints = x.size();
	}

	// broadcast total number of body points to all processes
	MPI_Bcast(&totalPoints, 1, MPIU_INT, 0, PETSC_COMM_WORLD);

	MPI_Barrier(PETSC_COMM_WORLD);

	// allocate memory
	x.resize(totalPoints);
	y.resize(totalPoints);
	z.resize(totalPoints);

	// broadcast vectors to all processes
	MPI_Bcast(&x.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&y.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&z.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD);

	bodyGlobalIndices.resize(x.size());

	PetscPrintf(PETSC_COMM_WORLD, "Number of body points: %d\n", x.size());
}
/***************************************************************************//**
 * \file TairaColoniusSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class `TairaColoniusSolver`.
 */


#if !defined(TAIRA_COLONIUS_SOLVER_H)
#define TAIRA_COLONIUS_SOLVER_H

#include "navierStokes/NavierStokesSolver.h"

#include "Body.h"


/**
 * \class TairaColoniusSolver
 * \brief Solves the Navier-Stokes equations 
 *        with immersed boundary projection method (Taira and Colonius, 2007).
 */
template <PetscInt dim>
class TairaColoniusSolver : public NavierStokesSolver<dim>
{
public:
  DM bda; ///< DMDA object

  PetscInt numBodies;
  std::vector<Body<dim> > bodies; ///< info about each body immersed

  PetscInt numBoundaryPoints; ///< total numbers of Lagrangian body points

  PetscInt  startGlobalIndex;
  Mat       ET;
  PetscReal force[3];
  Vec       nullSpaceVec, regularizedForce;

  std::ofstream forcesFile;

  std::vector<PetscInt>  globalIndexMapping;
  std::vector<PetscInt>  numBoundaryPointsOnProcess;
  std::vector<PetscInt>  numPhiOnProcess;
  std::vector< std::vector<PetscInt> > boundaryPointIndices;
  
  PetscErrorCode initializeBodies();
  PetscErrorCode generateBodyInfo();
  PetscErrorCode calculateCellIndices();
  PetscErrorCode createDMs();
  PetscErrorCode createVecs();
  PetscErrorCode setNullSpace();
  PetscErrorCode generateBNQ();
  PetscErrorCode generateR2();
  PetscErrorCode createGlobalMappingBodies();
  PetscErrorCode calculateForces();
  PetscErrorCode calculateForcesTC();
  
  PetscErrorCode readLambda(std::string directory);
  PetscErrorCode writeData();
  PetscErrorCode writeLambda(std::string directory);
  PetscErrorCode writeForces();

  PetscReal dhRoma(PetscReal x, PetscReal h);
  PetscReal delta(PetscReal x, PetscReal y, PetscReal h);
  PetscReal delta(PetscReal x, PetscReal y, PetscReal z, PetscReal h);
  PetscBool isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal xBody, PetscReal yBody, PetscReal radius, PetscReal *delta);
  PetscBool isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal zGrid, PetscReal xBody, PetscReal yBody, PetscReal zBody, PetscReal radius, PetscReal *delta);

public:
  // constructors
  TairaColoniusSolver(){ };
  TairaColoniusSolver(CartesianMesh *cartesianMesh, 
                      FlowDescription<dim> *flowDescription, 
                      SimulationParameters *simulationParameters);
  // destructor
  ~TairaColoniusSolver(){ };
  PetscErrorCode initialize();
  PetscErrorCode finalize();

}; // TairaColoniusSolver

#endif
/*! Definition of the class `TairaColoniusSolver`.
 * \file TairaColoniusSolver.h
 */


#if !defined(TAIRA_COLONIUS_SOLVER_H)
#define TAIRA_COLONIUS_SOLVER_H

#include "navierStokes/NavierStokesSolver.h"

#include "Body.h"


/*!
 * \class TairaColoniusSolver
 * \brief Solves the Navier-Stokes equations 
 *        with immersed boundary projection method (Taira and Colonius, 2007).
 */
template <PetscInt dim>
class TairaColoniusSolver : public NavierStokesSolver<dim>
{
public:
  PetscInt numBodies;  ///< number of immersed boundaries
  PetscInt numLagrangianPoints; ///< total number of Lagrangian points
  std::vector<Body<dim> > bodies;  ///< info about each immersed boundary
  
  DM bda;  ///< DMDA for all immersed boundaries

  std::vector<PetscInt> localNumPhiPoints;
  std::vector<PetscInt>  globalIndexMapping;

  Vec nullSpaceVec;

  std::ofstream forcesFile;
  
  PetscErrorCode initializeBodies();
  PetscErrorCode setLocalIndexPointsBodies();
  PetscErrorCode setLocalNumPhiPoints();
  PetscErrorCode createDMs();
  PetscErrorCode createVecs();
  PetscErrorCode setNullSpace();
  PetscErrorCode generateBNQ();
  PetscErrorCode generateR2();
  PetscErrorCode createGlobalMappingBodies();
  PetscErrorCode calculateForces();
  
  PetscErrorCode readLambda(std::string directory);
  PetscErrorCode writeData();
  PetscErrorCode writeLambda(std::string directory);
  PetscErrorCode writeForces();

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

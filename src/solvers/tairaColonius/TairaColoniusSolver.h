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
  PetscInt numBodies; ///< number of immersed boundaries
  std::vector<Body<dim> > bodies; ///< info about each immersed boundary
  
  DM bda; ///< DMDA object for all immersed boundaries

  Vec nullSpaceVec; ///< nullspace object to attach to the matrix QTBNQ

  std::ofstream forcesFile; ///< stream the file containing the forces acting on each immersed boundary
  
  PetscErrorCode initializeBodies();
  PetscErrorCode getNumLagPoints(PetscInt &n);
  PetscErrorCode getNumLagPointsOnProcess(std::vector<PetscInt> &numOnProcess);
  PetscErrorCode registerLagPointsOnProcess();
  PetscErrorCode createDMs();
  PetscErrorCode createVecs();
  PetscErrorCode createGlobalMappingBodies();
  PetscErrorCode generateBNQ();
  PetscErrorCode generateR2();
  PetscErrorCode setNullSpace();
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

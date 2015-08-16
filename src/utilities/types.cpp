/***************************************************************************//**
 * \file types.cpp
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the functions to map enums to strings.
 */


#include "types.h"

#include <iostream>


/**
 * \brief Returns the boundary location as an enum.
 */
BoundaryLocation stringToBoundaryLocation(std::string s)
{
  if (s == "xMinus" || s == "left")
    return XMINUS;
  if (s == "xPlus" || s == "right")
    return XPLUS;
  if (s == "yMinus" || s == "bottom")
    return YMINUS;
  if (s == "yPlus" || s == "top")
    return YPLUS;
  if (s == "zMinus" || s == "back")
    return ZMINUS;
  if (s == "zPlus" || s == "front")
    return ZPLUS;
  std::cout << "\nERROR: " << s << " - unknown boundary location.\n";
  std::cout << "Boundary locations available:\n";
  std::cout << "\txMinus (or left)\n";
  std::cout << "\txPlus (or right)\n";
  std::cout << "\tyMinus (or bottom)\n";
  std::cout << "\tyPlus (or top)\n";
  std::cout << "\tzMinus (or back)\n";
  std::cout << "\tzPlus (or front)\n" << std::endl;
  exit(0);
} // stringToBoundaryLocation


/**
 * \brief Returns the boundary location as a string.
 */
std::string stringFromBoundaryLocation(BoundaryLocation location)
{
  switch(location)
  {
    case XMINUS:
      return "xMinus (left)";
      break;
    case XPLUS:
      return "xPlus (right)";
      break;
    case YMINUS:
      return "yMinus (bottom)";
      break;
    case YPLUS:
      return "yPlus (top)";
      break;
    case ZMINUS:
      return "zMinus (back)";
      break;
    case ZPLUS:
      return "zPlus (front)";
      break;
    default:
      break;
  }
} // stringFromBoundaryLocation


/**
 * \brief Returns the boundary condition type as a enum.
 */
BoundaryType stringToBoundaryType(std::string s)
{
  if (s == "DIRICHLET")
    return DIRICHLET;
  if (s == "NEUMANN")
    return NEUMANN;
  if (s == "CONVECTIVE")
    return CONVECTIVE;
  if (s == "PERIODIC")
    return PERIODIC;
  std::cout << "\nERROR: " << s << " - unknown boundary condition type.\n";
  std::cout << "Boundary condition types available:\n";
  std::cout << "\tDIRICHLET\n";
  std::cout << "\tNEUMANN\n";
  std::cout << "\tCONVECTIVE\n";
  std::cout << "\tPERIODIC\n" << std::endl;
  exit(0);
} // stringToBoundaryType


/**
 * \brief Returns the boundary condition type as a string.
 */
std::string stringFromBoundaryType(BoundaryType type)
{
  switch(type)
  {
    case DIRICHLET:
      return "DIRICHLET";
      break;
    case NEUMANN:
      return "NEUMANN";
      break;
    case CONVECTIVE:
      return "CONVECTIVE";
      break;
    case PERIODIC:
      return "PERIODIC";
      break;
    default:
    break;
  }
} // stringFromBoundaryType


/**
 * \brief Returns the immersed-boundary method as an enum.
 */
IBMScheme stringToIBMScheme(std::string s)
{
  if (s == "NAVIER_STOKES")
    return NAVIER_STOKES;
  if (s == "TAIRA_COLONIUS")
    return TAIRA_COLONIUS;
  std::cout << "\nERROR: " << s << " - unknown IBM.\n";
  std::cout << "IBMs available:\n";
  std::cout << "\tNAVIER_STOKES\n";
  std::cout << "\tTAIRA_COLONIUS\n" << std::endl;
  exit(0);
} // stringToIBMScheme


/**
 * \brief Returns the immersed-boundary method as a string.
 */
std::string stringFromIBMScheme(IBMScheme ibmScheme)
{
  switch(ibmScheme)
  {
    case NAVIER_STOKES:
      return "Navier-Stokes (Perot, 1993)";
      break;
    case TAIRA_COLONIUS:
      return "Immersed-Boundary Projection method (Taira and Colonius, 2007)";
      break;
    default:
      break;
  }
} // stringFromIBMScheme


/**
 * \brief Returns the iterative method to use as an enum.
 */
IterativeMethod stringToIterativeMethod(std::string s)
{
  if (s == "CG")
    return CG;
  if (s == "BCGS")
    return BCGS;
  if (s == "GMRES")
    return GMRES;
  std::cout << "\nERROR: " << s << " - unknown iterative method.\n";
  std::cout << "Iterative methods implemented:\n";
  std::cout << "\tCG\n";
  std::cout << "\tBCGS\n";
  std::cout << "\tGMRES\n" << std::endl;
  exit(0);
} // stringToIterativeMethod


/**
 * \brief Returns the iterative method as a string.
 */
std::string stringFromIterativeMethod(IterativeMethod method)
{
  switch(method)
  {
    case CG:
      return "Conjugate-Gradient";
      break;
    case BCGS:
      return "biCGStab";
      break;
    case GMRES:
      return "GMRES";
      break;
    default:
      break;
  }
} // stringFromIterativeMethod


/**
 * \brief Returns the type of preconditioner as an enum.
 */
PreconditionerType stringToPreconditionerType(std::string s)
{
  if (s == "NONE")
    return NONE;
  if (s == "DIAGONAL")
    return DIAGONAL;
  if (s == "GAMG")
    return GAMG;
  std::cout << "\nERROR: " << s << " - unknown type of preconditioner.\n";
  std::cout << "Preconditioners available:\n";
  std::cout << "\tNONE\n";
  std::cout << "\tDIAGONAL\n";
  std::cout << "\tGAMG\n" << std::endl;
  exit(0);
} // stringToPreconditionerType


/**
 * \brief Returns the type of preconditioner as a string.
 */
std::string stringFromPreconditionerType(PreconditionerType preconditioner)
{
  switch(preconditioner)
  {
    case NONE:
      return "no preconditioning";
      break;
    case DIAGONAL:
      return "diagonal";
      break;
    case GAMG:
      return "GAMG - Smooth aggregation";
      break;
    default:
      break;
  }
} // stringFromPreconditionerType


/**
 * \brief Returns the time-integration scheme as an enum.
 */
TimeScheme stringToTimeScheme(std::string s)
{
  if (s == "EULER_EXPLICIT")
    return EULER_EXPLICIT;
  if (s == "EULER_IMPLICIT")
    return EULER_IMPLICIT;
  if (s == "ADAMS_BASHFORTH_2")
    return ADAMS_BASHFORTH_2;
  if (s == "CRANK_NICOLSON")
    return CRANK_NICOLSON;
  std::cout << "\nERROR: " << s << " - unknown time-integration scheme.\n";
  std::cout << "Time-integration schemes available:\n";
  std::cout << "\tEULER_EXPLICIT\n";
  std::cout << "\tEULER_IMPLICIT\n";
  std::cout << "\tADAMS_BASHFORTH_2\n";
  std::cout << "\tCRANK_NICOLSON\n" << std::endl;
  exit(0);
} // stringToTimeScheme


/**
 * \brief Returns the time-integration scheme as a string.
 */
std::string stringFromTimeScheme(TimeScheme scheme)
{
  switch(scheme)
  {
    case EULER_EXPLICIT:
      return "Euler-explicit";
      break;
    case EULER_IMPLICIT:
      return "Euler-implicit";
      break;
    case ADAMS_BASHFORTH_2:
      return "second-order Adams-Bashforth";
      break;
    case CRANK_NICOLSON:
      return "second-order Crank-Nicolson";
      break;
    default:
      break;
  }
} // stringFromTimeScheme
/***************************************************************************//**
 * \file types.cpp
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the functions to map enums to strings.
 */


#include "types.h"

#include <iostream>
#include <stdlib.h>


/**
 * \brief Returns the boundary location as an enum.
 *
 * \param s string that describes the boundary location.
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
  exit(1);
} // stringToBoundaryLocation


/**
 * \brief Returns the boundary location as a string.
 *
 * \param location boundary location as an enum.
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
      return "ERROR";
      break;
  }
} // stringFromBoundaryLocation


/**
 * \brief Returns the boundary condition type as a enum.
 *
 * \param s string that describes the type of boundary condition.
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
  exit(1);
} // stringToBoundaryType


/**
 * \brief Returns the boundary condition type as a string.
 *
 * \param type boundary condition as an enum.
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
      return "ERROR";
      break;
  }
} // stringFromBoundaryType


/**
 * \brief Returns the immersed-boundary method as an enum.
 *
 * \param s string that describes the immersed-boundary method.
 */
IBMethod stringToIBMethod(std::string s)
{
  if (s == "NONE")
    return NAVIER_STOKES;
  if (s == "TAIRA_COLONIUS")
    return TAIRA_COLONIUS;
  std::cout << "\nERROR: " << s << " - unknown Immersed Boudnary Method.\n";
  std::cout << "Immersed Boundary methods implemented:\n";
  std::cout << "\tTAIRA_COLONIUS\n";
  std::cout << "\tNONE\n" << std::endl;
  exit(1);
} // stringToIBMethod


/**
 * \brief Returns the immersed-boundary method as a string.
 *
 * \param method immersed-boundary method as an enum.
 */
std::string stringFromIBMethod(IBMethod method)
{
  switch(method)
  {
    case TAIRA_COLONIUS:
      return "Immersed-Boundary Projection method (Taira and Colonius, 2007)";
      break;
    case NAVIER_STOKES:
      return "Navier-Stokes solver (Perot, 1993)";
      break;
    default:
      return "ERROR";
      break;
  }
} // stringFromIBMethod


/**
 * \brief Returns the time-integration scheme as an enum.
 *
 * \param s string that describes the time-integration scheme.
 */
TimeScheme stringToTimeScheme(std::string s)
{
  if (s == "NONE")
    return NONE;
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
  std::cout << "\tNONE\n";
  std::cout << "\tEULER_EXPLICIT\n";
  std::cout << "\tEULER_IMPLICIT\n";
  std::cout << "\tADAMS_BASHFORTH_2\n";
  std::cout << "\tCRANK_NICOLSON\n" << std::endl;
  exit(1);
} // stringToTimeScheme


/**
 * \brief Returns the time-integration scheme as a string.
 * 
 * \param scheme time-integration scheme as an enum.
 */
std::string stringFromTimeScheme(TimeScheme scheme)
{
  switch(scheme)
  {
    case NONE:
      return "none";
      break;
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
      return "ERROR";
      break;
  }
} // stringFromTimeScheme


/**
 * \brief Returns the executing space as an enum.
 *
 * \param s string that describes the executing space of solvers.
 */
ExecuteType stringToExecuteType(std::string s)
{
  if (s == "GPU")
    return GPU;
  if (s == "CPU")
    return CPU;
  std::cout << "\nERROR: " << s << " - unknown executing space.\n";
  std::cout << "Acceptable executing spaces:\n";
  std::cout << "\tGPU\n";
  std::cout << "\tCPU\n" << std::endl;
  exit(EXIT_FAILURE);
} // stringToExecuteType


/**
 * \brief Returns the executing space as a string.
 *
 * \param exeType executing space as an enum.
 */
std::string stringFromExecuteType(ExecuteType exeType)
{
  switch(exeType)
  {
    case GPU:
      return "GPU-based (currently, AmgX)";
      break;
    case CPU:
      return "CPU-based (currently, PETSc KSPs)";
      break;
    default:
      return "ERROR";
      break;
  }
} // stringFromExecuteType

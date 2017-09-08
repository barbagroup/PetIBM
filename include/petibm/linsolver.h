/*! Implementation of the class `LinSolver`.
 * \file linsolver.h
 */


# pragma once

// STL
# include <string>
# include <memory>

// PETSc
# include <petscsys.h>
# include <petscvec.h>
# include <petscmat.h>
# include <petscviewer.h>

// PetIBM
# include "types.h"


namespace petibm
{
namespace linsolvers
{

/**
 * \class LinSolver.
 * \brief Super class for an iterative solver.
 */
class LinSolver
{
public:

    /** \brief default constructor. */
    LinSolver() = default;

    /**
     * \brief constructor.
     *
     * \param solverName the name of the solver.
     * \param file the full path to configuration file of the solver.
     */
    LinSolver(const std::string &solverName, const std::string &file):
        name(solverName), options(file) {};

    /** \brief virtual destruction function. */
    virtual ~LinSolver() = default;

    /**
     * \brief the function to set the coefficient matrix A in Ax=b.
     *
     * \param A the coefficient matrix.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode setMatrix(const Mat &A) = 0;

    /**
     * \brief the function to solve Ax=b.
     *
     * \param x the vector of unknowns.
     * \param b the RHS vector.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode solve(Vec &x, Vec &b) = 0;

    /**
     * \brief the function to get the number of iterations.
     *
     * \param iters the returned number of iterations.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getIters(PetscInt &iters) = 0;

    /**
     * \brief print info.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode printInfo(
            PetscViewer viewer=PETSC_VIEWER_STDOUT_WORLD) = 0;

protected:

    /** \brief the name of this solver.
     *
     * This name will also be used as the prefix foe PETSc solver options. For
     * example, if name is "velocity", the prefix in PETSc solver options will
     * be "velocity_".
     */
    std::string     name;

    /** \brief the full path to the configuration file for this solver. */
    std::string     options;

    /**
     * \brief the private initialization function.
     *
     * The design of linear solver classes is that only the factory function can
     * initialize the instances, and it doesn't allow user to declare instance
     * first and initialize the instance later. So it seems not necessary to 
     * make initialization function public.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode init() = 0;

}; // LinSolver


/**
 * \brief a factory function for creating a ponter to linear solver instance.
 *
 * \param solverName the name for the solver.
 * \param configFile the full path to the configuration file for the sovler.
 * \param type execution type.
 * \param solver a pointer to the solver instance.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createLinSolver(const std::string &solverName, 
        const std::string &configFile,
        const petibm::utilities::types::ExecuteType &type,
        std::shared_ptr<LinSolver> &solver);

} // end of namespace linsolvers
} // end of namespace petibm

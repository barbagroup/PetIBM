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

// YAML CPP
# include <yaml-cpp/yaml.h>

// PetIBM
# include <petibm/type.h>


namespace petibm
{
namespace linsolver
{
/**
 * \class LinSolverBase.
 * \brief Super class for an iterative solver.
 */
class LinSolverBase
{
public:

    /** \brief default constructor. */
    LinSolverBase() = default;

    /**
     * \brief constructor.
     *
     * \param solverName the name of the solver.
     * \param file the full path to configuration file of the solver.
     */
    LinSolverBase(const std::string &solverName, const std::string &file):
        name(solverName), config(file) {};

    /** \brief virtual destruction function. */
    virtual ~LinSolverBase() = default;

    /**
     * \brief Manually destroy the instance.
     *
     * \return  PetscErrorCode.
     */
    virtual PetscErrorCode destroy();
    
    /**
     * \brief print information to standard output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;
    

    /**
     * \brief get the type of this instance.
     * 
     * \param type [out] a string representing the type.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getType(std::string &type) const;

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
     * \param iters [out] the returned number of iterations.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getIters(PetscInt &iters) = 0;

    /**
     * \brief the function to get final residual.
     *
     * \param res[out] the returned residual (norm).
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getResidual(PetscReal &res) = 0;

protected:

    /** \brief the name of this solver.
     *
     * This name will also be used as the prefix foe PETSc solver options. For
     * example, if name is "velocity", the prefix in PETSc solver options will
     * be "velocity_".
     */
    std::string     name;

    /** \brief the full path to the configuration file for this solver. */
    std::string     config;
    
    /** \brief type of this linear solver. */
    std::string     type;

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
}


namespace type
{
    /** \brief type definition of LinSolver. */
    typedef std::shared_ptr<linsolver::LinSolverBase> LinSolver;
}


namespace linsolver
{
    /**
     * \brief a factory function for creating LinSolver.
     *
     * \param name [in] the name for the solver.
     * \param config [in] the full path to the configuration file of the sovler.
     * \param type [in] execution type.
     * \param solver [out] a LinSolver.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createLinSolver(const std::string &solverName,
            const YAML::Node &node, type::LinSolver &solver);
}

} // end of namespace petibm

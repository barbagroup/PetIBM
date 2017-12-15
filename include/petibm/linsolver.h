/** 
 * \file linsolver.h
 * \brief Def. of LinSolverBase, LinSolver, and factory function.
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
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

/** 
 * \defgroup linsolver Linear solvers
 * 
 * The design of PetIBM is to use shared pointer (std::shared_ptr) to hold 
 * instances. The class LinSolverBase is the abstract class for all kinds of 
 * linear solver classes. With an abstract class, we can unified the interfaces
 * of different linear solver classes. Users should use the factory function
 * petibm::linsolver::createLinSolver to create an instance, instead of 
 * initializing the instance directly.
 * 
 * Currently, there are two different linear solvers: PETSc KSP and NVIDIA AmgX.
 * Please see petibm::linsolver::createLinSolver for how to create different
 * types of linear solver instances.
 * 
 * The shared pointer holding a linear solver instance is defined as 
 * petibm::type::LinSolver.
 * 
 * \ingroup petibm
 */

namespace petibm
{
    
/** 
 * \namespace petibm::linsolver
 * \brief Collection of linear solvers from different libraries. 
 * 
 * \see petibm::type::LinSolver, petibm::linsolver::createLinSolver.
 * \ingroup linsolver
 */
namespace linsolver
{
    
/**
 * \class LinSolverBase
 * \brief The abstract (base) class for different iterative solvers.
 * \see petibm::type::LinSolver, petibm::linsolver::createLinSolver.
 * \ingroup linsolver
 */
class LinSolverBase
{
public:

    /** \brief Default constructor. */
    LinSolverBase() = default;

    /**
     * \brief Constructor.
     *
     * \param solverName [in] the name of the solver.
     * \param file [in] the path to the configuration file of the underlying solver.
     */
    LinSolverBase(const std::string &solverName, const std::string &file):
        name(solverName), config(file) {};

    /** \brief Default destructor. */
    virtual ~LinSolverBase() = default;

    /**
     * \brief Manually destroy the data in the current instance.
     *
     * \return  PetscErrorCode.
     */
    virtual PetscErrorCode destroy();
    
    /**
     * \brief Print information to standard output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;
    

    /**
     * \brief Get the type of this instance.
     * 
     * \param type [out] a string representing the type.
     * 
     * The possible return of `type` is either `NVIDIA AmgX` or `PETSc KSP`.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getType(std::string &type) const;

    /**
     * \brief Set the coefficient matrix A as in Ax=b.
     *
     * \param A [in] a PETSc Mat; the coefficient matrix.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode setMatrix(const Mat &A) = 0;

    /**
     * \brief Solve Ax=b.
     *
     * \param x [in, out] a PETSc Vec; the vector of unknowns.
     * \param b [in] a PETSc Vec; the RHS vector.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode solve(Vec &x, Vec &b) = 0;

    /**
     * \brief Get the number of iterations of the last solve.
     *
     * \param iters [out] number of iterations.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getIters(PetscInt &iters) = 0;

    /**
     * \brief Get the final residual of the last solve.
     *
     * \param res [out] residual (norm).
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getResidual(PetscReal &res) = 0;

protected:

    /** 
     * \brief The name of this solver.
     *
     * This name will also be used as the prefix for PETSc solver options. For
     * example, if name is "velocity", the prefix in PETSc solver options will
     * be "velocity_".
     */
    std::string     name;

    /** \brief The path to the configuration file for the underlying solver. */
    std::string     config;
    
    /** \brief The type of the underlying linear solver. */
    std::string     type;

    /**
     * \brief Private initialization function.
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
    /** \brief Type definition of LinSolver. 
     * 
     * Please use petibm::linsolver::createLinSolver to create a LinSolver 
     * instance.
     * 
     * Exmaple Usage:
     * \code
     * PetscErrorCode ierr;
     * Mat A;
     * Vec lhs, rhs;
     * petibm::type::LinSolver solver;
     * 
     * // create LinSolver with petibm::linsolver::createLinSolver
     * // create matrix A, vectors lhs and rhs. See PETSc manual.
     * 
     * ierr = solver->setMatrix(A); CKERRQ(ierr);
     * ierr = solver->solve(lhs, rhs); CHKERRQ(ierr);
     * \endcode
     * 
     * \see petibm::linsolver::createLinSolver, petibm::linsolver::LinSolverBase.
     * \ingroup linsolver
     */
    typedef std::shared_ptr<linsolver::LinSolverBase> LinSolver;
}


namespace linsolver
{
    /**
     * \brief A factory function for creating a LinSolver.
     *
     * \param solverName [in] the name for the solver.
     * \param node [in] the YAML::Node containing all configuration.
     * \param solver [out] a LinSolver.
     * 
     * This function will look into the `parameters` and find the child key 
     * `<solverName>Solver`. Under `<solverName>Solver`, there should at least
     * be two child keys: `type` and `config`. `type` indicate the type of the
     * solver, and `config` provides the path to the the solver configuration 
     * file.
     * 
     * For example, if the name of the target LinSolver is **velocity**, this
     * factory function will use `node[parameters][velocitySolver][type]`and
     * `node[parameters][velocitySolver][config]` to determine the type and the
     * path to the configuration file for the underlying linear solver.
     * 
     * Currently in PetIBM, the key `type` only accepts `CPU` and `GPU`.
     *  
     * An example of creating a LinSolver instance with KSP:
     * \code
     * YAML::Node node;
     * node["parameters"]["velocitySolver"]["type"] = "CPU";
     * node["parameters"]["velocitySolver"]["config"] = "solversPetscOptions.info";
     * 
     * PetscErrorCode ierr;
     * petibm::type::LinSolver solver;
     * ierr = petibm::linsolver::createLinSolver("velocity", node, solver); CHKERRQ(ierr);
     * \endcode
     *
     * \return PetscErrorCode.
     * 
     * \see petibm::type::LinSolver
     * \ingroup linsolver
     */
    PetscErrorCode createLinSolver(const std::string &solverName,
            const YAML::Node &node, type::LinSolver &solver);
}

} // end of namespace petibm

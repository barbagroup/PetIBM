/**
 * \file linsolver.h
 * \brief Def. of LinSolverBase, LinSolver, and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <memory>
#include <string>

#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>
#include <yaml-cpp/yaml.h>

#include <petibm/type.h>

/**
 * \defgroup linsolver Linear solvers
 * \brief Interfaces to different libraries of linear solvers
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
 *
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
     * \param solverName [in] Name of the solver.
     * \param file [in] Path of the configuration file for the solver.
     */
    LinSolverBase(const std::string &solverName, const std::string &file)
        : name(solverName), config(file){};

    /** \brief Default destructor. */
    virtual ~LinSolverBase() = default;

    /** \brief Manually destroy the data in the current instance. */
    virtual PetscErrorCode destroy();

    /** \brief Print information to standard output. */
    PetscErrorCode printInfo() const;

    /**
     * \brief Get the type of this instance.
     *
     * \param _type [out] String representing the type.
     *
     * Possible returns for `_type` are `NVIDIA AmgX` or `PETSc KSP`.
     */
    PetscErrorCode getType(std::string &_type) const;

    /**
     * \brief Set the coefficient matrix of the linear system.
     *
     * \param A [in] Coefficient matrix.
     */
    virtual PetscErrorCode setMatrix(const Mat &A) = 0;

    /**
     * \brief Solve the linear system.
     *
     * \param x [in, out] Vector of unknowns.
     * \param b [in] Right-hand side of the linear system.
     */
    virtual PetscErrorCode solve(Vec &x, Vec &b) = 0;

    /**
     * \brief Get the number of iterations of the solver.
     *
     * \param iters [out] Number of iterations.
     */
    virtual PetscErrorCode getIters(PetscInt &iters) = 0;

    /**
     * \brief Get the final residual of the solver.
     *
     * \param res [out] Residual.
     */
    virtual PetscErrorCode getResidual(PetscReal &res) = 0;

protected:
    /**
     * \brief Name of the linear solver.
     *
     * The name should match the prefix used for to configure the solver
     * from the command line or from a configuration file.
     * For example, if the name is "velocity", the prefix should be "_velocity".
     */
    std::string name;

    /** \brief Path of the solver configuration file to read. */
    std::string config;

    /** \brief Type of the linear solver. */
    std::string type;

    /**
     * \brief Private initialization function.
     *
     * The design of linear solver classes is that only the factory function can
     * initialize the instances, and it doesn't allow user to declare instance
     * first and initialize the instance later. So it seems not necessary to
     * make initialization function public.
     */
    virtual PetscErrorCode init() = 0;

};  // LinSolver

}  // end of namespace linsolver

namespace type
{
/** \brief Type definition of LinSolver.
 *
 * Use petibm::linsolver::createLinSolver to create a LinSolver object.
 *
 * Example Usage:
 * \code
 * PetscErrorCode ierr;
 * Mat A;
 * Vec lhs, rhs;
 * petibm::type::LinSolver solver;
 *.
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

}  // end of namespace type

namespace linsolver
{
/**
 * \brief A factory function for creating a LinSolver.
 *
 * \param solverName [in] Name of the linear solver.
 * \param node [in] YAML configuration node.
 * \param solver [out] The linear solver.
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
 * ierr = petibm::linsolver::createLinSolver(
 *     "velocity", node, solver); CHKERRQ(ierr);
 * \endcode
 *
 * \return PetscErrorCode.
 *
 * \see petibm::type::LinSolver
 * \ingroup linsolver
 */
PetscErrorCode createLinSolver(const std::string &solverName,
                               const YAML::Node &node, type::LinSolver &solver);

}  // end of namespace linsolver

}  // end of namespace petibm

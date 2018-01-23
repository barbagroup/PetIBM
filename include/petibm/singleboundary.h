/**
 * \file singleboundary.h
 * \brief Definition of the class `SingleBoundaryBase`.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


# pragma once

// C++ STL
# include <memory>

// here goes PETSc headers
# include <petscsys.h>
# include <petscvec.h>

// here goes headers from our PetIBM
# include <petibm/type.h>
# include <petibm/mesh.h>


namespace petibm
{
namespace boundary
{

/**
 * \brief Base (abstract) class for ghost points & BC on a single boundary.
 * \see boundaryModule, petibm::type::SingleBoundary
 * \ingroup boundaryModule
 */
class SingleBoundaryBase
{
public:

    /** \brief Dimension. */
    PetscInt            dim;

    /** \brief The location of this boundary. */
    type::BCLoc         loc;
    
    /** \brief The field of which the ghost points represent. */
    type::Field         field;
    
    /** \brief The type of boundary conditions. */
    type::BCType        type;

    /** \brief A constant value representing BC value. */
    PetscReal           value;

    /** \brief The direction of normal vector. */
    PetscReal           normal;

    /** \brief The list of ghost points on this boundary and at this field. */
    type::GhostPointsList   points;

    /** \brief Indicate if this process holds part of this boundary. */
    PetscBool           onThisProc;



    /**
     * \brief Constructor.
     * \param mesh [in] a Mesh instance.
     * \param loc [in] the location of the target boundary.
     * \param field [in] the target field.
     * \param type [in] the type of BC.
     * \param value [in] BC value.
     */
    SingleBoundaryBase(const type::Mesh &mesh,
            const type::BCLoc &loc, const type::Field &field,
            const type::BCType &type, const PetscReal &value); 

    /** \brief Default constructor. */
    SingleBoundaryBase() = default;

    /** \brief Default destructor. */
    virtual ~SingleBoundaryBase();


    /**
     * \brief Manually destroy data.
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode destroy();


    /**
     * \brief Set up the initial values of the ghost points.
     * \param vec [in] a packed solution Vec containing initial values.
     * \return PetscErrorCode.
     */
    PetscErrorCode setGhostICs(const Vec &vec);

    /**
     * \brief Modify the coefficients in the equation of ghost points.
     * \param vec [in] a packed solution Vec at current time step.
     * \param dt [in] a PetscReal representing time-step size.
     * \return PetscErrorCode.
     * 
     * The equation of ghost points means the relationship between boundary 
     * points and ghost points. The equation has a form 
     * u_ghost = a0 * u_boundary + a1. This function changes a1 and a0 according
     * to the type of BC.
     */
    PetscErrorCode updateEqs(const Vec &vec, const PetscReal &dt);

    /**
     * \brief Update the values of ghost points using the equation.
     * \param vec [in] a packed solution Vec at current time step.
     * \return PetscErrorCode.
     */
    PetscErrorCode updateGhostValues(const Vec &vec);

    /**
     * \brief Copy the values of ghost points to a local Vec.
     * \param lclVec [in, out] a local Vec with memory allocations for ghost points.
     * \return PetscErrorCode.
     * 
     * In PetIBM, we use a global packed Vec for velocity fields. But in some
     * occasions, we may need local Vecs that have ghost points in them. So we
     * have to copy the ghost values from this instance to the local Vecs.
     */
    PetscErrorCode copyValues2LocalVec(Vec &lclVec);

protected:

    /**
     * \brief Underlying initialization function.
     * \param mesh [in] a Mesh instance.
     * \param loc [in] the location of the target boundary.
     * \param field [in] the target field.
     * \param type [in] the type of BC.
     * \param bcValue [in] BC value.
     * \return  PetscErrorCode.
     */
    PetscErrorCode init(const type::Mesh &mesh,
            const type::BCLoc &loc, const type::Field &field,
            const type::BCType &type, const PetscReal &bcValue); 

    /**
     * \brief The underlying kernel for setting initial values and equations.
     * \param targetValue [in] the value of the corresponding boundary point.
     * \param p [in, out] the target ghost point.
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode setGhostICsKernel(
            const PetscReal &targetValue, type::GhostPointInfo &p) = 0;

    /**
     * \brief Underlying kernel for updating the coefficients of the equation.
     * \param targetValue [in] the value of corresponding boundary point.
     * \param dt [in] the size of a time step.
     * \param p [in, out] the target ghost point.
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p) = 0;


    /** \brief MPI communicator. */
    MPI_Comm        comm;

    /** \brief The size of the MPI communicator. */
    PetscMPIInt     mpiSize;
    
    /** \brief The rank of this process. */
    PetscMPIInt     mpiRank;


    /** \brief The corresponding Mesh object. */
    type::Mesh      mesh;

}; // SingleBoundaryBase

} // end of namespace boundary


namespace type
{
    /**
     * \brief Definition of type petibm::type::SingleBoundary.
     * \see boundaryModule, petibm::boundary::SingleBoundaryBase, petibm::boundary::createSingleBoundary
     * \ingroup boundaryModule
     * 
     * Please use petibm::boundary::createSingleBoundary to create a Mesh object.
     */
    typedef std::shared_ptr<boundary::SingleBoundaryBase> SingleBoundary;
} // end of namespace type

namespace boundary
{
    /**
     * \brief Factory function for creating a SingleBoundary object.
     * \param mesh [in] a Mesh instance.
     * \param loc [in] the location of the target boundary.
     * \param field [in] the target field.
     * \param value [in] BC value.
     * \param bcType [in] Type of the BC.
     * \param singleBd [out] resulting SingleBoundary object.
     * \return PetscErrorCode.
     * \see boundaryModule, petibm::type::SingleBoundarym petibm::boundary::createBoundary
     * \ingroup boundaryModule
     * 
     * This is generally not used by API users. API users normally only need to
     * use petibm::boundary::createBoundary.
     */
    PetscErrorCode createSingleBoundary(
            const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &value,
            const type::BCType &bcType,
            type::SingleBoundary &singleBd);
} // end of namespace boundary

} // end of namespace petibm

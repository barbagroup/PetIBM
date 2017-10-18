/***************************************************************************//**
 * \file singleboundary.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SingleBoundaryBase`.
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

/** \brief abstract class for ghost points & BC on a single boundary. */
class SingleBoundaryBase
{
public:

    /** \brief dimension. */
    PetscInt            dim;

    /** \brief the location of this boundary. */
    type::BCLoc         loc;
    
    /** \brief the field of which the ghost points represent. */
    type::Field         field;

    /** \brief a constant value representing BC value. */
    PetscReal           value;

    /** \brief the direction of normal vector. */
    PetscReal           normal;

    /** \brief the list of ghost points on this boundary and at this field. */
    type::GhostPointsList   points;

    /** \brief indicate if this process holds part of this boundary. */
    PetscBool           onThisProc;



    /** \brief default constructor. */
    SingleBoundaryBase() = default;

    /**
     * \brief constructor.
     *
     * \param mesh [in] a Mesh instance.
     * \param loc [in] the location of the target boundary.
     * \param field [in] the target field.
     * \param value [in] BC value.
     */
    SingleBoundaryBase(const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &value); 

    /** \brief default destructor. */
    virtual ~SingleBoundaryBase() = default;


    /**
     * \brief set up the initial values of the ghost points.
     *
     * \param vec [in] a packed solution Vec containing initial values.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode setGhostICs(const Vec &vec);

    /**
     * \brief modify the coefficients in the equation of ghost points.
     * 
     * The equation of ghost points means the relationship between boundary 
     * points and ghost points. The equation has a form 
     * u_ghost = a1 * u_boundary + a0. This function changes a1 and a0 according
     * to the type of BC.
     *
     * \param vec [in] a packed solution Vec at current time step.
     * \param dt [in] a PetscReal representing time-step size.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode updateEqs(const Vec &vec, const PetscReal &dt);

    /**
     * \brief update the values of ghost points using the equation.
     *
     * \param vec [in] a packed solution Vec at current time step.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode updateGhostValues(const Vec &vec);

    /**
     * \brief copy the values of ghost points to a local Vec.
     * 
     * In PetIBM, we use a global packed Vec for velocity fields. But in some
     * occasions, we may need local Vecs that have ghost points in them. So we
     * have to copy the ghost values from this instance to the local Vesc.
     *
     * \param lclVec [in, out] a local Vec with memory allocations for ghost points.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode copyValues2LocalVec(Vec &lclVec);

protected:

    /**
     * \brief underlying initialization function.
     *
     * \param mesh [in] a Mesh instance.
     * \param loc [in] the location of the target boundary.
     * \param field [in] the target field.
     * \param value [in] BC value.
     *
     * \return  PetscErrorCode.
     */
    PetscErrorCode init(const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &bcValue); 

    /**
     * \brief the underlying kernel for setting initial values and equations.
     *
     * \param targetValuea [in] the value of the corresponding boundary point.
     * \param p [in, out] the target ghost point.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode setGhostICsKernel(
            const PetscReal &targetValue, type::GhostPointInfo &p) = 0;

    /**
     * \brief underlying kernel for updating the coefficients of the equation.
     *
     * \param targetValue [in] the value of corresponding boundary point.
     * \param dt [in] the size of a time step.
     * \param p [in, out] the target ghost point.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p) = 0;


    /** \brief MPI communicator. */
    MPI_Comm        comm;

    /** \brief the size of the MPI communicator. */
    PetscMPIInt     mpiSize;
    
    /** \brief the rank of this process. */
    PetscMPIInt     mpiRank;


    /** \brief the corresponding Mesh object. */
    type::Mesh      mesh;

};

} // end of namespace boundary


namespace type
{
    /** \brief definition of type SingleBoundary. */
    typedef std::shared_ptr<boundary::SingleBoundaryBase> SingleBoundary;
}

namespace boundary
{
    /**
     * \brief factory function for creating a SingleBoundary object.
     *
     * \param mesh [in] a Mesh instance.
     * \param loc [in] the location of the target boundary.
     * \param field [in] the target field.
     * \param value [in] BC value.
     * \param singleBd [out] resulting SingleBoundary object.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createSingleBoundary(
            const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &value,
            type::SingleBoundary &singleBd);
}

} // end of namespace petibm

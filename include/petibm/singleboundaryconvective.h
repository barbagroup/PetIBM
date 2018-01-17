/**
 * \file singleboundaryconvective.h
 * \brief Definition of the class `SingleBoundaryConvective`.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


# pragma once

// here goes headers from our PetIBM
# include <petibm/singleboundary.h>


namespace petibm
{
namespace boundary
{

/**
 * \brief An implementation of SingleBoundaryBase for convective BC.
 * \see boundaryModule, petibm::type::SingleBoundary, petibm::boundary::createSingleBoundary
 * \ingroup boundaryModule
 */
class SingleBoundaryConvective : public SingleBoundaryBase
{
public:

    /**
     * \brief Constructor.
     * \param mesh [in] a Mesh instance.
     * \param loc [in] the location of the target boundary.
     * \param field [in] the target field.
     * \param value [in] BC value.
     */
    SingleBoundaryConvective(const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &value); 

    /** \copydoc SingleBoundaryBase::~SingleBoundaryBase */
    virtual ~SingleBoundaryConvective() = default;

    // re-implementation of SingleBoundaryBase::destroy
    virtual PetscErrorCode destroy();

protected:

    // implementation of SingleBoundaryBase::setGhostICsKernel
    virtual PetscErrorCode setGhostICsKernel(
            const PetscReal &targetValue, type::GhostPointInfo &p);

    // implementation of SingleBoundaryBase::updateEqsKernel
    virtual PetscErrorCode updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p);

    /**
     * \brief Underlying kernel that will determined during runtime according to
     *        the location of the boundary and the target field in a staggered grid.
     * \param targetValue [in] The value at time \f$t-1\f$ of the corresponding 
     *        boundary point. 
     * \param dt [in] Time-step size.
     * \param normal [in] Normal of the boundary: \f$+1\f$ represent *PLUS faces, 
     *        while \f$-1\f$ means *MINUS faces.
     * \param value [in] The value of the boundary condition.
     * \param p [in, out] A petibm::type::GhostPointInfo.
     * \return PetscErrorCode.
     */
    PetscErrorCode (*kernel)(
            const PetscReal &targetValue, const PetscReal &dt, 
            const PetscReal &normal, const PetscReal &value,
            type::GhostPointInfo &p);
};

} // end of namespace boundary
} // end of namespace petibm

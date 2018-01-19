/**
 * \file singlebodypoints.h
 * \brief Definition of body::SingleBodyPoints.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


# pragma once

// PetIBM
# include <petibm/singlebody.h>


namespace petibm
{
namespace body
{
/**
 * \brief An implementation of body::SingleBodyBase that uses point data as input.
 * \see bodyModule, petibm::type::SingleBody, petibm::body::SingleBodyBase
 * \ingroup bodyModule 
 * 
 * This implementation uses an ASCII file of coordinates of Lagrangian points as
 * its input.
 * 
 * Users should not initialize an instance of this class directly. They should
 * use petibm::body::createSingleBody. Actually, the design of PetIBM is to
 * handle multiple bodies, users should use petibm::type::BodyPack to handle
 * their bodies even if there is only one body.
 */
class SingleBodyPoints : public SingleBodyBase
{
public:


    /**
     * \brief Constructor using CartesainMesh and input file.
     *
     * \param mesh [in] an instance of type::Mesh.
     * \param name [in] the name of this body.
     * \param file [in] the ASCII file containing coordinates of Lagrangian points.
     */
    SingleBodyPoints(const type::Mesh &mesh,
            const std::string &name, const std::string &file);


    /** \copydoc SingleBodyBase::~SingleBodyBase */
    virtual ~SingleBodyPoints() = default;


    // implementation of SingleBodyBase::findProc
    virtual PetscErrorCode findProc(const PetscInt &i, PetscMPIInt &p) const;


    // implementation of SingleBodyBase::getGlobalIndex
    virtual PetscErrorCode getGlobalIndex(
            const PetscInt &i, const PetscInt &dof, PetscInt &idx) const;


    // implementation of SingleBodyBase::getGlobalIndex
    virtual PetscErrorCode getGlobalIndex(
            const MatStencil &s, PetscInt &idx) const;


    // implementation of SingleBodyBase::calculateAvgForces
    virtual PetscErrorCode calculateAvgForces(
            const Vec &f, type::RealVec1D &fAvg) const;
    

    // implementation of SingleBodyBase::updateMeshIdx
    virtual PetscErrorCode updateMeshIdx();
    

protected:


    /**
     * \brief Underlying initialization function.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode init(const type::Mesh &mesh,
            const std::string &name, const std::string &file);
    
    
    /**
     * \brief Find the indices of pressure cells that own local Lagrangian points.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode findCellIdx();


    /**
     * \brief Create a parallel 1D DMDA for this body.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createDMDA();


    /**
     * \brief Create a string for printing information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();
}; // SingleBodyPoints

} // end of namespace body
} // end of namespace petibm

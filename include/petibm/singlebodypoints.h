/***************************************************************************//**
 * \file singlebody.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SingleBody`.
 */


# pragma once

// PetIBM
# include <petibm/singlebody.h>


namespace petibm
{
namespace body
{
/** \brief class for a single body. */
class SingleBodyPoints : public SingleBodyBase
{
public:


    SingleBodyPoints(const type::Mesh &mesh,
            const std::string &name, const std::string &file);


    /** \brief the default destructor. */
    virtual ~SingleBodyPoints() = default;


    virtual PetscErrorCode findProc(const PetscInt &i, PetscMPIInt &p) const;


    virtual PetscErrorCode getGlobalIndex(
            const PetscInt &i, const PetscInt &dof, PetscInt &idx) const;


    virtual PetscErrorCode getGlobalIndex(
            const MatStencil &s, PetscInt &idx) const;


    virtual PetscErrorCode calculateAvgForces(
            const Vec &f, type::RealVec1D &fAvg) const;
    

    virtual PetscErrorCode updateMeshIdx();
    

protected:


    /**
     * \brief underlying initialization function.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode init(const type::Mesh &mesh,
            const std::string &name, const std::string &file);
    
    
    /**
     * \brief find the indices of presure cells that own local Lagrangian points.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode findCellIdx();


    /**
     * \brief create a parallel 1D DMDA for this body.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createDMDA();


    /**
     * \brief create a string for printing information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();
};

} // end of namespace body
} // end of namespace petibm

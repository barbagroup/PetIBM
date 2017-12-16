/***************************************************************************//**
 * \file bodypack.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `BodyPack`.
 */


# pragma once


// STL
# include <string>
# include <vector>
# include <memory>

// PETSc
# include <petscsys.h>
# include <petscdm.h>

// YAML
# include <yaml-cpp/yaml.h>

// PetIBM
# include <petibm/singlebody.h>


namespace petibm
{
namespace body
{

/** \brief class for a pack of multiple bodies. */
class BodyPackBase
{
public:

    /** \brief dimension. */
    PetscInt                    dim;

    /** \brief number of bodies in this pack. */
    PetscInt                    nBodies;

    /** \brief total number of Lagrangian points. */
    PetscInt                    nPts;

    /** \brief total number of local Lagrangian points. */
    PetscInt                    nLclPts;

    /** \brief a vector of SingleBody instances. */
    std::vector<type::SingleBody>     bodies;

    /** \brief a DMComposite of DMs of all `SingleBody`s. */
    DM                          dmPack;

    /** \brief a string for printing information. */
    std::string                 info;


    /** \brief default constructor. */
    BodyPackBase() = default;


    /**
     * \brief constructor of using a type::Mesh and a YAML node.
     *
     * \param mesh [in] an instance of type::Mesh for background mesh.
     * \param node [in] a YAML node specifying information of all bodies.
     */
    BodyPackBase(const type::Mesh &mesh, const YAML::Node &node);


    /** \brief default destructor. */
    virtual ~BodyPackBase();


    /**
     * \brief manually destroy data.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode destroy();


    /**
     * \brief print information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;


    /**
     * \brief find which process owns the target Lagrangian point of target body.
     *
     * \param bIdx [in] index of target body.
     * \param ptIdx [in] index of target point.
     * \param proc [in] returned process index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode findProc(const PetscInt &bIdx, 
            const PetscInt &ptIdx, PetscMPIInt &proc) const;


    /**
     * \brief find un-packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx [in] index of target body.
     * \param ptIdx [in] index of target point.
     * \param dof [in] index of target DoF.
     * \param idx [in] returned un-packed global index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getGlobalIndex(const PetscInt &bIdx, 
            const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const;


    /**
     * \brief find un-packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx [in] index of target body.
     * \param s [in] MatStencil of target point.
     * \param idx [in] returned un-packed global index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getGlobalIndex(const PetscInt &bIdx, 
            const MatStencil &s, PetscInt &idx) const;



    /**
     * \brief find packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx [in] index of target body.
     * \param ptIdx [in] index of target point.
     * \param dof [in] index of target DoF.
     * \param idx [in] returned packed global index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getPackedGlobalIndex(const PetscInt &bIdx, 
            const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const;



    /**
     * \brief find packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx [in] index of target body.
     * \param s [in] MatStencil of target point.
     * \param idx [in] returned packed global index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getPackedGlobalIndex(const PetscInt &bIdx, 
            const MatStencil &s, PetscInt &idx) const;


    /**
     * \brief calculate the averaged force of each body.
     *
     * \param f [in] packed force Vec of Lagrangian points.
     * \param fAvg [in] return averaged force for each body.
     *
     * Note: fAvg doesn't have to have correct size. This function will resize
     * the fAvg.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode calculateAvgForces(const Vec &f, type::RealVec2D &fAvg);
    
    
    /**
     * \brief update the indices of backgrounf Pressure cells for all body.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode updateMeshIdx();
    

protected:


    /**
     * \brief constructor of using a type::Mesh and a YAML node.
     *
     * \param mesh [in] an instance of type::Mesh for background mesh.
     * \param node [in] a YAML node specifying information of all bodies.
     */
    PetscErrorCode init(const type::Mesh &mesh, const YAML::Node &node);
    

    /** \brief reference to backgrounf mesh. */
    type::Mesh          mesh;


    /** \brief reference to the MPI communicator. */
    MPI_Comm            comm;


    /** \brief the total number of processes. */
    PetscMPIInt         mpiSize;


    /** \brief the rank of this process. */
    PetscMPIInt         mpiRank;


    /** \brief number of local packed variables of all processes. */
    type::IntVec1D      nLclAllProcs;


    /** \brief offsets of packed variables of all processes. */
    type::IntVec1D      offsetsAllProcs;


    /**
     * \brief create DMComposite.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createDmPack();


    /**
     * \brief create a string for printing information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();
};

} // end of namespace body


namespace type
{
    /** \brief definition of type::BodyPack. */
    typedef std::shared_ptr<body::BodyPackBase> BodyPack;
} // end of namespace type


namespace body
{
    /**
     * \brief factory for creating type::BodyPack
     *
     * \param mesh [in] background mesh.
     * \param node [in] YMAL::Node.
     * \param bodies [in] output type::BodyPack instance.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createBodyPack(const type::Mesh &mesh, 
            const YAML::Node &node, type::BodyPack &bodies);
} // end of namespace body
} // end of namespace petibm

/***************************************************************************//**
 * \file singlebody.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SingleBody`.
 */


# pragma once

// STL
# include <memory>

// PETSc
# include <petscsys.h>
# include <petscdm.h>

// PetIBM
# include <petibm/type.h>
# include <petibm/mesh.h>


namespace petibm
{
namespace body
{

// TODO: the way we store coords and meshIdx is different from CartesianMesh. 
//       Will this cause confusion?
/** \brief the abstract class for a single body. */
class SingleBodyBase
{
    friend class BodyPackBase;
    
public:

    /** \brief dimension. */
    PetscInt            dim;

    /** \brief the name of this body. */
    std::string         name;
    
    /** \brief the path/name of the mesh file. */
    std::string         file;

    /** \brief the total number of Lagrangian points. */
    PetscInt            nPts;

    /** \brief coordinates of ALL Lagrangian points. */
    type::RealVec2D     coords;

    /** \brief number of Lagrangain points owned locally. */
    PetscInt            nLclPts;

    /** \brief indices of the background pressure cells that own LOCAL 
     *         Lagrangian points.*/
    type::IntVec2D      meshIdx;

    /** \brief the beginning index of local Lagrangian points in all points. */
    PetscInt            bgPt;

    /** \brief the ending index of local Lagrangian points in all points. */
    PetscInt            edPt;

    /** \brief the underlying parallel 1D DMDA object. */
    DM                  da;

    /** \brief a string for printing information. */
    std::string         info;


    /** \brief the default constructor. */
    SingleBodyBase() = default;


    /**
     * \brief constructor using CartesainMesh and input file.
     *
     * \param mesh [in] an instance of type::Mesh.
     * \param name [in] the name of this body.
     * \param file [in] the ASCII file containing coordinates of Lagrangian points.
     */
    SingleBodyBase(const type::Mesh &mesh,
            const std::string &name, const std::string &file);


    /** \brief the default destructor. */
    virtual ~SingleBodyBase();


    /**
     * \brief manually destroy data.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode destroy();
    
    
    /**
     * \brief print information of this body to standard output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;


    /**
     * \brief find which process owns the Lagrangian point with index i.
     *
     * \param i [in] index of target Lagrangian point of this body.
     * \param p [out] returned process id.
     *
     * Note: all degree of freedom of a Lagrangian point are on the same 
     * process, so we don't need to know which degree of freedom users is asking.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode findProc(const PetscInt &i, PetscMPIInt &p) const = 0;


    /**
     * \brief find the global index in un-packed DM of a specified DoF of a 
     *        Lagrangian point.
     *
     * \param i [in] index of target Lagrangian point of this body.
     * \param dof [in] target degree of freedom.
     * \param idx [out] returned index.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getGlobalIndex(const PetscInt &i, 
            const PetscInt &dof, PetscInt &idx) const = 0;


    /**
     * \brief find the global index in un-packed DM of specified DoF of a 
     *        Lagrangian point.
     *
     * \param s [in] MatStencil of target point and DoF.
     * \param idx [out] returned index.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getGlobalIndex(
            const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief calculate the averaged force of this body.
     *
     * \param f [in] Vec of forces on each Lagrangian point of this body.
     * \param fAvg [out] return averaged force with length equal to dimension.
     *
     * Note: fAvg should have correct size. This function won't check if fAvg
     * has been allocated correctly.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode calculateAvgForces(
            const Vec &f, type::RealVec1D &fAvg) const = 0;
    
    
    /**
     * \brief update the indices of backgrounf Pressure cells.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode updateMeshIdx() = 0;
    

protected:


    /** \brief MPI communicator. */
    MPI_Comm            comm;

    /** \brief the total number of processes. */
    PetscMPIInt         mpiSize;

    /** \brief the rank of this process. */
    PetscMPIInt         mpiRank;
    

    /** \brief reference to backgrounf mesh. */
    type::Mesh          mesh;


    /** \brief number of varaibles (nLcLPts x dim) on each process. */
    type::IntVec1D      nLclAllProcs;


    /** \brief offset on each process. */
    type::IntVec1D      offsetsAllProcs;
};

} // end of namespace body


namespace type
{
    /** \brief definition of type::SingleBody. */
    typedef std::shared_ptr<body::SingleBodyBase> SingleBody;
} // end of namespace type


namespace body
{
    /**
     * \brief factory for creating a SingleBody object.
     *
     * \param mesh [in] the background Eulerian mesh.
     * \param type [in] the type of mesh file (currently only accept "point").
     * \param name [in] the name of this body.
     * \param file [in] the path to the mesh file.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createSingleBody(
            const type::Mesh &mesh, const std::string &type, 
            const std::string &name, const std::string &file,
            type::SingleBody &body);
} // end of namespace body
} // end of namespace petibm

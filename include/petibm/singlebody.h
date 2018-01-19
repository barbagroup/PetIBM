/**
 * \file singlebody.h
 * \brief body::SingleBodyBase, type::SingleBody factory function.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
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
/**
 * \brief Base (abstract) class for a single body.
 * \see bodyModule, petibm::type::SingleBody, petibm::body::createSingleBody
 * \ingroup bodyModule 
 * 
 * There is currently only one implementation for this abstract class: body::SingleBodyPoints.
 */
class SingleBodyBase
{
    friend class BodyPackBase;
    
public:

    /** \brief Dimension. */
    PetscInt            dim;

    /** \brief The name of this body. */
    std::string         name;
    
    /** \brief The path/name of the mesh file. */
    std::string         file;

    /** \brief The total number of Lagrangian points. */
    PetscInt            nPts;

    /** \brief Coordinates of ALL Lagrangian points. */
    type::RealVec2D     coords;

    /** \brief Number of Lagrangian points owned locally. */
    PetscInt            nLclPts;

    /**
     * \brief Indices of the background pressure cells that own LOCAL 
     *        Lagrangian points.
     */
    type::IntVec2D      meshIdx;

    /** \brief The beginning index of local Lagrangian points in all points. */
    PetscInt            bgPt;

    /** \brief The ending index of local Lagrangian points in all points. */
    PetscInt            edPt;

    /** \brief The underlying parallel 1D DMDA object. */
    DM                  da;

    /** \brief A string for printing information. */
    std::string         info;


    /**
     * \brief Constructor using CartesainMesh and input file.
     *
     * \param mesh [in] an instance of type::Mesh.
     * \param name [in] the name of this body.
     * \param file [in] the input file containing necessary information for this body.
     * 
     * The information carried by `file` will be different according to different
     * implementations. For example, the `file` of the current only one 
     * implementation, body::SingleBodyPoints, is a file containing the 
     * coordinates of all Lagrangian points for this body.
     */
    SingleBodyBase(const type::Mesh &mesh,
            const std::string &name, const std::string &file);


    /** \brief The default constructor. */
    SingleBodyBase() = default;


    /** \brief Default destructor. */
    virtual ~SingleBodyBase();


    /**
     * \brief Manually destroy data.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode destroy();
    
    
    /**
     * \brief Print information of this body to standard output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;


    /**
     * \brief Find which process owns the Lagrangian point with index i.
     *
     * \param i [in] index of target Lagrangian point of this body.
     * \param p [out] returned process id.
     *
     * Note: all degree of freedoms of a Lagrangian point are on the same 
     * process, so we don't need to know which degree of freedom users is asking.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode findProc(const PetscInt &i, PetscMPIInt &p) const = 0;


    /**
     * \brief Find the global index in unpacked DM of a specified DoF of a 
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
     * \brief Find the global index in unpacked DM of specified DoF of a 
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
     * \brief Calculate the averaged force of this body.
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
     * \brief Update the indices of background Pressure cells.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode updateMeshIdx() = 0;
    

protected:


    /** \brief MPI communicator. */
    MPI_Comm            comm;

    /** \brief The total number of processes. */
    PetscMPIInt         mpiSize;

    /** \brief The rank of this process. */
    PetscMPIInt         mpiRank;
    

    /** \brief Reference to background mesh. */
    type::Mesh          mesh;


    /** \brief Number of variables (nLcLPts x dim) on each process. */
    type::IntVec1D      nLclAllProcs;


    /** \brief Offset on each process. */
    type::IntVec1D      offsetsAllProcs;
}; // SingleBodyBase

} // end of namespace body


namespace type
{
    /**
     * \brief Definition of type::SingleBody.
     * \see bodyModule, petibm::body::SingleBodyBase, petibm::body::createSingleBody
     * \ingroup bodyModule
     */
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
     * \param body [out] returned type::SingleBody instance.
     *
     * \return PetscErrorCode.
     * \see bodyModule, petibm::type::SingleBody
     * \ingroup bodyModule
     */
    PetscErrorCode createSingleBody(
            const type::Mesh &mesh, const std::string &type, 
            const std::string &name, const std::string &file,
            type::SingleBody &body);
} // end of namespace body
} // end of namespace petibm

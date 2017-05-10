/***************************************************************************//**
 * \file SingleBody.h
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
# include "types.h"
# include "CartesianMesh.h"


// TODO: the way we store coords and meshIdx is different from CartesianMesh. 
//       Will this cause confusion?
/** \brief class for a single body. */
class SingleBody
{
    friend class BodyPack;

public:

    /** \brief dimension. */
    PetscInt            dim;

    /** \brief the name of this body. */
    std::string         name;

    /** \brief the total number of Lagrangian points. */
    PetscInt            nPts;

    /** \brief coordinates of ALL Lagrangian points. */
    types::RealVec2D    coords;

    /** \brief number of Lagrangain points owned locally. */
    PetscInt            nLclPts;

    /** \brief indices of the cells in background pressure mesh that own LOCAL 
     *         Lagrangian points.*/
    types::IntVec2D     meshIdx;

    /** \brief the beginning index of local Lagrangian points in all points. */
    PetscInt            bgPt;

    /** \brief the ending index of local Lagrangian points in all points. */
    PetscInt            edPt;

    /** \brief the underlying parallel 1D DMDA onject. */
    DM                  da;

    /** \brief a string for printing information. */
    std::string         info;


    /** \brief the default constructor. */
    SingleBody();


    /**
     * \brief constructor using CartesainMesh and input file.
     *
     * \param mesh an instance of CartesianMesh.
     * \param file the ASCII file containing coordinates of Lagrangian points.
     * \param name the name of this body.
     */
    SingleBody(const CartesianMesh &mesh, 
            const std::string &file, const std::string &name);


    /**
     * \brief constructor using CartesainMesh and a coordinate vector.
     *
     * \param mesh an instance of CartesianMesh.
     * \param coords a 3 by nPts vector of doubles.
     * \param name the name of this body.
     *
     * Note: Even with 2D problem, this function still requires coordinates in 
     * z-direction. They are simply zeros in this case.
     */
    SingleBody(const CartesianMesh &mesh, 
            const types::RealVec2D &coords, const std::string &name);


    /** \brief the default destructor. */
    ~SingleBody();


    /**
     * \brief initialize the instance using CartesainMesh and input file.
     *
     * \param mesh an instance of CartesianMesh.
     * \param file the ASCII file containing coordinates of Lagrangian points.
     * \param name the name of this body.
     */
    PetscErrorCode init(const CartesianMesh &mesh, 
            const std::string &file, const std::string &name);


    /**
     * \brief initialize the instance using CartesainMesh and a coordinate vector.
     *
     * \param mesh an instance of CartesianMesh.
     * \param coords a 3 by nPts vector of doubles.
     * \param name the name of this body.
     *
     * Note: Even with 2D problem, this function still requires coordinates in 
     * z-direction. They are simply zeros in this case.
     */
    PetscErrorCode init(const CartesianMesh &mesh, 
            const types::RealVec2D &coords, const std::string &name);


    /**
     * \brief print the information of this body.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;


    /**
     * \brief find which process owns the Lagrangian point with index i.
     *
     * \param i index of target Lagrangian point of this body.
     * \param p return process id.
     *
     * Note: all degree of freedom of a Lagrangian point are on the same 
     * process, so we don't need to know which degree of freedom users is asking.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode findProc(const PetscInt &i, PetscMPIInt &p) const;


    /**
     * \brief find the global index in un-packed DM of specified DoF of a 
     *        Lagrangian point.
     *
     * \param i index of target Lagrangian point of this body.
     * \param dof target degree of freedom.
     * \param idx returned index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getGlobalIndex(
            const PetscInt &i, const PetscInt &dof, PetscInt &idx) const;


    /**
     * \brief find the global index in un-packed DM of specified DoF of a 
     *        Lagrangian point.
     *
     * \param s MatStencil of target point and DoF.
     * \param idx returned index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getGlobalIndex(const MatStencil &s, PetscInt &idx) const;


    /**
     * \brief calculate the averaged force of this body.
     *
     * \param f Vec of forces on each Lagrangian point of this body.
     * \param fAvg return averaged force.
     *
     * Note: fAvg should have correct size. This function won't check if fAvg
     * has been allocated correctly.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode calculateAvgForces(const Vec &f, types::RealVec1D &fAvg);

protected:

    /** \brief reference to backgrounf CartesianMesh. */
    std::shared_ptr<const CartesianMesh>    mesh;


    /** \brief reference to the MPI communicator. */
    std::shared_ptr<const MPI_Comm>         comm;


    /** \brief the total number of processes. */
    PetscMPIInt                             mpiSize;


    /** \brief the rank of this process. */
    PetscMPIInt                             mpiRank;


    /** \brief number of varaibles (nLcLPts x dim) on each process. */
    types::IntVec1D                         nLclAllProcs;


    /** \brief offset on each process. */
    types::IntVec1D                         offsetsAllProcs;


    /**
     * \brief part of initialization that initialize mesh and MPI information.
     *
     * \param mesh an instance of CartesianMesh.
     * \param name the name of this body.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode preInit(const CartesianMesh &mesh, const std::string &name);


    /**
     * \brief part of initialization that initialize remaining members.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode postInit();


    /**
     * \brief read coordinates and number of points from an ASCII file.
     *
     * \param file the ASCII file containing coordinates of Lagrangian points.
     *
     * \return 
     */
    PetscErrorCode readFromFile(const std::string &file);


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

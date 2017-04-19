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


/** \brief class for a single body. */
class SingleBody
{
public:

    /** \brief dimension. */
    PetscInt            dim;

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
    PetscInt            bg;

    /** \brief the ending index of local Lagrangian points in all points. */
    PetscInt            ed;

    /** \brief the underlying parallel 1D DMDA onject. */
    DM                  da;


    /** \brief the default constructor. */
    SingleBody();


    /**
     * \brief constructor using CartesainMesh and input file.
     *
     * \param mesh an instance of CartesianMesh.
     * \param file the ASCII file containing coordinates of Lagrangian points.
     */
    SingleBody(const CartesianMesh &mesh, const std::string &file);


    /**
     * \brief constructor using CartesainMesh and a coordinate vector.
     *
     * \param mesh an instance of CartesianMesh.
     * \param coords a 3 by nPts vector of doubles.
     *
     * Note: Even with 2D problem, this function still requires coordinates in 
     * z-direction. They are simply zeros in this case.
     */
    SingleBody(const CartesianMesh &mesh, const types::RealVec2D &coords);


    /** \brief the default destructor. */
    ~SingleBody();


    /**
     * \brief initialize the instance using CartesainMesh and input file.
     *
     * \param mesh an instance of CartesianMesh.
     * \param file the ASCII file containing coordinates of Lagrangian points.
     */
    PetscErrorCode init(
            const CartesianMesh &mesh, const std::string &file);


    /**
     * \brief initialize the instance using CartesainMesh and a coordinate vector.
     *
     * \param mesh an instance of CartesianMesh.
     * \param coords a 3 by nPts vector of doubles.
     *
     * Note: Even with 2D problem, this function still requires coordinates in 
     * z-direction. They are simply zeros in this case.
     */
    PetscErrorCode init(
            const CartesianMesh &mesh, const types::RealVec2D &coords);

protected:

    /** \brief reference to backgrounf CartesianMesh. */
    std::shared_ptr<const CartesianMesh>    mesh;


    /** \brief reference to the MPI communicator. */
    std::shared_ptr<const MPI_Comm>         comm;


    /** \brief the total number of processes. */
    PetscMPIInt                             mpiSize;


    /** \brief the rank of this process. */
    PetscMPIInt                             mpiRank;


    /**
     * \brief part of initialization that initialize mesh and MPI information.
     *
     * \param mesh an instance of CartesianMesh.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode preInit(const CartesianMesh &mesh);


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
};

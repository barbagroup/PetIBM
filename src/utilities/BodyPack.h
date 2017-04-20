/***************************************************************************//**
 * \file BodyPack.h
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
# include "SingleBody.h"


/** \brief class for a pack of multiple bodies. */
class BodyPack
{
public:

    /** \brief dimension. */
    PetscInt        dim;

    /** \brief number of bodies in this pack. */
    PetscInt        nBodies;

    /** \brief a vector of SingleBody instances. */
    std::vector<SingleBody>     bodies;

    /** \brief a DMComposite of DMs of all `SingleBody`s. */
    DM              dmPack;

    /** \brief a string for printing information. */
    std::string     info;


    /** \brief default constructor. */
    BodyPack();


    /**
     * \brief constructor of using CartesianMesh and a YAML node.
     *
     * \param mesh an instance of CartesianMesh for background mesh.
     * \param node a YAML node specifying information of all bodies.
     */
    BodyPack(const CartesianMesh &mesh, const YAML::Node &node);


    /** \brief default destructor. */
    ~BodyPack();


    /**
     * \brief initialization using CartesianMesh and a YAML node.
     *
     * \param mesh an instance of CartesianMesh for background mesh.
     * \param node a YAML node specifying information of all bodies.
     */
    PetscErrorCode init(const CartesianMesh &mesh, const YAML::Node &node);


    /**
     * \brief print information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo();

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

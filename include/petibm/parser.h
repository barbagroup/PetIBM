/**
 * \file parser.h
 * \brief Prototypes of parser functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

# pragma once

// here goes C++ STL
# include <string>

// here goes PETSc
# include <petscsys.h>

// here goes yaml-cpp
# include <yaml-cpp/yaml.h>

// here goes our owen headers
# include <petibm/type.h>


namespace petibm
{
/**
 * \brief A collection of YAML node parsers.
 * \ingroup miscModule
 */
namespace parser
{
    /**
     * \brief Get settings from command line arguments and read YAML files.
     * \param node [out] a YAML node containing all settings.
     * \return PetscErrorCode.
     * \ingroup miscModule
     *
     * The function will look for the following command-line arguments 
     *
     * -# `-directory`: the working directory. Default is the current 
     *    directory if not found.
     * -# `-config`: location of config.yaml. If not provided, the default is
     *    [working directory]/config.yaml.
     * -# `-mesh`: location of mesh.yaml. This provides a way to overwrite the
     *    mesh section inside config.yaml.
     * -# `-flow`: location of flow.yaml. This provides a way to overwrite the
     *    flow section inside config.yaml.
     * -# `-parameters`: location of parameters.yaml. This provides a way to 
     *    overwrite the parameters section inside config.yaml.
     * -# `-bodies`: location of bodies.yaml. This provides a way to overwrite
     *    the bodies section inside config.yaml.
     *
     * If users provide non-empty YAML node as input, the data inside the node
     * will be discarded.
     */
    PetscErrorCode getSettings(YAML::Node &node);

    /** 
     * \brief Parse a YAML node of cartesianMesh.
     * \param meshNode [in] the YAML node holding mesh settings.
     * \param dim [out] returned dimension
     * \param bg [out] a vector of the starting boundary in each direction.
     * \param ed [out] a vector of the ending boundary in each direction.
     * \param nTotal [out] a vector for number of pressure cell in each direction.
     * \param dL [out] a nested vector for pressure cell sizes in each direction.
     * \return PetscErrorCode
     * \ingroup miscModule
     * 
     * This function will retrieve information from and do necessary calculations
     * based on the data provided under the key `mesh` in the YAML node.
     */
    PetscErrorCode parseMesh(
            const YAML::Node &meshNode, PetscInt &dim, type::RealVec1D &bg, 
            type::RealVec1D &ed, type::IntVec1D &nTotal, type::RealVec2D &dL);

    /** 
     * \brief Parse the info of only one direction from YAML node.
     * \param axis [in] the YAML node.
     * \param dir [out] returned direction label.
     * \param bg [out] returned starting boundary in this direction.
     * \param ed [out] returned ending boundary in this direction.
     * \param nTotal [out] returned total number of pressure cells in this direction.
     * \param dL [out] returned 1D vector for the size of each pressure cell in this direction.
     * \return  PetscErrorCode
     * \ingroup miscModule
     * 
     * This function only parse information and do calculation for only one axis.
     */
    PetscErrorCode parseOneAxis(
            const YAML::Node &axis, PetscInt &dir, PetscReal &bg, 
            PetscReal &ed, PetscInt &nTotal, type::RealVec1D &dL);

    /**
     * \brief Parse all sub-domains in a direction.
     * \param subs [in] the YAML node
     * \param bg [in] an input value providing starting of this direction.
     * \param nTotal [out] returned total number of pressure cells in this direction.
     * \param ed [out] returned ending of this direction.
     * \param dL [out] returned pressure cell sizes in this direction.
     * \return PetscErrorCode
     * \ingroup miscModule
     * 
     * This function parse the `subDomains` section of a direction.
     */
    PetscErrorCode parseSubDomains(
            const YAML::Node &subs, const PetscReal bg,
            PetscInt &nTotal, PetscReal &ed, type::RealVec1D &dL);

    /**
     * \brief Parse only one sub-domain
     * \param sub [in] the YAML node
     * \param bg [in] an input providing the starting of this sub-domain.
     * \param n [out] returned number of pressure cells in this sub-domain.
     * \param ed [out] returned ending of this sub-domain.
     * \param dL [out] returned 1D vector for sizes of pressure cells in this sub-domain.
     * \return PetscErrorCode
     * \ingroup miscModule
     */
    PetscErrorCode parseOneSubDomain(
            const YAML::Node &sub, const PetscReal bg,
            PetscInt &n, PetscReal &ed, type::RealVec1D &dL);
    
    
    /**
     * \brief Parse boundary conditions from a YAML node.
     * \param node [in] a YAML node containing boundary conditions.
     * \param bcTypes [out] BC types of different fields on different boundaries.
     * \param bcValues [out] BC values of different fields on different boundaries.
     * \return PetscErrorCode.
     * \ingroup miscModule
     * 
     * This function will look into the key `boundaryConditions` under the key 
     * `flow` in the provided YAML node.
     */
    PetscErrorCode parseBCs(const YAML::Node &node, 
            type::IntVec2D &bcTypes, type::RealVec2D &bcValues);
    
    
    /**
     * \brief Parse initial conditions from a YAML node.
     * \param node [in] a YAML node.
     * \param icValues [out] IC values of different fields.
     * \return PetscErrorCode.
     * \ingroup miscModule
     * 
     * This function will look into the key `initialVelocity` under the key 
     * `flow` in the provided YAML node.
     */
    PetscErrorCode parseICs(const YAML::Node &node, type::RealVec1D &icValues);

} // end of namespace parser
} // end of namespace petibm

/***************************************************************************//**
 * \file parser.h
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Declaration of functions under namespace `parser`.
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
/** \brief YAML node parsers for PetIBM components. */
namespace parser
{
    /**
     * \brief get settings from command line arguments and read YAML files.
     *
     * The function will look for the following command-line arguments 
     *
     *     1. -directory: the working directory. Default is the current 
     *        directory if not found.
     *     2. -config: location of config.yaml. If not provided, the default is
     *        [working directory]/config.yaml.
     *     3. -mesh: location of mesh.yaml. This provides a way to overwrite the
     *        mesh section indide config.yaml.
     *     4. -flow: location of flow.yaml. This provides a way to overwrite the
     *        flow section indide config.yaml.
     *     5. -parameters: location of parameters.yaml. This provides a way to 
     *        overwrite the parameters section indide config.yaml.
     *     6. -bodies: location of bodies.yaml. This provides a way to overwrite
     *        the bodies section indide config.yaml.
     *
     * \param node [out] a YAML node containing all settings.
     *
     * If users provide non-empty YAML node as input, the data inside the node
     * will be discarded.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getSettings(YAML::Node &node);

    /** \brief parse a YAML node of cartesianMesh.
     *
     * \param meshNode the YAML node holding info for creating `cartesianMesh`.
     * \param dim returned dimension
     * \param bg a vector of the starting boundary in each direction.
     * \param ed a vector of the ending boundary in each direction.
     * \param nTotal a vector for number of pressure cell in each direction.
     * \param dL a nested vector for pressure cell sizes in each direction.
     *
     * \return PetscErrorCode
     */
    PetscErrorCode parseMesh(
            const YAML::Node &meshNode, PetscInt &dim, type::RealVec1D &bg, 
            type::RealVec1D &ed, type::IntVec1D &nTotal, type::RealVec2D &dL);

    /** \brief parse the info of only one direction from YAML node.
     *
     * \param axis the YAML node.
     * \param dir returned direction label.
     * \param bg returned starting boundary in this direction.
     * \param ed returned ending boundary in this direction.
     * \param nTotal returned total number of pressure cells in this direction.
     * \param dL returned 1D vector for the size of each pressure cell in this direction.
     *
     * \return  PetscErrorCode
     */
    PetscErrorCode parseOneAxis(
            const YAML::Node &axis, PetscInt &dir, PetscReal &bg, 
            PetscReal &ed, PetscInt &nTotal, type::RealVec1D &dL);

    /** \brief parse all subdomains in a direction.
     *
     * \param subs the YAML node
     * \param bg an input value providing starting of this direction.
     * \param nTotal returned total number of pressure cells in this direction.
     * \param ed returned ending of this direction.
     * \param dL returned pressure cell sizes in this direction.
     *
     * \return PetscErrorCode
     */
    PetscErrorCode parseSubDomains(
            const YAML::Node &subs, const PetscReal bg,
            PetscInt &nTotal, PetscReal &ed, type::RealVec1D &dL);

    /** \brief parse only one subdomain
     *
     * \param sub the YAML node
     * \param bg an input providing the starting of this subdomain.
     * \param n returned number of pressure cells in this subdomain.
     * \param ed returned ending of this subdomain.
     * \param dL returned 1D vector for sizes of pressure cells in this subdomain.
     *
     * \return PetscErrorCode
     */
    PetscErrorCode parseOneSubDomain(
            const YAML::Node &sub, const PetscReal bg,
            PetscInt &n, PetscReal &ed, type::RealVec1D &dL);
    
    
    /**
     * \brief parse boundary conditions from a YAML node.
     *
     * \param node [in] a YAML node.
     * \param bcTypes [out] BC types of different fields on different boundaries.
     * \param bcValues [out] BC values of different fields on different boundaries.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode parseBCs(const YAML::Node &node, 
            type::IntVec2D &bcTypes, type::RealVec2D &bcValues);
    
    
    /**
     * \brief parse initial conditions from a YAML node.
     *
     * \param node [in] a YAML node.
     * \param icValues [out] IC values of different fileds.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode parseICs(const YAML::Node &node, type::RealVec1D &icValues);

} // end of namespace parser
} // end of namespace petibm

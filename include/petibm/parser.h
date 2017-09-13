/***************************************************************************//**
 * \file parser.h
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Declaration of functions under namespace `parser`.
 */

# pragma once

// here goes C++ STL
# include <string>
# include <experimental/filesystem>

// here goes PETSc
# include <petscsys.h>

// here goes yaml-cpp
# include <yaml-cpp/yaml.h>

// here goes our owen headers
# include "types.h"


namespace petibm
{
namespace utilities
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
	 *     	  [working directory]/config.yaml.
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

    /** \brief parse a YAML node and obtained info for class `FlowDescription`.
     *
     * \param flowNode a YAML node.
     * \param dim dimension will be returned here.
     * \param nu kinematic viscosity will be returned here.
     * \param customIC a flog indicating if using customized IC.
     * \param IC a vector holding constant IC value for each component will be returned here.
     * \param pertb perturbation info will be returned here.
     * \param nBC number of total boundaries will be returned here
     * \param BCInfo a `BCInfoHolder` will be returned here 
     *
     * \return PetscErrorCode
     */
    PetscErrorCode parseFlowDescription(
            const YAML::Node &flowNode, PetscInt &dim, PetscReal &nu, 
            PetscBool &customIC, types::RealVec1D &IC, types::Perturbation &pertb, 
            PetscInt &nBC, types::BCInfoHolder &BCInfo);

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
            const YAML::Node &meshNode, PetscInt &dim, types::RealVec1D &bg, 
            types::RealVec1D &ed, types::IntVec1D &nTotal, types::RealVec2D &dL);

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
            PetscReal &ed, PetscInt &nTotal, types::RealVec1D &dL);

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
            PetscInt &nTotal, PetscReal &ed, types::RealVec1D &dL);

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
            PetscInt &n, PetscReal &ed, types::RealVec1D &dL);

} // end of namespace parser
} // end of namespace utilities
} // end of namespace petibm

/**
 * \file parser.h
 * \brief Prototypes of parser functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <string>

#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/type.h>

namespace petibm
{
/**
 * \brief A collection of YAML node parsers.
 *
 * \ingroup miscModule
 */
namespace parser
{
/**
 * \brief Get settings from command line arguments and read YAML files.
 *
 * \param node [out] YAML configuration node.
 *
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
 *
 * \param meshNode [in] YAML configuration node.
 * \param dim [out] Number of dimensions
 * \param bg [out] Vector with starting coordinate in each direction.
 * \param ed [out] Vector with ending coordinate in each direction.
 * \param nTotal [out] Number of cells in each direction.
 * \param dL [out] Nested vector with grid cell sizes in each direction.
 *
 * \ingroup miscModule
 *
 * This function will retrieve information from and do necessary calculations
 * based on the data provided under the key `mesh` in the YAML node.
 */
PetscErrorCode parseMesh(const YAML::Node &meshNode, PetscInt &dim,
                         type::RealVec1D &bg, type::RealVec1D &ed,
                         type::IntVec1D &nTotal, type::RealVec2D &dL);

/**
 * \brief Parse the info of only one direction from YAML node.
 *
 * \param axis [in] YAML configuration node.
 * \param dir [out] Label of the direction.
 * \param bg [out] Starting coordinate.
 * \param ed [out] Ending coordinate.
 * \param nTotal [out] Number of cells along the direction.
 * \param dL [out] 1D vector with cell sizes along the direction.
 *
 * \ingroup miscModule
 *
 * This function only parse information and do calculation for only one axis.
 */
PetscErrorCode parseOneAxis(const YAML::Node &axis, PetscInt &dir,
                            PetscReal &bg, PetscReal &ed, PetscInt &nTotal,
                            type::RealVec1D &dL);

/**
 * \brief Parse all sub-domains in a direction.
 *
 * \param subs [in] YAML configuration node.
 * \param bg [in] Starting coordinate.
 * \param nTotal [out] Number of cell along the direction.
 * \param ed [out] Ending coordinate.
 * \param dL [out] Grid cell sizes along the direction.
 *
 * \ingroup miscModule
 *
 * This function parse the `subDomains` section of a direction.
 */
PetscErrorCode parseSubDomains(const YAML::Node &subs, const PetscReal bg,
                               PetscInt &nTotal, PetscReal &ed,
                               type::RealVec1D &dL);

/**
 * \brief Parse only one sub-domain.
 *
 * \param sub [in] YAML configuration node.
 * \param bg [in] Starting coordinate.
 * \param n [out] Number of cells along the direction of the sub-domain.
 * \param ed [out] Ending coordinate.
 * \param dL [out] Grid cell sizes along the direction of the sub-domain.
 *
 * \ingroup miscModule
 */
PetscErrorCode parseOneSubDomain(const YAML::Node &sub, const PetscReal bg,
                                 PetscInt &n, PetscReal &ed,
                                 type::RealVec1D &dL);

/**
 * \brief Parse boundary conditions from a YAML node.
 *
 * \param node [in] YAML configuration node.
 * \param bcTypes [out] BC types of different fields on different boundaries.
 * \param bcValues [out] BC values of different fields on different boundaries.
 *
 * \ingroup miscModule
 *
 * This function will look into the key `boundaryConditions` under the key
 * `flow` in the provided YAML node.
 */
PetscErrorCode parseBCs(const YAML::Node &node, type::IntVec2D &bcTypes,
                        type::RealVec2D &bcValues);

/**
 * \brief Parse initial conditions from a YAML node.
 *
 * \param node [in] YAML configuration node.
 * \param icValues [out] IC values of different fields.
 *
 * \ingroup miscModule
 *
 * This function will look into the key `initialVelocity` under the key
 * `flow` in the provided YAML node.
 */
PetscErrorCode parseICs(const YAML::Node &node, type::RealVec1D &icValues);

}  // end of namespace parser

}  // end of namespace petibm

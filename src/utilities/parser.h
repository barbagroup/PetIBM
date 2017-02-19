/***************************************************************************//**
 * \file parser.h
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Declaration of functions under namespace `parser` and `YAML`.
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
# include <types.h>


/** \brief YAML node parsers for PetIBM components. */
namespace parser
{
    /** \brief parse the config.yaml in a given case directory
     *
     * \param dir the directory of the targeting case
     * \param node returned YAML node
     *
     * \return PetscErrorCode
     */
    PetscErrorCode parseYAMLConfigFile(
            const std::string &dir, YAML::Node &node);

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
}


/** \brief YAML node converters used by yamp-cpp */
namespace YAML
{
    /** \brief converter for `types::Dir` */
    template <>
    struct convert<types::Dir>
    {
        static Node encode(const types::Dir &dir);
        static bool decode(const Node &node, types::Dir &dir);
    };

    /** \brief converter for `tpyes::VelocityComponent` */
    template <>
    struct convert<types::VelocityComponent>
    {
        static Node encode(const types::VelocityComponent &vc);
        static bool decode(const Node &node, types::VelocityComponent &vc);
    };

    /** \brief converter for `types::BCType` */
    template <>
    struct convert<types::BCType>
    {
        static Node encode(const types::BCType &bc);
        static bool decode(const Node &node, types::BCType &bc);
    };

    /** \brief converter for `types::BCLoc` */
    template <>
    struct convert<types::BCLoc>
    {
        static Node encode(const types::BCLoc &loc);
        static bool decode(const Node &node, types::BCLoc &loc);
    };

    /** \brief converter for `types::TimeScheme` */
    template <>
    struct convert<types::TimeScheme>
    {
        static Node encode(const types::TimeScheme &ts);
        static bool decode(const Node &node, types::TimeScheme &ts);
    };

    /** \brief converter for `types::IBMethod` */
    template <>
    struct convert<types::IBMethod>
    {
        static Node encode(const types::IBMethod &ibm);
        static bool decode(const Node &node, types::IBMethod &ibm);
    };

    /** \brief converter for `types::StaggeredMode` */
    template <>
    struct convert<types::StaggeredMode>
    {
        static Node encode(const types::StaggeredMode &sm);
        static bool decode(const Node &node, types::StaggeredMode &sm);
    };

    /** \brief converter for types::ExecuteType` */
    template <>
    struct convert<types::ExecuteType>
    {
        static Node encode(const types::ExecuteType &et);
        static bool decode(const Node &node, types::ExecuteType &et);
    };

    /** \brief converter for `types::BCTypeValuePair` */
    template<>
    struct convert<types::BCTypeValuePair>
    {
        static Node encode(const types::BCTypeValuePair &bcInfo);
        static bool decode(const Node &node, types::BCTypeValuePair &bcInfo);
    };

    /** \brief converter for `types::BCInfoHolder` */
    template<>
    struct convert<types::BCInfoHolder>
    {
        static Node encode(const types::BCInfoHolder &bcInfo);
        static bool decode(const Node &node, types::BCInfoHolder &bcInfo);
    };

    /** \brief converter for `types::Perturbation` */
    template<>
    struct convert<types::Perturbation>
    {
        static Node encode(const types::Perturbation &pertb);
        static bool decode(const Node &node, types::Perturbation &pertb);
    };

    /** \brief converter for `PetscBool` */
    template<>
    struct convert<PetscBool>
    {
        static Node encode(const PetscBool &b);
        static bool decode(const Node &node, PetscBool &b);
    };
}

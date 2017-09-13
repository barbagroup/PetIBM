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
# include "types.h"


/** \brief YAML node converters used by yamp-cpp */
namespace YAML
{

using namespace petibm::utilities::types;

    /** \brief converter for `types::Dir` */
    template <>
    struct convert<Dir>
    {
        static Node encode(const Dir &dir);
        static bool decode(const Node &node, Dir &dir);
    };

    /** \brief converter for `tpyes::Field` */
    template <>
    struct convert<Field>
    {
        static Node encode(const Field &vc);
        static bool decode(const Node &node, Field &vc);
    };

    /** \brief converter for `types::BCType` */
    template <>
    struct convert<BCType>
    {
        static Node encode(const BCType &bc);
        static bool decode(const Node &node, BCType &bc);
    };

    /** \brief converter for `types::BCLoc` */
    template <>
    struct convert<BCLoc>
    {
        static Node encode(const BCLoc &loc);
        static bool decode(const Node &node, BCLoc &loc);
    };

    /** \brief converter for `types::TimeScheme` */
    template <>
    struct convert<TimeScheme>
    {
        static Node encode(const TimeScheme &ts);
        static bool decode(const Node &node, TimeScheme &ts);
    };

    /** \brief converter for `types::IBMethod` */
    template <>
    struct convert<IBMethod>
    {
        static Node encode(const IBMethod &ibm);
        static bool decode(const Node &node, IBMethod &ibm);
    };

    /** \brief converter for `types::StaggeredMode` */
    template <>
    struct convert<StaggeredMode>
    {
        static Node encode(const StaggeredMode &sm);
        static bool decode(const Node &node, StaggeredMode &sm);
    };

    /** \brief converter for types::ExecuteType` */
    template <>
    struct convert<ExecuteType>
    {
        static Node encode(const ExecuteType &et);
        static bool decode(const Node &node, ExecuteType &et);
    };

    /** \brief converter for `types::BCTypeValuePair` */
    template<>
    struct convert<BCTypeValuePair>
    {
        static Node encode(const BCTypeValuePair &bcInfo);
        static bool decode(const Node &node, BCTypeValuePair &bcInfo);
    };

    /** \brief converter for `types::BCInfoHolder` */
    template<>
    struct convert<BCInfoHolder>
    {
        static Node encode(const BCInfoHolder &bcInfo);
        static bool decode(const Node &node, BCInfoHolder &bcInfo);
    };

    /** \brief converter for `types::Perturbation` */
    template<>
    struct convert<Perturbation>
    {
        static Node encode(const Perturbation &pertb);
        static bool decode(const Node &node, Perturbation &pertb);
    };

    /** \brief converter for `types::OutputType` */
    template<>
    struct convert<OutputType>
    {
        static Node encode(const OutputType &out);
        static bool decode(const Node &node, OutputType &out);
    };

    /** \brief converter for `PetscBool` */
    template<>
    struct convert<PetscBool>
    {
        static Node encode(const PetscBool &b);
        static bool decode(const Node &node, PetscBool &b);
    };

    /** \brief converter for `types::OutputInfo` */
    template<>
    struct convert<OutputInfo>
    {
        static Node encode(const OutputInfo &output);
        static bool decode(const Node &node, OutputInfo &output);
    };

    /** \brief converter for `types::LinSolverInfo` */
    template<>
    struct convert<LinSolverInfo>
    {
        static Node encode(const LinSolverInfo &solver);
        static bool decode(const Node &node, LinSolverInfo &solver);
    };

    /** \brief converter for `types::SchemeInfo` */
    template<>
    struct convert<SchemeInfo>
    {
        static Node encode(const SchemeInfo &scheme);
        static bool decode(const Node &node, SchemeInfo &scheme);
    };

    /** \brief converter for `types::SteppingInfo` */
    template<>
    struct convert<SteppingInfo>
    {
        static Node encode(const SteppingInfo &stepping);
        static bool decode(const Node &node, SteppingInfo &stepping);
    };

} // end of namespace YAML


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

    /**
     * \brief parse YAML node to get information of SimulationParameters
     *
     * \param param the YAML node.
     * \param output a struct holding output format settings.
     * \param velocity a struct holding velocity solver settings.
     * \param poisson a struct holding velocity solver settings.
     * \param methods a struct holding numerical schemes.
     * \param stepping a struct holding information of time stepping.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode parseSimulationParameters(const YAML::Node &param,
            types::OutputInfo  &output, types::LinSolverInfo &velocity, 
            types::LinSolverInfo &poisson, types::SchemeInfo  &methods, 
            types::SteppingInfo &stepping);

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

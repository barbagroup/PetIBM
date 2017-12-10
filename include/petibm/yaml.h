/***************************************************************************//**
 * \file yaml.h
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Declaration of functions under namespace `YAML`.
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


/** \brief YAML node converters used by yamp-cpp */
namespace YAML
{

using namespace petibm::type;

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

    /** \brief converter for `types::SteppingInfo` */
    template<>
    struct convert<SteppingInfo>
    {
        static Node encode(const SteppingInfo &stepping);
        static bool decode(const Node &node, SteppingInfo &stepping);
    };

} // end of namespace YAML

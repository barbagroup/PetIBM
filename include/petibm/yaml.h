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

    /** \brief converter for `PetscBool` */
    template<>
    struct convert<PetscBool>
    {
        static Node encode(const PetscBool &b);
        static bool decode(const Node &node, PetscBool &b);
    };

} // end of namespace YAML

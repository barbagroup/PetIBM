/**
 * \file yaml.h
 * \brief Prototypes of YAML converters.
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


/**
 * \brief A supplement to YAML-CPP that adds converters of our user-defined types.
 * \ingroup miscModule
 */
namespace YAML
{

    /** \brief converter for `types::Dir` */
    template <>
    struct convert<petibm::type::Dir>
    {
        static Node encode(const petibm::type::Dir &dir);
        static bool decode(const Node &node, petibm::type::Dir &dir);
    };

    /** \brief converter for `types::Field` */
    template <>
    struct convert<petibm::type::Field>
    {
        static Node encode(const petibm::type::Field &vc);
        static bool decode(const Node &node, petibm::type::Field &vc);
    };

    /** \brief converter for `types::BCType` */
    template <>
    struct convert<petibm::type::BCType>
    {
        static Node encode(const petibm::type::BCType &bc);
        static bool decode(const Node &node, petibm::type::BCType &bc);
    };

    /** \brief converter for `types::BCLoc` */
    template <>
    struct convert<petibm::type::BCLoc>
    {
        static Node encode(const petibm::type::BCLoc &loc);
        static bool decode(const Node &node, petibm::type::BCLoc &loc);
    };

    /** \brief converter for `PetscBool` */
    template<>
    struct convert<PetscBool>
    {
        static Node encode(const PetscBool &b);
        static bool decode(const Node &node, PetscBool &b);
    };

} // end of namespace YAML

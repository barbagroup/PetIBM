/**
 * \file yaml.cpp
 * \brief Implementations of YAML converters.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// here goes our own headers
# include "petibm/yaml.h"
# include "petibm/misc.h"


namespace YAML
{

    // for Dir
    Node convert<petibm::type::Dir>::encode(const petibm::type::Dir &dir)
    {
        Node node;
        node = petibm::type::dir2str[dir];
        return node;
    }

    bool convert<petibm::type::Dir>::decode(const Node &node, petibm::type::Dir &dir)
    {
        if (! node.IsDefined()) return false;

        dir = petibm::type::str2dir[node.as<std::string>()];
        return true;
    }

    // for Field
    Node convert<petibm::type::Field>::encode(const petibm::type::Field &vc)
    {
        Node node;
        node = petibm::type::fd2str[vc];
        return node;
    }

    bool convert<petibm::type::Field>::decode(const Node &node, petibm::type::Field &vc)
    {
        if (! node.IsDefined()) return false;

        vc = petibm::type::str2fd[node.as<std::string>()];
        return true;
    }

    // for BCType
    Node convert<petibm::type::BCType>::encode(const petibm::type::BCType &bc)
    {
        Node node;
        node = petibm::type::bt2str[bc];
        return node;
    }

    bool convert<petibm::type::BCType>::decode(const Node &node, petibm::type::BCType &bc)
    {
        if (! node.IsDefined()) return false;
        bc = petibm::type::str2bt[node.as<std::string>()];
        return true;
    }

    // for BCLoc
    Node convert<petibm::type::BCLoc>::encode(const petibm::type::BCLoc &loc)
    {
        Node node;
        node = petibm::type::bl2str[loc];
        return node;
    }

    bool convert<petibm::type::BCLoc>::decode(const Node &node, petibm::type::BCLoc &loc)
    {
        if (! node.IsDefined()) return false;

        loc = petibm::type::str2bl[node.as<std::string>()];
        return true;
    }

    // for PetscBool
    Node convert<PetscBool>::encode(const PetscBool &b)
    {
        YAML::Node  node;
        node = bool(b);
        return node;
    }

    bool convert<PetscBool>::decode(const Node &node, PetscBool &b)
    {
        if (! node.IsDefined()) return false;

        b = PetscBool(node.as<bool>());
        return true;
    }
    
} // end of namespace YAML

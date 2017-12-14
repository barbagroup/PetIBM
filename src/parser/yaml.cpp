/*
 * yaml.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */


// here goes our own headers
# include "petibm/yaml.h"
# include "petibm/misc.h"

/**
 * \brief user-defined YAML converter
 * \todo when converting our data types to YAML map, the order of keys are 
 * not retained. Need som extra work if we want to keep the order.
 */
namespace YAML
{

using namespace petibm::type;

    // for Dir
    Node convert<Dir>::encode(const Dir &dir)
    {
        Node node;
        node = dir2str[dir];
        return node;
    }

    bool convert<Dir>::decode(const Node &node, Dir &dir)
    {
        if (! node.IsDefined()) return false;

        dir = str2dir[node.as<std::string>()];
        return true;
    }

    // for Field
    Node convert<Field>::encode(const Field &vc)
    {
        Node node;
        node = fd2str[vc];
        return node;
    }

    bool convert<Field>::decode(const Node &node, Field &vc)
    {
        if (! node.IsDefined()) return false;

        vc = str2fd[node.as<std::string>()];
        return true;
    }

    // for BCType
    Node convert<BCType>::encode(const BCType &bc)
    {
        Node node;
        node = bt2str[bc];
        return node;
    }

    bool convert<BCType>::decode(const Node &node, BCType &bc)
    {
        if (! node.IsDefined()) return false;
        bc = str2bt[node.as<std::string>()];
        return true;
    }

    // for BCLoc
    Node convert<BCLoc>::encode(const BCLoc &loc)
    {
        Node node;
        node = bl2str[loc];
        return node;
    }

    bool convert<BCLoc>::decode(const Node &node, BCLoc &loc)
    {
        if (! node.IsDefined()) return false;

        loc = str2bl[node.as<std::string>()];
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

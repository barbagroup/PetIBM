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

using namespace petibm::utilities::types;

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

    // for TimeScheme
    Node convert<TimeScheme>::encode(const TimeScheme &ts)
    {
        Node node;
        node = ts2str[ts];
        return node;
    }

    bool convert<TimeScheme>::decode(const Node &node, TimeScheme &ts)
    {
        if (! node.IsDefined()) return false;
        ts = str2ts[node.as<std::string>()];
        return true;
    }

    // for IBMethod
    Node convert<IBMethod>::encode(const IBMethod &ibm)
    {
        Node node;
        node = ibm2str[ibm];
        return node;
    }

    bool convert<IBMethod>::decode(const Node &node, IBMethod &ibm)
    {
        if (! node.IsDefined()) return false;

        ibm = str2ibm[node.as<std::string>()];
        return true;
    }

    // for StaggeredMode
    Node convert<StaggeredMode>::encode(const StaggeredMode &sm)
    {
        Node node;
        node = sm2str[sm];
        return node;
    }

    bool convert<StaggeredMode>::decode(const Node &node, StaggeredMode &sm)
    {
        if (! node.IsDefined()) return false;

        sm = str2sm[node.as<std::string>()];
        return true;
    }

    // for ExecutiveType
    Node convert<ExecuteType>::encode(const ExecuteType &et)
    {
        Node node;
        node = et2str[et];
        return node;
    }

    bool convert<ExecuteType>::decode(const Node &node, ExecuteType &et)
    {
        if (! node.IsDefined()) return false;
        
        et = str2et[node.as<std::string>()];
        return true;
    }

    // for BCTypeValuePair
    Node convert<BCTypeValuePair>::encode(const BCTypeValuePair &bcInfo)
    {
        Emitter     out;
        out << Flow;
        out << BeginSeq << bt2str[bcInfo.type] << bcInfo.value << EndSeq;
        return Load(out.c_str());
    }

    bool convert<BCTypeValuePair>::decode(const Node &node,
                                          BCTypeValuePair &bcInfo)
    {
        if((! node[0].IsDefined()) || (! node[1].IsDefined())) return false;

        bcInfo.type = node[0].as<BCType>();
        bcInfo.value = node[1].as<PetscReal>();
        return true;
    }

    // for BCInfoHolder
    Node convert<BCInfoHolder>::encode(const BCInfoHolder &bc)
    {
        Node    node(NodeType::Sequence);
        for(auto loc: bc)
        {
            Node childNode(NodeType::Map);
            childNode["location"] = loc.first;
            for(auto comp: loc.second)
                childNode[fd2str[comp.first]] = comp.second;
            node.push_back(childNode);
        }
        return node;
    }

    bool convert<BCInfoHolder>::decode(const Node &node, BCInfoHolder &bc)
    {
        if (! node.IsDefined()) return false;

        for(auto loc: node)
            for(auto comp: loc)
                if (comp.first.as<std::string>() != "location")
                    bc[loc["location"].as<BCLoc>()][comp.first.as<Field>()]
                    		= comp.second.as<BCTypeValuePair>();
        return true;
    }

    // for Perturbation
    Node convert<Perturbation>::encode(const Perturbation &pertb)
    {
        YAML::Node  node(NodeType::Map);
        node["frequency"] = pertb.freq;
        node["amplitude"] = pertb.amp;
        return node;
    }

    bool convert<Perturbation>::decode(const Node &node, Perturbation &pertb)
    {
        if ((! node["frequency"].IsDefined()) || 
                (! node["amplitude"].IsDefined())) return false;

        pertb.freq = node["frequency"].as<PetscReal>();
        pertb.amp = node["amplitude"].as<PetscReal>();
        return true;
    }

    // for OutputType
    Node convert<OutputType>::encode(const OutputType &out)
    {
        YAML::Node  node;
        node = out2str[out];
        return node;
    }

    bool convert<OutputType>::decode(const Node &node, OutputType &out)
    {
        if (! node.IsDefined()) return false;

        out = str2out[node.as<std::string>()];
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

    // for OutputInfo
    Node convert<OutputInfo>::encode(const OutputInfo &output)
    {
        YAML::Node  node(NodeType::Map);
        node["outputFormat"] = output.format;
        node["outputFlux"] = output.outputFlux;
        node["outputVelocity"] = output.outputVelocity;
        return node;
    }

    bool convert<OutputInfo>::decode(const Node &node, OutputInfo &output)
    {
        output.format = 
        		node["outputFormat"].as<OutputType>(OutputType::Binary);
        output.outputFlux = node["outputFlux"].as<PetscBool>(PETSC_TRUE);
        output.outputVelocity = node["outputVelocity"].as<PetscBool>(PETSC_FALSE);
        return true;
    }

    // for LinSolverInfo
    Node convert<LinSolverInfo>::encode(const LinSolverInfo &solver)
    {
        YAML::Node  node(NodeType::Map);
        node["type"] = solver.type;
        node["config"] = solver.config.string();
        return node;
    }

    bool convert<LinSolverInfo>::decode(const Node &node, LinSolverInfo &solver)
    {
        solver.type = node["type"].as<ExecuteType>(ExecuteType::CPU);
        solver.config = node["config"].as<std::string>("");
        return true;
    }

    // for SchemeInfo
    Node convert<SchemeInfo>::encode(const SchemeInfo &scheme)
    {
        YAML::Node  node(NodeType::Map);
        node["ibm"] = scheme.ibm;
        node["convection"] = scheme.convection;
        node["diffusion"] = scheme.diffusion;
        return node;
    }

    bool convert<SchemeInfo>::decode(const Node &node, SchemeInfo &scheme)
    {
        scheme.ibm = node["ibm"].as<IBMethod>(IBMethod::NAVIER_STOKES);
        scheme.convection = node["convection"].as<TimeScheme>(TimeScheme::ADAMS_BASHFORTH_2);
        scheme.diffusion = node["diffusion"].as<TimeScheme>(TimeScheme::CRANK_NICOLSON);
        return true;
    }

    // for SteppingInfo
    Node convert<SteppingInfo>::encode(const SteppingInfo &stepping)
    {
        YAML::Node  node(NodeType::Map);
        node["dt"] = stepping.dt;
        node["startStep"] = stepping.nStart;
        node["nt"] = stepping.nTotal;
        node["nsave"] = stepping.nSave;
        node["nrestart"] = stepping.nRestart;
        return node;
    }

    bool convert<SteppingInfo>::decode(const Node &node, SteppingInfo &stepping)
    {
        stepping.dt = node["dt"].as<PetscReal>();
        stepping.nStart = node["startStep"].as<PetscInt>(0);
        stepping.nTotal = node["nt"].as<PetscInt>(1);
        stepping.nSave = node["nsave"].as<PetscInt>(1);
        stepping.nRestart = node["nrestart"].as<PetscInt>(1);
        return true;
    }
    
} // end of namespace YAML

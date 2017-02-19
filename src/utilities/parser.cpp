/*
 * parser.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

// here goes our own headers
# include "parser.h"
# include "misc.h"


/** \copydoc parser::parseYAMLConfigFile */
PetscErrorCode parser::parseYAMLConfigFile(
        const std::string &dir, YAML::Node &node)
{
    PetscFunctionBeginUser;
    namespace fs = std::experimental::filesystem;

    PetscErrorCode  ierr;

    fs::path p = fs::system_complete(dir);

    if (! fs::exists(p.append("config.yaml")))
        SETERRQ1(PETSC_COMM_WORLD, 65, 
                "Can not find config.yaml in %s !!", p.parent_path().c_str());

    node = YAML::LoadFile(p);

    node["caseDir"] = p.parent_path().c_str();
    node["configFile"] = p.filename().c_str();

    PetscFunctionReturn(0);
}


/** \copydoc parser::parseFlowDescription */
PetscErrorCode parser::parseFlowDescription(
        const YAML::Node &flowNode, PetscInt &dim, PetscReal &nu, 
        PetscBool &customIC, types::RealVec1D &IC, types::Perturbation &pertb, 
        PetscInt &nBC, types::BCInfoHolder &BCInfo)
{
    PetscFunctionBeginUser;

    // viscosity
    nu = flowNode["nu"].as<PetscReal>();

    // the flag for customized initial velocity
    customIC = flowNode["initialCustomField"].as<PetscBool>(PETSC_FALSE);

    // get dimension
    dim = flowNode["initialVelocity"].size();

    // reset the size of IC to dim
    IC.resize(dim);

    // asign values to IC
    for(unsigned int i=0; i<dim; ++i)
        IC[i] = flowNode["initialVelocity"][i].as<PetscReal>();

    // perturbation
    pertb = flowNode["perturbation"].as<
        types::Perturbation>(types::Perturbation());

    // boundary conditions
    BCInfo = flowNode["boundaryConditions"].as<types::BCInfoHolder>();

    // number of bcs
    nBC = BCInfo.size();

    // check the number of BCs
    if (nBC != (2*dim))
        SETERRQ2(PETSC_COMM_WORLD, 60, "The number of BCs does not match ! "
                "There should be %d BCs, while we only get %d from YAML file.",
                2 * dim, nBC);

    PetscFunctionReturn(0);
}


/** \copydoc parser::parseMesh */
PetscErrorCode parser::parseMesh(
        const YAML::Node &meshNode, PetscInt &dim, types::RealVec1D &bg, 
        types::RealVec1D &ed, types::IntVec1D &nTotal, types::RealVec2D &dL)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // get the dimension of the mesh; no checking here
    dim = meshNode.size();

    // initialize vectors according to the dimensions
    bg.resize(dim);
    ed.resize(dim);
    nTotal.resize(dim);
    dL.resize(dim);

    // loop through all dimensions
    for(auto ax: meshNode)
    {
        // note that the order of dimensions in the YAML file is not guaranteed,
        // so we have to use some temporary variables
        PetscInt    dir, nTotalAx;
        PetscReal   bgAx, edAx;
        types::RealVec1D    dLAx;

        // parse current dimension
        ierr = parseOneAxis(ax, dir, bgAx, edAx, nTotalAx, dLAx); CHKERRQ(ierr);

        // assign results back
        bg[dir] = bgAx;
        ed[dir] = edAx;
        nTotal[dir] = nTotalAx;
        dL[dir] = dLAx;
    }

    PetscFunctionReturn(0);
}


/** \copydoc parser::parseOneAxis */
PetscErrorCode parser::parseOneAxis(
        const YAML::Node &axis, PetscInt &dir, PetscReal &bg, 
        PetscReal &ed, PetscInt &nTotal, types::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // a map to transform a string to int

    // get the direction
    dir = PetscInt(axis["direction"].as<types::Dir>());

    // get the far left boundary
    bg = axis["start"].as<PetscReal>();

    // parse sub-domains
    ierr = parseSubDomains(axis["subDomains"], bg, nTotal, ed, dL); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc parser::parseSubDomains */
PetscErrorCode parser::parseSubDomains(
        const YAML::Node &subs, const PetscReal bg,
        PetscInt &nTotal, PetscReal &ed, types::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    // initialize nTotal
    nTotal = 0;

    // initialize max
    ed = bg;

    // initialize dL
    dL = types::RealVec1D();

    // loop through all subdomains
    for(auto sub: subs)
    {
        types::RealVec1D   dLSub;  // cell sizes of the subdomains
        PetscInt    nSub;  // number of the cells of the subdomains

        // the 1st ed is passed by value, while the 2nd one is by reference
        ierr = parseOneSubDomain(sub, ed, nSub, ed, dLSub); CHKERRQ(ierr);

        // add number of the subdomain to total number of cells
        nTotal += nSub;

        // append the subdomain dL to global dL
        dL.insert(dL.end(), dLSub.begin(), dLSub.end());
    }

    PetscFunctionReturn(0);
}


/** \copydoc parser::parseOneSubDomain */
PetscErrorCode parser::parseOneSubDomain(
        const YAML::Node &sub, const PetscReal bg,
        PetscInt &n, PetscReal &ed, types::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    // get the number of the cells in this subdomain
    n = sub["cells"].as<PetscInt>();

    // get the end of the subdomain
    ed = sub["end"].as<PetscReal>();
    
    // get the stretching ratio
    PetscReal   r = sub["stretchRatio"].as<PetscReal>();

    // obtain the cell sizes
    if (std::abs(r - 1.0) <= 1e-12)
        dL = types::RealVec1D(n, (ed - bg) / n);  // uniform grid
    else
        misc::stretchGrid(bg, ed, n, r, dL);  // stretch grid

    PetscFunctionReturn(0);
}


/**
 * \brief user-defined YAML converter
 * \todo when converting our data types to YAML map, the order of keys are 
 * not retained. Need som extra work if we want to keep the order.
 */
namespace YAML
{
    // for Dir
    Node convert<types::Dir>::encode(const types::Dir &dir)
    {
        Node node;
        node = types::dir2str[dir];
        return node;
    }

    bool convert<types::Dir>::decode(const Node &node, types::Dir &dir)
    {
        if (! node.IsDefined()) return false;

        dir = types::str2dir[node.as<std::string>()];
        return true;
    }

    // for VelocityComponent
    Node convert<types::VelocityComponent>::encode(const types::VelocityComponent &vc)
    {
        Node node;
        node = types::vc2str[vc];
        return node;
    }

    bool convert<types::VelocityComponent>::decode(const Node &node, types::VelocityComponent &vc)
    {
        if (! node.IsDefined()) return false;

        vc = types::str2vc[node.as<std::string>()];
        return true;
    }

    // for BCType
    Node convert<types::BCType>::encode(const types::BCType &bc)
    {
        Node node;
        node = types::bt2str[bc];
        return node;
    }

    bool convert<types::BCType>::decode(const Node &node, types::BCType &bc)
    {
        if (! node.IsDefined()) return false;
        bc = types::str2bt[node.as<std::string>()];
        return true;
    }

    // for BCLoc
    Node convert<types::BCLoc>::encode(const types::BCLoc &loc)
    {
        Node node;
        node = types::bl2str[loc];
        return node;
    }

    bool convert<types::BCLoc>::decode(const Node &node, types::BCLoc &loc)
    {
        if (! node.IsDefined()) return false;

        loc = types::str2bl[node.as<std::string>()];
        return true;
    }

    // for TimeScheme
    Node convert<types::TimeScheme>::encode(const types::TimeScheme &ts)
    {
        Node node;
        node = types::ts2str[ts];
        return node;
    }

    bool convert<types::TimeScheme>::decode(const Node &node, types::TimeScheme &ts)
    {
        if (! node.IsDefined()) return false;
        ts = types::str2ts[node.as<std::string>()];
        return true;
    }

    // for IBMethod
    Node convert<types::IBMethod>::encode(const types::IBMethod &ibm)
    {
        Node node;
        node = types::ibm2str[ibm];
        return node;
    }

    bool convert<types::IBMethod>::decode(const Node &node, types::IBMethod &ibm)
    {
        if (! node.IsDefined()) return false;

        ibm = types::str2ibm[node.as<std::string>()];
        return true;
    }

    // for StaggeredMode
    Node convert<types::StaggeredMode>::encode(const types::StaggeredMode &sm)
    {
        Node node;
        node = types::sm2str[sm];
        return node;
    }

    bool convert<types::StaggeredMode>::decode(const Node &node, types::StaggeredMode &sm)
    {
        if (! node.IsDefined()) return false;

        sm = types::str2sm[node.as<std::string>()];
        return true;
    }

    // for ExecutiveType
    Node convert<types::ExecuteType>::encode(const types::ExecuteType &et)
    {
        Node node;
        node = types::et2str[et];
        return node;
    }

    bool convert<types::ExecuteType>::decode(const Node &node, types::ExecuteType &et)
    {
        if (! node.IsDefined()) return false;
        
        et = types::str2et[node.as<std::string>()];
        return true;
    }

    // for BCTypeValuePair
    Node convert<types::BCTypeValuePair>::encode(const types::BCTypeValuePair &bcInfo)
    {
        Emitter     out;
        out << Flow;
        out << BeginSeq << types::bt2str[bcInfo.type] << bcInfo.value << EndSeq;
        return Load(out.c_str());
    }

    bool convert<types::BCTypeValuePair>::decode(const Node &node, types::BCTypeValuePair &bcInfo)
    {
        if((! node[0].IsDefined()) || (! node[1].IsDefined())) return false;

        bcInfo.type = node[0].as<types::BCType>();
        bcInfo.value = node[1].as<PetscReal>();
        return true;
    }

    // for BCInfoHolder
    Node convert<types::BCInfoHolder>::encode(const types::BCInfoHolder &bc)
    {
        Node    node(NodeType::Sequence);
        int     counter = 0;
        for(auto loc: bc)
        {
            Node childNode(NodeType::Map);
            childNode["location"] = loc.first;
            for(auto comp: loc.second)
                childNode[types::vc2str[comp.first]] = comp.second;
            node.push_back(childNode);
        }
        return node;
    }

    bool convert<types::BCInfoHolder>::decode(const Node &node, types::BCInfoHolder &bc)
    {
        if (! node.IsDefined()) return false;

        for(auto loc: node)
            for(auto comp: loc)
                if (comp.first.as<std::string>() != "location")
                    bc
                        [loc["location"].as<types::BCLoc>()]
                        [comp.first.as<types::VelocityComponent>()]
                            = comp.second.as<types::BCTypeValuePair>();
        return true;
    }

    // for Perturbation
    Node convert<types::Perturbation>::encode(const types::Perturbation &pertb)
    {
        YAML::Node  node(NodeType::Map);
        node["frequency"] = pertb.freq;
        node["amplitude"] = pertb.amp;
        return node;
    }

    bool convert<types::Perturbation>::decode(const Node &node, types::Perturbation &pertb)
    {
        if ((! node["frequency"].IsDefined()) || 
                (! node["amplitude"].IsDefined())) return false;

        pertb.freq = node["frequency"].as<PetscReal>();
        pertb.amp = node["amplitude"].as<PetscReal>();
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
    
}

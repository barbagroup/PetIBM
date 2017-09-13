/*
 * parser.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

// STL
# include <string>

// here goes our own headers
# include "petibm/parser.h"
# include "petibm/misc.h"

// private function. Read config.yaml and other files for overwriting
PetscErrorCode readYAMLs(YAML::Node &node);

// private function. Read a single YAML file overwriting part of config.yaml
PetscErrorCode readSingleYAML(YAML::Node &node, const std::string &s);


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


namespace petibm
{
namespace utilities
{
namespace parser
{

// get all settings into a sinfle YAML node
PetscErrorCode getSettings(YAML::Node &node)
{
	PetscFunctionBeginUser;

    PetscErrorCode  ierr;

	// use whole new YAML node
	node = YAML::Node();

	char 		s[PETSC_MAX_PATH_LEN];
	PetscBool 	flg;

    // directory: the working directory. Default is the current directory
    node["directory"] = "./";

	ierr = PetscOptionsGetString(nullptr, nullptr, "-directory", 
			s, sizeof(s), &flg); CHKERRQ(ierr);

	if (flg) node["directory"] = s;

    // config: location of config.yaml. Default is under worling directory.
    // TODO: what if users provide a relative path? Where should it relative to?
	node["config.yaml"] = node["directory"].as<std::string>() + "/config.yaml";

	ierr = PetscOptionsGetString(nullptr, nullptr, "-config",
			s, sizeof(s), &flg); CHKERRQ(ierr);
	
	if (flg) node["config.ymal"] = s;

	// the following four arguments will overwrite corresponding sections in
	// config.yaml, if users pass these argument through command line

    // mesh: mesh.yaml. No default value.
    // TODO: what if users provide a relative path? Where should it relative to?
	ierr = PetscOptionsGetString(nullptr, nullptr, "-mesh",
			s, sizeof(s), &flg); CHKERRQ(ierr);

	if (flg) node["mesh.yaml"] = s;

    // flow: flow.yaml. No default value.
    // TODO: what if users provide a relative path? Where should it relative to?
	ierr = PetscOptionsGetString(nullptr, nullptr, "-flow",
			s, sizeof(s), &flg); CHKERRQ(ierr);

	if (flg) node["flow.ymal"] = s;

    // parameters: parameters.yaml. No default value.
    // TODO: what if users provide a relative path? Where should it relative to?
	ierr = PetscOptionsGetString(nullptr, nullptr, "-parameters",
			s, sizeof(s), &flg); CHKERRQ(ierr);

	if (flg) node["parameters.ymal"] = s;

    // bodies: bodies.yaml. No default value.
    // TODO: what if users provide a relative path? Where should it relative to?
	ierr = PetscOptionsGetString(nullptr, nullptr, "-bodies",
			s, sizeof(s), &flg); CHKERRQ(ierr);

	if (flg) node["bodies.ymal"] = s;

    // solution: path to solution folder. Always under working directory.
   	node["solution"] = node["directory"].as<std::string>() + "/solution"; 

	// read setting from YAML files
	ierr = readYAMLs(node); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc parser::parseSimulationParameters */
PetscErrorCode parseSimulationParameters(const YAML::Node &param,
                                         types::OutputInfo  &output,
                                         types::LinSolverInfo &velocity,
                                         types::LinSolverInfo &poisson,
                                         types::SchemeInfo &methods,
                                         types::SteppingInfo &stepping)
{
    PetscFunctionBeginUser;

    output = param.as<types::OutputInfo>();
    velocity = param["velocitySolver"].as<types::LinSolverInfo>();
    poisson = param["poissonSolver"].as<types::LinSolverInfo>();
    methods = param.as<types::SchemeInfo>();
    stepping = param.as<types::SteppingInfo>();

    PetscFunctionReturn(0);
}



/** \copydoc parser::parseFlowDescription */
PetscErrorCode parseFlowDescription(const YAML::Node &flowNode,
                                    PetscInt &dim, PetscReal &nu,
                                    PetscBool &customIC, types::RealVec1D &IC,
                                    types::Perturbation &pertb, PetscInt &nBC,
                                    types::BCInfoHolder &BCInfo)
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
    for(int i=0; i<dim; ++i)
        IC[i] = flowNode["initialVelocity"][i].as<PetscReal>();

    // perturbation
    pertb = flowNode["perturbation"].as<types::Perturbation>(
    		types::Perturbation());

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
PetscErrorCode parseMesh(const YAML::Node &meshNode, PetscInt &dim,
                         types::RealVec1D &bg, types::RealVec1D &ed,
                         types::IntVec1D &nTotal, types::RealVec2D &dL)
{
    using namespace types;

    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // get the dimension of the mesh; no checking here
    dim = meshNode.size();

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
PetscErrorCode parseOneAxis(const YAML::Node &axis, PetscInt &dir,
                            PetscReal &bg, PetscReal &ed, PetscInt &nTotal,
                            types::RealVec1D &dL)
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
PetscErrorCode parseSubDomains(const YAML::Node &subs, const PetscReal bg,
                               PetscInt &nTotal, PetscReal &ed,
                               types::RealVec1D &dL)
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
        PetscInt           nSub;  // number of the cells of the subdomains

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
PetscErrorCode parseOneSubDomain(const YAML::Node &sub, const PetscReal bg,
                                 PetscInt &n, PetscReal &ed,
                                 types::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    // get the number of the cells in this subdomain
    n = sub["cells"].as<PetscInt>();

    // get the end of the subdomain
    ed = sub["end"].as<PetscReal>();
    
    // get the stretching ratio
    PetscReal r = sub["stretchRatio"].as<PetscReal>();

    // obtain the cell sizes
    if (std::abs(r - 1.0) <= 1e-12)
        dL = types::RealVec1D(n, (ed - bg) / n);  // uniform grid
    else
        misc::stretchGrid(bg, ed, n, r, dL);  // stretch grid

    PetscFunctionReturn(0);
}

} // end of namespace parser
} // end of namespace utilities
} // end of namespace petibm


// private function. Read config.yaml and other files for overwriting
PetscErrorCode readYAMLs(YAML::Node &node)
{
	PetscFunctionBeginUser;

	PetscErrorCode 	ierr;
	YAML::Node 	temp;

	// open config.yaml. If failed, return error through PETSc and terminate.
	try
	{
		temp = YAML::LoadFile(node["config.yaml"].as<std::string>());
	}
	catch(YAML::BadFile &err)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
				"Unable to open %s"
				"when reading settings from config.yaml\n",
				node["config.yaml"].as<std::string>().c_str());
	}

	// copy what we read from config.yaml to the `node`
	for(auto it: temp) node[it.first] = it.second;

	// overwrite mesh section if mesh.yaml is provided
	ierr = readSingleYAML(node, "mesh"); CHKERRQ(ierr);

	// overwrite flow section if flow.yaml is provided
	ierr = readSingleYAML(node, "flow"); CHKERRQ(ierr);

	// overwrite parameters section if parameters.yaml is provided
	ierr = readSingleYAML(node, "parameters"); CHKERRQ(ierr);

	// overwrite bodies section if bodies.yaml is provided
	ierr = readSingleYAML(node, "bodies"); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


// private function. Read from a single YAML file that will overwrite a part 
// of config.yaml
PetscErrorCode readSingleYAML(YAML::Node &node, const std::string &s)
{
	PetscFunctionBeginUser;

	YAML::Node 	temp;

	if (node[s+".yaml"].IsDefined())
	{
		try
		{
			temp = YAML::LoadFile(node[s+".yaml"].as<std::string>());
		}
		catch(YAML::BadFile &err)
		{
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
					"Unable to open %s"
					"when reading settings from %s.yaml\n",
					node[s+".yaml"].as<std::string>().c_str(), s.c_str());
		}

		if (! node[s].IsDefined()) node[s] = YAML::Node();
		for(auto it: temp) node[s][it.first] = it.second;
	}

	PetscFunctionReturn(0);
}

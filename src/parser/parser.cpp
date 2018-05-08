/**
 * \file parser.cpp
 * \brief Implementations of parser functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

// STL
# include <string>
# include <sys/stat.h>

// here goes our own headers
# include <petibm/parser.h>
# include <petibm/misc.h>


namespace // anonymous namespace for internal linkage
{

// private function. Read a single YAML file overwriting part of config.yaml
PetscErrorCode readSingleYAML(YAML::Node &node, const std::string &s)
{
    PetscFunctionBeginUser;

    YAML::Node  temp;

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
} // readSingleYAML


// private function. Read config.yaml and other files for overwriting
PetscErrorCode readYAMLs(YAML::Node &node)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;
    YAML::Node temp;

    // open config.yaml. If failed, return error through PETSc and terminate.
    try
    {
        temp = YAML::LoadFile(node["config.yaml"].as<std::string>());
    }
    catch(YAML::BadFile &err)
    {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
            "Unable to open %s "
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
} // readYAMLs


// private function. Check if a directory exists. Create one if not.
PetscErrorCode checkAndCreateFolder(const std::string &dir)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    PetscBool           flg;

    ierr = PetscTestDirectory(dir.c_str(), 'w', &flg); CHKERRQ(ierr);

    if (! flg)
    {
        PetscMPIInt     rank;
        ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

        if (rank == 0)
        {
            int status;
            status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

            if (status != 0)
                SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
                        "Could not create the folder \"%s\".\n",
                        dir.c_str());
        }

        ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
} // checkAndCreateFolder
} // end of anonymous namespace


// A supplement to YAML-CPP that adds converters of our user-defined types.
namespace YAML
{

    // converter for `types::Dir`
    template <>
    struct convert<petibm::type::Dir>
    {
        static Node encode(const petibm::type::Dir &dir)
        {
            Node node;
            node = petibm::type::dir2str[dir];
            return node;
        }

        static bool decode(const Node &node, petibm::type::Dir &dir)
        {
            if (! node.IsDefined()) return false;

            dir = petibm::type::str2dir[node.as<std::string>()];
            return true;
        }
    };

    // converter for `types::Field`
    template <>
    struct convert<petibm::type::Field>
    {
        static Node encode(const petibm::type::Field &vc)
        {
            Node node;
            node = petibm::type::fd2str[vc];
            return node;
        }

        static bool decode(const Node &node, petibm::type::Field &vc)
        {
            if (! node.IsDefined()) return false;

            vc = petibm::type::str2fd[node.as<std::string>()];
            return true;
        }
    };

    // converter for `types::BCType`
    template <>
    struct convert<petibm::type::BCType>
    {
        static Node encode(const petibm::type::BCType &bc)
        {
            Node node;
            node = petibm::type::bt2str[bc];
            return node;
        }

        static bool decode(const Node &node, petibm::type::BCType &bc)
        {
            if (! node.IsDefined()) return false;
            bc = petibm::type::str2bt[node.as<std::string>()];
            return true;
        }
    };

    // converter for `types::BCLoc`
    template <>
    struct convert<petibm::type::BCLoc>
    {
        static Node encode(const petibm::type::BCLoc &loc)
        {
            Node node;
            node = petibm::type::bl2str[loc];
            return node;
        }

        static bool decode(const Node &node, petibm::type::BCLoc &loc)
        {
            if (! node.IsDefined()) return false;

            loc = petibm::type::str2bl[node.as<std::string>()];
            return true;
        }
    };

    // converter for `PetscBool`
    template<>
    struct convert<PetscBool>
    {
        static Node encode(const PetscBool &b)
        {
            YAML::Node  node;
            node = bool(b);
            return node;
        }

        static bool decode(const Node &node, PetscBool &b)
        {
            if (! node.IsDefined()) return false;

            b = PetscBool(node.as<bool>());
            return true;
        }
    };

} // end of namespace YAML


namespace petibm
{
namespace parser
{

// get all settings into a single YAML node
PetscErrorCode getSettings(YAML::Node &node)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    // use whole new YAML node
    node = YAML::Node();

    char        s[PETSC_MAX_PATH_LEN];
    PetscBool   flg;

    // directory: the working directory. Default is the current directory
    node["directory"] = "./";

    ierr = PetscOptionsGetString(nullptr, nullptr, "-directory",
            s, sizeof(s), &flg); CHKERRQ(ierr);

    if (flg) node["directory"] = s;

    // config: location of config.yaml. Default is under working directory.
    // TODO: what if users provide a relative path? Where should it relative to?
    node["config.yaml"] = node["directory"].as<std::string>() + "/config.yaml";

    ierr = PetscOptionsGetString(nullptr, nullptr, "-config",
            s, sizeof(s), &flg); CHKERRQ(ierr);

    if (flg) node["config.yaml"] = s;

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

    if (flg) node["flow.yaml"] = s;

    // parameters: parameters.yaml. No default value.
    // TODO: what if users provide a relative path? Where should it relative to?
    ierr = PetscOptionsGetString(nullptr, nullptr, "-parameters",
            s, sizeof(s), &flg); CHKERRQ(ierr);

    if (flg) node["parameters.yaml"] = s;

    // bodies: bodies.yaml. No default value.
    // TODO: what if users provide a relative path? Where should it relative to?
    ierr = PetscOptionsGetString(nullptr, nullptr, "-bodies",
            s, sizeof(s), &flg); CHKERRQ(ierr);

    if (flg) node["bodies.yaml"] = s;

    // solution: path to solution folder. Always under working directory.
    node["solution"] = node["directory"].as<std::string>() + "/solution";
    ierr = checkAndCreateFolder(node["solution"].as<std::string>());

    // read setting from YAML files
    ierr = readYAMLs(node); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getSettings


PetscErrorCode parseMesh(const YAML::Node &meshNode, PetscInt &dim,
                         type::RealVec1D &bg, type::RealVec1D &ed,
                         type::IntVec1D &nTotal, type::RealVec2D &dL)
{
    using namespace type;

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
        type::RealVec1D    dLAx;

        // parse current dimension
        ierr = parseOneAxis(ax, dir, bgAx, edAx, nTotalAx, dLAx); CHKERRQ(ierr);

        // assign results back
        bg[dir] = bgAx;
        ed[dir] = edAx;
        nTotal[dir] = nTotalAx;
        dL[dir] = dLAx;
    }

    PetscFunctionReturn(0);
} // parseMesh


PetscErrorCode parseOneAxis(const YAML::Node &axis, PetscInt &dir,
                            PetscReal &bg, PetscReal &ed, PetscInt &nTotal,
                            type::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // a map to transform a string to int

    // get the direction
    dir = PetscInt(axis["direction"].as<type::Dir>());

    // get the far left boundary
    bg = axis["start"].as<PetscReal>();

    // parse sub-domains
    ierr = parseSubDomains(axis["subDomains"], bg, nTotal, ed, dL); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // parseOneAxis


PetscErrorCode parseSubDomains(const YAML::Node &subs, const PetscReal bg,
                               PetscInt &nTotal, PetscReal &ed,
                               type::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    // initialize nTotal
    nTotal = 0;

    // initialize max
    ed = bg;

    // initialize dL
    dL = type::RealVec1D();

    // loop through all sub-domains
    for(auto sub: subs)
    {
        type::RealVec1D   dLSub;  // cell sizes of the sub-domains
        PetscInt           nSub;  // number of the cells of the sub-domains

        // the 1st ed is passed by value, while the 2nd one is by reference
        ierr = parseOneSubDomain(sub, ed, nSub, ed, dLSub); CHKERRQ(ierr);

        // add number of the sub-domain to total number of cells
        nTotal += nSub;

        // append the sub-domain dL to global dL
        dL.insert(dL.end(), dLSub.begin(), dLSub.end());
    }

    PetscFunctionReturn(0);
} // parseSubDomains


PetscErrorCode parseOneSubDomain(const YAML::Node &sub, const PetscReal bg,
                                 PetscInt &n, PetscReal &ed,
                                 type::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    // get the number of the cells in this sub-domain
    n = sub["cells"].as<PetscInt>();

    // get the end of the sub-domain
    ed = sub["end"].as<PetscReal>();

    // get the stretching ratio
    PetscReal r = sub["stretchRatio"].as<PetscReal>();

    // obtain the cell sizes
    if (std::abs(r - 1.0) <= 1e-12)
        dL = type::RealVec1D(n, (ed - bg) / n);  // uniform grid
    else
        misc::stretchGrid(bg, ed, n, r, dL);  // stretch grid

    PetscFunctionReturn(0);
} // parseOneSubDomain


// get information about if we have periodic BCs from YAML node
PetscErrorCode parseBCs(const YAML::Node &node,
        type::IntVec2D &bcTypes, type::RealVec2D &bcValues)
{
    PetscFunctionBeginUser;

    if (! node["flow"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"flow\" in the YAML "
                "node passed to parseBCs.\n");

    if (! node["flow"]["boundaryConditions"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"boundaryConditions\" under the key "
                "\"flow\" in the YAML node passed to parseBCs.\n");

    bcTypes = type::IntVec2D(3, type::IntVec1D(6, int(type::BCType::NOBC)));
    bcValues = type::RealVec2D(3, type::RealVec1D(6, 0.0));

    for(auto sub: node["flow"]["boundaryConditions"])
    {
        int loc = int(sub["location"].as<type::BCLoc>());

        for(auto ssub: sub)
        {
            if (ssub.first.as<std::string>() != "location")
            {
                int f = int(ssub.first.as<type::Field>());
                bcTypes[f][loc] = int(ssub.second[0].as<type::BCType>());
                bcValues[f][loc] = ssub.second[1].as<PetscReal>();
            }
        }
    }

    PetscFunctionReturn(0);
} // parseBCs


// get initial conditions
PetscErrorCode parseICs(const YAML::Node &node, type::RealVec1D &icValues)
{
    PetscFunctionBeginUser;

    if (! node["flow"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"flow\" in the YAML "
                "node passed to parseICs.\n");

    if (! node["flow"]["initialVelocity"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"initialVelocity\" under the key "
                "\"flow\" in the YAML node passed to parseICs.\n");

    const YAML::Node &temp = node["flow"]["initialVelocity"];

    icValues = type::RealVec1D(3, 0.0);

    for(unsigned int i=0; i<temp.size(); i++) icValues[i] = temp[i].as<PetscReal>();

    PetscFunctionReturn(0);
} // parseICs

} // end of namespace parser
} // end of namespace petibm

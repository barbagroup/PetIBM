/**
 * \file parser.cpp
 * \brief Implementations of parser functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

// STL
#include <sys/stat.h>
#include <string>
#include <vector>

// symengine
#include <symengine/parser.h>
#include <symengine/expression.h>
#include <symengine/lambda_double.h>

// here goes our own headers
#include <petibm/misc.h>
#include <petibm/parser.h>

namespace  // anonymous namespace for internal linkage
{
// private function. Create a directory if not already existing.
PetscErrorCode createDirectory(const std::string &dir)
{
    PetscFunctionBeginUser;

    if (mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    {
        if (errno != EEXIST)  // if error != "File exists"
        {
            SETERRQ2(PETSC_COMM_WORLD, 1,
                     "Could not create the folder \"%s\" (%s).\n", dir.c_str(),
                     strerror(errno));
        }
    }

    PetscFunctionReturn(0);
}  // createDirectory
}  // end of anonymous namespace

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
        if (!node.IsDefined()) return false;

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
        if (!node.IsDefined()) return false;

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
        if (!node.IsDefined()) return false;
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
        if (!node.IsDefined()) return false;

        loc = petibm::type::str2bl[node.as<std::string>()];
        return true;
    }
};

// converter for `PetscBool`
template <>
struct convert<PetscBool>
{
    static Node encode(const PetscBool &b)
    {
        YAML::Node node;
        node = bool(b);
        return node;
    }

    static bool decode(const Node &node, PetscBool &b)
    {
        if (!node.IsDefined()) return false;

        b = PetscBool(node.as<bool>());
        return true;
    }
};

}  // end of namespace YAML

namespace petibm
{
namespace parser
{
// Load nodes from a given YAML file.
PetscErrorCode readYAMLFile(const std::string &filePath, YAML::Node &node)
{
    YAML::Node tmp;

    PetscFunctionBeginUser;

    // Load the content of the given YAML file.
    try
    {
        tmp = YAML::LoadFile(filePath);
    }
    catch (YAML::BadFile &err)
    {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Unable to open %s",
                 filePath.c_str());
    }
    // Add new nodes, overwrite existing ones.
    for (auto item : tmp) node[item.first.as<std::string>()] = item.second;

    PetscFunctionReturn(0);
}  // readYAMLFile

// Get all settings into a single YAML node.
PetscErrorCode getSettings(YAML::Node &node)
{
    PetscErrorCode ierr;
    char s[PETSC_MAX_PATH_LEN];
    PetscBool flag;
    std::string filePath;

    PetscFunctionBeginUser;

    // Get the working directory; default is the current working directory.
    node["directory"] = "./";
    ierr = PetscOptionsGetString(nullptr, nullptr, "-directory", s, sizeof(s),
                                 &flag); CHKERRQ(ierr);
    if (flag) node["directory"] = s;

    // Get the output directory where to save field solutions;
    // default is the `output` folder under the working directory.
    // If the user provides a relative path, it is relative to $PWD.
    node["output"] = node["directory"].as<std::string>() + "/output";
    ierr =
        PetscOptionsGetString(nullptr, nullptr, "-output", s, sizeof(s), &flag);
    CHKERRQ(ierr);
    if (flag) node["output"] = s;
    ierr = createDirectory(node["output"].as<std::string>()); CHKERRQ(ierr);

    // Get the directory where to save PETSc logging files;
    // default is the `logs` folder under the working directory.
    // TODO: if user provides a relative path, where should it be relative to?
    node["logs"] = node["output"].as<std::string>() + "/logs";
    ierr =
        PetscOptionsGetString(nullptr, nullptr, "-logs", s, sizeof(s), &flag);
    CHKERRQ(ierr);
    if (flag) node["logs"] = s;
    ierr = createDirectory(node["logs"].as<std::string>()); CHKERRQ(ierr);

    // Get the path of the global YAML configuration file;
    // default is the file `config.yaml` in the working directory.
    // TODO: if user provides a relative path, where should it be relative to?
    filePath = node["directory"].as<std::string>() + "/config.yaml";
    ierr = PetscOptionsGetString(
        nullptr, nullptr, "-config", s, sizeof(s), &flag); CHKERRQ(ierr);
    if (flag) filePath = s;
    // Load the settings into the YAML node.
    ierr = readYAMLFile(filePath, node); CHKERRQ(ierr);

    // Load YAML configuration files provided through the command-line.
    std::vector<const char *> cliOpts = {"-mesh", "-flow", "-parameters",
                                         "-bodies", "-probes"};
    for (auto cliOpt : cliOpts)
    {
        // Get the path of the YAML configuration file.
        // TODO: if relative path, where should it be relative to?
        ierr = PetscOptionsGetString(nullptr, nullptr, cliOpt, s, sizeof(s),
                                     &flag); CHKERRQ(ierr);
        if (flag)
        {
            // Load the YAML file and overwrite existing values into the node.
            ierr = readYAMLFile(std::string(s), node); CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}  // getSettings

PetscErrorCode parseMesh(const YAML::Node &meshNode, PetscInt &dim,
                         type::RealVec1D &bg, type::RealVec1D &ed,
                         type::IntVec1D &nTotal, type::RealVec2D &dL)
{
    using namespace type;

    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    // get the dimension of the mesh; no checking here
    dim = meshNode.size();

    // loop through all dimensions
    for (auto ax : meshNode)
    {
        // note that the order of dimensions in the YAML file is not guaranteed,
        // so we have to use some temporary variables
        PetscInt dir, nTotalAx;
        PetscReal bgAx, edAx;
        type::RealVec1D dLAx;

        // parse current dimension
        ierr = parseOneAxis(ax, dir, bgAx, edAx, nTotalAx, dLAx);
        CHKERRQ(ierr);

        // assign results back
        bg[dir] = bgAx;
        ed[dir] = edAx;
        nTotal[dir] = nTotalAx;
        dL[dir] = dLAx;
    }

    PetscFunctionReturn(0);
}  // parseMesh

PetscErrorCode parseOneAxis(const YAML::Node &axis, PetscInt &dir,
                            PetscReal &bg, PetscReal &ed, PetscInt &nTotal,
                            type::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    // a map to transform a string to int

    // get the direction
    dir = PetscInt(axis["direction"].as<type::Dir>());

    // get the far left boundary
    bg = axis["start"].as<PetscReal>();

    // parse sub-domains
    ierr = parseSubDomains(axis["subDomains"], bg, nTotal, ed, dL);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // parseOneAxis

PetscErrorCode parseSubDomains(const YAML::Node &subs, const PetscReal bg,
                               PetscInt &nTotal, PetscReal &ed,
                               type::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    // initialize nTotal
    nTotal = 0;

    // initialize max
    ed = bg;

    // initialize dL
    dL = type::RealVec1D();

    // loop through all sub-domains
    for (auto sub : subs)
    {
        type::RealVec1D dLSub;  // cell sizes of the sub-domains
        PetscInt nSub;          // number of the cells of the sub-domains

        // the 1st ed is passed by value, while the 2nd one is by reference
        ierr = parseOneSubDomain(sub, ed, nSub, ed, dLSub); CHKERRQ(ierr);

        // add number of the sub-domain to total number of cells
        nTotal += nSub;

        // append the sub-domain dL to global dL
        dL.insert(dL.end(), dLSub.begin(), dLSub.end());
    }

    PetscFunctionReturn(0);
}  // parseSubDomains

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
}  // parseOneSubDomain

// get information about if we have periodic BCs from YAML node
PetscErrorCode parseBCs(const YAML::Node &node, type::IntVec2D &bcTypes,
                        type::RealVec2D &bcValues)
{
    PetscFunctionBeginUser;

    if (!node["flow"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"flow\" in the YAML "
                "node passed to parseBCs.\n");

    if (!node["flow"]["boundaryConditions"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"boundaryConditions\" under the key "
                "\"flow\" in the YAML node passed to parseBCs.\n");

    bcTypes = type::IntVec2D(3, type::IntVec1D(6, int(type::BCType::NOBC)));
    bcValues = type::RealVec2D(3, type::RealVec1D(6, 0.0));

    for (auto sub : node["flow"]["boundaryConditions"])
    {
        int loc = int(sub["location"].as<type::BCLoc>());

        for (auto ssub : sub)
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
}  // parseBCs

// get initial conditions
PetscErrorCode parseICs(
        const YAML::Node &node,
        std::vector<SymEngine::LambdaRealDoubleVisitor> &lambdas)
{
    PetscFunctionBeginUser;

    if (!node["flow"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"flow\" in the YAML "
                "node passed to parseICs.\n");

    if (!node["flow"]["initialVelocity"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"initialVelocity\" under the key "
                "\"flow\" in the YAML node passed to parseICs.\n");

    const YAML::Node &temp = node["flow"]["initialVelocity"];

    lambdas = std::vector<SymEngine::LambdaRealDoubleVisitor>(3);

    // symbols for independent variables
    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");
    auto t = SymEngine::symbol("t");
    auto nu = SymEngine::symbol("nu");

    for (unsigned int i = 0; i < temp.size(); i++)
    {
        // get the pointer to the symbolic expression
        auto expr = SymEngine::parse(temp[i].as<std::string>());

        // assign the lambdified objects to the vector
        lambdas[i].init({x, y, z, t, nu}, *expr);
    }

    PetscFunctionReturn(0);
}  // parseICs

}  // end of namespace parser
}  // end of namespace petibm

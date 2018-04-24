/**
 * \file createxdmf/main.cpp
 * \brief An utility that generates XDMF files for visualization.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see createxdmf
 * \ingroup createxdmf
 */


# include <string>
# include <fstream>
# include <sstream>
# include <petscsys.h>
# include <yaml-cpp/yaml.h>
# include <petibm/type.h>
# include <petibm/parser.h>
# include <petibm/mesh.h>


/**
 * \defgroup createxdmf Post-processing utility: createxdmf
 * \brief A post-processing utility that creates XDMF files for visualization.
 *
 * This is a helper utility built with PetIBM components, and it creates XDMF
 * files based on the HDF5 files produced by
 * \ref nssolver "Navier-Stokes solver",
 * \ref tairacolonius "IBPM solver", and
 * \ref decoupledibpm "decoupled IBPM solver".
 * XDMF files can be used for visualization in
 * [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)
 *
 * If readers are interested in using this utility,
 * please refer to
 * \ref md_doc_markdowns_runpetibm "Running PetIBM",
 * \ref md_doc_markdowns_examples2d "2D Exmaples", and
 * \ref md_doc_markdowns_examples3d "3D Examples".
 *
 * \ingroup apps
 */

PetscErrorCode writeSingleXDMF(
        const std::string &directory, const std::string &name,
        const PetscInt &dim, const petibm::type::IntVec1D &n,
        const PetscInt &bg, const PetscInt &ed, const PetscInt &step);

int main(int argc, char **argv)
{

    PetscErrorCode      ierr;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);



    YAML::Node  setting;

    ierr = petibm::parser::getSettings(setting); CHKERRQ(ierr);

    petibm::type::Mesh  mesh;

    ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, setting, mesh); CHKERRQ(ierr);



    PetscInt                    bg;
    PetscInt                    ed;
    PetscInt                    step;
    PetscBool                   isSet;

    // get the range of solutions that are going to be calculated
    ierr = PetscOptionsGetInt(nullptr, nullptr, "-bg", &bg, &isSet); CHKERRQ(ierr);
    if (! isSet) bg = setting["parameters"]["startStep"].as<PetscInt>(0);

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-ed", &ed, &isSet); CHKERRQ(ierr);
    if (! isSet) ed = bg + setting["parameters"]["nt"].as<PetscInt>();

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-step", &step, &isSet); CHKERRQ(ierr);
    if (! isSet) step = setting["parameters"]["nsave"].as<PetscInt>();



    // u
    ierr = writeSingleXDMF(setting["directory"].as<std::string>(),
            "u", mesh->dim, mesh->n[0], bg, ed, step); CHKERRQ(ierr);

    // v
    ierr = writeSingleXDMF(setting["directory"].as<std::string>(),
            "v", mesh->dim, mesh->n[1], bg, ed, step); CHKERRQ(ierr);

    // p
    ierr = writeSingleXDMF(setting["directory"].as<std::string>(),
            "p", mesh->dim, mesh->n[3], bg, ed, step); CHKERRQ(ierr);

    // wz
    petibm::type::IntVec1D  wn(3);

    wn[0] = mesh->n[4][0]; wn[1] = mesh->n[4][1]; wn[2] = mesh->n[3][2];
    ierr = writeSingleXDMF(setting["directory"].as<std::string>(),
            "wz", mesh->dim, wn, bg, ed, step); CHKERRQ(ierr);

    if (mesh->dim == 3)
    {
        // w
        ierr = writeSingleXDMF(setting["directory"].as<std::string>(),
                "w", mesh->dim, mesh->n[2], bg, ed, step); CHKERRQ(ierr);

        // wx
        wn[0] = mesh->n[3][0]; wn[1] = mesh->n[4][1]; wn[2] = mesh->n[4][2];
        ierr = writeSingleXDMF(setting["directory"].as<std::string>(),
                "wx", mesh->dim, wn, bg, ed, step); CHKERRQ(ierr);

        // wy
        wn[0] = mesh->n[4][0]; wn[1] = mesh->n[3][1]; wn[2] = mesh->n[4][2];
        ierr = writeSingleXDMF(setting["directory"].as<std::string>(),
                "wy", mesh->dim, wn, bg, ed, step); CHKERRQ(ierr);
    }

    // manually destroy PETSc objects inside the mesh instance
    ierr = mesh->destroy(); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
} // main


PetscErrorCode writeSingleXDMF(
        const std::string &directory, const std::string &name,
        const PetscInt &dim, const petibm::type::IntVec1D &n,
        const PetscInt &bg, const PetscInt &ed, const PetscInt &step)
{
    using namespace petibm;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    PetscViewer         viewer;

    std::string         file = directory + "/" + name + ".xmf";

    ierr = PetscViewerASCIIOpen(
            PETSC_COMM_WORLD, file.c_str(), &viewer); CHKERRQ(ierr);

    // write header
    ierr = PetscViewerASCIIPrintf(viewer,
            "<?xml version=\'1.0\' ?>\n\n"); CHKERRQ(ierr);


    // write macro definitions
    ierr = PetscViewerASCIIPrintf(viewer,
            "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "\t"
            "<!ENTITY CaseDir \"%s\">\n", "./"); CHKERRQ(ierr);

    // always use 3D XDMF format, so both Visit and Paraview works
    for(int i=0; i<3; ++i)
    {
        ierr = PetscViewerASCIIPrintf(viewer, "\t"
                "<!ENTITY N%s \"%d\">\n",
                type::dir2str[type::Dir(i)].c_str(), n[i]); CHKERRQ(ierr);
    }

    // topology
    ierr = PetscViewerASCIIPrintf(viewer, "\t"
            "<!ENTITY Topo \"<Topology TopologyType=\'3DRectMesh\' "
            "Dimensions=\'&Nz; &Ny; &Nx;\'/>\">\n"); CHKERRQ(ierr);

    // geometry
    ierr = PetscViewerASCIIPrintf(viewer, "\t<!ENTITY Geo\n"); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "\t\t"
            "\"<Geometry GeometryType=\'VXVYVZ\'>\n"); CHKERRQ(ierr);

    for(int i=0; i<dim; ++i)
    {
        std::string dir = type::dir2str[type::Dir(i)];
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t"
                "<DataItem Dimensions=\'&N%s;\' Format=\'HDF\' Precision=\'8\'>\n"
                "\t\t\t\t&CaseDir;/grid.h5:/%s/%s\n"
                "\t\t\t</DataItem>\n",
                dir.c_str(), name.c_str(), dir.c_str()); CHKERRQ(ierr);
    }

    if (dim == 2) // if dim == 2, use dummy z-axis
    {
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t"
                "<DataItem Dimensions=\'&Nz;\' Format=\'XML\' Precision=\'8\'>\n"
                "\t\t\t\t0.0\n"
                "\t\t\t</DataItem>\n"); CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\t\t</Geometry>\"\n"); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "\t>\n"); CHKERRQ(ierr);

    // the end of enitiy definitions
    ierr = PetscViewerASCIIPrintf(viewer, "]>\n\n"); CHKERRQ(ierr);


    // write Xdmf block
    ierr = PetscViewerASCIIPrintf(viewer, "<Xdmf Version=\"3.0\">\n"); CHKERRQ(ierr);

    // write Domain block
    ierr = PetscViewerASCIIPrintf(viewer, "\t<Domain>\n"); CHKERRQ(ierr);


    //write temporal grid collection
    ierr = PetscViewerASCIIPrintf(viewer, "\t"
            "<Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n"); CHKERRQ(ierr);

    // write each step
    for(PetscInt t=bg; t<=ed; t+=step)
    {
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t"
                "<Grid GridType=\"Uniform\" Name=\"%s Grid\">\n", name.c_str()); CHKERRQ(ierr);

        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t"
                "<Time Value=\"%07d\" />\n", t); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t&Topo; &Geo;\n"); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t"
                "<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",
                name.c_str()); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t\t"
                "<DataItem Dimensions=\"&Nz; &Ny; &Nx;\" "); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,
                "Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">\n"); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t\t\t"
                "&CaseDir;/solution/%07d.h5:/%s\n", t, name.c_str()); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t\t"
                "</DataItem>\n"); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t\t"
                "</Attribute>\n"); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer, "\t\t" "</Grid>\n"); CHKERRQ(ierr);
    }

    ierr = PetscViewerASCIIPrintf(viewer, "\t</Grid>\n"); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "\t</Domain>\n"); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "</Xdmf>\n"); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // writeSingleXDMF

/**
 * \file vorticity/main.cpp
 * \brief An utility that calculates vorticity fields.
 * \see vorticity
 * \ingroup vorticity
 */

// STL
# include <sstream>
# include <iomanip>

// PETSc
# include <petsc.h>
# include <petscviewerhdf5.h>

// YAML-CPP
# include <yaml-cpp/yaml.h>

// PetIBM
# include <petibm/type.h>
# include <petibm/parser.h>
# include <petibm/mesh.h>
# include <petibm/boundary.h>
# include <petibm/solution.h>
# include <petibm/io.h>


/**
 * \defgroup vorticity Post-processing utility: vorticity
 * \brief A post-processing utility that calculates vorticity fields.
 * 
 * This is a helper utility built with PetIBM components, and it calculates
 * the vorticity field of simulation results from the
 * \ref nssolver "Navier-Stokes solver",
 * \ref tairacolonius "IBPM solver", and 
 * \ref decoupledibpm "decoupled IBPM solver".
 * 
 * If readers are interested in using this utility,
 * please refer to 
 * \ref md_runpetibm "Running PetIBM",
 * \ref md_examples2d "2D Examples", and
 * \ref md_examples3d "3D Examples".
 * 
 * \ingroup apps
 */


PetscErrorCode initVorticityMesh(const petibm::type::Mesh &mesh, 
        petibm::type::IntVec2D &n, petibm::type::RealVec3D &coord, 
        std::vector<DM> &wDMs, std::vector<std::string> &names, 
        std::vector<Vec> &vecs);


PetscErrorCode calculateVorticity2D(
        const petibm::type::Mesh &mesh, const petibm::type::Boundary &bds, 
        const petibm::type::Solution &soln, const std::vector<DM> &wDm, 
        std::vector<Vec> &w);


PetscErrorCode calculateVorticity3D(
        const petibm::type::Mesh &mesh, const petibm::type::Boundary &bds, 
        const petibm::type::Solution &soln, const std::vector<DM> &wDm, 
        std::vector<Vec> &w);


int main(int argc, char **argv)
{
    PetscErrorCode      ierr;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    YAML::Node                  setting;
    petibm::type::Mesh          mesh;
    petibm::type::Boundary      bc;
    petibm::type::Solution      soln;

    
    PetscInt                    bg;
    PetscInt                    ed;
    PetscInt                    step;

    petibm::type::IntVec2D      wN;
    petibm::type::RealVec3D     wCoord;
    std::vector<DM>             wDMs;
    std::vector<Vec>            w;
    std::vector<std::string>    wNames;

    PetscBool                   isSet;


    // read the config.yaml
    ierr = petibm::parser::getSettings(setting); CHKERRQ(ierr);
    
    // initialize mesh, bc, soln
    ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, setting, mesh); CHKERRQ(ierr);
    ierr = petibm::boundary::createBoundary(mesh, setting, bc); CHKERRQ(ierr);
    ierr = petibm::solution::createSolution(mesh, soln); CHKERRQ(ierr);


    // initialize vorticity mesh
    ierr = initVorticityMesh(mesh, wN, wCoord, wDMs, wNames, w); CHKERRQ(ierr);

    // output grid coordinates to existing gird.h5
    if (mesh->mpiRank == 0)
    {
        for(unsigned int i=0; i<wNames.size(); ++i)
        {
            ierr = petibm::io::writeHDF5Vecs(PETSC_COMM_SELF,
                    setting["directory"].as<std::string>() + "/grid",
                    wNames[i], {"x", "y", "z"}, wCoord[i], FILE_MODE_APPEND); 
            CHKERRQ(ierr);
        }
    }
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    // get the range of solutions that are going to be calculated
    ierr = PetscOptionsGetInt(nullptr, nullptr, "-bg", &bg, &isSet); CHKERRQ(ierr);
    if (! isSet) bg = setting["parameters"]["startStep"].as<PetscInt>(0);

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-ed", &ed, &isSet); CHKERRQ(ierr);
    if (! isSet) ed = bg + setting["parameters"]["nt"].as<PetscInt>(); 

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-step", &step, &isSet); CHKERRQ(ierr);
    if (! isSet) step = setting["parameters"]["nsave"].as<PetscInt>(); 


    // start calculating vorticity
    ierr = PetscPrintf(PETSC_COMM_WORLD, 
            "====================================\n"
            "Start to calculate vorticity fields!\n"
            "====================================\n\n"); CHKERRQ(ierr);
    
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Beginning: %d\n" 
            "End: %d\n" "Step: %d\n\n", bg, ed, step); CHKERRQ(ierr);


    for(int i=bg; i<=ed; i+=step)
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, 
                "Calculating vorticity fields for time step %d ... ", i); 
        CHKERRQ(ierr);

        // read solution
        std::stringstream   ss;
        ss << setting["solution"].as<std::string>() << "/" 
            << std::setfill('0') << std::setw(7) << i;
        ierr = soln->read(ss.str());
        CHKERRQ(ierr);

        // calculate values at ghost points based on the current solution
        ierr = bc->setGhostICs(soln); CHKERRQ(ierr);

        // calculate vorticity
        if (mesh->dim == 2)
        { ierr = calculateVorticity2D(mesh, bc, soln, wDMs, w); CHKERRQ(ierr); }
        else
        { ierr = calculateVorticity3D(mesh, bc, soln, wDMs, w); CHKERRQ(ierr); }

        // append to existing HDF5 file
        ierr = petibm::io::writeHDF5Vecs(mesh->comm, ss.str(), "/",
                wNames, w, FILE_MODE_APPEND); CHKERRQ(ierr);

        ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);
    }

    // destroy everything from PETSc
    for(unsigned int f=0; f<w.size(); ++f)
    {
        ierr = VecDestroy(&w[f]); CHKERRQ(ierr);
        ierr = DMDestroy(&wDMs[f]); CHKERRQ(ierr);
    }
    

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
} // main


PetscErrorCode calculateVorticity2D(const petibm::type::Mesh &mesh, 
        const petibm::type::Boundary &bds, const petibm::type::Solution &soln,
        const std::vector<DM> &wDm, std::vector<Vec> &w)
{
    PetscFunctionBeginUser;

    PetscErrorCode              ierr;

    std::vector<Vec>            uLocals;

    std::vector<PetscReal**>    uLocalArrys;

    PetscReal                   **wArray;

    DMDALocalInfo               info;


    // get local info for vorticity meshes
    ierr = DMDAGetLocalInfo(wDm[0], &info); CHKERRQ(ierr);

    // resize the STL vectors
    uLocals.resize(mesh->dim);
    uLocalArrys.resize(mesh->dim);

    // create local Vecs
    for(PetscInt f=0; f<mesh->dim; ++f)
    {
        ierr = DMCreateLocalVector(mesh->da[f], &uLocals[f]); CHKERRQ(ierr);
    }

    // scatter from global Vecs to local Vecs
    ierr = DMCompositeScatterArray(
            mesh->UPack, soln->UGlobal, uLocals.data()); CHKERRQ(ierr);

    // copy the values at ghost points from Boundary instance to local Vecs
    ierr = bds->copyValues2LocalVecs(uLocals); CHKERRQ(ierr);

    // get read-only raw array from local Vecs
    for(PetscInt f=0; f<mesh->dim; ++f)
    {
        ierr = DMDAVecGetArrayRead(
                mesh->da[f], uLocals[f], &uLocalArrys[f]); CHKERRQ(ierr);
    }

    // get raw array from vorticity Vecs
    ierr = DMDAVecGetArray(wDm[0], w[0], &wArray); CHKERRQ(ierr);


    // calculate z-vorticity
    for(PetscInt j=info.ys; j<info.ys+info.ym; ++j)
    {
        for(PetscInt i=info.xs; i<info.xs+info.xm; ++i)
        {
            wArray[j][i] = 
                (uLocalArrys[1][j-1][i] - uLocalArrys[1][j-1][i-1]) / 
                (mesh->coord[1][0][i] - mesh->coord[1][0][i-1]) -
                (uLocalArrys[0][j][i-1] - uLocalArrys[0][j-1][i-1]) / 
                (mesh->coord[0][1][j] - mesh->coord[0][1][j-1]);
        }
    }


    // return the ownership of vorticity raw array
    ierr = DMDAVecRestoreArray(wDm[0], w[0], &wArray); CHKERRQ(ierr);


    // return the ownership of local velocity array and destroy Vecs
    for(PetscInt f=0; f<mesh->dim; ++f)
    {
        ierr = DMDAVecRestoreArrayRead(
                mesh->da[f], uLocals[f], &uLocalArrys[f]); CHKERRQ(ierr);

        ierr = VecDestroy(&uLocals[f]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
} // calculateVorticity2D


PetscErrorCode calculateVorticity3D(const petibm::type::Mesh &mesh, 
        const petibm::type::Boundary &bds, const petibm::type::Solution &soln,
        const std::vector<DM> &wDm, std::vector<Vec> &w)
{
    PetscFunctionBeginUser;

    PetscErrorCode              ierr;

    std::vector<Vec>            uLocals;

    std::vector<PetscReal***>   uLocalArrys;

    std::vector<DMDALocalInfo>  info;

    PetscReal                   ***wArray;


    // resize the STL vectors
    info.resize(mesh->dim);
    uLocals.resize(mesh->dim);
    uLocalArrys.resize(mesh->dim);

    // create local Vecs and local infos
    for(PetscInt f=0; f<mesh->dim; ++f)
    {
        ierr = DMDAGetLocalInfo(wDm[f], &info[f]); CHKERRQ(ierr);
        ierr = DMCreateLocalVector(mesh->da[f], &uLocals[f]); CHKERRQ(ierr);
    }

    // scatter from global Vecs to local Vecs
    ierr = DMCompositeScatterArray(
            mesh->UPack, soln->UGlobal, uLocals.data()); CHKERRQ(ierr);

    // copy the values at ghost points from Boundary instance to local Vecs
    ierr = bds->copyValues2LocalVecs(uLocals); CHKERRQ(ierr);

    // get read-only raw array from local Vecs
    for(PetscInt f=0; f<mesh->dim; ++f)
    {
        ierr = DMDAVecGetArrayRead(
                mesh->da[f], uLocals[f], &uLocalArrys[f]); CHKERRQ(ierr);
    }


    // wx
    {
        // get raw array from vorticity Vecs
        ierr = DMDAVecGetArray(wDm[0], w[0], &wArray); CHKERRQ(ierr);

        // calculate x-vorticity
        for(PetscInt k=info[0].zs; k<info[0].zs+info[0].zm; ++k)
        {
            for(PetscInt j=info[0].ys; j<info[0].ys+info[0].ym; ++j)
            {
                for(PetscInt i=info[0].xs; i<info[0].xs+info[0].xm; ++i)
                {
                    wArray[k][j][i] = 

                        (uLocalArrys[2][k-1][j][i-1] - uLocalArrys[2][k-1][j-1][i-1]) / 
                        (mesh->coord[2][1][j] - mesh->coord[2][1][j-1]) -

                        (uLocalArrys[1][k][j-1][i-1] - uLocalArrys[1][k-1][j-1][i-1]) / 
                        (mesh->coord[1][2][k] - mesh->coord[1][2][k-1]);
                }
            }
        }

        // return the ownership of vorticity raw array
        ierr = DMDAVecRestoreArray(wDm[0], w[0], &wArray); CHKERRQ(ierr);
    }


    // wy
    {
        // get raw array from vorticity Vecs
        ierr = DMDAVecGetArray(wDm[1], w[1], &wArray); CHKERRQ(ierr);

        // calculate y-vorticity
        for(PetscInt k=info[1].zs; k<info[1].zs+info[1].zm; ++k)
        {
            for(PetscInt j=info[1].ys; j<info[1].ys+info[1].ym; ++j)
            {
                for(PetscInt i=info[1].xs; i<info[1].xs+info[1].xm; ++i)
                {
                    wArray[k][j][i] = 

                        (uLocalArrys[0][k][j-1][i-1] - uLocalArrys[0][k-1][j-1][i-1]) / 
                        (mesh->coord[0][2][k] - mesh->coord[0][2][k-1]) -

                        (uLocalArrys[2][k-1][j-1][i] - uLocalArrys[2][k-1][j-1][i-1]) / 
                        (mesh->coord[2][0][i] - mesh->coord[2][0][i-1]);
                }
            }
        }

        // return the ownership of vorticity raw array
        ierr = DMDAVecRestoreArray(wDm[1], w[1], &wArray); CHKERRQ(ierr);
    }


    // wz
    {
        // get raw array from vorticity Vecs
        ierr = DMDAVecGetArray(wDm[2], w[2], &wArray); CHKERRQ(ierr);

        // calculate z-vorticity
        for(PetscInt k=info[2].zs; k<info[2].zs+info[2].zm; ++k)
        {
            for(PetscInt j=info[2].ys; j<info[2].ys+info[2].ym; ++j)
            {
                for(PetscInt i=info[2].xs; i<info[2].xs+info[2].xm; ++i)
                {
                    wArray[k][j][i] = 

                        (uLocalArrys[1][k][j-1][i] - uLocalArrys[1][k][j-1][i-1]) / 
                        (mesh->coord[1][0][i] - mesh->coord[1][0][i-1]) -

                        (uLocalArrys[0][k][j][i-1] - uLocalArrys[0][k][j-1][i-1]) / 
                        (mesh->coord[0][1][j] - mesh->coord[0][1][j-1]);
                }
            }
        }

        // return the ownership of vorticity raw array
        ierr = DMDAVecRestoreArray(wDm[2], w[2], &wArray); CHKERRQ(ierr);
    }


    // return the ownership of local velocity array and destroy local Vecs
    for(PetscInt f=0; f<mesh->dim; ++f)
    {
        ierr = DMDAVecRestoreArrayRead(
                mesh->da[f], uLocals[f], &uLocalArrys[f]); CHKERRQ(ierr);

        ierr = VecDestroy(&uLocals[f]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
} // calculateVorticity3D


PetscErrorCode initVorticityMesh(const petibm::type::Mesh &mesh, 
        petibm::type::IntVec2D &n, petibm::type::RealVec3D &coord, 
        std::vector<DM> &wDMs, std::vector<std::string> &names, 
        std::vector<Vec> &vecs)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    if (mesh->dim == 2)
    {
        // initialize size of meshes
        n = petibm::type::IntVec2D(1, petibm::type::IntVec1D(3));

        // set values of n for wz
        n[0][0] = mesh->n[4][0];
        n[0][1] = mesh->n[4][1];
        n[0][2] = mesh->n[3][2]; // this value should be one

        // initialize vectors holding coordinates
        coord = petibm::type::RealVec3D(1, petibm::type::RealVec2D(3));

        // set coordinates for wz
        coord[0][0].assign(mesh->coord[4][0], mesh->coord[4][0]+n[0][0]);
        coord[0][1].assign(mesh->coord[4][1], mesh->coord[4][1]+n[0][1]);
        coord[0][2].assign(mesh->coord[3][2], mesh->coord[3][2]+n[0][2]);

        // name of vorticity fields
        names = std::vector<std::string>(1);

        // set name for wz
        names[0] = "wz";

        // initialize vectors holding DMs and Vecs
        wDMs = std::vector<DM>(1);
        vecs = std::vector<Vec>(1);

        // create DMs
        ierr = DMDACreate2d(mesh->comm,
                DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                n[0][0], n[0][1], mesh->nProc[0], mesh->nProc[1],
                1, 1, nullptr, nullptr, &wDMs[0]); CHKERRQ(ierr);
        
        ierr = DMSetUp(wDMs[0]); CHKERRQ(ierr);

        // create Vecs
        ierr = DMCreateGlobalVector(wDMs[0], &vecs[0]); CHKERRQ(ierr);

        // set names
        ierr = PetscObjectSetName(PetscObject(vecs[0]), "wz"); CHKERRQ(ierr);
    }
    else
    {
        // initialize size of meshes
        n = petibm::type::IntVec2D(3, petibm::type::IntVec1D(3));

        // set values of n for wx
        n[0][0] = mesh->n[3][0];
        n[0][1] = mesh->n[4][1];
        n[0][2] = mesh->n[4][2];
        // set values of n for wy
        n[1][0] = mesh->n[4][0];
        n[1][1] = mesh->n[3][1];
        n[1][2] = mesh->n[4][2];
        // set values of n for wz
        n[2][0] = mesh->n[4][0];
        n[2][1] = mesh->n[4][1];
        n[2][2] = mesh->n[3][2];

        // initialize vectors holding coordinates
        coord = petibm::type::RealVec3D(3, petibm::type::RealVec2D(3));

        // set coordinates for wx
        coord[0][0].assign(mesh->coord[3][0], mesh->coord[3][0]+n[0][0]);
        coord[0][1].assign(mesh->coord[4][1], mesh->coord[4][1]+n[0][1]);
        coord[0][2].assign(mesh->coord[4][2], mesh->coord[4][2]+n[0][2]);
        // set coordinates for wy
        coord[1][0].assign(mesh->coord[4][0], mesh->coord[4][0]+n[1][0]);
        coord[1][1].assign(mesh->coord[3][1], mesh->coord[3][1]+n[1][1]);
        coord[1][2].assign(mesh->coord[4][2], mesh->coord[4][2]+n[1][2]);
        // set coordinates for wz
        coord[2][0].assign(mesh->coord[4][0], mesh->coord[4][0]+n[2][0]);
        coord[2][1].assign(mesh->coord[4][1], mesh->coord[4][1]+n[2][1]);
        coord[2][2].assign(mesh->coord[3][2], mesh->coord[3][2]+n[2][2]);

        // name of vorticity fields
        names = std::vector<std::string>(3);

        // set name for wx
        names[0] = "wx";
        // set name for wy
        names[1] = "wy";
        // set name for wz
        names[2] = "wz";

        wDMs = std::vector<DM>(3);
        vecs = std::vector<Vec>(3);

        for(PetscInt f=0; f<3; ++f)
        {
            // create DMs
            ierr = DMDACreate3d(mesh->comm,
                    DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                    DMDA_STENCIL_BOX,
                    n[f][0], n[f][1], n[f][2], 
                    mesh->nProc[0], mesh->nProc[1], mesh->nProc[2],
                    1, 1, nullptr, nullptr, nullptr, &wDMs[f]); CHKERRQ(ierr);
        
            ierr = DMSetUp(wDMs[f]); CHKERRQ(ierr);

            // create Vecs
            ierr = DMCreateGlobalVector(wDMs[f], &vecs[f]); CHKERRQ(ierr);

            // set name
            ierr = PetscObjectSetName(
                    PetscObject(vecs[f]), names[f].c_str()); CHKERRQ(ierr);
        }
    }
    PetscFunctionReturn(0);
} // initVorticityMesh

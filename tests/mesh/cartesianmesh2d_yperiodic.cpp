/**
 * \file cartesianmesh_test.cpp
 * \brief Unit-tests for the class `CartesianMesh`.
 */

#include <petsc.h>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <petibm/mesh.h>


class CartesainMeshTest2D_YPeriodic : public ::testing::Test
{
protected:

    CartesainMeshTest2D_YPeriodic() {};

    virtual ~CartesainMeshTest2D_YPeriodic() {};

    static void SetUpTestCase()
    {
        using namespace YAML;
        
        Node config;
        
        config["mesh"].push_back(Node(NodeType::Map));
        
        config["mesh"][0]["direction"] = "x";
        config["mesh"][0]["start"] = "0.1";
        config["mesh"][0]["subDomains"].push_back(Node(NodeType::Map));
        config["mesh"][0]["subDomains"][0]["end"] = 1.6;
        config["mesh"][0]["subDomains"][0]["cells"] = 4;
        config["mesh"][0]["subDomains"][0]["stretchRatio"] = 0.5;
        config["mesh"][0]["subDomains"][1]["end"] = 1.9;
        config["mesh"][0]["subDomains"][1]["cells"] = 3;
        config["mesh"][0]["subDomains"][1]["stretchRatio"] = 1;
        config["mesh"][0]["subDomains"][2]["end"] = 5.0;
        config["mesh"][0]["subDomains"][2]["cells"] = 5;
        config["mesh"][0]["subDomains"][2]["stretchRatio"] = 2.0;
        
        config["mesh"][1]["direction"] = "y";
        config["mesh"][1]["start"] = "0.05";
        config["mesh"][1]["subDomains"].push_back(Node(NodeType::Map));
        config["mesh"][1]["subDomains"][0]["end"] = 0.8625;
        config["mesh"][1]["subDomains"][0]["cells"] = 4;
        config["mesh"][1]["subDomains"][0]["stretchRatio"] = 0.666666666666666;
        config["mesh"][1]["subDomains"][1]["end"] = 1.1625;
        config["mesh"][1]["subDomains"][1]["cells"] = 3;
        config["mesh"][1]["subDomains"][1]["stretchRatio"] = 1.0;
        config["mesh"][1]["subDomains"][2]["end"] = 1.975;
        config["mesh"][1]["subDomains"][2]["cells"] = 4;
        config["mesh"][1]["subDomains"][2]["stretchRatio"] = 1.5;
        
        config["flow"] = YAML::Node(NodeType::Map);
        config["flow"]["boundaryConditions"].push_back(Node(NodeType::Map));
        config["flow"]["boundaryConditions"][0]["location"] = "xMinus";
        config["flow"]["boundaryConditions"][1]["location"] = "xPlus";
        config["flow"]["boundaryConditions"][2]["location"] = "yMinus";
        config["flow"]["boundaryConditions"][3]["location"] = "yPlus";
        
        for(unsigned int i=0; i<2; ++i)
        {
            config["flow"]["boundaryConditions"][i]["u"][0] = "DIRICHLET";
            config["flow"]["boundaryConditions"][i]["u"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["v"][0] = "DIRICHLET";
            config["flow"]["boundaryConditions"][i]["v"][1] = 0.0;
        }
        
        for(unsigned int i=2; i<4; ++i)
        {
            config["flow"]["boundaryConditions"][i]["u"][0] = "PERIODIC";
            config["flow"]["boundaryConditions"][i]["u"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["v"][0] = "PERIODIC";
            config["flow"]["boundaryConditions"][i]["v"][1] = 0.0;
        }
        
        config["flow"]["boundaryConditions"][3]["u"][1] = 1.0;
        
        
        PetscErrorCode ierr;
        PetscInt nProc[3];
        
        ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh); ASSERT_FALSE(ierr);
        
        ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, 
                DMDA_STENCIL_BOX, 12, 11, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
                nullptr, nullptr, &da[3]); ASSERT_FALSE(ierr);
        ierr = DMSetUp(da[3]); ASSERT_FALSE(ierr);
        
        ierr = DMDAGetInfo(da[3], nullptr, nullptr, nullptr, nullptr, &nProc[0], 
                &nProc[1], &nProc[2], nullptr, nullptr, nullptr, nullptr, 
                nullptr, nullptr); ASSERT_FALSE(ierr);
        
        ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, 
                DMDA_STENCIL_BOX, 11, 11, nProc[0], nProc[1], 1, 1,
                nullptr, nullptr, &da[0]); ASSERT_FALSE(ierr);
        ierr = DMSetUp(da[0]); ASSERT_FALSE(ierr);
        
        ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, 
                DMDA_STENCIL_BOX, 12, 11, nProc[0], nProc[1], 1, 1,
                nullptr, nullptr, &da[1]); ASSERT_FALSE(ierr);
        ierr = DMSetUp(da[1]); ASSERT_FALSE(ierr);
    };

    virtual void SetUp() {};

    virtual void TearDown() {};

    static void TearDownTestCase()
    { 
        mesh.~shared_ptr(); 
        
        PetscErrorCode ierr;
        ierr = DMDestroy(&da[0]); ASSERT_FALSE(ierr);
        ierr = DMDestroy(&da[1]); ASSERT_FALSE(ierr);
        ierr = DMDestroy(&da[3]); ASSERT_FALSE(ierr);
    };

    static petibm::type::Mesh mesh;
    static DM da[5];
}; // CartesianMeshTest

petibm::type::Mesh CartesainMeshTest2D_YPeriodic::mesh = nullptr;
DM CartesainMeshTest2D_YPeriodic::da[5] {nullptr, nullptr, nullptr, nullptr, nullptr};

// test dimension
TEST_F(CartesainMeshTest2D_YPeriodic, check_dim)
{
    ASSERT_EQ(2, mesh->dim);
}

// check min
TEST_F(CartesainMeshTest2D_YPeriodic, check_min)
{
    ASSERT_DOUBLE_EQ(0.1, mesh->min[0]);
    ASSERT_DOUBLE_EQ(0.05, mesh->min[1]);
    ASSERT_DOUBLE_EQ(0.0, mesh->min[2]);
}

// check max
TEST_F(CartesainMeshTest2D_YPeriodic, check_max)
{
    ASSERT_DOUBLE_EQ(5.0, mesh->max[0]);
    ASSERT_DOUBLE_EQ(1.975, mesh->max[1]);
    ASSERT_DOUBLE_EQ(1.0, mesh->max[2]);
}

// check n
TEST_F(CartesainMeshTest2D_YPeriodic, check_n)
{
    petibm::type::RealVec2D n {
        {11, 11, 1}, {12, 11, 1}, {1, 1, 1}, {12, 11, 1}, {13, 12, 1}};
    
    for(unsigned int f=0; f<5; ++f)
        for(unsigned int d=0; d<3; ++d)
            ASSERT_EQ(n[f][d], mesh->n[f][d]);
}

// check periodic
TEST_F(CartesainMeshTest2D_YPeriodic, check_periodic)
{
    petibm::type::BoolVec2D periodic{{PETSC_FALSE, PETSC_TRUE, PETSC_FALSE},
        {PETSC_FALSE, PETSC_TRUE, PETSC_FALSE}, {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE}};
    
    for(unsigned int f=0; f<3; ++f)
        for(unsigned int d=0; d<3; ++d)
            ASSERT_EQ(periodic[f][d], mesh->periodic[f][d]);
}

// check coord
TEST_F(CartesainMeshTest2D_YPeriodic, check_coord)
{
    petibm::type::RealVec3D coordTrue {
            {{0.1, 0.9, 1.3, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.6, 3.4, 5.0},
             {-0.11875, 0.21875, 0.5, 0.6875, 0.8125, 0.9125, 1.0125, 1.1125, 1.2125, 1.3375, 1.525, 1.80625, 2.14375},
             {0.0}},
            {{-0.3, 0.5, 1.1, 1.4, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.4, 3.0, 4.2, 5.8},
             {0.05, 0.3875, 0.6125, 0.7625, 0.8625, 0.9625, 1.0625, 1.1625, 1.2625, 1.4125, 1.6375, 1.975, 2.3125},
             {0.0}},
            {{0.0},
             {0.0},
             {0.0}},
            {{0.5, 1.1, 1.4, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.4, 3.0, 4.2},
             {0.21875, 0.5, 0.6875, 0.8125, 0.9125, 1.0125, 1.1125, 1.2125, 1.3375, 1.525, 1.80625},
             {0.0}},
            {{0.1, 0.9, 1.3, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.6, 3.4, 5.0},
             {0.05, 0.3875, 0.6125, 0.7625, 0.8625, 0.9625, 1.0625, 1.1625, 1.2625, 1.4125, 1.6375, 1.975},
             {0.0}}};
    
    petibm::type::GhostedVec3D coord {
        {&coordTrue[0][0][1], &coordTrue[0][1][1], &coordTrue[0][2][0]},
        {&coordTrue[1][0][1], &coordTrue[1][1][1], &coordTrue[1][2][0]},
        {&coordTrue[2][0][0], &coordTrue[2][1][0], &coordTrue[2][2][0]},
        {&coordTrue[3][0][0], &coordTrue[3][1][0], &coordTrue[3][2][0]},
        {&coordTrue[4][0][0], &coordTrue[4][1][0], &coordTrue[4][2][0]},
        };
    
    for(unsigned int f=0; f<2; ++f)
    {
        for(unsigned int d=0; d<2; ++d)
            for(int i=-1; i<mesh->n[f][d]+1; ++i)
                ASSERT_NEAR(coord[f][d][i], mesh->coord[f][d][i], 1e-12)
                    << "dL incorrect at field " << f
                    << ", direction " << d
                    << ", index " << i << std::endl;
        
        ASSERT_NEAR(coord[f][2][0], mesh->coord[f][2][0], 1e-12)
            << "dL incorrect at field " << f
            << ", direction " << 2
            << ", index " << 0 << std::endl;
    }
    
    for(unsigned int f=2; f<5; ++f)
        for(unsigned int d=0; d<3; ++d)
            for(int i=0; i<mesh->n[f][d]; ++i)
                ASSERT_NEAR(coord[f][d][i], mesh->coord[f][d][i], 1e-12)
                    << "dL incorrect at field " << f
                    << ", direction " << d
                    << ", index " << i << std::endl;
}

// check dL
TEST_F(CartesainMeshTest2D_YPeriodic, check_dL)
{
    petibm::type::RealVec3D dLTrue {
            {{0.8, 0.6, 0.3, 0.15, 0.1, 0.1, 0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 1.6},
             {0.3375, 0.3375, 0.225, 0.15, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.225, 0.3375, 0.3375},
             {1.0}},
            {{0.8, 0.8, 0.4, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.4, 0.8, 1.6, 1.6},
             {0.3375, 0.28125, 0.1875, 0.125, 0.1, 0.1, 0.1, 0.1, 0.125, 0.1875, 0.28125, 0.3375, 0.28125},
             {1.0}},
            {{1.0},
             {1.0},
             {1.0}},
            {{0.8, 0.4, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.4, 0.8, 1.6},
             {0.3375, 0.225, 0.15, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.225, 0.3375},
             {1.0}},
            {{1.0},
             {1.0},
             {1.0}}};
    
    petibm::type::GhostedVec3D dL {
        {&dLTrue[0][0][1], &dLTrue[0][1][1], &dLTrue[0][2][0]},
        {&dLTrue[1][0][1], &dLTrue[1][1][1], &dLTrue[1][2][0]},
        {&dLTrue[2][0][0], &dLTrue[2][1][0], &dLTrue[2][2][0]},
        {&dLTrue[3][0][0], &dLTrue[3][1][0], &dLTrue[3][2][0]},
        {nullptr, nullptr, nullptr}};
    
    for(unsigned int f=0; f<2; ++f)
    {
        for(unsigned int d=0; d<2; ++d)
            for(int i=-1; i<mesh->n[f][d]+1; ++i)
                ASSERT_NEAR(dL[f][d][i], mesh->dL[f][d][i], 1e-12)
                    << "dL incorrect at field " << f
                    << ", direction " << d
                    << ", index " << i << std::endl;
        
        ASSERT_NEAR(dL[f][2][0], mesh->dL[f][2][0], 1e-12)
                    << "dL incorrect at field " << f
                    << ", direction " << 2
                    << ", index " << 0 << std::endl;
    }
    
    for(unsigned int f=2; f<4; ++f)
        for(unsigned int d=0; d<3; ++d)
            for(int i=0; i<mesh->n[f][d]; ++i)
                ASSERT_NEAR(dL[f][d][i], mesh->dL[f][d][i], 1e-12)
                    << "dL incorrect at field " << f
                    << ", direction " << d
                    << ", index " << i << std::endl;
    
    for(unsigned int d=0; d<3; ++d)
        ASSERT_FALSE(mesh->dL[4][d])
            << "The dL of the vertex mesh is not nullptr." << std::endl;
}

// check UN
TEST_F(CartesainMeshTest2D_YPeriodic, check_UN)
{
    ASSERT_EQ(253, mesh->UN)
        << "The total number of velocity points is not correct." << std::endl;
}

// check pN
TEST_F(CartesainMeshTest2D_YPeriodic, check_pN)
{
    ASSERT_EQ(132, mesh->pN) 
        << "The total number of pressure points is not correct." << std::endl;
}

// check da
TEST_F(CartesainMeshTest2D_YPeriodic, check_da)
{
    PetscErrorCode ierr;
    DMType  type;
    
    ierr = DMGetType(mesh->da[0], &type); ASSERT_FALSE(ierr);
    ASSERT_STREQ(type, DMDA) << "u-velocity mesh is " << type << ", not expected da" << std::endl;
    
    ierr = DMGetType(mesh->da[1], &type); ASSERT_FALSE(ierr);
    ASSERT_STREQ(type, DMDA) << "v-velocity mesh is " << type << ", not expected da" << std::endl;
    
    ierr = DMGetType(mesh->da[3], &type); ASSERT_FALSE(ierr);
    ASSERT_STREQ(type, DMDA) << "v-velocity mesh is " << type << ", not expected da" << std::endl;
    
    ASSERT_FALSE(mesh->da[2]) << "w-velocity mesh should not be initialized." << std::endl;
    ASSERT_FALSE(mesh->da[4]) << "vertex mesh should not be initialized." << std::endl;
}

// check nProc
TEST_F(CartesainMeshTest2D_YPeriodic, check_nProc)
{
    PetscErrorCode ierr;
    
    PetscInt nProc[3];
    
    ierr = DMDAGetInfo(da[3], nullptr, nullptr, nullptr, nullptr, &nProc[0], 
            &nProc[1], &nProc[2], nullptr, nullptr, nullptr, nullptr, 
            nullptr, nullptr); ASSERT_FALSE(ierr);
    
    ASSERT_EQ(nProc[0], mesh->nProc[0]);
    ASSERT_EQ(nProc[1], mesh->nProc[1]);
    ASSERT_EQ(nProc[2], mesh->nProc[2]);
}

// check bg
TEST_F(CartesainMeshTest2D_YPeriodic, check_bg)
{
    PetscErrorCode ierr;
    for(unsigned int f=0; f<5; f++)
    {
        if ((f == 2) || (f == 4)) continue;
        DMDALocalInfo info;
        ierr = DMDAGetLocalInfo(da[f], &info); ASSERT_FALSE(ierr); 
        ASSERT_EQ(info.xs, mesh->bg[f][0]) << "f = " << f << std::endl;
        ASSERT_EQ(info.ys, mesh->bg[f][1]) << "f = " << f << std::endl;
        ASSERT_EQ(info.zs, mesh->bg[f][2]) << "f = " << f << std::endl;
    }
    
    ASSERT_EQ(0, mesh->bg[2][0]) << "f = " << 2 << std::endl;
    ASSERT_EQ(0, mesh->bg[2][1]) << "f = " << 2 << std::endl;
    ASSERT_EQ(0, mesh->bg[2][2]) << "f = " << 2 << std::endl;
}

// check ed
TEST_F(CartesainMeshTest2D_YPeriodic, check_ed)
{
    PetscErrorCode ierr;
    for(unsigned int f=0; f<5; f++)
    {
        if ((f == 2) || (f == 4)) continue;
        DMDALocalInfo info;
        ierr = DMDAGetLocalInfo(da[f], &info); ASSERT_FALSE(ierr); 
        ASSERT_EQ(info.xs+info.xm, mesh->ed[f][0]) << "f = " << f << std::endl;
        ASSERT_EQ(info.ys+info.ym, mesh->ed[f][1]) << "f = " << f << std::endl;
        ASSERT_EQ(info.zs+info.zm, mesh->ed[f][2]) << "f = " << f << std::endl;
    }
    
    ASSERT_EQ(1, mesh->ed[2][0]) << "f = " << 2 << std::endl;
    ASSERT_EQ(1, mesh->ed[2][1]) << "f = " << 2 << std::endl;
    ASSERT_EQ(1, mesh->ed[2][2]) << "f = " << 2 << std::endl;
}

// check m
TEST_F(CartesainMeshTest2D_YPeriodic, check_m)
{
    PetscErrorCode ierr;
    for(unsigned int f=0; f<5; f++)
    {
        if ((f == 2) || (f == 4)) continue;
        DMDALocalInfo info;
        ierr = DMDAGetLocalInfo(da[f], &info); ASSERT_FALSE(ierr); 
        ASSERT_EQ(info.xm, mesh->m[f][0]) << "f = " << f << std::endl;
        ASSERT_EQ(info.ym, mesh->m[f][1]) << "f = " << f << std::endl;
        ASSERT_EQ(info.zm, mesh->m[f][2]) << "f = " << f << std::endl;
    }
    
    ASSERT_EQ(0, mesh->m[2][0]) << "f = " << 2 << std::endl;
    ASSERT_EQ(0, mesh->m[2][1]) << "f = " << 2 << std::endl;
    ASSERT_EQ(0, mesh->m[2][2]) << "f = " << 2 << std::endl;
}

// check UNLocal
TEST_F(CartesainMeshTest2D_YPeriodic, check_UNLocal)
{
    int N = mesh->m[0][0] * mesh->m[0][1] + 
        mesh->m[1][0] * mesh->m[1][1];
    ASSERT_EQ(N, mesh->UNLocal);
}

// check pNLocal
TEST_F(CartesainMeshTest2D_YPeriodic, check_pNLocal)
{
    int N = mesh->m[3][0] * mesh->m[3][1]; 
    ASSERT_EQ(N, mesh->pNLocal);
}

// check UPack
TEST_F(CartesainMeshTest2D_YPeriodic, check_UPack)
{
    PetscErrorCode ierr;
    DMType  type;
    
    ierr = DMGetType(mesh->UPack, &type); ASSERT_FALSE(ierr);
    ASSERT_STREQ(type, DMCOMPOSITE);
    
    DM dms[3] {nullptr, nullptr, nullptr};
    ierr = DMCompositeGetEntriesArray(mesh->UPack, dms); ASSERT_FALSE(ierr);
    ASSERT_EQ(dms[0], mesh->da[0]);
    ASSERT_EQ(dms[1], mesh->da[1]);
    ASSERT_EQ(dms[2], mesh->da[2]);
    
    ierr = DMDestroy(&dms[0]); ASSERT_FALSE(ierr);
    ierr = DMDestroy(&dms[1]); ASSERT_FALSE(ierr);
}

// check mpi
TEST_F(CartesainMeshTest2D_YPeriodic, check_mpi)
{
    PetscErrorCode  ierr;
    PetscMPIInt     rank, size;
    
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
    ASSERT_FALSE(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    ASSERT_FALSE(ierr);
    
    ASSERT_EQ(PETSC_COMM_WORLD, mesh->comm);
    ASSERT_EQ(size, mesh->mpiSize);
    ASSERT_EQ(rank, mesh->mpiRank);
}

// check getLocalIndex with i, j, k
TEST_F(CartesainMeshTest2D_YPeriodic, check_getLocalIndex_ijk)
{
    for(int f=0; f<5; ++f)
    {
        switch (f)
        {
            case 2: case 4: break;
            default:
                for(int k=mesh->bg[f][2]; k<mesh->ed[f][2]; k++)
                    for(int j=mesh->bg[f][1]-1; j<mesh->ed[f][1]+1; j++)
                        for(int i=mesh->bg[f][0]-1; i<mesh->ed[f][0]+1; i++)
                        {
                            PetscInt idx;
                            mesh->getLocalIndex(f, i, j, k, idx);
                            
                            PetscInt i2, j2, k2, idx2;
                            
                            i2 = i - mesh->bg[f][0] + 1;
                            j2 = j - mesh->bg[f][1] + 1;
                            k2 = k - mesh->bg[f][2];
                            
                            idx2 = i2 + j2 * (mesh->m[f][0] + 2) + 
                                k2 * ((mesh->m[f][0] + 2) * (mesh->m[f][1] + 2));
                            
                            ASSERT_EQ(idx2, idx)
                                << "getLocalIndex  wrong at field " << f <<" and (" 
                                << i << ", " << j << ", " << k << ") " << std::endl;
                        }
        }
    }
}

// check getNaturalIndex with i, j, k
TEST_F(CartesainMeshTest2D_YPeriodic, check_getNaturalIndex_ijk)
{
    for(int f=0; f<5; ++f)
    {
        // internal points
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int j=0; j<mesh->n[f][1]; j++)
                for(int i=0; i<mesh->n[f][0]; i++)
                {
                    PetscInt idx, idx2;
                    mesh->getNaturalIndex(f, i, j, k, idx);
                    
                    idx2 = i + j * mesh->n[f][0] + 
                        k * (mesh->n[f][0] * mesh->n[f][1]);
                    
                    ASSERT_EQ(idx2, idx)
                        << "getNatrualIndex wrong at field " << f 
                        << " and internal points (" 
                        << i << ", " << j << ", " << k << ") " << std::endl;
                }
        
        // ghost points at top and bottom
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int i=0; i<mesh->n[f][0]; i++)
            {
                PetscInt idx, idx2, j;
                
                j = mesh->n[f][1];
                mesh->getNaturalIndex(f, i, j, k, idx);
                idx2 = i + k * (mesh->n[f][0] * mesh->n[f][1]);
                ASSERT_EQ(idx2, idx)
                    << "getNatrualIndex wrong at field " << f 
                    << " and top ghost point ("
                    << i << ", " << j << ", " << k << ") " << std::endl;
                
                j = -1;
                mesh->getNaturalIndex(f, i, j, k, idx);
                idx2 = i + (mesh->n[f][1] - 1) * mesh->n[f][0] + 
                    k * (mesh->n[f][0] * mesh->n[f][1]);
                ASSERT_EQ(idx2, idx)
                    << "getNatrualIndex  wrong at field " << f 
                    << " and bottom ghost point ("
                    << i << ", " << j << ", " << k << ") " << std::endl;
            }
        
        // ghost points at right and left
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int j=0; j<mesh->n[f][1]; j++)
            {
                PetscInt idx, i;
                
                i = mesh->n[f][0];
                mesh->getNaturalIndex(f, i, j, k, idx);
                ASSERT_EQ(-1, idx)
                    << "getNatrualIndex wrong at field " << f  
                    << " and right ghost point ("
                    << i << ", " << j << ", " << k << ") " << std::endl;
                
                i = -1;
                mesh->getNaturalIndex(f, i, j, k, idx);
                ASSERT_EQ(-1, idx)
                    << "getNatrualIndex wrong at field " << f 
                    << " and left ghost point ("
                    << i << ", " << j << ", " << k << ") " << std::endl;
            }
        
        // corners
        for(int k=0; k<mesh->n[f][2]; k++)
        {
            PetscInt idx;
            
            mesh->getNaturalIndex(f, -1, -1, 0, idx);
            ASSERT_EQ(-1, idx)
                << "getNatrualIndex wrong at field " << f 
                << " and corner ghost point ("
                << -1 << ", " << -1 << ", " << k << ") " << std::endl;
            
            mesh->getNaturalIndex(f, -1, mesh->n[f][1], 0, idx);
            ASSERT_EQ(-1, idx)
                << "getNatrualIndex wrong at field " << f
                << " and corner ghost point ("
                << -1 << ", " << mesh->n[f][1] << ", " << k << ") " << std::endl;
            
            mesh->getNaturalIndex(f, mesh->n[f][0], -1, 0, idx);
            ASSERT_EQ(-1, idx)
                << "getNatrualIndex wrong at field " << f
                << " and corner ghost point ("
                << mesh->n[f][0] << ", " << -1 << ", " << k << ") " << std::endl;
            
            mesh->getNaturalIndex(f, mesh->n[f][0], mesh->n[f][1], 0, idx);
            ASSERT_EQ(-1, idx)
                << "getNatrualIndex wrong at field " << f
                << " and corner ghost point ("
                << mesh->n[f][0] << ", " << mesh->n[f][1] << ", " << k << ") " << std::endl;
        }
    }
}

// check getGlobalIndex with i, j, k of internal points
TEST_F(CartesainMeshTest2D_YPeriodic, check_getGlobalIndex_ijk_internal)
{
    for(int f=0; f<5; ++f)
    {
        if ((f == 2) || (f == 4)) continue;
        
        PetscErrorCode  ierr;
        
        ISLocalToGlobalMapping mapping;
        ierr = DMGetLocalToGlobalMapping(mesh->da[f], &mapping); ASSERT_FALSE(ierr);
        
        Vec vini, vinj, vink;
        ierr = DMGetGlobalVector(mesh->da[f], &vini); ASSERT_FALSE(ierr);
        ierr = DMGetGlobalVector(mesh->da[f], &vinj); ASSERT_FALSE(ierr);
        ierr = DMGetGlobalVector(mesh->da[f], &vink); ASSERT_FALSE(ierr);
        
        for(int k=mesh->bg[f][2]; k<mesh->ed[f][2]; k++)
            for(int j=mesh->bg[f][1]; j<mesh->ed[f][1]; j++)
                for(int i=mesh->bg[f][0]; i<mesh->ed[f][0]; i++)
                {
                    PetscInt idx;
                    
                    ierr = DMDAConvertToCell(mesh->da[f], {k, j, i, 0}, &idx); ASSERT_FALSE(ierr);
                    ierr = ISLocalToGlobalMappingApply(mapping, 1, &idx, &idx); ASSERT_FALSE(ierr);
                    
                    ierr = VecSetValue(vini, idx, i, INSERT_VALUES); ASSERT_FALSE(ierr);
                    ierr = VecSetValue(vinj, idx, j, INSERT_VALUES); ASSERT_FALSE(ierr);
                    ierr = VecSetValue(vink, idx, k, INSERT_VALUES); ASSERT_FALSE(ierr);
                }
        ierr = VecAssemblyBegin(vini); ASSERT_FALSE(ierr);
        ierr = VecAssemblyBegin(vinj); ASSERT_FALSE(ierr);
        ierr = VecAssemblyBegin(vink); ASSERT_FALSE(ierr);
        ierr = VecAssemblyEnd(vini); ASSERT_FALSE(ierr);
        ierr = VecAssemblyEnd(vinj); ASSERT_FALSE(ierr);
        ierr = VecAssemblyEnd(vink); ASSERT_FALSE(ierr);
        
        Vec vouti, voutj, voutk;
        VecScatter ctx;
        ierr = VecScatterCreateToAll(vini, &ctx, &vouti); ASSERT_FALSE(ierr);
        ierr = VecDuplicate(vouti, &voutj); ASSERT_FALSE(ierr);
        ierr = VecDuplicate(vouti, &voutk); ASSERT_FALSE(ierr);
        ierr = VecScatterBegin(ctx, vini, vouti, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
        ierr = VecScatterEnd(ctx, vini, vouti, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
        ierr = VecScatterBegin(ctx, vinj, voutj, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
        ierr = VecScatterEnd(ctx, vinj, voutj, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
        ierr = VecScatterBegin(ctx, vink, voutk, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
        ierr = VecScatterEnd(ctx, vink, voutk, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);

        // internal points
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int j=0; j<mesh->n[f][1]; j++)
                for(int i=0; i<mesh->n[f][0]; i++)
                {
                    PetscInt idx;
                    ierr = mesh->getGlobalIndex(f, i, j, k, idx); ASSERT_FALSE(ierr);
                    
                    PetscReal vi, vj, vk;
                    ierr = VecGetValues(vouti, 1, &idx, &vi); ASSERT_FALSE(ierr);
                    ierr = VecGetValues(voutj, 1, &idx, &vj); ASSERT_FALSE(ierr);
                    ierr = VecGetValues(voutk, 1, &idx, &vk); ASSERT_FALSE(ierr);
                    
                    ASSERT_EQ(vi, PetscReal(i));
                    ASSERT_EQ(vj, PetscReal(j));
                    ASSERT_EQ(vk, PetscReal(k));
                }
        
        ierr = VecScatterDestroy(&ctx); ASSERT_FALSE(ierr);
        ierr = VecDestroy(&voutk); ASSERT_FALSE(ierr);
        ierr = VecDestroy(&voutj); ASSERT_FALSE(ierr);
        ierr = VecDestroy(&vouti); ASSERT_FALSE(ierr);
        ierr = VecDestroy(&vink); ASSERT_FALSE(ierr);
        ierr = VecDestroy(&vinj); ASSERT_FALSE(ierr);
        ierr = VecDestroy(&vini); ASSERT_FALSE(ierr);
        
        //ierr = ISLocalToGlobalMappingDestroy(&mapping); ASSERT_FALSE(ierr);
        // Seem there's a bug in ISLocalToGlobalMappingDestroy. It should just
        // decrease the reference number by one. Instead, it now destroies the 
        // underlying object completely....
        mapping = PETSC_NULL;
    }
}

// check getGlobalIndex with i, j, k of ghost points
TEST_F(CartesainMeshTest2D_YPeriodic, check_getGlobalIndex_ijk_ghosts)
{
    for(int f=0; f<5; ++f)
    {
        if ((f == 2) || (f == 4)) continue;
        
        // ghost points at top and bottom
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int i=0; i<mesh->n[f][0]; i++)
            {
                PetscInt idx, idx2, j;
                
                j = mesh->n[f][1];
                mesh->getGlobalIndex(f, i, j, k, idx);
                
                // internal part should have already passed
                mesh->getGlobalIndex(f, i, 0, k, idx2);
                ASSERT_EQ(idx2, idx)
                    << "getGlobalIndex wrong at field " << f <<" and (" 
                    << i << ", " << j << ", " << k << ") " << std::endl;
                
                j = -1;
                mesh->getGlobalIndex(f, i, j, k, idx);
                
                // internal part should have already passed
                mesh->getGlobalIndex(f, i, mesh->n[f][1]-1, k, idx2);
                ASSERT_EQ(idx2, idx)
                    << "getGlobalIndex wrong at field " << f <<" and (" 
                    << i << ", " << j << ", " << k << ") " << std::endl;
            }
        
        // ghost points at right and left
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int j=0; j<mesh->n[f][1]; j++)
            {
                PetscInt idx, i;
                
                i = mesh->n[f][0];
                mesh->getGlobalIndex(f, i, j, k, idx);
                ASSERT_EQ(-1, idx)
                    << "getGlobalIndex wrong at field " << f <<" and (" 
                    << i << ", " << j << ", " << k << ") " << std::endl;
                
                i = -1;
                mesh->getGlobalIndex(f, i, j, k, idx);
                ASSERT_EQ(-1, idx)
                    << "getGlobalIndex wrong at field " << f <<" and (" 
                    << i << ", " << j << ", " << k << ") " << std::endl;
            }
        
        // corners
        for(int k=0; k<mesh->n[f][2]; k++)
        {
            PetscInt idx;
            
            mesh->getGlobalIndex(f, -1, -1, 0, idx);
            ASSERT_EQ(-1, idx)
                << "getGlobalIndex wrong at field " << f <<" and (" 
                << -1 << ", " << -1 << ", " << k << ") " << std::endl;
            
            mesh->getGlobalIndex(f, -1, mesh->n[f][1], 0, idx);
            ASSERT_EQ(-1, idx)
                << "getGlobalIndex wrong at field " << f <<" and (" 
                << -1 << ", " << mesh->n[f][1] << ", " << k << ") " << std::endl;
            
            mesh->getGlobalIndex(f, mesh->n[f][0], -1, 0, idx);
            ASSERT_EQ(-1, idx)
                << "getGlobalIndex wrong at field " << f <<" and (" 
                << mesh->n[f][0] << ", " << -1 << ", " << k << ") " << std::endl;
            
            mesh->getGlobalIndex(f, mesh->n[f][0], mesh->n[f][1], 0, idx);
            ASSERT_EQ(-1, idx)
                << "getGlobalIndex wrong at field " << f <<" and (" 
                << mesh->n[f][0] << ", " << mesh->n[f][1] << ", " << k << ") " << std::endl;
        }
    }
}

// check getPackedGlobalIndex with i, j, k of internal points
TEST_F(CartesainMeshTest2D_YPeriodic, check_getPackedGlobalIndex_ijk_internal)
{
    PetscErrorCode  ierr;
    
    DM UPack;
    
    ierr = DMCompositeCreate(PETSC_COMM_WORLD, &UPack); ASSERT_FALSE(ierr);
    ierr = DMCompositeAddDM(UPack, da[0]); ASSERT_FALSE(ierr);
    ierr = DMCompositeAddDM(UPack, da[1]); ASSERT_FALSE(ierr);
    ierr = DMSetUp(UPack); ASSERT_FALSE(ierr);

    Vec fp, ip, jp, kp;
    ierr = DMCreateGlobalVector(mesh->UPack, &fp); ASSERT_FALSE(ierr);
    ierr = DMCreateGlobalVector(mesh->UPack, &ip); ASSERT_FALSE(ierr);
    ierr = DMCreateGlobalVector(mesh->UPack, &jp); ASSERT_FALSE(ierr);
    ierr = DMCreateGlobalVector(mesh->UPack, &kp); ASSERT_FALSE(ierr);
    
    Vec fv[2], iv[2], jv[2], kv[2];
    ierr = DMCompositeGetAccessArray(mesh->UPack, fp, 2, nullptr, fv); ASSERT_FALSE(ierr);
    ierr = DMCompositeGetAccessArray(mesh->UPack, ip, 2, nullptr, iv); ASSERT_FALSE(ierr);
    ierr = DMCompositeGetAccessArray(mesh->UPack, jp, 2, nullptr, jv); ASSERT_FALSE(ierr);
    ierr = DMCompositeGetAccessArray(mesh->UPack, kp, 2, nullptr, kv); ASSERT_FALSE(ierr);
    
    for(int f=0; f<2; ++f)
    {
        ISLocalToGlobalMapping mapping;
        ierr = DMGetLocalToGlobalMapping(mesh->da[f], &mapping); ASSERT_FALSE(ierr);
        
        for(int k=mesh->bg[f][2]; k<mesh->ed[f][2]; ++k)
            for(int j=mesh->bg[f][1]; j<mesh->ed[f][1]; ++j)
                for(int i=mesh->bg[f][0]; i<mesh->ed[f][0]; ++i)
                {
                    PetscInt    idx;
                    ierr = DMDAConvertToCell(mesh->da[f], {k, j, i, 0}, &idx); ASSERT_FALSE(ierr);
                    ierr = ISLocalToGlobalMappingApply(mapping, 1, &idx, &idx); ASSERT_FALSE(ierr);
                    
                    ierr = VecSetValue(fv[f], idx, (PetscReal) f, INSERT_VALUES); ASSERT_FALSE(ierr);
                    ierr = VecSetValue(iv[f], idx, (PetscReal) i, INSERT_VALUES); ASSERT_FALSE(ierr);
                    ierr = VecSetValue(jv[f], idx, (PetscReal) j, INSERT_VALUES); ASSERT_FALSE(ierr);
                    ierr = VecSetValue(kv[f], idx, (PetscReal) k, INSERT_VALUES); ASSERT_FALSE(ierr);
                }
        
        ierr = VecAssemblyBegin(fv[f]); ASSERT_FALSE(ierr);
        ierr = VecAssemblyBegin(iv[f]); ASSERT_FALSE(ierr);
        ierr = VecAssemblyBegin(jv[f]); ASSERT_FALSE(ierr);
        ierr = VecAssemblyBegin(kv[f]); ASSERT_FALSE(ierr);
        ierr = VecAssemblyEnd(fv[f]); ASSERT_FALSE(ierr);
        ierr = VecAssemblyEnd(iv[f]); ASSERT_FALSE(ierr);
        ierr = VecAssemblyEnd(jv[f]); ASSERT_FALSE(ierr);
        ierr = VecAssemblyEnd(kv[f]); ASSERT_FALSE(ierr);
        
        ierr = ISLocalToGlobalMappingDestroy(&mapping); ASSERT_FALSE(ierr);
    }
    
    ierr = DMCompositeRestoreAccessArray(mesh->UPack, fp, 2, nullptr, fv); ASSERT_FALSE(ierr);
    ierr = DMCompositeRestoreAccessArray(mesh->UPack, ip, 2, nullptr, iv); ASSERT_FALSE(ierr);
    ierr = DMCompositeRestoreAccessArray(mesh->UPack, jp, 2, nullptr, jv); ASSERT_FALSE(ierr);
    ierr = DMCompositeRestoreAccessArray(mesh->UPack, kp, 2, nullptr, kv); ASSERT_FALSE(ierr);
    
    VecScatter  ctx;
    Vec fps, ips, jps, kps;
    
    ierr = VecScatterCreateToAll(fp, &ctx, &fps); ASSERT_FALSE(ierr);
    ierr = VecDuplicate(fps, &ips); ASSERT_FALSE(ierr);
    ierr = VecDuplicate(fps, &jps); ASSERT_FALSE(ierr);
    ierr = VecDuplicate(fps, &kps); ASSERT_FALSE(ierr);
    
    ierr = VecScatterBegin(ctx, fp, fps, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
    ierr = VecScatterEnd(ctx, fp, fps, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
    ierr = VecScatterBegin(ctx, ip, ips, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
    ierr = VecScatterEnd(ctx, ip, ips, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
    ierr = VecScatterBegin(ctx, jp, jps, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
    ierr = VecScatterEnd(ctx, jp, jps, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
    ierr = VecScatterBegin(ctx, kp, kps, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
    ierr = VecScatterEnd(ctx, kp, kps, INSERT_VALUES, SCATTER_FORWARD); ASSERT_FALSE(ierr);
    
    ierr = VecScatterDestroy(&ctx); ASSERT_FALSE(ierr);
    ierr = VecDestroy(&kp); ASSERT_FALSE(ierr);
    ierr = VecDestroy(&jp); ASSERT_FALSE(ierr);
    ierr = VecDestroy(&ip); ASSERT_FALSE(ierr);
    ierr = VecDestroy(&fp); ASSERT_FALSE(ierr);

    for(int f=0; f<2; ++f)
    {
        // internal points
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int j=0; j<mesh->n[f][1]; j++)
                for(int i=0; i<mesh->n[f][0]; i++)
                {
                    PetscInt idx;
                    ierr = mesh->getPackedGlobalIndex(f, i, j, k, idx); ASSERT_FALSE(ierr);
                    
                    PetscReal vf, vi, vj, vk;
                    ierr = VecGetValues(fps, 1, &idx, &vf); ASSERT_FALSE(ierr);
                    ierr = VecGetValues(ips, 1, &idx, &vi); ASSERT_FALSE(ierr);
                    ierr = VecGetValues(jps, 1, &idx, &vj); ASSERT_FALSE(ierr);
                    ierr = VecGetValues(kps, 1, &idx, &vk); ASSERT_FALSE(ierr);
                    
                    ASSERT_EQ(vf, f) << f << ".(" << i << ", " << j << ", " << k << ")" << std::endl;
                    ASSERT_EQ(vi, i) << f << ".(" << i << ", " << j << ", " << k << ")" << std::endl;
                    ASSERT_EQ(vj, j) << f << ".(" << i << ", " << j << ", " << k << ")" << std::endl;
                    ASSERT_EQ(vk, k) << f << ".(" << i << ", " << j << ", " << k << ")" << std::endl;
                }
    }
    
    ierr = VecDestroy(&kps); ASSERT_FALSE(ierr);
    ierr = VecDestroy(&jps); ASSERT_FALSE(ierr);
    ierr = VecDestroy(&ips); ASSERT_FALSE(ierr);
    ierr = VecDestroy(&fps); ASSERT_FALSE(ierr);
}

// check getPackedGlobalIndex with i, j, k of ghosts
TEST_F(CartesainMeshTest2D_YPeriodic, check_getPackedGlobalIndex_ijk_ghost)
{
    for(int f=0; f<2; ++f)
    {
        // ghost points at top and bottom
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int i=0; i<mesh->n[f][0]; i++)
            {
                PetscInt idx, idx2, j;
                
                j = mesh->n[f][1];
                mesh->getPackedGlobalIndex(f, i, j, k, idx);
                
                // internal part should have already passed
                mesh->getPackedGlobalIndex(f, i, 0, k, idx2);
                ASSERT_EQ(idx2, idx)
                    << "getPackedGlobalIndex wrong at field " << f <<" and (" 
                    << i << ", " << j << ", " << k << ") " << std::endl;
                
                j = -1;
                mesh->getPackedGlobalIndex(f, i, j, k, idx);
                
                // internal part should have already passed
                mesh->getPackedGlobalIndex(f, i, mesh->n[f][1]-1, k, idx2);
                ASSERT_EQ(idx2, idx)
                    << "getPackedGlobalIndex wrong at field " << f <<" and (" 
                    << i << ", " << j << ", " << k << ") " << std::endl;
            }
        
        // ghost points at right and left
        for(int k=0; k<mesh->n[f][2]; k++)
            for(int j=0; j<mesh->n[f][1]; j++)
            {
                PetscInt idx, i;
                
                i = mesh->n[f][0];
                mesh->getPackedGlobalIndex(f, i, j, k, idx);
                ASSERT_EQ(-1, idx)
                    << "getPackedGlobalIndex wrong at field " << f <<" and (" 
                    << i << ", " << j << ", " << k << ") " << std::endl;
                
                i = -1;
                mesh->getPackedGlobalIndex(f, i, j, k, idx);
                ASSERT_EQ(-1, idx)
                    << "getPackedGlobalIndex wrong at field " << f <<" and (" 
                    << i << ", " << j << ", " << k << ") " << std::endl;
            }
        
        // corners
        for(int k=0; k<mesh->n[f][2]; k++)
        {
            PetscInt idx;
            
            mesh->getPackedGlobalIndex(f, -1, -1, 0, idx);
            ASSERT_EQ(-1, idx)
                << "getPackedGlobalIndex wrong at field " << f <<" and (" 
                << -1 << ", " << -1 << ", " << k << ") " << std::endl;
            
            mesh->getPackedGlobalIndex(f, -1, mesh->n[f][1], 0, idx);
            ASSERT_EQ(-1, idx)
                << "getPackedGlobalIndex wrong at field " << f <<" and (" 
                << -1 << ", " << mesh->n[f][1] << ", " << k << ") " << std::endl;
            
            mesh->getPackedGlobalIndex(f, mesh->n[f][0], -1, 0, idx);
            ASSERT_EQ(-1, idx)
                << "getPackedGlobalIndex wrong at field " << f <<" and (" 
                << mesh->n[f][0] << ", " << -1 << ", " << k << ") " << std::endl;
            
            mesh->getPackedGlobalIndex(f, mesh->n[f][0], mesh->n[f][1], 0, idx);
            ASSERT_EQ(-1, idx)
                << "getPackedGlobalIndex wrong at field " << f <<" and (" 
                << mesh->n[f][0] << ", " << mesh->n[f][1] << ", " << k << ") " << std::endl;
        }
    }
}

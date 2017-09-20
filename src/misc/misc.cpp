/*
 * misc.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

#include <petibm/misc.h>


namespace petibm
{
namespace misc
{
    PetscErrorCode checkPeriodicBC(
            const type::IntVec2D &bcTypes, type::BoolVec2D &periodic)
    {
        using namespace petibm::type;
        
        PetscFunctionBeginUser;
        
        // dimension
        PetscInt dim = 3;
        
        // initialize periodic bools
        periodic = BoolVec2D(3, BoolVec1D(3, PETSC_FALSE));

        // if a boundary of a field is periodic, the counterpart should be, too
        for(int f=0; f<3; f++)
        {
            for(int d=0; d<3; d++)
            {
                bool p1 = (bcTypes[f][d*2] == int(BCType::PERIODIC));
                bool p2 = (bcTypes[f][d*2+1] == int(BCType::PERIODIC));
                
                if (p1 && (!p2))
                {
                    SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                            "Boundary %s for velocity field %s is periodic BC, "
                            "but its counterpart, %s, is not!\n",
                            bl2str[BCLoc(d*2)].c_str(), fd2str[Field(f)].c_str(),
                            bl2str[BCLoc(d*2+1)].c_str());
                }
                
                if ((!p1) && p2)
                {
                    SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                            "Boundary %s for velocity field %s is periodic BC, "
                            "but its counterpart, %s, is not!\n",
                            bl2str[BCLoc(d*2+1)].c_str(), fd2str[Field(f)].c_str(),
                            bl2str[BCLoc(d*2)].c_str());
                }
                
                if (p1 && p2) periodic[f][d] = PETSC_TRUE;
            }
        }
        
        // if all BCs of w velocity are NOBC, it's a 3D simulation
        if (std::all_of(bcTypes[2].begin(), bcTypes[2].end(),
                    [](int t){return t == int(BCType::NOBC);}))
            dim = 2;

        // if a field at a boundary is periodic BC, other fields should be, too
        for(int d=0; d<3; d++)
        {
            PetscInt c = 0;
            for(int f=0; f<3; f++)
            {
                if (periodic[f][d]) c += 1;
            }
            if ((c != 0) && (c != dim))
                SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                        "Not all velocity fields in direction %s are periodic!",
                        type::dir2str[Dir(d)].c_str());
        }

        PetscFunctionReturn(0);
    }
}
}

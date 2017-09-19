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
        PetscFunctionBeginUser;
        
        // dimension
        PetscInt dim = 0;
        
        // initialize periodic bools
        periodic = type::BoolVec2D(3, type::BoolVec1D(3, PETSC_FALSE));

        // if a boundary of a field is periodic, the counterpart should be, too
        for(int f=0; f<3; f++)
        {
            for(int d=0; d<3; d++)
            {
                bool p1 = (bcTypes[f][d*2] == int(type::BCType::PERIODIC));
                bool p2 = (bcTypes[f][d*2+1] == int(type::BCType::PERIODIC));
                
                if (p1 && (!p2))
                {
                    SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                            "Boundary %s for velocity field %s is periodic BC, "
                            "but its counterpart, %s, is not!\n",
                            type::bl2str[d*2].c_str(), type::fd2str[f].c_str(),
                            type::bl2str[d*2+1].c_str());
                }
                
                if ((!p1) && p2)
                {
                    SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                            "Boundary %s for velocity field %s is periodic BC, "
                            "but its counterpart, %s, is not!\n",
                            type::bl2str[d*2+1].c_str(), type::fd2str[f].c_str(),
                            type::bl2str[d*2].c_str());
                }
                
                if (p1 && p2) periodic[f][d] = PETSC_TRUE;
            }
        }
        
        // if all BCs of w velocity are NONE, it's a 3D simulation
        if (std::all_of(bcTypes[2].begin(), bcTypes[2].end(),
                    [](type::BCType t){return t == type::BCType::NONE;}))
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
                        type::dir2str[d].c_str());
        }

        PetscFunctionReturn(0);
    }
}
}

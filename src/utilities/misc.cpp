/*
 * misc.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */
# include "misc.h"

/** \copydoc misc::findCell1D */
PetscErrorCode misc::findCell1D(
        const PetscReal &p, const std::vector<PetscReal> &x, PetscInt &loc)
{
    PetscFunctionBeginUser;

    if (x.back() < p || x[0] > p)
        SETERRQ1(PETSC_COMM_WORLD, 56, 
                "body coordinate %e is outside domain !", p);

    auto it = std::lower_bound(x.begin(), x.end(), p);

    loc = it - x.begin() - 1;

    PetscFunctionReturn(0);
}


/** \copydoc misc::stretchGrid */
PetscErrorCode misc::stretchGrid(
        const PetscReal &bg, const PetscReal &ed, 
        const PetscInt &n, const PetscReal &r, types::RealVec1D &dL)
{
    PetscFunctionBeginUser;

    dL.resize(n);

    PetscReal                   h;

    // calculate the size of the first cell
    dL[0] = (ed - bg) * (r - 1.0) / (std::pow(r, n) - 1.0);

    // dL[i] = dL[i-1] * r
    for(auto it=dL.begin()+1; it<dL.end(); ++it) *it = *(it -1) * r;

    PetscFunctionReturn(0);
}

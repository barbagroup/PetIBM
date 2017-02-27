# include "headers.hpp"


/**
 * @brief Generate 3D Cartesian grid
 *
 * @param grid
 * @param Nx
 * @param Ny
 * @param Nz
 * @param Lx
 * @param Ly
 * @param Lz
 * @param dx
 * @param dy
 * @param dz
 * @param x
 * @param y
 * @param z
 *
 * @return 
 */
PetscErrorCode generateGrid(const DM &grid, 
        const PetscInt &Nx, const PetscInt &Ny, const PetscInt &Nz,
        const PetscScalar &Lx, const PetscScalar &Ly, const PetscScalar &Lz,
        PetscScalar &dx, PetscScalar &dy, PetscScalar &dz,
        Vec &x, Vec &y, Vec &z)
{
    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    PetscInt            xbg, xed,
                        ybg, yed,
                        zbg, zed;

    PetscScalar           ***x_arry,
                        ***y_arry,
                        ***z_arry;

    // calculate spacing
    dx = Lx / (PetscScalar) Nx;
    dy = Ly / (PetscScalar) Ny;
    dz = Lz / (PetscScalar) Nz;


    // get indices for left-bottom and right-top corner
    ierr = DMDAGetCorners(grid, &xbg, &ybg, &zbg, &xed, &yed, &zed); CHKERRQ(ierr);

    xed += xbg;
    yed += ybg;
    zed += zbg;

    // generate x coordinates
    ierr = DMDAVecGetArray(grid, x, &x_arry);                      CHKERRQ(ierr);
    for(int k=zbg; k<zed; ++k)
    {
        for(int j=ybg; j<yed; ++j)
        {
            for(int i=xbg; i<xed; ++i)
                x_arry[k][j][i] = (0.5 + (PetscScalar)i) * dx;
        }
    }
    ierr = DMDAVecRestoreArray(grid, x, &x_arry);                  CHKERRQ(ierr);

    // generate y coordinates
    ierr = DMDAVecGetArray(grid, y, &y_arry);                      CHKERRQ(ierr);
    for(int k=zbg; k<zed; ++k)
    {
        for(int j=ybg; j<yed; ++j)
        {
            for(int i=xbg; i<xed; ++i)
                y_arry[k][j][i] = (0.5 + (PetscScalar)j) * dy;
        }
    }
    ierr = DMDAVecRestoreArray(grid, y, &y_arry);                  CHKERRQ(ierr);

    // generate z coordinates
    ierr = DMDAVecGetArray(grid, z, &z_arry);                      CHKERRQ(ierr);
    for(int k=zbg; k<zed; ++k)
    {
        for(int j=ybg; j<yed; ++j)
        {
            for(int i=xbg; i<xed; ++i)
                z_arry[k][j][i] = (0.5 + (PetscScalar)k) * dz;
        }
    }
    ierr = DMDAVecRestoreArray(grid, z, &z_arry);                  CHKERRQ(ierr);

    return 0;
}


/**
 * @brief Generate 2D Cartesian grid
 *
 * @param grid
 * @param Nx
 * @param Ny
 * @param Lx
 * @param Ly
 * @param dx
 * @param dy
 * @param x
 * @param y
 *
 * @return 
 */
PetscErrorCode generateGrid(const DM &grid, 
        const PetscInt &Nx, const PetscInt &Ny,
        const PetscScalar &Lx, const PetscScalar &Ly,
        PetscScalar &dx, PetscScalar &dy, Vec &x, Vec &y)
{
    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    PetscInt            xbg, xed,
                        ybg, yed;

    PetscScalar           **x_arry,
                        **y_arry;

    // calculate spacing
    dx = Lx / (PetscScalar) Nx;
    dy = Ly / (PetscScalar) Ny;


    // get indices for left-bottom and right-top corner
    ierr = DMDAGetCorners(grid, &xbg, &ybg, nullptr, &xed, &yed, nullptr); CHKERRQ(ierr);

    xed += xbg;
    yed += ybg;

    // generate x coordinates
    ierr = DMDAVecGetArray(grid, x, &x_arry);                      CHKERRQ(ierr);
    for(int j=ybg; j<yed; ++j)
    {
        for(int i=xbg; i<xed; ++i)
            x_arry[j][i] = (0.5 + (PetscScalar)i) * dx;
    }
    ierr = DMDAVecRestoreArray(grid, x, &x_arry);                  CHKERRQ(ierr);

    // generate y coordinates
    ierr = DMDAVecGetArray(grid, y, &y_arry);                      CHKERRQ(ierr);
    for(int j=ybg; j<yed; ++j)
    {
        for(int i=xbg; i<xed; ++i)
            y_arry[j][i] = (0.5 + (PetscScalar)j) * dy;
    }
    ierr = DMDAVecRestoreArray(grid, y, &y_arry);                  CHKERRQ(ierr);

    return 0;
}



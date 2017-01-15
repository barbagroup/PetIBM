# include "headers.hpp"


/**
 * @brief Generte exact solution for the 3D Poisson equation
 *
 * @param grid
 * @param x
 * @param y
 * @param z
 * @param exact
 *
 * @return 
 */
PetscErrorCode generateExt(const DM &grid, 
        const Vec &x, const Vec &y, const Vec &z, Vec &exact)
{
    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    PetscInt            xbg, xed,
                        ybg, yed,
                        zbg, zed;

    PetscScalar           ***x_arry,
                        ***y_arry,
                        ***z_arry,
                        ***exact_arry;


    // get indices for left-bottom and right-top corner
    ierr = DMDAGetCorners(grid, &xbg, &ybg, &zbg, &xed, &yed, &zed); CHKERRQ(ierr);

    xed += xbg;
    yed += ybg;
    zed += zbg;

    // generate rhs
    ierr = DMDAVecGetArray(grid, x, &x_arry);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(grid, y, &y_arry);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(grid, z, &z_arry);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(grid, exact, &exact_arry);              CHKERRQ(ierr);
    for(int k=zbg; k<zed; ++k)
    {
        for(int j=ybg; j<yed; ++j)
        {
            for(int i=xbg; i<xed; ++i)
                exact_arry[k][j][i] = std::cos(c1 * x_arry[k][j][i]) * 
                    std::cos(c1 * y_arry[k][j][i]) * std::cos(c1 * z_arry[k][j][i]);
        }
    }
    ierr = DMDAVecRestoreArray(grid, x, &x_arry);                  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(grid, y, &y_arry);                  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(grid, z, &z_arry);                  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(grid, exact, &exact_arry);          CHKERRQ(ierr);

    return 0;
}


/**
 * @brief Generate exact solution for the 2D Poisson equation
 *
 * @param grid
 * @param x
 * @param y
 * @param exact
 *
 * @return 
 */
PetscErrorCode generateExt(const DM &grid, 
        const Vec &x, const Vec &y, Vec &exact)
{
    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    PetscInt            xbg, xed,
                        ybg, yed;

    PetscScalar           **x_arry,
                        **y_arry,
                        **exact_arry;


    // get indices for left-bottom and right-top corner
    ierr = DMDAGetCorners(grid, &xbg, &ybg, nullptr, &xed, &yed, nullptr); CHKERRQ(ierr);

    xed += xbg;
    yed += ybg;

    // generate rhs
    ierr = DMDAVecGetArray(grid, x, &x_arry);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(grid, y, &y_arry);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(grid, exact, &exact_arry);              CHKERRQ(ierr);
    for(int j=ybg; j<yed; ++j)
    {
        for(int i=xbg; i<xed; ++i)
            exact_arry[j][i] = 
                std::cos(c1 * x_arry[j][i]) * std::cos(c1 * y_arry[j][i]);
    }
    ierr = DMDAVecRestoreArray(grid, x, &x_arry);                  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(grid, y, &y_arry);                  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(grid, exact, &exact_arry);          CHKERRQ(ierr);

    return 0;
}

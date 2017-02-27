# include "headers.hpp"

/**
 * @brief Generate a matrix for 3D Poisson equation with all-Neumann BCs
 *
 * @param grid
 * @param dx
 * @param dy
 * @param dz
 * @param A
 *
 * @return 
 */
PetscErrorCode generateA(const DM &grid, 
        const PetscScalar &dx, const PetscScalar &dy, const PetscScalar &dz, Mat &A)
{
    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    PetscInt            xbg, xed,
                        ybg, yed,
                        zbg, zed;

    PetscInt            Nx, Ny, Nz;

    MatStencil          row;

    PetscScalar         Cx,
                        Cy,
                        Cz,
                        Cd;



    // get indices for left-bottom and right-top corner
    ierr = DMDAGetCorners(grid, &xbg, &ybg, &zbg, 
            &xed, &yed, &zed);                                     CHKERRQ(ierr);

    ierr = DMDAGetInfo(grid, nullptr, &Nx, &Ny, &Nz, 
            nullptr, nullptr, nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr);                   CHKERRQ(ierr);

    xed += xbg;
    yed += ybg;
    zed += zbg;

    Cx = 1.0 / dx / dx;
    Cy = 1.0 / dy / dy;
    Cz = 1.0 / dz / dz;
    Cd = -2.0 * (Cx + Cy + Cz);

    for(int k=zbg; k<zed; ++k)
    {
        for(int j=ybg; j<yed; ++j)
        {
            for(int i=xbg; i<xed; ++i)
            {
                row.i = i; row.j = j; row.k = k; row.c = 0;

                // Corners
                if (k == 0 && j == 0 && i == 0)
                {
                    MatStencil      col[4];
                    PetscScalar     values[4];

                    for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                                   col[2].j += 1;
                                   col[3].i += 1;

                    values[0] = Cd + Cz + Cy + Cx;
                                    values[1] = Cz;
                                    values[2] = Cy;
                                    values[3] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == 0 && j == 0 && i == Nx-1)
                {
                    MatStencil      col[4];
                    PetscScalar     values[4];

                    for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                                   col[2].j += 1;
                    col[3].i -= 1;               

                    values[0] = Cd + Cz + Cy + Cx;
                                    values[1] = Cz;
                                    values[2] = Cy;
                    values[3] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == 0 && j == Ny-1 && i == Nx-1)
                {
                    MatStencil      col[4];
                    PetscScalar     values[4];

                    for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                    col[2].j -= 1;               
                    col[3].i -= 1;               

                    values[0] = Cd + Cz + Cy + Cx;
                                    values[1] = Cz;
                    values[2] = Cy;                
                    values[3] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == 0 && j == Ny-1 && i == 0)
                {
                    MatStencil      col[4];
                    PetscScalar     values[4];

                    for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                    col[2].j -= 1;               
                                   col[3].i += 1;

                    values[0] = Cd + Cz + Cy + Cx;
                                    values[1] = Cz;
                    values[2] = Cy;                
                                    values[3] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1 && j == 0 && i == 0)
                {
                    MatStencil      col[4];
                    PetscScalar     values[4];

                    for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                    col[1].k -= 1;               
                                   col[2].j += 1;
                                   col[3].i += 1;

                    values[0] = Cd + Cz + Cy + Cx;
                    values[1] = Cz;                
                                    values[2] = Cy;
                                    values[3] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1 && j == 0 && i == Nx-1)
                {
                    MatStencil      col[4];
                    PetscScalar     values[4];

                    for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                    col[1].k -= 1;               
                                   col[2].j += 1;
                    col[3].i -= 1;               

                    values[0] = Cd + Cz + Cy + Cx;
                    values[1] = Cz;                
                                    values[2] = Cy;
                    values[3] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1 && j == Ny-1 && i == Nx-1)
                {
                    MatStencil      col[4];
                    PetscScalar     values[4];

                    for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                    col[1].k -= 1;               
                    col[2].j -= 1;               
                    col[3].i -= 1;               

                    values[0] = Cd + Cz + Cy + Cx;
                    values[1] = Cz;                
                    values[2] = Cy;                
                    values[3] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1 && j == Ny-1 && i == 0)
                {
                    MatStencil      col[4];
                    PetscScalar     values[4];

                    for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                    col[1].k -= 1;               
                    col[2].j -= 1;               
                                   col[3].i += 1;

                    values[0] = Cd + Cz + Cy + Cx;
                    values[1] = Cz;                
                    values[2] = Cy;                
                                    values[3] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }

                // Edges
                else if (k == 0 && j == 0)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                                   col[2].j += 1;
                    col[3].i -= 1; col[4].i += 1;

                    values[0] = Cd + Cz + Cy;
                                    values[1] = Cz;
                                    values[2] = Cy;
                    values[3] = Cx; values[4] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == 0 && j == Ny-1)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                    col[2].j -= 1;               
                    col[3].i -= 1; col[4].i += 1;

                    values[0] = Cd + Cz + Cy;
                                    values[1] = Cz;
                    values[2] = Cy;                
                    values[3] = Cx; values[4] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == 0 && i == 0)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                    col[2].j -= 1; col[3].j += 1;
                                   col[4].i += 1;

                    values[0] = Cd + Cz + Cx;
                                    values[1] = Cz;
                    values[2] = Cy; values[3] = Cy;
                                    values[4] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == 0 && i == Nx-1)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                    col[2].j -= 1; col[3].j += 1;
                    col[4].i -= 1;               

                    values[0] = Cd + Cz + Cx;
                                    values[1] = Cz;
                    values[2] = Cy; values[3] = Cy;
                    values[4] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1 && j == 0)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                    col[1].k -= 1; 
                                   col[2].j += 1;
                    col[3].i -= 1; col[4].i += 1;

                    values[0] = Cd + Cz + Cy;
                    values[1] = Cz;                
                                    values[2] = Cy;
                    values[3] = Cx; values[4] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1 && j == Ny-1)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                    col[1].k -= 1;               
                    col[2].j -= 1;               
                    col[3].i -= 1; col[4].i += 1;

                    values[0] = Cd + Cz + Cy;
                    values[1] = Cz; 
                    values[2] = Cy; 
                    values[3] = Cx; values[4] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1 && i == 0)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                    col[1].k -= 1;               
                    col[2].j -= 1; col[3].j += 1;
                                   col[4].i += 1;

                    values[0] = Cd + Cz + Cx;
                    values[1] = Cz;                
                    values[2] = Cy; values[3] = Cy;
                                  ; values[4] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1 && i == Nx-1)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                    col[1].k -= 1;               
                    col[2].j -= 1; col[3].j += 1;
                    col[4].i -= 1;               

                    values[0] = Cd + Cz + Cx;
                    values[1] = Cz;                
                    values[2] = Cy; values[3] = Cy;
                    values[4] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (j == 0 && i == 0)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                                   col[3].j += 1;
                                   col[4].i += 1;

                    values[0] = Cd + Cy + Cx;
                    values[1] = Cz; values[2] = Cz;
                                    values[3] = Cy;
                                    values[4] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (j == 0 && i == Nx-1)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                                   col[3].j += 1;
                    col[4].i -= 1;               

                    values[0] = Cd + Cy + Cx;
                    values[1] = Cz; values[2] = Cz;
                                    values[3] = Cy;
                    values[4] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (j == Ny-1 && i == 0)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                    col[3].j -= 1;               
                                   col[4].i += 1;

                    values[0] = Cd + Cy + Cx;
                    values[1] = Cz; values[2] = Cz;
                    values[3] = Cy;                
                                    values[4] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (j == Ny-1 && i == Nx-1)
                {
                    MatStencil      col[5];
                    PetscScalar     values[5];

                    for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                    col[3].j -= 1;               
                    col[4].i -= 1;               

                    values[0] = Cd + Cy + Cx;
                    values[1] = Cz; values[2] = Cz;
                    values[3] = Cy;                
                    values[4] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }

                // Faces
                else if (k == 0) // bottom face
                {
                    MatStencil      col[6];
                    PetscScalar     values[6];

                    for (int stcl=0; stcl<6; stcl++)  col[stcl] = row;

                                   col[1].k += 1;
                    col[2].j -= 1; col[3].j += 1;
                    col[4].i -= 1; col[5].i += 1;

                    values[0] = Cd + Cz;
                                  ; values[1] = Cz;
                    values[2] = Cy; values[3] = Cy;
                    values[4] = Cx; values[5] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 6, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (k == Nz-1) // top face
                {
                    MatStencil      col[6];
                    PetscScalar     values[6];

                    for (int stcl=0; stcl<6; stcl++)  col[stcl] = row;

                    col[1].k -= 1;               
                    col[2].j -= 1; col[3].j += 1;
                    col[4].i -= 1; col[5].i += 1;

                    values[0] = Cd + Cz;
                    values[1] = Cz;                
                    values[2] = Cy; values[3] = Cy;
                    values[4] = Cx; values[5] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 6, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (j == 0) // front face
                {
                    MatStencil      col[6];
                    PetscScalar     values[6];

                    for (int stcl=0; stcl<6; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                                   col[3].j += 1;
                    col[4].i -= 1; col[5].i += 1;

                    values[0] = Cd + Cy;
                    values[1] = Cz; values[2] = Cz;
                                    values[3] = Cy;
                    values[4] = Cx; values[5] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 6, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (j == Ny-1) // rear face
                {
                    MatStencil      col[6];
                    PetscScalar     values[6];

                    for (int stcl=0; stcl<6; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                    col[3].j -= 1; 
                    col[4].i -= 1; col[5].i += 1;

                    values[0] = Cd + Cy;
                    values[1] = Cz; values[2] = Cz;
                    values[3] = Cy; 
                    values[4] = Cx; values[5] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 6, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (i == 0) // left face
                {
                    MatStencil      col[6];
                    PetscScalar     values[6];

                    for (int stcl=0; stcl<6; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                    col[3].j -= 1; col[4].j += 1;
                                   col[5].i += 1;

                    values[0] = Cd + Cx;
                    values[1] = Cz; values[2] = Cz;
                    values[3] = Cy; values[4] = Cy;
                                    values[5] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 6, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
                else if (i == Nx-1) //right face
                {
                    MatStencil      col[6];
                    PetscScalar     values[6];

                    for (int stcl=0; stcl<6; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                    col[3].j -= 1; col[4].j += 1;
                    col[5].i -= 1;               

                    values[0] = Cd + Cx;
                    values[1] = Cz; values[2] = Cz;
                    values[3] = Cy; values[4] = Cy;
                    values[5] = Cx;                
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 6, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }

                // Interior
                else
                {
                    MatStencil      col[7];
                    PetscScalar     values[7];

                    for (int stcl=0; stcl<7; stcl++)  col[stcl] = row;

                    col[1].k -= 1; col[2].k += 1;
                    col[3].j -= 1; col[4].j += 1;
                    col[5].i -= 1; col[6].i += 1;

                    values[0] = Cd;
                    values[1] = Cz; values[2] = Cz;
                    values[3] = Cy; values[4] = Cy;
                    values[5] = Cx; values[6] = Cx;
               
                    ierr = MatSetValuesStencil(
                        A, 1, &row, 7, col, values, INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        }
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);                CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);                  CHKERRQ(ierr);

    return 0;
}

/**
 * @brief Generate a matrix for 2D Poisson equation with all-Neumann BCs
 *
 * @param grid
 * @param dx
 * @param dy
 * @param A
 *
 * @return 
 */
PetscErrorCode generateA(const DM &grid, 
        const PetscScalar &dx, const PetscScalar &dy, Mat &A)
{
    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    PetscInt            xbg, xed,
                        ybg, yed;

    PetscInt            Nx, Ny;

    MatStencil          row;

    PetscScalar         Cx,
                        Cy,
                        Cd;



    // get indices for left-bottom and right-top corner
    ierr = DMDAGetCorners(grid, &xbg, &ybg, nullptr, &xed, &yed, nullptr); CHKERRQ(ierr);

    ierr = DMDAGetInfo(grid, nullptr, &Nx, &Ny, nullptr, 
            nullptr, nullptr, nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr);                   CHKERRQ(ierr);

    xed += xbg;
    yed += ybg;

    Cx = 1.0 / dx / dx;
    Cy = 1.0 / dy / dy;
    Cd = -2.0 * (Cx + Cy);

    for(int j=ybg; j<yed; ++j)
    {
        for(int i=xbg; i<xed; ++i)
        {
            row.i = i; row.j = j; row.k = 0; row.c = 0;

            // Corners
            if (j == 0 && i == 0)
            {
                MatStencil      col[3];
                PetscScalar     values[3];

                for (int stcl=0; stcl<3; stcl++)  col[stcl] = row;

                               col[1].j += 1;
                               col[2].i += 1;

                values[0] = Cd + Cy + Cx;
                                values[1] = Cy;
                                values[2] = Cx;
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 3, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (j == 0 && i == Nx-1)
            {
                MatStencil      col[3];
                PetscScalar     values[3];

                for (int stcl=0; stcl<3; stcl++)  col[stcl] = row;

                               col[1].j += 1;
                col[2].i -= 1;               

                values[0] = Cd + Cy + Cx;
                                values[1] = Cy;
                values[2] = Cx;                
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 3, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (j == Ny-1 && i == Nx-1)
            {
                MatStencil      col[3];
                PetscScalar     values[3];

                for (int stcl=0; stcl<3; stcl++)  col[stcl] = row;

                col[1].j -= 1;               
                col[2].i -= 1;               

                values[0] = Cd +  Cy + Cx;
                values[1] = Cy;                
                values[2] = Cx;                
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 3, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (j == Ny-1 && i == 0)
            {
                MatStencil      col[3];
                PetscScalar     values[3];

                for (int stcl=0; stcl<3; stcl++)  col[stcl] = row;

                col[1].j -= 1;               
                               col[2].i += 1;

                values[0] = Cd +  Cy + Cx;
                values[1] = Cy;                
                                values[2] = Cx;
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 3, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }


            // Edges
            else if (j == 0)
            {
                MatStencil      col[4];
                PetscScalar     values[4];

                for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                               col[1].j += 1;
                col[2].i -= 1; col[3].i += 1;

                values[0] = Cd +  Cy;
                                values[1] = Cy;
                values[2] = Cx; values[3] = Cx;
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (j == Ny-1)
            {
                MatStencil      col[4];
                PetscScalar     values[4];

                for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                col[1].j -= 1;               
                col[2].i -= 1; col[3].i += 1;

                values[0] = Cd +  Cy;
                values[1] = Cy;                
                values[2] = Cx; values[3] = Cx;
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (i == 0)
            {
                MatStencil      col[4];
                PetscScalar     values[4];

                for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                col[1].j -= 1; col[2].j += 1;
                               col[3].i += 1;

                values[0] = Cd +  Cx;
                values[1] = Cy; values[2] = Cy;
                                values[3] = Cx;
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (i == Nx-1)
            {
                MatStencil      col[4];
                PetscScalar     values[4];

                for (int stcl=0; stcl<4; stcl++)  col[stcl] = row;

                col[1].j -= 1; col[2].j += 1;
                col[3].i -= 1;               

                values[0] = Cd +  Cx;
                values[1] = Cy; values[2] = Cy;
                values[3] = Cx;                
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 4, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }


            // Interior
            else
            {
                MatStencil      col[5];
                PetscScalar     values[5];

                for (int stcl=0; stcl<5; stcl++)  col[stcl] = row;

                col[1].j -= 1; col[2].j += 1;
                col[3].i -= 1; col[4].i += 1;

                values[0] = Cd;
                values[1] = Cy; values[2] = Cy;
                values[3] = Cx; values[4] = Cx;
           
                ierr = MatSetValuesStencil(
                    A, 1, &row, 5, col, values, INSERT_VALUES); CHKERRQ(ierr);
            }
        }
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);                CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);                  CHKERRQ(ierr);

    return 0;
}

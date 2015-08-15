/***************************************************************************//**
 * \file calculateExplicitTerms.inl
 * \author Anush Krishnan (anuhs@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the method `calculateExplicitTerms`
 *        of the class `NavierStokesSolver`.
 */


/**
 * \brief Calculates the second derivative on a non-uniform grid using a central
 *        difference scheme.
 */
inline PetscReal d2udx2(PetscReal uMinus, PetscReal uCenter, PetscReal uPlus, 
                        PetscReal dxMinus, PetscReal dxPlus)
{
  return (dxPlus*uMinus + dxMinus*uPlus - (dxPlus+dxMinus)*uCenter)*2.0/dxMinus/dxPlus/(dxMinus+dxPlus);
} // du2dx2


/**
 * \brief Calculates the explicit terms in the discretized Navier-Stokes equations. 
 * This includes the convection term, and the explicit portion of the diffusion
 * term. The velocity value at the previous time step that appears in the time 
 * discretization is also added to the explicit terms.
 *
 * The discretization scheme used to calculate the convection term is the
 * second-order conservative scheme proposed by Morinishi et al. (1998).
 *
 * A central difference scheme on a non-uniform grid is used to calculate the
 * diffusion term.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::calculateExplicitTerms()
{
  return 0;
} // calculateExplicitTerms


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::calculateExplicitTerms()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices 
           M, N,           // global number of nodes along each direction
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  PetscReal HnMinus1; // convection term at previous time-step
  PetscReal u_W, u_P, u_E, u_S, u_N, // West, Point, East, South, North x-velocities
            v_W, v_P, v_E, v_S, v_N, // West, Point, East, South, North y-velocities
            u_w, u_e, u_s, u_n,      // interpolated x-velocities at cell-vertices
            v_w, v_e, v_s, v_n;      // interpolated y-velocities at cell-vertices
  PetscReal convectionTerm, // value of convection term
            diffusionTerm;  // value of diffusion term
  
  PetscReal *dx = &mesh->dx[0],
            *dy = &mesh->dy[0];
  PetscInt nx = mesh->nx,
           ny = mesh->ny;

  PetscReal dt = parameters->dt, // time-increment
            nu = flow->nu, // viscosity
            alpha = parameters->diffusion.coefficients[1],  // explicit (n) diffusion coefficient
            gamma = parameters->convection.coefficients[1], // explicit (n) convection coefficient
            zeta  = parameters->convection.coefficients[2]; // explicit (n-1) convection coefficient

  // copy global fluxes vector to local vectors
  ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal); CHKERRQ(ierr);
  
  Vec HxGlobal, HyGlobal;
  ierr = DMCompositeGetAccess(qPack, H,  &HxGlobal, &HyGlobal); CHKERRQ(ierr);
  Vec rxGlobal, ryGlobal;
  ierr = DMCompositeGetAccess(qPack, rn, &rxGlobal, &ryGlobal); CHKERRQ(ierr);
  
  // access local vectors through multi-dimensional pointers
  PetscReal **qx, **qy;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  
  // compute explicit terms: x-component
  PetscReal **Hx, **rx;
  ierr = DMDAVecGetArray(uda, HxGlobal, &Hx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(uda, rxGlobal, &rx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      // velocity values
      u_W = qx[j][i-1]/dy[j];
      u_P = qx[j][i]/dy[j];
      u_E = qx[j][i+1]/dy[j];
      u_S = qx[j-1][i];
      if (j == 0 && flow->boundaries[YMINUS][0].type == PERIODIC) // bottom boundary and y-periodic
        u_S /= dy[ny-1];
      else if (j > 0) // inside domain
        u_S /= dy[j-1];
      u_N = qx[j+1][i];
      if (j == N-1 && flow->boundaries[YPLUS][0].type == PERIODIC) // top boundary and y-periodic
        u_N /= dy[0];
      else if (j < N-1) // inside domain
        u_N /= dy[j+1];

      // convection term
      u_w = 0.5 * (u_W + u_P);
      u_e = 0.5 * (u_P + u_E);
      u_s = 0.5 * (u_S + u_P);
      u_n = 0.5 * (u_P + u_N);
      v_s = 0.5 * (qy[j-1][i]/dx[i] + qy[j-1][i+1]/dx[i+1]);
      v_n = 0.5 * (qy[j][i]/dx[i] + qy[j][i+1]/dx[i+1]);
      // Hx = d(u^2)/dx + d(uv)/dy
      HnMinus1 = Hx[j][i];
      Hx[j][i] = (  (u_e*u_e - u_w*u_w)/(0.5*(dx[i]+dx[i+1])) 
                  + (u_n*v_n - u_s*v_s)/dy[j] );
      convectionTerm = gamma*Hx[j][i] + zeta*HnMinus1;

      // diffusion term
      diffusionTerm = alpha*nu * (  d2udx2(u_W, u_P, u_E, dx[i], dx[i+1])
                                  + d2udx2(u_S, u_P, u_N, 0.5*(dy[j-1]+dy[j]), 0.5*(dy[j]+dy[j+1])) );

      // explicit term
      rx[j][i] = u_P/dt - convectionTerm + diffusionTerm;
    }
  }
  ierr = DMDAVecRestoreArray(uda, HxGlobal, &Hx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, rxGlobal, &rx); CHKERRQ(ierr);
  
  // compute explicit terms: y-component
  PetscReal **Hy, **ry;
  ierr = DMDAVecGetArray(vda, HyGlobal, &Hy); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(vda, ryGlobal, &ry); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      // velocity values
      v_W = qy[j][i-1];
      if (i == 0 && flow->boundaries[XMINUS][1].type == PERIODIC) // left boundary and x-periodic
        v_W /= dx[nx-1];
      else if (i > 0) // inside domain
        v_W /= dx[i-1];
      v_P = qy[j][i]/dx[i];
      v_E = qy[j][i+1];
      if (i == M-1 && flow->boundaries[XPLUS][1].type == PERIODIC) // right boundary and x-periodic
        v_E /= dx[0];
      else if (i < M-1) // inside domain
        v_E /= dx[i+1];
      v_S = qy[j-1][i]/dx[i];
      v_N = qy[j+1][i]/dx[i];

      // convection term
      u_w = 0.5 * (qx[j][i]/dy[j] + qx[j+1][i]/dy[j+1]);
      u_e = 0.5 * (qx[j][i+1]/dy[j] + qx[j+1][i+1]/dy[j+1]);
      v_w = 0.5 * (v_W + v_P); 
      v_e = 0.5 * (v_P + v_E);
      v_s = 0.5 * (v_S + v_P);
      v_n = 0.5 * (v_P + v_N);
      // Hy = d(uv)/dx + d(v^2)/dy
      HnMinus1 = Hy[j][i];
      Hy[j][i] = (  (u_e*v_e - u_w*v_w)/dx[i] 
                  + (v_n*v_n - v_s*v_s)/(0.5*(dy[j]+dy[j+1])) );
      convectionTerm = gamma*Hy[j][i] + zeta*HnMinus1;

      // diffusion term
      diffusionTerm = alpha*nu * (  d2udx2(v_W, v_P, v_E, 0.5*(dx[i-1]+dx[i]), 0.5*(dx[i]+dx[i+1]))
                                  + d2udx2(v_S, v_P, v_N, dy[j], dy[j+1]) );

      // explicit term
      ry[j][i] = v_P/dt - convectionTerm + diffusionTerm;
    }
  }
  ierr = DMDAVecRestoreArray(vda, HyGlobal, &Hy); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, ryGlobal, &ry); CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  
  ierr = DMCompositeRestoreAccess(qPack, H,  &HxGlobal, &HyGlobal); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(qPack, rn, &rxGlobal, &ryGlobal); CHKERRQ(ierr);

  return 0;
} // calculateExplicitTerms


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::calculateExplicitTerms()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           M, N, P,                // global number of nodes along each direction
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices

  PetscReal HnMinus1; // convection term at previous time-step
  PetscReal u_W, u_P, u_E, u_S, u_N, u_B, u_F, // West, Point, East, South, North, Back, Front x-velocities
            v_W, v_P, v_E, v_S, v_N, v_B, v_F, // West, Point, East, South, North, Back, Front y-velocities
            w_W, w_P, w_E, w_S, w_N, w_B, w_F, // West, Point, East, South, North, Back, Front z-velocities
            u_w, u_e, u_s, u_n, u_b, u_f, // interpolated x-velocities
            v_w, v_e, v_s, v_n, v_b, v_f, // interpolated y-velocities
            w_w, w_e, w_s, w_n, w_b, w_f; // interpolated z-velocities
  PetscReal convectionTerm, // value of convection term
            diffusionTerm;  // value of diffusion term

  PetscReal *dx = &mesh->dx[0],
            *dy = &mesh->dy[0],
            *dz = &mesh->dz[0];
  PetscInt nx = mesh->nx,
           ny = mesh->ny,
           nz = mesh->nz;

  PetscReal nu = flow->nu, // viscosity
            dt = parameters->dt, // time-increment
            alpha = parameters->diffusion.coefficients[1],  // explicit (n) diffusion coefficient
            gamma = parameters->convection.coefficients[1], // explicit (n) convection coefficient
            zeta  = parameters->convection.coefficients[2]; // explicit (n-1) convection coefficient

  // copy global fluxes vector to local vectors
  ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal, qzLocal); CHKERRQ(ierr);

  Vec HxGlobal, HyGlobal, HzGlobal;  
  ierr = DMCompositeGetAccess(qPack, H,  &HxGlobal, &HyGlobal, &HzGlobal); CHKERRQ(ierr);
  Vec rxGlobal, ryGlobal, rzGlobal;
  ierr = DMCompositeGetAccess(qPack, rn, &rxGlobal, &ryGlobal, &rzGlobal); CHKERRQ(ierr);
  
  // access local vectors through multi-dimensional pointers
  PetscReal ***qx, ***qy, ***qz;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRQ(ierr);
  
  // compute explicit terms: x-component
  PetscReal ***Hx, ***rx;
  ierr = DMDAVecGetArray(uda, HxGlobal, &Hx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(uda, rxGlobal, &rx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // velocity values
        u_W = qx[k][j][i-1]/(dy[j]*dz[k]);
        u_P = qx[k][j][i]/(dy[j]*dz[k]);
        u_E = qx[k][j][i+1]/(dy[j]*dz[k]);
        u_S = qx[k][j-1][i];
        if (j == 0 && flow->boundaries[YMINUS][0].type == PERIODIC) // bottom boundary and y-periodic
          u_S /= dy[ny-1]*dz[k];
        else if (j > 0) // inside domain
          u_S /= dy[j-1]*dz[k];
        u_N = qx[k][j+1][i];
        if (j == N-1 && flow->boundaries[YPLUS][0].type == PERIODIC) // top boundary and y-periodic
          u_N /= dy[0]*dz[k];
        else if (j < N-1) // inside domain
          u_N /= dy[j+1]*dz[k];
        u_B = qx[k-1][j][i];
        if (k == 0 && flow->boundaries[ZMINUS][0].type == PERIODIC) // back boundary and z-periodic
          u_B /= dy[j]*dz[nz-1];
        else if (k < P-1) // inside domain
          u_B /= dy[j]*dz[k];
        u_F = qx[k+1][j][i];
        if (k == 0 && flow->boundaries[ZPLUS][0].type == PERIODIC) // front boundary and z-periodic
          u_F /= dy[j]*dz[0];
        else if (k < P-1) // inside domain
          u_B /= dy[j]*dz[k];

        // convection term
        u_w = 0.5 * (u_W + u_P);
        u_e = 0.5 * (u_P + u_E);
        u_s = 0.5 * (u_S + u_P);
        u_n = 0.5 * (u_P + u_N);
        u_b = 0.5 * (u_B + u_P);
        u_f = 0.5 * (u_P + u_F);
        v_s = 0.5 * (qy[k][j-1][i]/(dx[i]*dz[k]) + qy[k][j-1][i+1]/(dx[i+1]*dz[k]));
        v_n = 0.5 * (qy[k][j][i]/(dx[i]*dz[k]) + qy[k][j][i+1]/(dx[i+1]*dz[k]));
        w_b = 0.5 * (qz[k-1][j][i]/(dx[i]*dy[j]) + qz[k-1][j][i+1]/(dx[i]*dy[j]));
        w_f = 0.5 * (qz[k][j][i]/(dx[i]*dy[j]) + qz[k][j][i+1]/(dx[i]*dy[j]));
        // Hx = d(u^2)/dx + d(uv)/dy + d(uw)/dz
        HnMinus1 = Hx[j][j][i];
        Hx[k][j][i] = (  (u_e*u_e - u_w*u_w)/(0.5*(dx[i]+dx[i+1])) 
                       + (u_n*v_n - u_s*v_s)/dy[j]
                       + (u_f*w_f - u_b*w_b)/dz[k] );
        convectionTerm = gamma*Hx[k][j][i] + zeta*HnMinus1;

        // diffusion term
        diffusionTerm = alpha*nu * (  d2udx2(u_W, u_P, u_E, dx[i], dx[i+1])
                                    + d2udx2(u_S, u_P, u_N, 0.5*(dy[j-1]+dy[j]), 0.5*(dy[j]+dy[j+1]))
                                    + d2udx2(u_B, u_P, u_F, 0.5*(dz[k-1]+dz[k]), 0.5*(dz[k]+dz[k+1])) );

        // explicit term
        rx[k][j][i] = u_P/dt - convectionTerm + diffusionTerm;
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, HxGlobal, &Hx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, rxGlobal, &rx); CHKERRQ(ierr);

  // y-component
  PetscReal ***Hy, ***ry;
  ierr = DMDAVecGetArray(vda, HyGlobal, &Hy); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(vda, ryGlobal, &ry); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // velocity values
        v_W = qy[k][j][i-1];
        if (i == 0 && flow->boundaries[XMINUS][1].type == PERIODIC) // left boundary and x-periodic
          v_W /= dx[nx-1]*dz[k];
        else if (i > 0) // inside domain
          v_W /= dx[i-1]*dz[k];
        v_P = qy[k][j][i]/(dx[i]*dz[k]);
        v_E = qy[k][j][i+1];
        if (i == M-1 && flow->boundaries[XPLUS][1].type == PERIODIC) // right boundary and x-periodic
          v_E /= dx[0]*dz[k];
        else if (i < M-1) // inside domain
          v_E /= dx[i+1]*dz[k];
        v_S = qy[k][j-1][i]/(dx[i]*dz[k]);
        v_N = qy[k][j+1][i]/(dx[i]*dz[k]);
        v_B = qy[k-1][j][i];
        if (k == 0 && flow->boundaries[ZMINUS][1].type == PERIODIC) // back boundary and z-periodic
          v_B /= dx[i]*dz[nz-1];
        else if (k > 0) // inside domain
          v_B /= dx[i]*dz[k-1];
        v_F = qy[k+1][j][i];
        if (k == P-1 && flow->boundaries[ZPLUS][1].type == PERIODIC) // front boundary and z-periodic
          v_F /= dx[i]*dz[0];
        else if (k < P-1) // inside domain
          v_F /= dx[i]*dz[k+1];

        // convection term
        u_w = 0.5 * (qx[k][j][i]/(dy[j]*dz[k]) + qx[k][j+1][i]/(dy[j+1]*dz[k]));
        u_e = 0.5 * (qx[k][j][i+1]/(dy[j]*dz[k]) + qx[k][j+1][i+1]/(dy[j+1]*dz[k]));
        v_w = 0.5 * (v_W + v_P); 
        v_e = 0.5 * (v_P + v_E);
        v_s = 0.5 * (v_S + v_P);
        v_n = 0.5 * (v_P + v_N);
        v_b = 0.5 * (v_B + v_P);
        v_f = 0.5 * (v_P + v_F);
        w_b = 0.5 * (qz[k-1][j][i]/(dx[i]*dy[j]) + qz[k-1][j][i+1]/(dx[i]*dy[j]));
        w_f = 0.5 * (qz[k][j][i]/(dx[i]*dy[j]) + qz[k][j][i+1]/(dx[i]*dy[j]));
        // Hy = d(vu)/dx + d(v^2)/dy + d(vw)/dz
        HnMinus1 = Hy[k][j][i];
        Hy[k][j][i] = (  (u_e*v_e - u_w*v_w)/(0.5*(dx[i]+dx[i+1])) 
                       + (v_n*v_n - v_s*v_s)/dy[j]
                       + (v_f*w_f - v_b*w_b)/(0.5*(dz[k]+dz[k+1])) );
        convectionTerm = gamma*Hy[k][j][i] + zeta*HnMinus1;

        // diffusion term
        diffusionTerm = alpha*nu * (  d2udx2(v_W, v_P, v_E, 0.5*(dx[i-1]+dx[i]), 0.5*(dx[i]+dx[i+1]))
                                    + d2udx2(v_S, v_P, v_N, dy[j], dy[j+1])
                                    + d2udx2(v_B, v_P, v_F, 0.5*(dz[k-1]+dz[k]), 0.5*(dz[k]+dz[k+1])) );

        // explicit term
        ry[k][j][i] = v_P/dt - convectionTerm + diffusionTerm;
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, HyGlobal, &Hy); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, ryGlobal, &ry); CHKERRQ(ierr);

  // z-component
  PetscReal ***Hz, ***rz;
  ierr = DMDAVecGetArray(wda, HzGlobal, &Hz); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(wda, rzGlobal, &rz); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(wda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // velocity values
        w_P = qz[k][j][i]/(dx[i]*dy[j]);
        w_W = qz[k][j][i-1];
        if (i == 0 && flow->boundaries[XMINUS][2].type == PERIODIC) // left boundary and x-periodic
          w_W /= dx[nx-1]*dy[j];
        else if (i > 0) // inside domain
          w_W /= dx[i-1]*dy[j];
        w_E = qz[k][j][i+1];
        if (i == M-1 && flow->boundaries[XPLUS][2].type == PERIODIC) // right boundary and x-periodic
          w_E /= dx[0]*dy[j];
        else if (i < M-1) // inside domain
          w_E /= dx[i+1]*dy[j];
        w_S = qz[k][j-1][i];
        if (j == 0 && flow->boundaries[YMINUS][2].type == PERIODIC) // bottom boundary and y-periodic
          w_S /= dx[i]*dy[ny-1];
        else if (j > 0) // inside domain
          w_S /= dx[i]*dy[j-1];
        w_N = qz[k][j+1][i];
        if (j == N-1 && flow->boundaries[YPLUS][2].type == PERIODIC) // top boundary and y-periodic
          w_N /= dx[i]*dy[0];
        else if (j < N-1) // inside domain
          w_N /= dx[i]*dy[j+1];
        w_B = qz[k-1][j][i]/(dx[i]*dy[j]);
        w_F = qz[k+1][j][i]/(dx[i]*dy[j]);

        // convection term
        u_w = 0.5 * (qx[k][j][i]/(dy[j]*dz[k]) + qx[k+1][j][i]/(dy[j]*dz[k+1]));
        u_e = 0.5 * (qx[k][j][i+1]/(dy[j]*dz[k]) + qx[k+1][j][i+1]/(dy[j]*dz[k+1]));
        v_s = 0.5 * (qy[k][j][i]/(dx[i]*dz[k]) + qy[k+1][j][i]/(dx[i]*dz[k+1]));
        v_n = 0.5 * (qy[k][j+1][i]/(dx[i]*dz[k]) + qy[k+1][j+1][i]/(dx[i]*dz[k+1]));
        w_w = 0.5 * (w_W + w_P);
        w_e = 0.5 * (w_P + w_E);
        w_s = 0.5 * (w_S + w_P);
        w_n = 0.5 * (w_P + w_N);
        w_b = 0.5 * (w_B + w_P);
        w_f = 0.5 * (w_P + w_F);
        // Hz = d(wu)/dx + d(wv)/dy + d(w^2)/dz
        HnMinus1 = Hz[k][j][i];
        Hz[k][j][i] = (  (u_e*w_e - u_w*w_w)/(0.5*(dx[i]+dx[i+1])) 
                       + (v_n*w_n - v_s*w_s)/(0.5*(dy[j]+dy[j+1]))
                       + (w_f*w_f - w_b*w_b)/dz[k] );
        convectionTerm = gamma*Hz[k][j][i] + zeta*HnMinus1;

        // diffusion term
        diffusionTerm = alpha*nu * (  d2udx2(w_W, w_P, w_E, 0.5*(dx[i-1]+dx[i]), 0.5*(dx[i]+dx[i+1]))
                                    + d2udx2(w_S, w_P, w_N, 0.5*(dy[j-1]+dy[j]), 0.5*(dy[j]+dy[j+1]))
                                    + d2udx2(v_B, v_P, v_F, dz[k], dz[k+1]) );

        // explicit term
        rz[k][j][i] = w_P/dt - convectionTerm + diffusionTerm;
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, HzGlobal, &Hz); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(wda, rzGlobal, &rz); CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRQ(ierr);
  
  ierr = DMCompositeRestoreAccess(qPack, H,  &HxGlobal, &HyGlobal, &HzGlobal); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(qPack, rn, &rxGlobal, &ryGlobal, &rzGlobal); CHKERRQ(ierr);

  return 0;
} // calculateExplicitTerms
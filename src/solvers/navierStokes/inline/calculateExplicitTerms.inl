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
} // d2udx2


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
  PetscReal dxMinus, dxPlus, dyMinus, dyPlus;
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

  PetscBool periodicX = (flow->boundaries[XMINUS][0].type == PERIODIC) ? PETSC_TRUE : PETSC_FALSE,
            periodicY = (flow->boundaries[YMINUS][0].type == PERIODIC) ? PETSC_TRUE : PETSC_FALSE;

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
      // velocity value at nodes
      u_P = qx[j][i]/dy[j];
      u_W = qx[j][i-1]/dy[j];
      u_E = qx[j][i+1]/dy[j];
      u_S = (j > 0) ? qx[j-1][i]/dy[j-1] : (periodicY) ? qx[j-1][i]/dy[ny-1] : qx[j-1][i];
      u_N = (j < N-1) ? qx[j+1][i]/dy[j+1] : (periodicY) ? qx[j+1][i]/dy[0] : qx[j+1][i];
      // interpolated velocity values
      u_w = 0.5*(u_W + u_P);
      u_e = 0.5*(u_P + u_E);
      u_s = (j == 0 && !periodicY) ? u_S : 0.5*(u_S + u_P);
      u_n = (j == N-1 && !periodicY) ? u_N : 0.5*(u_P + u_N);      
      dxMinus = dx[i];
      dxPlus = (i == M-1 && periodicX) ? dx[0] : dx[i+1];
      v_s = 0.5*(qy[j-1][i]/dxMinus + qy[j-1][i+1]/dxPlus);
      v_n = 0.5*(qy[j][i]/dxMinus + qy[j][i+1]/dxPlus);
      // convection term: Hx = d(u^2)/dx + d(uv)/dy
      HnMinus1 = Hx[j][i];
      Hx[j][i] = (  (u_e*u_e - u_w*u_w)/(0.5*(dxMinus+dxPlus)) 
                  + (v_n*u_n - v_s*u_s)/dy[j] );
      convectionTerm = gamma*Hx[j][i] + zeta*HnMinus1;
      // diffusion term: d^2u/dx^2 + d^2u/dy^2
      dyMinus = (j > 0) ? 0.5*(dy[j-1] + dy[j]) : (periodicY) ? 0.5*(dy[ny-1] + dy[j]) : 0.5*dy[j];
      dyPlus = (j < N-1) ? 0.5*(dy[j] + dy[j+1]) : (periodicY) ? 0.5*(dy[j] + dy[0]) : 0.5*dy[j];
      diffusionTerm = alpha*nu * (  d2udx2(u_W, u_P, u_E, dxMinus, dxPlus)
                                  + d2udx2(u_S, u_P, u_N, dyMinus, dyPlus) );
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
      // velocity value at nodes
      v_P = qy[j][i]/dx[i];
      v_W = (i > 0) ? qy[j][i-1]/dx[i-1] : (periodicX) ? qy[j][i-1]/dx[nx-1] : qy[j][i-1];
      v_E = (i < M-1) ? qy[j][i+1]/dx[i+1] : (periodicX) ? qy[j][i+1]/dx[0] : qy[j][i+1];
      v_S = qy[j-1][i]/dx[i];
      v_N = qy[j+1][i]/dx[i];
      // interpolated velocity values
      v_w = (i == 0 && !periodicX) ? v_W : 0.5*(v_W + v_P);
      v_e = (i == M-1 && !periodicX) ? v_E : 0.5*(v_P + v_E);
      v_s = 0.5*(v_S + v_P);
      v_n = 0.5*(v_P + v_N);
      dyMinus = dy[j];
      dyPlus = (j == N-1 && periodicY) ? dy[0] : dy[j+1];
      u_w = 0.5*(qx[j][i-1]/dyMinus + qx[j+1][i-1]/dyPlus);
      u_e = 0.5*(qx[j][i]/dyMinus + qx[j+1][i]/dyPlus);
      // convection term: Hy = d(uv)/dx + d(v^2)/dy
      HnMinus1 = Hy[j][i];
      Hy[j][i] = (  (u_e*v_e - u_w*v_w)/dx[i] 
                  + (v_n*v_n - v_s*v_s)/(0.5*(dyMinus+dyPlus)) );
      convectionTerm = gamma*Hy[j][i] + zeta*HnMinus1;
      // diffusion term: d^2v/dx^2 + d^2v/dy^2
      dxMinus = (i > 0) ? 0.5*(dx[i-1] + dx[i]) : (periodicX) ? 0.5*(dx[nx-1] + dx[i]) : 0.5*dx[i];
      dxPlus = (i < M-1) ? 0.5*(dx[i] + dx[i+1]) : (periodicX) ? 0.5*(dx[i] + dx[0]) : 0.5*dx[i];
      diffusionTerm = alpha*nu * (  d2udx2(v_W, v_P, v_E, dxMinus, dxPlus)
                                  + d2udx2(v_S, v_P, v_N, dyMinus, dyPlus) );
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

  PetscReal dxMinus, dxPlus, dyMinus, dyPlus, dzMinus, dzPlus;

  PetscBool periodicX = (flow->boundaries[XMINUS][0].type == PERIODIC) ? PETSC_TRUE : PETSC_FALSE,
            periodicY = (flow->boundaries[YMINUS][0].type == PERIODIC) ? PETSC_TRUE : PETSC_FALSE,
            periodicZ = (flow->boundaries[ZMINUS][0].type == PERIODIC) ? PETSC_TRUE : PETSC_FALSE;

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
        // velocity value at nodes
        u_P = qx[k][j][i]/(dy[j]*dz[k]);
        u_W = qx[k][j][i-1]/(dy[j]*dz[k]);
        u_E = qx[k][j][i+1]/(dy[j]*dz[k]);
        u_S = (j > 0)   ? qx[k][j-1][i]/(dy[j-1]*dz[k]) : (periodicY) ? qx[k][j-1][i]/(dy[ny-1]*dz[k]) : qx[k][j-1][i];
        u_N = (j < N-1) ? qx[k][j+1][i]/(dy[j+1]*dz[k]) : (periodicY) ? qx[k][j+1][i]/(dy[0]*dz[k])    : qx[k][j+1][i];
        u_B = (k > 0)   ? qx[k-1][j][i]/(dy[j]*dz[k-1]) : (periodicZ) ? qx[k-1][j][i]/(dy[j]*dz[nz-1]) : qx[k-1][j][i];
        u_F = (k < P-1) ? qx[k+1][j][i]/(dy[j]*dz[k+1]) : (periodicZ) ? qx[k+1][j][i]/(dy[j]*dz[0])    : qx[k+1][j][i];
        // interpolated velocity values
        u_w = 0.5*(u_W + u_P);
        u_e = 0.5*(u_P + u_E);
        u_s = (j == 0 && !periodicY)   ? u_S : 0.5*(u_S + u_P);
        u_n = (j == N-1 && !periodicY) ? u_N : 0.5*(u_P + u_N);
        u_b = (k == 0 && !periodicZ)   ? u_B : 0.5*(u_B + u_P);
        u_f = (k == P-1 && !periodicZ) ? u_F : 0.5*(u_P + u_F);
        dxMinus = dx[i];
        dxPlus = (i == M-1 && periodicX) ? dx[0] : dx[i+1];
        v_s = 0.5*(qy[k][j-1][i]/(dxMinus*dz[k]) + qy[k][j-1][i+1]/(dxPlus*dz[k]));
        v_n = 0.5*(qy[k][j][i]/(dxMinus*dz[k])   + qy[k][j][i+1]/(dxPlus*dz[k]));
        w_b = 0.5*(qz[k-1][j][i]/(dxMinus*dy[j]) + qz[k-1][j][i+1]/(dxPlus*dy[j]));
        w_f = 0.5*(qz[k][j][i]/(dxMinus*dy[j])   + qz[k][j][i+1]/(dxPlus*dy[j]));
        // convection term: Hx = d(u^2)/dx + d(uv)/dy + d(uw)/dz
        HnMinus1 = Hx[k][j][i];
        Hx[k][j][i] = (  (u_e*u_e - u_w*u_w)/(0.5*(dxMinus + dxPlus)) 
                       + (v_n*u_n - v_s*u_s)/dy[j]
                       + (w_f*u_f - w_b*u_b)/dz[k] );
        convectionTerm = gamma*Hx[k][j][i] + zeta*HnMinus1;
        // diffusion term: d^2u/dx^2 + d^2u/dy^2 + d^2u/dz^2
        dyMinus = (j > 0)  ? 0.5*(dy[j-1] + dy[j]) : (periodicY) ? 0.5*(dy[ny-1] + dy[j]) : 0.5*dy[j];
        dyPlus = (j < N-1) ? 0.5*(dy[j] + dy[j+1]) : (periodicY) ? 0.5*(dy[j] + dy[0])    : 0.5*dy[j];
        dzMinus = (k > 0)  ? 0.5*(dz[k-1] + dz[k]) : (periodicZ) ? 0.5*(dz[nz-1] + dz[k]) : 0.5*dz[k];
        dzPlus = (k < P-1) ? 0.5*(dz[k] + dz[k+1]) : (periodicZ) ? 0.5*(dz[k] + dz[0])    : 0.5*dz[k];
        diffusionTerm = alpha*nu * (  d2udx2(u_W, u_P, u_E, dxMinus, dxPlus)
                                    + d2udx2(u_S, u_P, u_N, dyMinus, dyPlus)
                                    + d2udx2(u_B, u_P, u_F, dzMinus, dzPlus) );
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
        // velocity value at nodes
        v_P = qy[k][j][i]/(dx[i]*dz[k]);
        v_W = (i > 0)   ? qy[k][j][i-1]/(dx[i-1]*dz[k]) : (periodicX) ? qy[k][j][i-1]/(dx[nx-1]*dz[k]) : qy[k][j][i-1];
        v_E = (i < M-1) ? qy[k][j][i+1]/(dx[i+1]*dz[k]) : (periodicX) ? qy[k][j][i+1]/(dx[0]*dz[k])    : qy[k][j][i+1];
        v_S = qy[k][j-1][i]/(dx[i]*dz[k]);
        v_N = qy[k][j+1][i]/(dx[i]*dz[k]);
        v_B = (k > 0)   ? qy[k-1][j][i]/(dx[i]*dz[k-1]) : (periodicZ) ? qy[k-1][j][i]/(dx[i]*dz[nz-1]) : qy[k-1][j][i];
        v_F = (k < P-1) ? qy[k+1][j][i]/(dx[i]*dz[k+1]) : (periodicZ) ? qy[k+1][j][i]/(dx[i]*dz[0])    : qy[k+1][j][i];
        // interpolated velocity values
        v_w = (i == 0 && !periodicX)   ? v_W : 0.5*(v_W + v_P);
        v_e = (i == M-1 && !periodicX) ? v_E : 0.5*(v_P + v_E);
        v_s = 0.5*(v_S + v_P);
        v_n = 0.5*(v_P + v_N);
        v_b = (k == 0 && !periodicZ)   ? v_B : 0.5*(v_B + v_P);
        v_f = (k == P-1 && !periodicZ) ? v_F : 0.5*(v_P + v_F);
        dyMinus = dy[j];
        dyPlus = (j == N-1 && periodicY) ? dy[0] : dy[j+1];
        u_w = 0.5*(qx[k][j][i-1]/(dyMinus*dz[k]) + qx[k][j+1][i-1]/(dyPlus*dz[k]));
        u_e = 0.5*(qx[k][j][i]/(dyMinus*dz[k]) + qx[k][j+1][i]/(dyPlus*dz[k]));
        w_b = 0.5*(qz[k-1][j][i]/(dx[i]*dyMinus) + qz[k-1][j+1][i]/(dx[i]*dyPlus));
        w_f = 0.5*(qz[k][j][i]/(dx[i]*dyMinus) + qz[k][j+1][i]/(dx[i]*dyPlus));
        // convection term: Hy = d(vu)/dx + d(v^2)/dy + d(vw)/dz
        HnMinus1 = Hy[k][j][i];
        Hy[k][j][i] = (  (u_e*v_e - u_w*v_w)/dx[i] 
                       + (v_n*v_n - v_s*v_s)/(0.5*(dyMinus + dyPlus))
                       + (w_f*v_f - w_b*v_b)/dz[k] );
        convectionTerm = gamma*Hy[k][j][i] + zeta*HnMinus1;
        // diffusion term: d^2v/dx^2 + d^2v/dy^2 + d^2v/dz^2
        dxMinus = (i > 0)  ? 0.5*(dx[i-1] + dx[i]) : (periodicX) ? 0.5*(dx[nx-1] + dx[i]) : 0.5*dx[i];
        dxPlus = (i < M-1) ? 0.5*(dx[i] + dx[i+1]) : (periodicX) ? 0.5*(dx[i] + dx[0])    : 0.5*dx[i];
        dzMinus = (k > 0)  ? 0.5*(dz[k-1] + dz[k]) : (periodicZ) ? 0.5*(dz[nz-1] + dz[k]) : 0.5*dz[k];
        dzPlus = (k < P-1) ? 0.5*(dz[k] + dz[k+1]) : (periodicZ) ? 0.5*(dz[k] + dz[0])    : 0.5*dz[k];
        diffusionTerm = alpha*nu * (  d2udx2(v_W, v_P, v_E, dxMinus, dxPlus)
                                    + d2udx2(v_S, v_P, v_N, dyMinus, dyPlus)
                                    + d2udx2(v_B, v_P, v_F, dzMinus, dzPlus) );
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
        // velocity value at nodes
        w_P = qz[k][j][i]/(dx[i]*dy[j]);
        w_W = (i > 0)   ? qz[k][j][i-1]/(dx[i-1]*dy[j]) : (periodicX) ? qz[k][j][i-1]/(dx[nx-1]*dy[j]) : qz[k][j][i-1];
        w_E = (i < M-1) ? qz[k][j][i+1]/(dx[i+1]*dy[j]) : (periodicX) ? qz[k][j][i+1]/(dx[0]*dy[j])    : qz[k][j][i+1];
        w_S = (j > 0)   ? qz[k][j-1][i]/(dx[i]*dy[j-1]) : (periodicY) ? qz[k][j-1][i]/(dx[i]*dy[ny-1]) : qz[k][j-1][i];
        w_N = (j < N-1) ? qz[k][j+1][i]/(dx[i]*dy[j+1]) : (periodicY) ? qz[k][j+1][i]/(dx[i]*dy[0])    : qz[k][j+1][i];
        w_B = qz[k-1][j][i]/(dx[i]*dy[j]);
        w_F = qz[k+1][j][i]/(dx[i]*dy[j]);
        // interpolated velocity values
        w_w = (i == 0 && !periodicX)   ? w_W : 0.5*(w_W + w_P);
        w_e = (i == M-1 && !periodicX) ? w_E : 0.5*(w_P + w_E);
        w_s = (j == 0 && !periodicY)   ? w_S : 0.5*(w_S + w_P);
        w_n = (j == N-1 && !periodicY) ? w_N : 0.5*(w_P + w_N);
        w_b = 0.5*(w_B + w_P);
        w_f = 0.5*(w_P + w_F);
        dzMinus = dz[k];
        dzPlus = (k == P-1 && periodicZ) ? dz[0] : dz[k+1];
        u_w = 0.5*(qx[k][j][i-1]/(dy[j]*dzMinus) + qx[k+1][j][i-1]/(dy[j]*dzPlus));
        u_e = 0.5*(qx[k][j][i]/(dy[j]*dzMinus) + qx[k+1][j][i]/(dy[j]*dzPlus));
        v_s = 0.5*(qy[k][j-1][i]/(dx[i]*dzMinus) + qy[k+1][j-1][i]/(dx[i]*dzPlus));
        v_n = 0.5*(qy[k][j][i]/(dx[i]*dzMinus) + qy[k+1][j][i]/(dx[i]*dzPlus));
        // convection term: Hz = d(wu)/dx + d(wv)/dy + d(w^2)/dz
        HnMinus1 = Hz[k][j][i];
        Hz[k][j][i] = (  (u_e*w_e - u_w*w_w)/dx[i] 
                       + (v_n*w_n - v_s*w_s)/dy[j]
                       + (w_f*w_f - w_b*w_b)/(0.5*(dzMinus + dzPlus)) );
        convectionTerm = gamma*Hz[k][j][i] + zeta*HnMinus1;
        // diffusion term: d^2w/dx^2 + d^2w/dy^2 + d^2w/dz^2
        dxMinus = (i > 0)  ? 0.5*(dx[i-1] + dx[i]) : (periodicX) ? 0.5*(dx[nx-1] + dx[i]) : 0.5*dx[i];
        dxPlus = (i < M-1) ? 0.5*(dx[i] + dx[i+1]) : (periodicX) ? 0.5*(dx[i] + dx[0])    : 0.5*dx[i];
        dyMinus = (j > 0)  ? 0.5*(dy[j-1] + dy[j]) : (periodicY) ? 0.5*(dy[ny-1] + dy[j]) : 0.5*dy[j];
        dyPlus = (j < N-1) ? 0.5*(dy[j] + dy[j+1]) : (periodicY) ? 0.5*(dy[j] + dy[0])    : 0.5*dy[j];
        diffusionTerm = alpha*nu * (  d2udx2(w_W, w_P, w_E, dxMinus, dxPlus)
                                    + d2udx2(w_S, w_P, w_N, dyMinus, dyPlus)
                                    + d2udx2(v_B, v_P, v_F, dzMinus, dzPlus) );
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
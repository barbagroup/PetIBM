#include <petscdmcomposite.h>

inline PetscReal du2dx2(PetscReal uMinus, PetscReal uCenter, PetscReal uPlus, PetscReal dxMinus, PetscReal dxPlus)
{
	return (dxPlus*uMinus + dxMinus*uPlus - (dxPlus+dxMinus)*uCenter)*2.0/dxMinus/dxPlus/(dxMinus+dxPlus);
}

template<>
void NavierStokesSolver<2>::calculateExplicitTerms()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j, M, N;
	Vec            HxGlobal, HyGlobal;
	Vec            rxGlobal, ryGlobal;
	PetscReal      **qx, **qy;
	PetscReal      **Hx, **Hy, **rx, **ry;
	PetscReal      HnMinus1, u, v;
	PetscReal      uNorth, uEast, uWest, uSouth;
	PetscReal      vNorth, vEast, vWest, vSouth;
	PetscReal      convectionTerm, diffusionTerm;
	PetscReal      nu = flowDesc->nu;
	PetscReal      alphaExplicit = simParams->alphaExplicit,
	               gamma = simParams->gamma,
	               zeta  = simParams->zeta;
	PetscReal      dt = simParams->dt;

	// copy fluxes to local vectors
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal); CHKERRV(ierr);
	
	ierr = DMCompositeGetAccess(qPack, H,  &HxGlobal, &HyGlobal); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(qPack, rn, &rxGlobal, &ryGlobal); CHKERRV(ierr);
	
	// access local vectors through multi-dimensional pointers
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	
	// x-component
	ierr = DMDAVecGetArray(uda, HxGlobal, &Hx); CHKERRV(ierr);
	ierr = DMDAVecGetArray(uda, rxGlobal, &rx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			// convection term
			u      = qx[j][i]/mesh->dy[j];
			uWest  = 0.5*(u + qx[j][i-1]/mesh->dy[j]);
			uEast  = 0.5*(u + qx[j][i+1]/mesh->dy[j]);
			// first check if the node is adjacent to the -Y or +Y boundaries
			// then check if the boundary condition in the y-direction is periodic
			uSouth = (j > 0)?   0.5*(u + qx[j-1][i]/mesh->dy[j-1]) : (flowDesc->bc[0][YMINUS].type!=PERIODIC)? qx[j-1][i] : 0.5*(u + qx[j-1][i]/mesh->dy[mesh->ny-1]);
			uNorth = (j < N-1)? 0.5*(u + qx[j+1][i]/mesh->dy[j+1]) : (flowDesc->bc[0][YPLUS].type !=PERIODIC)? qx[j+1][i] : 0.5*(u + qx[j+1][i]/mesh->dy[0]);
			vSouth = 0.5*(qy[j-1][i]/dxU[i] + qy[j-1][i+1]/dxU[i+1]);
			vNorth = 0.5*(  qy[j][i]/dxU[i] +   qy[j][i+1]/dxU[i+1]);
			// Hx = d(u^2)/dx + d(uv)/dy
			HnMinus1 = Hx[j][i];
			Hx[j][i] = (uEast*uEast - uWest*uWest)/(0.5*(dxU[i] + dxU[i+1]))
			           + (uNorth*vNorth - uSouth*vSouth)/mesh->dy[j];
			convectionTerm = gamma*Hx[j][i] + zeta*HnMinus1;
			
			// diffusion term
			// reuse the above variable names to calculate the diffusion term
			// their meanings change
			uWest  = qx[j][i-1]/mesh->dy[j];
			uEast  = qx[j][i+1]/mesh->dy[j];
			uSouth = (j > 0)?   qx[j-1][i]/mesh->dy[j-1] : (flowDesc->bc[0][YMINUS].type!=PERIODIC)? qx[j-1][i] : qx[j-1][i]/mesh->dy[mesh->ny-1];
			uNorth = (j < N-1)? qx[j+1][i]/mesh->dy[j+1] : (flowDesc->bc[0][YPLUS].type !=PERIODIC)? qx[j+1][i] : qx[j+1][i]/mesh->dy[0];
			// Dx = d^2(u)/dx^2 + d^2(u)/dy^2
			diffusionTerm = alphaExplicit*nu*(   du2dx2(uWest,  u, uEast,  dxU[i], dxU[i+1])
			                                   + du2dx2(uSouth, u, uNorth, dyU[j], dyU[j+1])
			                                 );
			rx[j][i] = (u/dt - convectionTerm + diffusionTerm);
		}
	}
	ierr = DMDAVecRestoreArray(uda, HxGlobal, &Hx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(uda, rxGlobal, &rx); CHKERRV(ierr);
	
	// y-component
	ierr = DMDAVecGetArray(vda, HyGlobal, &Hy); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, ryGlobal, &ry); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			// convection term
			v      = qy[j][i]/mesh->dx[i];
			vSouth = 0.5*(v + qy[j-1][i]/mesh->dx[i]);
			vNorth = 0.5*(v + qy[j+1][i]/mesh->dx[i]);
			uWest  = 0.5*(qx[j][i-1]/dyV[j] + qx[j+1][i-1]/dyV[j+1]);
			uEast  = 0.5*(qx[j][i]/dyV[j]   + qx[j+1][i]/dyV[j+1]);
			// first check if the node is adjacent to the -X or +X boundaries
			// then check if the boundary condition in the x-direction is periodic
			vWest  = (i > 0)?   0.5*(v + qy[j][i-1]/mesh->dx[i-1]) : (flowDesc->bc[1][XMINUS].type!=PERIODIC)? qy[j][i-1] : 0.5*(v + qy[j][i-1]/mesh->dx[mesh->nx-1]);
			vEast  = (i < M-1)? 0.5*(v + qy[j][i+1]/mesh->dx[i+1]) : (flowDesc->bc[1][XPLUS].type !=PERIODIC)? qy[j][i+1] : 0.5*(v + qy[j][i+1]/mesh->dx[0]);
			// Hx = d(uv)/dx + d(v^2)/dy
			HnMinus1 = Hy[j][i];
			Hy[j][i] = (uEast*vEast - uWest*vWest)/mesh->dx[i]
			           + (vNorth*vNorth - vSouth*vSouth)/(0.5*(dyV[j] + dyV[j+1]));
			convectionTerm = gamma*Hy[j][i] + zeta*HnMinus1;
			
			// diffusion term			
			// reuse the above variable names to calculate the diffusion term
			// their meanings change
			vSouth = qy[j-1][i]/mesh->dx[i];
			vNorth = qy[j+1][i]/mesh->dx[i];
			vWest  = (i > 0)?   qy[j][i-1]/mesh->dx[i-1] : (flowDesc->bc[1][XMINUS].type!=PERIODIC)? qy[j][i-1] : qy[j][i-1]/mesh->dx[mesh->nx-1];
			vEast  = (i < M-1)? qy[j][i+1]/mesh->dx[i+1] : (flowDesc->bc[1][XPLUS].type !=PERIODIC)? qy[j][i+1] : qy[j][i+1]/mesh->dx[0];
			// Dy = d^2(v)/dx^2 + d^2(v)/dy^2
			diffusionTerm = alphaExplicit*nu*(   du2dx2(vWest,  v, vEast,  dxV[i], dxV[i+1])
			                                   + du2dx2(vSouth, v, vNorth, dyV[j], dyV[j+1])
			                                 );
			
			ry[j][i] = (v/dt - convectionTerm + diffusionTerm);
		}
	}
	ierr = DMDAVecRestoreArray(vda, HyGlobal, &Hy); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, ryGlobal, &ry); CHKERRV(ierr);
	
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);
	
	ierr = DMCompositeRestoreAccess(qPack, H,  &HxGlobal, &HyGlobal); CHKERRV(ierr);
	ierr = DMCompositeRestoreAccess(qPack, rn, &rxGlobal, &ryGlobal); CHKERRV(ierr);
}

template<>
void NavierStokesSolver<3>::calculateExplicitTerms()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p, i, j, k, M, N, P;
	Vec            HxGlobal, HyGlobal, HzGlobal;
	Vec            rxGlobal, ryGlobal, rzGlobal;
	PetscReal      ***qx, ***qy, ***qz;
	PetscReal      ***Hx, ***Hy, ***Hz, ***rx, ***ry, ***rz;
	PetscReal      HnMinus1, u, v, w;
	PetscReal      uNorth, uEast, uWest, uSouth, uNadir, uZenith;
	PetscReal      vNorth, vEast, vWest, vSouth, vNadir, vZenith;
	PetscReal      wNorth, wEast, wWest, wSouth, wNadir, wZenith;
	PetscReal      convectionTerm, diffusionTerm;
	PetscReal      nu = flowDesc->nu;
	PetscReal      alphaExplicit = simParams->alphaExplicit,
	               gamma = simParams->gamma,
	               zeta  = simParams->zeta;
	PetscReal      dt = simParams->dt;

	// copy fluxes to local vectors
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal, qzLocal); CHKERRV(ierr);
	
	ierr = DMCompositeGetAccess(qPack, H,  &HxGlobal, &HyGlobal, &HzGlobal); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(qPack, rn, &rxGlobal, &ryGlobal, &rzGlobal); CHKERRV(ierr);
	
	// access local vectors through multi-dimensional pointers
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRV(ierr);
	
	// x-component
	ierr = DMDAVecGetArray(uda, HxGlobal, &Hx); CHKERRV(ierr);
	ierr = DMDAVecGetArray(uda, rxGlobal, &rx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// convection term
				u = qx[k][j][i]/(mesh->dy[j]*mesh->dz[k]);
				// x
				uWest   = 0.5*(u + qx[k][j][i-1]/(mesh->dy[j]*mesh->dz[k]));
				uEast   = 0.5*(u + qx[k][j][i+1]/(mesh->dy[j]*mesh->dz[k]));
				// y
				// first check if the node is adjacent to the -Y or +Y boundaries
				// then check if the boundary condition in the y-direction is periodic
				uSouth  = (j > 0)?   0.5*(u + qx[k][j-1][i]/(mesh->dy[j-1]*mesh->dz[k])) : (flowDesc->bc[0][YMINUS].type!=PERIODIC)? qx[k][j-1][i] : 0.5*(u + qx[k][j-1][i]/(mesh->dy[mesh->ny-1]*mesh->dz[k]));
				uNorth  = (j < N-1)? 0.5*(u + qx[k][j+1][i]/(mesh->dy[j+1]*mesh->dz[k])) : (flowDesc->bc[0][YPLUS].type !=PERIODIC)? qx[k][j+1][i] : 0.5*(u + qx[k][j+1][i]/(mesh->dy[0]*mesh->dz[k]));
				vSouth  = 0.5*(qy[k][j-1][i]/(mesh->dz[k]*dxU[i]) + qy[k][j-1][i+1]/(mesh->dz[k]*dxU[i+1]));
				vNorth  = 0.5*(  qy[k][j][i]/(mesh->dz[k]*dxU[i]) +   qy[k][j][i+1]/(mesh->dz[k]*dxU[i+1]));
				// z
				// first check if the node is adjacent to the -Z or +Z boundaries
				// then check if the boundary condition in the z-direction is periodic
				uNadir  = (k > 0)?   0.5*(u + qx[k-1][j][i]/(mesh->dy[j]*mesh->dz[k-1])) : (flowDesc->bc[0][ZMINUS].type!=PERIODIC)? qx[k-1][j][i] : 0.5*(u + qx[k-1][j][i]/(mesh->dy[j]*mesh->dz[mesh->nz-1]));
				uZenith = (k < P-1)? 0.5*(u + qx[k+1][j][i]/(mesh->dy[j]*mesh->dz[k+1])) : (flowDesc->bc[0][ZPLUS].type !=PERIODIC)? qx[k+1][j][i] : 0.5*(u + qx[k+1][j][i]/(mesh->dy[j]*mesh->dz[0]));
				wNadir  = 0.5*(qz[k-1][j][i]/(dxU[i]*mesh->dy[j]) + qz[k-1][j][i+1]/(dxU[i+1]*mesh->dy[j]));
				wZenith = 0.5*(  qz[k][j][i]/(dxU[i]*mesh->dy[j]) +   qz[k][j][i+1]/(dxU[i+1]*mesh->dy[j]));
				// Hx = d(u^2)/dx + d(uv)/dy + d(uw)/dz
				HnMinus1    = Hx[k][j][i];
				Hx[k][j][i] = (uEast*uEast - uWest*uWest)/(0.5*(dxU[i] + dxU[i+1]))
				              + (uNorth*vNorth - uSouth*vSouth)/mesh->dy[j]
				              + (uZenith*wZenith - uNadir*wNadir)/mesh->dz[k];
				convectionTerm = gamma*Hx[k][j][i] + zeta*HnMinus1;
				
				// diffusion term
				// reuse the above variable names to calculate the diffusion term
				// their meanings change
				uWest   = qx[k][j][i-1]/(mesh->dy[j]*mesh->dz[k]);
				uEast   = qx[k][j][i+1]/(mesh->dy[j]*mesh->dz[k]);
				uSouth  = (j > 0)?   qx[k][j-1][i]/(mesh->dy[j-1]*mesh->dz[k]) : (flowDesc->bc[0][YMINUS].type!=PERIODIC)? qx[k][j-1][i] : qx[k][j-1][i]/(mesh->dy[mesh->ny-1]*mesh->dz[k]);
				uNorth  = (j < N-1)? qx[k][j+1][i]/(mesh->dy[j+1]*mesh->dz[k]) : (flowDesc->bc[0][YPLUS].type !=PERIODIC)? qx[k][j+1][i] : qx[k][j+1][i]/(mesh->dy[0]*mesh->dz[k]);
				uNadir  = (k > 0)?   qx[k-1][j][i]/(mesh->dy[j]*mesh->dz[k-1]) : (flowDesc->bc[0][ZMINUS].type!=PERIODIC)? qx[k-1][j][i] : qx[k-1][j][i]/(mesh->dy[j]*mesh->dz[mesh->nz-1]);
				uZenith = (k < P-1)? qx[k+1][j][i]/(mesh->dy[j]*mesh->dz[k+1]) : (flowDesc->bc[0][ZPLUS].type !=PERIODIC)? qx[k+1][j][i] : qx[k+1][j][i]/(mesh->dy[j]*mesh->dz[0]);
				// Dx = d^2(u)/dx^2 + d^2(u)/dy^2 + d^2(u)/dz^2
				diffusionTerm = alphaExplicit*nu*(   du2dx2(uWest,  u, uEast,   dxU[i], dxU[i+1])
				                                   + du2dx2(uSouth, u, uNorth,  dyU[j], dyU[j+1])
				                                   + du2dx2(uNadir, u, uZenith, dzU[k], dzU[k+1])
				                                 );
				rx[k][j][i] = (u/dt - convectionTerm + diffusionTerm);
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, HxGlobal, &Hx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(uda, rxGlobal, &rx); CHKERRV(ierr);

	// y-component
	ierr = DMDAVecGetArray(vda, HyGlobal, &Hy); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, ryGlobal, &ry); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(vda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// convection term
				v = qy[k][j][i]/(mesh->dz[k]*mesh->dx[i]);
				// x
				// first check if the node is adjacent to the -X or +X boundaries
				// then check if the boundary condition in the x-direction is periodic
				vWest   = (i > 0)?   0.5*(v + qy[k][j][i-1]/(mesh->dz[k]*mesh->dx[i-1])) : (flowDesc->bc[1][XMINUS].type!=PERIODIC)? qy[k][j][i-1] : 0.5*(u + qy[k][j][i-1]/(mesh->dz[k]*mesh->dx[mesh->nx-1]));
				vEast   = (i < M-1)? 0.5*(v + qy[k][j][i+1]/(mesh->dz[k]*mesh->dx[i+1])) : (flowDesc->bc[1][XPLUS].type !=PERIODIC)? qy[k][j][i+1] : 0.5*(u + qy[k][j][i+1]/(mesh->dz[k]*mesh->dx[0]));
				uWest   = 0.5*(qx[k][j][i-1]/(dyV[j]*mesh->dz[k]) + qx[k][j+1][i-1]/(dyV[j+1]*mesh->dz[k]));
				uEast   = 0.5*(  qx[k][j][i]/(dyV[j]*mesh->dz[k]) +   qx[k][j+1][i]/(dyV[j+1]*mesh->dz[k]));
				// y
				vSouth  = 0.5*(v + qy[k][j-1][i]/(mesh->dz[k]*mesh->dx[i]));
				vNorth  = 0.5*(v + qy[k][j+1][i]/(mesh->dz[k]*mesh->dx[i]));
				// z
				// first check if the node is adjacent to the -Z or +Z boundaries
				// then check if the boundary condition in the z-direction is periodic
				vNadir  = (k > 0)?   0.5*(v + qy[k-1][j][i]/(mesh->dz[k-1]*mesh->dx[i])) : (flowDesc->bc[1][ZMINUS].type!=PERIODIC)? qy[k-1][j][i] : 0.5*(u + qy[k-1][j][i]/(mesh->dz[mesh->nz-1]*mesh->dx[i]));
				vZenith = (k < P-1)? 0.5*(v + qy[k+1][j][i]/(mesh->dz[k+1]*mesh->dx[i])) : (flowDesc->bc[1][ZPLUS].type !=PERIODIC)? qy[k+1][j][i] : 0.5*(u + qy[k+1][j][i]/(mesh->dz[0]*mesh->dx[i]));
				wNadir  = 0.5*(qz[k-1][j][i]/(mesh->dx[i]*dyV[j]) + qz[k-1][j+1][i]/(mesh->dx[i]*dyV[j+1]));
				wZenith = 0.5*(  qz[k][j][i]/(mesh->dx[i]*dyV[j]) +   qz[k][j+1][i]/(mesh->dx[i]*dyV[j+1]));
				// Hx = d(vu)/dx + d(v^2)/dy + d(vw)/dz
				HnMinus1 = Hy[k][j][i];
				Hy[k][j][i] = (vEast*uEast - vWest*uWest)/mesh->dx[i]
				              + (vNorth*vNorth - vSouth*vSouth)/(0.5*(dyV[j] + dyV[j+1]))
				              + (vZenith*wZenith - vNadir*wNadir)/mesh->dz[k];
				convectionTerm = gamma*Hy[k][j][i] + zeta*HnMinus1;

				// diffusion term			
				// reuse the above variables to calculate the diffusion term
				// their meanings change
				vWest   = (i > 0)?   qy[k][j][i-1]/(mesh->dz[k]*mesh->dx[i-1]) : (flowDesc->bc[1][XMINUS].type!=PERIODIC)? qy[k][j][i-1] : qy[k][j][i-1]/(mesh->dz[k]*mesh->dx[mesh->nx-1]);
				vEast   = (i < M-1)? qy[k][j][i+1]/(mesh->dz[k]*mesh->dx[i+1]) : (flowDesc->bc[1][XPLUS].type !=PERIODIC)? qy[k][j][i+1] : qy[k][j][i+1]/(mesh->dz[k]*mesh->dx[0]);
				vSouth  = qy[k][j-1][i]/(mesh->dz[k]*mesh->dx[i]);
				vNorth  = qy[k][j+1][i]/(mesh->dz[k]*mesh->dx[i]);
				vNadir  = (k > 0)?   qy[k-1][j][i]/(mesh->dz[k-1]*mesh->dx[i]) : (flowDesc->bc[1][ZMINUS].type!=PERIODIC)? qy[k-1][j][i] : qy[k-1][j][i]/(mesh->dz[mesh->nz-1]*mesh->dx[i]);
				vZenith = (k < P-1)? qy[k+1][j][i]/(mesh->dz[k+1]*mesh->dx[i]) : (flowDesc->bc[1][ZPLUS].type !=PERIODIC)? qy[k+1][j][i] : qy[k+1][j][i]/(mesh->dz[0]*mesh->dx[i]);
				// Dy = d^2(v)/dx^2 + d^2(v)/dy^2 + d^2(v)/dz^2
				diffusionTerm = alphaExplicit*nu*(   du2dx2(vWest,  v, vEast,   dxV[i], dxV[i+1])
				                                   + du2dx2(vSouth, v, vNorth,  dyV[j], dyV[j+1])
				                                   + du2dx2(vNadir, v, vZenith, dzV[k], dzV[k+1])
				                                 );
				ry[k][j][i] = (v/dt - convectionTerm + diffusionTerm);
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, HyGlobal, &Hy); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, ryGlobal, &ry); CHKERRV(ierr);

	// z-component
	ierr = DMDAVecGetArray(wda, HzGlobal, &Hz); CHKERRV(ierr);
	ierr = DMDAVecGetArray(wda, rzGlobal, &rz); CHKERRV(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(wda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// convection term
				w = qz[k][j][i]/(mesh->dx[i]*mesh->dy[j]);
				// x
				// first check if the node is adjacent to the -X or +X boundaries
				// then check if the boundary condition in the x-direction is periodic
				wWest   = (i > 0)?   0.5*(w + qz[k][j][i-1]/(mesh->dx[i-1]*mesh->dy[j])) : (flowDesc->bc[2][XMINUS].type!=PERIODIC)? qz[k][j][i-1] : 0.5*(u + qz[k][j][i-1]/(mesh->dx[mesh->nx-1]*mesh->dy[j]));
				wEast   = (i < M-1)? 0.5*(w + qz[k][j][i+1]/(mesh->dx[i+1]*mesh->dy[j])) : (flowDesc->bc[2][XPLUS].type !=PERIODIC)? qz[k][j][i+1] : 0.5*(u + qz[k][j][i+1]/(mesh->dx[0]*mesh->dy[j]));
				uWest   = 0.5*(qx[k][j][i-1]/(mesh->dy[j]*dzW[k]) + qx[k+1][j][i-1]/(mesh->dy[j]*dzW[k+1]));
				uEast   = 0.5*(  qx[k][j][i]/(mesh->dy[j]*dzW[k]) +   qx[k+1][j][i]/(mesh->dy[j]*dzW[k+1]));
				// y
				// first check if the node is adjacent to the -Y or +Y boundaries
				// then check if the boundary condition in the y-direction is periodic
				wSouth  = (j > 0)?   0.5*(w + qz[k][j-1][i]/(mesh->dx[i]*mesh->dy[j-1])) : (flowDesc->bc[2][YMINUS].type!=PERIODIC)? qz[k][j-1][i] : 0.5*(u + qz[k][j-1][i]/(mesh->dx[i]*mesh->dy[mesh->ny-1]));
				wNorth  = (j < N-1)? 0.5*(w + qz[k][j+1][i]/(mesh->dx[i]*mesh->dy[j+1])) : (flowDesc->bc[2][YPLUS].type !=PERIODIC)? qz[k][j+1][i] : 0.5*(u + qz[k][j+1][i]/(mesh->dx[i]*mesh->dy[0]));
				vSouth  = 0.5*(qy[k][j-1][i]/(dzW[k]*mesh->dx[i]) + qy[k+1][j-1][i]/(dzW[k+1]*mesh->dx[i]));
				vNorth  = 0.5*(  qy[k][j][i]/(dzW[k]*mesh->dx[i]) +   qy[k+1][j][i]/(dzW[k+1]*mesh->dx[i]));
				// z
				wNadir  = 0.5*(w + qz[k-1][j][i]/(mesh->dx[i]*mesh->dy[j]));
				wZenith = 0.5*(w + qz[k+1][j][i]/(mesh->dx[i]*mesh->dy[j]));
				// Hx = d(wu)/dx + d(wv)/dy + d(w^2)/dz
				HnMinus1 = Hz[k][j][i];
				Hz[k][j][i] = (wEast*uEast - wWest*uWest)/mesh->dx[i]
				              + (wNorth*vNorth - wSouth*vSouth)/mesh->dy[j]
				              + (wZenith*wZenith - wNadir*wNadir)/(0.5*(dzW[k] + dzW[k+1]));
				convectionTerm = gamma*Hz[k][j][i] + zeta*HnMinus1;

				// diffusion term			
				// reuse the above variables to calculate the diffusion term
				// their meanings change
				wWest   = (i > 0)?   qz[k][j][i-1]/(mesh->dx[i-1]*mesh->dy[j]) : (flowDesc->bc[2][XMINUS].type!=PERIODIC)? qz[k][j][i-1] : qz[k][j][i-1]/(mesh->dx[mesh->nx-1]*mesh->dy[j]);
				wEast   = (i < M-1)? qz[k][j][i+1]/(mesh->dx[i+1]*mesh->dy[j]) : (flowDesc->bc[2][XPLUS].type !=PERIODIC)? qz[k][j][i+1] : qz[k][j][i+1]/(mesh->dx[0]*mesh->dy[j]);
				wSouth  = (j > 0)?   qz[k][j-1][i]/(mesh->dx[i]*mesh->dy[j-1]) : (flowDesc->bc[2][YMINUS].type!=PERIODIC)? qz[k][j-1][i] : qz[k][j-1][i]/(mesh->dx[i]*mesh->dy[mesh->ny-1]);
				wNorth  = (j < N-1)? qz[k][j+1][i]/(mesh->dx[i]*mesh->dy[j+1]) : (flowDesc->bc[2][YPLUS].type !=PERIODIC)? qz[k][j+1][i] : qz[k][j+1][i]/(mesh->dx[i]*mesh->dy[0]);
				wNadir  = qz[k-1][j][i]/(mesh->dx[i]*mesh->dy[j]);
				wZenith = qz[k+1][j][i]/(mesh->dx[i]*mesh->dy[j]);
				// Dz = d^2(w)/dx^2 + d^2(w)/dy^2 + d^2(w)/dz^2
				diffusionTerm = alphaExplicit*nu*(   du2dx2(wWest,  w, wEast,   dxW[i], dxW[i+1])
				                                   + du2dx2(wSouth, w, wNorth,  dyW[j], dyW[j+1])
				                                   + du2dx2(wNadir, w, wZenith, dzW[k], dzW[k+1])
				                                 );
				rz[k][j][i] = (w/dt - convectionTerm + diffusionTerm);
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, HzGlobal, &Hz); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(wda, rzGlobal, &rz); CHKERRV(ierr);
	
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRV(ierr);
	
	ierr = DMCompositeRestoreAccess(qPack, H,  &HxGlobal, &HyGlobal, &HzGlobal); CHKERRV(ierr);
	ierr = DMCompositeRestoreAccess(qPack, rn, &rxGlobal, &ryGlobal, &rzGlobal); CHKERRV(ierr);
}

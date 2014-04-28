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
	PetscReal      dxMinus, dxPlus, dyMinus, dyPlus;
	PetscReal      convectionTerm, diffusionTerm;
	PetscReal      alphaExplicit = simParams->alphaExplicit,
	               gamma = simParams->gamma,
	               zeta  = simParams->zeta;
	PetscReal      dt = simParams->dt;
	
	ierr = DMCompositeGetAccess(pack, H, &HxGlobal, &HyGlobal); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(pack, rn, &rxGlobal, &ryGlobal); CHKERRV(ierr);
	
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
			
			HnMinus1 = Hx[j][i];
			
			if(i<M-1)
				dxPlus = mesh->dx[i+1];
			else
				dxPlus = (flowDesc->bc[0][XPLUS].type!=PERIODIC)? mesh->dx[i+1] : mesh->dx[0];
			
			u  = qx[j][i]/mesh->dy[j];
			
			uWest  = 0.5*(u + qx[j][i-1]/mesh->dy[j]);
			uEast  = 0.5*(u + qx[j][i+1]/mesh->dy[j]);

			if(j>0)
				uSouth = 0.5*(u + qx[j-1][i]/mesh->dy[j-1]);
			else
				uSouth = (flowDesc->bc[0][YPLUS].type!=PERIODIC)? qx[j-1][i] : 0.5*(u + qx[j-1][i]/mesh->dy[mesh->ny-1]);
			
			if(j<N-1)
				uNorth = 0.5*(u + qx[j+1][i]/mesh->dy[j+1]);
			else
				uNorth = (flowDesc->bc[0][YPLUS].type!=PERIODIC)? qx[j+1][i] : 0.5*(u + qx[j+1][i]/mesh->dy[0]);

			vSouth = 0.5*(qy[j-1][i]/mesh->dx[i] + qy[j-1][i+1]/mesh->dx[i+1]);
			vNorth = 0.5*(  qy[j][i]/mesh->dx[i] +   qy[j][i+1]/mesh->dx[i+1]);
			
			// Hx = d(u^2)/dx + d(uv)/dy
			Hx[j][i] = (uEast*uEast - uWest*uWest)/(0.5*(mesh->dx[i] + dxPlus)) + (uNorth*vNorth - uSouth*vSouth)/mesh->dy[j];
			
			convectionTerm = gamma*Hx[j][i] + zeta*HnMinus1;
			
			// diffusion term
			
			// reuse the above variables to calculate the diffusion term
			// their meanings change
			
			uWest = qx[j][i-1]/mesh->dy[j];
			uEast = qx[j][i+1]/mesh->dy[j];
			
			if(j>0)
			{
				uSouth  = qx[j-1][i]/mesh->dy[j-1];
				dyMinus = 0.5*(mesh->dy[j] + mesh->dy[j-1]);
			}
			else
			{
				uSouth  = (flowDesc->bc[0][YPLUS].type!=PERIODIC)? qx[j-1][i] : qx[j-1][i]/mesh->dy[mesh->ny-1];
				dyMinus = (flowDesc->bc[0][YPLUS].type!=PERIODIC)? 0.5*mesh->dy[j] : 0.5*(mesh->dy[j] + mesh->dy[mesh->ny-1]); 
			}

			if(j<N-1)
			{
				uNorth = qx[j+1][i]/mesh->dy[j+1];
				dyPlus = 0.5*(mesh->dy[j] + mesh->dy[j+1]);
			}
			else
			{
				uNorth = (flowDesc->bc[0][YPLUS].type!=PERIODIC)? qx[j+1][i] : qx[j+1][i]/mesh->dy[0];
				dyPlus = (flowDesc->bc[0][YPLUS].type!=PERIODIC)? 0.5*mesh->dy[j] : 0.5*(mesh->dy[j] + mesh->dy[0]);
			}
			
			diffusionTerm = alphaExplicit*(  du2dx2(uWest, u, uEast, mesh->dx[i], dxPlus)
			                       + du2dx2(uSouth, u, uNorth, dyMinus, dyPlus)
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
			
			HnMinus1 = Hy[j][i];
			
			if(j<N-1)
				dyPlus = mesh->dy[j+1];
			else
				dyPlus = (flowDesc->bc[1][YPLUS].type!=PERIODIC)? mesh->dy[j+1] : mesh->dy[0];
			
			v  = qy[j][i]/mesh->dx[i];
			
			vSouth = 0.5*(v + qy[j-1][i]/mesh->dx[i]);
			vNorth = 0.5*(v + qy[j+1][i]/mesh->dx[i]);
			
			uWest  = 0.5*(qx[j][i-1]/mesh->dy[j] + qx[j+1][i-1]/mesh->dy[j+1]);
			uEast  = 0.5*(  qx[j][i]/mesh->dy[j] +   qx[j+1][i]/mesh->dy[j+1]);
			
			if(i>0)
				vWest  = 0.5*(v + qy[j][i-1]/mesh->dx[i-1]);
			else
				vWest  = (flowDesc->bc[1][XPLUS].type!=PERIODIC)? qy[j][i-1] : 0.5*(v + qy[j][i-1]/mesh->dx[mesh->nx-1]);
			
			if(i<M-1)
				vEast  = 0.5*(v + qy[j][i+1]/mesh->dx[i+1]);
			else
				vEast  = (flowDesc->bc[1][XPLUS].type!=PERIODIC)? qy[j][i+1] : 0.5*(v + qy[j][i+1]/mesh->dx[0]);
			
			// Hx = d(uv)/dx + d(v^2)/dy
			Hy[j][i] = (uEast*vEast - uWest*vWest)/mesh->dx[i] + (vNorth*vNorth - vSouth*vSouth)/(0.5*(mesh->dy[j] + dyPlus));
			
			convectionTerm = gamma*Hy[j][i] + zeta*HnMinus1;
			
			// diffusion term
			
			// reuse the above variables to calculate the diffusion term
			// their meanings change
			
			vSouth = qy[j-1][i]/mesh->dx[i];
			vNorth = qy[j+1][i]/mesh->dx[i];
			
			if(i>0)
			{
				vWest   = qy[j][i-1]/mesh->dx[i-1];
				dxMinus = 0.5*(mesh->dx[i] + mesh->dx[i-1]);
			}
			else
			{
				vWest   = (flowDesc->bc[1][XPLUS].type!=PERIODIC)? qy[j][i-1] : qy[j][i-1]/mesh->dx[mesh->nx-1];
				dxMinus = (flowDesc->bc[1][XPLUS].type!=PERIODIC)? 0.5*mesh->dx[i] : 0.5*(mesh->dx[i] + mesh->dx[mesh->nx-1]);
			}
			
			if(i<M-1)
			{
				vEast  = qy[j][i+1]/mesh->dx[i+1];
				dxPlus = 0.5*(mesh->dx[i] + mesh->dx[i+1]);
			}
			else
			{
				vEast  = (flowDesc->bc[1][XPLUS].type!=PERIODIC)? qy[j][i+1] : qy[j][i+1]/mesh->dx[0];
				dxPlus = (flowDesc->bc[1][XPLUS].type!=PERIODIC)? 0.5*mesh->dx[i] : 0.5*(mesh->dx[i] + mesh->dx[0]);
			}
			
			diffusionTerm = alphaExplicit*(  du2dx2(vWest, v, vEast, dxMinus, dxPlus)
			                       + du2dx2(vSouth, u, vNorth, mesh->dy[j], dyPlus)
			                      );
			
			ry[j][i] = (v/dt - convectionTerm + diffusionTerm);
		}
	}
	ierr = DMDAVecRestoreArray(vda, HyGlobal, &Hy); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, ryGlobal, &ry); CHKERRV(ierr);
	
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);
	
	ierr = DMCompositeRestoreAccess(pack, H, &HxGlobal, &HyGlobal); CHKERRV(ierr);
	ierr = DMCompositeRestoreAccess(pack, rn, &rxGlobal, &ryGlobal); CHKERRV(ierr);
}

template<>
void NavierStokesSolver<3>::calculateExplicitTerms()
{
}

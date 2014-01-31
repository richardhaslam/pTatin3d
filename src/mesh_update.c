/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      mesh_update.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#include "stdio.h"
#include "stdlib.h"
#include "petsc.h"
#include "petscdm.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_element_q2p1.h"
#include "quadrature.h"
#include "element_utils_q2.h"
#include "dmda_update_coords.h"
#include "element_utils_q1.h"
#include "dmda_remesh.h"
#include "dmda_iterator.h"
#include "dmda_redundant.h"
#include "material_point_point_location.h"

#include "mesh_update.h"


#undef __FUNCT__
#define __FUNCT__ "DMDABilinearizeQ2Elements"
PetscErrorCode DMDABilinearizeQ2Elements(DM dau)
{
	DM cda;
	Vec gcoords;
	PetscScalar *LA_gcoords;
	PetscInt nel,nen,e,n,k,ii,jj,kk;
	const PetscInt *elnidx;
	PetscScalar elcoordsQ2[3*Q2_NODES_PER_EL_3D];
	PetscScalar elcoordsQ1[3*Q1_NODES_PER_EL_3D],Ni[Q1_NODES_PER_EL_3D];
	PetscInt cnt;
	PetscScalar xi_nodal_coordsQ2[3*Q2_NODES_PER_EL_3D],x_new[3];
	PetscInt vel_el_lidx[3*Q2_NODES_PER_EL_3D];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* define some xi coords to interpolate to - these correspond to the nodal location of the q2 basis functions */
	/* this loop should be ordered the same way as the natural nodal ordering of the Q2 nodes in element space */
	cnt = 0;
	for (kk=0; kk<3; kk++) {
		for (jj=0; jj<3; jj++) {
			for (ii=0; ii<3; ii++) {
				xi_nodal_coordsQ2[3*cnt+0] = -1.0 + (double)ii;
				xi_nodal_coordsQ2[3*cnt+1] = -1.0 + (double)jj;
				xi_nodal_coordsQ2[3*cnt+2] = -1.0 + (double)kk;
				cnt++;
			}
		}
	}	
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dau,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dau,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen,&elnidx);CHKERRQ(ierr);
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoordsQ2,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		for (k=0; k<3; k++) {
			elcoordsQ1[3*0+k] = elcoordsQ2[3*0+k];
			elcoordsQ1[3*1+k] = elcoordsQ2[3*2+k];
			elcoordsQ1[3*2+k] = elcoordsQ2[3*6+k];
			elcoordsQ1[3*3+k] = elcoordsQ2[3*8+k];
			
			elcoordsQ1[3*4+k] = elcoordsQ2[3*18+k];
			elcoordsQ1[3*5+k] = elcoordsQ2[3*20+k];
			elcoordsQ1[3*6+k] = elcoordsQ2[3*24+k];
			elcoordsQ1[3*7+k] = elcoordsQ2[3*26+k];
		}
		
		/* for each interior point */
		for (k=0; k<27; k++) {
			P3D_ConstructNi_Q1_3D(&xi_nodal_coordsQ2[3*k],Ni);
			
			/* inpterpolate */
			x_new[0] = 0.0;
			x_new[1] = 0.0;
			x_new[2] = 0.0;
			for (n=0; n<Q1_NODES_PER_EL_3D; n++) {
				x_new[0] += elcoordsQ1[3*n+0] * Ni[n];
				x_new[1] += elcoordsQ1[3*n+1] * Ni[n];
				x_new[2] += elcoordsQ1[3*n+2] * Ni[n];
			}
			
			elcoordsQ2[3*k+0] = x_new[0];
			elcoordsQ2[3*k+1] = x_new[1];
			elcoordsQ2[3*k+2] = x_new[2];
			//printf("e-%d: %1.4e  %1.4e  %1.4e \n", e,x_new[0],x_new[1],x_new[2] );
		}
		
		/* push modification */
		ierr = DMDASetValuesLocalStencil_InsertValues_Stokes_Velocity(LA_gcoords, vel_el_lidx,elcoordsQ2);CHKERRQ(ierr);
		
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	/* send ghostes values into global vector */
	ierr = DMDASetCoordinatesFromLocalVector(dau,gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateMeshGeometry_FullLagrangian"
PetscErrorCode UpdateMeshGeometry_FullLagrangian(DM dav,Vec velocity,PetscReal step)
{
	Vec            coordinates;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
	ierr = VecAXPY(coordinates,step,velocity);CHKERRQ(ierr); /* x = x + dt.vel_advect_mesh */
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh"
PetscErrorCode UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(DM dav,Vec velocity,PetscReal step)
{
	Vec            velocity_ale;
	Vec            coordinates;
	PetscInt       M,N,P;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

	ierr = DMGetGlobalVector(dav,&velocity_ale);CHKERRQ(ierr);
	ierr = VecCopy(velocity,velocity_ale);CHKERRQ(ierr);
	
	ierr = VecStrideSet(velocity_ale,0,0.0);CHKERRQ(ierr); /* zero x component */
	ierr = VecStrideSet(velocity_ale,2,0.0);CHKERRQ(ierr); /* zero y component */

	ierr = DMDAVecTraverseIJK(dav,velocity_ale,1,DMDAVecTraverseIJK_ZeroInteriorMinusNmax,(void*)&N);CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
	ierr = VecAXPY(coordinates,step,velocity_ale);CHKERRQ(ierr); /* x = x + dt.vel_ale */
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	
	ierr = DMRestoreGlobalVector(dav,&velocity_ale);CHKERRQ(ierr);

	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,N);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateMeshGeometry_FullLagrangianWithVerticalSurfaceRemesh"
PetscErrorCode UpdateMeshGeometry_FullLagrangianWithVerticalSurfaceRemesh(DM dav,Vec velocity,PetscReal step)
{
	PetscInt       M,N,P;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* take a full lagrangian step */
	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);
	
	{
		PetscReal dx;
		PetscReal gmin[3],gmax[3];
		PetscInt i,k,si,sj,sk,nx,ny,nz;
		Vec coordinates;
		DM cda;
		DMDACoor3d ***LA_coords;
		
		
		ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
		ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
		dx = (gmax[0]-gmin[0])/( (PetscReal)((M-1)) );
		
		/* resample bottom nodes so they are equi-distance in x */
		ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
		
		ierr = DMDAGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
		ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);

		if (sj == 0) { /* bottom layer nodes */
			for (k=sk; k<sk+nz; k++) {
				for (i=si; i<si+nx; i++) {
					LA_coords[k][0][i].x = gmin[0] + i * dx;
				}
			}
		}
		
		ierr = DMDAVecRestoreArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);

		ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	}
	
	
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,N);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateMeshGeometry_DecoupledHorizontalVerticalMeshMovement"
PetscErrorCode UpdateMeshGeometry_DecoupledHorizontalVerticalMeshMovement(DM dav,Vec velocity,PetscReal step)
{
	PetscInt       M,N,P;
	Vec            velocity_ale;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMGetGlobalVector(dav,&velocity_ale);CHKERRQ(ierr);
	ierr = VecCopy(velocity,velocity_ale);CHKERRQ(ierr);

	/* kill y componenet */
	ierr = VecStrideSet(velocity_ale,1,0.0);CHKERRQ(ierr); /* zero y component */
	/* take a full lagrangian step in x-z  */
	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity_ale,step);CHKERRQ(ierr);

	
	/* resample bottom nodes so they are equi-distance in x-z */
	{
		PetscReal dx,dz;
		PetscReal gmin[3],gmax[3];
		PetscInt i,j,k,si,sj,sk,nx,ny,nz;
		Vec coordinates;
		DM cda;
		DMDACoor3d ***LA_coords;
		
		
		ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
		ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
		dx = (gmax[0]-gmin[0])/( (PetscReal)((M-1)) );
		dz = (gmax[2]-gmin[2])/( (PetscReal)((P-1)) );
		
		ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
		
		ierr = DMDAGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
		ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
		
		if (sj == 0) { /* bottom layer nodes */
			j = 0;
			for (k=sk; k<sk+nz; k++) {
				for (i=si; i<si+nx; i++) {
					LA_coords[k][j][i].x = gmin[0] + i * dx;
					LA_coords[k][j][i].z = gmin[2] + k * dz;
				}
			}
		}

		if (sj + ny == N) { /* top layer nodes */
			j = N - 1;
			for (k=sk; k<sk+nz; k++) {
				for (i=si; i<si+nx; i++) {
					LA_coords[k][j][i].x = gmin[0] + i * dx;
					LA_coords[k][j][i].z = gmin[2] + k * dz;
				}
			}
		}
		
		
		ierr = DMDAVecRestoreArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
		
		ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	}
	

	ierr = VecCopy(velocity,velocity_ale);CHKERRQ(ierr);
	ierr = VecStrideSet(velocity_ale,0,0.0);CHKERRQ(ierr); /* zero x component */
	ierr = VecStrideSet(velocity_ale,2,0.0);CHKERRQ(ierr); /* zero z component */
	

	/* take a full lagrangian step in y  */
	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity_ale,step);CHKERRQ(ierr);
	
	/* move exterior based on interior */
	ierr = DMDARemeshJMAX_UpdateHeightsFromInterior(dav);CHKERRQ(ierr);
	
	/* clean up spacing inside */
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,N);CHKERRQ(ierr);
	
	
	ierr = DMRestoreGlobalVector(dav,&velocity_ale);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX_InterpolateMSurfToVolSurf"
PetscErrorCode UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX_InterpolateMSurfToVolSurf(DM dm_msurf,DM dm_vol)
{
	PetscErrorCode ierr;
	PetscInt       M,N,P,si,sj,sk,nx,ny,nz,ni,nj,nk,msurf_mi,msurf_mk;
	PetscInt       i,j,k,n,ei,ek;
	PetscReal      H_p,el_H[Q2_NODES_PER_EL_2D],Ni[Q2_NODES_PER_EL_2D];
	Vec            msurf_coords,v_coords;
	PetscScalar    *LA_msurf_coords;
	DMDACoor3d     ***LA_v_coords;
	PetscReal      *msurf_coords_2d;
	PetscInt       *msurf_element,elcnt,*el,nid[Q2_NODES_PER_EL_2D];
	DM             v_cda;
	
	PetscFunctionBegin;
	/* generate an element->node mapping for dm_surf for Q2 */
	ierr = DMDAGetInfo(dm_msurf,0,&ni,&nj,&nk,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*ni*nk*Q2_NODES_PER_EL_2D,&msurf_element);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscReal)*2*ni*nk,&msurf_coords_2d);CHKERRQ(ierr);
	
	msurf_mi = (ni - 1)/2;
	msurf_mk = (nk - 1)/2;
	elcnt = 0;
	for (ek=0; ek<msurf_mk; ek++) {
		k = 2*ek;
		for (ei=0; ei<msurf_mi; ei++) {
			i = 2*ei;
			
			el = &msurf_element[Q2_NODES_PER_EL_2D*elcnt];
			
			nid[0] = (i  ) + (k  ) *ni;
			nid[1] = (i+1) + (k  ) *ni;
			nid[2] = (i+2) + (k  ) *ni;
			
			nid[3] = (i  ) + (k+1) *ni;
			nid[4] = (i+1) + (k+1) *ni;
			nid[5] = (i+2) + (k+1) *ni;
			
			nid[6] = (i  ) + (k+2) *ni;
			nid[7] = (i+1) + (k+2) *ni;
			nid[8] = (i+2) + (k+2) *ni;
			
			for (n=0; n<Q2_NODES_PER_EL_2D; n++) {
				if (nid[n] > ni*nj*nk) { 
					SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of global nodes");
				}
				el[n] = nid[n];
			}
			elcnt++;
			
		}
	}
	
	ierr = DMDAGetCoordinates(dm_msurf,&msurf_coords);CHKERRQ(ierr);
	ierr = VecGetArray(msurf_coords,&LA_msurf_coords);CHKERRQ(ierr);
	
	for (k=0; k<ni*nk; k++) {
		msurf_coords_2d[2*k+0] = LA_msurf_coords[3*k+0]; // x
		msurf_coords_2d[2*k+1] = LA_msurf_coords[3*k+2]; // z
	}


	ierr = DMDAGetInfo(dm_vol,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm_vol,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinates(dm_vol,&v_coords);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(dm_vol,&v_cda);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(v_cda,v_coords,&LA_v_coords);CHKERRQ(ierr);
	
	/* a) Insert coords of vol->surf into material point structure for interpolation purposes */
	/* b) Map vol->surf.x -> mp.x; vol->surf.z -> mp.y */
	/* c) perform point location with msurf */
	/* d) perform interpolation from msurf -> vol->surf.y */

	if (sj + ny == N) {
		j = N - 1;

		for (k=sk; k<sk+nz; k++) {
			for (i=si; i<si+nx; i++) {
				MPntStd mp;
				
				/* a),b) */
				mp.coor[0] = LA_v_coords[k][j][i].x;
				mp.coor[1] = LA_v_coords[k][j][i].z;
				mp.coor[2] = 0.0;
				/* c) */
				InverseMappingDomain_2dQ2(1.0e-10,20,PETSC_FALSE,PETSC_FALSE,
																	msurf_coords_2d,msurf_mi,msurf_mk,msurf_element,
																	1,&mp);
				/*
				if (mp.wil == -1) {
					PetscPrintf(PETSC_COMM_SELF,"Error: point %d not found\n",n);
					PetscPrintf(PETSC_COMM_SELF,"  pos(x,z) %1.8e %1.8e\n",mp.coor[0],mp.coor[1]);
					SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"  location of node in interpolated surface could not be determined");
				}
				 */

				if (mp.wil == -1) {
					PetscReal sep2,min_sep2 = 1.0e32;
					PetscInt index;
					
					//PetscPrintf(PETSC_COMM_SELF,"Warning: point %d not found\n",(i-si)+(k-sk)*nx);
					//PetscPrintf(PETSC_COMM_SELF,"  pos(x,z) %1.8e %1.8e\n",mp.coor[0],mp.coor[1]);

					index = -1;
					sep2 = 0.0;
					for (n=0; n<ni*nk; n++) {
						sep2 = sep2 + (msurf_coords_2d[2*n+0]-mp.coor[0])*(msurf_coords_2d[2*n+0]-mp.coor[0]);
						sep2 = sep2 + (msurf_coords_2d[2*n+1]-mp.coor[1])*(msurf_coords_2d[2*n+1]-mp.coor[1]);
						if (sep2 < min_sep2) {
							min_sep2 = sep2;
							index = n;
						}
					}
					
					H_p = LA_msurf_coords[3*index+1];
					
				} else {
					
					/* d) */
					/* get element heights */
					for (n=0; n<Q2_NODES_PER_EL_2D; n++) {
						PetscInt nidx;
						
						nidx = msurf_element[ Q2_NODES_PER_EL_2D*mp.wil + n ];
						el_H[n] = LA_msurf_coords[3*nidx+1]; /* fetch the J component from the coordinate array */
					}
					
					/* evaluate basis and interplate */
					P3D_ConstructNi_Q2_2D(mp.xi,Ni);
					H_p = 0.0;
					for (n=0; n<Q2_NODES_PER_EL_2D; n++) {
						H_p += Ni[n] * el_H[n];
					}
					
				}
				
				/* assign interplated height */
				LA_v_coords[k][j][i].y = H_p;
			}
		}
		
	}
	
	ierr = DMDAVecRestoreArray(v_cda,v_coords,&LA_v_coords);CHKERRQ(ierr);
	ierr = VecRestoreArray(msurf_coords,&LA_msurf_coords);CHKERRQ(ierr);
	ierr = PetscFree(msurf_element);CHKERRQ(ierr);
	ierr = PetscFree(msurf_coords_2d);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX"
PetscErrorCode UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX(DM dav,Vec vel_fluid,Vec vel_mesh,PetscReal step)
{
	PetscInt       M,N,P;
	PetscErrorCode ierr;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz;
	Vec coordinates;
	DM cda,dm_surf_max_j;
	DMDACoor3d ***LA_coords,***LA_vel;
	PetscInt surf_si,surf_ei,surf_sk,surf_ek;
	
	
	PetscFunctionBegin;
	
	/* [a] Surface advection only: take a full lagrangian step */
	//PetscPrintf(PETSC_COMM_WORLD,"[a] --- \n");
	ierr = DMDAGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dav,vel_fluid,&LA_vel);CHKERRQ(ierr);

	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	if (sj + ny == N) { /* top layer nodes */
		j = N - 1;
		for (k=sk; k<sk+nz; k++) {
			for (i=si; i<si+nx; i++) {
				LA_coords[k][j][i].x = LA_coords[k][j][i].x + step * LA_vel[k][j][i].x;
				LA_coords[k][j][i].y = LA_coords[k][j][i].y + step * LA_vel[k][j][i].y;
				LA_coords[k][j][i].z = LA_coords[k][j][i].z + step * LA_vel[k][j][i].z;
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(dav,vel_fluid,&LA_vel);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
	
	/* [b] extract the advected surface */
	//PetscPrintf(PETSC_COMM_WORLD,"[b] --- \n");
/*
	{
		surf_si = si - 4;      if (surf_si < 0) { surf_si = 0; }
		surf_ei = si + nx + 4; if (surf_ei > M) { surf_ei = M; }
		surf_sk = sk - 4;      if (surf_sk < 0) { surf_sk = 0; }
		surf_ek = sk + nz + 4; if (surf_ek > P) { surf_ek = P; }
		
	}
*/
	{
		PetscInt esi,esj,esk,mx,my,mz,mo;
		
		ierr = DMDAGetCornersElementQ2(dav,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
		mo = 2; /* number of overlapping elements */
		surf_si = esi - 2*mo;      if (surf_si < 0) { surf_si = 0; }
		surf_ei = esi + nx + 2*mo; if (surf_ei > M) { surf_ei = M; }
		surf_sk = esk - 2*mo;      if (surf_sk < 0) { surf_sk = 0; }
		surf_ek = esk + nz + 2*mo; if (surf_ek > P) { surf_ek = P; }
	}
	ierr = DMDACreate3dRedundant(dav,surf_si,surf_ei,N-1,N,surf_sk,surf_ek, 1, &dm_surf_max_j);CHKERRQ(ierr);


	/* push advected surface back - it's a local operation so there are no concerns regarding efficiency */
	ierr = DMDAGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dav,vel_fluid,&LA_vel);CHKERRQ(ierr);
	if (sj + ny == N) { /* top layer nodes */
		j = N - 1;
		for (k=sk; k<sk+nz; k++) {
			for (i=si; i<si+nx; i++) {
				LA_coords[k][j][i].x = LA_coords[k][j][i].x - step * LA_vel[k][j][i].x;
				LA_coords[k][j][i].y = LA_coords[k][j][i].y - step * LA_vel[k][j][i].y;
				LA_coords[k][j][i].z = LA_coords[k][j][i].z - step * LA_vel[k][j][i].z;
			}
		}
	}
	ierr = DMDAVecRestoreArray(dav,vel_fluid,&LA_vel);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
	
	
	/* 
	 [c] Update other mesh directions, if there is a vector defined for moving the mesh
	 */
	/*
	// x = x + dt.vx
	// z = z + dt.vz
	if (vel_mesh) {
		Vec vel_mesh_copy;
		
		ierr = VecDuplicate(vel_mesh,&vel_mesh_copy);CHKERRQ(ierr);
		ierr = VecCopy(vel_mesh,vel_mesh_copy);CHKERRQ(ierr);
		ierr = VecStrideSet(vel_mesh_copy,1,0.0);CHKERRQ(ierr); // zero y component //
		
		ierr = DMDAGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
		ierr = VecAXPY(coordinates,step,vel_mesh_copy);CHKERRQ(ierr); // x = x + dt.vel_mesh[i]; z = z + dt.vel_mesh[k] //
		ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
		
		ierr = VecDestroy(&vel_mesh_copy);CHKERRQ(ierr);
	}	
	*/
	if (vel_mesh) {
		//PetscPrintf(PETSC_COMM_WORLD,"[c] --- \n");
		ierr = DMDAGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
		ierr = VecAXPY(coordinates,step,vel_mesh);CHKERRQ(ierr); // (x,y,z) = (x,y,z) + dt.vel_mesh(i,j,k) //
		ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	}		
	
	
	/* [d] interpolate surface geometry onto new mesh coordinates */
	//PetscPrintf(PETSC_COMM_WORLD,"[d] --- \n");
	ierr = UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX_InterpolateMSurfToVolSurf(dm_surf_max_j,dav);CHKERRQ(ierr);
	
	
	/* [e] clean up between top-bottom of mesh */
	//PetscPrintf(PETSC_COMM_WORLD,"[e] --- \n");
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,N);CHKERRQ(ierr);
	
	
	ierr = DMDestroy(&dm_surf_max_j);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateMeshGeometry_ComputeSurfaceCourantTimestep"
PetscErrorCode UpdateMeshGeometry_ComputeSurfaceCourantTimestep(DM dav,Vec velocity,PetscReal vert_displacement_max,PetscReal *step)
{
	PetscReal vy;
	DMDACoor3d ***LA_vel;
	PetscReal dt,dt_min_local,dt_min;
	PetscInt M,N,P,i,k,si,sj,sk,nx,ny,nz;
	MPI_Comm comm;
	PetscErrorCode ierr;

	dt_min_local = 1.0e32;
	
	ierr = DMDAGetInfo(dav,0,&M,&N,&P, 0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dav,velocity,&LA_vel);CHKERRQ(ierr);

	if (sj + ny == N) {
		
		for (k=sk; k<sk+nz; k++) {
			for (i=si; i<si+nx; i++) {
				vy = PetscRealPart(LA_vel[k][N-1][i].y);
				
				dt = vert_displacement_max / PetscAbsReal(vy);
				if (dt < dt_min_local) {
					dt_min_local = dt;
				}
			}
		}
		
	}
	ierr = DMDAVecRestoreArray(dav,velocity,&LA_vel);CHKERRQ(ierr);

	ierr = PetscObjectGetComm((PetscObject)dav,&comm);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&dt_min_local,&dt_min,1,MPIU_REAL,MPI_MIN,comm);CHKERRQ(ierr);
	
	*step = dt_min;
	
	PetscFunctionReturn(0);
}

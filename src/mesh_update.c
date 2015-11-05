/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   mesh_update.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/


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
#include "dmda_view_petscvtk.h"
#include "dmda_element_q1.h"
#include "model_utils.h"

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
	ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
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
	
	ierr = DMGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
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
	
	ierr = DMGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
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
		
		ierr = DMGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
		ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
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
		
		ierr = DMGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
		ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
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
	
	ierr = DMGetCoordinates(dm_msurf,&msurf_coords);CHKERRQ(ierr);
	ierr = VecGetArray(msurf_coords,&LA_msurf_coords);CHKERRQ(ierr);
	
	for (k=0; k<ni*nk; k++) {
		msurf_coords_2d[2*k+0] = LA_msurf_coords[3*k+0]; // x
		msurf_coords_2d[2*k+1] = LA_msurf_coords[3*k+2]; // z
	}


	ierr = DMDAGetInfo(dm_vol,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm_vol,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	
	ierr = DMGetCoordinates(dm_vol,&v_coords);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dm_vol,&v_cda);CHKERRQ(ierr);
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
					PetscPrintf(PETSC_COMM_SELF,"Error: point %D not found\n",n);
					PetscPrintf(PETSC_COMM_SELF,"  pos(x,z) %1.8e %1.8e\n",mp.coor[0],mp.coor[1]);
					SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"  location of node in interpolated surface could not be determined");
				}
				 */

				if (mp.wil == -1) {
					PetscReal sep2,min_sep2 = 1.0e32;
					PetscInt index;
					
					//PetscPrintf(PETSC_COMM_SELF,"Warning: point %D not found\n",(i-si)+(k-sk)*nx);
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
	ierr = DMGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
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
	ierr = DMGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
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
		
		ierr = DMGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
		ierr = VecAXPY(coordinates,step,vel_mesh_copy);CHKERRQ(ierr); // x = x + dt.vel_mesh[i]; z = z + dt.vel_mesh[k] //
		ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
		
		ierr = VecDestroy(&vel_mesh_copy);CHKERRQ(ierr);
	}	
	*/
	if (vel_mesh) {
		//PetscPrintf(PETSC_COMM_WORLD,"[c] --- \n");
		ierr = DMGetCoordinates(dav,&coordinates);CHKERRQ(ierr);
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
	ierr = MPI_Allreduce(&dt_min_local,&dt_min,1,MPIU_REAL,MPIU_MIN,comm);CHKERRQ(ierr);
	
	*step = dt_min;
	
	PetscFunctionReturn(0);
}

/*

 Applies a diffusion operator on the surface, j = JMAX-1
 
 dh(x,z,t)/dt = -kappa \nabla . (\nabla h(x,z,t))
 
 where h(x,z) is the height of the nodes along the surface j = JMAX-1.
 \nabla = ( \partial / \partial x, \partial / \partial z )

 The solution is discretised in space using Q1 finite elements. 
 Solution is updated using a first order explicit time integration.
 
 The user must specify at least two dirichlet boundary conditions.
 These will fix the hieght of the surface h(x,z) at walls {east,west,front,back}
 
 M h^{k+1} = M h^k + dt K h^k
 which is approximated as
 h^{k+1} = h^k + dt inv(diag(M)) K h^k
 One point quadrature rule is used to evaluate all matrices.
 
 
*/
#undef __FUNCT__
#define __FUNCT__ "UpdateMeshGeometry_ApplyDiffusionJMAX"
PetscErrorCode UpdateMeshGeometry_ApplyDiffusionJMAX(DM dav,PetscReal diffusivity,PetscReal timespan,
                                                     PetscBool dirichlet_east,PetscBool dirichlet_west,PetscBool dirichlet_front,PetscBool dirichlet_back,PetscBool only_update_surface)
{
	PetscErrorCode ierr;
    DM daH,cda;
    PetscInt nsteps,s,nM;
    PetscReal dt_explicit;
	Vec Hinit,H,local_H,coords,gcoords,diagM,rhs,local_rhs,local_M;
    PetscScalar *LA_H,*LA_gcoords,*LA_rhs,*LA_M;
    PetscInt i,j,k,e,nel,nen,MX,MY,MZ;
    const PetscInt *elnidx;
    PetscInt ge_eqnums[NODES_PER_EL_Q1_3D];
    PetscScalar el_coords[NSD*NODES_PER_EL_Q1_3D];
	PetscScalar el_h[NODES_PER_EL_Q1_3D];
    PetscScalar ***LA_H3,***LA_Hinit3;
    PetscInt NI,NJ,NK,si,sj,sk,ni,nj,nk;
    PetscReal MeshMin[3],MeshMax[3],ds;
    PetscBool view_surface_diffusion = PETSC_FALSE;
    
	PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD,"Applying surface diffision. Dirichlet sides: { east=%D , west=%D , front=%D , back=%D } \n",dirichlet_east,dirichlet_west,dirichlet_front,dirichlet_back);
    if (only_update_surface) {
        PetscPrintf(PETSC_COMM_WORLD,"  [Only applying diffusion operator to JMAX surface] \n");
    }
    if (view_surface_diffusion) {
        PetscPrintf(PETSC_COMM_WORLD,"  [Surface diffusion vizualisation on <debugging>] \n");
    }
    
    /* clone the Q2 DMDA and make a Q1 representation */
    /* duplicate DMDA for a scalar field */
    ierr = DMDACreateOverlappingQ1FromQ2(dav,1,&daH);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(daH,&H);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(daH,&Hinit);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(daH,&diagM);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(daH,&rhs);CHKERRQ(ierr);
    
    ierr = VecGetLocalSize(diagM,&nM);CHKERRQ(ierr);
    
	ierr = DMDAGetInfo(daH,0,&NI,&NJ,&NK,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
    ierr = DMDAGetCorners(daH,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
    ierr = DMDAGetBoundingBox(daH,MeshMin,MeshMax);CHKERRQ(ierr);
    
    ds = (MeshMax[0] - MeshMin[0])/((PetscReal)(NI-1));
    ds = PetscMin(ds, (MeshMax[2] - MeshMin[2])/((PetscReal)(NK-1)));
    PetscPrintf(PETSC_COMM_WORLD,"  [approximate min. element size] ds = %1.4e \n",ds);
    
    /* compute a time step based on x,z cell sizes and a CFL of 0.5 */
    dt_explicit = 0.1 *(ds*ds) / diffusivity;
    PetscPrintf(PETSC_COMM_WORLD,"  [time step] dt_explicit = %1.4e \n",dt_explicit);
    
    nsteps = timespan / dt_explicit;
    nsteps++;
    if (nsteps < 1) { nsteps = 1; }
    dt_explicit = timespan/((PetscReal)nsteps);

    
    PetscPrintf(PETSC_COMM_WORLD,"  [time period of diffusion] timespan = %1.4e \n",timespan);
    PetscPrintf(PETSC_COMM_WORLD,"  [number of explicit time steps] nsteps = %D \n",nsteps);
    PetscPrintf(PETSC_COMM_WORLD,"  [time step size] dt_explicit = %1.4e \n",dt_explicit);
    
    /* extract height */
    ierr = DMGetCoordinateDM(daH,&cda);CHKERRQ(ierr);
    ierr = DMGetCoordinates(daH,&coords);CHKERRQ(ierr);
    ierr = VecStrideGather(coords,1,H,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecCopy(H,Hinit);CHKERRQ(ierr);
    
	/* setup for coords */
    ierr = DMGetCoordinateDM(daH,&cda);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(daH,&gcoords);CHKERRQ(ierr);
    ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

    ierr = DMGetLocalVector(daH,&local_H);CHKERRQ(ierr);
    ierr = DMGetLocalVector(daH,&local_rhs);CHKERRQ(ierr);
    ierr = DMGetLocalVector(daH,&local_M);CHKERRQ(ierr);

    ierr = DMDAGetSizeElementQ2(dav,&MX,&MY,&MZ);CHKERRQ(ierr);
	ierr = DMDAGetElementsQ1(daH,&nel,&nen,&elnidx);CHKERRQ(ierr);

    for (s=0; s<nsteps; s++) {
    
        /* get acces to the vector V */
        ierr = DMGlobalToLocalBegin(daH,H,INSERT_VALUES,local_H);CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(  daH,H,INSERT_VALUES,local_H);CHKERRQ(ierr);
        ierr = VecGetArray(local_H,&LA_H);CHKERRQ(ierr);

        ierr = VecZeroEntries(rhs);CHKERRQ(ierr);
        ierr = VecZeroEntries(local_rhs);CHKERRQ(ierr);
        ierr = VecGetArray(local_rhs,&LA_rhs);CHKERRQ(ierr);

        ierr = VecZeroEntries(diagM);CHKERRQ(ierr);
        ierr = VecZeroEntries(local_M);CHKERRQ(ierr);
        ierr = VecGetArray(local_M,&LA_M);CHKERRQ(ierr);
        
        /* compute rhs */
        for (e=0; e<nel; e++) {
            PetscInt idx2d[] = { 2, 3, 6, 7 };
            PetscInt ii;
            PetscReal xi[] = {0.0,0.0};
            PetscReal Ni[4],GNix[4],GNiz[4],dNdx[4],dNdz[4];
            PetscReal dhdx,dhdz,el_rhs[4],el_M[4],detJ,el_h2d[4],el_coords2d[2*4];
            PetscInt gI,gJ,gK;
            PetscInt ge_eqnums2d[4];
            PetscReal el_rhs3d[8],el_M3d[8];
            
            /* extract equation numbers for element */
            ierr = DMDAEQ1_GetElementLocalIndicesDOF(ge_eqnums,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
            
            /* extract x,z coords for element */
            ierr = DMDAEQ1_GetVectorElementField_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);

            /* extract height for element */
            ierr = DMDAEQ1_GetScalarElementField_3D(el_h,(PetscInt*)&elnidx[nen*e],LA_H);CHKERRQ(ierr);

            /* extract face values */
            for (i=0; i<4; i++) {
                PetscInt idx3;
                
                idx3 = idx2d[i];
                
                el_h2d[i] = el_h[idx3];

                el_coords2d[2*i  ] = el_coords[3*idx3    ]; // x coordinate
                el_coords2d[2*i+1] = el_coords[3*idx3 + 2]; // z coordinate
                                 
                ge_eqnums2d[i] = ge_eqnums[idx3]; // face indices
            }
            
            /* determine if element is located on surface */
            if (only_update_surface) {
                PetscBool surface_face = PETSC_FALSE;
                
                for (i=0; i<4; i++) {
                    ierr = DMDAConvertLocalGhostNodeIndex2GlobalnInJnK(daH,ge_eqnums2d[i],&gI,&gJ,&gK);CHKERRQ(ierr);
                    if (gJ == NJ-1) {
                        surface_face = PETSC_TRUE;
                        break;
                    }
                }
                if (!surface_face) {
                    continue;
                }
            }
            
            
            /* compute grad(h) at xi = (0 , 0) */
            /* compute grad(N) at xi = (0 , 0) */
            /* Ni[0] = 0.25 * (1-x)*(1-z) ; GNix[] = -0.25*(1-z) ; GNiz[] = -0.25*(1-x) */
            /* Ni[1] = 0.25 * (1+x)*(1-z) ; GNix[] =  0.25*(1-z) ; GNiz[] = -0.25*(1+x) */
            /* Ni[2] = 0.25 * (1-x)*(1+z) ; GNix[] =  0.25*(1+z) ; GNiz[] =  0.25*(1-x) */
            /* Ni[3] = 0.25 * (1+x)*(1+z) ; GNix[] = -0.25*(1+z) ; GNiz[] =  0.25*(1+x) */
            
            P3D_ConstructNi_Q1_2D(xi,Ni);
            P3D_ConstructGNi_Q1_2D(xi,GNix,GNiz);
            P3D_evaluate_geometry_elementQ1_2D(el_coords2d,GNix,GNiz,&detJ,dNdx,dNdz);
            
            dhdx = dNdx[0]*el_h2d[0] + dNdx[1]*el_h2d[1] + dNdx[2]*el_h2d[2] + dNdx[3]*el_h2d[3];
            dhdz = dNdz[0]*el_h2d[0] + dNdz[1]*el_h2d[1] + dNdz[2]*el_h2d[2] + dNdz[3]*el_h2d[3];
            
            // 1 pnt quadrature rule, so w_q = 4.0 and there is no summation to perform
            //PetscMemzero(el_rhs,sizeof(PetscReal)*4);
            //PetscMemzero(el_M,sizeof(PetscReal)*4);
            for (ii=0; ii<4; ii++) {
                /* Kij = GNix[i].GNix[j] + GNiz[i].GNiz[j] */
                el_rhs[ii] = -4.0 * diffusivity * (dNdx[ii] * dhdx + dNdz[ii] * dhdz) * detJ;
                
                /* compute diag(M) */
                el_M[ii]   = 4.0 * (Ni[ii] * Ni[ii]) * detJ;
            }
            //printf("el_rhs %1.4e %1.4e %1.4e %1.4e \n",el_rhs[0],el_rhs[1],el_rhs[2],el_rhs[3]);
            //printf("el_M %1.4e %1.4e %1.4e %1.4e \n",el_M[0],el_M[1],el_M[2],el_M[3]);
            
            PetscMemzero(el_rhs3d,sizeof(PetscReal)*8);
            PetscMemzero(el_M3d,sizeof(PetscReal)*8);
            /* I drop in 1.0 here as later, post assembly we will invert this vector */
            /*
            for (i=0; i<8; i++) {
                el_M3d[i] = 1.0;
            }
            */
            /* map 2d diffusion into 3d dmda vector */
            for (i=0; i<4; i++) {
                el_rhs3d[ idx2d[i] ] = el_rhs[i];
                el_M3d[   idx2d[i] ] = el_M[i];
            }
            
            ierr = DMDAEQ1_SetValuesLocalStencil_AddValues_DOF(LA_rhs,1,ge_eqnums,el_rhs3d);CHKERRQ(ierr);
            ierr = DMDAEQ1_SetValuesLocalStencil_AddValues_DOF(LA_M,1,ge_eqnums,el_M3d);CHKERRQ(ierr);
        }

        ierr = VecRestoreArray(local_H,&LA_H);CHKERRQ(ierr);
        ierr = VecRestoreArray(local_M,&LA_M);CHKERRQ(ierr);
        ierr = VecRestoreArray(local_rhs,&LA_rhs);CHKERRQ(ierr);

        /* assemble rhs, diag(M) */
        ierr = DMLocalToGlobalBegin(daH,local_rhs,ADD_VALUES,rhs);CHKERRQ(ierr);
        ierr = DMLocalToGlobalEnd(  daH,local_rhs,ADD_VALUES,rhs);CHKERRQ(ierr);

        ierr = DMLocalToGlobalBegin(daH,local_M,ADD_VALUES,diagM);CHKERRQ(ierr);
        ierr = DMLocalToGlobalEnd(  daH,local_M,ADD_VALUES,diagM);CHKERRQ(ierr);

        /* update */
        //ierr = VecPointwiseDivide(rhs,rhs,diagM);CHKERRQ(ierr);
        ierr = VecGetArray(diagM,&LA_M);CHKERRQ(ierr);
        for (i=0; i<nM; i++) {
            if (fabs(PetscRealPart(LA_M[i])) > 1.0e-20) {
                LA_M[i] = 1.0/LA_M[i];
            }
        }
        ierr = VecRestoreArray(diagM,&LA_M);CHKERRQ(ierr);
        
        ierr = VecPointwiseMult(rhs,rhs,diagM);CHKERRQ(ierr);
        
        ierr = VecAXPY(H,dt_explicit,rhs);CHKERRQ(ierr);
        
        /* force boundary conditions */
        ierr = DMDAVecGetArray(daH,H,&LA_H3);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(daH,Hinit,&LA_Hinit3);CHKERRQ(ierr);

        for (k=sk; k<sk+nk; k++) {
            for (j=sj; j<sj+nj; j++) {
                for (i=si; i<si+ni; i++) {
                    
                    if (dirichlet_west) {
                        if (i == 0) { LA_H3[k][j][i]    = LA_Hinit3[k][j][i]; }
                    }
                    if (dirichlet_east) {
                        if (i == NI-1) { LA_H3[k][j][i] = LA_Hinit3[k][j][i]; }
                    }
                    if (dirichlet_front) {
                        if (k == 0) { LA_H3[k][j][i]    = LA_Hinit3[k][j][i]; }
                    }
                    if (dirichlet_back) {
                        if (k == NK-1) { LA_H3[k][j][i] = LA_Hinit3[k][j][i]; }
                    }
                }
            }
        }
        
        ierr = DMDAVecRestoreArray(daH,Hinit,&LA_Hinit3);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(daH,H,&LA_H3);CHKERRQ(ierr);

        // output
        if (view_surface_diffusion){
            char name[1000];
            Vec coordsH;
            DM cdaH;
            int s32;
            
            ierr = DMGetCoordinateDM(daH,&cdaH);CHKERRQ(ierr);
            ierr = DMGetCoordinates(daH,&coordsH);CHKERRQ(ierr);
            ierr = VecStrideScatter(H,1,coordsH,INSERT_VALUES);CHKERRQ(ierr);

            PetscMPIIntCast(s,&s32);
            sprintf(name,"surface_diffusion_%.4d.vtk",s32);
            ierr = DMDAViewPetscVTK(daH,H,name);CHKERRQ(ierr);
        }
        //
        
    }
    ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
    
    ierr = DMRestoreLocalVector(daH,&local_rhs);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(daH,&local_M);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(daH,&local_H);CHKERRQ(ierr);
    
    
    ierr = DMGetLocalVector(daH,&local_H);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(daH,H,INSERT_VALUES,local_H);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(  daH,H,INSERT_VALUES,local_H);CHKERRQ(ierr);
    
    {
        DMDACoor3d ***LA_coords3;
        PetscInt siv,niv,sjv,njv,skv,nkv;
        PetscInt sig,nig,sjg,njg,skg,nkg;
        PetscInt NIv,NJv,NKv;
        
        ierr = DMDAGetInfo(dav,0,&NIv,&NJv,&NKv,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
        
        ierr = DMDAGetGhostCorners(daH,&sig,&sjg,&skg,&nig,&njg,&nkg);CHKERRQ(ierr);
        
        ierr = DMDAVecGetArray(daH,local_H,&LA_H3);CHKERRQ(ierr);
        
        ierr = DMDAGetCorners(dav,&siv,&sjv,&skv,&niv,&njv,&nkv);CHKERRQ(ierr);
        ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
        ierr = DMGetCoordinates(dav,&coords);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(cda,coords,&LA_coords3);CHKERRQ(ierr);

        if (only_update_surface) {
            /* copy only heights from JMAX surface onto Q2 mesh */
            if ( NJv == sjv+njv ) {
                j = NJv - 1;
                for (k=skv; k<skv+nkv; k++) {
                    for (i=siv; i<siv+niv; i++) {
                        if (i%2 != 0) { continue; }
                        if (j%2 != 0) { continue; }
                        if (k%2 != 0) { continue; }
                        
                        if (i/2 < sig) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"i min bound failed"); }
                        if (i/2 >= sig+nig) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"i max bound failed"); }

                        if (j/2 < sjg) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"j min bound failed"); }
                        if (j/2 >= sjg+njg) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"j max bound failed"); }

                        if (k/2 < skg) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"k min bound failed"); }
                        if (k/2 >= skg+nkg) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"k max bound failed"); }

                        LA_coords3[k][j][i].y = LA_H3[k/2][j/2][i/2];
                    }
                }
            }
        } else {
            /* copy all heights onto Q2 mesh */
            for (k=skv; k<skv+nkv; k++) {
                for (j=sjv; j<sjv+njv; j++) {
                    for (i=siv; i<siv+niv; i++) {
                        if (i%2 != 0) { continue; }
                        if (j%2 != 0) { continue; }
                        if (k%2 != 0) { continue; }
                        
                        if (i/2 < sig) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"i min bound failed"); }
                        if (i/2 >= sig+nig) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"i max bound failed"); }
                        
                        if (j/2 < sjg) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"j min bound failed"); }
                        if (j/2 >= sjg+njg) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"j max bound failed"); }
                        
                        if (k/2 < skg) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"k min bound failed"); }
                        if (k/2 >= skg+nkg) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"k max bound failed"); }
                        
                        LA_coords3[k][j][i].y = LA_H3[k/2][j/2][i/2];
                    }
                }
            }
        }
        
        ierr = DMDAVecRestoreArray(cda,coords,&LA_coords3);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(daH,local_H,&LA_H3);CHKERRQ(ierr);

        ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
    }
    
    ierr = DMRestoreLocalVector(daH,&local_H);CHKERRQ(ierr);
    
    
    /* billinearize Q2 mesh */
    ierr = DMDABilinearizeQ2Elements(dav);CHKERRQ(ierr);
    
    ierr = VecDestroy(&rhs);CHKERRQ(ierr);
    ierr = VecDestroy(&diagM);CHKERRQ(ierr);
    ierr = VecDestroy(&H);CHKERRQ(ierr);
    ierr = VecDestroy(&Hinit);CHKERRQ(ierr);
    
    ierr = DMDestroy(&daH);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

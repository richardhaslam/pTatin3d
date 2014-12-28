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
 **    Filename:      spm_utils.c
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

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_redundant.h"
#include "ptatin3d_defs.h"
#include "element_type_Q2.h"
#include "element_utils_q2.h"
#include "MPntStd_def.h"
#include "material_point_point_location.h"
#include "dmda_update_coords.h"
#include "dmda_remesh.h"
#include "mesh_update.h"
#include "dmda_update_coords.h"
#include "dmda_view_petscvtk.h"
#include "spm_utils.h"

/* 
Contains routines for assisting with:
 
 * Parallel(mechanical model) + sequential (spm)
     DMDAGatherIKRedundantSurfaceDMDA:  Gather a redundant surface mesh defined in the I-K plane on rank 0 (MSurf0) from a parallel mechanical model (MSurf)
     InterpolateMSurf0ToSPMSurfIKGrid:  Perform interpolation between the MSurf0 and a user defined surface to be used in the SPM (SPMSurf)
     InterpolateSPMSurfIKGridToMSurf0:  Perform interpolation between SPMSurf and MSurf0
     DMDAScatterIKRedundantSurfaceDMDA: scatter between the MSurf0 and the mechanical model MSurf

 * Parallel (mechanical model) + parallel, semi-redundant (spm) using only cores defined at the surface of the domain
     DMDAGatherIKSemiRedundantSurfaceDMDA:  Gather a semi redundant parallel surface mesh defined in the I-K plane
                                            using on all surface ranks (MSurfSR) from a parallel mechanical model (MSurf)
     DMDAScatterIKSemiRedundantSurfaceDMDA: Scatter between the MSurfSR and the mechanical model MSurf

 * Parallel (mechanical model) + parallel (spm) using all cores 
     DMDAGatherIKSurfaceDMDA:  Gather a parallel surface mesh defined in the I-K plane
                                            using on all ranks (MSurf) from a parallel mechanical model (MSurf)
     DMDAScatterIKSurfaceDMDA: Scatter between the MSurf and the mechanical model MSurf
     Note that these functions create an spm DM which overlaps the global node layout of the mechanical DM
*/

/* Parallel mechanical -> sequential spm */
#undef __FUNCT__
#define __FUNCT__ "DMDAGatherIKRedundantSurfaceDMDA"
PetscErrorCode DMDAGatherIKRedundantSurfaceDMDA(DM dm_mech,DM *dm_msurf0)
{
	PetscErrorCode ierr;
	DM             dm_surf;
	PetscInt       si,sj,sk,ni,nj,nk,M,N,P,startI,endI,startK,endK,endJ;
	MPI_Comm       comm;
	PetscMPIInt    rank;
	
	
	PetscFunctionBegin;
	PetscObjectGetComm((PetscObject)dm_mech,&comm);
	ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	
	/* collect top and bottom surface da on rank0 */
	ierr = DMDAGetInfo(dm_mech,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm_mech,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
	
	/* collective call
	 rank0 fetchs N, all others fetch something local and tiny 
	*/
	if (rank != 0) {
		startI = si;
		endI   = si+1;
		startK = sk;
		endK   = sk+1;

		endJ   = sj+nj;
	} else {
		startI = 0;
		endI   = M;
		startK = 0;
		endK   = P;

		endJ   = N;
	}
	ierr = DMDACreate3dRedundant(dm_mech,startI,endI,endJ-1,endJ,startK,endK,1,&dm_surf);CHKERRQ(ierr);
	
	*dm_msurf0 = NULL;
	if (rank == 0) {
		*dm_msurf0 = dm_surf;
	}
	
	if (rank != 0) {
		ierr = DMDestroy(&dm_surf);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InterpolateMSurf0ToSPMSurfIKGrid"
PetscErrorCode InterpolateMSurf0ToSPMSurfIKGrid(DM dm_msurf0,PetscInt spm_mi,PetscInt spm_mj,PetscReal *spm_coords,PetscReal *spm_H)
{
	PetscErrorCode ierr;
	PetscInt       ei,ek,ni,nj,nk,msurf_mi,msurf_mk;
	PetscInt       i,k,n,spm_nodes;
	PetscReal      H_p,el_H[Q2_NODES_PER_EL_2D],Ni[Q2_NODES_PER_EL_2D];
	Vec            coords;
	PetscScalar    *msurf_coords;
	PetscReal      *_msurf_coords;
	PetscInt       *msurf_element,elcnt,*el,nid[Q2_NODES_PER_EL_2D];
	
	
	PetscFunctionBegin;
	/* generate an element->node mapping for msurf0 for Q2 */
	ierr = DMDAGetInfo(dm_msurf0,0,&ni,&nj,&nk,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*ni*nk*Q2_NODES_PER_EL_2D,&msurf_element);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscReal)*2*ni*nk,&_msurf_coords);CHKERRQ(ierr);
	
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
	
	ierr = DMGetCoordinates(dm_msurf0,&coords);CHKERRQ(ierr);
	ierr = VecGetArray(coords,&msurf_coords);CHKERRQ(ierr);
	
	for (k=0; k<ni*nk; k++) {
		_msurf_coords[2*k+0] = msurf_coords[3*k+0]; // x
		_msurf_coords[2*k+1] = msurf_coords[3*k+2]; // z
	}
	
	/* a) Insert coords of spmsurf into material point structure for interpolation purposes */
	/* b) Map spmsurf.x -> mp.x; spmsurf.z -> mp.y */
	/* c) perform point location */
	/* d) perform interpolation from msurf0 -> spm_H */
	spm_nodes = (spm_mi+1)*(spm_mj+1);
	for (n=0; n<spm_nodes; n++) {
		MPntStd mp;

		/* a),b) */
		mp.coor[0] = spm_coords[2*n+0];
		mp.coor[1] = spm_coords[2*n+1];
		mp.coor[2] = 0.0;
		/* c) */
		InverseMappingDomain_2dQ2(1.0e-10,10,PETSC_FALSE,PETSC_FALSE,
															_msurf_coords,msurf_mi,msurf_mk,msurf_element,
															1,&mp);
		if (mp.wil == -1) {
			PetscPrintf(PETSC_COMM_SELF,"Error: point %d not found\n",n);
			PetscPrintf(PETSC_COMM_SELF,"  pos(x,z) %1.4e %1.4e\n",mp.coor[0],mp.coor[1]);
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"  location of spm node in mechanical surface could not be determined");
		}
		
		/* d) */
		/* get element heights */
		for (i=0; i<Q2_NODES_PER_EL_2D; i++) {
			PetscInt nidx;
			
			nidx = msurf_element[ Q2_NODES_PER_EL_2D*mp.wil + i ];
			el_H[i] = msurf_coords[3*nidx+1]; /* fetch the J component from the coordinate array */
		}
		
		/* evaluate basis and interplate */
		P3D_ConstructNi_Q2_2D(mp.xi,Ni);
		H_p = 0.0;
		for (i=0; i<Q2_NODES_PER_EL_2D; i++) {
			H_p += Ni[i] * el_H[i];
		}
		
		/* assign interplated height */
		spm_H[n] = H_p;
	}

	ierr = VecRestoreArray(coords,&msurf_coords);CHKERRQ(ierr);
	ierr = PetscFree(msurf_element);CHKERRQ(ierr);
	ierr = PetscFree(_msurf_coords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
	Treat nodal spm data as "cell-centered" with a mesh which extends from x0-0.5*dx 
*/
#undef __FUNCT__
#define __FUNCT__ "InterpolateSPMSurfIKGridToMSurf0"
PetscErrorCode InterpolateSPMSurfIKGridToMSurf0(PetscInt spm_mi,PetscInt spm_mj,PetscReal *spm_coords,PetscReal *spm_H,DM dm_msurf0)
{
	PetscErrorCode ierr;
	PetscInt       ii,jj,i,ni,nj,nk,msurf_nodes,D;
	PetscInt       spm_i,spm_j,spm_ni,spm_nj,spm_nodes,spm_nidx;
	Vec            coords;
	PetscScalar    *msurf_coords;
	PetscReal      spmsurf_x0,spmsurf_y0,spmsurf_x1,spmsurf_y1,spmsurf_dx,spmsurf_dy;
	
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(dm_msurf0,0,&ni,&nj,&nk,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	msurf_nodes = ni*nj*nk;
	ierr = DMGetCoordinates(dm_msurf0,&coords);CHKERRQ(ierr);
	ierr = VecGetArray(coords,&msurf_coords);CHKERRQ(ierr);

	spm_ni = spm_mi + 1;
	spm_nj = spm_mj + 1;
	spm_nodes = spm_ni * spm_nj;
	D = 8;
	
	spmsurf_x0 = spm_coords[2*0+0];
	spmsurf_y0 = spm_coords[2*0+1];
	spmsurf_x1 = spm_coords[2*(spm_nodes-1)+0];
	spmsurf_y1 = spm_coords[2*(spm_nodes-1)+1];
	spmsurf_dx = (spmsurf_x1 -spmsurf_x0)/(PetscReal)spm_mi;
	spmsurf_dy = (spmsurf_y1 -spmsurf_y0)/(PetscReal)spm_mj;

	spmsurf_x0 = spmsurf_x0 - 0.5 * spmsurf_dx;
	spmsurf_y0 = spmsurf_y0 - 0.5 * spmsurf_dy;

	for (i=0; i<msurf_nodes; i++) {
		PetscReal msurf_x,msurf_z,spm_surface_Hi;
		PetscInt nsamples;
		
		msurf_x = msurf_coords[3*i+0];
		msurf_z = msurf_coords[3*i+2];
		
		/* locate surf_x,z within the spm "cell-centered" representation */
		spm_i = ( msurf_x - spmsurf_x0 ) / spmsurf_dx;
		spm_j = ( msurf_z - spmsurf_y0 ) / spmsurf_dy;
		
		if (spm_i == spm_ni) { spm_i--; }
		if (spm_j == spm_nj) { spm_j--; }
		
		/* scan over neighbourhood and average */
		spm_surface_Hi = 0.0;
		nsamples = 0;
		for (jj=spm_j-D; jj<spm_j+D; jj++) {
			for (ii=spm_i-D; ii<spm_i+D; ii++) {
				
				if (ii < 0) { continue; }
				if (jj < 0) { continue; }
				if (ii >= spm_ni) { continue; }
				if (jj >= spm_nj) { continue; }
				
				/* node index in spm data */
				spm_nidx = ii + jj * spm_ni;
				if (spm_nidx >= spm_nodes) {
					SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"spm node index %D larger than max %D",spm_nidx,spm_nodes);
				}
				if (spm_nidx < 0) {
					SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"spm node index %D < 0",spm_nidx);
				}
				
				spm_surface_Hi += spm_H[spm_nidx];
				nsamples++;
			}
		}
		spm_surface_Hi = spm_surface_Hi / (PetscReal)nsamples;
		
		/* set new coordinate value into J */
		msurf_coords[3*i+1] = spm_surface_Hi;
	}
	
	ierr = VecRestoreArray(coords,&msurf_coords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
 Copy the parts we need for each sub-domain on the surface from spm_H.
 MPI_Send
 Insert on each sub-domain 
*/
#undef __FUNCT__
#define __FUNCT__ "DMDAScatterIKRedundantSurfaceDMDA"
PetscErrorCode DMDAScatterIKRedundantSurfaceDMDA(DM dm_msurf0,DM dm_mech)
{
	PetscErrorCode ierr;
	PetscInt       pID,pI,pJ,pK,rI,rJ,rK,rIJ,rank;
	PetscInt       i,k,ii,kk,si,ei,sk,ek,sj,ej,M,N,P,pi,pk;
	PetscMPIInt    _rank;
	PetscBool      surface_rank;
	PetscReal      *surface_chunk;
	PetscInt       nchunk;
	DM             cda_mech;
	const PetscInt *lni,*lnj,*lnk;
	PetscInt       *lsi,*lsk;
	Vec            coords;
	DMDACoor3d     ***LA_coords;
	MPI_Comm       comm;
	MPI_Status     status;
	PetscInt       nl;
	
	
	PetscFunctionBegin;
	/* --------------------------------------------------- */
	/* determine ranks living on the free surface j=JMAX-1 */
	PetscObjectGetComm((PetscObject)dm_mech,&comm);
	ierr = MPI_Comm_rank(comm,&_rank);CHKERRQ(ierr);
	rank = (PetscInt)_rank;
	ierr = DMDAGetInfo(dm_mech,0,&M,&N,&P,&pI,&pJ,&pK,0,0,0,0,0,0);CHKERRQ(ierr);

	/* convert rank to rI,rJ,rK */
	rK  = rank / (pI*pJ);
	rIJ = rank - rK * pI*pJ;  
	rJ = rIJ / pI;
	rI = rIJ - rJ*pI;
	
	surface_rank = PETSC_FALSE;
	if (rJ == (pJ-1)) {
		surface_rank = PETSC_TRUE;
	}
	//nsurface_ranks = pI * pK;
	
	
	/* --------------------------------------------- */
	/* need corner info for all ranks on rank 0      */
	ierr = DMDAGetOwnershipRanges(dm_mech,&lni,&lnj,&lnk);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*pI,&lsi);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*pK,&lsk);CHKERRQ(ierr);
	lsi[0] = 0;
	for (pi=1; pi<pI; pi++) {
		lsi[pi] = lsi[pi-1] + lni[pi-1]; 
	}
	lsk[0] = 0;
	for (pk=1; pk<pK; pk++) {
		lsk[pk] = lsk[pk-1] + lnk[pk-1]; 
	}
		
	/* allocate space for receiving elevations */
	ierr = DMDAGetCorners(dm_mech,&si,&sj,&sk,&ei,&ej,&ek);CHKERRQ(ierr);
	nchunk = ei * ek;
	ierr = PetscMalloc(sizeof(PetscReal)*nchunk,&surface_chunk);CHKERRQ(ierr);
			
	//if (surface_rank) {
	//	printf("rank[%d] surface rank %d [%d-%d]\n",rank,surface_rank,sj,sj+ej);
	//}
	
	if (rank == 0) {
		DM cda_spmsurf;
		
		
		ierr = DMGetCoordinateDM(dm_msurf0,&cda_spmsurf);CHKERRQ(ierr);
		ierr = DMGetCoordinates(dm_msurf0,&coords);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda_spmsurf,coords,&LA_coords);CHKERRQ(ierr);
		
		for (pk=0; pk<pK; pk++) {
			for (pi=0; pi<pI; pi++) {
				PetscReal *sendchunk;
				
				pID = pi + (pJ-1)*pI + pk*pI*pJ;

				//printf("I: sending chunks lsi=%d to lei=%d to rank %d\n",lsi[pi],lsi[pi]+lni[pi],pID);
				//printf("K: sending chunks lsk=%d to lek=%d to rank %d\n",lsk[pk],lsk[pk]+lnk[pk],pID);
				
				/* alloc space */
				nl = lni[pi] * lnk[pk];
				//printf("...sending a total of %d values to rank %d\n",nl,pID);
				ierr = PetscMalloc(sizeof(PetscReal)*nl,&sendchunk);CHKERRQ(ierr);
				
				/* copy chunks - order of insertition is important when extracting result */
				for (k=lsk[pk]; k<lsk[pk]+lnk[pk]; k++) {
					for (i=lsi[pi]; i<lsi[pi]+lni[pi]; i++) {
						ii = i - lsi[pi];
						kk = k - lsk[pk];
						
						sendchunk[ii + kk*lni[pi]*1] = LA_coords[k][0][i].y;
					}
				}
				
				/* send */
				ierr = MPI_Send(sendchunk,nl,MPIU_REAL,pID,0,comm);CHKERRQ(ierr);
				
				ierr = PetscFree(sendchunk);CHKERRQ(ierr);
			}
		}

		ierr = DMDAVecRestoreArray(cda_spmsurf,coords,&LA_coords);CHKERRQ(ierr);
	}

	/* collect results */
	ierr = DMGetCoordinateDM(dm_mech,&cda_mech);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dm_mech,&coords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_mech,coords,&LA_coords);CHKERRQ(ierr);
	if (surface_rank) {
		if (rJ != (pJ-1)) {
			SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"Rank of surface sub-domain must define rJ = %D - detected %D",pJ-1,rJ);
		}
		
		/* receive */
		pID = rI + rJ*pI + rK*pI*pJ;
		nl = lni[rI] * lnk[rK];
		//printf("[%d]...reciving a total of %d values to rank [allocated space %d]\n",rank,nl,nchunk);
		ierr = MPI_Recv(surface_chunk,nl,MPIU_REAL,0,MPI_ANY_TAG,comm,&status);CHKERRQ(ierr);
		
		/* insert */
		for (k=sk; k<sk+ek; k++) {
			for (i=si; i<si+ei; i++) {
				PetscInt local_idx;
				
				local_idx = (i-si) + (k-sk)*ei;
				LA_coords[k][N-1][i].y = surface_chunk[ local_idx ];
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda_mech,coords,&LA_coords);CHKERRQ(ierr);
	
	ierr = PetscFree(surface_chunk);CHKERRQ(ierr);
	ierr = PetscFree(lsi);CHKERRQ(ierr);
	ierr = PetscFree(lsk);CHKERRQ(ierr);
	
	/* update */
	ierr = DMDAUpdateGhostedCoordinates(dm_mech);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test_spm_utils_MPItoSEQ"
PetscErrorCode test_spm_utils_MPItoSEQ(DM dav)
{
	DM             dm_spmsurf0;
	PetscInt       JMAX;
	PetscErrorCode ierr;
	
	
	ierr = DMDAGatherIKRedundantSurfaceDMDA(dav,&dm_spmsurf0);CHKERRQ(ierr);
	if (dm_spmsurf0) {
		ierr = DMDAViewPetscVTS(dm_spmsurf0,NULL,"surf_extraction_ic.vtk");CHKERRQ(ierr);
	}
	if (dm_spmsurf0) {
		PetscInt smx,smy,ii,jj;
		PetscReal *sc,*sH;
		
		smx = 128;
		smy = 128;
		
		PetscMalloc(sizeof(PetscReal)*2*(smx+1)*(smy+1),&sc);
		PetscMalloc(sizeof(PetscReal)*(smx+1)*(smy+1),&sH);
		
		for (jj=0; jj<smy+1; jj++) {
			for (ii=0; ii<smx+1; ii++) {
				sc[2*(ii+(smx+1)*jj)+0] = 0.0 + ii * 1.0/((PetscReal)smx);
				sc[2*(ii+(smx+1)*jj)+1] = 0.0 + jj * 1.0/((PetscReal)smy);
				sH[ii+(smx+1)*jj] = 0.0;
			}
		}
		
		ierr = InterpolateMSurf0ToSPMSurfIKGrid(dm_spmsurf0,smx,smy,sc,sH);CHKERRQ(ierr);
		// dump
		{
			FILE *fp = fopen("smesh.gp","w");
			for (jj=0; jj<smy+1; jj++) {
				for (ii=0; ii<smx+1; ii++) {
					fprintf(fp,"%1.4e %1.4e %1.4e\n",sc[2*(ii+(smx+1)*jj)+0],sc[2*(ii+(smx+1)*jj)+1],sH[ii+(smx+1)*jj]);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);
		}
		// perturb and dump
		{
			FILE *fp = fopen("smesh2.gp","w");
			for (jj=0; jj<smy+1; jj++) {
				for (ii=0; ii<smx+1; ii++) {
					double x,y;
					
					x = sc[2*(ii+(smx+1)*jj)+0];
					y = sc[2*(ii+(smx+1)*jj)+1];
					
					sH[ii+(smx+1)*jj] += 0.1 * sin(3.2*x*M_PI)*cos(6.6*y*M_PI);
					
					fprintf(fp,"%1.4e %1.4e %1.4e\n",sc[2*(ii+(smx+1)*jj)+0],sc[2*(ii+(smx+1)*jj)+1],sH[ii+(smx+1)*jj]);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);
		}
		
		ierr = InterpolateSPMSurfIKGridToMSurf0(smx,smy,sc,sH,dm_spmsurf0);CHKERRQ(ierr);
		ierr = DMDAViewPetscVTS(dm_spmsurf0,NULL,"surf_extraction_interp.vtk");CHKERRQ(ierr);
		
		PetscFree(sc);
		PetscFree(sH);
	}
	
	ierr = DMDAScatterIKRedundantSurfaceDMDA(dm_spmsurf0,dav);CHKERRQ(ierr);
	
	if (dm_spmsurf0) {
		ierr = DMDestroy(&dm_spmsurf0);CHKERRQ(ierr);
	}
	
	ierr = DMDAGetInfo(dav,0,0,&JMAX,0,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,JMAX);CHKERRQ(ierr);
	ierr = DMDABilinearizeQ2Elements(dav);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* Parallel mechanical -> semi-redundant spm */

/*
  Create communicator for the surface ranks
  Extract partitioning from dm_mechanical, using this to create semi-redundant 2D DMDA defined in x,y
  Fetch part of coordinates need from dm_mechanical
*/
#undef __FUNCT__
#define __FUNCT__ "DMDAGatherIKSemiRedundantSurfaceDMDA"
PetscErrorCode DMDAGatherIKSemiRedundantSurfaceDMDA(DM dm_mech,DM *_dm_msurf)
{
	PetscFunctionBegin;
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"SemiRedundantSurfaceDMDA is not supported - only parallel or sequential are permitted");
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAScatterIKSemiRedundantSurfaceDMDA"
PetscErrorCode DMDAScatterIKSemiRedundantSurfaceDMDA(DM dm_msurf,DM dm_mech)
{
	PetscFunctionBegin;
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"SemiRedundantSurfaceDMDA is not supported - only parallel or sequential are permitted");
	PetscFunctionReturn(0);
}

/* Parallel mechanical -> parallel spm (with overlapping global nodes) */

#undef __FUNCT__
#define __FUNCT__ "DMDAGatherIKSurfaceDMDA"
PetscErrorCode DMDAGatherIKSurfaceDMDA(DM dm_mech,DM *_dm_msurf,Vec *_elevation)
{
	PetscErrorCode ierr;
	PetscInt       si,sj,sk,nx,ny,nz,M,N,P;
	PetscInt       si2d,sj2d,nx2d,ny2d,i,k;
	MPI_Comm       comm;
	DM             dm_red_spm,dm_surf;
	DM             dm_red_spm_coords,dm_surf_coords;
	Vec            coords_surf,coords_red_spm;
	DMDACoor3d     ***LA_coords_red_spm;
	DMDACoor2d     **LA_coords_surf;
	PetscScalar    **LA_height;
	PetscMPIInt    rank;
	Vec            height;

	
	PetscFunctionBegin;	
	PetscObjectGetComm((PetscObject)dm_mech,&comm);
	ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dm_mech,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDACreate2d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,M,P, PETSC_DECIDE,PETSC_DECIDE, 1,1, 0,0,&dm_surf);CHKERRQ(ierr);	
	ierr = DMDASetUniformCoordinates(dm_surf, 0.0,1.0, 0.0,1.0, 0.,0.);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm_surf,&si2d,&sj2d,NULL,&nx2d,&ny2d,NULL);CHKERRQ(ierr);
	
	/* fetch coordinates i need from the mechanical domain */
	ierr = DMDAGetCorners(dm_mech,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	
	//printf("rank %d: vol  si,ei=%.4d-%.4d: sj,ej=%.4d-%.4d: sk,ek=%.4d-%.4d \n", rank, si,si+nx,sj,sj+ny,sk,sk+nz);
	//printf("rank %d: surf si,ei=%.4d-%.4d:                  sj,ej=%.4d-%.4d \n", rank, si2d,si2d+nx2d,sj2d,sj2d+ny2d);
	
	ierr = DMDACreate3dRedundant(dm_mech,si2d,si2d+nx2d,N-1,N,sj2d,sj2d+ny2d, 1, &dm_red_spm);CHKERRQ(ierr);
	//sprintf(name,"surf%d.vtk",rank);
	//ierr = DMDAViewPetscVTS(da_red_spm,0,name);CHKERRQ(ierr);
	
	/* copy these values into my parallel surface mesh x,y,z (vol) => x,y (surf) */
	ierr = DMGetCoordinates(dm_surf,&coords_surf);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dm_red_spm,&coords_red_spm);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(dm_surf,&dm_surf_coords);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dm_red_spm,&dm_red_spm_coords);CHKERRQ(ierr);
	
	
	ierr = DMDAVecGetArray(dm_surf_coords,coords_surf,&LA_coords_surf);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dm_red_spm_coords,coords_red_spm,&LA_coords_red_spm);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(dm_surf,&si,&sk,NULL,&nx,&nz,NULL);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (i=si; i<si+nx; i++) {
			LA_coords_surf[k][i].x = LA_coords_red_spm[k-sk][0][i-si].x;
			LA_coords_surf[k][i].y = LA_coords_red_spm[k-sk][0][i-si].z;
		}
	}
	
	ierr = DMCreateGlobalVector(dm_surf,&height);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dm_surf,height,&LA_height);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (i=si; i<si+nx; i++) {
			LA_height[k][i] = LA_coords_red_spm[k-sk][0][i-si].y;
		}
	}
	ierr = DMDAVecRestoreArray(dm_surf,height,&LA_height);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(dm_red_spm_coords,coords_red_spm,&LA_coords_red_spm);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(dm_surf_coords,coords_surf,&LA_coords_surf);CHKERRQ(ierr);
	
	ierr = DMDestroy(&dm_red_spm);CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(dm_surf);CHKERRQ(ierr);
	
	*_elevation = height;
	*_dm_msurf  = dm_surf;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAScatterIKSurfaceDMDA"
PetscErrorCode DMDAScatterIKSurfaceDMDA(DM dm_msurf,Vec height,DM dm_mech)
{
	PetscErrorCode ierr;
	MPI_Comm       comm;
	PetscMPIInt    _rank;
	PetscInt       rank;
	PetscInt       M,N,P,pI,pJ,pK,rJ,rK,rIJ,pi,pk;
	PetscBool      surface_rank;
	const PetscInt *lni,*lnj,*lnk;
	PetscInt       *lsi,*lsk;
	PetscInt       si,sj,sk,nx,ny,nz,i,j,k,c,nsurf_nodes,*surf_indices,M2d,P2d;
	Vec            height_self,height_natural;
	PetscScalar    *LA_height_self;
	IS             is_surf,is_local;
	VecScatter     sctx;
	Vec            coords_mech;
	DM             cda_mech;
	DMDACoor3d     ***LA_coords_mech;
	
	
	PetscFunctionBegin;
	/* --------------------------------------------------- */
	/* determine ranks living on the free surface j=JMAX-1 */
	PetscObjectGetComm((PetscObject)dm_mech,&comm);
	ierr = MPI_Comm_rank(comm,&_rank);CHKERRQ(ierr);
	rank = (PetscInt)_rank;
	ierr = DMDAGetInfo(dm_mech,0,&M,&N,&P,&pI,&pJ,&pK,0,0,0,0,0,0);CHKERRQ(ierr);
	
	/* convert rank to rI,rJ,rK */
	rK  = rank / (pI*pJ);
	rIJ = rank - rK * pI*pJ;  
	rJ = rIJ / pI;
	//rI = rIJ - rJ*pI;
	
	surface_rank = PETSC_FALSE;
	if (rJ == (pJ-1)) {
		surface_rank = PETSC_TRUE;
	}
	//nsurface_ranks = pI * pK;

	/* check surface (I,J) and volume have matching (I,K) counts */
	{
		PetscInt sM,sN,sP;
		ierr = DMDAGetInfo(dm_msurf,0,&sM,&sN,&sP,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
		
		if (sP != 1) { SETERRQ1(comm,PETSC_ERR_SUP,"Surface DM (P=%D) must have P = 1 in K",sN); }
		if (sM != M) { SETERRQ2(comm,PETSC_ERR_SUP,"Surface DM (I:M=%D) doesn't overlap with volume DM (I:M=%D) in I",sM,M); }
		if (sN != P) { SETERRQ2(comm,PETSC_ERR_SUP,"Surface DM (J:N=%D) doesn't overlap with volume DM (K:P=%D) in J",sN,P); }
	}
	
	/* ----------------------------- */
	/* Get corner info for all ranks */
	ierr = DMDAGetOwnershipRanges(dm_mech,&lni,&lnj,&lnk);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*pI,&lsi);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*pK,&lsk);CHKERRQ(ierr);
	lsi[0] = 0;
	for (pi=1; pi<pI; pi++) {
		lsi[pi] = lsi[pi-1] + lni[pi-1]; 
	}
	lsk[0] = 0;
	for (pk=1; pk<pK; pk++) {
		lsk[pk] = lsk[pk-1] + lnk[pk-1]; 
	}
	
	/* Determine pieces i need on the volume mesh */
	ierr = DMDAGetCorners(dm_mech,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	nsurf_nodes = nx * nz;
	if (!surface_rank) {
		nsurf_nodes = 1;
	}
	ierr = PetscMalloc(sizeof(PetscInt)*nsurf_nodes*1,&surf_indices);CHKERRQ(ierr);

	ierr = DMDAGetInfo(dm_msurf,0,&M2d,&P2d,0,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	if (surface_rank) {
		j = N - 1;
		c = 0;
		for (k=sk; k<sk+nz; k++){
			for (i=si; i<si+nx; i++){
				PetscInt nidx;
				
				nidx = i + k*M2d*1;
				surf_indices[c] = nidx; 
				c++;
			}
		}
	} else {
		/* fetch something local */
		surf_indices[0] = si + sk*nx;
	}
	
	/* Create local comm_self vector to store height */
	ierr = VecCreate(PETSC_COMM_SELF,&height_self);CHKERRQ(ierr);
	ierr = VecSetSizes(height_self,PETSC_DECIDE,nsurf_nodes);CHKERRQ(ierr);
	ierr = VecSetFromOptions(height_self);CHKERRQ(ierr);
	
	/* Create natural ordering vector */
	/* create natural vector in i,j ordering for height */
	ierr = DMDACreateNaturalVector(dm_msurf,&height_natural);CHKERRQ(ierr);
	ierr = DMDAGlobalToNaturalBegin( dm_msurf,height,INSERT_VALUES, height_natural );CHKERRQ(ierr);
	ierr = DMDAGlobalToNaturalEnd(   dm_msurf,height,INSERT_VALUES, height_natural );CHKERRQ(ierr);
	
	/* Create scatter */
	ierr = ISCreateGeneral( comm,             nsurf_nodes, surf_indices, PETSC_USE_POINTER, &is_surf );CHKERRQ(ierr);
	ierr = ISCreateStride(  PETSC_COMM_SELF,  nsurf_nodes,0,1,&is_local);CHKERRQ(ierr);
	ierr = VecScatterCreate( height_natural,is_surf, height_self,is_local, &sctx );CHKERRQ(ierr);
	
	ierr = VecScatterBegin( sctx, height_natural,height_self,INSERT_VALUES, SCATTER_FORWARD );CHKERRQ(ierr);
	ierr = VecScatterEnd(   sctx, height_natural,height_self,INSERT_VALUES, SCATTER_FORWARD );CHKERRQ(ierr);
	
	/* update volume mesh */
	ierr = DMGetCoordinates(dm_mech,&coords_mech);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dm_mech,&cda_mech);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_mech,coords_mech,&LA_coords_mech);CHKERRQ(ierr);
	ierr = VecGetArray(height_self,&LA_height_self);CHKERRQ(ierr);
	
	/* identify the nodes on the surface of the mechanical domain */
	if (surface_rank) {
		j = N - 1;
		if ((j>=sj) && (j<sj+ny)) {
			for (k=sk; k<sk+nz; k++) {
				for (i=si; i<si+nx; i++) {
					LA_coords_mech[k][j][i].y = LA_height_self[(i-si)+(k-sk)*nx];
				}
			}
		}	
	}
	
	ierr = VecRestoreArray(height_self,&LA_height_self);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda_mech,coords_mech,&LA_coords_mech);CHKERRQ(ierr);
	
	/* insert surface values */
	ierr = DMDAUpdateGhostedCoordinates(dm_mech);CHKERRQ(ierr);
	
	/* clean up */
	ierr = VecScatterDestroy(&sctx);CHKERRQ(ierr);
	ierr = ISDestroy(&is_surf);CHKERRQ(ierr);
	ierr = ISDestroy(&is_local);CHKERRQ(ierr);
	
	ierr = VecDestroy(&height_natural);CHKERRQ(ierr);	
	ierr = VecDestroy(&height_self);CHKERRQ(ierr);
	ierr = PetscFree(surf_indices);CHKERRQ(ierr);
	ierr = PetscFree(lsi);CHKERRQ(ierr);
	ierr = PetscFree(lsk);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test_spm_utils_MPItoMPI"
PetscErrorCode test_spm_utils_MPItoMPI(DM dav)
{
	DM             dm_spmsurf;
	Vec            height;
	PetscInt       JMAX;
	PetscErrorCode ierr;
	
	
	ierr = DMDAGatherIKSurfaceDMDA(dav,&dm_spmsurf,&height);CHKERRQ(ierr);
	ierr = DMDAViewPetscVTS(dm_spmsurf,height,"surf_extraction_ic.vtk");CHKERRQ(ierr);

	ierr = VecShift(height,0.5);CHKERRQ(ierr);
	
	ierr = DMDAScatterIKSurfaceDMDA(dm_spmsurf,height,dav);CHKERRQ(ierr);
	
	ierr = VecDestroy(&height);CHKERRQ(ierr);
	ierr = DMDestroy(&dm_spmsurf);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dav,0,0,&JMAX,0,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,JMAX);CHKERRQ(ierr);
	ierr = DMDABilinearizeQ2Elements(dav);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
 Keep surface DM which overlaps with the volume DM (_dm_msurf)
 The surface DM used for the SPM should be _dm_msurf_ref.
*/
#undef __FUNCT__
#define __FUNCT__ "DMDAGatherIKNestedSurfaceDMDA"
PetscErrorCode DMDAGatherIKNestedSurfaceDMDA(DM dm_mech,PetscInt ref[],DM *_dm_msurf,Mat *_R,DM *_dm_msurf_ref,Vec *_elevation_ref)
{
	PetscErrorCode ierr;
	PetscInt       si,sj,sk,nx,ny,nz,M,N,P;
	PetscInt       si2d,sj2d,nx2d,ny2d,i,k;
	MPI_Comm       comm;
	DM             dm_red_spm,dm_surf,dm_surf_ref;
	DM             dm_red_spm_coords,dm_surf_coords;
	Vec            coords_surf,coords_red_spm;
	DMDACoor3d     ***LA_coords_red_spm;
	DMDACoor2d     **LA_coords_surf;
	PetscScalar    **LA_height;
	PetscMPIInt    rank;
	Vec            height,height_ref;
	Mat            R;
	
	
	PetscFunctionBegin;
	PetscObjectGetComm((PetscObject)dm_mech,&comm);
	ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dm_mech,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDACreate2d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,M,P, PETSC_DECIDE,PETSC_DECIDE, 1,1, 0,0,&dm_surf);CHKERRQ(ierr);	
	ierr = DMDASetUniformCoordinates(dm_surf, 0.0,1.0, 0.0,1.0, 0.,0.);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dm_surf,&si2d,&sj2d,NULL,&nx2d,&ny2d,NULL);CHKERRQ(ierr);
	
	/* fetch coordinates i need from the mechanical domain */
	ierr = DMDAGetCorners(dm_mech,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	
	//printf("rank %d: vol  si,ei=%.4d-%.4d: sj,ej=%.4d-%.4d: sk,ek=%.4d-%.4d \n", rank, si,si+nx,sj,sj+ny,sk,sk+nz);
	//printf("rank %d: surf si,ei=%.4d-%.4d:                  sj,ej=%.4d-%.4d \n", rank, si2d,si2d+nx2d,sj2d,sj2d+ny2d);
	
	ierr = DMDACreate3dRedundant(dm_mech,si2d,si2d+nx2d,N-1,N,sj2d,sj2d+ny2d, 1, &dm_red_spm);CHKERRQ(ierr);
	//sprintf(name,"surf%d.vtk",rank);
	//ierr = DMDAViewPetscVTS(da_red_spm,0,name);CHKERRQ(ierr);
	
	/* copy these values into my parallel surface mesh x,y,z (vol) => x,y (surf) */
	ierr = DMGetCoordinates(dm_surf,&coords_surf);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dm_red_spm,&coords_red_spm);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(dm_surf,&dm_surf_coords);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dm_red_spm,&dm_red_spm_coords);CHKERRQ(ierr);
	
	
	ierr = DMDAVecGetArray(dm_surf_coords,coords_surf,&LA_coords_surf);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dm_red_spm_coords,coords_red_spm,&LA_coords_red_spm);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(dm_surf,&si,&sk,NULL,&nx,&nz,NULL);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (i=si; i<si+nx; i++) {
			LA_coords_surf[k][i].x = LA_coords_red_spm[k-sk][0][i-si].x;
			LA_coords_surf[k][i].y = LA_coords_red_spm[k-sk][0][i-si].z;
		}
	}
	
	ierr = DMGetGlobalVector(dm_surf,&height);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(dm_surf,height,&LA_height);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (i=si; i<si+nx; i++) {
			LA_height[k][i] = LA_coords_red_spm[k-sk][0][i-si].y;
		}
	}
	ierr = DMDAVecRestoreArray(dm_surf,height,&LA_height);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(dm_red_spm_coords,coords_red_spm,&LA_coords_red_spm);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(dm_surf_coords,coords_surf,&LA_coords_surf);CHKERRQ(ierr);
	
	ierr = DMDestroy(&dm_red_spm);CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(dm_surf);CHKERRQ(ierr);
	
	
	/* refine DM */
	ierr = DMDASetRefinementFactor(dm_surf,ref[0],ref[1],1);CHKERRQ(ierr);
	dm_surf_ref = NULL;
	ierr = DMRefine(dm_surf,comm,&dm_surf_ref);CHKERRQ(ierr);
	if (!dm_surf_ref) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"No refinment occurred on surface DM");
	}
	
	/* refine height vector */
	ierr = DMCreateInterpolation(dm_surf,dm_surf_ref,&R,NULL);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(dm_surf_ref,&height_ref);CHKERRQ(ierr);
	ierr = MatInterpolate(R,height,height_ref);CHKERRQ(ierr);
	
	ierr = DMRestoreGlobalVector(dm_surf,&height);CHKERRQ(ierr);
	
	*_R             = R;
	*_elevation_ref = height_ref;
	*_dm_msurf      = dm_surf;
	*_dm_msurf_ref  = dm_surf_ref;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAScatterIKNestedSurfaceDMDA"
PetscErrorCode DMDAScatterIKNestedSurfaceDMDA(DM dm_msurf,DM dm_msurf_overlap,Mat R,Vec height,DM dm_mech)
{
	Vec            scale,height_vol;
	PetscErrorCode ierr;

	
	PetscFunctionBegin;
	ierr = DMGetGlobalVector(dm_msurf_overlap,&height_vol);CHKERRQ(ierr);
	
	/* Interpolate height to height_mech */
	ierr = DMCreateInterpolationScale(dm_msurf_overlap,dm_msurf,R,&scale);CHKERRQ(ierr);
	ierr = MatRestrict(R,height,height_vol);CHKERRQ(ierr);
	ierr = VecPointwiseMult(height_vol,height_vol,scale);CHKERRQ(ierr);
	ierr = VecDestroy(&scale);CHKERRQ(ierr);
	//ierr = DMDAViewPetscVTS(dm_msurf_overlap,height_vol,"surf_interp_ic.vtk");CHKERRQ(ierr);

	/* Use overlapping routine to scatter from surface -> volume */
	ierr = DMDAScatterIKSurfaceDMDA(dm_msurf_overlap,height_vol,dm_mech);CHKERRQ(ierr);
	
	ierr = DMRestoreGlobalVector(dm_msurf_overlap,&height_vol);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test_nested_spm_utils_MPItoMPI"
PetscErrorCode test_nested_spm_utils_MPItoMPI(DM dav)
{
	DM             dm_spmsurf,dm_spmsurf_overlap;
	Vec            height;
	PetscInt       JMAX;
	PetscInt       ref[] = { 10, 6 };
	Mat            R;
	PetscErrorCode ierr;
	
	
	ierr = DMDAGatherIKNestedSurfaceDMDA(dav,ref,&dm_spmsurf_overlap,&R,&dm_spmsurf,&height);CHKERRQ(ierr);
	ierr = DMDAViewPetscVTS(dm_spmsurf,height,"surf_extraction_ic.vtk");CHKERRQ(ierr);
	
	ierr = VecShift(height,0.5);CHKERRQ(ierr);
	
	ierr = DMDAScatterIKNestedSurfaceDMDA(dm_spmsurf,dm_spmsurf_overlap,R,height,dav);CHKERRQ(ierr);
	
	ierr = VecDestroy(&height);CHKERRQ(ierr);
	ierr = MatDestroy(&R);CHKERRQ(ierr);
	ierr = DMDestroy(&dm_spmsurf);CHKERRQ(ierr);
	ierr = DMDestroy(&dm_spmsurf_overlap);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dav,0,0,&JMAX,0,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,JMAX);CHKERRQ(ierr);
	ierr = DMDABilinearizeQ2Elements(dav);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}

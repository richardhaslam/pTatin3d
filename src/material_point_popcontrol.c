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
 **    filename:   material_point_popcontrol.c
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


#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "MPntStd_def.h"

#include "data_bucket.h"
#include "data_exchanger.h"
#include "element_type_Q2.h"
#include "element_utils_q2.h"
#include "dmda_element_q2p1.h"
#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "material_point_popcontrol.h"


#define MPPC_LOG_LEVEL 0 /* 0 - no logging; 1 - logging per mesh; 2 - logging per cell */

PetscLogEvent PTATIN_MaterialPointPopulationControlInsert;
PetscLogEvent PTATIN_MaterialPointPopulationControlRemove;


int sort_ComparePSortCtx(const void *dataA,const void *dataB)
{
	PSortCtx *pointA;
	PSortCtx *pointB;
	
	pointA = (PSortCtx*)dataA;
	pointB = (PSortCtx*)dataB;
	
	if (pointA->cell_index < pointB->cell_index) {
		return -1;
	} else if (pointA->cell_index > pointB->cell_index) {
		return 1;
	} else {
		return 0;
	}
}
void sort_PSortCx(const int _np, PSortCtx list[])
{
	PetscLogDouble t0,t1;
    size_t np;
	
    np = (size_t)_np;
	PetscTime(&t0);
	qsort( list, np, sizeof(PSortCtx), sort_ComparePSortCtx );
	PetscTime(&t1);
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PETSC_COMM_WORLD,"[LOG]  sort_PSortCx [npoints = %d] -> qsort %1.4e (sec)\n", _np,t1-t0 );
#endif
}


typedef struct _p_NNSortCtx {
	PetscInt point_index;
	double   coor[3];
	double   sep;
} NNSortCtx;

int sort_CompareNNSortCtx(const void *dataA,const void *dataB)
{
	NNSortCtx *pointA;
	NNSortCtx *pointB;
	
	pointA = (NNSortCtx*)dataA;
	pointB = (NNSortCtx*)dataB;
	
	if (pointA->sep < pointB->sep) {
		return -1;
	} else if (pointA->sep > pointB->sep) {
		return 1;
	} else {
		return 0;
	}
}
void sort_NNSortCx(const PetscInt _np, NNSortCtx list[])
{
    size_t np;
    np = (size_t)_np;
	qsort( list, np, sizeof(NNSortCtx), sort_CompareNNSortCtx );
}



#undef __FUNCT__  
#define __FUNCT__ "_find_min"
PetscErrorCode _find_min(double pos[],PetscInt point_count,double patch_point_coords[],PetscInt *idx)
{
	PetscInt p;
	double   dx,dy,dz;
	double   sep,min_sep;
	
	PetscFunctionBegin;

	min_sep = 1.0e32;
	*idx = -1;

	for (p=0; p<point_count; p++) {
		dx = (pos[0] - patch_point_coords[3*p+0]);
		dy = (pos[1] - patch_point_coords[3*p+1]);
		dz = (pos[2] - patch_point_coords[3*p+2]);
		sep = dx*dx + dy*dy + dz*dz;
		if (sep < min_sep) {
			min_sep = sep;
			*idx = p;
		}
	}
	
	PetscFunctionReturn(0);
}

/* 
 dist = AVD3dDistanceTest(p0,p1,p2);
 if dist > 0
   p2 is closer to p0 than p1
*/
static inline double AVD3dDistanceTest(double p0[],double p1[],double p2[])
{
	return (p1[0]+p2[0]-p0[0]-p0[0])*(p1[0]-p2[0]) + (p1[1]+p2[1]-p0[1]-p0[1])*(p1[1]-p2[1]) + (p1[2]+p2[2]-p0[2]-p0[2])*(p1[2]-p2[2]);
}

#undef __FUNCT__  
#define __FUNCT__ "_find_min_fast"
PetscErrorCode _find_min_fast(double pos[],PetscInt point_count,double patch_point_coords[],PetscInt *idx)
{
	PetscInt p,closest;
	double   dist;
	double   *p1,*p2;
	
	PetscFunctionBegin;
	
	closest = 0;
	p1 = &patch_point_coords[3*closest];
	
	for (p=1; p<point_count; p++) {
		
		p2 = &patch_point_coords[3*p];
		dist = AVD3dDistanceTest(pos,p1,p2);
		
		if (dist > 0.0) { /* p2 is closer than p1 */
			closest = p;
			//p1 = &patch_point_coords[3*closest];
			p1 = p2;
		}
	}

	*idx = closest;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "_find_min_sep_brute_force"
PetscErrorCode _find_min_sep_brute_force(double pos[],PetscInt point_count,NNSortCtx patch_points[],PetscInt *idx)
{
	PetscInt p;
	double   dx,dy,dz;
	double   sep,min_sep;
	
	PetscFunctionBegin;
	
	for (p=0; p<point_count; p++) {
		dx = (pos[0] - patch_points[p].coor[0]);
		dy = (pos[1] - patch_points[p].coor[1]);
		dz = (pos[2] - patch_points[p].coor[2]);
		sep = dx*dx + dy*dy + dz*dz;
		patch_points[p].sep = sep;
	}

	min_sep = 1.0e32;
	*idx = -1;
	for (p=0; p<point_count; p++) {
		sep = patch_points[p].sep;
		if (sep < min_sep) {
			min_sep = sep;
			*idx = p;
		}
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "_find_min_sep_qsort"
PetscErrorCode _find_min_sep_qsort(double pos[],PetscInt point_count,NNSortCtx patch_points[],PetscInt *idx)
{
	PetscInt p;
	double   dx,dy,dz;
	double   sep;
	
	PetscFunctionBegin;
	
	for (p=0; p<point_count; p++) {
		dx = (pos[0] - patch_points[p].coor[0]);
		dy = (pos[1] - patch_points[p].coor[1]);
		dz = (pos[2] - patch_points[p].coor[2]);
		sep = dx*dx + dy*dy + dz*dz;
		patch_points[p].sep = sep;
	}

	sort_NNSortCx(point_count,patch_points);
	
	*idx = 0;
	
	PetscFunctionReturn(0);
}


/*

 -std=gnu99 -O2 -Wall -Wno-unused-variable -Wstrict-aliasing -fstrict-aliasing -funroll-loops
 -mx 32 -my 32 -mz 32   -lattice_layout_Nx 4 -lattice_layout_Ny 4 -lattice_layout_Nz 4 -mp_popctrl_np_lower 100
 
_find_min     :  time_nn           = 1.2912e+00 (sec)
_find_min_fast:  time_nn           = 2.9685e+00 (sec) 
 
 apply_mppc_nn_patch(_find_min)
 
 [LOG]  sort_PSortCx [npoints = 2097152] -> qsort 1.5965e-01 (sec)
 [LOG]  cells with points < np_lower (32768) 
 [LOG]  cells with points > np_upper (0) 
 [LOG]  np_per_patch_max = 1728 
 [LOG]  time_nn           = 1.2912e+00 (sec)
 [LOG]  npoints_init      = 2097152 
 [LOG]  npoints_current-1 = 2359296 
 [LOG]  npoints_current-2 = 2359296 
 [LOG]  time(apply_mppc_nn_patch): 2.5713e+00 (sec)
 
 -mx 32 -my 32 -mz 32
 -lattice_layout_Nx 4 -lattice_layout_Ny 4 -lattice_layout_Nz 4 
 -mp_popctrl_np_lower 100 
 -mp_popctrl_nxp 4 -mp_popctrl_nyp 4 -mp_popctrl_nzp 4

 [LOG] MPPC_NearestNeighbourPatch: 
 [LOG]  sort_PSortCx [npoints = 2097152] -> qsort 1.5914e-01 (sec)
 [LOG]  cells with points < np_lower (32768) 
 [LOG]  cells with points > np_upper (0) 
 [LOG]  np_per_patch_max = 1728 
 [LOG]  time_nn           = 1.0283e+01 (sec)
 [LOG]  npoints_init      = 2097152 
 [LOG]  npoints_current-1 = 4194304 
 [LOG]  npoints_current-2 = 4194304 
 [LOG]  time(apply_mppc_nn_patch): 1.2710e+01 (sec)

 */
 
#undef __FUNCT__  
#define __FUNCT__ "apply_mppc_nn_patch"
PetscErrorCode apply_mppc_nn_patch(
														PetscInt ncells, PetscInt pcell_list[],
														PSortCtx plist[],
														PetscInt np_lower,
														PetscInt patch_extend,PetscInt nxp,PetscInt nyp,PetscInt nzp,PetscScalar perturb,DM da,DataBucket db)
{
	PetscInt        np_per_cell_max,mx,my,mz;
	PetscInt        c,i,j,k,cell_index_i,cell_index_j,cell_index_k,cidx2d,point_count,points_per_patch;
	PetscInt        p;
	const PetscInt  *elnidx;
	PetscInt        nel,nen;
	Vec             gcoords;
	PetscScalar     *LA_coords;
	PetscScalar     el_coords[Q2_NODES_PER_EL_3D*NSD];
	DM              cda;
	DataField       PField;
	int             Lnew,nx3,npoints_current,npoints_init;
	long int        cells_needing_new_points,cells_needing_new_points_g;
	double          *patch_point_coords;
	PetscInt        *patch_point_idx;
	PetscLogDouble  t0_nn,t1_nn,time_nn = 0.0;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	/* get mx,my from the da */
	ierr = DMDAGetLocalSizeElementQ2(da,&mx,&my,&mz);CHKERRQ(ierr);
	
    DataBucketGetSizes(db,&npoints_init,NULL,NULL);
	npoints_current = npoints_init;
    
	/* find max np_per_cell I will need */
    nx3 = (int)(nxp * nyp * nzp);
	np_per_cell_max = 0;
	cells_needing_new_points = 0;
	DataBucketGetSizes(db,&Lnew,NULL,NULL);
	for (c=0; c<nel; c++) {
		PetscInt points_per_cell;
		PetscInt points_per_patch;
		
		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		if (points_per_cell > np_lower) { continue; }

		cell_index_k = c / (mx*my);
		cidx2d = c - cell_index_k*(mx*my);
		cell_index_j = cidx2d / mx;
		cell_index_i = cidx2d - cell_index_j * mx;
		
		points_per_patch = 0;
		for ( k=cell_index_k - patch_extend; k<=cell_index_k + patch_extend; k++ ) {
			for ( j=cell_index_j - patch_extend; j<=cell_index_j + patch_extend; j++ ) {
				for ( i=cell_index_i - patch_extend; i<=cell_index_i + patch_extend; i++ ) {
					PetscInt patch_cell_id;
					
					if (i>=mx) { continue; }
					if (j>=my) { continue; }
					if (k>=mz) { continue; }
					if (i<0) { continue; }
					if (j<0) { continue; }
					if (k<0) { continue; }
					
					patch_cell_id = i + j * mx + k * mx*my;
					
					points_per_patch = points_per_patch + (pcell_list[patch_cell_id+1] - pcell_list[patch_cell_id]);
				}
			}
		}
		
		if (points_per_patch > np_per_cell_max) {
			np_per_cell_max = points_per_patch;
		}
		
		Lnew = Lnew + nx3;
		cells_needing_new_points++;
	}
	
#if (MPPC_LOG_LEVEL >= 1)
	printf("[LOG]  np_per_patch_max = %d \n", np_per_cell_max );
#endif	
	ierr = MPI_Allreduce( &cells_needing_new_points, &cells_needing_new_points_g, 1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD );CHKERRQ(ierr);
	if (cells_needing_new_points_g == 0) {
		//		PetscPrintf(PETSC_COMM_WORLD,"!! No population control required <global>!!\n");
		PetscFunctionReturn(0);
	}
	//PetscPrintf(PETSC_COMM_WORLD,"!! Population control required <global>!!\n");
	
	DataBucketSetSizes(db,Lnew,-1);

	ierr = PetscMalloc(sizeof(double)*3*np_per_cell_max,&patch_point_coords);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*np_per_cell_max,&patch_point_idx);CHKERRQ(ierr);
	
  /* setup for coords */
  ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	
	for (c=0; c<nel; c++) {
		PetscInt points_per_cell;
		
		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		if (points_per_cell > np_lower) { continue; }
		
		cell_index_k = c / (mx*my);
		cidx2d = c - cell_index_k*(mx*my);
		cell_index_j = cidx2d / mx;
		cell_index_i = cidx2d - cell_index_j * mx;
		
		/* load points */
		ierr = PetscMemzero( patch_point_coords, sizeof(double)*3*np_per_cell_max );CHKERRQ(ierr);
		ierr = PetscMemzero( patch_point_idx, sizeof(PetscInt)*np_per_cell_max );CHKERRQ(ierr);
		
		point_count = 0;
		
		DataFieldGetAccess(PField);
		for ( k=cell_index_k - patch_extend; k<=cell_index_k + patch_extend; k++ ) {
			for ( j=cell_index_j - patch_extend; j<=cell_index_j + patch_extend; j++ ) {
				for ( i=cell_index_i - patch_extend; i<=cell_index_i + patch_extend; i++ ) {
					PetscInt patch_cell_id;
					
					if (i>=mx) { continue; }
					if (j>=my) { continue; }
					if (k>=mz) { continue; }
					if (i<0) { continue; }
					if (j<0) { continue; }
					if (k<0) { continue; }
					
					patch_cell_id = i + j * mx + k * mx*my;
					points_per_patch = (pcell_list[patch_cell_id+1] - pcell_list[patch_cell_id]);
#if (MPPC_LOG_LEVEL >= 2)
					printf("[LOG]     patch(%d)-(%d,%d,%d) cell(%d)-(%d,%d,%d)  : ppcell = %d \n", c, cell_index_i,cell_index_j,cell_index_k, patch_cell_id,i,j,k,points_per_patch);
#endif					
					for (p=0; p<points_per_patch; p++) {
						MPntStd *marker_p;
						PetscInt pid, pid_unsorted;
						
						pid = pcell_list[patch_cell_id] + p;
						pid_unsorted = plist[pid].point_index;
						
						DataFieldAccessPoint(PField, pid_unsorted ,(void**)&marker_p);
						
						patch_point_coords[3*point_count+0] = marker_p->coor[0];
						patch_point_coords[3*point_count+1] = marker_p->coor[1];
						patch_point_coords[3*point_count+2] = marker_p->coor[2];
						patch_point_idx[point_count]        = pid_unsorted;
#if (MPPC_LOG_LEVEL >= 2)
						printf("[LOG]       patch(%d)/cell(%d) -> p(%d):p->wil,x,y,z = %d %1.4e %1.4e %1.4e \n", c, patch_cell_id, p,marker_p->wil, marker_p->coor[0],marker_p->coor[1],marker_p->coor[2] );
#endif
						point_count++;
					}
					
				}
			}
		}
		DataFieldRestoreAccess(PField);
#if (MPPC_LOG_LEVEL >= 2)
		printf("[LOG]  cell = %d: total points per patch = %d \n", c,point_count);
#endif
		
		/* create trial coordinates - find closest point */
		{
			PetscInt Nxp[NSD],pi,pj,pk,k;
			PetscScalar dxi,deta,dzeta;
			PetscInt marker_index;
			
			Nxp[0] = nxp;
			Nxp[1] = nyp;
			Nxp[2] = nzp;
			
			dxi    = 2.0/(PetscReal)Nxp[0];
			deta   = 2.0/(PetscReal)Nxp[1];
			dzeta  = 2.0/(PetscReal)Nxp[2];
			
			ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx[nen*c],LA_coords);CHKERRQ(ierr);
			for (pk=0; pk<Nxp[2]; pk++) {
				for (pj=0; pj<Nxp[1]; pj++) {
					for (pi=0; pi<Nxp[0]; pi++) {
						PetscInt idx;
						double   xip[NSD],xip_shift[NSD],xip_rand[NSD],xp_rand[NSD],Ni[Q2_NODES_PER_EL_3D];
						MPntStd  *marker_p,*marker_nearest;
						
						xip[0] = -1.0 + dxi    * (pi + 0.5);
						xip[1] = -1.0 + deta   * (pj + 0.5);
						xip[2] = -1.0 + dzeta  * (pk + 0.5);
						
						/* random between -0.5 <= shift <= 0.5 */
						xip_shift[0] = 1.0*(rand()/(RAND_MAX+1.0)) - 0.5;
						xip_shift[1] = 1.0*(rand()/(RAND_MAX+1.0)) - 0.5;
						xip_shift[2] = 1.0*(rand()/(RAND_MAX+1.0)) - 0.5;
						
						xip_rand[0] = xip[0] + perturb * dxi    * xip_shift[0];
						xip_rand[1] = xip[1] + perturb * deta   * xip_shift[1];
						xip_rand[2] = xip[2] + perturb * dzeta  * xip_shift[2];
						
						if (fabs(xip_rand[0]) > 1.0) {
							SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"fabs(x-point coord) greater than 1.0");
						}
						if (fabs(xip_rand[1]) > 1.0) {
							SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"fabs(y-point coord) greater than 1.0");
						}
						if (fabs(xip_rand[2]) > 1.0) {
							SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"fabs(z-point coord) greater than 1.0");
						}
						
						P3D_ConstructNi_Q2_3D(xip_rand,Ni);
						
						xp_rand[0] = xp_rand[1] = xp_rand[2] = 0.0;
						for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
							xp_rand[0] += Ni[k] * el_coords[NSD*k+0];
							xp_rand[1] += Ni[k] * el_coords[NSD*k+1];
							xp_rand[2] += Ni[k] * el_coords[NSD*k+2];
						}
						
						/* locate nearest point */
						PetscTime(&t0_nn);
						ierr = _find_min(xp_rand,point_count,patch_point_coords,&idx);CHKERRQ(ierr);
						//ierr = _find_min_fast(xp_rand,point_count,patch_point_coords,&idx);CHKERRQ(ierr);
						PetscTime(&t1_nn);
						time_nn += (t1_nn - t0_nn);
						
						marker_index = patch_point_idx[ idx ];
						
						DataBucketCopyPoint(db,(int)marker_index, db,npoints_current);
						
						DataFieldGetAccess(PField);
						
						DataFieldAccessPoint(PField,npoints_current,(void**)&marker_p);
						DataFieldAccessPoint(PField,(int)marker_index,(void**)&marker_nearest);
						
						marker_p->phase   = marker_nearest->phase;
						
						marker_p->coor[0] = xp_rand[0];
						marker_p->coor[1] = xp_rand[1];
						marker_p->coor[2] = xp_rand[2];
						marker_p->wil     = (int)c;
						marker_p->xi[0]   = xip_rand[0];
						marker_p->xi[1]   = xip_rand[1];
						marker_p->xi[2]   = xip_rand[2];
						
						DataFieldRestoreAccess(PField);
						
						
						npoints_current++;
					}
				}
			}
			
		}
		
		
	}
  ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
	ierr = PetscFree(patch_point_coords);CHKERRQ(ierr);
	ierr = PetscFree(patch_point_idx);CHKERRQ(ierr);
	
	DataBucketGetSizes(db,&Lnew,NULL,NULL);
	
	ierr = SwarmMPntStd_AssignUniquePointIdentifiers(PetscObjectComm((PetscObject)da),db,npoints_init,Lnew);CHKERRQ(ierr);
	
#if (MPPC_LOG_LEVEL >= 1)
	printf("[LOG]  time_nn           = %1.4e (sec)\n", time_nn);
	printf("[LOG]  npoints_init      = %d \n", npoints_init);
	printf("[LOG]  npoints_current-1 = %d \n", npoints_current);
	printf("[LOG]  npoints_current-2 = %d \n", Lnew);
#endif
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "apply_mppc_nn_patch2"
PetscErrorCode apply_mppc_nn_patch2(
																	 PetscInt ncells, PetscInt pcell_list[],
																	 PetscInt np, PSortCtx plist[],
																	 PetscInt np_lower,
																	 PetscInt patch_extend,PetscInt nxp,PetscInt nyp,PetscInt nzp,PetscScalar perturb,DM da,DataBucket db)
{
	PetscInt        np_per_cell_max,mx,my,mz;
	PetscInt        c,i,j,k,cell_index_i,cell_index_j,cell_index_k,cidx2d,point_count,points_per_patch;
	PetscInt        p;
	const PetscInt  *elnidx;
	PetscInt        nel,nen;
	Vec             gcoords;
	PetscScalar     *LA_coords;
	PetscScalar     el_coords[Q2_NODES_PER_EL_3D*NSD];
	DM              cda;
	DataField       PField;
	int             Lnew,npoints_current,npoints_init;
	long int        cells_needing_new_points,cells_needing_new_points_g;
	NNSortCtx       *patch_points;
	PetscLogDouble  t0_nn,t1_nn,time_nn = 0.0;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	/* get mx,my from the da */
	ierr = DMDAGetLocalSizeElementQ2(da,&mx,&my,&mz);CHKERRQ(ierr);
	
	/* find max np_per_cell I will need */
	np_per_cell_max = 0;
	cells_needing_new_points = 0;
	DataBucketGetSizes(db,&Lnew,NULL,NULL);
	for (c=0; c<nel; c++) {
		PetscInt points_per_cell;
		PetscInt points_per_patch;
		
		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		if (points_per_cell > np_lower) { continue; }
		
		cell_index_k = c / (mx*my);
		cidx2d = c - cell_index_k*(mx*my);
		cell_index_j = cidx2d / mx;
		cell_index_i = cidx2d - cell_index_j * mx;
		
		points_per_patch = 0;
		for ( k=cell_index_k - patch_extend; k<=cell_index_k + patch_extend; k++ ) {
			for ( j=cell_index_j - patch_extend; j<=cell_index_j + patch_extend; j++ ) {
				for ( i=cell_index_i - patch_extend; i<=cell_index_i + patch_extend; i++ ) {
					PetscInt patch_cell_id;
					
					if (i>=mx) { continue; }
					if (j>=my) { continue; }
					if (k>=mz) { continue; }
					if (i<0) { continue; }
					if (j<0) { continue; }
					if (k<0) { continue; }
					
					patch_cell_id = i + j * mx + k * mx*my;
					
					points_per_patch = points_per_patch + (pcell_list[patch_cell_id+1] - pcell_list[patch_cell_id]);
				}
			}
		}
		
		if (points_per_patch > np_per_cell_max) {
			np_per_cell_max = points_per_patch;
		}
		
		Lnew = Lnew + nxp * nyp * nzp;
		cells_needing_new_points++;
	}
	
#if (MPPC_LOG_LEVEL >= 1)
	printf("[LOG]  np_per_patch_max = %d \n", np_per_cell_max );
#endif	
	ierr = MPI_Allreduce( &cells_needing_new_points, &cells_needing_new_points_g, 1, MPI_LONG, MPI_SUM, PetscObjectComm((PetscObject)da) );CHKERRQ(ierr);
	if (cells_needing_new_points_g==0) {
		//		PetscPrintf(PETSC_COMM_WORLD,"!! No population control required <global>!!\n");
		PetscFunctionReturn(0);
	}
	//PetscPrintf(PETSC_COMM_WORLD,"!! Population control required <global>!!\n");
	
	DataBucketSetSizes(db,Lnew,-1);
	
	ierr = PetscMalloc(sizeof(NNSortCtx)*np_per_cell_max,&patch_points);CHKERRQ(ierr);
	
  /* setup for coords */
  ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	
	npoints_init    = np;
	npoints_current = np;
	for (c=0; c<nel; c++) {
		PetscInt points_per_cell;
		
		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		if (points_per_cell > np_lower) { continue; }
		
		cell_index_k = c / (mx*my);
		cidx2d = c - cell_index_k*(mx*my);
		cell_index_j = cidx2d / mx;
		cell_index_i = cidx2d - cell_index_j * mx;
		
		/* load points */
		ierr = PetscMemzero( patch_points, sizeof(NNSortCtx)*np_per_cell_max );CHKERRQ(ierr);
		
		point_count = 0;
		
		DataFieldGetAccess(PField);
		for ( k=cell_index_k - patch_extend; k<=cell_index_k + patch_extend; k++ ) {
			for ( j=cell_index_j - patch_extend; j<=cell_index_j + patch_extend; j++ ) {
				for ( i=cell_index_i - patch_extend; i<=cell_index_i + patch_extend; i++ ) {
					PetscInt patch_cell_id;
					
					if (i>=mx) { continue; }
					if (j>=my) { continue; }
					if (k>=mz) { continue; }
					if (i<0) { continue; }
					if (j<0) { continue; }
					if (k<0) { continue; }
					
					patch_cell_id = i + j * mx + k * mx*my;
					points_per_patch = (pcell_list[patch_cell_id+1] - pcell_list[patch_cell_id]);
#if (MPPC_LOG_LEVEL >= 2)
					printf("[LOG]     patch(%d)-(%d,%d,%d) cell(%d)-(%d,%d,%d)  : ppcell = %d \n", c, cell_index_i,cell_index_j,cell_index_k, patch_cell_id,i,j,k,points_per_patch);
#endif					
					for (p=0; p<points_per_patch; p++) {
						MPntStd *marker_p;
						PetscInt pid, pid_unsorted;
						
						pid = pcell_list[patch_cell_id] + p;
						pid_unsorted = plist[pid].point_index;
						
						DataFieldAccessPoint(PField, pid_unsorted ,(void**)&marker_p);
						
						patch_points[point_count].coor[0] = marker_p->coor[0];
						patch_points[point_count].coor[0] = marker_p->coor[1];
						patch_points[point_count].coor[0] = marker_p->coor[2];
						patch_points[point_count].point_index = pid_unsorted;
#if (MPPC_LOG_LEVEL >= 2)
						printf("[LOG]       patch(%d)/cell(%d) -> p(%d):p->wil,x,y,z = %d %1.4e %1.4e %1.4e \n", c, patch_cell_id, p,marker_p->wil, marker_p->coor[0],marker_p->coor[1],marker_p->coor[2] );
#endif
						point_count++;
					}
					
				}
			}
		}
		DataFieldRestoreAccess(PField);
#if (MPPC_LOG_LEVEL >= 2)
		printf("[LOG]  cell = %d: total points per patch = %d \n", c,point_count);
#endif
		
		/* create trial coordinates - find closest point */
		{
			PetscInt Nxp[NSD],pi,pj,pk,k;
			PetscScalar dxi,deta,dzeta;
			PetscInt marker_index;
			
			Nxp[0] = nxp;
			Nxp[1] = nyp;
			Nxp[2] = nzp;
			
			dxi    = 2.0/(PetscReal)Nxp[0];
			deta   = 2.0/(PetscReal)Nxp[1];
			dzeta  = 2.0/(PetscReal)Nxp[2];
			
			ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx[nen*c],LA_coords);CHKERRQ(ierr);
			for (pk=0; pk<Nxp[2]; pk++) {
				for (pj=0; pj<Nxp[1]; pj++) {
					for (pi=0; pi<Nxp[0]; pi++) {
						int idx;
						double xip[NSD],xip_shift[NSD],xip_rand[NSD],xp_rand[NSD],Ni[Q2_NODES_PER_EL_3D];
						MPntStd *marker_p,*marker_nearest;
						
						xip[0] = -1.0 + dxi    * (pi + 0.5);
						xip[1] = -1.0 + deta   * (pj + 0.5);
						xip[2] = -1.0 + dzeta  * (pk + 0.5);
						
						/* random between -0.5 <= shift <= 0.5 */
						xip_shift[0] = 1.0*(rand()/(RAND_MAX+1.0)) - 0.5;
						xip_shift[1] = 1.0*(rand()/(RAND_MAX+1.0)) - 0.5;
						xip_shift[2] = 1.0*(rand()/(RAND_MAX+1.0)) - 0.5;
						
						xip_rand[0] = xip[0] + perturb * dxi    * xip_shift[0];
						xip_rand[1] = xip[1] + perturb * deta   * xip_shift[1];
						xip_rand[2] = xip[2] + perturb * dzeta  * xip_shift[2];
						
						if (fabs(xip_rand[0]) > 1.0) {
							SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"fabs(x-point coord) greater than 1.0");
						}
						if (fabs(xip_rand[1]) > 1.0) {
							SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"fabs(y-point coord) greater than 1.0");
						}
						if (fabs(xip_rand[2]) > 1.0) {
							SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"fabs(z-point coord) greater than 1.0");
						}
						
						P3D_ConstructNi_Q2_3D(xip_rand,Ni);
						
						xp_rand[0] = xp_rand[1] = xp_rand[2] = 0.0;
						for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
							xp_rand[0] += Ni[k] * el_coords[NSD*k+0];
							xp_rand[1] += Ni[k] * el_coords[NSD*k+1];
							xp_rand[2] += Ni[k] * el_coords[NSD*k+2];
						}
						
						/* locate nearest point */
						PetscTime(&t0_nn);
						//ierr = _find_min_sep_brute_force(xp_rand,point_count,patch_points,&idx);CHKERRQ(ierr);
						ierr = _find_min_sep_qsort(xp_rand,point_count,patch_points,&idx);CHKERRQ(ierr);
						PetscTime(&t1_nn);
						time_nn += (t1_nn - t0_nn);
						
						marker_index = patch_points[ idx ].point_index;
						
						DataBucketCopyPoint(db,marker_index, db,npoints_current);
						
						DataFieldGetAccess(PField);
						
						DataFieldAccessPoint(PField,npoints_current,(void**)&marker_p);
						DataFieldAccessPoint(PField,marker_index,(void**)&marker_nearest);
						
						marker_p->phase   = marker_nearest->phase;
						
						marker_p->coor[0] = xp_rand[0];
						marker_p->coor[1] = xp_rand[1];
						marker_p->coor[2] = xp_rand[2];
						marker_p->wil     = c;
						marker_p->xi[0]   = xip_rand[0];
						marker_p->xi[1]   = xip_rand[1];
						marker_p->xi[2]   = xip_rand[2];
						
						DataFieldRestoreAccess(PField);
						
						
						npoints_current++;
					}
				}
			}
			
		}
		
		
	}
  ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
	ierr = PetscFree(patch_points);CHKERRQ(ierr);
	
	DataBucketGetSizes(db,&Lnew,NULL,NULL);
	
	ierr = SwarmMPntStd_AssignUniquePointIdentifiers(PetscObjectComm((PetscObject)da),db,npoints_init,Lnew);CHKERRQ(ierr);
	
#if (MPPC_LOG_LEVEL >= 1)
	printf("[LOG]  time_nn           = %1.4e (sec)\n", time_nn);
	printf("[LOG]  npoints_init      = %d \n", npoints_init);
	printf("[LOG]  npoints_current-1 = %d \n", npoints_current);
	printf("[LOG]  npoints_current-2 = %d \n", Lnew);
#endif
	
	PetscFunctionReturn(0);
}

/*
 
 1) create list of markers, sorted by element
 2) for any cell c, which has less than np_lower points
 3)  get cell patch bounds using nodal coordinates
 4)  assemble avd using marker global coordinates
 5)  inject nxp x nyp into cell
 
 */
#undef __FUNCT__  
#define __FUNCT__ "MPPC_NearestNeighbourPatch"
PetscErrorCode MPPC_NearestNeighbourPatch(PetscInt np_lower,PetscInt np_upper,PetscInt patch_extend,PetscInt nxp,PetscInt nyp,PetscInt nzp,PetscScalar pertub,DM da,DataBucket db)
{
	PetscInt        *pcell_list;
	PSortCtx        *plist;
	int             p,npoints;
	PetscInt        tmp,c,count,cells_np_lower,cells_np_upper;
	const PetscInt  *elnidx;
	PetscInt        nel,nen;
	DataField       PField;
	PetscLogDouble  t0,t1;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_MaterialPointPopulationControlInsert,0,0,0,0);CHKERRQ(ierr);
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG] %s: \n", __FUNCTION__);
#endif	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	ierr = PetscMalloc( sizeof(PetscInt)*(nel+1),&pcell_list );CHKERRQ(ierr);
	
	DataBucketGetSizes(db,&npoints,NULL,NULL);
	ierr = PetscMalloc( sizeof(PSortCtx)*(npoints), &plist);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p;
		
		DataFieldAccessPoint(PField,p,(void**)&marker_p);
		plist[p].point_index = (PetscInt)p;
		plist[p].cell_index  = (PetscInt)marker_p->wil;
	}
	DataFieldRestoreAccess(PField);
	
	sort_PSortCx(npoints,plist);
	
	/* sum points per cell */
	ierr = PetscMemzero( pcell_list,sizeof(PetscInt)*(nel+1) );CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		pcell_list[ plist[p].cell_index ]++;
	}
	
	/* create offset list */
	count = 0;
	for (c=0; c<nel; c++) {
		tmp = pcell_list[c];
		pcell_list[c] = count;
		count = count + tmp;
	}
	pcell_list[c] = count;
	
	cells_np_lower = 0;
	cells_np_upper = 0;
	for (c=0; c<nel; c++) {
		PetscInt points_per_cell = pcell_list[c+1] - pcell_list[c];
        
		if (points_per_cell <= np_lower) { cells_np_lower++; }
		if (points_per_cell > np_upper) { cells_np_upper++; }
	}
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  %d cells with points < np_lower (%d) \n", cells_np_lower,np_lower );
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  %d cells with points > np_upper (%d) \n", cells_np_upper,np_upper);
#endif
	
	/* apply point injection routine */
	PetscTime(&t0);
	ierr = apply_mppc_nn_patch(
											nel, pcell_list,
											plist,
											np_lower,
											patch_extend, nxp,nyp,nzp, pertub, da,db);CHKERRQ(ierr);
	PetscTime(&t1);
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  time(apply_mppc_nn_patch): %1.4e (sec)\n", t1-t0);
#endif
	
	
	ierr = PetscFree(plist);CHKERRQ(ierr);
	ierr = PetscFree(pcell_list);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(PTATIN_MaterialPointPopulationControlInsert,0,0,0,0);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MPPC_SimpleRemoval"
PetscErrorCode MPPC_SimpleRemoval(PetscInt np_upper,DM da,DataBucket db,PetscBool reverse_order_removal)
{
	PetscInt        *cell_count,count;
	int             p,npoints;
	PetscInt        c,nel,nen;
	const PetscInt  *elnidx;
	DataField       PField;
	PetscLogDouble  t0,t1;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	ierr = PetscLogEventBegin(PTATIN_MaterialPointPopulationControlRemove,0,0,0,0);CHKERRQ(ierr);
	
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG] %s: \n", __FUNCTION__);
#endif	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	ierr = PetscMalloc( sizeof(PetscInt)*(nel),&cell_count );CHKERRQ(ierr);
	ierr = PetscMemzero( cell_count, sizeof(PetscInt)*(nel) );CHKERRQ(ierr);
	
	DataBucketGetSizes(db,&npoints,NULL,NULL);

	/* compute number of points per cell */
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	DataFieldGetAccess(PField);
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p;

		DataFieldAccessPoint(PField,p,(void**)&marker_p);
		if (marker_p->wil < 0) { continue; }
		
		cell_count[ marker_p->wil ]++;
	}
	DataFieldRestoreAccess(PField);

	
	count = 0;
	for (c=0; c<nel; c++) {
		if (cell_count[c] > np_upper) {
			count++;
		}
	}
	
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  %d cells with points > np_upper (%d) \n", count, np_upper);
#endif
	
	
	if (count == 0) {
		ierr = PetscFree(cell_count);CHKERRQ(ierr);
		ierr = PetscLogEventEnd(PTATIN_MaterialPointPopulationControlRemove,0,0,0,0);CHKERRQ(ierr);
		PetscFunctionReturn(0);
	}

	PetscTime(&t0);

	if (!reverse_order_removal) {
		/* remove points from cells with excessive number */
		DataFieldGetAccess(PField);
		for (p=0; p<npoints; p++) {
			MPntStd *marker_p;
			int wil;
			
			DataFieldAccessPoint(PField,p,(void**)&marker_p);
			wil = marker_p->wil;
			
			if (cell_count[wil] > np_upper) {
				DataBucketRemovePointAtIndex(db,p);
				
				DataBucketGetSizes(db,&npoints,0,0); /* you need to update npoints as the list size decreases! */
				p--; /* check replacement point */
				cell_count[wil]--;
			}
		}
		DataFieldRestoreAccess(PField);
	}
	
	/* scan in reverse order so that most recent points added to list will be removed as a priority */
/*
	if (reverse_order_removal) {
		DataFieldGetAccess(PField);
		for (p=npoints-1; p>=0; p--) {
			MPntStd *marker_p;
			int wil;
			
			DataFieldAccessPoint(PField,p,(void**)&marker_p);
			wil = marker_p->wil;
			
			if (cell_count[wil] > np_upper) {
				DataBucketRemovePointAtIndex(db,p);
				
				DataBucketGetSizes(db,&npoints,0,0); // you need to update npoints as the list size decreases! //
				cell_count[wil]--;
			}
		}
		DataFieldRestoreAccess(PField);	
	}
*/
	if (reverse_order_removal) {
		int remove;

		remove = 1;
		while (remove == 1) {
		
			remove = 0;

			DataBucketGetSizes(db,&npoints,0,0); // you need to update npoints as the list size decreases! //

			DataFieldGetAccess(PField);
			for (p=npoints-1; p>=0; p--) {
				MPntStd *marker_p;
				int wil;
				
				DataFieldAccessPoint(PField,p,(void**)&marker_p);
				wil = marker_p->wil;
				
				if (cell_count[wil] > np_upper) {
					DataBucketRemovePointAtIndex(db,p);
					cell_count[wil]--;
					
					remove = 1;
					break;
				}
			}
			DataFieldRestoreAccess(PField);	
		}
		
	}
	
	PetscTime(&t1);
	
	
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  time(MPPC_SimpleRemoval): %1.4e (sec)\n", t1-t0);
#endif
	
	ierr = PetscFree(cell_count);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(PTATIN_MaterialPointPopulationControlRemove,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MaterialPointPopulationControl_v1"
PetscErrorCode MaterialPointPopulationControl_v1(pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscInt       np_lower,np_upper,patch_extent,nxp,nyp,nzp;
	PetscScalar    perturb;
	PetscBool      flg;
	DataBucket     db;
	PetscBool      reverse_order_removal;
    MPI_Comm       comm;
	
	PetscFunctionBegin;

	/* options for control number of points per cell */
	np_lower = 0;
	np_upper = 60;
	ierr = PetscOptionsGetInt(NULL,"-mp_popctrl_np_lower",&np_lower,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,"-mp_popctrl_np_upper",&np_upper,&flg);CHKERRQ(ierr);

	/* options for injection of markers */
	nxp = 2;
	nyp = 2;
	nzp = 2;
	ierr = PetscOptionsGetInt(NULL,"-mp_popctrl_nxp",&nxp,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,"-mp_popctrl_nyp",&nyp,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,"-mp_popctrl_nzp",&nzp,&flg);CHKERRQ(ierr);
	
	perturb = 0.1;
	ierr = PetscOptionsGetScalar(NULL,"-mp_popctrl_perturb",&perturb,&flg);CHKERRQ(ierr);
	patch_extent = 1;
	ierr = PetscOptionsGetInt(NULL,"-mp_popctrl_patch_extent",&patch_extent,&flg);CHKERRQ(ierr);

	ierr = pTatinGetMaterialPoints(ctx,&db,NULL);CHKERRQ(ierr);
    ierr = PetscObjectGetComm((PetscObject)ctx->stokes_ctx->dav,&comm);CHKERRQ(ierr);
	
	/* insertion */
#if (MPPC_LOG_LEVEL >= 1)
    {
        long int np_g;
        DataBucketGetGlobalSizes(comm,db,&np_g,NULL,NULL);
        PetscPrintf(comm,"[LOG]  total markers before population control (%ld) \n", np_g );
    }
#endif
	
	if (np_lower != -1) {
        ierr = MPPC_NearestNeighbourPatch(np_lower,np_upper,patch_extent,nxp,nyp,nzp,perturb,ctx->stokes_ctx->dav,db);CHKERRQ(ierr);
	}
#if (MPPC_LOG_LEVEL >= 1)
    {
        long int np_g;
        DataBucketGetGlobalSizes(comm,db,&np_g,NULL,NULL);
        PetscPrintf(comm,"[LOG]  total markers after INJECTION (%ld) \n", np_g );
    }
#endif
	
	/* removal */
	if (np_upper != -1) {
		reverse_order_removal = PETSC_TRUE;
		ierr = MPPC_SimpleRemoval(np_upper,ctx->stokes_ctx->dav,db,reverse_order_removal);CHKERRQ(ierr);
	}	

#if (MPPC_LOG_LEVEL >= 1)
    {
        long int np_g;
        DataBucketGetGlobalSizes(comm,db,&np_g,NULL,NULL);
        PetscPrintf(comm,"[LOG]  total markers after DELETION (%ld) \n", np_g );
    }
#endif
	
	PetscFunctionReturn(0);
}


/*
 Assign all markers with phase = MATERIAL_POINT_PHASE_UNASSIGNED to closest phase
*/
#undef __FUNCT__  
#define __FUNCT__ "apply_mppc_region_assignment"
PetscErrorCode apply_mppc_region_assignment(
																	 PetscInt nel, PetscInt cell_count[], PetscInt pcell_list[],
																	 PetscInt np, PSortCtx plist[],
																	 PetscInt patch_extend,DM da,DataBucket db)
{
	PetscInt        np_per_cell_max,mx,my,mz;
	PetscInt        c,i,j,k,cell_index_i,cell_index_j,cell_index_k,cidx2d,point_count;
	PetscInt        points_per_cell,points_per_patch;
	PetscInt        p;
	DataField       PField;
	double          *patch_point_coords;
	int             *patch_point_idx;
	PetscLogDouble  t0_nn,t1_nn,time_nn = 0.0;
	int             points_assigned = 0;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	
	/* get mx,my from the da */
	ierr = DMDAGetLocalSizeElementQ2(da,&mx,&my,&mz);CHKERRQ(ierr);

	/* find max np_per_cell I will need */
	np_per_cell_max = 0;
	for (c=0; c<nel; c++) {
		
		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		if (cell_count[c] == 0) { continue; }
		
		cell_index_k = c / (mx*my);
		cidx2d = c - cell_index_k*(mx*my);
		cell_index_j = cidx2d / mx;
		cell_index_i = cidx2d - cell_index_j * mx;
		
		points_per_patch = 0;
		for ( k=cell_index_k - patch_extend; k<=cell_index_k + patch_extend; k++ ) {
			for ( j=cell_index_j - patch_extend; j<=cell_index_j + patch_extend; j++ ) {
				for ( i=cell_index_i - patch_extend; i<=cell_index_i + patch_extend; i++ ) {
					PetscInt patch_cell_id;
					
					if (i >= mx) { continue; }
					if (j >= my) { continue; }
					if (k >= mz) { continue; }
					if (i < 0) { continue; }
					if (j < 0) { continue; }
					if (k < 0) { continue; }
					
					patch_cell_id = i + j * mx + k * mx*my;
					
					points_per_patch = points_per_patch + (pcell_list[patch_cell_id+1] - pcell_list[patch_cell_id]);
				}
			}
		}
		
		if (points_per_patch > np_per_cell_max) {
			np_per_cell_max = points_per_patch;
		}
	}
	
#if (MPPC_LOG_LEVEL >= 1)
	printf("[LOG]  np_per_patch_max = %d \n", np_per_cell_max );
#endif	
	
	ierr = PetscMalloc(sizeof(double)*3*np_per_cell_max,&patch_point_coords);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(int)*np_per_cell_max,&patch_point_idx);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	
	for (c=0; c<nel; c++) {
		
		/* if cell doesn't contain any points which need re-assignment - skip */
		if (cell_count[c] == 0) { continue; }

		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		cell_index_k = c / (mx*my);
		cidx2d = c - cell_index_k*(mx*my);
		cell_index_j = cidx2d / mx;
		cell_index_i = cidx2d - cell_index_j * mx;
		
		/* load points from the patch into a list - only load points with an assigned phase index */
		ierr = PetscMemzero( patch_point_coords, sizeof(double)*3*np_per_cell_max );CHKERRQ(ierr);
		ierr = PetscMemzero( patch_point_idx, sizeof(int)*np_per_cell_max );CHKERRQ(ierr);
		
		point_count = 0;
		
		DataFieldGetAccess(PField);
		for ( k=cell_index_k - patch_extend; k<=cell_index_k + patch_extend; k++ ) {
			for ( j=cell_index_j - patch_extend; j<=cell_index_j + patch_extend; j++ ) {
				for ( i=cell_index_i - patch_extend; i<=cell_index_i + patch_extend; i++ ) {
					PetscInt patch_cell_id;
					
					if (i >= mx) { continue; }
					if (j >= my) { continue; }
					if (k >= mz) { continue; }
					if (i < 0) { continue; }
					if (j < 0) { continue; }
					if (k < 0) { continue; }
					
					patch_cell_id = i + j * mx + k * mx*my;
					points_per_patch = (pcell_list[patch_cell_id+1] - pcell_list[patch_cell_id]);
#if (MPPC_LOG_LEVEL >= 2)
					printf("[LOG]     patch(%d)-(%d,%d,%d) cell(%d)-(%d,%d,%d)  : ppcell = %d \n", c, cell_index_i,cell_index_j,cell_index_k, patch_cell_id,i,j,k,points_per_patch);
#endif					
					for (p=0; p<points_per_patch; p++) {
						MPntStd *marker_p;
						PetscInt pid, pid_unsorted;
						
						pid = pcell_list[patch_cell_id] + p;
						pid_unsorted = plist[pid].point_index;
						
						DataFieldAccessPoint(PField, pid_unsorted ,(void**)&marker_p);
						
						/* skip markers from patch which need to be assigned */
						if (marker_p->phase == MATERIAL_POINT_PHASE_UNASSIGNED) { continue; }
						
						patch_point_coords[3*point_count+0] = marker_p->coor[0];
						patch_point_coords[3*point_count+1] = marker_p->coor[1];
						patch_point_coords[3*point_count+2] = marker_p->coor[2];
						patch_point_idx[point_count]        = pid_unsorted;
#if (MPPC_LOG_LEVEL >= 2)
						printf("[LOG]       patch(%d)/cell(%d) -> p(%d):p->wil,x,y,z = %d %1.4e %1.4e %1.4e \n", c, patch_cell_id, p,marker_p->wil, marker_p->coor[0],marker_p->coor[1],marker_p->coor[2] );
#endif
						point_count++;
					}
					
				}
			}
		}
		DataFieldRestoreAccess(PField);
#if (MPPC_LOG_LEVEL >= 2)
		printf("[LOG]  cell = %d: total points per patch = %d \n", c,point_count);
#endif
		
		/* traverse points in this cell with phase = MATERIAL_POINT_PHASE_UNASSIGNED and find closest point */
		points_per_cell = pcell_list[c+1] - pcell_list[c];

		for (p=0; p<points_per_cell; p++) {
			MPntStd  *marker_p,*marker_nearest;
			double   *pos_p;
			int      pid,pid_unsorted,nearest_idx,marker_index;
			double   xp_orig[3],xip_orig[3];
			long int pid_orig;
			
			pid = pcell_list[c] + p;
			pid_unsorted = plist[pid].point_index;
			
			DataFieldGetAccess(PField);
			DataFieldAccessPoint(PField,pid_unsorted,(void**)&marker_p);
			
			/* if marker is assigned - skip */
			if (marker_p->phase != MATERIAL_POINT_PHASE_UNASSIGNED) { 
				DataFieldRestoreAccess(PField);
				continue; 
			}

			pid_orig    = marker_p->pid;
			xp_orig[0]  = marker_p->coor[0];
			xp_orig[1]  = marker_p->coor[1];
			xp_orig[2]  = marker_p->coor[2];
			xip_orig[0] = marker_p->xi[0];
			xip_orig[1] = marker_p->xi[1];
			xip_orig[2] = marker_p->xi[2];
			
			
#if (MPPC_LOG_LEVEL >= 2)
			PetscPrintf(PETSC_COMM_SELF,"[LOG]  cell(%d) point(%d) is un-assigned\n",c,pid_unsorted);
#endif
			
			pos_p = marker_p->coor;
			
			/* locate nearest point */
			PetscTime(&t0_nn);
			ierr = _find_min(pos_p, point_count,patch_point_coords, &nearest_idx);CHKERRQ(ierr);
			PetscTime(&t1_nn);
			time_nn += (t1_nn - t0_nn);

			/* marker index of nearest point */
			marker_index = patch_point_idx[ nearest_idx ];

			/* fetch nearest with index "marker_index" */
			DataFieldAccessPoint(PField,marker_index,(void**)&marker_nearest);
			DataFieldRestoreAccess(PField);

			/* set phase to match the nearest */
			//marker_p->phase = marker_nearest->phase;

			DataBucketCopyPoint(db,marker_index, db,pid_unsorted);

			/* override unique values */
			DataFieldGetAccess(PField);
			DataFieldAccessPoint(PField,pid_unsorted,(void**)&marker_p);
			marker_p->pid     = pid_orig;
			marker_p->coor[0] = xp_orig[0];
			marker_p->coor[1] = xp_orig[1];
			marker_p->coor[2] = xp_orig[2];
			marker_p->wil     = c;
			marker_p->xi[0]   = xip_orig[0];
			marker_p->xi[1]   = xip_orig[1];
			marker_p->xi[2]   = xip_orig[2];
			DataFieldRestoreAccess(PField);
			
			points_assigned++;
			
#if (MPPC_LOG_LEVEL >= 2)
			PetscPrintf(PETSC_COMM_SELF,"[LOG]  point(%d) nearest neighbour(%d) -> phase %d\n",pid_unsorted,marker_index,marker_nearest->phase);
#endif
			
			if (marker_p->phase == MATERIAL_POINT_PHASE_UNASSIGNED) { 
				SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Assigned a region index which is itself un-assigned");
			}
			
		}
		
	} /* end loop on elements */
	
	ierr = PetscFree(patch_point_coords);CHKERRQ(ierr);
	ierr = PetscFree(patch_point_idx);CHKERRQ(ierr);
	
	
#if (MPPC_LOG_LEVEL >= 1)
	printf("[LOG]  points assigned   = %d\n", points_assigned);
	printf("[LOG]  time_nn           = %1.4e (sec)\n", time_nn);
#endif
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "MaterialPointRegionAssignment_v1"
PetscErrorCode MaterialPointRegionAssignment_v1(DataBucket db,DM da)
{
	PetscInt       *pcell_list;
	PSortCtx       *plist;
	int            p,npoints;
	PetscInt       tmp,c,count;
	const PetscInt *elnidx;
	PetscInt       nel,nen;
	DataField      PField;
	PetscLogDouble t0,t1;
	long int       cells_needing_reassignment,cells_needing_reassignment_g;
	long int       points_needing_reassignment;
	PetscInt       *cell_count;
	PetscInt       patch_extend;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	DataBucketGetSizes(db,&npoints,NULL,NULL);

#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG] %s: \n", __FUNCTION__);
#endif	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	
	/* compute number of cells with unassigned region index */
	ierr = PetscMalloc(sizeof(PetscInt)*nel,&cell_count);CHKERRQ(ierr);
	ierr = PetscMemzero(cell_count,sizeof(PetscInt)*nel);CHKERRQ(ierr);
	
	/* count number of points in each cell with phase = MATERIAL_POINT_PHASE_UNASSIGNED */
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));

	for (p=0; p<npoints; p++) {
		MPntStd *marker_p;
		
		DataFieldAccessPoint(PField,p,(void**)&marker_p);
		if (marker_p->phase == MATERIAL_POINT_PHASE_UNASSIGNED) {
			cell_count[ marker_p->wil ]++;
		}
	}

	DataFieldRestoreAccess(PField);

	/* scan number of elements need to be re-assigned */
	points_needing_reassignment = 0;
	cells_needing_reassignment = 0;
	for (c=0; c<nel; c++) {
		points_needing_reassignment += (long int)cell_count[c];
		if (cell_count[c] != 0) {
			cells_needing_reassignment++;	
		}
	}

	/* check if we can exit early */
	ierr = MPI_Allreduce( &cells_needing_reassignment, &cells_needing_reassignment_g, 1, MPI_LONG, MPI_SUM, PetscObjectComm((PetscObject)da) );CHKERRQ(ierr);
	if (cells_needing_reassignment_g == 0) {
#if (MPPC_LOG_LEVEL >= 1)
		PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  !! No region re-assignment equired <global> !!\n");
#endif
		ierr = PetscFree(cell_count);CHKERRQ(ierr);
		PetscFunctionReturn(0);
	}
#if (MPPC_LOG_LEVEL >= 1)
    {
        long int points_needing_reassignment_g;
        PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  !! Region re-assignment required for %D cells <global> !!\n",cells_needing_reassignment_g);
        ierr = MPI_Allreduce( &points_needing_reassignment, &points_needing_reassignment_g, 1, MPI_LONG, MPI_SUM, PetscObjectComm((PetscObject)da) );CHKERRQ(ierr);
        PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  !! Region re-assignment required for %D points <global> !!\n",points_needing_reassignment_g);
    }
#endif

	
	
	/* create sorted list */
	ierr = PetscMalloc(sizeof(PetscInt)*(nel+1),&pcell_list);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PSortCtx)*(npoints),&plist);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db, MPntStd_classname,&PField);
	DataFieldGetAccess(PField);
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p;
		
		DataFieldAccessPoint(PField,p,(void**)&marker_p);
		plist[p].point_index = p;
		plist[p].cell_index  = marker_p->wil;
	}
	DataFieldRestoreAccess(PField);
	
	sort_PSortCx(npoints,plist);
	
	/* sum points per cell */
	ierr = PetscMemzero( pcell_list,sizeof(PetscInt)*(nel+1) );CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		pcell_list[ plist[p].cell_index ]++;
	}
	
	/* create offset list */
	count = 0;
	for (c=0; c<nel; c++) {
		tmp = pcell_list[c];
		pcell_list[c] = count;
		count = count + tmp;
	}
	pcell_list[c] = count;
	
	patch_extend = 1;
	
	PetscTime(&t0);
	ierr = apply_mppc_region_assignment(
														 nel, cell_count, pcell_list,
														 npoints, plist,
														 patch_extend, da,db);CHKERRQ(ierr);
	PetscTime(&t1);
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  time(apply_mppc_region_assignment): %1.4e (sec)\n", t1-t0);
#endif
	
	
	ierr = PetscFree(cell_count);CHKERRQ(ierr);
	ierr = PetscFree(plist);CHKERRQ(ierr);
	ierr = PetscFree(pcell_list);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "apply_mppc_region_assignment_v2"
PetscErrorCode apply_mppc_region_assignment_v2(
																						PetscInt nel, PetscInt cell_count[], PetscInt pcell_list[],
																						PetscInt np, PSortCtx plist[],
																						PetscInt patch_extend,DM da,DataBucket db)
{
	PetscInt        np_per_cell_max,mx,my,mz;
	PetscInt        c,i,j,k,cell_index_i,cell_index_j,cell_index_k,cidx2d,point_count;
	PetscInt        points_per_cell,points_per_patch;
	PetscInt        p;
	DataField       PField;
	double          *patch_point_coords;
	int             *patch_point_idx;
	PetscLogDouble  t0_nn,t1_nn,time_nn = 0.0;
	int             points_assigned = 0;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	
	/* get mx,my from the da */
	ierr = DMDAGetLocalSizeElementQ2(da,&mx,&my,&mz);CHKERRQ(ierr);
	
	/* find max np_per_cell I will need */
	np_per_cell_max = 0;
	for (c=0; c<nel; c++) {
		
		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		if (cell_count[c] == 0) { continue; }
		
		cell_index_k = c / (mx*my);
		cidx2d = c - cell_index_k*(mx*my);
		cell_index_j = cidx2d / mx;
		cell_index_i = cidx2d - cell_index_j * mx;
		
		points_per_patch = 0;
		for ( k=cell_index_k - patch_extend; k<=cell_index_k + patch_extend; k++ ) {
			for ( j=cell_index_j - patch_extend; j<=cell_index_j + patch_extend; j++ ) {
				for ( i=cell_index_i - patch_extend; i<=cell_index_i + patch_extend; i++ ) {
					PetscInt patch_cell_id;
					
					if (i >= mx) { continue; }
					if (j >= my) { continue; }
					if (k >= mz) { continue; }
					if (i < 0) { continue; }
					if (j < 0) { continue; }
					if (k < 0) { continue; }
					
					patch_cell_id = i + j * mx + k * mx*my;
					
					points_per_patch = points_per_patch + (pcell_list[patch_cell_id+1] - pcell_list[patch_cell_id]);
				}
			}
		}
		
		if (points_per_patch > np_per_cell_max) {
			np_per_cell_max = points_per_patch;
		}
	}
	
#if (MPPC_LOG_LEVEL >= 1)
	printf("[LOG]  np_per_patch_max = %d \n", np_per_cell_max );
#endif	
	
	ierr = PetscMalloc(sizeof(double)*3*np_per_cell_max,&patch_point_coords);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(int)*np_per_cell_max,&patch_point_idx);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	
	for (c=0; c<nel; c++) {
		
		/* if cell doesn't contain any points which need re-assignment - skip */
		if (cell_count[c] == 0) { continue; }
		
		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		cell_index_k = c / (mx*my);
		cidx2d = c - cell_index_k*(mx*my);
		cell_index_j = cidx2d / mx;
		cell_index_i = cidx2d - cell_index_j * mx;
		
		/* load points from the patch into a list - only load points with an assigned phase index */
		ierr = PetscMemzero( patch_point_coords, sizeof(double)*3*np_per_cell_max );CHKERRQ(ierr);
		ierr = PetscMemzero( patch_point_idx, sizeof(int)*np_per_cell_max );CHKERRQ(ierr);
		
		point_count = 0;
		
		DataFieldGetAccess(PField);
		for ( k=cell_index_k - patch_extend; k<=cell_index_k + patch_extend; k++ ) {
			for ( j=cell_index_j - patch_extend; j<=cell_index_j + patch_extend; j++ ) {
				for ( i=cell_index_i - patch_extend; i<=cell_index_i + patch_extend; i++ ) {
					PetscInt patch_cell_id;
					
					if (i >= mx) { continue; }
					if (j >= my) { continue; }
					if (k >= mz) { continue; }
					if (i < 0) { continue; }
					if (j < 0) { continue; }
					if (k < 0) { continue; }
					
					patch_cell_id = i + j * mx + k * mx*my;
					points_per_patch = (pcell_list[patch_cell_id+1] - pcell_list[patch_cell_id]);
#if (MPPC_LOG_LEVEL >= 2)
					printf("[LOG]     patch(%d)-(%d,%d,%d) cell(%d)-(%d,%d,%d)  : ppcell = %d \n", c, cell_index_i,cell_index_j,cell_index_k, patch_cell_id,i,j,k,points_per_patch);
#endif					
					for (p=0; p<points_per_patch; p++) {
						MPntStd *marker_p;
						PetscInt pid, pid_unsorted;
						
						pid = pcell_list[patch_cell_id] + p;
						pid_unsorted = plist[pid].point_index;
						
						DataFieldAccessPoint(PField, pid_unsorted ,(void**)&marker_p);
						
						/* skip markers from patch which need to be assigned */
						if (marker_p->phase == MATERIAL_POINT_PHASE_UNASSIGNED) { continue; }
						
						patch_point_coords[3*point_count+0] = marker_p->coor[0];
						patch_point_coords[3*point_count+1] = marker_p->coor[1];
						patch_point_coords[3*point_count+2] = marker_p->coor[2];
						patch_point_idx[point_count]        = pid_unsorted;
#if (MPPC_LOG_LEVEL >= 2)
						printf("[LOG]       patch(%d)/cell(%d) -> p(%d):p->wil,x,y,z = %d %1.4e %1.4e %1.4e \n", c, patch_cell_id, p,marker_p->wil, marker_p->coor[0],marker_p->coor[1],marker_p->coor[2] );
#endif
						point_count++;
					}
					
				}
			}
		}
		DataFieldRestoreAccess(PField);
#if (MPPC_LOG_LEVEL >= 2)
		printf("[LOG]  cell = %d: total points per patch = %d \n", c,point_count);
#endif
		
		/* traverse points in this cell with phase = MATERIAL_POINT_PHASE_UNASSIGNED and find closest point */
		points_per_cell = pcell_list[c+1] - pcell_list[c];
		
		DataFieldGetAccess(PField);

		for (p=0; p<points_per_cell; p++) {
			MPntStd  *marker_p,*marker_nearest_p;
			double   *pos_p;
			int      pid,pid_unsorted,nearest_idx,marker_index;
			
			pid = pcell_list[c] + p;
			pid_unsorted = plist[pid].point_index;
			
			DataFieldAccessPoint(PField,pid_unsorted,(void**)&marker_p);
			
			/* if marker is assigned - skip */
			if (marker_p->phase != MATERIAL_POINT_PHASE_UNASSIGNED) { 
				continue; 
			}
			
			
#if (MPPC_LOG_LEVEL >= 2)
			PetscPrintf(PETSC_COMM_SELF,"[LOG]  cell(%d) point(%d) is un-assigned\n",c,pid_unsorted);
#endif
			
			pos_p = marker_p->coor;
			
			/* locate nearest point */
			PetscTime(&t0_nn);
			ierr = _find_min(pos_p, point_count,patch_point_coords, &nearest_idx);CHKERRQ(ierr);
			PetscTime(&t1_nn);
			time_nn += (t1_nn - t0_nn);
			
			/* marker index of nearest point */
			marker_index = patch_point_idx[ nearest_idx ];
			
			/* fetch nearest with index "marker_index" */
			DataFieldAccessPoint(PField,marker_index,(void**)&marker_nearest_p);
			
			/* set phase to match the nearest */
			marker_p->phase = marker_nearest_p->phase;
			
			//MPntStdSetField_phase_index(marker_p,marker_nearest->phase);
			
			points_assigned++;
			
#if (MPPC_LOG_LEVEL >= 2)
			PetscPrintf(PETSC_COMM_SELF,"[LOG]  point(%d) nearest neighbour(%d) -> phase %d\n",pid_unsorted,marker_index,marker_nearest_p->phase);
#endif
			
			if (marker_p->phase == MATERIAL_POINT_PHASE_UNASSIGNED) { 
				SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Assigned a region index which is itself un-assigned");
			}
			
		}
		
	} /* end loop on elements */

	DataFieldRestoreAccess(PField);
	
	ierr = PetscFree(patch_point_coords);CHKERRQ(ierr);
	ierr = PetscFree(patch_point_idx);CHKERRQ(ierr);
	
	
#if (MPPC_LOG_LEVEL >= 1)
	printf("[LOG]  points assigned   = %d\n", points_assigned);
	printf("[LOG]  time_nn           = %1.4e (sec)\n", time_nn);
#endif
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "MaterialPointRegionAssignment_v2"
PetscErrorCode MaterialPointRegionAssignment_v2(DataBucket db,DM da)
{
	PetscInt       *pcell_list;
	PSortCtx       *plist;
	int            p,npoints;
	PetscInt       tmp,c,count;
	const PetscInt *elnidx;
	PetscInt       nel,nen;
	DataField      PField;
	PetscLogDouble t0,t1;
	long int       cells_needing_reassignment,cells_needing_reassignment_g;
	long int       points_needing_reassignment;
	PetscInt       *cell_count;
	PetscInt       patch_extend;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(db,&npoints,NULL,NULL);
	
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG] %s: \n", __FUNCTION__);
#endif	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	
	/* compute number of cells with unassigned region index */
	ierr = PetscMalloc(sizeof(PetscInt)*nel,&cell_count);CHKERRQ(ierr);
	ierr = PetscMemzero(cell_count,sizeof(PetscInt)*nel);CHKERRQ(ierr);
	
	/* count number of points in each cell with phase = MATERIAL_POINT_PHASE_UNASSIGNED */
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p;
		
		DataFieldAccessPoint(PField,p,(void**)&marker_p);
		if (marker_p->phase == MATERIAL_POINT_PHASE_UNASSIGNED) {
			cell_count[ marker_p->wil ]++;
		}
	}
	
	DataFieldRestoreAccess(PField);
	
	/* scan number of elements need to be re-assigned */
	points_needing_reassignment = 0;
	cells_needing_reassignment = 0;
	for (c=0; c<nel; c++) {
		points_needing_reassignment += (long int)cell_count[c];
		if (cell_count[c] != 0) {
			cells_needing_reassignment++;	
		}
	}
	
	/* check if we can exit early */
	ierr = MPI_Allreduce( &cells_needing_reassignment, &cells_needing_reassignment_g, 1, MPI_LONG, MPI_SUM, PetscObjectComm((PetscObject)da) );CHKERRQ(ierr);
	if (cells_needing_reassignment_g == 0) {
#if (MPPC_LOG_LEVEL >= 1)
		PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  !! No region re-assignment equired <global> !!\n");
#endif
		ierr = PetscFree(cell_count);CHKERRQ(ierr);
		PetscFunctionReturn(0);
	}
#if (MPPC_LOG_LEVEL >= 1)
    {
        long int points_needing_reassignment_g;
        PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  !! Region re-assignment required for %D cells <global> !!\n",cells_needing_reassignment_g);
        ierr = MPI_Allreduce( &points_needing_reassignment, &points_needing_reassignment_g, 1, MPI_LONG, MPI_SUM, PetscObjectComm((PetscObject)da) );CHKERRQ(ierr);
        PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  !! Region re-assignment required for %D points <global> !!\n",points_needing_reassignment_g);
    }
#endif
	
	
	
	/* create sorted list */
	ierr = PetscMalloc(sizeof(PetscInt)*(nel+1),&pcell_list);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PSortCtx)*(npoints),&plist);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db, MPntStd_classname,&PField);
	DataFieldGetAccess(PField);
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p;
		
		DataFieldAccessPoint(PField,p,(void**)&marker_p);
		plist[p].point_index = p;
		plist[p].cell_index  = marker_p->wil;
	}
	DataFieldRestoreAccess(PField);
	
	sort_PSortCx(npoints,plist);
	
	/* sum points per cell */
	ierr = PetscMemzero( pcell_list,sizeof(PetscInt)*(nel+1) );CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		pcell_list[ plist[p].cell_index ]++;
	}
	
	/* create offset list */
	count = 0;
	for (c=0; c<nel; c++) {
		tmp = pcell_list[c];
		pcell_list[c] = count;
		count = count + tmp;
	}
	pcell_list[c] = count;
	
	patch_extend = 1;
	
	PetscTime(&t0);
	ierr = apply_mppc_region_assignment_v2(
																			nel, cell_count, pcell_list,
																			npoints, plist,
																			patch_extend, da,db);CHKERRQ(ierr);
	PetscTime(&t1);
#if (MPPC_LOG_LEVEL >= 1)
	PetscPrintf(PetscObjectComm((PetscObject)da),"[LOG]  time(apply_mppc_region_assignment): %1.4e (sec)\n", t1-t0);
#endif
	
	
	ierr = PetscFree(cell_count);CHKERRQ(ierr);
	ierr = PetscFree(plist);CHKERRQ(ierr);
	ierr = PetscFree(pcell_list);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPPCCreateSortedCtx"
PetscErrorCode MPPCCreateSortedCtx(DataBucket db,DM da,PetscInt *_np,PetscInt *_nc,PSortCtx **_plist,PetscInt **_pcell_list)
{
	PetscInt        *pcell_list;
	PSortCtx        *plist;
	PetscInt        p,npoints;
    int             npoints32;
	PetscInt        tmp,c,count;
	const PetscInt  *elnidx;
	PetscInt        nel,nen;
	DataField       PField;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	ierr = PetscMalloc(sizeof(PetscInt)*(nel+1),&pcell_list);CHKERRQ(ierr);
	
	DataBucketGetSizes(db,&npoints32,NULL,NULL);
    npoints = (PetscInt)npoints32;
	ierr = PetscMalloc(sizeof(PSortCtx)*(npoints),&plist);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p;
		
		DataFieldAccessPoint(PField,p,(void**)&marker_p);
		plist[p].point_index = p;
		plist[p].cell_index  = marker_p->wil;
	}
	DataFieldRestoreAccess(PField);
	
	sort_PSortCx(npoints,plist);
	
	/* sum points per cell */
	ierr = PetscMemzero( pcell_list,sizeof(PetscInt)*(nel+1) );CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		pcell_list[ plist[p].cell_index ]++;
	}
	
	/* create offset list */
	count = 0;
	for (c=0; c<nel; c++) {
		tmp = pcell_list[c];
		pcell_list[c] = count;
		count = count + tmp;
	}
	pcell_list[c] = count;
	
    *_np = (PetscInt)npoints;
    *_nc = (PetscInt)nel;
	*_plist      = plist;
    *_pcell_list = pcell_list;
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPPCDestroySortedCtx"
PetscErrorCode MPPCDestroySortedCtx(DataBucket db,DM da,PSortCtx **_plist,PetscInt **_pcell_list)
{
	PetscInt        *pcell_list;
	PSortCtx        *plist;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
    if (_plist)      {
        plist      = *_plist;
        ierr = PetscFree(plist);CHKERRQ(ierr);
        *_plist      = NULL;
    }
    if (_pcell_list) {
        pcell_list = *_pcell_list;
        ierr = PetscFree(pcell_list);CHKERRQ(ierr);
        *_pcell_list = NULL;
    }
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPPCSortedCtxGetNumberOfPointsPerCell"
PetscErrorCode MPPCSortedCtxGetNumberOfPointsPerCell(DataBucket db,PetscInt cell_idx,PetscInt pcell_list[],PetscInt *np)
{
    PetscInt       points_per_cell;
    
    points_per_cell = pcell_list[cell_idx+1] - pcell_list[cell_idx];
    *np = points_per_cell;
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPPCSortedCtxGetPointByCell"
PetscErrorCode MPPCSortedCtxGetPointByCell(DataBucket db,PetscInt cell_idx,PetscInt pidx,PSortCtx plist[],PetscInt pcell_list[],MPntStd **point)
{
    PetscInt       points_per_cell;
	DataField      PField;
    MPntStd        *mp_std,*marker_p;
    PetscInt       pid,pid_unsorted;
    
    points_per_cell = pcell_list[cell_idx+1] - pcell_list[cell_idx];
    if (pidx >= points_per_cell) {
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Requesting a marker index which is larger than the number of points per cell");
    }
    
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
    mp_std = PField->data;
    
    pid = pcell_list[cell_idx] + pidx;
    pid_unsorted = plist[pid].point_index;

    marker_p = &mp_std[pid_unsorted];
    *point = marker_p;
    
	PetscFunctionReturn(0);
}

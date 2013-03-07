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
 **    Filename:      dmda_iterator.c
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

#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
#include "dmda_iterator.h"

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d"
PetscErrorCode DMDAVecTraverse3d(DM da,Vec X,PetscInt dof_idx,PetscBool (*eval)(PetscScalar*,PetscScalar*,void*),void *ctx)
{
	PetscInt i,j,k,si,sj,sk,m,n,p,M,N,P,ndof;
	DM cda;
	Vec coords;
	DMDACoor3d ***LA_coords;	
	PetscScalar pos[3];
	PetscScalar val;
	PetscBool impose_value;
	PetscScalar ****LA_X;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;	
	ierr = DMDAGetInfo(da,0, &M,&N,&P, 0,0,0, &ndof,0, 0,0,0, 0);CHKERRQ(ierr);
	if (dof_idx >= ndof) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"dof_index >= dm->blocksize"); }
	
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coords);CHKERRQ(ierr);
	if (!coords) { 
		PetscPrintf(PETSC_COMM_WORLD,"WARNING: coordinates not set\n"); 
	} else {
		ierr = DMDAVecGetArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	}
	ierr = DMDAVecGetArrayDOF(da,X,&LA_X);CHKERRQ(ierr);
	
	for (k=sk; k<sk+p; k++) {
		for (j=sj; j<sj+n; j++) {
			for (i=si; i<si+m; i++) {
				
				if (coords) {
					pos[0] = LA_coords[k][j][i].x;
					pos[1] = LA_coords[k][j][i].y;
					pos[2] = LA_coords[k][j][i].z;
				} else {
					pos[0] = pos[1] = pos[2] = 0.0;
				}
				
				impose_value = PETSC_FALSE;
				impose_value = eval(pos,&val,ctx);
				if (impose_value) {
					LA_X[k][j][i][dof_idx] = val;
				}
				
			}
		}
	}
	ierr = DMDAVecGetArrayDOF(da,X,&LA_X);CHKERRQ(ierr);
	if (coords) {
		ierr = DMDAVecRestoreArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverseIJK"
PetscErrorCode DMDAVecTraverseIJK(DM da,Vec X,PetscInt dof_idx,PetscBool (*eval)(PetscScalar*,PetscInt*,PetscInt*,PetscScalar*,void*),void *ctx)
{
	PetscInt i,j,k,si,sj,sk,m,n,p,M,N,P,ndof;
	DM cda;
	Vec coords;
	DMDACoor3d ***LA_coords;	
	PetscScalar pos[3];
	PetscScalar val;
	PetscInt G[3],L[3];
	PetscBool impose_value;
	PetscScalar ****LA_X;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;	
	ierr = DMDAGetInfo(da,0, &M,&N,&P, 0,0,0, &ndof,0, 0,0,0, 0);CHKERRQ(ierr);
	if (dof_idx >= ndof) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"dof_index >= dm->blocksize"); }
	
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coords);CHKERRQ(ierr);
	if (!coords) { 
		PetscPrintf(PETSC_COMM_WORLD,"WARNING: coordinates not set\n"); 
	} else {
		ierr = DMDAVecGetArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	}
	
	ierr = DMDAVecGetArrayDOF(da,X,&LA_X);CHKERRQ(ierr);
	
	for (k=sk; k<sk+p; k++) {
		for (j=sj; j<sj+n; j++) {
			for (i=si; i<si+m; i++) {
				
				if (coords) {
					pos[0] = LA_coords[k][j][i].x;
					pos[1] = LA_coords[k][j][i].y;
					pos[2] = LA_coords[k][j][i].z;
				} else {
					pos[0] = pos[1] = pos[2] = 0.0;
				}
				
				G[0] = i;
				G[1] = j;
				G[2] = k;
				
				L[0] = i-si;
				L[1] = j-sj;
				L[2] = k-sk;
				
				impose_value = PETSC_FALSE;
				impose_value = eval(pos,G,L,&val,ctx);
				if (impose_value) {
					LA_X[k][j][i][dof_idx] = val;
				}
				
			}
		}
	}
	ierr = DMDAVecGetArrayDOF(da,X,&LA_X);CHKERRQ(ierr);
	if (coords) {
		ierr = DMDAVecRestoreArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	}	
	PetscFunctionReturn(0);
}


/* HELPERS */
/* constant evaluator */
#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_Constant"
PetscBool DMDAVecTraverse3d_Constant(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	PetscScalar c = *( (PetscScalar*)ctx );
	PetscBool impose;
	*val = c;
	impose = PETSC_TRUE;
	return impose;
}
#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverseIJK_Constant"
PetscBool DMDAVecTraverseIJK_Constant(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx)
{
	PetscScalar c = *( (PetscScalar*)ctx );
	PetscBool impose;
	*val = c;
	impose = PETSC_TRUE;
	return impose;
}

/* linearly interpolate evaluator */
#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_InterpCtxSetUp_X"
PetscErrorCode DMDAVecTraverse3d_InterpCtxSetUp_X(DMDAVecTraverse3d_InterpCtx *c,PetscScalar a, PetscScalar b,PetscScalar ox)
{
	PetscFunctionBegin;	
	
	c->A[0] = a;
	c->B[0] = b;
	c->Ox[0] = ox;
	
	c->Ox[1] = c->Ox[2] = 0.0;
	c->A[1] = c->A[2] = 0.0;
	c->B[1] = c->B[2] = 0.0;
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_InterpCtxSetUp_Y"
PetscErrorCode DMDAVecTraverse3d_InterpCtxSetUp_Y(DMDAVecTraverse3d_InterpCtx *c,PetscScalar a, PetscScalar b,PetscScalar oy)
{
	PetscFunctionBegin;	
	
	c->A[1] = a;
	c->B[1] = b;
	c->Ox[1] = oy;
	
	c->Ox[0] = c->Ox[2] = 0.0;
	c->A[0] = c->A[2] = 0.0;
	c->B[0] = c->B[2] = 0.0;
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_InterpCtxSetUp_Z"
PetscErrorCode DMDAVecTraverse3d_InterpCtxSetUp_Z(DMDAVecTraverse3d_InterpCtx *c,PetscScalar a, PetscScalar b,PetscScalar oz)
{
	PetscFunctionBegin;	
	
	c->A[2] = a;
	c->B[2] = b;
	c->Ox[2] = oz;
	
	c->Ox[0] = c->Ox[1] = 0.0;
	c->A[0] = c->A[1] = 0.0;
	c->B[0] = c->B[1] = 0.0;
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_InterpCtxSetUp"
PetscErrorCode DMDAVecTraverse3d_InterpCtxSetUp(DMDAVecTraverse3d_InterpCtx *c,PetscInt dir,PetscScalar a, PetscScalar b,PetscScalar ox)
{
	PetscFunctionBegin;	
	
	if (dir < 0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"0 <= dir <=2");
	}
	if (dir > 2) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"0 <= dir <=2");
	}
	
	c->A[dir] = a;
	c->B[dir] = b;
	c->Ox[dir] = ox;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_Interp"
PetscBool DMDAVecTraverse3d_Interp(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	DMDAVecTraverse3d_InterpCtx *c = (DMDAVecTraverse3d_InterpCtx*)ctx;
	PetscBool impose;
	
	*val = (c->A[0]*(pos[0]-c->Ox[0]) + c->B[0]) + (c->A[1]*(pos[1]-c->Ox[1]) + c->B[1]) + (c->A[2]*(pos[2]-c->Ox[2]) + c->B[2]);
	
	impose = PETSC_TRUE;
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverseIJK_Interp"
PetscBool DMDAVecTraverseIJK_Interp(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx)
{
	DMDAVecTraverse3d_InterpCtx *c = (DMDAVecTraverse3d_InterpCtx*)ctx;
	PetscBool impose;
	
	*val = (c->A[0]*(pos[0]-c->Ox[0]) + c->B[0]) + (c->A[1]*(pos[1]-c->Ox[1]) + c->B[1]) + (c->A[2]*(pos[2]-c->Ox[2]) + c->B[2]);
	
	impose = PETSC_TRUE;
	return impose;
}

/* hydrostatic pressure evaluator */
#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverseIJK_HydroStaticPressure_v1"
PetscBool DMDAVecTraverseIJK_HydroStaticPressure_v1(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx)
{
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx *c = (DMDAVecTraverse3d_HydrostaticPressureCalcCtx*)ctx;
	PetscScalar z,P;
	PetscBool impose;
	PetscErrorCode ierr;
	
	z = c->ref_height - pos[1];
	P = c->rho * c->grav * z;
	
	*val = P;
	impose = PETSC_TRUE;
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverseIJK_HydroStaticPressure_v2"
PetscBool DMDAVecTraverseIJK_HydroStaticPressure_v2(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx)
{
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx *c = (DMDAVecTraverse3d_HydrostaticPressureCalcCtx*)ctx;
	PetscScalar z,P;
	PetscReal dz;
	PetscBool impose;
	PetscErrorCode ierr;
	
	dz = c->ref_height / ((PetscReal)( c->ref_N+1 ) );
	z = (PetscScalar)(c->ref_N - global_index[1]);
	z = z / ((PetscScalar)( c->ref_N ));
	z = z * c->ref_height + 0.5*dz;
	//printf("z = %1.4e \n",z);
	P = c->surface_pressure + c->rho * c->grav * (z);
	
	*val = P;
	impose = PETSC_TRUE;
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverseIJK_HydroStaticPressure_dpdy_v2"
PetscBool DMDAVecTraverseIJK_HydroStaticPressure_dpdy_v2(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx)
{
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx *c = (DMDAVecTraverse3d_HydrostaticPressureCalcCtx*)ctx;
	PetscScalar dPdy;
	PetscBool impose;
	PetscErrorCode ierr;
	
	dPdy = -c->rho * c->grav;
	
	*val = 0.5*dPdy;
	impose = PETSC_TRUE;
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverseIJK_HydroStaticPressure"
PetscBool DMDAVecTraverseIJK_HydroStaticPressure(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx)
{
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx *c = (DMDAVecTraverse3d_HydrostaticPressureCalcCtx*)ctx;
	PetscScalar z,P;
	PetscScalar ****LA_surface_coords;
	PetscBool impose;
	DM cda;
	Vec surface_coords;
	PetscErrorCode ierr;
	
	ierr = DMDAGetCoordinateDA(c->surface_da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(c->surface_da,&surface_coords);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(cda,surface_coords,&LA_surface_coords);CHKERRQ(ierr);
	
	z = LA_surface_coords[local_index[2]][0][local_index[0]][1] - pos[1];
	P = c->rho * c->grav * z;
	
	ierr = DMDAVecRestoreArrayDOF(cda,surface_coords,&LA_surface_coords);CHKERRQ(ierr);

	*val = P;
	impose = PETSC_TRUE;
	return impose;
}


#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_GaussianXY"
PetscBool DMDAVecTraverse3d_GaussianXY(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	PetscBool impose;
	PetscScalar x,y,z;
	
	impose = PETSC_TRUE;
	x = pos[0];
	y = pos[1];
	z = pos[2];
	*val = exp( -100.0*( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) )); 

	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_GaussianXYZ"
PetscBool DMDAVecTraverse3d_GaussianXYZ(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	PetscBool impose;
	PetscScalar x,y,z;
	
	impose = PETSC_TRUE;
	x = pos[0];
	y = pos[1];
	z = pos[2];
	*val = exp( -100.0*( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5) )); 
	
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_StepX"
PetscBool DMDAVecTraverse3d_StepX(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	PetscBool impose;
	PetscScalar x,y,z;
	
	impose = PETSC_TRUE;
	x = pos[0];
	y = pos[1];
	z = pos[2];
	if (x < 0.5) {
		*val = 2.0;
	} else {
		*val = 1.0;
	}
	
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_StepXY"
PetscBool DMDAVecTraverse3d_StepXY(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	PetscBool impose;
	PetscScalar x,y,z;
	
	impose = PETSC_TRUE;
	x = pos[0];
	y = pos[1];
	z = pos[2];
	if ( (x < 0.5) && (y < 0.5) ) {
		*val = 2.0;
	} else {
		*val = 1.0;
	}
	
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_StepXYZ"
PetscBool DMDAVecTraverse3d_StepXYZ(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	PetscBool impose;
	PetscScalar x,y,z;
	
	impose = PETSC_TRUE;
	x = pos[0];
	y = pos[1];
	z = pos[2];
	if ( (x < 0.5) && (y < 0.5) && (z < 0.5) ) {
		*val = 2.0;
	} else {
		*val = 1.0;
	}
	
	return impose;
}


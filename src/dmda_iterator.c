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
#include "dmda_update_coords.h"
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
#define __FUNCT__ "DMDAVecTraverse3d_StepWithDirection"
PetscBool DMDAVecTraverse3d_StepWithDirection(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	int dir = *((int*)ctx);
	PetscBool impose;
	PetscScalar x,y,z;
	
	impose = PETSC_TRUE;
	x = pos[0];
	y = pos[1];
	z = pos[2];
	
	switch (dir) {
		case 0:
			if (x < 0.5) {
				*val = 2.0;
			} else {
				*val = 1.0;
			}
			break;

		case 1:
			if (y < 0.5) {
				*val = 2.0;
			} else {
				*val = 1.0;
			}
			break;

		case 2:
			if (z < 0.5) {
				*val = 2.0;
			} else {
				*val = 1.0;
			}
			break;
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

/*
 Example usage:
 
 PetscReal vals[4];
 
 // will create gradients T = -933.0 * y //
 vals[0] = 0.0;    // offset
 vals[1] = 0.0;    // x grad
 vals[2] = -933.0; // y grad
 vals[3] = 0.0;    // z grad
 
 ierr = DMDAVecTraverse3d(daT,temperature,0, DMDAVecTraverse3d_LinearFunctionXYZ, (void*)vals);CHKERRQ(ierr);
 
 */
#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_LinearFunctionXYZ"
PetscBool DMDAVecTraverse3d_LinearFunctionXYZ(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	PetscScalar x,y,z;
	PetscReal  *coeffA;
	PetscBool  impose;
	
	/* get coordinates */
	x = pos[0];
	y = pos[1];
	z = pos[2];

	/* fetch user data */
	coeffA = (PetscReal*)ctx;

	/* define value you want to set into the vector */
	*val = coeffA[0] + coeffA[1]*x + coeffA[2]*y + coeffA[3]*z; 
	
	/* indicate you want to set this value into the vector */
	impose = PETSC_TRUE;

	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_ERFC3DFunctionXYZ"
PetscBool DMDAVecTraverse3d_ERFC3DFunctionXYZ(PetscScalar pos[],PetscScalar *val,void *ctx)
{
	PetscScalar x,y,z;
	PetscReal  *coeffs;
	PetscBool  impose;
	PetscReal  xc,zc,age0,age_anom,L_bar,Tbot,wx,wz,age;
	/* get coordinates */
	x = pos[0];
	y = pos[1];
	z = pos[2];
	/* fetch user data */
	coeffs = (PetscReal*)ctx;
    xc       = coeffs[0]; 
    zc       = coeffs[1]; 
    age0     = coeffs[2];
    age_anom = coeffs[3];
    L_bar    = coeffs[4];
    Tbot     = coeffs[5];
    wx       = coeffs[6];
    wz       = coeffs[7];
    
	/* define value you want to set into the vector */ 
    age = (1-exp(-pow((x-xc)*wx,2))*exp(-pow((z-zc)*wz,2))*(age0-age_anom)/age0)*age0;
    *val = -Tbot*erfc((-y*L_bar)/sqrt(1e-6*age*3.14e13))+Tbot;
	/* indicate you want to set this value into the vector */
	impose = PETSC_TRUE;
    
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverseIJK_ZeroInteriorMinusNmax"
PetscBool DMDAVecTraverseIJK_ZeroInteriorMinusNmax(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx)
{
	PetscInt  upper_limit = *((PetscInt*)ctx);
	PetscBool impose;
	
	*val = 0.0;
	impose = PETSC_TRUE;
	if (global_index[1] == (upper_limit-1) ) {
		impose = PETSC_FALSE;
	}
	
	return impose;
}

/*
 PetscInt plane: The plane (I,J or K) we wish to sweep over. Valid values are {I=0, J=1, K=1}
								 If plane = 0, then we sweep of J,K, taking the i-index value from the variable "index"
 PetscInt index: The I,J or K value which defines the plane in the mesh we will sweep over.
 PetscInt coord_dof: The degree of freedom we wish to assign in the coordinate vector. Valid values are {0,1,2}.
 
 The function pointer "eval" is defined the following way;
   PetscBool eval(PetscScalar position[],PetscInt global_ijk[],PetscInt local_ijk[],PetscScalar *val,void *ctx);
 
 Input: 
   PetscScalar position[]: the coordinate vector for node i,j,k
   PetscInt global_ijk[]: the current value of i,j,k in the mesh in the global (natural) numbering
   PetscInt local_ijk[]: the current value of i,j,k in the mesh in the local numbering
   void *ctx: pointer to used defined context
 Output:
   PetscScalar *val: the value of the new coodinate
   PetscBool flag: flag indicates whether the output "val" should be set into the coordinate vector for the DMDA   
*/
#undef __FUNCT__
#define __FUNCT__ "DMDACoordTraverseIJK"
PetscErrorCode DMDACoordTraverseIJK(DM da,PetscInt plane,PetscInt index,PetscInt coord_dof,PetscBool (*eval)(PetscScalar*,PetscInt*,PetscInt*,PetscScalar*,void*),void *ctx)
{
	PetscInt       i,j,k,si,sj,sk,m,n,p;
	DM             cda;
	Vec            coords;
	DMDACoor3d     ***LA_coords;	
	PetscScalar    pos[3];
	PetscInt       L[3],G[3];
	PetscScalar    val;
	PetscBool      impose_value;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;	

	if (plane >= 3) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"0 <= plane <= 2"); }
	if (coord_dof >= 3) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"0 <= coord_dof <= 2"); }

	ierr = DMDAGetCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coords);CHKERRQ(ierr);
	if (!coords) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"coordinates not set"); }
	
	ierr = DMDAVecGetArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	
	for (k=sk; k<sk+p; k++) {
		if ((plane == 2) && (k != index)) { continue; }

		for (j=sj; j<sj+n; j++) {
			if ((plane == 1) && (j != index)) { continue; }

			for (i=si; i<si+m; i++) {
				if ((plane == 0) && (i != index)) { continue; }
				
				pos[0] = LA_coords[k][j][i].x;
				pos[1] = LA_coords[k][j][i].y;
				pos[2] = LA_coords[k][j][i].z;
				
				G[0] = i;
				G[1] = j;
				G[2] = k;
				
				L[0] = i-si;
				L[1] = j-sj;
				L[2] = k-sk;
				
				impose_value = PETSC_FALSE;
				impose_value = eval(pos,G,L,&val,ctx);
				if (impose_value) {
					pos[ coord_dof ] = val;

					LA_coords[k][j][i].x = pos[0];
					LA_coords[k][j][i].y = pos[1];
					LA_coords[k][j][i].z = pos[2];
				}
				
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

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
 **    Filename:      ptatin_std_dirichlet_boundary_conditions.c
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
/*

 Collection of standard boundary conditions for hex meshes.

 Definitions:
 Regarding directions, we will adopt the following nomenclecture
 NX,NY,NZ are the number of nodes in x,y,z direction respectively

 NorthFace: Implies boundary defined by all nodes with ny = NY-1
 
 EastFace:  Implies boundary defined by all nodes with nx = NX-1
 
 SouthFace: Implies boundary defined by all nodes with ny = 0
 
 WestFace:  Implies boundary defined by all nodes with nx = 0

 FrontFace: Implies boundary defined by all nodes with nz = NZ-1
 
 BackFace:  Implies boundary defined by all nodes with nz = 0
 
*/

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_bcs.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"


#undef __FUNCT__
#define __FUNCT__ "DirichletBC_FreeSlip"
PetscErrorCode DirichletBC_FreeSlip(BCList bclist,DM dav,BoundaryFaceType face)
{
	PetscScalar    zero = 0.0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	switch (face) {

		case EAST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case WEST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case NORTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;
			
		case SOUTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case FRONT_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case BACK_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Face must N,E,S,W,F,B");
			break;
	}
		
	PetscFunctionReturn(0);	
}

#undef __FUNCT__
#define __FUNCT__ "DirichletBC_SetConstant"
PetscErrorCode DirichletBC_SetConstant(BCList bclist,DM dav,BoundaryFaceType face,PetscInt dof,PetscScalar value)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
    
	switch (face) {
            
		case EAST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,dof,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
            
		case WEST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,dof,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
            
		case NORTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,dof,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case SOUTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,dof,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
            
		case FRONT_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,dof,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
            
		case BACK_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,dof,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
            
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Face must N,E,S,W,F,B");
			break;
	}
    
	PetscFunctionReturn(0);	
}

/* 
 Note: This is not general and will not work on deformed boundaries.
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyNormalVelocity"
PetscErrorCode DirichletBC_ApplyNormalVelocity(BCList bclist,DM dav,BoundaryFaceType face,PetscReal v_normal)
{
	PetscScalar    value = v_normal;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	switch (face) {
			
		case EAST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case WEST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case NORTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case SOUTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case FRONT_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case BACK_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;

		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Face must N,E,S,W,F,B");
			break;
	}
	
	PetscFunctionReturn(0);	
}

/* 
 Extend mesh in x direction with a specified strain rate, Exx.
 The applied velocity has the same magnitude in east and west directions.

 Compute size of mesh (Lx) in x direction.
 Exx = du/dx = (vx_east - vx_west)/Lx = 2 Vx_bc / Lx
 Vx_bc = 0.5 . Exx . Lx 
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyStrainRateExx"
PetscErrorCode DirichletBC_ApplyStrainRateExx(BCList bclist,DM dav,PetscReal exx_bc)
{
	PetscScalar value;
	PetscReal MeshMin[3],MeshMax[3];
	PetscReal origin[3];
	PetscReal Lx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	origin[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	origin[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	origin[2] = 0.5*(MeshMax[2] + MeshMin[2]);

	Lx = (MeshMax[0] - MeshMin[0]);
	value = 0.5 * exx_bc * (Lx);

	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);

	value = -value;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);

	PetscFunctionReturn(0);	
}

/* 
 Apply tangential velocity along the north/south faces via a specified shear strain rate, Exy.
 Normal velocity on north/south face is set to be zero
 The applied velocity has the same magnitude in east and west directions.
 
 Compute size of mesh (Ly) in y direction.
 Exy = 0.5.(du/dy + dv/dx)= 0.5.(vx_east - vx_west)/Lx = Vx_bc / Ly
 Vx_bc = Exy . Ly
*/
/*
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyStrainRateExy"
PetscErrorCode DirichletBC_ApplyStrainRateExy(BCList bclist,DM dav,PetscReal exy_bc)
{
	PetscScalar value;
	PetscReal MeshMin[3],MeshMax[3];
	PetscReal origin[3];
	PetscReal Ly;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	origin[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	origin[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	origin[2] = 0.5*(MeshMax[2] + MeshMin[2]);

	Ly = (MeshMax[1] - MeshMin[1]);

	value = exy_bc * (Ly);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	value = -value;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	value = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}
*/

/* 
 Extend/compress in x direction, compress/extend in y direction.

 v_shear(x,y) = value.(x.i - y.j)  
 
 Apply x component of v_shear() along east-west boundary faces
 Apply y component of v_shear() along north-south boundary faces 
*/
/*
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningY"
PetscErrorCode DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningY(BCList list,DM dav,PetscReal vx_bc)
{
	PetscScalar v_normal;
	PetscReal MeshMin[3],MeshMax[3];
	PetscReal origin[3];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	origin[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	origin[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	origin[2] = 0.5*(MeshMax[2] + MeshMin[2]);
		
	v_normal = vx_bc * (MeshMax[0]-origin[0]);
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,EAST_FACE,v_normal);CHKERRQ(ierr);
	v_normal = vx_bc * (MeshMin[0]-origin[0]);
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,WEST_FACE,v_normal);CHKERRQ(ierr);

	v_normal = -vx_bc * (MeshMax[1]-origin[1]);
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,NORTH_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc * (MeshMin[1]-origin[1]);
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,SOUTH_FACE,v_normal);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}
*/
/* 
 Extend/compress in x direction, compress/extend in z direction.
 
 v_shear(x,z) = value.(x.i - z.k)  
 
 Apply x component of v_shear() along east-west boundary faces
 Apply z component of v_shear() along front-back boundary faces 

 Actually, i prefer to implement it like this
 2.vx/Lx = 2.vz/Lz
 
 vz = Lz.vx/Lx
 
*/
/*
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ"
PetscErrorCode DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ(BCList list,DM dav,PetscReal vx_bc)
{
	PetscScalar v_normal;
	PetscReal MeshMin[3],MeshMax[3];
	PetscReal origin[3],Lx,Lz;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	origin[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	origin[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	origin[2] = 0.5*(MeshMax[2] + MeshMin[2]);	
	
	Lx = (MeshMax[0] - MeshMin[0]);
	Lz = (MeshMax[2] - MeshMin[2]);

	v_normal =  vx_bc ;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,EAST_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc ;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,WEST_FACE,v_normal);CHKERRQ(ierr);
	
	v_normal =  vx_bc * Lz/Lx;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,FRONT_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc * Lz/Lx;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,BACK_FACE,v_normal);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}
*/

/* 
 Extend/compress in x direction, compress/extend in y,z direction.
 
 v_shear(x,y,z) = value.(x.i - 0.5.y.j - 0.5.z.k)  
 
 Apply x component of v_shear() along east-west boundary faces
 Apply y component of v_shear() along north-south boundary faces 
 Apply z component of v_shear() along front-back boundary faces 
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantVolumeDomain_ExtensionX"
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionX(BCList list,DM dav,PetscReal vx_bc)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(list,dav,0.5,vx_bc);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

/* 
 Extend/compress in x direction, compress/extend in y,z direction.
 
 v_shear(x,y,z) = value.(x.i - beta.y.j - (1-beta).z.k)  
 
 Apply x component of v_shear() along east-west boundary faces
 Apply y component of v_shear() along north-south boundary faces 
 Apply z component of v_shear() along front-back boundary faces
 
 du/dx + dv/dy + dw/dz = 0
 
 (u_1 - u_2)/Lx + beta (v_1 - v_2)/Ly + (1-beta) (v_1 - v_2)/Lz = 0

 beta (v_1 - v_2)/Ly + (1-beta) (v_1 - v_2)/Lz = -2.0.u_1/Lx
 
*/
/*
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening"
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(BCList list,DM dav,PetscReal beta,PetscReal vx_bc)
{
	PetscScalar v_normal;
	PetscReal MeshMin[3],MeshMax[3];
	PetscReal origin[3];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	origin[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	origin[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	origin[2] = 0.5*(MeshMax[2] + MeshMin[2]);
	
	v_normal = vx_bc;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,EAST_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,WEST_FACE,v_normal);CHKERRQ(ierr);
	
	Exx = 2.0 * vx_bc / Lx;
	v_normal = -0.5 * Exx * beta;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,NORTH_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -v_normal;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,SOUTH_FACE,v_normal);CHKERRQ(ierr);
	
	v_normal = -vx_bc * (1.0-beta) * (MeshMax[2]-origin[2]);
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,FRONT_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc * (1.0-beta) * (MeshMin[2]-origin[2]);
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,BACK_FACE,v_normal);CHKERRQ(ierr);

	PetscFunctionReturn(0);	
}
*/

typedef struct {
	PetscInt  dof_idx;
	PetscReal alpha,beta,gamma;
	PetscReal Ox[3];
} StrainRateBCCtx;

#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_StrainRate"
PetscBool BCListEvaluator_StrainRate(PetscScalar position[],PetscScalar *value,void *ctx) 
{
	PetscReal       V[3];
	PetscBool       impose_dirichlet = PETSC_TRUE;
	StrainRateBCCtx *strainrrate_ctx;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	strainrrate_ctx = (StrainRateBCCtx*)ctx;
	
	V[0] = strainrrate_ctx->alpha * ( position[0] - strainrrate_ctx->Ox[0] );
	V[1] = strainrrate_ctx->beta  * ( position[1] - strainrrate_ctx->Ox[1] );
	V[2] = strainrrate_ctx->gamma * ( position[2] - strainrrate_ctx->Ox[2] );
	
	switch (strainrrate_ctx->dof_idx) {
		case 0:
			*value = V[0];
			break;
		case 1:
			*value = V[1];
			break;
		case 2:
			*value = V[2];
			break;

		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"dof_idx must be 0,1,2");
			break;
	}
	return impose_dirichlet;
}

#undef __FUNCT__
#define __FUNCT__ "DirichletBC_DefineExx"
PetscReal DirichletBC_DefineExx(StrainRateBCCtx *ctx,PetscReal Exx,DM dav)
{
	PetscReal MeshMin[3],MeshMax[3];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ctx->alpha = Exx;
	ctx->beta  = 0.0;
	ctx->gamma = 0.0;

	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	ctx->Ox[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	ctx->Ox[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	ctx->Ox[2] = 0.5*(MeshMax[2] + MeshMin[2]);	

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DirichletBC_DefineEyy"
PetscReal DirichletBC_DefineEyy(StrainRateBCCtx *ctx,PetscReal Eyy,DM dav)
{
	PetscReal MeshMin[3],MeshMax[3];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ctx->alpha = 0.0;
	ctx->beta  = Eyy;
	ctx->gamma = 0.0;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	ctx->Ox[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	ctx->Ox[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	ctx->Ox[2] = 0.5*(MeshMax[2] + MeshMin[2]);	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DirichletBC_DefineEzz"
PetscReal DirichletBC_DefineEzz(StrainRateBCCtx *ctx,PetscReal Ezz,DM dav)
{
	PetscReal MeshMin[3],MeshMax[3];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ctx->alpha = 0.0;
	ctx->beta  = 0.0;
	ctx->gamma = Ezz;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	ctx->Ox[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	ctx->Ox[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	ctx->Ox[2] = 0.5*(MeshMax[2] + MeshMin[2]);	
	
	PetscFunctionReturn(0);
}

/* 
 Extend mesh in x direction with a specified strain rate, Exx.
 The applied velocity has the same magnitude in east and west directions.
 
 Compute size of mesh (Lx) in x direction.
 Exx = du/dx = (vx_east - vx_west)/Lx = 2 Vx_bc / Lx
 Vx_bc = 0.5 . Exx . Lx 
 */
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyDirectStrainRate"
PetscErrorCode DirichletBC_ApplyDirectStrainRate(BCList bclist,DM dav,PetscReal Evalue,PetscInt direction)
{
	StrainRateBCCtx ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	

	ctx.dof_idx = direction;
	switch (direction) {
		case 0:
			DirichletBC_DefineExx(&ctx,Evalue,dav);
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
			break;

		case 1:
			DirichletBC_DefineEyy(&ctx,Evalue,dav);
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
			break;

		case 2:
			DirichletBC_DefineEzz(&ctx,Evalue,dav);
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
			break;
			
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"direction must be 0,1,2");
			break;

	}
	
	
	PetscFunctionReturn(0);	
}

/* 
 Apply tangential velocity (vx=0,vz) along the east/west faces via a specified shear strain rate, Exz.
 Normal velocity on east/west is set to be zero
 Along front/back faces we prescribe a normal velocity in (vx=0,vz) which is consistent with the sense of shear.
 x velocity compnent on front/bac is set to be zero
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyStrainRateExz"
PetscErrorCode DirichletBC_ApplyStrainRateExz(BCList bclist,DM dav,PetscReal exz_bc)
{
	PetscReal       MeshMin[3],MeshMax[3];
	PetscScalar     value;
	StrainRateBCCtx ctx;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	ctx.Ox[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	ctx.Ox[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	ctx.Ox[2] = 0.5*(MeshMax[2] + MeshMin[2]);
	
	value = exz_bc * (MeshMax[0] - MeshMin[0]);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	value = -value;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	value = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);

	
	ctx.dof_idx = 0;
	value       = exz_bc * (MeshMax[0] - MeshMin[0]);
	ctx.alpha   = 2.0 * value / (MeshMax[0] - MeshMin[0]);
	ctx.beta    = 0.0;
	ctx.gamma   = 0.0;	
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	
	value = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

/* 
 Apply tangential velocity (vx,vz=0) along the front/back faces via a specified shear strain rate, Exz.
 Normal velocity on front/back is set to be zero
 Along east/west faces we prescribe a normal velocity in (vx,vz=0) which is consistent with the sense of shear.
 z velocity compnent on east/west is set to be zero
 */
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyStrainRateExz_b"
PetscErrorCode DirichletBC_ApplyStrainRateExz_b(BCList bclist,DM dav,PetscReal exz_bc)
{
	PetscReal       MeshMin[3],MeshMax[3];
	PetscScalar     value;
	StrainRateBCCtx ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	ctx.Ox[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	ctx.Ox[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	ctx.Ox[2] = 0.5*(MeshMax[2] + MeshMin[2]);
	
	value = exz_bc * (MeshMax[0] - MeshMin[0]);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	value = -value;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	value = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	ctx.dof_idx = 2;
	value       = exz_bc * (MeshMax[0] - MeshMin[0]);
	ctx.alpha   = 0.0;
	ctx.beta    = 0.0;
	ctx.gamma   = 2.0 * value / (MeshMax[0] - MeshMin[0]);	
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	
	value = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

/* 
 Along east/west faces we prescribe a normal velocity in (vx,vz=0) which is consistent with the sense of shear.
 z velocity compnent on east/west is set to be zero
 */
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyStrainRateExz_c"
PetscErrorCode DirichletBC_ApplyStrainRateExz_c(BCList bclist,DM dav,PetscReal exz_bc)
{
	PetscReal       MeshMin[3],MeshMax[3];
	PetscScalar     value;
	StrainRateBCCtx ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	ctx.Ox[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	ctx.Ox[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	ctx.Ox[2] = 0.5*(MeshMax[2] + MeshMin[2]);
	
	ctx.dof_idx = 2;
	value       = exz_bc * (MeshMax[0] - MeshMin[0]);
	ctx.alpha   = 0.0;
	ctx.beta    = 0.0;
	ctx.gamma   = 2.0 * value / (MeshMax[0] - MeshMin[0]);	
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	
	value = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ"
PetscErrorCode DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ(BCList list,DM dav,PetscReal vx_bc)
{
	PetscReal       MeshMin[3],MeshMax[3];
	PetscScalar     v_normal;
	StrainRateBCCtx ctx;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	ctx.Ox[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	ctx.Ox[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	ctx.Ox[2] = 0.5*(MeshMax[2] + MeshMin[2]);
	
	v_normal =  vx_bc;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,EAST_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,WEST_FACE,v_normal);CHKERRQ(ierr);

	
	ctx.dof_idx = 2;
	ctx.alpha   = 0.0;
	ctx.beta    = 0.0;
	ctx.gamma   = -2.0 * vx_bc/(MeshMax[0] - MeshMin[0]);
	
	ierr = DMDABCListTraverse3d(list,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(list,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);	
}

#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening"
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(BCList list,DM dav,PetscReal factor,PetscReal vx_bc)
{
	PetscScalar     v_normal;
	StrainRateBCCtx ctx;
	PetscReal       Exx,MeshMin[3],MeshMax[3];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	ctx.Ox[0] = 0.5*(MeshMax[0] + MeshMin[0]);
	ctx.Ox[1] = 0.5*(MeshMax[1] + MeshMin[1]);
	ctx.Ox[2] = 0.5*(MeshMax[2] + MeshMin[2]);
	
	if (factor<0.0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"factor must be >= 0.0");
	}
	if (factor>1.0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"factor must be <= 1.0");
	}
	
	v_normal =  vx_bc;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,EAST_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc;
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,WEST_FACE,v_normal);CHKERRQ(ierr);
	
	Exx = 2.0 * vx_bc / (MeshMax[0]-MeshMin[0]);
	ctx.alpha   = 0.0;
	ctx.beta    = -factor * Exx;
	ctx.gamma   = -(1.0-factor)*Exx;
	
	ctx.dof_idx = 1;
	ierr = DMDABCListTraverse3d(list,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(list,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	
	ctx.dof_idx = 2;
	ierr = DMDABCListTraverse3d(list,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(list,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_StrainRate,(void*)&ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

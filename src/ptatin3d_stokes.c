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
 **    Filename:      ptatin3d_stokes.c
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
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_element_q2p1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"

#include "QPntVolCoefStokes_def.h"


#undef __FUNCT__
#define __FUNCT__ "StokesVelocity_GetElementLocalIndices"
PetscErrorCode StokesVelocity_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		el_localIndices[3*n  ] = 3*elnid[n]  ;
		el_localIndices[3*n+1] = 3*elnid[n]+1;
		el_localIndices[3*n+2] = 3*elnid[n]+2;
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "StokesPressure_GetElementLocalIndices"
PetscErrorCode StokesPressure_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<P_BASIS_FUNCTIONS; n++) {
		el_localIndices[n] = elnid[n];
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "StokesVelocityScalar_GetElementLocalIndices"
PetscErrorCode StokesVelocityScalar_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		el_localIndices[n] = elnid[n] ;
	}
	PetscFunctionReturn(0);
}


/* physics component loader */
#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreate_Stokes"
PetscErrorCode PhysCompCreate_Stokes(PhysCompStokes *ctx)
{
	PetscErrorCode ierr;
	PhysCompStokes stokes;
	
	PetscFunctionBegin;
	ierr = PetscMalloc(sizeof(struct _p_PhysCompStokes),&stokes);CHKERRQ(ierr);
	*ctx = stokes;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompDestroy_Stokes"
PetscErrorCode PhysCompDestroy_Stokes(PhysCompStokes *ctx)
{
	PetscErrorCode ierr;
	PhysCompStokes user;
	PetscInt e;
	
	PetscFunctionBegin;
	
	if (!ctx) {PetscFunctionReturn(0);}
	user = *ctx;
	
//	for (e=0; e<QUAD_EDGES; e++) {
//		if (user->surfQ[e]) { ierr = SurfaceQuadratureStokesDestroy(&user->surfQ[e]);CHKERRQ(ierr); }
//	}
//	if (user->Q) { ierr = QuadratureStokesDestroy(&user->Q);CHKERRQ(ierr); }
	if (user->p_bclist) { ierr = BCListDestroy(&user->p_bclist);CHKERRQ(ierr); }
	if (user->u_bclist) { ierr = BCListDestroy(&user->u_bclist);CHKERRQ(ierr); }
  if (user->stokes_pack) { ierr = DMDestroy(&user->stokes_pack);CHKERRQ(ierr); }
	if (user->dap) { ierr = DMDestroy(&user->dap);CHKERRQ(ierr); }
  if (user->dav) { ierr = DMDestroy(&user->dav);CHKERRQ(ierr); }
	if (user) { ierr = PetscFree(user);CHKERRQ(ierr); }
	
	*ctx = PETSC_NULL;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateMesh_Stokes3d"
PetscErrorCode PhysCompCreateMesh_Stokes3d(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx)
{
	DM dav,dap,multipys_pack;
	PetscInt vbasis_dofs;
	PetscInt pbasis_dofs;
	const PetscInt *lxp,*lyp,*lzp;
	PetscInt MX,MY,MZ,p,Mp,Np,Pp,*lxv,*lyv,*lzv,i;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	MX = mx;
	MY = my;
	MZ = mz;

#if 0
	/* pressure */
	pbasis_dofs = P_BASIS_FUNCTIONS;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX,MY,MZ, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, pbasis_dofs,0, 0,0,0, &dap );CHKERRQ(ierr);
	ierr = DMDASetElementType_P1(dap);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(dap,&lxp,&lyp,&lzp);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dap,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Mp+1),&lxv);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Np+1),&lyv);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Pp+1),&lzv);CHKERRQ(ierr);
	for (p=0; p<Mp; p++) {
		lxv[p] = lxp[p] * 2;
	} lxv[Mp-1]++;
	for (p=0; p<Np; p++) {
		lyv[p] = lyp[p] * 2;
	} lyv[Np-1]++;
	for (p=0; p<Pp; p++) {
		lzv[p] = lzp[p] * 2;
	} lzv[Pp-1]++;
	
	/* velocity */
	vbasis_dofs = 3;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, 2*MX+1,2*MY+1,2*MZ+1, Mp,Np,Pp, vbasis_dofs,2, lxv,lyv,lzv, &dav );CHKERRQ(ierr);
	ierr = DMDASetElementType_Q2(dav);CHKERRQ(ierr);
	ierr = PetscFree(lxv);CHKERRQ(ierr);
	ierr = PetscFree(lyv);CHKERRQ(ierr);
	ierr = PetscFree(lzv);CHKERRQ(ierr);
#endif


	/* velocity */
	vbasis_dofs = 3;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, 2*MX+1,2*MY+1,2*MZ+1, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, vbasis_dofs,2, PETSC_NULL,PETSC_NULL,PETSC_NULL, &dav );CHKERRQ(ierr);
	ierr = DMDASetElementType_Q2(dav);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dav,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRangesElementQ2(dav, 0,0,0, 0,0,0, &lxv,&lyv,&lzv);CHKERRQ(ierr);
	
	/* pressure */
	pbasis_dofs = P_BASIS_FUNCTIONS;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX,MY,MZ, Mp,Np,Pp, pbasis_dofs,0, lxv,lyv,lzv, &dap );CHKERRQ(ierr);
	ierr = DMDASetElementType_P1(dap);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(dap,&lxp,&lyp,&lzp);CHKERRQ(ierr);
	
	/* set an initial geometry */
	ierr = DMDASetUniformCoordinates(dav,0.0,1.0, 0.0,1.0, 0.0,1.0);CHKERRQ(ierr);
	
	/* stokes */
	ierr = DMCompositeCreate(PETSC_COMM_WORLD,&multipys_pack);CHKERRQ(ierr);
	ierr = DMCompositeAddDM(multipys_pack,dav);CHKERRQ(ierr);	
	ierr = DMCompositeAddDM(multipys_pack,dap);CHKERRQ(ierr);	
	
	ierr = DMDASetFieldName(dav,0,"ux");CHKERRQ(ierr);
	ierr = DMDASetFieldName(dav,1,"uy");CHKERRQ(ierr);
	ierr = DMDASetFieldName(dav,2,"uz");CHKERRQ(ierr);
	switch (P_BASIS_FUNCTIONS) {
		case 1:
			ierr = DMDASetFieldName(dap,0,"P0_p");CHKERRQ(ierr);
			break;
		case 4:
			ierr = DMDASetFieldName(dap,0,"P1_p");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,1,"P1_dpdx");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,2,"P1_dpdy");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,3,"P1_dpdz");CHKERRQ(ierr);
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Pressure space may ONLY contain 1 or 4 basis functions");
			break;
	}
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dav,"stk_velocity_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dap,"stk_pressure_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)multipys_pack,"stk_pack_");CHKERRQ(ierr);
	
	ctx->dav  = dav;
	ctx->dap  = dap;
	ctx->stokes_pack = multipys_pack;

	ierr = PetscFree(lxv);CHKERRQ(ierr);
	ierr = PetscFree(lyv);CHKERRQ(ierr);
	ierr = PetscFree(lzv);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateBoundaryList_Stokes"
PetscErrorCode PhysCompCreateBoundaryList_Stokes(PhysCompStokes ctx)
{
	DM dav;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* vel bc's */
	dav = ctx->dav;
	if (!dav) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"dav must be set"); }
	
	ierr = DMDABCListCreate(dav,&ctx->u_bclist);CHKERRQ(ierr);
	
	/* pressure bc's */
	ctx->p_bclist = PETSC_NULL;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateVolumeQuadrature_Stokes"
PetscErrorCode PhysCompCreateVolumeQuadrature_Stokes(PhysCompStokes ctx)
{
	DM dav;
	PetscInt lmx,lmy,lmz;
	PetscInt np_per_dim,ncells;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	dav = ctx->dav;

	np_per_dim = 3;
  ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	ncells = lmx * lmy * lmz;
	ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&ctx->volQ);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetValuesLocalStencil_AddValues_Stokes_Velocity"
PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[])
{
  PetscInt n,idx;
	
  PetscFunctionBegin;
  for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
		idx = u_eqn[3*n  ];
    fields_F[idx] += Fe_u[NSD*n  ];
		
		idx = u_eqn[3*n+1];
    fields_F[idx] += Fe_u[NSD*n+1];
		
		idx = u_eqn[3*n+2];
    fields_F[idx] += Fe_u[NSD*n+2];
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetValuesLocalStencil_InsertValues_Stokes_Velocity"
PetscErrorCode DMDASetValuesLocalStencil_InsertValues_Stokes_Velocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[])
{
  PetscInt n,idx;
	
  PetscFunctionBegin;
  for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
		idx = u_eqn[3*n  ];
    fields_F[idx] = Fe_u[NSD*n  ];
		
		idx = u_eqn[3*n+1];
    fields_F[idx] = Fe_u[NSD*n+1];
		
		idx = u_eqn[3*n+2];
    fields_F[idx] = Fe_u[NSD*n+2];
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetValuesLocalStencil_AddValues_Stokes_Pressure"
PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(PetscScalar *fields_F,PetscInt p_eqn[],PetscScalar Fe_p[])
{
  PetscInt n,idx;
	
  PetscFunctionBegin;
  for (n = 0; n<P_BASIS_FUNCTIONS; n++) {
		idx = p_eqn[n];
    fields_F[idx] += Fe_p[n];
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetValuesLocalStencil_AddValues_Stokes_ScalarVelocity"
PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_ScalarVelocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[])
{
  PetscInt n,idx;
	
  PetscFunctionBegin;
  for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
		idx = u_eqn[n];
    fields_F[idx] += Fe_u[n];
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompLoadMesh_Stokes3d"
PetscErrorCode PhysCompLoadMesh_Stokes3d(PhysCompStokes ctx,const char fname_vel[],const char fname_p[])
{
	DM dav,dap,multipys_pack;
	PetscInt vbasis_dofs;
	PetscInt pbasis_dofs;
	PetscInt p,Mp,Np,Pp,*lxv,*lyv,*lzv,i,MX,MY,MZ;
	const PetscInt *lxp,*lyp,*lzp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

#if 0	
	/* pressure */
	ierr = DMDACreateFromPackDataToFile(PETSC_COMM_WORLD,fname_p,&dap);CHKERRQ(ierr);
	ierr = DMDASetElementType_P1(dap);CHKERRQ(ierr);

	ierr = DMDAGetOwnershipRanges(dap,&lxp,&lyp,&lzp);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dap,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Mp+1),&lxv);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Np+1),&lyv);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Pp+1),&lzv);CHKERRQ(ierr);
	for (p=0; p<Mp; p++) {
		lxv[p] = lxp[p] * 2;
	} lxv[Mp-1]++;
	for (p=0; p<Np; p++) {
		lyv[p] = lyp[p] * 2;
	} lyv[Np-1]++;
	for (p=0; p<Pp; p++) {
		lzv[p] = lzp[p] * 2;
	} lzv[Pp-1]++;

	
	/* velocity */
	/* FUCK ME, the layout changes and the Q2-P1 elements won't line up anymore...lame */
//	ierr = DMDACreateFromPackDataToFile(PETSC_COMM_WORLD,fname_vel,&dav);CHKERRQ(ierr);
//	ierr = DMDASetElementType_Q2(dav);CHKERRQ(ierr);

	MX = ctx->mx;
	MY = ctx->my;
	MZ = ctx->mz;
	
	vbasis_dofs = 3;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, 2*MX+1,2*MY+1,2*MZ+1, Mp,Np,Pp, vbasis_dofs,2, lxv,lyv,lzv, &dav );CHKERRQ(ierr);
	ierr = DMDASetElementType_Q2(dav);CHKERRQ(ierr);
	ierr = PetscFree(lxv);CHKERRQ(ierr);
	ierr = PetscFree(lyv);CHKERRQ(ierr);
	ierr = PetscFree(lzv);CHKERRQ(ierr);
#endif
	
	
	/* velocity */
	vbasis_dofs = 3;
	ierr = DMDACreateFromPackDataToFile(PETSC_COMM_WORLD,fname_vel,&dav);CHKERRQ(ierr);
	/* the above function call will load the initial geometry */
	ierr = DMDASetElementType_Q2(dav);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dav,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRangesElementQ2(dav, 0,0,0, 0,0,0, &lxv,&lyv,&lzv);CHKERRQ(ierr);
	
	ierr = DMDAGetSizeElementQ2(dav,&ctx->mx,&ctx->my,&ctx->mz);CHKERRQ(ierr);
	
	/* pressure */
	MX = ctx->mx;
	MY = ctx->my;
	MZ = ctx->mz;
	
	pbasis_dofs = P_BASIS_FUNCTIONS;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX,MY,MZ, Mp,Np,Pp, pbasis_dofs,0, lxv,lyv,lzv, &dap );CHKERRQ(ierr);
	ierr = DMDASetElementType_P1(dap);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(dap,&lxp,&lyp,&lzp);CHKERRQ(ierr);
	
	ierr = PetscFree(lxv);CHKERRQ(ierr);
	ierr = PetscFree(lyv);CHKERRQ(ierr);
	ierr = PetscFree(lzv);CHKERRQ(ierr);
	
	
	/* stokes */
	ierr = DMCompositeCreate(PETSC_COMM_WORLD,&multipys_pack);CHKERRQ(ierr);
	ierr = DMCompositeAddDM(multipys_pack,dav);CHKERRQ(ierr);	
	ierr = DMCompositeAddDM(multipys_pack,dap);CHKERRQ(ierr);	
	
	ierr = DMDASetFieldName(dav,0,"ux");CHKERRQ(ierr);
	ierr = DMDASetFieldName(dav,1,"uy");CHKERRQ(ierr);
	ierr = DMDASetFieldName(dav,2,"uz");CHKERRQ(ierr);
	switch (P_BASIS_FUNCTIONS) {
		case 1:
			ierr = DMDASetFieldName(dap,0,"P0_p");CHKERRQ(ierr);
			break;
		case 4:
			ierr = DMDASetFieldName(dap,0,"P1_p");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,1,"P1_dpdx");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,2,"P1_dpdy");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,3,"P1_dpdz");CHKERRQ(ierr);
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Pressure space may ONLY contain 1 or 4 basis functions");
			break;
	}
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dav,"stk_velocity_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dap,"stk_pressure_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)multipys_pack,"stk_pack_");CHKERRQ(ierr);
	
	ctx->dav  = dav;
	ctx->dap  = dap;
	ctx->stokes_pack = multipys_pack;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompSaveMesh_Stokes3d"
PetscErrorCode PhysCompSaveMesh_Stokes3d(PhysCompStokes ctx,const char fname_vel[],const char fname_p[],const char fname_coors[])
{
	DM dav,dap;
	Vec coords;
	PetscViewer viewer;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	dav = ctx->dav;
	dap = ctx->dap;

	/* velocity */
	ierr = DMDAPackDataToFile(dav,fname_vel);CHKERRQ(ierr);

	/* pressure */
	ierr = DMDAPackDataToFile(dap,fname_p);CHKERRQ(ierr);
	
	/* coords */
	if (fname_coors) {
		ierr = DMDAGetCoordinates(dav,&coords);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen( ((PetscObject)dav)->comm,fname_coors,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
		ierr = VecView(coords,viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinStokesKSPMonitorBlocks"
PetscErrorCode pTatinStokesKSPMonitorBlocks(KSP ksp,PetscInt n,PetscReal rnorm,void *data)
{
	PetscErrorCode ierr;
	pTatinCtx ctx;
	PetscReal norms[4];
	Vec X,Xu,Xp,v,w;
	Mat A;
	
	PetscFunctionBegin;
	ctx = (pTatinCtx)data;
	ierr = KSPGetOperators(ksp,&A,0,0);CHKERRQ(ierr);
	ierr = MatGetVecs(A,&w,&v);CHKERRQ(ierr);
	
	ierr = KSPBuildResidual(ksp,v,w,&X);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
	
	ierr = VecStrideNorm(Xu,0,NORM_2,&norms[0]);CHKERRQ(ierr);
	ierr = VecStrideNorm(Xu,1,NORM_2,&norms[1]);CHKERRQ(ierr);
	ierr = VecStrideNorm(Xu,2,NORM_2,&norms[2]);CHKERRQ(ierr);
	ierr = VecNorm(Xp,NORM_2,&norms[3]);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
	ierr = VecDestroy(&v);CHKERRQ(ierr);
	ierr = VecDestroy(&w);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"%3D KSP Component U,V,W,P residual norm [ %1.12e, %1.12e, %1.12e, %1.12e ]\n",n,norms[0],norms[1],norms[2],norms[3]);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureCreate_GaussLegendreStokes"
PetscErrorCode VolumeQuadratureCreate_GaussLegendreStokes(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature)
{
	Quadrature Q;
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
	
	ierr = QuadratureCreate(&Q);CHKERRQ(ierr);
	Q->dim  = nsd;
	Q->type = VOLUME_QUAD;
	
	PetscPrintf(PETSC_COMM_WORLD,"VolumeQuadratureCreate_GaussLegendreStokes:\n");
	switch (np_per_dim) {
		case 1:
			PetscPrintf(PETSC_COMM_WORLD,"\tUsing 1 pnt Gauss Legendre quadrature\n");
			//QuadratureCreateGauss_1pnt_3D(&ngp,gp_xi,gp_weight);
			break;
			
		case 2:
			PetscPrintf(PETSC_COMM_WORLD,"\tUsing 2x2 pnt Gauss Legendre quadrature\n");
			QuadratureCreateGauss_2pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
			
		case 3:
			PetscPrintf(PETSC_COMM_WORLD,"\tUsing 3x3 pnt Gauss Legendre quadrature\n");
			QuadratureCreateGauss_3pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
			
		default:
			PetscPrintf(PETSC_COMM_WORLD,"\tUsing 3x3 pnt Gauss Legendre quadrature\n");
			QuadratureCreateGauss_3pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
	}
	
	Q->n_elements = ncells;
	if (ncells!=0) {
		
		DataBucketCreate(&Q->properties_db);
		DataBucketRegisterField(Q->properties_db,QPntVolCoefStokes_classname, sizeof(QPntVolCoefStokes),PETSC_NULL);
		DataBucketFinalize(Q->properties_db);
		
		DataBucketSetInitialSizes(Q->properties_db,Q->npoints*ncells,1);
		
		DataBucketView(PETSC_COMM_WORLD, Q->properties_db,"GaussLegendre StokesCoefficients",DATABUCKET_VIEW_STDOUT);
	}
	
	*quadrature = Q;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureGetAllCellData_Stokes"
PetscErrorCode VolumeQuadratureGetAllCellData_Stokes(Quadrature Q,QPntVolCoefStokes *coeffs[])
{
	QPntVolCoefStokes *quadraturepoint_data;
  DataField          PField;
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(Q->properties_db, QPntVolCoefStokes_classname ,&PField);
	quadraturepoint_data = PField->data;
	*coeffs = quadraturepoint_data;
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureGetCellData_Stokes"
PetscErrorCode VolumeQuadratureGetCellData_Stokes(Quadrature Q,QPntVolCoefStokes coeffs[],PetscInt cidx,QPntVolCoefStokes *cell[])
{
  PetscFunctionBegin;
	if (cidx>=Q->n_elements) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_SIZ,"cidx > max cells");
	}
	
	*cell = &coeffs[cidx*Q->npoints];
  PetscFunctionReturn(0);
}



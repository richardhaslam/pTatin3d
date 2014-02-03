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
 **    Filename:      monitors.c
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



#include "petsc-private/dmdaimpl.h" 
#include "petscdm.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "output_paraview.h"


#undef __FUNCT__
#define __FUNCT__ "pTatin_KSPMonitor_StdoutStokesResiduals3d"
PetscErrorCode pTatin_KSPMonitor_StdoutStokesResiduals3d(KSP ksp,PetscInt n,PetscReal rnorm,void *data)
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
#define __FUNCT__ "pTatin_KSPMonitor_ParaviewStokesResiduals3d"
PetscErrorCode pTatin_KSPMonitor_ParaviewStokesResiduals3d(KSP ksp,PetscInt n,PetscReal rnorm,void *data)
{
	PetscErrorCode ierr;
	pTatinCtx ctx;
	DM stokes_pack;
	Vec X,v,w,UP;
	Mat A;
	static char pvdfilename[1000];
	char vtkfilename[1000];
	PetscInt its;
	
	PetscFunctionBegin;
	ctx = (pTatinCtx)data;
	ierr = KSPGetOperators(ksp,&A,0,0);CHKERRQ(ierr);
	ierr = MatGetVecs(A,&w,&v);CHKERRQ(ierr);
	
	ierr = KSPBuildResidual(ksp,v,w,&X);CHKERRQ(ierr);

	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

	if (its==0) {
		sprintf(pvdfilename,"%s/residualseries_vp_residuals_step%d.pvd",ctx->outputpath,ctx->step);
		PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename );
		ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
	}
	sprintf(vtkfilename, "iteration%d_vp_residuals_step%d.pvts",its,ctx->step);
	ierr = ParaviewPVDAppend(pvdfilename,its, vtkfilename, "");CHKERRQ(ierr);
	
	
	// PVTS + VTS
	sprintf(vtkfilename, "iteration%d_vp_residuals_step%d",its,ctx->step);
	
	stokes_pack = ctx->stokes_ctx->stokes_pack;
	UP = X;
	ierr = pTatinOutputParaViewMeshVelocityPressure(stokes_pack,UP,ctx->outputpath,vtkfilename);CHKERRQ(ierr);

	
	ierr = VecDestroy(&v);CHKERRQ(ierr);
	ierr = VecDestroy(&w);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"%3D KSP Residual: ptatin linear solution viewer wrote file \n",its);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin_SNESMonitor_ParaviewStokesResiduals3d"
PetscErrorCode pTatin_SNESMonitor_ParaviewStokesResiduals3d(SNES snes,PetscInt n,PetscReal rnorm,void *data)
{
	PetscErrorCode ierr;
	pTatinCtx ctx;
	DM stokes_pack;
	Vec X,UP;
	const char *prefix;
	static char pvdfilename[1000];
	char vtkfilename[1000];
	PetscInt its;
	
	PetscFunctionBegin;
	ctx = (pTatinCtx)data;

	ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);
	ierr = SNESGetSolution(snes,&X);CHKERRQ(ierr);
	ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
	
	if (its==0) {
		if (!prefix) {
			sprintf(pvdfilename,"%s/residualseries_vp_snes_residuals_step%d.pvd",ctx->outputpath,ctx->step);
		}else {
			sprintf(pvdfilename,"%s/residualseries_vp_snes_%sresiduals_step%d.pvd",ctx->outputpath,prefix,ctx->step);
		}
		PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename );
		ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
	}
	if (!prefix) {
		sprintf(vtkfilename, "snes_iteration%d_vp_residuals_step%d.pvts",its,ctx->step);
	} else {
		sprintf(vtkfilename, "snes_%siteration%d_vp_residuals_step%d.pvts",prefix,its,ctx->step);
	}
	ierr = ParaviewPVDAppend(pvdfilename,its, vtkfilename, "");CHKERRQ(ierr);
	
	
	// PVTS + VTS
	if (!prefix) {
		sprintf(vtkfilename, "snes_iteration%d_vp_residuals_step%d",its,ctx->step);
	} else {
		sprintf(vtkfilename, "snes_%siteration%d_vp_residuals_step%d",prefix,its,ctx->step);
	}
	
	stokes_pack = ctx->stokes_ctx->stokes_pack;
	UP = X;
	ierr = pTatinOutputParaViewMeshVelocityPressure(stokes_pack,UP,ctx->outputpath,vtkfilename);CHKERRQ(ierr);

	if (!prefix) {
		PetscPrintf(PETSC_COMM_WORLD,"%3D SNES Residual: ptatin non-linear solution viewer wrote file \n",its);
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"%3D SNES(%s) Residual: ptatin non-linear solution viewer wrote file \n",its,prefix);
	}	
	PetscFunctionReturn(0);
}

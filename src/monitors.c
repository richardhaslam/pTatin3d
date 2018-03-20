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
 **    filename:   monitors.c
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

#include "petscdm.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "petsc_utils.h"
#include "output_paraview.h"


PetscErrorCode pTatin_SNESMonitor_StdoutStokesResiduals3d(SNES snes,PetscInt n,PetscReal fnorm,void *data)
{
  PetscErrorCode ierr;
  pTatinCtx      ctx;
  PetscReal      norms[4];
  Vec            F,Fu,Fp;
  MPI_Comm       comm;

  PetscFunctionBegin;
  ctx = (pTatinCtx)data;
  ierr = SNESGetFunction(snes,&F,NULL,NULL);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(ctx->stokes_ctx->stokes_pack,F,&Fu,&Fp);CHKERRQ(ierr);

  ierr = VecStrideNormLocal(Fu,0,NORM_2,&norms[0]);CHKERRQ(ierr);
  ierr = VecStrideNormLocal(Fu,1,NORM_2,&norms[1]);CHKERRQ(ierr);
  ierr = VecStrideNormLocal(Fu,2,NORM_2,&norms[2]);CHKERRQ(ierr);
  ierr = VecNormLocal(Fp,NORM_2,&norms[3]);CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(ctx->stokes_ctx->stokes_pack,F,&Fu,&Fp);CHKERRQ(ierr);

  ierr = PetscObjectGetComm((PetscObject)F,&comm);CHKERRQ(ierr);
  PetscPrintf(comm,"%3D SNES Component Fu,Fv,Fw,Fp function norm [ %1.12e, %1.12e, %1.12e, %1.12e ]\n",n,norms[0],norms[1],norms[2],norms[3]);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_KSPMonitor_StdoutStokesResiduals3d(KSP ksp,PetscInt n,PetscReal rnorm,void *data)
{
  PetscErrorCode ierr;
  pTatinCtx      ctx;
  PetscReal      norms[4];
  Vec            X,Xu,Xp,v,w;
  Mat            A;
  MPI_Comm       comm;

  PetscFunctionBegin;
  ctx = (pTatinCtx)data;
  ierr = KSPGetOperators(ksp,&A,0);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&w,&v);CHKERRQ(ierr);

  ierr = KSPBuildResidual(ksp,v,w,&X);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);

  ierr = VecStrideNorm(Xu,0,NORM_2,&norms[0]);CHKERRQ(ierr);
  ierr = VecStrideNorm(Xu,1,NORM_2,&norms[1]);CHKERRQ(ierr);
  ierr = VecStrideNorm(Xu,2,NORM_2,&norms[2]);CHKERRQ(ierr);
  ierr = VecNorm(Xp,NORM_2,&norms[3]);CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);

  ierr = PetscObjectGetComm((PetscObject)ksp,&comm);CHKERRQ(ierr);
  PetscPrintf(comm,"%3D KSP Component U,V,W,P residual norm [ %1.12e, %1.12e, %1.12e, %1.12e ]\n",n,norms[0],norms[1],norms[2],norms[3]);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_KSPMonitor_ParaviewStokesResiduals3d(KSP ksp,PetscInt n,PetscReal rnorm,void *data)
{
  PetscErrorCode ierr;
  pTatinCtx      ctx;
  DM             stokes_pack;
  Vec            X,v,w,UP;
  Mat            A;
  static char    pvdfilename[PETSC_MAX_PATH_LEN];
  char           vtkfilename[PETSC_MAX_PATH_LEN];
  PetscInt       its;
  MPI_Comm       comm;

  PetscFunctionBegin;
  ctx = (pTatinCtx)data;
  ierr = KSPGetOperators(ksp,&A,0);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&w,&v);CHKERRQ(ierr);

  ierr = KSPBuildResidual(ksp,v,w,&X);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)ksp,&comm);CHKERRQ(ierr);

  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

  if (its == 0) {
    PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/stokes_ksp_r_step%.6d.pvd",ctx->outputpath,ctx->step);
    PetscPrintf(comm,"  writing pvdfilename %s \n",pvdfilename );
    ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
  }
  PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"stokes_ksp_r_it%.4d_step%.6d.pvts",its,ctx->step);
  ierr = ParaviewPVDAppend(pvdfilename,its,vtkfilename,"");CHKERRQ(ierr);

  // PVTS + VTS
  PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"stokes_ksp_r_it%.4d_step%.6d",its,ctx->step);

  stokes_pack = ctx->stokes_ctx->stokes_pack;
  UP = X;
  ierr = pTatinOutputParaViewMeshVelocityPressure(stokes_pack,UP,ctx->outputpath,vtkfilename);CHKERRQ(ierr);

  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);

  PetscPrintf(comm,"%3D KSP Residual: ptatin linear solution viewer wrote file \n",its);

  PetscFunctionReturn(0);
}

PetscErrorCode _pTatin_SNESMonitorStokes_Paraview(SNES snes,pTatinCtx ctx,Vec X,const char field[])
{
  PetscErrorCode ierr;
  DM             stokes_pack;
  const char     *prefix;
  static char    pvdfilename[PETSC_MAX_PATH_LEN];
  char           vtkfilename[PETSC_MAX_PATH_LEN];
  PetscInt       its;
  MPI_Comm       comm;

  PetscFunctionBegin;
  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

  if (its == 0) {
    if (!prefix) {
      PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/stokes_snes_%s_step%.6d.pvd",ctx->outputpath,field,ctx->step);
    } else {
      PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/stokes_snes_%s%s_step%.6d.pvd",ctx->outputpath,prefix,field,ctx->step);
    }
    PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n",pvdfilename );
    ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
  }
  if (!prefix) {
    PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"stokes_snes_%s_it%.4d_step%.6d.pvts",field,its,ctx->step);
  } else {
    PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"stokes_snes_%s%s_it%.4d_step%.6d.pvts",prefix,field,its,ctx->step);
  }
  ierr = ParaviewPVDAppend(pvdfilename,(double)its,vtkfilename,"");CHKERRQ(ierr);

  // PVTS + VTS
  if (!prefix) {
    PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"stokes_snes_%s_it%.4d_step%.6d",field,its,ctx->step);
  } else {
    PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"stokes_snes_%s%s_it%.4d_step%.6d",prefix,field,its,ctx->step);
  }

  stokes_pack = ctx->stokes_ctx->stokes_pack;
  ierr = pTatinOutputParaViewMeshVelocityPressure(stokes_pack,X,ctx->outputpath,vtkfilename);CHKERRQ(ierr);

  ierr = PetscObjectGetComm((PetscObject)X,&comm);CHKERRQ(ierr);
  if (!prefix) {
    PetscPrintf(comm,"%3D SNES [field %s]: ptatin non-linear solution viewer wrote file \n",its,field);
  } else {
    PetscPrintf(comm,"%3D SNES(%s) [field %s]: ptatin non-linear solution viewer wrote file \n",its,field,prefix);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_SNESMonitorStokes_Solution_Paraview(SNES snes,PetscInt n,PetscReal rnorm,void *data)
{
  PetscErrorCode ierr;
  pTatinCtx      ctx;
  Vec            X;

  PetscFunctionBegin;
  ctx = (pTatinCtx)data;

  ierr = SNESGetSolution(snes,&X);CHKERRQ(ierr);
  ierr = _pTatin_SNESMonitorStokes_Paraview(snes,ctx,X,"X");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_SNESMonitorStokes_SolutionUpdate_Paraview(SNES snes,PetscInt n,PetscReal rnorm,void *data)
{
  PetscErrorCode ierr;
  pTatinCtx      ctx;
  Vec            X;

  PetscFunctionBegin;
  ctx = (pTatinCtx)data;

  ierr = SNESGetSolutionUpdate(snes,&X);CHKERRQ(ierr);
  ierr = _pTatin_SNESMonitorStokes_Paraview(snes,ctx,X,"dX");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_SNESMonitorStokes_Residual_Paraview(SNES snes,PetscInt n,PetscReal rnorm,void *data)
{
  PetscErrorCode ierr;
  pTatinCtx      ctx;
  Vec            X;

  PetscFunctionBegin;
  ctx = (pTatinCtx)data;

  ierr = SNESGetFunction(snes,&X,0,0);CHKERRQ(ierr);
  ierr = _pTatin_SNESMonitorStokes_Paraview(snes,ctx,X,"F");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* monitors */
PetscErrorCode pTatin_Stokes_ActivateMonitors(pTatinCtx user,SNES snes)
{
  PetscErrorCode ierr;
  PetscBool      monitor_set;
  KSP            ksp;

  PetscFunctionBegin;
  monitor_set = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-stokes_snes_monitor_up",&monitor_set,NULL);CHKERRQ(ierr);
  if (monitor_set) { ierr = SNESMonitorSet(snes,pTatin_SNESMonitor_StdoutStokesResiduals3d,(void*)user,NULL);CHKERRQ(ierr); }

  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);

  monitor_set = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-stokes_ksp_monitor_up",&monitor_set,NULL);CHKERRQ(ierr);
  if (monitor_set) { ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_StdoutStokesResiduals3d,(void*)user,NULL);CHKERRQ(ierr); }

  monitor_set = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-stokes_ksp_monitor_paraview",&monitor_set,NULL);CHKERRQ(ierr);
  if (monitor_set) { ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_ParaviewStokesResiduals3d,(void*)user,NULL);CHKERRQ(ierr); }

  monitor_set = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-stokes_snes_monitor_F_paraview",&monitor_set,NULL);CHKERRQ(ierr);
  if (monitor_set) { ierr = SNESMonitorSet(snes,pTatin_SNESMonitorStokes_Residual_Paraview,(void*)user,NULL);CHKERRQ(ierr); }

  monitor_set = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-stokes_snes_monitor_dX_paraview",&monitor_set,NULL);CHKERRQ(ierr);
  if (monitor_set) { ierr = SNESMonitorSet(snes,pTatin_SNESMonitorStokes_SolutionUpdate_Paraview,(void*)user,NULL);CHKERRQ(ierr); }

  monitor_set = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-stokes_snes_monitor_X_paraview",&monitor_set,NULL);CHKERRQ(ierr);
  if (monitor_set) { ierr = SNESMonitorSet(snes,pTatin_SNESMonitorStokes_Solution_Paraview,(void*)user,NULL);CHKERRQ(ierr); }

  PetscFunctionReturn(0);
}

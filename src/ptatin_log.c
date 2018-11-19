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
 **    filename:   ptatin_log.c
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

#include "petsc.h"
#include "petscksp.h"
#include "petscsnes.h"

#include "ptatin3d.h"
#include "ptatin_version_info.h"
#include "ptatin_utils.h"
#include "private/ptatin_impl.h"


PetscErrorCode pTatinLogOpenFile(pTatinCtx ctx)
{
  PetscBool      stdout = PETSC_FALSE;
  char           name[PETSC_MAX_PATH_LEN];
  char           date_time[1024];
  PetscErrorCode ierr;

  pTatinGenerateFormattedTimestamp(date_time);
  ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/ptatin.log-%s",ctx->outputpath,date_time);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,NULL,"-ptatin_log_stdout",&stdout,NULL);CHKERRQ(ierr);
  if (!stdout) {
    if (ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is already open");
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,name,&ctx->log);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"[pTatin] Created log file: %s\n",name);
  }
  else {
    ctx->log = PETSC_VIEWER_STDOUT_WORLD;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogCloseFile(pTatinCtx ctx)
{
  PetscBool      stdout = PETSC_FALSE;
  PetscErrorCode ierr;

  ierr = PetscOptionsGetBool(NULL,NULL,"-ptatin_log_stdout",&stdout,NULL);CHKERRQ(ierr);
  if (!stdout) {
    ierr = PetscViewerDestroy(&ctx->log);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogHeader(pTatinCtx ctx)
{
  char username[PETSC_MAX_PATH_LEN];
  char date[PETSC_MAX_PATH_LEN];
  char machine[PETSC_MAX_PATH_LEN];
  char prgname[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  /* write header into options file */
  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  ierr = PetscGetUserName(username,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
  ierr = PetscGetDate(date,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
  ierr = PetscGetHostName(machine,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
  ierr = PetscGetProgramName(prgname,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"** ===================================================================================\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**   pTatin3d Log File\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**     %s\n", PTATIN_VERSION_CNTR_REPO);
#ifdef PTATIN_DEVELOPMENT_VERSION
  ierr = PetscViewerASCIIPrintf(ctx->log,"**     %s\n", PTATIN_VERSION_CNTR_REVISION);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**     %s\n", PTATIN_VERSION_CNTR_LOG);
  #ifdef PTATIN_GIT_REPO_STATUS
  ierr = PetscViewerASCIIPrintf(ctx->log,"**     %s\n", PTATIN_GIT_REPO_STATUS);
  #endif
#endif
#ifdef PTATIN_RELEASE
  ierr = PetscViewerASCIIPrintf(ctx->log,"**     %s\n", PTATIN_VERSION_CNTR_REVISION);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**     Release v%d.%d-p%d\n", PTATIN_VERSION_MAJOR,PTATIN_VERSION_MINOR,PTATIN_VERSION_PATCH);
#endif
  ierr = PetscViewerASCIIPrintf(ctx->log,"**\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**   Generated by user: %s\n",username);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**   Date             : %s\n",date);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**   Machine          : %s\n",machine);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**   Driver           : %s\n",prgname);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"**\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"** ===================================================================================\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ctx->log,"\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogBasic(pTatinCtx ctx)
{
  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  PetscViewerASCIIPrintf(ctx->log,"------------------------------------------------------------------------------------------\n");
  PetscViewerASCIIPrintf(ctx->log,"  time step %6d\n", ctx->step);
  PetscViewerASCIIPrintf(ctx->log,"  time %1.4e;  dt %1.4e;  dt_min %1.4e;  dt_max %1.4e\n", ctx->time, ctx->dt, ctx->dt_min, ctx->dt_max);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogBasicKSP(pTatinCtx ctx,const char kspname[],KSP ksp)
{
  PetscReal rnorm;
  PetscInt its;
  const char *prefix;
  KSPConvergedReason reason;
  PetscErrorCode ierr;

  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  ierr = KSPGetOptionsPrefix(ksp,&prefix);CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp,&rnorm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);

  PetscViewerASCIIPrintf(ctx->log,"  KSP: (%18.18s)[prefix %8.8s] residual %1.4e;  iterations %1.4d;  reason %s;\n", kspname,prefix,rnorm,its,KSPConvergedReasons[reason]);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogBasicSNES(pTatinCtx ctx,const char snesname[],SNES snes)
{
  Vec       r;
  PetscReal rnorm;
  PetscInt its,lits,nkspfails,nFevals,nstepfails;
  const char *prefix;
  SNESConvergedReason reason;
  PetscBool same;
  PetscErrorCode ierr;

  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);
  ierr = SNESGetFunction(snes,&r,NULL,NULL);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&rnorm);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetLinearSolveIterations(snes,&lits);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);

  ierr = PetscObjectTypeCompare((PetscObject)snes,SNESKSPONLY,&same);CHKERRQ(ierr);

  if (!same) {
    PetscViewerASCIIPrintf(ctx->log,"  SNES: (%18.18s)[prefix %8.8s] residual %1.4e;  iterations %1.4d;  total linear its. %1.4d;  reason %s;\n", snesname,prefix,rnorm,its,lits,SNESConvergedReasons[reason]);

    ierr = SNESGetLinearSolveFailures(snes,&nkspfails);CHKERRQ(ierr);
    ierr = SNESGetNonlinearStepFailures(snes,&nstepfails);CHKERRQ(ierr);
    ierr = SNESGetNumberFunctionEvals(snes,&nFevals);CHKERRQ(ierr);
    PetscViewerASCIIPrintf(ctx->log,"                                              function evals %1.4d;  ksp failures %1.4d;  step failures %1.4d\n", nFevals, nkspfails, nstepfails);
  }
  else if (same) {
    char kspname[] = "snes->ksp";
    KSP ksp;

    PetscViewerASCIIPrintf(ctx->log,"  SNES: (%18.18s)[prefix %8.8s] -> ksponly\n", snesname,prefix);

    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = pTatinLogBasicKSP(ctx,kspname,ksp);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogBasicStokesSolution(pTatinCtx ctx,DM pack,Vec X)
{
  char fieldname[120];
  Vec velocity,pressure;
  PetscReal minval,maxval,norm2,norm1;
  PetscInt loc;
  PetscErrorCode ierr;

  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  ierr = DMCompositeGetAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);

  PetscViewerASCIIPrintf( ctx->log, "  Field         min          max          norm_2       norm_1\n" );
  VecMin( velocity, &loc, &minval );
  VecMax( velocity, &loc, &maxval );
  VecNorm( velocity, NORM_2, &norm2 );
  VecNorm( velocity, NORM_1, &norm1 );

  sprintf(fieldname,"velocity");
  PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e\n", fieldname,minval,maxval,norm2,norm1 );

  VecMin( pressure, &loc, &minval );
  VecMax( pressure, &loc, &maxval );
  VecNorm( pressure, NORM_2, &norm2 );
  VecNorm( pressure, NORM_1, &norm1 );
  sprintf(fieldname,"pressure");
  PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e\n", fieldname,minval,maxval,norm2,norm1 );

  VecMin( X, &loc, &minval );
  VecMax( X, &loc, &maxval );
  VecNorm( X, NORM_2, &norm2 );
  VecNorm( X, NORM_1, &norm1 );
  sprintf(fieldname,"X");
  PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e\n", fieldname,minval,maxval,norm2,norm1 );

  ierr = DMCompositeRestoreAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinViewBasicStokesSolution(pTatinCtx ctx,DM pack,Vec X)
{
  PetscViewer  tmp;
  PetscErrorCode ierr;

  tmp = ctx->log;
  ctx->log  = PETSC_VIEWER_STDOUT_WORLD;
  PetscPrintf(PETSC_COMM_WORLD,"pTatinViewBasicStokesSolution:\n");
  ierr = pTatinLogBasicStokesSolution(ctx,pack,X);CHKERRQ(ierr);
  ctx->log = tmp;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogBasicStokesSolutionResiduals(pTatinCtx ctx,SNES snes,DM pack,Vec X)
{
  char fieldname[120];
  Vec F,Fu,Fp;
  PetscReal minval,maxval,norm2,norm1;
  PetscInt loc;
  PetscErrorCode ierr;

  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
//  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
  ierr = SNESGetFunction(snes,&F,NULL,NULL);CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);

  ierr = DMCompositeGetAccess(pack,F,&Fu,&Fp);CHKERRQ(ierr);

  PetscViewerASCIIPrintf( ctx->log, "  Field         min          max          norm_2       norm_1\n" );
  VecMin( Fu, &loc, &minval );
  VecMax( Fu, &loc, &maxval );
  VecNorm( Fu, NORM_2, &norm2 );
  VecNorm( Fu, NORM_1, &norm1 );

  sprintf(fieldname,"Fu");
  PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e\n", fieldname,minval,maxval,norm2,norm1 );

  VecMin( Fp, &loc, &minval );
  VecMax( Fp, &loc, &maxval );
  VecNorm( Fp, NORM_2, &norm2 );
  VecNorm( Fp, NORM_1, &norm1 );
  sprintf(fieldname,"Fp");
  PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e\n", fieldname,minval,maxval,norm2,norm1 );

  VecMin( F, &loc, &minval );
  VecMax( F, &loc, &maxval );
  VecNorm( F, NORM_2, &norm2 );
  VecNorm( F, NORM_1, &norm1 );
  sprintf(fieldname,"FX");
  PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e\n", fieldname,minval,maxval,norm2,norm1 );

  ierr = DMCompositeRestoreAccess(pack,F,&Fu,&Fp);CHKERRQ(ierr);
//  ierr = VecDestroy(&F);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinViewBasicStokesSolutionResiduals(pTatinCtx ctx,SNES snes,DM pack,Vec X)
{
  PetscViewer  tmp;
  PetscErrorCode ierr;

  tmp = ctx->log;
  ctx->log  = PETSC_VIEWER_STDOUT_WORLD;
  PetscPrintf(PETSC_COMM_WORLD,"pTatinViewBasicStokesSolutionResiduals:\n");
  ierr = pTatinLogBasicStokesSolutionResiduals(ctx,snes,pack,X);CHKERRQ(ierr);
  ctx->log = tmp;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogBasicDMDA(pTatinCtx ctx,const char dmname[],DM dm)
{
  PetscReal min[3],max[3];
  PetscInt M,N,P,m,n,p,k;
  const PetscInt *lx,*ly,*lz;
  Vec coords;
  const char *prefix;
  PetscErrorCode ierr;

  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  ierr = DMGetOptionsPrefix(dm,&prefix);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dm,0,&M,&N,&P,&m,&n,&p, 0,0, 0,0,0, 0);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRanges(dm,&lx,&ly,&lz);CHKERRQ(ierr);

  PetscViewerASCIIPrintf(ctx->log,"  DMDA: (%18.18s)[prefix %s]\n    node [ %1.4d x %1.4d x %1.4d ]\n", dmname,prefix,M,N,P);

  coords = NULL;
  ierr = DMGetCoordinates(dm,&coords);CHKERRQ(ierr);
  if (coords) {
    ierr = DMDAGetBoundingBox(dm,min,max);CHKERRQ(ierr);
    PetscViewerASCIIPrintf(ctx->log,"    span [ %1.4e , %1.4e ] x [ %1.4e , %1.4e ] x [ %1.4e , %1.4e ]\n", min[0],max[0],min[1],max[1],min[2],max[2]);
  }
  PetscViewerASCIIPrintf(ctx->log,"    processor decomp [ %1.4d x %1.4d x %1.4d ]\n", m,n,p);

  PetscViewerASCIIPrintf(ctx->log,"    node decomp : i [");
  for (k=0; k<m-1; k++) {
    PetscViewerASCIIPrintf(ctx->log," %1.4d x",lx[k]);
  } PetscViewerASCIIPrintf(ctx->log," %1.4d ]\n",lx[m-1]);

  PetscViewerASCIIPrintf(ctx->log,"    node decomp : j [");
  for (k=0; k<n-1; k++) {
    PetscViewerASCIIPrintf(ctx->log," %1.4d x",ly[k]);
  } PetscViewerASCIIPrintf(ctx->log," %1.4d ]\n",ly[n-1]);

  PetscViewerASCIIPrintf(ctx->log,"    node decomp : k [");
  for (k=0; k<p-1; k++) {
    PetscViewerASCIIPrintf(ctx->log," %1.4d x",lz[k]);
  } PetscViewerASCIIPrintf(ctx->log," %1.4d ]\n",lz[p-1]);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogBasicMaterialPoints(pTatinCtx ctx,const char mpname[],DataBucket db)
{
  int npoints,buffer,allocated;

  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  DataBucketGetSizes(db,&npoints,&buffer,&allocated);
  PetscViewerASCIIPrintf(ctx->log,"  MaterialPoints: (%8.8s)  current %1.4d;  buffer size %1.4d; total allocated %1.4d\n", mpname,npoints,buffer,allocated);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogBasicCPUtime(pTatinCtx ctx,const char component_description[],double time)
{
  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  PetscViewerASCIIPrintf(ctx->log,"  CPU time (%s):  %1.4e (sec);\n", component_description,time);
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogNote(pTatinCtx ctx,const char comment[])
{
  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  PetscViewerASCIIPrintf(ctx->log,"  Note: %s\n",comment);
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogNote2(pTatinCtx ctx,const char comment1[],const char comment2[])
{
  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  PetscViewerASCIIPrintf(ctx->log,"  Note: %s %s\n",comment1,comment2);
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinLogPetscLog(pTatinCtx ctx,const char comment[])
{
  PetscErrorCode ierr;

  if (!ctx->log) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"pTatinCtx->log is NULL");
  PetscViewerASCIIPrintf(ctx->log,">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  if (comment) {
    PetscViewerASCIIPrintf(ctx->log,">>>>>>>>>>>  %s\n",comment);
  }
  ierr = PetscLogView(ctx->log);CHKERRQ(ierr);
  PetscViewerASCIIPrintf(ctx->log,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
  PetscFunctionReturn(0);
}

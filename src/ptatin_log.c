
#include "petsc.h"
#include "petscksp.h"
#include "petscsnes.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"



#undef __FUNCT__  
#define __FUNCT__ "pTatinLogOpenFile"
PetscErrorCode pTatinLogOpenFile(pTatinCtx ctx)
{
	PetscBool stdout = PETSC_FALSE;
	char name[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;

	sprintf(name,"%s/ptatin.log",ctx->outputpath);
	
	ierr = PetscOptionsGetBool(PETSC_NULL,"-ptatin_log_stdout",&stdout,PETSC_NULL);CHKERRQ(ierr);
	if (!stdout) {
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,name,&ctx->log);CHKERRQ(ierr);
	}
	else {
		ctx->log = PETSC_VIEWER_STDOUT_WORLD;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogCloseFile"
PetscErrorCode pTatinLogCloseFile(pTatinCtx ctx)
{
	PetscBool stdout = PETSC_FALSE;
	PetscErrorCode ierr;

	ierr = PetscOptionsGetBool(PETSC_NULL,"-ptatin_log_stdout",&stdout,PETSC_NULL);CHKERRQ(ierr);
	if (!stdout) {
		ierr = PetscViewerDestroy(&ctx->log);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogHeader"
PetscErrorCode pTatinLogHeader(pTatinCtx ctx)
{
	char username[PETSC_MAX_PATH_LEN];
	char date[PETSC_MAX_PATH_LEN];
	char machine[PETSC_MAX_PATH_LEN];
	char prgname[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	
	/* write header into options file */
	ierr = PetscGetUserName(username,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
	ierr = PetscGetDate(date,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
	ierr = PetscGetHostName(machine,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
	ierr = PetscGetProgramName(prgname,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"## =============================================================== \n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"##\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"##   pTatin Log File\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"##\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"##   Generated by user: %s\n",username);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"##   Date             : %s\n",date);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"##   Machine          : %s\n",machine);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"##   Using driver     : %s\n",prgname);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"##\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"## =============================================================== \n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(ctx->log,"\n");CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogBasic"
PetscErrorCode pTatinLogBasic(pTatinCtx ctx)
{
	PetscErrorCode ierr;

	PetscViewerASCIIPrintf(ctx->log,"------------------------------------------------------------------------------------------\n");
	PetscViewerASCIIPrintf(ctx->log,"  time step %6d \n", ctx->step);
	PetscViewerASCIIPrintf(ctx->log,"  time %1.4e;  dt %1.4e;  dt_min %1.4e;  dt_max %1.4e \n", ctx->time, ctx->dt, ctx->dt_min, ctx->dt_max);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogBasicKSP"
PetscErrorCode pTatinLogBasicKSP(pTatinCtx ctx,const char kspname[],KSP ksp)
{
	PetscReal rnorm;
	PetscInt its;
	const char *prefix;
	KSPConvergedReason reason;
	PetscErrorCode ierr;

	ierr = KSPGetOptionsPrefix(ksp,&prefix);CHKERRQ(ierr);
	ierr = KSPGetResidualNorm(ksp,&rnorm);CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);

	PetscViewerASCIIPrintf(ctx->log,"  KSP: (%8.8s)[prefix %8.8s] residual %1.4e;  iterations %1.4d;  reason %s;  \n", kspname,prefix,rnorm,its,KSPConvergedReasons[reason]);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogBasicSNES"
PetscErrorCode pTatinLogBasicSNES(pTatinCtx ctx,const char snesname[],SNES snes)
{
	PetscReal rnorm;
	PetscInt its,lits,nkspfails,nFevals,nstepfails;
	const char *prefix;
	SNESConvergedReason reason;
	PetscErrorCode ierr;
	
	ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);
	ierr = SNESGetFunctionNorm(snes,&rnorm);CHKERRQ(ierr);
	ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
	ierr = SNESGetLinearSolveIterations(snes,&lits);CHKERRQ(ierr);
	ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);

	PetscViewerASCIIPrintf(ctx->log,"  SNES: (%8.8s)[prefix %8.8s] residual %1.4e;  iterations %1.4d;  total linear its. %1.4d;  reason %s;  \n", snesname,prefix,rnorm,its,lits,SNESConvergedReasons[reason]);

	ierr = SNESGetLinearSolveFailures(snes,&nkspfails);CHKERRQ(ierr);
	ierr = SNESGetNonlinearStepFailures(snes,&nstepfails);CHKERRQ(ierr);
	ierr = SNESGetNumberFunctionEvals(snes,&nFevals);CHKERRQ(ierr);
	PetscViewerASCIIPrintf(ctx->log,"        function evals %1.4d;  ksp failures %1.4d;  step failures %1.4d\n", nFevals, nkspfails, nstepfails);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogBasicStokesSolution"
PetscErrorCode pTatinLogBasicStokesSolution(pTatinCtx ctx,DM pack,Vec X)
{
	char fieldname[120];
	Vec velocity,pressure;
	PetscReal minval,maxval,norm2,norm1;
	PetscInt loc;
	PetscErrorCode ierr;
	
	ierr = DMCompositeGetAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);

	PetscViewerASCIIPrintf( ctx->log, "  Field         min          max          norm_2       norm_1 \n" );
	VecMin( velocity, &loc, &minval );
	VecMax( velocity, &loc, &maxval );
	VecNorm( velocity, NORM_2, &norm2 );
	VecNorm( velocity, NORM_1, &norm1 );

	sprintf(fieldname,"velocity");
	PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e \n", fieldname,minval,maxval,norm2,norm1 );
	
	VecMin( pressure, &loc, &minval );
	VecMax( pressure, &loc, &maxval );
	VecNorm( pressure, NORM_2, &norm2 );
	VecNorm( pressure, NORM_1, &norm1 );
	sprintf(fieldname,"pressure");
	PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e \n", fieldname,minval,maxval,norm2,norm1 );
	
	VecMin( X, &loc, &minval );
	VecMax( X, &loc, &maxval );
	VecNorm( X, NORM_2, &norm2 );
	VecNorm( X, NORM_1, &norm1 );
	sprintf(fieldname,"X");
	PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e \n", fieldname,minval,maxval,norm2,norm1 );
	
	ierr = DMCompositeRestoreAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogBasicStokesSolutionResiduals"
PetscErrorCode pTatinLogBasicStokesSolutionResiduals(pTatinCtx ctx,SNES snes,DM pack,Vec X)
{
	char fieldname[120];
	Vec F,Fu,Fp;
	PetscReal minval,maxval,norm2,norm1;
	PetscInt loc;
	PetscErrorCode ierr;
	
//	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
//	ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
	
	ierr = SNESGetFunction(snes,&F,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

	ierr = DMCompositeGetAccess(pack,F,&Fu,&Fp);CHKERRQ(ierr);
	
	PetscViewerASCIIPrintf( ctx->log, "  Field         min          max          norm_2       norm_1 \n" );
	VecMin( Fu, &loc, &minval );
	VecMax( Fu, &loc, &maxval );
	VecNorm( Fu, NORM_2, &norm2 );
	VecNorm( Fu, NORM_1, &norm1 );
	
	sprintf(fieldname,"Fu");
	PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e \n", fieldname,minval,maxval,norm2,norm1 );
	
	VecMin( Fp, &loc, &minval );
	VecMax( Fp, &loc, &maxval );
	VecNorm( Fp, NORM_2, &norm2 );
	VecNorm( Fp, NORM_1, &norm1 );
	sprintf(fieldname,"Fp");
	PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e \n", fieldname,minval,maxval,norm2,norm1 );
	
	VecMin( F, &loc, &minval );
	VecMax( F, &loc, &maxval );
	VecNorm( F, NORM_2, &norm2 );
	VecNorm( F, NORM_1, &norm1 );
	sprintf(fieldname,"FX");
	PetscViewerASCIIPrintf( ctx->log, "  %8.8s      %+1.4e  %+1.4e  %+1.4e  %+1.4e \n", fieldname,minval,maxval,norm2,norm1 );
	
	ierr = DMCompositeRestoreAccess(pack,F,&Fu,&Fp);CHKERRQ(ierr);
//	ierr = VecDestroy(&F);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogBasicDMDA"
PetscErrorCode pTatinLogBasicDMDA(pTatinCtx ctx,const char dmname[],DM dm)
{
	const char *prefix;
	PetscReal min[3],max[3];
	PetscInt M,N,P;
	PetscErrorCode ierr;
	
	//ierr = DMGetOptionsPrefix(dm,&prefix);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dm,0,&M,&N,&P,0,0,0, 0,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(dm,min,max);CHKERRQ(ierr);
	
	PetscViewerASCIIPrintf(ctx->log,"  DMDA: (%8.8s)[prefix %8.8s]  node [ %1.4d x %1.4d x %1.4d ] \n", dmname,PETSC_NULL,M,N,P);
	PetscViewerASCIIPrintf(ctx->log,"                                     span [ %1.4e , %1.4e ] x [ %1.4e , %1.4e ] x [ %1.4e , %1.4e ] \n", min[0],max[0],min[1],max[1],min[2],max[2]);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogBasicMaterialPoints"
PetscErrorCode pTatinLogBasicMaterialPoints(pTatinCtx ctx,const char mpname[],DataBucket db)
{
	int npoints,buffer,allocated;
	
	DataBucketGetSizes(db,&npoints,&buffer,&allocated);
	PetscViewerASCIIPrintf(ctx->log,"  MaterialPoints: (%8.8s)  current %1.4d;  buffer size %1.4d; total allocated %1.4d \n", mpname,npoints,buffer,allocated);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinLogBasicCPUtime"
PetscErrorCode pTatinLogBasicCPUtime(pTatinCtx ctx,const char component_description[],double time)
{
	PetscViewerASCIIPrintf(ctx->log,"  CPU time (%s):  %1.4e (sec);\n", component_description,time);
	
	PetscFunctionReturn(0);
}



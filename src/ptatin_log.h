
#ifndef __ptatin_log_h__
#define __ptatin_log_h__

PetscErrorCode pTatinLogOpenFile(pTatinCtx ctx);
PetscErrorCode pTatinLogCloseFile(pTatinCtx ctx);
PetscErrorCode pTatinLogHeader(pTatinCtx ctx);

PetscErrorCode pTatinLogBasic(pTatinCtx ctx);
PetscErrorCode pTatinLogBasicKSP(pTatinCtx ctx,const char kspname[],KSP ksp);
PetscErrorCode pTatinLogBasicSNES(pTatinCtx ctx,const char snesname[],SNES snes);
PetscErrorCode pTatinLogBasicStokesSolution(pTatinCtx ctx,DM pack,Vec X);
PetscErrorCode pTatinLogBasicStokesSolutionResiduals(pTatinCtx ctx,SNES snes,DM pack,Vec X);
PetscErrorCode pTatinLogBasicDMDA(pTatinCtx ctx,const char dmname[],DM dm);
PetscErrorCode pTatinLogBasicMaterialPoints(pTatinCtx ctx,const char mpname[],DataBucket db);
PetscErrorCode pTatinLogBasicCPUtime(pTatinCtx ctx,const char component_description[],double time);


#endif

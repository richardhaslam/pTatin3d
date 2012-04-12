
#ifndef __pTatin_models_h__
#define __pTatin_models_h__

#include "petsc.h"
#include "petscvec.h"
#include "ptatin3d.h"
#include "dmda_bcs.h"

typedef enum {
	PTATIN_MODEL_INIT=0,
	PTATIN_MODEL_APPLY_BC,
	PTATIN_MODEL_APPLY_BCMG,
	PTATIN_MODEL_APPLY_MAT_BC,
	PTATIN_MODEL_APPLY_INIT_SOLUTION,
	PTATIN_MODEL_APPLY_INIT_MESH_GEOM,
	PTATIN_MODEL_APPLY_INIT_MAT_GEOM,
	PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,
	PTATIN_MODEL_OUTPUT,
	PTATIN_MODEL_DESTROY
} pTatinModelOperation;



extern pTatinModel *registered_model_list;

PetscErrorCode pTatinModelCreate(pTatinModel *model);
PetscErrorCode pTatinModelSetName(pTatinModel model,const char name[]);
PetscErrorCode pTatinModelSetUserData(pTatinModel model,void *data);
PetscErrorCode pTatinModelSetFunctionPointer(pTatinModel model,pTatinModelOperation type,void(*func)(void));
PetscErrorCode pTatinModelRegister(pTatinModel model);
PetscErrorCode pTatinModelGetByName(const char name[],pTatinModel *model);

PetscErrorCode pTatinModelGetUserData(pTatinModel model,void **data);

PetscErrorCode pTatinModelRegisterAll(void);


PetscErrorCode pTatinModel_Initialize(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_Destroy(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_Output(pTatinModel model,pTatinCtx ctx,Vec X,const char name[]);
PetscErrorCode pTatinModel_UpdateMeshGeometry(pTatinModel model,pTatinCtx ctx,Vec X);
PetscErrorCode pTatinModel_ApplyInitialMaterialGeometry(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyInitialMeshGeometry(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_Initialize(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyBoundaryCondition(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyBoundaryConditionMG(PetscInt nl,BCList bclist[],DM dav[],pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyMaterialBoundaryCondition(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyInitialSolution(pTatinModel model,pTatinCtx ctx,Vec X);


#endif

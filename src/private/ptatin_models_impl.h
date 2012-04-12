
#ifndef __private_ptatin_models_impl_h__
#define __private_ptatin_models_impl_h__

#include "petsc.h"
#include "ptatin3d.h"

struct _p_pTatinModel {
	pTatinCtx ptat_ctx;
	char *model_name;
	void *model_data;
	PetscErrorCode (*FP_pTatinModel_Initialize)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyInitialSolution)(pTatinCtx,Vec,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyBoundaryCondition)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyBoundaryConditionMG)(PetscInt,BCList*,DM*,pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyMaterialBoundaryCondition)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyInitialMeshGeometry)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyInitialMaterialGeometry)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_UpdateMeshGeometry)(pTatinCtx,Vec,void*);
	PetscErrorCode (*FP_pTatinModel_Output)(pTatinCtx,Vec,const char*,void*);
	PetscErrorCode (*FP_pTatinModel_Destroy)(pTatinCtx,void*);
};

#endif

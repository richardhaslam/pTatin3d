

#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "ptatin_models.h"

#include "template_modelctx.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Template"
PetscErrorCode ModelInitialize_Template(pTatinCtx c,void *ctx)
{
	TemplateModelCtx *data = (TemplateModelCtx*)ctx;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_template_param1",&data->param1,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-model_template_param2",&data->param2,&flg);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Template"
PetscErrorCode ModelApplyBoundaryCondition_Template(pTatinCtx c,void *ctx)
{
	TemplateModelCtx *data = (TemplateModelCtx*)ctx;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"param1 = %lf \n", data->param1 );
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Template"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Template(pTatinCtx c,void *ctx)
{
	TemplateModelCtx *data = (TemplateModelCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Template"
PetscErrorCode ModelApplyInitialMeshGeometry_Template(pTatinCtx c,void *ctx)
{
	TemplateModelCtx *data = (TemplateModelCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"param2 = %d \n", data->param2 );

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Template"
PetscErrorCode ModelApplyInitialMaterialGeometry_Template(pTatinCtx c,void *ctx)
{
	TemplateModelCtx *data = (TemplateModelCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Template"
PetscErrorCode ModelApplyUpdateMeshGeometry_Template(pTatinCtx c,void *ctx)
{
	TemplateModelCtx *data = (TemplateModelCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Template"
PetscErrorCode ModelOutput_Template(pTatinCtx c,void *ctx)
{
	TemplateModelCtx *data = (TemplateModelCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Template"
PetscErrorCode ModelDestroy_Template(pTatinCtx c,void *ctx)
{
	TemplateModelCtx *data = (TemplateModelCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Template"
PetscErrorCode pTatinModelRegister_Template(void)
{
	TemplateModelCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(TemplateModelCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(TemplateModelCtx));CHKERRQ(ierr);
	
	/* set initial values for model parameters */
	data->param1 = 0.0;
	data->param2 = 0;
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"template");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Template);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Template);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Template);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Template);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Template);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Template);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Template);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Template);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

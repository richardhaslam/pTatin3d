
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_models.h"

pTatinModel *registered_model_list = NULL;



#undef __FUNCT__
#define __FUNCT__ "pTatinModelCreate"
PetscErrorCode pTatinModelCreate(pTatinModel *model)
{
	PetscErrorCode ierr;
	pTatinModel m;
	PetscFunctionBegin;
	
	PetscMalloc(sizeof(struct _p_pTatinModel),&m);
	PetscMemzero(m,sizeof(struct _p_pTatinModel));
	*model = m;
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "pTatinModelSetName"
PetscErrorCode pTatinModelSetName(pTatinModel model,const char name[])
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	asprintf(&model->model_name,"%s",name);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelSetUserData"
PetscErrorCode pTatinModelSetUserData(pTatinModel model,void *data)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	model->model_data = data;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelGetUserData"
PetscErrorCode pTatinModelGetUserData(pTatinModel model,void **data)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	if (data) {
		*data = model->model_data;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinModelGetModelData"
PetscErrorCode pTatinModelGetModelData(pTatinModel ctx,const char name[],void **data)
{
	PetscErrorCode ierr;
	PetscContainer container;
	
  PetscFunctionBegin;
	ierr = PetscObjectQuery((PetscObject)ctx->data,name,(PetscObject*)&container);CHKERRQ(ierr);
	if (!container) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"No data with name \"%s\" was composed with model->model_data",name);
	ierr = PetscContainerGetPointer(container,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinModelSetModelData"
PetscErrorCode pTatinModelSetModelData(pTatinModel ctx,const char name[],void *data)
{
	PetscContainer container;
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
  ierr = PetscContainerCreate(PETSC_COMM_WORLD,&container);CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(container,(void*)data);CHKERRQ(ierr);
	
	ierr = PetscObjectCompose((PetscObject)ctx->data,name,(PetscObject)container);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelSetFunctionPointer"
PetscErrorCode pTatinModelSetFunctionPointer(pTatinModel model,pTatinModelOperation type,void(*func)(void))
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	switch (type) {
		case PTATIN_MODEL_INIT:
			model->FP_pTatinModel_Initialize = ( PetscErrorCode(*)(pTatinCtx,void*) )func;
			break;
		case PTATIN_MODEL_APPLY_BC:
			model->FP_pTatinModel_ApplyBoundaryCondition = ( PetscErrorCode(*)(pTatinCtx,void*) )func;
			break;
		case PTATIN_MODEL_APPLY_BCMG:
			model->FP_pTatinModel_ApplyBoundaryConditionMG = ( PetscErrorCode(*)(PetscInt,BCList*,DM*,pTatinCtx,void*) )func;
			break;
		case PTATIN_MODEL_APPLY_MAT_BC:
			model->FP_pTatinModel_ApplyMaterialBoundaryCondition = ( PetscErrorCode(*)(pTatinCtx,void*) )func;
			break;
		case PTATIN_MODEL_APPLY_INIT_SOLUTION:
			model->FP_pTatinModel_ApplyInitialSolution = ( PetscErrorCode(*)(pTatinCtx,Vec,void*) )func;
			break;
		case PTATIN_MODEL_APPLY_INIT_MESH_GEOM:
			model->FP_pTatinModel_ApplyInitialMeshGeometry = ( PetscErrorCode(*)(pTatinCtx,void*) )func;
			break;
		case PTATIN_MODEL_APPLY_INIT_MAT_GEOM:
			model->FP_pTatinModel_ApplyInitialMaterialGeometry = ( PetscErrorCode(*)(pTatinCtx,void*) )func;
			break;
		case PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM:
			model->FP_pTatinModel_UpdateMeshGeometry = ( PetscErrorCode(*)(pTatinCtx,Vec,void*) )func;
			break;
		case PTATIN_MODEL_OUTPUT:
			model->FP_pTatinModel_Output = ( PetscErrorCode(*)(pTatinCtx,Vec,const char*,void*) )func;
			break;
		case PTATIN_MODEL_DESTROY:
			model->FP_pTatinModel_Destroy = ( PetscErrorCode(*)(pTatinCtx,void*) )func;
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"  [pTatinModel]: unknown operation specified");
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ptatin_match_model_index"
PetscErrorCode ptatin_match_model_index(const char modelname[],PetscInt *index)
{
	PetscErrorCode ierr;
	PetscBool match;
	pTatinModel item;	
	int cnt;
	
	*index = -1;
	cnt = 0;
	item = registered_model_list[0];
	while (item!=PETSC_NULL) {
		match = PETSC_FALSE;
		PetscStrcmp(modelname,item->model_name,&match);
		if (match) {
			*index = cnt;
			break;
		}
		cnt++;
		item = registered_model_list[cnt];
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister"
PetscErrorCode pTatinModelRegister(pTatinModel model)
{
	PetscInt cnt,list_length,index;
	pTatinModel ii;
	pTatinModel *tmp;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (registered_model_list==PETSC_NULL) {
		registered_model_list = malloc( sizeof(pTatinModel) );
		registered_model_list[0] = PETSC_NULL;
	} 
	
	/* find list size */
	cnt = 0;
	ii = registered_model_list[0];
	while (ii!=PETSC_NULL) {
		cnt++;
		ii = registered_model_list[cnt];
	}
	list_length = cnt;
	
	/* check model not already loaded with same name */
	ierr = ptatin_match_model_index(model->model_name,&index);CHKERRQ(ierr);
	if (index!=-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"  [pTatinModel]: Model with name \"%s\" has already been registered",model->model_name );
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"  [pTatinModel]: Registering model [%d] with name \"%s\"\n",list_length, model->model_name );		
	}
	
	tmp = realloc( registered_model_list,sizeof(pTatinModel)*(list_length+2) );
	registered_model_list = tmp;
	
	registered_model_list[list_length  ] = model;
	registered_model_list[list_length+1] = PETSC_NULL;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelGetByName"
PetscErrorCode pTatinModelGetByName(const char name[],pTatinModel *model)
{	
	PetscErrorCode ierr;
	PetscInt index;
	PetscFunctionBegin;
	
	*model = PETSC_NULL;
	ierr = ptatin_match_model_index(name,&index);CHKERRQ(ierr);
	if ( index==-1 ) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"  [pTatinModel]: -ptatin_model \"%s\" wasn't identified in list registered_model_list[]",name );
	}
	(*model) = registered_model_list[index];
	
	PetscFunctionReturn(0);
}

/* wrappers */
/*
PetscErrorCode (*FP_pTatin2d_ModelOutput)(void*);
PetscErrorCode (*FP_pTatin2d_ModelDestroy)(void*);
*/

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_Destroy"
PetscErrorCode pTatinModel_Destroy(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_Destroy) {
		ierr = model->FP_pTatinModel_Destroy(ctx,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_Output"
PetscErrorCode pTatinModel_Output(pTatinModel model,pTatinCtx ctx,Vec X,const char name[])
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_Output) {
		ierr = model->FP_pTatinModel_Output(ctx,X,name,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyInitialSolution"
PetscErrorCode pTatinModel_ApplyInitialSolution(pTatinModel model,pTatinCtx ctx,Vec X)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_ApplyInitialSolution) {
		ierr = model->FP_pTatinModel_ApplyInitialSolution(ctx,X,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_UpdateMeshGeometry"
PetscErrorCode pTatinModel_UpdateMeshGeometry(pTatinModel model,pTatinCtx ctx,Vec X)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_UpdateMeshGeometry) {
		ierr = model->FP_pTatinModel_UpdateMeshGeometry(ctx,X,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyInitialMaterialGeometry"
PetscErrorCode pTatinModel_ApplyInitialMaterialGeometry(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_ApplyInitialMaterialGeometry) {
		ierr = model->FP_pTatinModel_ApplyInitialMaterialGeometry(ctx,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyInitialMeshGeometry"
PetscErrorCode pTatinModel_ApplyInitialMeshGeometry(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_ApplyInitialMeshGeometry) {
		ierr = model->FP_pTatinModel_ApplyInitialMeshGeometry(ctx,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_Initialize"
PetscErrorCode pTatinModel_Initialize(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	if (model->FP_pTatinModel_Initialize) {
		ierr = model->FP_pTatinModel_Initialize(ctx,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyBoundaryCondition"
PetscErrorCode pTatinModel_ApplyBoundaryCondition(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_ApplyBoundaryCondition) {
		ierr = model->FP_pTatinModel_ApplyBoundaryCondition(ctx,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyBoundaryConditionMG"
PetscErrorCode pTatinModel_ApplyBoundaryConditionMG(PetscInt nl,BCList bclist[],DM dav[],pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_ApplyBoundaryConditionMG) {
		ierr = model->FP_pTatinModel_ApplyBoundaryConditionMG(nl,bclist,dav,ctx,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyMaterialBoundaryCondition"
PetscErrorCode pTatinModel_ApplyMaterialBoundaryCondition(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (model->FP_pTatinModel_ApplyMaterialBoundaryCondition) {
		ierr = model->FP_pTatinModel_ApplyMaterialBoundaryCondition(ctx,model->model_data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}


/* Users add prototypes here */
extern PetscErrorCode pTatinModelRegister_Template(void);
extern PetscErrorCode pTatinModelRegister_ViscousSinker(void);
extern PetscErrorCode pTatinModelRegister_GENE3D(void);
extern PetscErrorCode pTatinModelRegister_Indentor(void);

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegisterAll"
PetscErrorCode pTatinModelRegisterAll(void)
{
	PetscErrorCode ierr;

	PetscFunctionBegin;
	
	/* call registration functions for each model here */
	ierr = pTatinModelRegister_Template();CHKERRQ(ierr);
	ierr = pTatinModelRegister_ViscousSinker();CHKERRQ(ierr);
	ierr = pTatinModelRegister_GENE3D();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Indentor();CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}

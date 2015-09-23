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
 **    filename:   ptatin_models.c
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

#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_models.h"

pTatinModel *registered_model_list = NULL;

PetscLogEvent PTATIN_ModelInitialize;
PetscLogEvent PTATIN_ModelApplyInitialSolution;
PetscLogEvent PTATIN_ModelApplyInitialMeshGeometry;
PetscLogEvent PTATIN_ModelApplyInitialMaterialGeometry;
PetscLogEvent PTATIN_ModelApplyInitialStokesVariableMarkers;
PetscLogEvent PTATIN_ModelApplyBoundaryCondition;
PetscLogEvent PTATIN_ModelApplyBoundaryConditionMG;
PetscLogEvent PTATIN_ModelApplyMaterialBoundaryCondition;
PetscLogEvent PTATIN_ModelUpdateMeshGeometry;
PetscLogEvent PTATIN_ModelOutput;


#undef __FUNCT__
#define __FUNCT__ "pTatinModelCreate"
PetscErrorCode pTatinModelCreate(pTatinModel *model)
{
	PetscErrorCode ierr;
	pTatinModel m;
	PetscFunctionBegin;
	
	ierr = PetscMalloc(sizeof(struct _p_pTatinModel),&m);CHKERRQ(ierr);
	ierr = PetscMemzero(m,sizeof(struct _p_pTatinModel));CHKERRQ(ierr);

	m->disable_initial_solution          = PETSC_FALSE;
	m->disable_initial_stokes_variables  = PETSC_FALSE;
	m->disable_apply_bc                  = PETSC_FALSE;
	m->disable_apply_bc_mg               = PETSC_FALSE;
	m->disable_apply_material_bc         = PETSC_FALSE;
	m->disable_initial_mesh_geometry     = PETSC_FALSE;
	m->disable_initial_material_geometry = PETSC_FALSE;
	m->disable_update_mesh_geometry      = PETSC_FALSE;
	m->disable_output                    = PETSC_FALSE;

	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_intial_solution_disable",&m->disable_initial_solution,0);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_initial_stokes_variables_disable",&m->disable_initial_stokes_variables,0);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_apply_bc_disable",&m->disable_apply_bc,0);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_apply_bc_mg_disable",&m->disable_apply_bc_mg,0);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_apply_material_bc_disable",&m->disable_apply_material_bc,0);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_initial_mesh_geometry_disable",&m->disable_initial_mesh_geometry,0);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_initial_material_geometry_disable",&m->disable_initial_material_geometry,0);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_update_mesh_geometry_disable",&m->disable_update_mesh_geometry,0);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-ptatin_model_output_disable",&m->disable_output,0);CHKERRQ(ierr);
	
	*model = m;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelDestroy"
PetscErrorCode pTatinModelDestroy(pTatinModel *m)
{
  pTatinModel    model;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  if (!m) PetscFunctionReturn(0);
  model = *m;

  if (model->model_name) { free(model->model_name); }
  ierr = pTatinModel_Destroy(model,NULL);CHKERRQ(ierr);
  PetscFree(model);
  *m = NULL;
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelSetName"
PetscErrorCode pTatinModelSetName(pTatinModel model,const char name[])
{
	PetscFunctionBegin;
	asprintf(&model->model_name,"%s",name);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelSetUserData"
PetscErrorCode pTatinModelSetUserData(pTatinModel model,void *data)
{
	PetscFunctionBegin;
	model->model_data = data;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelGetUserData"
PetscErrorCode pTatinModelGetUserData(pTatinModel model,void **data)
{
	PetscFunctionBegin;
	if (data) {
		*data = model->model_data;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelGetName"
PetscErrorCode pTatinModelGetName(pTatinModel model,char **name)
{
	PetscFunctionBegin;
	if (name) {
		*name = model->model_name;
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
		case PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS:
			model->FP_pTatinModel_ApplyInitialStokesVariableMarkers = ( PetscErrorCode(*)(pTatinCtx,Vec,void*) )func;
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
	PetscBool match;
	pTatinModel item;	
	int cnt;
	
	*index = -1;
	cnt = 0;
	item = registered_model_list[0];
	while (item!=NULL) {
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
	
	if (registered_model_list==NULL) {
		registered_model_list = malloc( sizeof(pTatinModel) );
		registered_model_list[0] = NULL;
	} 
	
	/* find list size */
	cnt = 0;
	ii = registered_model_list[0];
	while (ii!=NULL) {
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
	registered_model_list[list_length+1] = NULL;
	
	if (model->disable_output) {
		PetscPrintf(PETSC_COMM_WORLD,"  [pTatinModel]: Output functionality for \"%s\" has been deactivated\n",model->model_name);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelGetByName"
PetscErrorCode pTatinModelGetByName(const char name[],pTatinModel *model)
{	
	PetscErrorCode ierr;
	PetscInt index;
	PetscFunctionBegin;
	
	*model = NULL;
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
	
	ierr = PetscLogEventBegin(PTATIN_ModelOutput,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_output) {
		if (model->FP_pTatinModel_Output) {
			ierr = model->FP_pTatinModel_Output(ctx,X,name,model->model_data);CHKERRQ(ierr);
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelOutput,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyInitialSolution"
PetscErrorCode pTatinModel_ApplyInitialSolution(pTatinModel model,pTatinCtx ctx,Vec X)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_ModelApplyInitialSolution,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_initial_solution) {
		if (model->FP_pTatinModel_ApplyInitialSolution) {
			ierr = model->FP_pTatinModel_ApplyInitialSolution(ctx,X,model->model_data);CHKERRQ(ierr);
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelApplyInitialSolution,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyInitialStokesVariableMarkers"
PetscErrorCode pTatinModel_ApplyInitialStokesVariableMarkers(pTatinModel model,pTatinCtx ctx,Vec X)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_ModelApplyInitialStokesVariableMarkers,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_initial_stokes_variables) {
		if (model->FP_pTatinModel_ApplyInitialStokesVariableMarkers) {
			ierr = model->FP_pTatinModel_ApplyInitialStokesVariableMarkers(ctx,X,model->model_data);CHKERRQ(ierr);
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelApplyInitialStokesVariableMarkers,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_UpdateMeshGeometry"
PetscErrorCode pTatinModel_UpdateMeshGeometry(pTatinModel model,pTatinCtx ctx,Vec X)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_ModelUpdateMeshGeometry,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_update_mesh_geometry) {
		if (model->FP_pTatinModel_UpdateMeshGeometry) {
			ierr = model->FP_pTatinModel_UpdateMeshGeometry(ctx,X,model->model_data);CHKERRQ(ierr);
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelUpdateMeshGeometry,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyInitialMaterialGeometry"
PetscErrorCode pTatinModel_ApplyInitialMaterialGeometry(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_ModelApplyInitialMaterialGeometry,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_initial_material_geometry) {
		if (model->FP_pTatinModel_ApplyInitialMaterialGeometry) {
			ierr = model->FP_pTatinModel_ApplyInitialMaterialGeometry(ctx,model->model_data);CHKERRQ(ierr);
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelApplyInitialMaterialGeometry,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyInitialMeshGeometry"
PetscErrorCode pTatinModel_ApplyInitialMeshGeometry(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_ModelApplyInitialMeshGeometry,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_initial_mesh_geometry) {
		if (model->FP_pTatinModel_ApplyInitialMeshGeometry) {
			ierr = model->FP_pTatinModel_ApplyInitialMeshGeometry(ctx,model->model_data);CHKERRQ(ierr);
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelApplyInitialMeshGeometry,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_Initialize"
PetscErrorCode pTatinModel_Initialize(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscLogEventBegin(PTATIN_ModelInitialize,0,0,0,0);CHKERRQ(ierr);
	if (model->FP_pTatinModel_Initialize) {
		ierr = model->FP_pTatinModel_Initialize(ctx,model->model_data);CHKERRQ(ierr);
	}
	ierr = PetscLogEventEnd(PTATIN_ModelInitialize,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyBoundaryCondition"
PetscErrorCode pTatinModel_ApplyBoundaryCondition(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_ModelApplyBoundaryCondition,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_apply_bc) {
		if (model->FP_pTatinModel_ApplyBoundaryCondition) {
			ierr = model->FP_pTatinModel_ApplyBoundaryCondition(ctx,model->model_data);CHKERRQ(ierr);
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelApplyBoundaryCondition,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyBoundaryConditionMG"
PetscErrorCode pTatinModel_ApplyBoundaryConditionMG(PetscInt nl,BCList bclist[],DM dav[],pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_ModelApplyBoundaryConditionMG,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_apply_bc_mg) {
		if (model->FP_pTatinModel_ApplyBoundaryConditionMG) {
			ierr = model->FP_pTatinModel_ApplyBoundaryConditionMG(nl,bclist,dav,ctx,model->model_data);CHKERRQ(ierr);
		} else {
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"  [pTatinModel]: -ptatin_model \"%s\" wasn't prodivided with the operation \"ApplyBoundaryConditionMG\"",model->model_name );
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelApplyBoundaryConditionMG,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModel_ApplyMaterialBoundaryCondition"
PetscErrorCode pTatinModel_ApplyMaterialBoundaryCondition(pTatinModel model,pTatinCtx ctx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscLogEventBegin(PTATIN_ModelApplyMaterialBoundaryCondition,0,0,0,0);CHKERRQ(ierr);
	if (!model->disable_apply_material_bc) {
		if (model->FP_pTatinModel_ApplyMaterialBoundaryCondition) {
			ierr = model->FP_pTatinModel_ApplyMaterialBoundaryCondition(ctx,model->model_data);CHKERRQ(ierr);
		}
	}
	ierr = PetscLogEventEnd(PTATIN_ModelApplyMaterialBoundaryCondition,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelDeRegisterAll"
PetscErrorCode pTatinModelDeRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscInt i;
  pTatinModel item;
  
  i = 0;
  item = registered_model_list[0];
  while (item) {
    ierr = pTatinModelDestroy(&item);CHKERRQ(ierr);
    i++;
    item = registered_model_list[i];
  }
  free(registered_model_list);
  
	PetscFunctionReturn(0);
}

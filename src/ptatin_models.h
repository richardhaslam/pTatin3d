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
 **    filename:   ptatin_models.h
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
    PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,
	PTATIN_MODEL_APPLY_INIT_MESH_GEOM,
	PTATIN_MODEL_APPLY_INIT_MAT_GEOM,
	PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,
	PTATIN_MODEL_OUTPUT,
	PTATIN_MODEL_DESTROY
} pTatinModelOperation;



extern pTatinModel *registered_model_list;

PetscErrorCode pTatinModelCreate(pTatinModel *model);
PetscErrorCode pTatinModelDestroy(pTatinModel *m);
PetscErrorCode pTatinModelSetName(pTatinModel model,const char name[]);
PetscErrorCode pTatinModelSetFunctionPointer(pTatinModel model,pTatinModelOperation type,void(*func)(void));
PetscErrorCode pTatinModelRegister(pTatinModel model);
PetscErrorCode pTatinModelGetByName(const char name[],pTatinModel *model);

PetscErrorCode pTatinModelSetUserData(pTatinModel model,void *data);
PetscErrorCode pTatinModelGetUserData(pTatinModel model,void **data);

PetscErrorCode pTatinModelGetModelData(pTatinModel ctx,const char name[],void **data);
PetscErrorCode pTatinModelSetModelData(pTatinModel ctx,const char name[],void *data);

PetscErrorCode pTatinModelGetName(pTatinModel model,char **name);

PetscErrorCode pTatinModelRegisterAll(void);
PetscErrorCode pTatinModelDeRegisterAll(void);

PetscErrorCode pTatinModel_Initialize(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_Destroy(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_Output(pTatinModel model,pTatinCtx ctx,Vec X,const char name[]);
PetscErrorCode pTatinModel_UpdateMeshGeometry(pTatinModel model,pTatinCtx ctx,Vec X);
PetscErrorCode pTatinModel_ApplyInitialMaterialGeometry(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyInitialMeshGeometry(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyBoundaryCondition(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyBoundaryConditionMG(PetscInt nl,BCList bclist[],DM dav[],pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyMaterialBoundaryCondition(pTatinModel model,pTatinCtx ctx);
PetscErrorCode pTatinModel_ApplyInitialSolution(pTatinModel model,pTatinCtx ctx,Vec X);
PetscErrorCode pTatinModel_ApplyInitialStokesVariableMarkers(pTatinModel model,pTatinCtx ctx,Vec X);

#endif


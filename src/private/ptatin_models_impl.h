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
 **    filename:   ptatin_models_impl.h
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

#ifndef __private_ptatin_models_impl_h__
#define __private_ptatin_models_impl_h__

#include <petsc.h>
#include "../ptatin3d.h"

struct _p_pTatinModel {
	pTatinCtx ptat_ctx;
	char *model_name;
	void *model_data;
	PetscContainer data;
	PetscErrorCode (*FP_pTatinModel_Initialize)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyInitialSolution)(pTatinCtx,Vec,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyInitialStokesVariableMarkers)(pTatinCtx,Vec,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyBoundaryCondition)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyBoundaryConditionMG)(PetscInt,BCList*,DM*,pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyMaterialBoundaryCondition)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyInitialMeshGeometry)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ApplyInitialMaterialGeometry)(pTatinCtx,void*);
	PetscErrorCode (*FP_pTatinModel_ExportInnerMesh)(pTatinCtx,DM, Vec,void*);
	PetscErrorCode (*FP_pTatinModel_UpdateMeshGeometry)(pTatinCtx,Vec,void*);
	PetscErrorCode (*FP_pTatinModel_Output)(pTatinCtx,Vec,const char*,void*);
	PetscErrorCode (*FP_pTatinModel_Destroy)(pTatinCtx,void*);
	PetscBool disable_initial_solution;
	PetscBool disable_initial_stokes_variables;
	PetscBool disable_apply_bc;
	PetscBool disable_apply_bc_mg;
	PetscBool disable_apply_material_bc;
	PetscBool disable_initial_mesh_geometry;
	PetscBool disable_initial_material_geometry;
	PetscBool disable_export_inner_mesh;
	PetscBool disable_update_mesh_geometry;
	PetscBool disable_output;
};

#endif

/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      model_ops_advdiff_example.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_models.h"

#include "dmda_iterator.h"
#include "dmda_view_petscvtk.h"
#include "energy_output.h"

#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "MPntPEnergy_def.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_AdvDiffExample"
PetscErrorCode ModelInitialize_AdvDiffExample(pTatinCtx c,void *ctx)
{
	PetscBool      flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_AdvDiffExample"
PetscErrorCode ModelApplyBoundaryCondition_AdvDiffExample(pTatinCtx c,void *ctx)
{
	PetscErrorCode ierr;
	PetscScalar    val;
	PetscBool      active_energy;
	PhysCompEnergy energy;
	DM             daT;
	BCList         bclist;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* ignore all velocity boundary conditions */
	
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"energy must be activated");
	}
	
	/* a) fetch the energy context */
	ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
	daT    = energy->daT;
	bclist = energy->T_bclist;
	
	val = 2.0;
	ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	
	//ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	//ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_AdvDiffExample"
PetscErrorCode ModelApplyBoundaryConditionMG_AdvDiffExample(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	PetscInt       n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_AdvDiffExample"
PetscErrorCode ModelApplyMaterialBoundaryCondition_AdvDiffExample(pTatinCtx c,void *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_AdvDiffExample"
PetscErrorCode ModelApplyInitialMeshGeometry_AdvDiffExample(pTatinCtx c,void *ctx)
{
	PetscErrorCode ierr;
	PhysCompStokes stokes;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	ierr = DMDASetUniformCoordinates(stokes->dav,0.0,1.0, 0.0,1.0, 0.0,1.0);CHKERRQ(ierr);
	
	/* note - Don't access the energy mesh here, its not yet created */
	/* note - The initial velocity mesh geometry will be copied into the energy mesh */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_AdvDiffExample"
PetscErrorCode ModelApplyInitialMaterialGeometry_AdvDiffExample(pTatinCtx c,void *ctx)
{
	PetscErrorCode ierr;
	DataBucket     materialpoint_db;
	DataField      PField_std,PField_stokes, PField_energy;
	int            p,n_mp_points;
	PetscBool      active_energy;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);


	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mp_points,0,0);
	
	/* define viscous params */
	DataBucketGetDataFieldByName(materialpoint_db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	
	DataBucketGetDataFieldByName(materialpoint_db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		double      *position;
		
		DataFieldAccessPoint(PField_std,    p, (void**)&material_point);
		DataFieldAccessPoint(PField_stokes, p, (void**)&mpprop_stokes);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,0);
		MPntPStokesSetField_eta_effective(mpprop_stokes,1.0);
		MPntPStokesSetField_density(mpprop_stokes,1.0);
	}
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	
	
	
	/* define thermal params */
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	
	if (active_energy) {
		DataBucketGetDataFieldByName(materialpoint_db,MPntStd_classname,&PField_std);
		DataFieldGetAccess(PField_std);

		DataBucketGetDataFieldByName(materialpoint_db,MPntPEnergy_classname,&PField_energy);
		DataFieldGetAccess(PField_energy);
		
		for (p=0; p<n_mp_points; p++) {
			MPntStd     *material_point;
			MPntPEnergy *mpprop_energy;
			double      *position;
			
			DataFieldAccessPoint(PField_std,   p, (void**)&material_point);
			DataFieldAccessPoint(PField_energy,p, (void**)&mpprop_energy);
			
			/* Access using the getter function provided for you (recommeneded for beginner user) */
			MPntStdGetField_global_coord(material_point,&position);
			
			/* user the setters provided for you */
			MPntPEnergySetField_diffusivity(mpprop_energy,1.0);
			MPntPEnergySetField_heat_source(mpprop_energy,0.0);
		}
		DataFieldRestoreAccess(PField_std);
		DataFieldRestoreAccess(PField_energy);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_AdvDiffExample"
PetscErrorCode ModelApplyInitialSolution_AdvDiffExample(pTatinCtx c,Vec X,void *ctx)
{
	PhysCompStokes stokes;
	PhysCompEnergy energy;
	PetscBool      active_energy;
	Vec            velocity,pressure,temperature;
	DM             daT,dav,dap,multipys_pack;
	PetscScalar    vx_const,vy_const,vz_const,p_init;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* note - X contains only the velocity and pressure */
	
	/* define velocity field for advection test */
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	multipys_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(multipys_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	vx_const = 0.0;
	vy_const = 0.0;
	vz_const = 0.0;
	ierr = DMDAVecTraverse3d(dav,velocity,0,DMDAVecTraverse3d_Constant,(void*)&vx_const);CHKERRQ(ierr);
	ierr = DMDAVecTraverse3d(dav,velocity,1,DMDAVecTraverse3d_Constant,(void*)&vy_const);CHKERRQ(ierr);
	ierr = DMDAVecTraverse3d(dav,velocity,2,DMDAVecTraverse3d_Constant,(void*)&vz_const);CHKERRQ(ierr);
	
	/* note - pressure has no coordinates set, so use ijk traversal */
	p_init = 0.0;
	ierr = DMDAVecTraverseIJK(dap,pressure,0,DMDAVecTraverseIJK_Constant,(void*)&p_init);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"energy must be activated");
	}
	

	/* define initial temperature field for advection test */
	/* a) fetch the energy context */
	ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
	/* b) get the DM */
	daT = energy->daT;
	/* c) fetch the temp vector */
	ierr = pTatinPhysCompGetData_Energy(c,&temperature,PETSC_NULL);CHKERRQ(ierr);


	ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_GaussianXY,PETSC_NULL);CHKERRQ(ierr);
	//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_GaussianXYZ,PETSC_NULL);CHKERRQ(ierr);
	//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepX,PETSC_NULL);CHKERRQ(ierr);
	//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepXY,PETSC_NULL);CHKERRQ(ierr);
	//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepXYZ,PETSC_NULL);CHKERRQ(ierr);
	ierr = VecShift(temperature,2.0);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_AdvDiffExample"
PetscErrorCode ModelApplyUpdateMeshGeometry_AdvDiffExample(pTatinCtx c,Vec X,void *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_AdvDiffExample"
PetscErrorCode ModelOutput_AdvDiffExample(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	PetscErrorCode ierr;
	DataBucket     materialpoint_db;
	PhysCompEnergy energy;
	Vec temperature;
	char name[256];
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* stokes variables */
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	
	/* energy variables */
	ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
	
	ierr = pTatinPhysCompGetData_Energy(c,&temperature,PETSC_NULL);CHKERRQ(ierr);

	/* standard viewer */
	ierr = pTatin3d_ModelOutput_Temperature_Energy(c,temperature,prefix);CHKERRQ(ierr);

	/* debugging adv-diff solver */
#if 0	
	sprintf(name,"%s/%s_T.vtk",c->outputpath,prefix);
	ierr = DMDAViewPetscVTK(energy->daT,temperature,name);CHKERRQ(ierr);
	
	sprintf(name,"%s/%s_Tlast.vtk",c->outputpath,prefix);
	ierr = DMDAViewPetscVTK(energy->daT,energy->Told,name);CHKERRQ(ierr);
	
	sprintf(name,"%s/%s_advdiff_u.vtk",c->outputpath,prefix);
	ierr = DMDAViewPetscVTK(energy->daT,energy->u_minus_V,name);CHKERRQ(ierr);
#endif
	
	/* stokes + energy material points */
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
	{
		const int                  nf = 2;
		const MaterialPointField   mp_prop_list[] = { MPField_Std, MPField_Energy }; 

		
		sprintf(name,"%s_advdiff_mpoints.vtk",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,name);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_AdvDiffExample"
PetscErrorCode ModelDestroy_AdvDiffExample(pTatinCtx c,void *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_AdvDiffExample"
PetscErrorCode pTatinModelRegister_AdvDiffExample(void)
{
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"advdiff_example");CHKERRQ(ierr);

	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_AdvDiffExample);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_AdvDiffExample);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

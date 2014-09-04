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
#include "output_material_points.h"

#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "MPntPEnergy_def.h"

PetscInt model_setup = 0;
PetscReal vel_scale       = 1.0;
PetscReal diffusion_scale = 1.0;

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_AdvDiffExample"
PetscErrorCode ModelInitialize_AdvDiffExample(pTatinCtx c,void *ctx)
{	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscOptionsGetInt(NULL,"-advdiff_setup",&model_setup,0);
	PetscOptionsGetReal(NULL,"-advdiff_vel_scale",&vel_scale,0);
	PetscOptionsGetReal(NULL,"-advdiff_diff_scale",&diffusion_scale,0);

	switch (model_setup) {
		case 0:
			PetscPrintf(PETSC_COMM_WORLD,"AdvDiff test: Pure diffusion of Gaussian hill.\n");
			break;
		case 1:
			PetscPrintf(PETSC_COMM_WORLD,"AdvDiff test: Advection of Gaussian hill in x-direction.\n");
			break;
		case 2:
			PetscPrintf(PETSC_COMM_WORLD,"AdvDiff test: Skew advection test of Hughes.\n");
			break;
		case 3:
			PetscPrintf(PETSC_COMM_WORLD,"AdvDiff test: Step advection in x-direction.\n");
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_HughesSkewTest"
PetscBool BCListEvaluator_HughesSkewTest( PetscScalar position[], PetscScalar *value, void *ctx ) 
{
	PetscBool impose_dirichlet = PETSC_TRUE;
	PetscScalar dv = *((PetscScalar*)ctx);

	if (position[1] < 0.2) {
		*value = dv;
	} else {
		*value = 0.0;
	}
	return impose_dirichlet;
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
	
	if (model_setup == 0) {
		val = 2.0;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	}	

	if (model_setup == 1) {
		val = 2.0;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
		// outflow
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	}

	if (model_setup == 2) {
		val = 1.0;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);

		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMIN_LOC,0,BCListEvaluator_HughesSkewTest,(void*)&val);CHKERRQ(ierr);
		// outflow IMAX,JMAX		
	}

	if (model_setup == 3) {
		/*
		val = 2.0;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
		//val = 1.0;
		//ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
		*/
		val = 2.0;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_AdvDiffExample"
PetscErrorCode ModelApplyBoundaryConditionMG_AdvDiffExample(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	PetscInt       n;
	
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

	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
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
		
		if (position[0] < 0.5) {
			MPntStdSetField_phase_index(material_point,1);
		}
		if (position[1] < 0.7) {
			MPntStdSetField_phase_index(material_point,2);
		}
		if (position[2] < 0.9) {
			MPntStdSetField_phase_index(material_point,3);
		}
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
			if (model_setup == 0) {
				MPntPEnergySetField_diffusivity(mpprop_energy,1.0);
				MPntPEnergySetField_heat_source(mpprop_energy,0.0);
			}
			if (model_setup == 1) {
				//MPntPEnergySetField_diffusivity(mpprop_energy,position[0]*position[1]+position[2]*position[2]+1.0);
				MPntPEnergySetField_diffusivity(mpprop_energy,1.0e-6);
				MPntPEnergySetField_heat_source(mpprop_energy,0.0);
			}
			if (model_setup == 2) {
				MPntPEnergySetField_diffusivity(mpprop_energy,1.0e-6);
				MPntPEnergySetField_heat_source(mpprop_energy,0.0);
			}
			if (model_setup == 3) {
				MPntPEnergySetField_diffusivity(mpprop_energy,1.0e-6*diffusion_scale);
				//MPntPEnergySetField_diffusivity(mpprop_energy,position[0]*position[1]+position[2]*position[2]+1.0);
				MPntPEnergySetField_heat_source(mpprop_energy,0.0);
			}
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
	
	if (model_setup == 0) {
		vx_const = 0.0;
		vy_const = 0.0;
		vz_const = 0.0;
		ierr = DMDAVecTraverse3d(dav,velocity,0,DMDAVecTraverse3d_Constant,(void*)&vx_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,1,DMDAVecTraverse3d_Constant,(void*)&vy_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,2,DMDAVecTraverse3d_Constant,(void*)&vz_const);CHKERRQ(ierr);
	}
	if (model_setup == 1) {
		vx_const = 1.0;
		vy_const = 0.0;
		vz_const = 0.0;
		ierr = DMDAVecTraverse3d(dav,velocity,0,DMDAVecTraverse3d_Constant,(void*)&vx_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,1,DMDAVecTraverse3d_Constant,(void*)&vy_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,2,DMDAVecTraverse3d_Constant,(void*)&vz_const);CHKERRQ(ierr);
	}
	if (model_setup == 2) {
		vx_const = 1.0;
		vy_const = 1.0;
		vz_const = 0.0;
		ierr = DMDAVecTraverse3d(dav,velocity,0,DMDAVecTraverse3d_Constant,(void*)&vx_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,1,DMDAVecTraverse3d_Constant,(void*)&vy_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,2,DMDAVecTraverse3d_Constant,(void*)&vz_const);CHKERRQ(ierr);
	}
	if (model_setup == 3) {
		/*
		vx_const = 1.0 * vel_scale;
		vy_const = 0.0;
		vz_const = 0.0;
		ierr = DMDAVecTraverse3d(dav,velocity,0,DMDAVecTraverse3d_Constant,(void*)&vx_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,1,DMDAVecTraverse3d_Constant,(void*)&vy_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,2,DMDAVecTraverse3d_Constant,(void*)&vz_const);CHKERRQ(ierr);
		 */
		vx_const = 0.0;
		vy_const = 1.0 * vel_scale;
		vz_const = 0.0;
		ierr = DMDAVecTraverse3d(dav,velocity,0,DMDAVecTraverse3d_Constant,(void*)&vx_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,1,DMDAVecTraverse3d_Constant,(void*)&vy_const);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dav,velocity,2,DMDAVecTraverse3d_Constant,(void*)&vz_const);CHKERRQ(ierr);
		
	}
	
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
	ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);

	if (model_setup == 0) {
		ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_GaussianXY,NULL);CHKERRQ(ierr);
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_GaussianXYZ,NULL);CHKERRQ(ierr);
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepX,NULL);CHKERRQ(ierr);
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepXY,NULL);CHKERRQ(ierr);
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepXYZ,NULL);CHKERRQ(ierr);
		ierr = VecShift(temperature,2.0);CHKERRQ(ierr);
	}	
	if (model_setup == 1) {
		ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_GaussianXY,NULL);CHKERRQ(ierr);
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_GaussianXYZ,NULL);CHKERRQ(ierr);
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepX,NULL);CHKERRQ(ierr);
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepXY,NULL);CHKERRQ(ierr);
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepXYZ,NULL);CHKERRQ(ierr);
		ierr = VecShift(temperature,2.0);CHKERRQ(ierr);
	}	
	if (model_setup == 2) {
		ierr = VecZeroEntries(temperature);CHKERRQ(ierr);
	}	
	if (model_setup == 3) {
		int direction;
		
		//ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepX,NULL);CHKERRQ(ierr);
		/*
		direction = 0;
		ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepWithDirection,(void*)&direction);CHKERRQ(ierr);
		*/
		//
		direction = 1;
		ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepWithDirection,(void*)&direction);CHKERRQ(ierr);
		//
		/*
		direction = 2;
		ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_StepWithDirection,(void*)&direction);CHKERRQ(ierr);
		 */
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_AdvDiffExample"
PetscErrorCode ModelApplyUpdateMeshGeometry_AdvDiffExample(pTatinCtx c,Vec X,void *ctx)
{
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
	Vec            temperature;
	char           name[256];
	PetscBool      write_markers,output_markers_once = PETSC_FALSE;
	static int     been_here = 0;
	
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* stokes variables */
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	
	/* energy variables */
	ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
	
	ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);

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
	ierr = PetscOptionsGetBool(NULL,"-advdiff_output_markers_once",&output_markers_once,NULL);CHKERRQ(ierr);
	write_markers = PETSC_TRUE;
	if ((output_markers_once) && (been_here != 0)) {
		write_markers = PETSC_FALSE;
	}
	if (write_markers) {
		const int                  nf = 2;
		const MaterialPointField   mp_prop_list[] = { MPField_Std, MPField_Energy }; 

		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
		sprintf(name,"%s_advdiff_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,name);CHKERRQ(ierr);
	}

	if (write_markers) {
		const int                   nf = -1;
		const MaterialPointVariable mp_prop_list[] = { MPV_viscosity, MPV_density }; 
		
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
		sprintf(name,"%s_mpoints_cell",prefix);
		ierr = pTatinOutputParaViewMarkerFields(c->stokes_ctx->stokes_pack,materialpoint_db,nf,mp_prop_list,c->outputpath,name);CHKERRQ(ierr);
	}
	
	been_here = 1;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_AdvDiffExample"
PetscErrorCode ModelDestroy_AdvDiffExample(pTatinCtx c,void *ctx)
{
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_AdvDiffExample"
PetscErrorCode pTatinModelRegister_AdvDiffExample(void)
{
	pTatinModel m;
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

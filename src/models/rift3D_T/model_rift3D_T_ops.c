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
 **    Filename:      model_rift3D_T_ops.c
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
 **    $Id: model_rift3D_T_ops.c 3960 2013-03-06 15:28:23Z dmay $
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

/*  
 Developed by Laetitia Le Pourhiet [laetitia.le_pourhiet@upmc.fr] 
*/


#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPStokesPl_def.h"
#include "MPntPEnergy_def.h"
#include "stokes_form_function.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_iterator.h"
#include "mesh_update.h"
#include "output_material_points.h"
#include "material_point_std_utils.h"
#include "material_point_utils.h"
#include "material_point_popcontrol.h"
#include "energy_output.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"

#include "rift3D_T_ctx.h"

PetscErrorCode ModelApplyUpdateMeshGeometry_Rift3D_T_semi_eulerian(pTatinCtx c,Vec X,void *ctx);
PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D_T_semi_eulerian(pTatinCtx c,void *ctx);

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Rift3D_T"
PetscErrorCode ModelInitialize_Rift3D_T(pTatinCtx c,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	RheologyConstants      *rheology;
	DataBucket materialconstants = c->material_constants;
	PetscBool nondim;
	PetscScalar vx,vy,vz,Sx,Sy,Sz;
	PetscInt regionidx;
    PetscReal cm_per_yer2m_per_sec = 1.0e-2 / ( 365.0 * 24.0 * 60.0 * 60.0 ) ;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscPrintf(PETSC_COMM_WORLD,"Rift model expects the following dimensions for input\n");
	PetscPrintf(PETSC_COMM_WORLD," Box geometry: [m] \n");
	PetscPrintf(PETSC_COMM_WORLD," Viscosity:    [Pa.s] \n");
	PetscPrintf(PETSC_COMM_WORLD," Velocity:     [m/sec] \n");
	PetscPrintf(PETSC_COMM_WORLD," Density:      [kg/m^3] \n");	
	
	PetscPrintf(PETSC_COMM_WORLD,"if you wish to use non dimensional input you must add -model_rift3D_T_dimensional \n");
	rheology                = &c->rheology_constants;
	rheology->rheology_type = RHEOLOGY_VP_STD;
	/* I REALLY DONT LIKE THE FOLLOWING ONE, SHOULD BE  in model data */
	rheology->nphases_active = 5;
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global = 1.e+25;
	rheology->eta_lower_cutoff_global = 1.e+19;
	data->runmises = PETSC_FALSE;	
	/* set the deffault values of the material constant for this particular model */
	/*scaling */ 
	data->length_bar    = 100.0 * 1.0e3;
	data->viscosity_bar = 1e25;
	data->velocity_bar  = 1.0e-10;
	data->dimensional   = PETSC_TRUE;
	/* box geometry, m */
	data->Lx =  6.0e5;
	data->Ly =  0.0e5;
	data->Lz =  6.0e5;
    //data->Ox =  -6.0e5;
	data->Ox =  0.0e5;
	data->Oy =  -1.5e5;
	data->Oz =  0.0e5;
	/* velocity cm/y */
	vx = 0.5*cm_per_yer2m_per_sec;
	vz = 0.1*cm_per_yer2m_per_sec;
	/* rho0 for initial pressure*/ 
	data->rho0 = 3140.0;
    /*Temperature */ 
    data->Tbottom = 1400.0;
    data->Ttop    = 0.0;
    data->thermal_age0 = 200;
    data->thermal_age_anom = 25;
    data->wx_anom  = 1;
    data->wz_anom  = 0.5;
    data->cx_anom  = 3.0;
    data->cz_anom  = 0.0;
	/* Material constant */
	MaterialConstantsSetDefaults(materialconstants);

	MaterialConstantsSetValues_MaterialType(materialconstants,0,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
    //MaterialConstantsSetValues_MaterialType(materialconstants,0,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
	MaterialConstantsSetValues_ViscosityFK(materialconstants,0,1.0e27,0.020);
    
    //MaterialConstantsSetValues_MaterialType(materialconstants,0,VISCOUS_ARRHENIUS,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	//MaterialConstantsSetValues_ViscosityArrh(materialconstants,0,1.0e27,0.05357142857145);
    
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,0,2700,2.e-5,0.0);
    MaterialConstantsSetValues_DensityConst(materialconstants,0,2700);
	MaterialConstantsSetValues_PlasticDP(materialconstants,0,0.6,0.1,2.e7,2.e7,1.e7,2.e8);
	MaterialConstantsSetValues_PlasticMises(materialconstants,0,1.e8,1.e8);
    MaterialConstantsSetValues_SoftLin(materialconstants,0,0.0,0.1);
    
    //MaterialConstantsSetValues_MaterialType(materialconstants,1,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
	MaterialConstantsSetValues_MaterialType(materialconstants,1,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);    
	MaterialConstantsSetValues_ViscosityFK(materialconstants,1,1.0e27,0.020);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,1,2700,5.e-6,3.e-12);
    MaterialConstantsSetValues_DensityConst(materialconstants,1,2800);
	MaterialConstantsSetValues_PlasticDP(materialconstants,1,0.6,0.1,2.e7,2.e7,1.e7,2.e8);
	MaterialConstantsSetValues_PlasticMises(materialconstants,1,1.e8,1.e8);
    MaterialConstantsSetValues_SoftLin(materialconstants,1,0.0,0.1);    
    
    //MaterialConstantsSetValues_MaterialType(materialconstants,2,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
	MaterialConstantsSetValues_MaterialType(materialconstants,2,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);    
	MaterialConstantsSetValues_ViscosityFK(materialconstants,2,1.0e27,0.012);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,2,3300,5.e-6,3.e-12);
	MaterialConstantsSetValues_DensityConst(materialconstants,2,3300);
    MaterialConstantsSetValues_PlasticDP(materialconstants,2,0.0,0.0,3.e8,2.e8,2.e7,3.e8);
	MaterialConstantsSetValues_PlasticMises(materialconstants,2,3.e8,3.e8);
    MaterialConstantsSetValues_SoftLin(materialconstants,2,0.0,0.1);
    
    //MaterialConstantsSetValues_MaterialType(materialconstants,3,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
	MaterialConstantsSetValues_MaterialType(materialconstants,3,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	MaterialConstantsSetValues_ViscosityFK(materialconstants,3,1.0e27,0.012);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,3,3300,5.e-6,3.e-12);
	MaterialConstantsSetValues_DensityConst(materialconstants,3,3300);
    MaterialConstantsSetValues_PlasticDP(materialconstants,3,0.3,0.3,2.e7,2.e7,1.e7,2.e8);
    MaterialConstantsSetValues_PlasticMises(materialconstants,3,3.e8,3.e8);
    MaterialConstantsSetValues_SoftLin(materialconstants,3,0.0,0.1);
    
    //MaterialConstantsSetValues_MaterialType(materialconstants,4,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
    MaterialConstantsSetValues_MaterialType(materialconstants,4,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);    
	MaterialConstantsSetValues_ViscosityFK(materialconstants,4,1.0e27,0.02);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,4,3300,5.e-6,3.e-12);
	MaterialConstantsSetValues_DensityConst(materialconstants,4,3300);
    MaterialConstantsSetValues_PlasticDP(materialconstants,4,0.6,0.6,2.e7,2.e7,1.e7,3.e8);
	MaterialConstantsSetValues_PlasticMises(materialconstants,4,3.e8,3.e8);
    MaterialConstantsSetValues_SoftLin(materialconstants,4,0.0,0.1);
    
	/* Read the options */
	/*cutoff */
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_rift3D_T_apply_viscosity_cutoff_global",&rheology->apply_viscosity_cutoff_global,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_eta_lower_cutoff_global",&rheology->eta_lower_cutoff_global,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_eta_upper_cutoff_global",&rheology->eta_upper_cutoff_global,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_rift3D_T_runwithmises",&data->runmises,PETSC_NULL);CHKERRQ(ierr);
	/*scaling */     
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_rift3D_T_nondimensional",&nondim,PETSC_NULL);CHKERRQ(ierr);
	if (nondim){
		data->dimensional = PETSC_FALSE;
	} else {
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_vis_bar",&data->viscosity_bar,PETSC_NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_vel_bar",&data->velocity_bar,PETSC_NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_length_bar",&data->length_bar,PETSC_NULL);CHKERRQ(ierr); 
	}
	/* box geometry, m */	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_Lx",&data->Lx,PETSC_NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_Ly",&data->Ly,PETSC_NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_Lz",&data->Lz,PETSC_NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_Ox",&data->Ox,PETSC_NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_Oy",&data->Oy,PETSC_NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_Oz",&data->Oz,PETSC_NULL);CHKERRQ(ierr); 
	
	/* velocity cm/y */    
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_vx",&vx,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_vz",&vz,PETSC_NULL);CHKERRQ(ierr); 
	
	/* rho0 for initial pressure */
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_rho0",&data->rho0,PETSC_NULL);CHKERRQ(ierr); 
	
    /* temperature initial condition */ 
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_Tbot",&data->Tbottom,PETSC_NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_Ttop",&data->Ttop,PETSC_NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_age0",&data->thermal_age0,PETSC_NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_ageAnom",&data->thermal_age_anom,PETSC_NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_wx",&data->wx_anom,PETSC_NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_wz",&data->wz_anom,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_cx",&data->cx_anom,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_T_cz",&data->cz_anom,PETSC_NULL);CHKERRQ(ierr);
 
    
	/* Material constant */
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
		PetscPrintf(PETSC_COMM_WORLD,"reading options");
		ierr= MaterialConstantsSetFromOptions(materialconstants,"model_rift3D_T",regionidx,PETSC_FALSE);CHKERRQ(ierr);
	}
	
	/*Compute velocity at bottom*/
	Sx = (data->Ly - data->Oy)*(data->Lz - data->Oz);
	Sz = (data->Ly - data->Oy)*(data->Lx - data->Ox);
	Sy = (data->Lx - data->Ox)*(data->Lz - data->Oz);
	vy = (vx*Sx-vz*Sz)/Sy;
	
	/* reports before scaling */
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_T_Ox %+1.4e [SI] -model_rift3D_T_Lx : %+1.4e [SI]\n", data->Ox ,data->Lx );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_T_Oy %+1.4e [SI] -model_rift3D_T_Ly : %+1.4e [SI]\n", data->Oy ,data->Ly );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_T_Oz %+1.4e [SI] -model_rift3D_T_Lz : %+1.4e [SI]\n", data->Oz ,data->Lz );
	PetscPrintf(PETSC_COMM_WORLD,"  -model_rift3D_T_vx [m/s]:  %+1.4e  -model_rift3D_T_vz [m/s]:  %+1.4e : computed vy [m/s]:  %+1.4e \n", vx,vz,vy);
	PetscPrintf(PETSC_COMM_WORLD,"-model_rift3D_T_rho0 [kg/m^3] :%+1.4e \n", data->rho0 );
    PetscPrintf(PETSC_COMM_WORLD,"-model_rift3D_T_Tbot:%+1.4e \t -model_rift3D_T_Ttop:%+1.4e \t -model_rift3D_T_age0:%+1.4e \n",data->Tbottom,data->Ttop,  data->thermal_age0);
    PetscPrintf(PETSC_COMM_WORLD,"ageAnom:%+1.4e \t wx:%+1.4e \t wz:%+1.4e cx:%+1.4e \t cz:%+1.4e \n",data->thermal_age_anom,data->wx_anom,data->wz_anom,data->cx_anom,data->cz_anom);
    
    
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) { 
		MaterialConstantsPrintAll(materialconstants,regionidx);
	} 
	
	if (data->dimensional) {
		/*Compute additional scaling parameters*/
		data->time_bar      = data->length_bar / data->velocity_bar;
		data->pressure_bar  = data->viscosity_bar/data->time_bar;
		data->density_bar   = data->pressure_bar / data->length_bar;
		
        
        
		PetscPrintf(PETSC_COMM_WORLD,"[rift3D_T]:  during the solve scaling will be done using \n");
		PetscPrintf(PETSC_COMM_WORLD,"  L*    : %1.4e [m]\n", data->length_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  U*    : %1.4e [m.s^-1]\n", data->velocity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  t*    : %1.4e [s]\n", data->time_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  eta*  : %1.4e [Pa.s]\n", data->viscosity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  rho*  : %1.4e [kg.m^-3]\n", data->density_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  P*    : %1.4e [Pa]\n", data->pressure_bar );
		//scale viscosity cutoff
		rheology->eta_lower_cutoff_global = rheology->eta_lower_cutoff_global / data->viscosity_bar;
		rheology->eta_upper_cutoff_global = rheology->eta_upper_cutoff_global / data->viscosity_bar;
		//scale length 
		data->Lx = data->Lx / data->length_bar;
		data->Ly = data->Ly / data->length_bar;
		data->Lz = data->Lz / data->length_bar;
		data->Ox = data->Ox / data->length_bar;
		data->Oy = data->Oy / data->length_bar;
		data->Oz = data->Oz / data->length_bar; 
		
		//scale velocity
		data->vx = vx/data->velocity_bar;
		data->vy = vy/data->velocity_bar;
		data->vz = vz/data->velocity_bar;
		//scale rho0
		data->rho0 = data->rho0/data->density_bar;
		
		// scale material properties
		for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
			MaterialConstantsScaleAll(materialconstants,regionidx,data->length_bar,data->velocity_bar,data->time_bar,data->viscosity_bar,data->density_bar,data->pressure_bar);
		}
		

		/*Reports scaled values*/ 
		
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_T_Ox   :  %+1.4e    -model_rift3D_T_Lx   :  %+1.4e  \n", data->Ox ,data->Lx );
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_T_Oy   :  %+1.4e    -model_rift3D_T_Ly   :  %+1.4e \n", data->Oy, data->Ly );
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_T_Oz   :  %+1.4e    -model_rift3D_T_Lz   :  %+1.4e\n", data->Oz , data->Lz );
		
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_T_Vx:%+1.4e    -model_rift3D_T_vy:%+1.4e    -model_rift3D_T_vz:  %+1.4e \n", data->vx ,data->vy, data->vz);
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_T_rho0:%+1.4e \n", data->rho0 );
		PetscPrintf(PETSC_COMM_WORLD,"scaled value for material parameters\n");
		for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) { 
			MaterialConstantsPrintAll(materialconstants,regionidx);
		}    
	}
	
	/* force energy equation to be introduced */
	ierr = PetscOptionsInsertString("-activate_energy");CHKERRQ(ierr);

	data->use_semi_eulerian_mesh = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_rift3D_T_use_semi_eulerian",&data->use_semi_eulerian_mesh,PETSC_NULL);CHKERRQ(ierr);
	if (data->use_semi_eulerian_mesh) {
		pTatinModel model;
		
		PetscPrintf(PETSC_COMM_WORLD,"rift3D_T: activating semi Eulerian mesh advection\n");
		ierr = pTatinGetModel(c,&model);CHKERRQ(ierr);
		ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Rift3D_T_semi_eulerian);CHKERRQ(ierr);
		ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Rift3D_T_semi_eulerian);CHKERRQ(ierr);
	}	

	data->output_markers = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_rift3D_T_output_markers",&data->output_markers,PETSC_NULL);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


/* 
 
 1/ Define boundary conditions in one place for this model.
 
 2/ Calling pattern should always be
 PetscErrorCode ModelRift3D_T_DefineBCList(BCList bclist,DM dav,pTatinCtx user,ModelRift3D_TCtx data)
 where ModelRift3D_TCtx data is a different type for each model.
 
 3/ Re-use this function in 
 ModelApplyBoundaryCondition_Rift3D_T();
 ModelApplyBoundaryConditionMG_Rift3D_T();
 
 */
#undef __FUNCT__
#define __FUNCT__ "ModelRift3D_T_DefineBCList"
PetscErrorCode ModelRift3D_T_DefineBCList(BCList bclist,DM dav,pTatinCtx user,ModelRift3D_TCtx *data)
{
	PetscScalar    vxl,vxr,vzf,vzb,vy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	vxl = 0.0;//-data->vx;
	vxr =  data->vx;
	vy  =  data->vy;
	vzf = -data->vz;
	vzb =  0.0;//data->vz;
	
	/* infilling free slip base */
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&vy);CHKERRQ(ierr);
	
	/* free surface top*/
	
	/*extension along face of normal x */ 
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&(vxl));CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&(vxr));CHKERRQ(ierr);
	
	/*compression along face of normal z */ 
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&(vzb));CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&(vzf));CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Rift3D_T"
PetscErrorCode ModelApplyBoundaryCondition_Rift3D_T(pTatinCtx user,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	PetscBool active_energy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = ModelRift3D_T_DefineBCList(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,user,data);CHKERRQ(ierr);
	
	/* set boundary conditions for temperature */
	ierr = pTatinContextValid_Energy(user,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PetscReal      val_T;
		PhysCompEnergy energy;
		BCList         bclist;
		DM             daT;
		
		ierr   = pTatinGetContext_Energy(user,&energy);CHKERRQ(ierr);
		daT    = energy->daT;
		bclist = energy->T_bclist;
		
		
		val_T = data->Tbottom;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);
		val_T = data->Ttop;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);		
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Rift3D_T"
PetscErrorCode ModelApplyBoundaryConditionMG_Rift3D_T(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	PetscInt       n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = ModelRift3D_T_DefineBCList(bclist[n],dav[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Rift3D_T"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D_T(pTatinCtx c,void *ctx)
{
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]] - Not implemented \n", __FUNCT__);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Rift3D_T_semi_eulerian"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D_T_semi_eulerian(pTatinCtx c,void *ctx)
{
	ModelRift3D_TCtx   *data = (ModelRift3D_TCtx*)ctx;
	PhysCompStokes     stokes;
	DM                 stokes_pack,dav,dap;
	PetscInt           Nxp[2];
	PetscReal          perturb;
	DataBucket         material_point_db,material_point_face_db;
	PetscInt           f, n_face_list=3, face_list[] = { 3, 4, 5 };
	int                p,n_mp_points;
	MPAccess           mpX;
	PetscErrorCode     ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

	
	ierr = pTatinGetMaterialPoints(c,&material_point_db,PETSC_NULL);CHKERRQ(ierr);

#if 0
	{
		const int nf = 2;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes };
		char mp_file_prefix[1024];
		
		sprintf(mp_file_prefix,"mpoints_remesh_vol0");
		ierr = SwarmViewGeneric_ParaView(material_point_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
#endif
	
	/* create face storage for markers */
	DataBucketDuplicateFields(material_point_db,&material_point_face_db);
	
	for (f=0; f<n_face_list; f++) {
		
		/* traverse */
		/* [0,1/east,west] ; [2,3/north,south] ; [4,5/front,back] */
		Nxp[0]  = 2;
		Nxp[1]  = 2;
		perturb = 0.1;

		/* reset size */
		DataBucketSetSizes(material_point_face_db,0,-1);

		/* assign coords */
		ierr = SwarmMPntStd_CoordAssignment_FaceLatticeLayout3d(dav,Nxp,perturb, face_list[f], material_point_face_db);CHKERRQ(ierr);

		/* assign values */
		DataBucketGetSizes(material_point_face_db,&n_mp_points,0,0);
		ierr = MaterialPointGetAccess(material_point_face_db,&mpX);CHKERRQ(ierr);
		for (p=0; p<n_mp_points; p++) {
			ierr = MaterialPointSet_phase_index(mpX,p,MATERIAL_POINT_PHASE_UNASSIGNED);CHKERRQ(ierr);
		}
		ierr = MaterialPointRestoreAccess(material_point_face_db,&mpX);CHKERRQ(ierr);
		
		/* output */
#if 0
		{
			const int nf = 2;
			const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes };
			char mp_file_prefix[1024];
			
			sprintf(mp_file_prefix,"mpoints_remesh_face%d",face_list[f]);
			ierr = SwarmViewGeneric_ParaView(material_point_face_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
		}
#endif
		
		/* insert into volume bucket */
		DataBucketInsertValues(material_point_db,material_point_face_db);
	}	

	/* Copy ALL values from nearest markers to newly inserted markers expect (xi,xip,pid) */
	ierr = MaterialPointRegionAssignment_v1(material_point_db,dav);CHKERRQ(ierr);
	
	/* reset any variables */
	DataBucketGetSizes(material_point_face_db,&n_mp_points,0,0);
	ierr = MaterialPointGetAccess(material_point_face_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mp_points; p++) {
		ierr = MaterialPointSet_plastic_strain(mpX,p,0.0);CHKERRQ(ierr);
		ierr = MaterialPointSet_yield_indicator(mpX,p,0);CHKERRQ(ierr);
	}
	ierr = MaterialPointRestoreAccess(material_point_face_db,&mpX);CHKERRQ(ierr);
	
	/* re-assign pid's for new particles such that they are consistent with the original volume marker set */
	
	
#if 1
	{
		const int nf = 2;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes };
		char mp_file_prefix[1024];
		
		sprintf(mp_file_prefix,"mpoints_remesh_vol1");
		ierr = SwarmViewGeneric_ParaView(material_point_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
#endif	
	
	/* delete */
	DataBucketDestroy(&material_point_face_db);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Rift3D_T"
PetscErrorCode ModelApplyInitialMeshGeometry_Rift3D_T(pTatinCtx c,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,data->Ox,data->Lx,data->Oy,data->Ly,data->Oz,data->Lz);CHKERRQ(ierr);
	
	/* note - Don't access the energy mesh here, its not yet created */
	/* note - The initial velocity mesh geometry will be copied into the energy mesh */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Rift3D_T"
PetscErrorCode ModelApplyInitialMaterialGeometry_Rift3D_T(pTatinCtx c,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	int                    e,p,n_mp_points;
	PetscScalar            y_lab,y_moho,y_midcrust,notch_l,notch_w2,xc;
	DataBucket             db;
	DataField              PField_std,PField_pls;
	int                    phase;
	MPAccess               mpX;
	PetscErrorCode         ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	/* define properties on material points */
	db = c->materialpoint_db;
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	
	DataBucketGetDataFieldByName(db,MPntPStokesPl_classname,&PField_pls);
	DataFieldGetAccess(PField_pls);
	DataFieldVerifyAccess(PField_pls,sizeof(MPntPStokesPl));
	
	
	/* m */
	y_lab      = -120.0e3; 
	y_moho     = -40.0e3;
	y_midcrust = -20.0e3;
	notch_w2   = 50.e3;
    notch_l    = 150.e3;
    xc         = (data->Lx + data->Ox)/2.0* data->length_bar;
    xc         = 0.0; 
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	srand(0);
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd       *material_point;
		MPntPStokes   *mpprop_stokes;
		MPntPStokesPl *mpprop_pls;
		double        *position,ycoord,xcoord,zcoord;
		float         pls;
		char          yield;
		PetscBool     norandomiseplastic;
        
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_pls,p,(void**)&mpprop_pls);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		/* convert to scaled units */
		ycoord = position[1] * data->length_bar;
		xcoord=position[0] * data->length_bar;
		zcoord=position[2] * data->length_bar;
		
		if (ycoord<y_lab) {
			phase = 3;
		} else if (ycoord<y_moho) {
			phase = 2;
		} else if (ycoord<y_midcrust) {
			phase = 1;
		} else {
			phase = 0;
		}
        norandomiseplastic = PETSC_FALSE;
        ierr = PetscOptionsGetBool(PETSC_NULL,"-model_rift3D_T_norandom",&norandomiseplastic,PETSC_NULL);CHKERRQ(ierr);
		if (norandomiseplastic) {
            
            pls   = 0.0;
            if (xcoord>250.e3 && xcoord < 350.e3 &&  zcoord < 150.e3){
                pls = 0.05;   
            }    
        }else{
            pls   = rand()/(RAND_MAX+0.5)*0.001;
            if (abs(xcoord - xc) < notch_w2 &&  zcoord < notch_l){
                pls = rand()/(RAND_MAX+0.5)*0.1;
            }
            
        }
        
        yield = 0; 
			/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,phase);
		MPntPStokesPlSetField_yield_indicator(mpprop_pls,yield);
		MPntPStokesPlSetField_plastic_strain(mpprop_pls,pls);
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_pls);

	
	ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mp_points; p++) {
		double kappa,H;
		
		ierr = MaterialPointGet_phase_index(mpX,p,&phase);CHKERRQ(ierr);

		kappa = 1.0e-6/data->length_bar/data->length_bar*data->time_bar;
		H     = 0.0;
		ierr = MaterialPointSet_diffusivity(mpX,p,kappa);CHKERRQ(ierr);
		ierr = MaterialPointSet_heat_source(mpX,p,H);CHKERRQ(ierr);
	}
	ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);
    
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Rift3D_T"
PetscErrorCode ModelApplyUpdateMeshGeometry_Rift3D_T(pTatinCtx c,Vec X,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	PetscReal        step;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* fully lagrangian update */
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Rift3D_T_semi_eulerian"
PetscErrorCode ModelApplyUpdateMeshGeometry_Rift3D_T_semi_eulerian(pTatinCtx c,Vec X,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	PetscReal        step;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* fully lagrangian update */
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Rift3D_T_CheckScales"
PetscErrorCode ModelOutput_Rift3D_T_CheckScales(pTatinCtx c,Vec X)
{
	Vec Xcopy,velocity,pressure,F,RHS;
	PetscReal fu,fp;
	PetscInt Nu,Np;
	PhysCompStokes    stokes;
	DM                stokes_pack;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	
	
	ierr = VecDuplicate(X,&Xcopy);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&RHS);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"[rift3D_T]: check scales \n");
	
	ierr = VecCopy(X,Xcopy);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,Xcopy,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = VecStrideMin(pressure,0,PETSC_NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," min|P0|   = %+1.4e \n",fp);
	ierr = VecStrideMax(pressure,0,PETSC_NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," max|P0|   = %+1.4e \n",fp);
	
	ierr = VecStrideMin(pressure,1,PETSC_NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," min|dPdx| = %+1.4e \n",fp);
	ierr = VecStrideMax(pressure,1,PETSC_NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," max|dPdx| = %+1.4e \n",fp);
	
	ierr = VecStrideMin(pressure,2,PETSC_NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," min|dPdy| = %+1.4e \n",fp);
	ierr = VecStrideMax(pressure,2,PETSC_NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," max|dPdy| = %+1.4e \n",fp);
	
	ierr = VecStrideMin(pressure,3,PETSC_NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," min|dPdz| = %+1.4e \n",fp);
	ierr = VecStrideMax(pressure,3,PETSC_NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," max|dPdz| = %+1.4e \n",fp);
	
	
	ierr = DMCompositeRestoreAccess(stokes_pack,Xcopy,&velocity,&pressure);CHKERRQ(ierr);
	
	
	
	ierr = VecZeroEntries(Xcopy);CHKERRQ(ierr);
	ierr = FormFunction_Stokes(PETSC_NULL,Xcopy,RHS,(void*)c);CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(stokes_pack,RHS,&velocity,&pressure);CHKERRQ(ierr);
	ierr = BCListInsertZero(stokes->u_bclist,velocity);CHKERRQ(ierr);
	ierr = VecGetSize(velocity,&Nu);CHKERRQ(ierr);
	ierr = VecGetSize(pressure,&Np);CHKERRQ(ierr);
	
	ierr = VecNorm(velocity,NORM_2,&fu);CHKERRQ(ierr);
	ierr = VecNorm(pressure,NORM_2,&fp);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD," |rho.g|    = %1.4e \n",fu/sqrt(Nu));
	PetscPrintf(PETSC_COMM_WORLD," |cont_rhs| = %1.4e \n",fp/sqrt(Np));
	ierr = DMCompositeRestoreAccess(stokes_pack,RHS,&velocity,&pressure);CHKERRQ(ierr);
	
	
	
	ierr = VecCopy(X,Xcopy);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,Xcopy,&velocity,&pressure);CHKERRQ(ierr);
	ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(stokes_pack,Xcopy,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = FormFunction_Stokes(PETSC_NULL,Xcopy,F,(void*)c);CHKERRQ(ierr);
	ierr = VecAXPY(F,1.0,RHS);CHKERRQ(ierr); /* F = F - RHS */
	
	ierr = DMCompositeGetAccess(stokes_pack,F,&velocity,&pressure);CHKERRQ(ierr);
	ierr = BCListInsertZero(stokes->u_bclist,velocity);CHKERRQ(ierr);
	ierr = VecNorm(velocity,NORM_2,&fu);CHKERRQ(ierr);
	ierr = VecNorm(pressure,NORM_2,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," |div(sigma_ij)| = %1.4e \n",fu/sqrt(Nu));
	PetscPrintf(PETSC_COMM_WORLD," |div(u_i)|      = %1.4e \n",fp/sqrt(Np));
	ierr = DMCompositeRestoreAccess(stokes_pack,F,&velocity,&pressure);CHKERRQ(ierr);
	
	
	ierr = VecCopy(X,Xcopy);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,Xcopy,&velocity,&pressure);CHKERRQ(ierr);
	ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(stokes_pack,Xcopy,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = FormFunction_Stokes(PETSC_NULL,Xcopy,F,(void*)c);CHKERRQ(ierr);
	ierr = VecAXPY(F,1.0,RHS);CHKERRQ(ierr); /* F = F - RHS */
	
	ierr = DMCompositeGetAccess(stokes_pack,F,&velocity,&pressure);CHKERRQ(ierr);
	ierr = BCListInsertZero(stokes->u_bclist,velocity);CHKERRQ(ierr);
	ierr = VecNorm(velocity,NORM_2,&fu);CHKERRQ(ierr);
	ierr = VecNorm(pressure,NORM_2,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," |grad(P)| = %1.4e \n",fu/sqrt(Nu));
	ierr = DMCompositeRestoreAccess(stokes_pack,F,&velocity,&pressure);CHKERRQ(ierr);
	
	
	ierr = VecDestroy(&Xcopy);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = VecDestroy(&RHS);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Rift3D_T"
PetscErrorCode ModelOutput_Rift3D_T(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelRift3D_TCtx  *data = (ModelRift3D_TCtx*)ctx;
	PetscBool         active_energy;
	DataBucket        materialpoint_db;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	//ierr = ModelOutput_Rift3D_T_CheckScales(c,X);CHKERRQ(ierr);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	
	if (data->output_markers)
	{
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
		//  Write out just the stokes variable?
		//  const int nf = 1;
		//  const MaterialPointField mp_prop_list[] = { MPField_Stokes };
		//
		//  Write out just std, stokes and plastic variables
		const int nf = 4;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
		char mp_file_prefix[256];
		
		sprintf(mp_file_prefix,"%s_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
	
		
	{
		const int                   nf = 3;
		const MaterialPointVariable mp_prop_list[] = { MPV_viscosity, MPV_density, MPV_plastic_strain }; 
		
		ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
	}	
	
	
	/* standard viewer */
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PhysCompEnergy energy;
		Vec            temperature;
		
		ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(c,&temperature,PETSC_NULL);CHKERRQ(ierr);

		ierr = pTatin3d_ModelOutput_Temperature_Energy(c,temperature,prefix);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Rift3D_T"
PetscErrorCode ModelDestroy_Rift3D_T(pTatinCtx c,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialStokesVariableMarkers_Rift3D_T"
PetscErrorCode ModelApplyInitialStokesVariableMarkers_Rift3D_T(pTatinCtx user,Vec X,void *ctx)
{
	
	DM                stokes_pack,dau,dap;
	PhysCompStokes    stokes;
	Vec               Uloc,Ploc;
	PetscScalar       *LA_Uloc,*LA_Ploc;
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
    DataField                    PField;
	MaterialConst_MaterialType   *truc;
	PetscErrorCode    ierr;
	PetscInt regionidx;	
	PetscFunctionBegin;
	
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	if (!data->runmises) {
        DataBucketGetDataFieldByName(user->material_constants,MaterialConst_MaterialType_classname,&PField);
        DataFieldGetAccess(PField);
		for (regionidx=0; regionidx<user->rheology_constants.nphases_active;regionidx++) {
            DataFieldAccessPoint(PField,regionidx,(void**)&truc);
            MaterialConst_MaterialTypeSetField_plastic_type(truc,PLASTIC_MISES);
		}
        DataFieldRestoreAccess(PField);
	}
	
	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	ierr = pTatin_EvaluateRheologyNonlinearities(user,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	if (!data->runmises) {
		for (regionidx=0; regionidx<user->rheology_constants.nphases_active;regionidx++) {
            DataBucketGetDataFieldByName(user->material_constants,MaterialConst_MaterialType_classname,&PField);
            DataFieldGetAccess(PField);
            for (regionidx=0; regionidx<user->rheology_constants.nphases_active;regionidx++) {
                DataFieldAccessPoint(PField,regionidx,(void**)&truc);
                MaterialConst_MaterialTypeSetField_plastic_type(truc,PLASTIC_DP);
            }
            DataFieldRestoreAccess(PField);

		}
	}
	
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialCondition_Rift3D_T"
PetscErrorCode ModelApplyInitialCondition_Rift3D_T(pTatinCtx c,Vec X,void *ctx)
{
	ModelRift3D_TCtx *data = (ModelRift3D_TCtx*)ctx;
	DM stokes_pack,dau,dap;
	Vec velocity,pressure;
	PetscReal vxl,vxr,vzb,vzf,vy;
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	DMDAVecTraverse3d_InterpCtx IntpCtx;
	PetscReal MeshMin[3],MeshMax[3],domain_height;
	PetscBool active_energy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	stokes_pack = c->stokes_ctx->stokes_pack;
	
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
    ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    
    vxl = -data->vx;
	vxr =  data->vx;
	vy  =  data->vy;
	vzf = -data->vz;
	vzb =  0.0;//data->vz;
	
	ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
	/* apply -5 < vx 5 across the domain x \in [0,1] */
	
	ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vxr-vxl)/(data->Lx-data->Ox),vxl,0.0);CHKERRQ(ierr);
	ierr = DMDAVecTraverse3d(dau,velocity,0,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
	ierr = DMDAVecTraverse3d_InterpCtxSetUp_Z(&IntpCtx,(vzf-vzb)/(data->Lz-data->Oz),vzb,0.0);CHKERRQ(ierr);
	ierr = DMDAVecTraverse3d(dau,velocity,2,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
	
	ierr = DMDAVecTraverse3d_InterpCtxSetUp_Y(&IntpCtx,-vy/(data->Ly-data->Oy),0.0,0.0);CHKERRQ(ierr);
	ierr = DMDAVecTraverse3d(dau,velocity,1,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
	
	
	ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
	
	
	ierr = DMDAGetBoundingBox(dau,MeshMin,MeshMax);CHKERRQ(ierr);
	domain_height = MeshMax[1] - MeshMin[1];
	
	HPctx.surface_pressure = 0.0;
	HPctx.ref_height = domain_height;
	HPctx.ref_N      = c->stokes_ctx->my-1;
	HPctx.grav       = 10.0;
	HPctx.rho        = data->rho0;
	
	
  ierr = DMDAVecTraverseIJK(dap,pressure,0,DMDAVecTraverseIJK_HydroStaticPressure_v2,     (void*)&HPctx);CHKERRQ(ierr); /* P = P0 + a.x + b.y + c.z, modify P0 (idx=0) */
  ierr = DMDAVecTraverseIJK(dap,pressure,2,DMDAVecTraverseIJK_HydroStaticPressure_dpdy_v2,(void*)&HPctx);CHKERRQ(ierr); /* P = P0 + a.x + b.y + c.z, modify b  (idx=2) */
  ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,"testHP");CHKERRQ(ierr);
	
	
	/* initial condition for temperature */
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PhysCompEnergy energy;
		Vec            temperature;
		DM             daT;
		PetscReal      coeffs[8];
		
		ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(c,&temperature,PETSC_NULL);CHKERRQ(ierr);
		daT  = energy->daT;
    

        coeffs[0] = data->cx_anom ; 
        coeffs[1] = data->cz_anom ; 
        coeffs[2] = data->thermal_age0;
        coeffs[3] = data->thermal_age_anom;
        coeffs[4] = data->length_bar;
        coeffs[5] = data->Tbottom;
        coeffs[6] = data->wx_anom;
        coeffs[7] = data->wz_anom;
		ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse3d_ERFC3DFunctionXYZ,(void*)coeffs);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Rift3D_T"
PetscErrorCode pTatinModelRegister_Rift3D_T(void)
{
	ModelRift3D_TCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelRift3D_TCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelRift3D_TCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"rift3D_T");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Rift3D_T);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialCondition_Rift3D_T);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,   (void (*)(void))ModelApplyInitialStokesVariableMarkers_Rift3D_T);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Rift3D_T);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Rift3D_T);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Rift3D_T);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Rift3D_T);CHKERRQ(ierr);

	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Rift3D_T);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Rift3D_T);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Rift3D_T);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Rift3D_T);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

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
 **    filename:   model_rift3D_ops.c
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

/*
 Developed by Laetitia Le Pourhiet [laetitia.le_pourhiet@upmc.fr] 
*/


#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "data_bucket.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "stokes_form_function.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_iterator.h"

#include "rift3D_ctx.h"

PetscReal cm_per_yer2m_per_sec = 1.0e-2 / ( 365.0 * 24.0 * 60.0 * 60.0 ) ;

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Rift3D"
PetscErrorCode ModelInitialize_Rift3D(pTatinCtx c,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	RheologyConstants      *rheology;
	DataBucket materialconstants = c->material_constants;
	PetscBool nondim;
	PetscScalar vx,vy,vz,Sx,Sy,Sz;
	PetscInt regionidx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscPrintf(PETSC_COMM_WORLD,"Rift model expects the following dimensions for input\n");
	PetscPrintf(PETSC_COMM_WORLD," Box geometry: [m] \n");
	PetscPrintf(PETSC_COMM_WORLD," Viscosity:    [Pa.s] \n");
	PetscPrintf(PETSC_COMM_WORLD," Velocity:     [m/sec] \n");
	PetscPrintf(PETSC_COMM_WORLD," Density:      [kg/m^3] \n");	
	
	PetscPrintf(PETSC_COMM_WORLD,"if you wish to use non dimensional input you must add -model_rift3D_dimensional \n");
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
	data->Ox =  0.0e5;
	data->Oy =  -1.5e5;
	data->Oz =  0.0e5;
	/* velocity cm/y */
	vx = 0.5*cm_per_yer2m_per_sec;
	vz = 0.1*cm_per_yer2m_per_sec;
	/* rho0 for initial pressure*/ 
	data->rho0 = 3300.0;
	/* Material constant */
	MaterialConstantsSetDefaults(materialconstants);

	MaterialConstantsSetValues_MaterialType(materialconstants,0,VISCOUS_Z,PLASTIC_DP,SOFTENING_NONE,DENSITY_CONSTANT);
	MaterialConstantsSetValues_ViscosityZ(materialconstants,0,1.0e28,2000.,10000.);
	MaterialConstantsSetValues_DensityConst(materialconstants,0,2700);
	MaterialConstantsSetValues_PlasticDP(materialconstants,0,0.6,0.6,2.e7,2.e7,1.e7,3e8);
	MaterialConstantsSetValues_PlasticMises(materialconstants,0,1.e8,1.e8);
    
	MaterialConstantsSetValues_MaterialType(materialconstants,1,VISCOUS_Z,PLASTIC_DP,SOFTENING_NONE,DENSITY_CONSTANT);    
	MaterialConstantsSetValues_ViscosityZ(materialconstants,1,1.0e28,2000.,10000.);
	MaterialConstantsSetValues_DensityConst(materialconstants,1,2700);
	MaterialConstantsSetValues_PlasticDP(materialconstants,1,0.6,0.6,2.e7,2.e7,1.e7,3e8);
	MaterialConstantsSetValues_PlasticMises(materialconstants,1,1.e8,1.e8);
    
	MaterialConstantsSetValues_MaterialType(materialconstants,2,VISCOUS_Z,PLASTIC_DP,SOFTENING_NONE,DENSITY_CONSTANT);    
	MaterialConstantsSetValues_ViscosityZ(materialconstants,2,1.0e26,12000.,10000.);
	MaterialConstantsSetValues_DensityConst(materialconstants,2,3300);
	MaterialConstantsSetValues_PlasticDP(materialconstants,2,0.6,0.6,2.e7,2.e7,1.e7,3e8);
	MaterialConstantsSetValues_PlasticMises(materialconstants,2,3.e8,3.e8);
    
	MaterialConstantsSetValues_MaterialType(materialconstants,3,VISCOUS_Z,PLASTIC_DP,SOFTENING_NONE,DENSITY_CONSTANT);
	MaterialConstantsSetValues_ViscosityZ(materialconstants,3,1.0e26,12000.,10000.);
	MaterialConstantsSetValues_DensityConst(materialconstants,3,3280);
	MaterialConstantsSetValues_PlasticDP(materialconstants,3,0.6,0.6,2.e7,2.e7,1.e7,3e8);
    MaterialConstantsSetValues_PlasticMises(materialconstants,3,3.e8,3.e8);
    
    MaterialConstantsSetValues_MaterialType(materialconstants,4,VISCOUS_Z,PLASTIC_DP,SOFTENING_NONE,DENSITY_CONSTANT);    
	MaterialConstantsSetValues_ViscosityZ(materialconstants,4,1.0e28,2000.,10000.);
	MaterialConstantsSetValues_DensityConst(materialconstants,4,3300);
	MaterialConstantsSetValues_PlasticDP(materialconstants,4,0.6,0.6,2.e7,2.e7,1.e7,3e8);
	MaterialConstantsSetValues_PlasticMises(materialconstants,4,3.e8,3.e8);
 
	/* Read the options */
	/*cutoff */
	ierr = PetscOptionsGetBool(NULL,"-model_rift3D_apply_viscosity_cutoff_global",&rheology->apply_viscosity_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_eta_lower_cutoff_global",&rheology->eta_lower_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_eta_upper_cutoff_global",&rheology->eta_upper_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-model_rift3D_runwithmises",&data->runmises,NULL);CHKERRQ(ierr);
	/*scaling */     
	nondim = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-model_rift3D_nondimensional",&nondim,NULL);CHKERRQ(ierr);
	if (nondim){
		data->dimensional = PETSC_FALSE;
	} else {
    ierr = PetscOptionsGetReal(NULL,"-model_rift3D_vis_bar",&data->viscosity_bar,NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(NULL,"-model_rift3D_vel_bar",&data->velocity_bar,NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetReal(NULL,"-model_rift3D_length_bar",&data->length_bar,NULL);CHKERRQ(ierr); 
	}
	/* box geometry, m */	
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_Lx",&data->Lx,NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_Ly",&data->Ly,NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_Lz",&data->Lz,NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_Ox",&data->Ox,NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_Oy",&data->Oy,NULL);CHKERRQ(ierr); 
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_Oz",&data->Oz,NULL);CHKERRQ(ierr); 
	
	/* velocity cm/y */    
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_vx",&vx,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_vz",&vz,NULL);CHKERRQ(ierr); 
	
	/* rho0 for initial pressure */
	ierr = PetscOptionsGetReal(NULL,"-model_rift3D_rho0",&data->rho0,NULL);CHKERRQ(ierr); 
	
	/* Material constant */
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
		PetscPrintf(PETSC_COMM_WORLD,"reading options");
		ierr= MaterialConstantsSetFromOptions(materialconstants,"model_rift3D",regionidx,PETSC_FALSE);CHKERRQ(ierr);
	}
	
	/*Compute velocity at bottom*/
	Sx = (data->Ly - data->Oy)*(data->Lz - data->Oz);
	Sz = (data->Ly - data->Oy)*(data->Lx - data->Ox);
	Sy = (data->Lx - data->Ox)*(data->Lz - data->Oz);
	vy = (vx*Sx+vz*Sz)*2/Sy;
	
	/* reports before scaling */
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_Ox %+1.4e [SI] -model_rift3D_Lx : %+1.4e [SI]\n", data->Ox ,data->Lx );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_Oy %+1.4e [SI] -model_rift3D_Ly : %+1.4e [SI]\n", data->Oy ,data->Ly );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_Oz %+1.4e [SI] -model_rift3D_Lz : %+1.4e [SI]\n", data->Oz ,data->Lz );
	PetscPrintf(PETSC_COMM_WORLD,"  -model_rift3D_vx [m/s]:  %+1.4e  -model_rift3D_vz [m/s]:  %+1.4e : computed vy [m/s]:  %+1.4e \n", vx,vz,vy);
	PetscPrintf(PETSC_COMM_WORLD,"-model_rift3D_rho0 [kg/m^3] :%+1.4e \n", data->rho0 );
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) { 
		MaterialConstantsPrintAll(materialconstants,regionidx);
	} 
	
	if (data->dimensional) {
		/*Compute additional scaling parameters*/
		data->time_bar      = data->length_bar / data->velocity_bar;
		data->pressure_bar  = data->viscosity_bar*data->velocity_bar / data->length_bar;
		data->density_bar   = data->viscosity_bar*data->velocity_bar / ( data->length_bar * data->length_bar );
		
		PetscPrintf(PETSC_COMM_WORLD,"[rift3D]:  during the solve scaling will be done using \n");
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
		
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_Ox   :  %+1.4e    -model_rift3D_Lx   :  %+1.4e  \n", data->Ox ,data->Lx );
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_Oy   :  %+1.4e    -model_rift3D_Ly   :  %+1.4e \n", data->Oy, data->Ly );
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_Oz   :  %+1.4e    -model_rift3D_Lz   :  %+1.4e\n", data->Oz , data->Lz );
		
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_Vx:%+1.4e    -model_rift3D_vy:%+1.4e    -model_rift3D_vz:  %+1.4e \n", data->vx ,data->vy, data->vz);
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_rho0:%+1.4e \n", data->rho0 );
		PetscPrintf(PETSC_COMM_WORLD,"scaled value for material parameters\n");
		for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) { 
			MaterialConstantsPrintAll(materialconstants,regionidx);
		}    
	}
	
	PetscFunctionReturn(0);
}


/* 
 
 1/ Define boundary conditions in one place for this model.
 
 2/ Calling pattern should always be
 PetscErrorCode ModelRift3D_DefineBCList(BCList bclist,DM dav,pTatinCtx user,ModelRift3DCtx data)
 where ModelRift3DCtx data is a different type for each model.
 
 3/ Re-use this function in 
 ModelApplyBoundaryCondition_Rift3D();
 ModelApplyBoundaryConditionMG_Rift3D();
 
 */
#undef __FUNCT__
#define __FUNCT__ "ModelRift3D_DefineBCList"
PetscErrorCode ModelRift3D_DefineBCList(BCList bclist,DM dav,pTatinCtx user,ModelRift3DCtx *data)
{
	PetscScalar    vxl,vxr,vzf,vzb,vy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	vxl = -data->vx;
	vxr =  data->vx;
	vy  =  data->vy;
	vzf = -data->vz;
	vzb =  data->vz;
	
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
#define __FUNCT__ "ModelApplyBoundaryCondition_Rift3D"
PetscErrorCode ModelApplyBoundaryCondition_Rift3D(pTatinCtx user,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = ModelRift3D_DefineBCList(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,user,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Rift3D"
PetscErrorCode ModelApplyBoundaryConditionMG_Rift3D(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscInt       n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = ModelRift3D_DefineBCList(bclist[n],dav[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Rift3D"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D(pTatinCtx c,void *ctx)
{
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Rift3D"
PetscErrorCode ModelApplyInitialMeshGeometry_Rift3D(pTatinCtx c,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,data->Ox,data->Lx,data->Oy,data->Ly,data->Oz,data->Lz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Rift3D"
PetscErrorCode ModelApplyInitialMaterialGeometry_Rift3D(pTatinCtx c,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	int                    p,n_mp_points;
	PetscScalar            y_lab,y_moho,y_midcrust;
	DataBucket             db;
	DataField              PField_std,PField_stokes,PField_DensConst,PField_pls,PField_ViscZ;
	int                    phase;
	DataBucket materialconstants = c->material_constants;
	MaterialConst_DensityConst *DensConst_data;
	MaterialConst_ViscosityZ *ViscZ_data;
			
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* define properties on material points */
	db = c->materialpoint_db;
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));
	
	DataBucketGetDataFieldByName(db,MPntPStokesPl_classname,&PField_pls);
	DataFieldGetAccess(PField_pls);
	DataFieldVerifyAccess(PField_pls,sizeof(MPntPStokesPl));
	
	DataBucketGetDataFieldByName(materialconstants,MaterialConst_DensityConst_classname,  &PField_DensConst);
	DensConst_data    = (MaterialConst_DensityConst*)PField_DensConst->data;
	
	DataBucketGetDataFieldByName(materialconstants,MaterialConst_ViscosityZ_classname,  &PField_ViscZ);
	ViscZ_data    = (MaterialConst_ViscosityZ*)PField_ViscZ->data;
	
	/* m */
	y_lab      = -120.0e3; 
	y_moho     = -35.0e3;
	y_midcrust = -20.0e3;
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd       *material_point;
		MPntPStokes   *mpprop_stokes;
		MPntPStokesPl *mpprop_pls;
		double        *position,ycoord,xcoord,zcoord;
		double        rho;
		float         pls;
		char          yield;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		DataFieldAccessPoint(PField_pls,p,(void**)&mpprop_pls);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		/* convert to scaled units */
		ycoord = position[1] * data->length_bar;
		xcoord=position[0] * data->length_bar;
		zcoord=position[2] * data->length_bar;
		y_moho     = -30.0e3;
		/*if (abs(data->Lx * data->length_bar/2 - xcoord)<40.e3  & abs(data->Lz * data->length_bar/2 - zcoord) < 40e3){ 
			y_moho     = -50.0e3;
		}*/
		if (ycoord<y_lab) {
			phase = 3;
			rho   = DensConst_data[3].density;
		} else if (ycoord<y_moho) {
			phase = 2;
			rho   = DensConst_data[2].density;
            //if (abs(data->Lx * data->length_bar/2 - xcoord)<40.e3  & abs(data->Lz * data->length_bar/2 - zcoord  ) < 40e3 & ycoord > -50.0e3){
            if ((zcoord<150.e3)  && (abs(data->Lx * data->length_bar/2 - xcoord  ) < 40e3) && (ycoord > -50.0e3)){    
            phase = 4;
            rho   = DensConst_data[4].density;   
            }
		} else if (ycoord<y_midcrust) {
			phase = 1;
			rho   = DensConst_data[1].density;
		} else {
			phase = 0;
			rho   = DensConst_data[0].density;
		}
		
		//eta = 1.0; 
		pls = 0.0;
		yield = 0; 
		rho = -rho * 10.0;
		
		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,phase);
		MPntPStokesPlSetField_yield_indicator(mpprop_pls,yield);
		MPntPStokesPlSetField_plastic_strain(mpprop_pls,pls);
		MPntPStokesSetField_density(mpprop_stokes,rho);
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	DataFieldRestoreAccess(PField_pls);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Rift3D"
PetscErrorCode ModelApplyUpdateMeshGeometry_Rift3D(pTatinCtx c,Vec X,void *ctx)
{
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Rift3D_CheckScales"
PetscErrorCode ModelOutput_Rift3D_CheckScales(pTatinCtx c,Vec X)
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
	
	PetscPrintf(PETSC_COMM_WORLD,"[rift3D]: check scales \n");
	
	ierr = VecCopy(X,Xcopy);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,Xcopy,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = VecStrideMin(pressure,0,NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," min|P0|   = %+1.4e \n",fp);
	ierr = VecStrideMax(pressure,0,NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," max|P0|   = %+1.4e \n",fp);
	
	ierr = VecStrideMin(pressure,1,NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," min|dPdx| = %+1.4e \n",fp);
	ierr = VecStrideMax(pressure,1,NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," max|dPdx| = %+1.4e \n",fp);
	
	ierr = VecStrideMin(pressure,2,NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," min|dPdy| = %+1.4e \n",fp);
	ierr = VecStrideMax(pressure,2,NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," max|dPdy| = %+1.4e \n",fp);
	
	ierr = VecStrideMin(pressure,3,NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," min|dPdz| = %+1.4e \n",fp);
	ierr = VecStrideMax(pressure,3,NULL,&fp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," max|dPdz| = %+1.4e \n",fp);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,Xcopy,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(Xcopy);CHKERRQ(ierr);
	ierr = FormFunction_Stokes(NULL,Xcopy,RHS,(void*)c);CHKERRQ(ierr);
	
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
	
	ierr = FormFunction_Stokes(NULL,Xcopy,F,(void*)c);CHKERRQ(ierr);
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
	
	ierr = FormFunction_Stokes(NULL,Xcopy,F,(void*)c);CHKERRQ(ierr);
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
#define __FUNCT__ "ModelOutput_Rift3D"
PetscErrorCode ModelOutput_Rift3D(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	//ierr = ModelOutput_Rift3D_CheckScales(c,X);CHKERRQ(ierr);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	
	{
		//  Write out just the stokes variable?
		//  const int nf = 1;
		//  const MaterialPointField mp_prop_list[] = { MPField_Stokes };
		//
		//  Write out just std, stokes and plastic variables
		const int nf = 3;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
		char mp_file_prefix[256];
		
		sprintf(mp_file_prefix,"%s_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(c->materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Rift3D"
PetscErrorCode ModelDestroy_Rift3D(pTatinCtx c,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialStokesVariableMarkers_Rift3D"
PetscErrorCode ModelApplyInitialStokesVariableMarkers_Rift3D(pTatinCtx user,Vec X,void *ctx)
{
    DM                stokes_pack,dau,dap;
	PhysCompStokes    stokes;
    Vec               Uloc,Ploc;
    PetscScalar       *LA_Uloc,*LA_Ploc;
    ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscErrorCode    ierr;
    PetscInt regionidx;	
    PetscFunctionBegin;
    
    
    PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    if (~data->runmises) {
        for (regionidx=0; regionidx<user->rheology_constants.nphases_active;regionidx++) {
            MaterialConstantsSetValues_MaterialType(user->material_constants,regionidx,VISCOUS_Z,PLASTIC_MISES,SOFTENING_NONE,DENSITY_CONSTANT);
        }
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

    if (~data->runmises) {
        for (regionidx=0; regionidx<user->rheology_constants.nphases_active;regionidx++) {
            MaterialConstantsSetValues_MaterialType(user->material_constants,regionidx,VISCOUS_Z,PLASTIC_DP,SOFTENING_NONE,DENSITY_CONSTANT);
        }
    }
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialCondition_Rift3D"
PetscErrorCode ModelApplyInitialCondition_Rift3D(pTatinCtx c,Vec X,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	DM stokes_pack,dau,dap;
	Vec velocity,pressure;
	PetscReal vxl,vxr,vzb,vzf,vy;
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	DMDAVecTraverse3d_InterpCtx IntpCtx;
	PetscReal MeshMin[3],MeshMax[3],domain_height;
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
	vzb =  data->vz;
	
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
	
    ierr = DMDAVecTraverseIJK(dap,pressure,0,DMDAVecTraverseIJK_HydroStaticPressure_v2,     (void*)&HPctx); /* P = P0 + a.x + b.y + c.z, modify P0 (idx=0) */
    ierr = DMDAVecTraverseIJK(dap,pressure,2,DMDAVecTraverseIJK_HydroStaticPressure_dpdy_v2,(void*)&HPctx); /* P = P0 + a.x + b.y + c.z, modify b  (idx=2) */
    ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,"testHP");CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Rift3D"
PetscErrorCode pTatinModelRegister_Rift3D(void)
{
	ModelRift3DCtx *data;
	pTatinModel m;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelRift3DCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelRift3DCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"rift3D");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialCondition_Rift3D);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,   (void (*)(void))ModelApplyInitialStokesVariableMarkers_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Rift3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Rift3D);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

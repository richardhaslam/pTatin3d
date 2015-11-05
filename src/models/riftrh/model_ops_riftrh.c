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
 **    filename:   model_ops_riftrh.c
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
#include "material_point_std_utils.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPStokesPl_def.h"
#include "ptatin3d_energy.h"
#include "energy_output.h"
#include "model_utils.h"
#include "dmda_iterator.h"
#include "mesh_update.h"

#include "ptatin_models.h"

#include "model_riftrh_ctx.h"

PetscErrorCode ModelApplyUpdateMeshGeometry_Riftrh_semi_eulerian(pTatinCtx c,Vec X,void *ctx);
PetscErrorCode ModelApplyMaterialBoundaryCondition_Riftrh_semi_eulerian(pTatinCtx c,void *ctx);
PetscBool BCListEvaluator_riftrhl( PetscScalar position[], PetscScalar *value, void *ctx );
PetscBool BCListEvaluator_riftrhr( PetscScalar position[], PetscScalar *value, void *ctx );

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Riftrh"
PetscErrorCode ModelInitialize_Riftrh(pTatinCtx c,void *ctx)
{
	ModelRiftrhCtx *data = (ModelRiftrhCtx*)ctx;
	PetscBool flg;
	RheologyConstants *rheology;
	DataBucket materialconstants = c->material_constants;
	PetscBool nondim,use_energy;
	PetscScalar dh_vx;
	PetscInt regionidx;
	PetscReal cm_per_yer2m_per_sec = 1.0e-2 / ( 365.0 * 24.0 * 60.0 * 60.0 ),phi1_rad,phi2_rad ;
	PetscReal preexpA,Ascale,entalpy,Vmol,nexp,Tref;
	PetscInt nlayers;
	
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Temperature flag */
	use_energy=PETSC_TRUE;
	//    use_energy=PETSC_FALSE;
	
	/* set the deffault values of the material constant for this particular model */
	/*scaling */
	data->length_bar    = 100.0 * 1.0e3;
	data->viscosity_bar = 1e26;
	data->velocity_bar  = 1.0e-10;
	data->dimensional   = PETSC_TRUE;
	
	/* Box GEOMETRY */
	data->Lx = 1200.0e3;
	data->Ly = 250.0e3;
	data->Lz = 1200.0e3;
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_Lx",&data->Lx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_Ly",&data->Ly,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_Lz",&data->Lz,&flg);CHKERRQ(ierr);
	
	data->hc = 60.0e3;
	data->hm = 60.0e3;
	data->ha = 130.0e3;
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_hc",&data->hc,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_hm",&data->hm,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_ha",&data->ha,&flg);CHKERRQ(ierr);
	
	/* seed geometry */
	data->dxs = 12.0e3;
	data->dys = 6.0e3;
	data->dzs = 1200.0e3;
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_dxs",&data->dxs,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_dys",&data->dys,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_dzs",&data->dzs,&flg);CHKERRQ(ierr);
	
	/* velocity boundary condition geometry */
	data->hvbx1 = 125.0e3;
	data->hvbx2 = 115.0e3;
	data->vx_up = 0.5*cm_per_yer2m_per_sec;
	/* VELOCITY BOUNDARY CONDITION GEOMETRY */
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_hvbx1",&data->hvbx1,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_hvbx2",&data->hvbx2,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_vx_up",&data->vx_up,&flg);CHKERRQ(ierr);
	
	/* material propoerties */
	data->rhoc = 2800.0;
	data->rhom = 3300.0;
	data->rhoa = 3250.0;
	data->etac = 1.0e28;
	data->etam = 1.0e21;
	data->etaa = 1.0e19;
	/* MATERIAL PARAMETERS */
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_rhoc",&data->rhoc,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_rhom",&data->rhom,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_rhoa",&data->rhoa,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_etac",&data->etac,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_etam",&data->etam,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_etaa",&data->etaa,&flg);CHKERRQ(ierr);
	
	/*Temperature */
	if(use_energy) {
		data->thermalparams.nlayers  = 3;
		data->thermalparams.lscale   = data->length_bar;
		
		/*         data->thermalparams.ytop[0]  = 250.0e3;
		 data->thermalparams.thick[0] = 120.0e3;
		 data->thermalparams.ttop[0]  = 0.0;
		 data->thermalparams.cond[0]  = 2.25;
		 data->thermalparams.hp[0]    = 0.0;
		 data->thermalparams.qbase[0] = 22.5e-3;
		 
		 data->thermalparams.ytop[1]  = 130.0e3;
		 data->thermalparams.thick[1] = 130.0e3;
		 data->thermalparams.ttop[1]  = 1200.0;
		 data->thermalparams.cond[1]  = 2.25;
		 data->thermalparams.hp[1]    = 0.0;
		 data->thermalparams.qbase[1] = 0.e-3;
		 */
		data->thermalparams.ytop[0]  = 250.0e3;
		data->thermalparams.thick[0] = 60.0e3;
		data->thermalparams.ttop[0]  = 0.0;
		data->thermalparams.cond[0]  = 2.25;
		data->thermalparams.hp[0]    = 0.0;
		data->thermalparams.qbase[0] = 22.5e-3;
		
		data->thermalparams.ytop[1]  = 190.0e3;
		data->thermalparams.thick[1] = 60.0e3;
		data->thermalparams.ttop[1]  = 600.0;
		data->thermalparams.cond[1]  = 2.25;
		data->thermalparams.hp[1]    = 0.0;
		data->thermalparams.qbase[1] = 22.5e-3;
		
		data->thermalparams.ytop[2]  = 130.0e3;
		data->thermalparams.thick[2] = 130.0e3;
		data->thermalparams.ttop[2]  = 1200.0;
		data->thermalparams.cond[2]  = 2.25;
		data->thermalparams.hp[2]    = 0.0;
		data->thermalparams.qbase[2] = 0.0e-3;
		
		nlayers = data->thermalparams.nlayers;
		data->Ttop = data->thermalparams.ttop[nlayers-1];
		data->Tbottom = data->thermalparams.ttop[nlayers-1] + ( data->thermalparams.hp[nlayers-1]* pow(data->thermalparams.thick[nlayers-1],2) / 2.0 + data->thermalparams.qbase[nlayers-1] * data->thermalparams.thick[nlayers-1] ) / data->thermalparams.cond[nlayers-1];
		data->Ttop=0;
		data->Tbottom=1200;
	}
	
	/* rheology parameters */
	rheology                = &c->rheology_constants;
	rheology->rheology_type = RHEOLOGY_VP_STD;
	rheology->nphases_active = 3;
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global = 1.e+28;
	rheology->eta_lower_cutoff_global = 1.e+19;
	data->runmises = PETSC_FALSE;
	
	//    ierr = PetscOptionsGetInt(NULL,"-model_Riftrh_param2",&data->param2,&flg);CHKERRQ(ierr);
	
	/* Material constants */
	MaterialConstantsSetDefaults(materialconstants);
	
	/* PHASE 0, ASTHENOSPHERE */
	regionidx = 0;
	if(use_energy == PETSC_FALSE){
		//    if(PETSC_TRUE){
		MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,
																						DENSITY_CONSTANT);
	}
	else {
		MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	}
	
	//VISCOSITY PARAMETERS - wet olivine karato
	preexpA  = 1.393e-14;
	Ascale  = 1.;
	entalpy = 430.e3;
	Vmol    = 15.e-6;
	nexp    = 3.;
	Tref    = 273.;
	MaterialConstantsSetValues_ViscosityConst(materialconstants,regionidx,data->etaa);
	MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexpA,Ascale,entalpy,Vmol,nexp,Tref);
	
	//DENSITY PARAMETERS
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,regionidx,data->rhom,2.e-5,0.);
	MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,data->rhoa);
	
	//PLASTICITY PARAMETERS
	phi1_rad=M_PI*15./180.;
	phi2_rad=M_PI*2./180.;
	MaterialConstantsSetValues_PlasticDP(materialconstants,regionidx,phi1_rad,phi2_rad,2.e7,2.e7,1.e7,1e20);
	MaterialConstantsSetValues_PlasticMises(materialconstants,regionidx,1.e8,1.e8);
	MaterialConstantsSetValues_SoftLin(materialconstants,regionidx,0.1,0.5);
	
	/* PHASE 1, MANTLE LITHOSPHERE */
	regionidx = 1;
	if(use_energy == PETSC_FALSE){
		//    if(PETSC_TRUE){
		MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,
																						DENSITY_CONSTANT);
	}
	else {
		MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	}
	
	//VISCOSITY PARAMETERS - wet olivine x 5, karato
	preexpA = 1.393e-14;
	Ascale  = 5.;
	entalpy = 430.e3;
	Vmol    = 0.e-6;
	nexp    = 3.;
	Tref    = 273.;
	MaterialConstantsSetValues_ViscosityConst(materialconstants,regionidx,data->etam);
	MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexpA,Ascale,entalpy,Vmol,nexp,Tref);
	
	//DENSITY PARAMETERS
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,regionidx,data->rhom,2.e-5,0.0);
	MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,data->rhom);
	
	//PLASTICITY PARAMETERS
	phi1_rad=M_PI*15./180.;
	phi2_rad=M_PI*2./180.;
	MaterialConstantsSetValues_PlasticDP(materialconstants,regionidx,phi1_rad,phi2_rad,2.e7,2.e7,1.e7,1e20);
	MaterialConstantsSetValues_PlasticMises(materialconstants,regionidx,1.e8,1.e8);
	MaterialConstantsSetValues_SoftLin(materialconstants,regionidx,0.1,0.5);
	
	/* PHASE 2, CRUST / UPPER LITHOSPHERE */
	regionidx = 2;
	if(use_energy == PETSC_FALSE){
		//    if(PETSC_TRUE){
		MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,
																						DENSITY_CONSTANT);
	}
	else {
		MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	}
	
	//VISCOSITY PARAMETERS
	preexpA = 1.393e-14;
	Ascale  = 5.;
	entalpy = 430.e3;
	Vmol    = 0.e-6;
	nexp    = 3.;
	Tref    = 273.;
	MaterialConstantsSetValues_ViscosityConst(materialconstants,regionidx,data->etac);
	MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexpA,Ascale,entalpy,Vmol,nexp,Tref);
	
	//DENSITY PARAMETERS
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,regionidx,data->rhoc,2.e-5,0.0);
	MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,data->rhoc);
	
	//PLASTICITY PARAMETERS
	phi1_rad=M_PI*15./180.;
	phi2_rad=M_PI*2./180.;
	MaterialConstantsSetValues_PlasticDP(materialconstants,regionidx,phi1_rad,phi2_rad,2.e7,2.e7,1.e7,1e20);
	MaterialConstantsSetValues_PlasticMises(materialconstants,regionidx,1.e8,1.e8);
	MaterialConstantsSetValues_SoftLin(materialconstants,regionidx,0.1,0.5);
	
	
	/* Read the options */
	/* cutoff */
	ierr = PetscOptionsGetBool(NULL,"-model_Riftrh_apply_viscosity_cutoff_global",&rheology->apply_viscosity_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_eta_lower_cutoff_global",&rheology->eta_lower_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_eta_upper_cutoff_global",&rheology->eta_upper_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-model_Riftrh_runwithmises",&data->runmises,NULL);CHKERRQ(ierr);
	/* scaling */
	nondim = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-model_Riftrh_nondimensional",&nondim,NULL);CHKERRQ(ierr);
	if (nondim){
		data->dimensional = PETSC_FALSE;
	} else {
		ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_vis_bar",&data->viscosity_bar,NULL);CHKERRQ(ierr);
		ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_vel_bar",&data->velocity_bar,NULL);CHKERRQ(ierr);
		ierr = PetscOptionsGetReal(NULL,"-model_Riftrh_length_bar",&data->length_bar,NULL);CHKERRQ(ierr);
	}
	
	/* compute vxdown */
	dh_vx = data->hvbx1 - data->hvbx2;
	data->vx_down = -data->vx_up * (data->hc + data->hm +dh_vx/2) / (data->ha - dh_vx/2);
	
	/* reports before scaling */
	PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  input: -model_Riftrh_Lx : %+1.4e [SI]\n", data->Lx );
	PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  input: -model_Riftrh_Ly : %+1.4e [SI]\n", data->Ly );
	PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  input: -model_Riftrh_Lz : %+1.4e [SI]\n", data->Lz );
	PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  -model_Riftrh_vx_up [m/s]:  %+1.4e \n", data->vx_up);
	PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  -model_Riftrh_rhoc [kg/m^3] :%+1.4e \n", data->rhoc );
	PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  -model_Riftrh_rhom [kg/m^3] :%+1.4e \n", data->rhom );
	PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  -model_Riftrh_rhoa [kg/m^3] :%+1.4e \n", data->rhoa );
	
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
		MaterialConstantsPrintAll(materialconstants,regionidx);
	}
	
	if (data->dimensional) {
		/* Compute additional scaling parameters */
		data->time_bar      = data->length_bar / data->velocity_bar;
		data->pressure_bar  = data->viscosity_bar/data->time_bar;
		data->density_bar   = data->pressure_bar / data->length_bar;
		
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh]:  during the solve scaling will be done using \n");
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  L*    : %1.4e [m]\n", data->length_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  U*    : %1.4e [m.s^-1]\n", data->velocity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  t*    : %1.4e [s]\n", data->time_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  eta*  : %1.4e [Pa.s]\n", data->viscosity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  rho*  : %1.4e [kg.m^-3]\n", data->density_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh]  P*    : %1.4e [Pa]\n", data->pressure_bar );
		//scale viscosity cutoff
		rheology->eta_lower_cutoff_global = rheology->eta_lower_cutoff_global / data->viscosity_bar;
		rheology->eta_upper_cutoff_global = rheology->eta_upper_cutoff_global / data->viscosity_bar;
		//scale length
		data->Lx = data->Lx / data->length_bar;
		data->Ly = data->Ly / data->length_bar;
		data->Lz = data->Lz / data->length_bar;
		data->hc = data->hc / data->length_bar;
		data->hm = data->hm / data->length_bar;
		data->ha = data->ha / data->length_bar;
		/* seed geometry */
		data->dxs = data->dxs / data->length_bar;
		data->dys = data->dys / data->length_bar;
		data->dzs = data->dzs / data->length_bar;
		//scale velocity
		data->vx_up   = data->vx_up  /data->velocity_bar;
		data->vx_down = data->vx_down/data->velocity_bar;
		data->hvbx1   = data->hvbx1/data->length_bar;
		data->hvbx2   = data->hvbx2/data->length_bar;
		//scale rho
		data->rhoc = data->rhoc/data->density_bar;
		data->rhom = data->rhom/data->density_bar;
		data->rhoa = data->rhoa/data->density_bar;
		
		// scale material properties
		for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
			MaterialConstantsScaleAll(materialconstants,regionidx,data->length_bar,data->velocity_bar,data->time_bar,data->viscosity_bar,data->density_bar,data->pressure_bar);
		}
		
		/* Reports scaled values */		
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] scaled value    -model_Riftrh_Lx   :  %+1.4e \n", data->Lx );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] scaled value    -model_Riftrh_Ly   :  %+1.4e \n", data->Ly );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] scaled value    -model_Riftrh_Lz   :  %+1.4e \n", data->Lz );
		
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] scaled value   -model_Riftrh_vx_up:%+1.4e    -model_Riftrh_vx_down:%+1.4e \n", data->vx_up ,data->vx_down);
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] scaled value   -model_Riftrh_rhoc:%+1.4e \n", data->rhoc );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] scaled value   -model_Riftrh_rhom:%+1.4e \n", data->rhom );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] scaled value   -model_Riftrh_rhoa:%+1.4e \n", data->rhoa );
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] scaled value for material parameters\n");
		for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
			MaterialConstantsPrintAll(materialconstants,regionidx);
		}
	}
	
	data->output_markers = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-model_Riftrh_output_markers",&data->output_markers,NULL);CHKERRQ(ierr);
	
	/* USE ENERGY EQUATION */
	if(use_energy){
		ierr = PetscOptionsInsertString("-activate_energy");CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

/* SET AN INITIAL BACK GROUND STRAIN RATE, TEMPERATURE, PRESSURE */
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_Riftrh"
PetscErrorCode ModelApplyInitialSolution_Riftrh(pTatinCtx c,Vec X,void *ctx)
{
	ModelRiftrhCtx *data = (ModelRiftrhCtx*)ctx;
	PetscErrorCode ierr;
	PetscBool use_energy;
	PhysCompStokes     stokes;
	DM stokes_pack,dau,dap;
	DMDAVecTraverse3d_InterpCtx IntpCtx;
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	Vec velocity,pressure;
	PetscReal MeshMax[3],MeshMin[3],height,length,vxl,vxr,vybottom;
	PetscBool use_initial_up_field = PETSC_FALSE;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscOptionsGetBool(NULL,"-model_Riftrh_use_initial_up_field",&use_initial_up_field,NULL);
	
	if (use_initial_up_field) {
		PetscPrintf(PETSC_COMM_WORLD,"[riftrh] Using velocity from boundary condition and mantle hydrostatic pressure as initial condition\n");
		
		ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
		stokes_pack = stokes->stokes_pack;
		ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
		ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		
		/* velocity intial condition - background strain */
		ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
		vxl = -data->vx_up;
		vxr =  data->vx_up;
		ierr = DMDAGetBoundingBox(dau,MeshMin,MeshMax);CHKERRQ(ierr);
		height = MeshMax[1] - MeshMin[1];
		length = MeshMax[0] - MeshMin[0];
		vybottom = 2.*data->vx_up * height / length;
		
		ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vxr-vxl)/(length),vxl,0.0);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dau,velocity,0,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
		
		ierr = DMDAVecTraverse3d_InterpCtxSetUp_Y(&IntpCtx,(0.0-vybottom)/(height),vybottom,0.0);CHKERRQ(ierr);
		ierr = DMDAVecTraverse3d(dau,velocity,1,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
		
		/* pressure initial condition - hydrostatic */
		ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
		HPctx.surface_pressure = 0.0;
		HPctx.ref_height = height;
		HPctx.ref_N      = c->stokes_ctx->my-1;
		HPctx.grav       = 10.0;
		HPctx.rho        = data->rhom;
		ierr = DMDAVecTraverseIJK(dap,pressure,0,DMDAVecTraverseIJK_HydroStaticPressure_v2,     (void*)&HPctx);CHKERRQ(ierr); /* P = P0 + a.x + b.y + c.z, modify P0 (idx=0) */
		ierr = DMDAVecTraverseIJK(dap,pressure,2,DMDAVecTraverseIJK_HydroStaticPressure_dpdy_v2,(void*)&HPctx);CHKERRQ(ierr); /* P = P0 + a.x + b.y + c.z, modify b  (idx=2) */
		
		ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	}
	
	ierr = pTatinContextValid_Energy(c,&use_energy);CHKERRQ(ierr);
	if (use_energy) {
		PhysCompEnergy energy;
		Vec            temperature;
		DM             daT;
		
		ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);
		daT  = energy->daT;
		
		ierr = DMDAVecTraverse3d(daT,temperature,0,DMDAVecTraverse_InitialThermalField3D,(void*)&data->thermalparams);CHKERRQ(ierr);
	}
	//    exit(1);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelRiftrh_DefineBCList"
PetscErrorCode ModelRiftrh_DefineBCList(BCList bclist,DM dav,pTatinCtx user,ModelRiftrhCtx *data)
{
	PetscScalar    vxl,vxr,vybottom,zero,height,length;
	PetscInt       vbc_type;
	PetscErrorCode ierr;
	PetscReal MeshMin[3],MeshMax[3];
	DM stokes_pack,dau,dap;
	
	
	PetscFunctionBegin;
	
	stokes_pack = user->stokes_ctx->stokes_pack;
	
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	
	zero=0.;
	//    vbc_type = 1; /* in / out flow condition on the sides */
	vbc_type = 2; /* outflow condition on the sides, inflow condition on the base */
	
	if (vbc_type == 1) {
		
		/* infilling free slip base */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		/* free surface top*/
		
		/*extension along face of normal x */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_riftrhl,(void*)data);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_riftrhr,(void*)data);CHKERRQ(ierr);
		
		/*free slip base*/
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		/* no flow in z*/
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}
	
	if (vbc_type == 2) {
		vxl = -data->vx_up;
		vxr = data->vx_up;
		ierr = DMDAGetBoundingBox(dau,MeshMin,MeshMax);CHKERRQ(ierr);
		height = MeshMax[1] - MeshMin[1];
		length = MeshMax[0] - MeshMin[0];
		vybottom = 2.*data->vx_up * height / length;
		
		/* infilling free slip base */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		/* free surface top*/
		
		/*extension along face of normal x */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&vxl);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&vxr);CHKERRQ(ierr);
		
		/*free slip base*/
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&vybottom);CHKERRQ(ierr);
		
		/* no flow in z*/
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_riftrhl"
PetscBool BCListEvaluator_riftrhl( PetscScalar position[], PetscScalar *value, void *data )
{
	PetscBool impose_dirichlet = PETSC_TRUE;
	PetscScalar vx,dh_vx;
	ModelRiftrhCtx *datal = (ModelRiftrhCtx*)data;
	
	PetscFunctionBegin;
	
	dh_vx = datal->hvbx1 - datal->hvbx2;
	if ( position[1] >= datal->ha ) {
		vx = -datal->vx_up;
	} else if ( position[1] < datal->ha - dh_vx) {
		vx = -datal->vx_down;
	} else {
		vx = -(datal->vx_up + (datal->vx_down - datal->vx_up) * (datal->hvbx1 - position[1]) / (datal->hvbx1 - datal->hvbx2) );
	}
	*value = vx;
	return impose_dirichlet;
}

#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_riftrhr"
PetscBool BCListEvaluator_riftrhr( PetscScalar position[], PetscScalar *value, void *data )
{
	PetscBool impose_dirichlet = PETSC_TRUE;
	PetscScalar vx,dh_vx;
	ModelRiftrhCtx *datal = (ModelRiftrhCtx*)data;
	
	PetscFunctionBegin;
	
	dh_vx = datal->hvbx1 - datal->hvbx2;
	if ( position[1] >= datal->ha ) {
		vx = datal->vx_up;
	} else if ( position[1] < datal->ha - dh_vx) {
		vx = datal->vx_down;
	} else {
		vx = (datal->vx_up + (datal->vx_down - datal->vx_up) * (datal->hvbx1 - position[1]) / (datal->hvbx1 - datal->hvbx2) );
	}
	*value = vx;
	return impose_dirichlet;
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Riftrh"
PetscErrorCode ModelApplyBoundaryCondition_Riftrh(pTatinCtx c,void *ctx)
{
	ModelRiftrhCtx *data = (ModelRiftrhCtx*)ctx;
	PetscBool use_energy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = ModelRiftrh_DefineBCList(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,c,data);CHKERRQ(ierr);
	
	ierr = pTatinContextValid_Energy(c,&use_energy);CHKERRQ(ierr);
	if (use_energy) {
		PetscReal      val_T;
		PhysCompEnergy energy;
		BCList         bclist;
		DM             daT;
		
		ierr   = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		daT    = energy->daT;
		bclist = energy->T_bclist;
		
		val_T = data->Ttop;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);
		
		val_T = data->Tbottom;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Riftrh"
PetscErrorCode ModelApplyBoundaryConditionMG_Riftrh(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	PetscInt n;
	ModelRiftrhCtx *data = (ModelRiftrhCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = ModelRiftrh_DefineBCList(bclist[n],dav[n],user,data);CHKERRQ(ierr);
		
	}	
	
	PetscFunctionReturn(0);
}

// adding particles on the lower boundary to accommodate inflow
// adding particles on the left and right boundary to accommodate inflow
#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Riftrh_semi_eulerian"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Riftrh_semi_eulerian(pTatinCtx c,void *ctx)
{
#if 0
	ModelRiftrhCtx     *data = (ModelRiftrhCtx*)ctx;
	PhysCompStokes     stokes;
	DM                 stokes_pack,dav,dap;
	PetscInt           Nxp[2];
	PetscReal          perturb;
	DataBucket         material_point_db,material_point_face_db;
	PetscInt           f, n_face_list=3, face_list[] = { 0, 1, 3 }; // xmin, xmax, ybase //
	int                p,n_mp_points;
	MPAccess           mpX;
	PetscErrorCode     ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
		
	ierr = pTatinGetMaterialPoints(c,&material_point_db,NULL);CHKERRQ(ierr);
		
	/* create face storage for markers */
	DataBucketDuplicateFields(material_point_db,&material_point_face_db);
	
	for (f=0; f<n_face_list; f++) {
		
		/* traverse */
		/* [0,1/east,west] ; [2,3/north,south] ; [4,5/front,back] */
		Nxp[0]  = 1;
		Nxp[1]  = 1;
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
	
	/* delete */
	DataBucketDestroy(&material_point_face_db);
#endif	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Riftrh"
PetscErrorCode ModelApplyInitialMeshGeometry_Riftrh(pTatinCtx c,void *ctx)
{
	ModelRiftrhCtx *data = (ModelRiftrhCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,0.0,data->Lx,0.0,data->Ly,0.0,data->Lz);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"[riftrh] Lx = %1.4e \n", data->Lx );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Riftrh"
PetscErrorCode ModelApplyInitialMaterialGeometry_Riftrh(pTatinCtx c,void *ctx)
{
	ModelRiftrhCtx *data = (ModelRiftrhCtx*)ctx;
	int                    p,n_mp_points;
	PetscInt               notch_type;
	DataBucket             db;
	DataField              PField_std,PField_stokes,PField_pls;
	int                    phase;
	PetscScalar            ha_dimensional,hm_dimensional,hc_dimensional,notch_height,notch_width,x_center,y_center,z_center;
	PetscScalar            xp_dimensional,yp_dimensional,zp_dimensional;
	PetscErrorCode ierr;
	MPAccess           mpX;
	PetscBool use_energy;
	
	
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
	
	ha_dimensional = data->ha * data->length_bar;
	hm_dimensional = data->hm * data->length_bar;
	hc_dimensional = data->hc * data->length_bar;
	x_center       = 0.5 * data->Lx * data->length_bar;
	y_center       = 0.5 * data->Ly * data->length_bar;
	z_center       = 0.5 * data->Lz * data->length_bar;
	notch_width    = 50.0e3;
	notch_height   = 30.0e3;
	
	/* marker loop */
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	/*define notch type, 1: one notch; 2 two notches/weak seeds */
	notch_type = 2; 
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point_std;
		MPntPStokes *material_point_stokes;
		MPntPStokesPl *mpprop_pls;
		double      *position;
		double      eta,rho;
		int         phase;
		float       plastic_strain;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point_std);
		DataFieldAccessPoint(PField_stokes,p,(void**)&material_point_stokes);
		DataFieldAccessPoint(PField_pls,p,(void**)&mpprop_pls);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point_std,&position);
		
		phase = 0;
		eta=data->etaa;
		rho=data->rhoa;
		if (position[1] > data->ha && position[1] < (data->ha+data->hm) ) {
			phase = 1;
			eta=data->etam;
			rho=data->rhom;
		}
		if (position[1]>data->ha+data->hm ) {
			phase = 2;
			eta=data->etac;
			rho=data->rhoc;
		}
		
		xp_dimensional = position[0] * data->length_bar;
		yp_dimensional = position[1] * data->length_bar;
		zp_dimensional = position[2] * data->length_bar;
		plastic_strain = 0.0;
		
		if (notch_type == 1 ) {
			if (xp_dimensional >= x_center - 0.5*notch_width && xp_dimensional <= x_center + 0.5*notch_width) {
				if(yp_dimensional>= ha_dimensional+hm_dimensional && yp_dimensional< ha_dimensional+hm_dimensional +notch_height) {
					plastic_strain = 1.0;
				}
			}
		}
		
		if (notch_type == 2 ) {
			if (xp_dimensional >= 1.15*x_center - 0.5*notch_width && xp_dimensional <= 1.15*x_center + 0.5*notch_width) {
				if(yp_dimensional>= ha_dimensional+hm_dimensional && yp_dimensional< ha_dimensional+hm_dimensional +notch_height) {
					if (zp_dimensional <= 0.8*z_center ) {
						plastic_strain = 1.0;
					}
				}
			}
			if (xp_dimensional >= 0.85*x_center - 0.5*notch_width && xp_dimensional <= 0.85*x_center + 0.5*notch_width) {
				if(yp_dimensional>= ha_dimensional+hm_dimensional && yp_dimensional< ha_dimensional+hm_dimensional +notch_height) {
					if (zp_dimensional >= 1.2*z_center) {
						plastic_strain = 1.0;
					}
				}
			}
		}
		
		/* user the setters provided for you */
		MPntPStokesPlSetField_plastic_strain(mpprop_pls,plastic_strain);
		MPntStdSetField_phase_index(material_point_std,         phase);
		MPntPStokesSetField_eta_effective(material_point_stokes,eta);
		MPntPStokesSetField_density(material_point_stokes,      rho);
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_pls);
	DataFieldRestoreAccess(PField_stokes);
	
	
	ierr = pTatinContextValid_Energy(c,&use_energy);CHKERRQ(ierr);
	if (use_energy) {
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
	}    
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialStokesVariableMarkers_Riftrh"
/* ASSIGN INITIAL VISCOSITY BASED ON INITIAL STRAIN RATE PRESSURE */
PetscErrorCode ModelApplyInitialStokesVariableMarkers_Riftrh(pTatinCtx user,Vec X,void *ctx)
{
	DM                           stokes_pack,dau,dap;
	PhysCompStokes               stokes;
	Vec                          Uloc,Ploc;
	PetscScalar                  *LA_Uloc,*LA_Ploc;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
	/* interpolate mesh velocity and pressure to markers and then sets viscosity */
	ierr = pTatin_EvaluateRheologyNonlinearities(user,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Riftrh_semi_eulerian"
/* DEFINE ALE */
PetscErrorCode ModelApplyUpdateMeshGeometry_Riftrh_semi_eulerian(pTatinCtx c,Vec X,void *ctx)
{
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
	
	/* ONLY VERTICAL REMESHING */
	ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Riftrh"
PetscErrorCode ModelOutput_Riftrh(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelRiftrhCtx  *data = (ModelRiftrhCtx*)ctx;
	PetscBool        active_energy;
	DataBucket       materialpoint_db;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	if (data->output_markers) { 
		ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	}
	
	if (data->output_markers) { 
		//  Write out just the stokes variable?
		//  const int nf = 1;
		//  const MaterialPointField mp_prop_list[] = { MPField_Stokes };
		//
		//  Write out just std, stokes and plastic variables
		const int nf = 3;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
		char mp_file_prefix[256];
		
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
		sprintf(mp_file_prefix,"%s_all_mp_data",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
	
	/* ENERGY OUTPUT */
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PhysCompEnergy energy;
		Vec            temperature;
		
		ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);
		
		ierr = pTatin3d_ModelOutput_Temperature_Energy(c,temperature,prefix);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Riftrh"
PetscErrorCode ModelDestroy_Riftrh(pTatinCtx c,void *ctx)
{
	ModelRiftrhCtx *data = (ModelRiftrhCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Riftrh"
/* ASSIGN ALL FUNCTIONS DEFINED ABOVE */
PetscErrorCode pTatinModelRegister_Riftrh(void)
{
	ModelRiftrhCtx *data;
	pTatinModel    m;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelRiftrhCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelRiftrhCtx));CHKERRQ(ierr);
	
	/* set initial values for model parameters */
	//	data->Lx = 0.0;
	//	data->param2 = 0;
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"Riftrh");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Riftrh);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Riftrh);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Riftrh);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Riftrh_semi_eulerian);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Riftrh);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Riftrh);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_Riftrh);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,   (void (*)(void))ModelApplyInitialStokesVariableMarkers_Riftrh);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Riftrh_semi_eulerian);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Riftrh);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Riftrh);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

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
 **    Filename:      rift_fastscape_3D.c
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

/*
	Input file for 3D rifting and surface processes using fastscape
*/


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

// added includes
#include "output_material_points.h"
#include "output_material_points_p0.h"
#include "pswarm.h"


#include "ptatin_models.h"

#include "rift_fastscape_3D_ctx.h"

#include "material_constants_energy.h"
#include "dmda_remesh.h"


PetscErrorCode ModelApplyUpdateMeshGeometry_rift_fastscape_3D_semi_eulerian(pTatinCtx c,Vec X,void *ctx);
PetscErrorCode ModelApplyMaterialBoundaryCondition_rift_fastscape_3D_semi_eulerian(pTatinCtx c,void *ctx);
PetscBool BCListEvaluator_rift_fastscape_3Dl( PetscScalar position[], PetscScalar *value, void *ctx );
PetscBool BCListEvaluator_rift_fastscape_3Dr( PetscScalar position[], PetscScalar *value, void *ctx );

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_rift_fastscape_3D"
PetscErrorCode ModelInitialize_rift_fastscape_3D(pTatinCtx c,void *ctx)
{
	Modelrift_fastscape_3DCtx  *data = (Modelrift_fastscape_3DCtx*)ctx;
	RheologyConstants       *rheology;
	EnergyMaterialConstants *matconstants_e;
	EnergyConductivityThreshold *matconstants_cond;
	EnergyConductivityConst *matconstants_cond_cst;
	EnergySourceDecay *matconstants_source_decay;
	EnergySourceAdiabaticAdvection *matconstants_source_adi_adv;
	PetscBool      flg;
	DataBucket     materialconstants;
	DataField      PField;
	PSwarm			pswarm;
	PetscBool      nondim,use_energy;
	//PetscScalar    dh_vx;
	PetscInt       regionidx;
	PetscReal      cm_per_yer2m_per_sec = 1.0e-2 / ( 365.0 * 24.0 * 60.0 * 60.0 ),phi1_rad,phi2_rad ;
	PetscReal      preexpA,Ascale,entalpy,Vmol,nexp,Tref;
	int 		   conductivity_type, density_type;
	double		   alpha, beta, rho_ref, Cp, k0, k1, T_threshold, dT, dTdy;

	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Temperature flag */
	use_energy = PETSC_TRUE;
	//use_energy=PETSC_FALSE;
	
	/* Particle Swarm */

	// PSwarm create
	ierr = PSwarmCreate(PETSC_COMM_WORLD,&pswarm);CHKERRQ(ierr);
	ierr = PSwarmSetOptionsPrefix(pswarm,"passive_");CHKERRQ(ierr);
	ierr = PSwarmSetPtatinCtx(pswarm,c);CHKERRQ(ierr);
	ierr = PSwarmSetTransportModeType(pswarm,PSWARM_TM_LAGRANGIAN);CHKERRQ(ierr);

	ierr = PSwarmSetFromOptions(pswarm);CHKERRQ(ierr);

	  /* Copy reference into model data for later use in different functions */
	  data->pswarm = pswarm;
	  
	/* set the default values of the material constant for this particular model */
	/*scaling */
	data->length_bar    = 100.0 * 1.0e3;
	data->viscosity_bar = 1.0e26;
	data->velocity_bar  = 1.0e-10;
	data->dimensional   = PETSC_TRUE;
	
	/* Box GEOMETRY */
	data->Lx = 600.0e3;
	data->Ly = 250.0e3;
	data->Lz = 600.0e3;
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_Lx",&data->Lx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_Ly",&data->Ly,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_Lz",&data->Lz,&flg);CHKERRQ(ierr);
	
	data->hc = 35.0e3;
	data->hm = 85.0e3;
	data->ha = 130.0e3;
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_hc",&data->hc,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_hm",&data->hm,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_ha",&data->ha,&flg);CHKERRQ(ierr);
	
	/* notch geometry */
	data->notch_type   = 1;
	data->notch_width  = 2.0e3;
	data->notch_height = 15.0e3;
	data->notch_base = data->Ly - data->notch_height;
	data->beta   = 60.; //Beta is the fault angle
	data->fault_width   = 5e3;
 	ierr = PetscOptionsGetInt(NULL,"-model_rift_fastscape_3D_notch_type",&data->notch_type,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_notch_width",&data->notch_width,&flg);CHKERRQ(ierr);
 	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_notch_height",&data->notch_height,&flg);CHKERRQ(ierr);
 	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_notch_base",&data->notch_base,&flg);CHKERRQ(ierr);
 	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_beta",&data->beta,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,"-model_rift_fastscape_3D_fault_width",&data->fault_width,&flg);CHKERRQ(ierr);

	/* velocity boundary condition geometry */
	//data->hvbx1 = 125.0e3;
	//data->hvbx2 = 115.0e3;
	data->vx_up = 0.5 * cm_per_yer2m_per_sec;
	//data->vybottom = 2.0 * data->vx_up * data->Ly / data->Lx;
	/* VELOCITY BOUNDARY CONDITION GEOMETRY */
	//ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_hvbx1",&data->hvbx1,&flg);CHKERRQ(ierr);
	//ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_hvbx2",&data->hvbx2,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_vx_up",&data->vx_up,&flg);CHKERRQ(ierr);
	
	/* material properties */
	data->rhoc = 2800.0;
	data->rhom = 3300.0;
	data->rhoa = 3250.0;
	data->etac = 1.0e28;
	data->etam = 1.0e21;
	data->etaa = 1.0e19;
	data->eps1 = 0.1; // initiation threshold of strain softening
	data->eps2 = 0.5; // threshold for full strain softening
	data->phi1 = 15.0; // friction angle before weakening
	data->phi2 = 2.0;
	data->coe1 = 2e7; // cohesion before weakening
	data->coe2 = 2e7;
	
	/* MATERIAL PARAMETERS */
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_rhoc",&data->rhoc,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_rhom",&data->rhom,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_rhoa",&data->rhoa,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_etac",&data->etac,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_etam",&data->etam,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_etaa",&data->etaa,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_eps1",&data->eps1,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_eps2",&data->eps2,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_phi1",&data->phi1,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_phi2",&data->phi2,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_coe1",&data->coe1,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_coe2",&data->coe2,&flg);CHKERRQ(ierr);
	
	/* Initial Temperature profile */
	if(use_energy) {
		data->thermalparams.nlayers  = 3;
		data->thermalparams.lscale   = data->length_bar;
		
		data->thermalparams.ytop[0]  = 250.0e3;
		data->thermalparams.thick[0] = 35.0e3;
		data->thermalparams.ttop[0]  = 0.0;
		data->thermalparams.cond[0]  = 2.25;
		data->thermalparams.hp[0]    = 0.9e-6;
		data->thermalparams.qbase[0] = 19.5e-3;
		
		data->thermalparams.ytop[1]  = 215.0e3;
		data->thermalparams.thick[1] = 85.0e3;
		data->thermalparams.ttop[1]  = 550.0;
		data->thermalparams.cond[1]  = 2.25;
		data->thermalparams.hp[1]    = 0.0;
		data->thermalparams.qbase[1] = 19.5e-3;
		
		data->thermalparams.ytop[2]  = 130.0e3;
		data->thermalparams.thick[2] = 130.0e3;
		data->thermalparams.ttop[2]  = 1330.0;
		data->thermalparams.cond[2]  = 48.75;
		data->thermalparams.hp[2]    = 0.0;
		data->thermalparams.qbase[2] = 19.5e-3;
		
		data->Ttop    = 0.0;
		data->Tbottom = 1382.0;
	}
	
	/* rheology parameters */
    ierr = pTatinGetRheology(c,&rheology);CHKERRQ(ierr);
	rheology->rheology_type  = RHEOLOGY_VP_STD;
	rheology->nphases_active = 6;
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global       = 1.0e+28;
	rheology->eta_lower_cutoff_global       = 1.0e+19;
	data->runmises = PETSC_FALSE;
	
	//    ierr = PetscOptionsGetInt(NULL,"-model_rift_fastscape_3D_param2",&data->param2,&flg);CHKERRQ(ierr);
	
	/* Material constants */
    ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
    ierr = MaterialConstantsSetDefaults(materialconstants);CHKERRQ(ierr);


    // ENERGY //
    // get fields entries for the various energy law methods

	/* Get the energy data fields, and the field entries */
	DataBucketGetDataFieldByName(materialconstants,EnergyMaterialConstants_classname,&PField);
	DataFieldGetEntries(PField,(void**)&matconstants_e);

	// Conductivity threshold //
	/* Get the conductivity threshold data fields, and the field entries */
	DataBucketGetDataFieldByName(materialconstants,EnergyConductivityThreshold_classname,&PField);
	DataFieldGetEntries(PField,(void**)&matconstants_cond);

	// Conductivity constant //
	/* Get the conductivity threshold data fields, and the field entries */
	DataBucketGetDataFieldByName(materialconstants,EnergyConductivityConst_classname,&PField);
	DataFieldGetEntries(PField,(void**)&matconstants_cond_cst);

	// Source Decay //
	/* Get the Energy source constant data fields, and the field entries */
	DataBucketGetDataFieldByName(materialconstants,EnergySourceDecay_classname,&PField);
	DataFieldGetEntries(PField,(void**)&matconstants_source_decay);

	// Source Adiabatic Advection //
	/* Get the Source data fields, and the field entries */
	DataBucketGetDataFieldByName(materialconstants,EnergySourceAdiabaticAdvection_classname,&PField);
	DataFieldGetEntries(PField,(void**)&matconstants_source_adi_adv);


	//-------------------------//
	/* PHASE 0, ASTHENOSPHERE */
    //-------------------------//
	regionidx = 0;
	alpha = 2e-5;
	beta = 0;
	rho_ref = data->rhoa;
	Cp = 1000;
	density_type = ENERGYDENSITY_CONSTANT;
	conductivity_type = ENERGYCONDUCTIVITY_TEMP_DEP_THRESHOLD;
	k0 = 2.25 ; //standard conductivity when T < T_threshold-dT
	k1 = 48.75 ; //conductivity for pseudo-adiabat, when T > T_threshold == Qm*dTdy
	T_threshold = 1380.0 ;
	dT = 50.0 ;
	dTdy = -0.4e-3;
	

	ierr = MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS_2,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);CHKERRQ(ierr);

	//VISCOSITY PARAMETERS - wet olivine karato
	preexpA  = 1.393e-14;
	Ascale   = 1.0;
	entalpy  = 430.e3;
	Vmol     = 15.0e-6;
	nexp     = 3.0;
	Tref     = 273.0;
	MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexpA,Ascale,entalpy,Vmol,nexp,Tref);
	
	//DENSITY PARAMETERS FOR STOKES
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,regionidx,data->rhom,2.0e-5,0.0);
	
	//PLASTICITY PARAMETERS
	phi1_rad = M_PI * data->phi1/180.0;
	phi2_rad = M_PI * data->phi2/180.0;
	MaterialConstantsSetValues_PlasticDP(materialconstants,regionidx,phi1_rad,phi2_rad,data->coe1,data->coe2,1.0e7,1.0e20);
	MaterialConstantsSetValues_PlasticMises(materialconstants,regionidx,1.0e8,1.0e8);
	MaterialConstantsSetValues_SoftLin(materialconstants,regionidx,data->eps1,data->eps2);
	
	// ENERGY //
	//Conductivity
	MaterialConstantsSetValues_ConductivityThreshold(regionidx,matconstants_cond, k0, k1, T_threshold, dT);
	//Source method: set all to NONE, then update the first entry of the array to ADIABATIC_ADVECTION
	EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[regionidx],ENERGYSOURCE_NONE);
	EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[regionidx],0,ENERGYSOURCE_ADIABATIC_ADVECTION);
	MaterialConstantsSetValues_SourceAdiabaticAdv(regionidx, matconstants_source_adi_adv, dTdy);

	MaterialConstantsSetValues_EnergyMaterialConstants(regionidx,matconstants_e,alpha,beta,rho_ref,Cp,density_type,conductivity_type,NULL);


	//------------------------------//
	/* PHASE 1, MANTLE LITHOSPHERE */
	//------------------------------//

	regionidx = 1;
	alpha = 2e-5;
	beta = 0;
	rho_ref = data->rhom;
	Cp = 1000;
	density_type = ENERGYDENSITY_CONSTANT;
	conductivity_type = ENERGYCONDUCTIVITY_CONSTANT;
	k0 = 2.25 ; //standard conductivity

	MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS_2,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);

	
	//VISCOSITY PARAMETERS - wet olivine x 5, karato
	preexpA = 1.393e-14;
	Ascale  = 5.0;
	entalpy = 430.0e3;
	Vmol    = 15.0e-6;
	nexp    = 3.0;
	Tref    = 273.0;
	MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexpA,Ascale,entalpy,Vmol,nexp,Tref);
	
	//DENSITY PARAMETERS FOR STOKES
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,regionidx,data->rhom,2.e-5,0.0);
	
	//PLASTICITY PARAMETERS
	phi1_rad = M_PI * data->phi1/180.0;
	phi2_rad = M_PI * data->phi2/180.0;
	MaterialConstantsSetValues_PlasticDP(materialconstants,regionidx,phi1_rad,phi2_rad,data->coe1,data->coe2,1.0e7,1.0e20);
	MaterialConstantsSetValues_PlasticMises(materialconstants,regionidx,1.0e8,1.0e8);
	MaterialConstantsSetValues_SoftLin(materialconstants,regionidx,data->eps1,data->eps2);
	

	//ENERGY//
	//MaterialConstantsSetValues_ConductivityConst(regionidx,matconstants_cond_cst,k0);
	EnergyConductivityConstSetField_k0(&matconstants_cond_cst[regionidx],k0);
	//Source method: set all energy source to NONE
	EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[regionidx],ENERGYSOURCE_NONE);
	MaterialConstantsSetValues_EnergyMaterialConstants(regionidx,matconstants_e,alpha,beta,rho_ref,Cp,density_type,conductivity_type,NULL);



	//-------------------------------------//
	/* PHASE 2, LOWER CRUST       */
	//-------------------------------------//

	regionidx = 2;
	alpha = 2e-5;
	beta = 0;
	rho_ref = data->rhoc;
	Cp = 1000;
	density_type = ENERGYDENSITY_CONSTANT;
	conductivity_type = ENERGYCONDUCTIVITY_CONSTANT;
	k0 = 2.25 ; //standard conductivity

	
	MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS_2,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);

	//VISCOSITY PARAMETERS
	preexpA = 8.5737e-28;
	Ascale  = 1.0;
	entalpy = 223.0e3;
	Vmol    = 0.0e-6;
	nexp    = 4.0;
	Tref    = 273.0;
	MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexpA,Ascale,entalpy,Vmol,nexp,Tref);
	
	//DENSITY PARAMETERS STOKES
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,regionidx,data->rhoc,2.0e-5,0.0);
	
	//PLASTICITY PARAMETERS
	phi1_rad = M_PI * data->phi1/180.0;
	phi2_rad = M_PI * data->phi2/180.0;
	MaterialConstantsSetValues_PlasticDP(materialconstants,regionidx,phi1_rad,phi2_rad,data->coe1,data->coe2,1.0e7,1.0e20);
	MaterialConstantsSetValues_PlasticMises(materialconstants,regionidx,1.0e8,1.0e8);
	MaterialConstantsSetValues_SoftLin(materialconstants,regionidx,data->eps1,data->eps2);



  /* ENERGY FLOW LAW DEFS */
  {
    DataField      PField;
    EnergyMaterialConstants        *matconstants_e;
    EnergyConductivityConst        *matconstants_k_const;
    EnergySourceConst              *matconstants_h_const;
    EnergyConductivityThreshold	   *matconstants_cond;
    EnergySourceAdiabaticAdvection *matconstants_source_adi_adv;
    double alpha,beta,rho_ref,Cp,k,k0,k1,dT,T_threshold,dTdy;    
    int 	 conductivity_type, density_type, conductivity_type_asthen;
    
    
    /* Get the energy data fields, and the field entries */
    DataBucketGetDataFieldByName(materialconstants,EnergyMaterialConstants_classname,&PField);
    DataFieldGetEntries(PField,(void**)&matconstants_e);

    DataBucketGetDataFieldByName(materialconstants,EnergyConductivityConst_classname,&PField);
    DataFieldGetEntries(PField,(void**)&matconstants_k_const);

    DataBucketGetDataFieldByName(materialconstants,EnergyConductivityThreshold_classname,&PField);
    DataFieldGetEntries(PField,(void**)&matconstants_cond);

    DataBucketGetDataFieldByName(materialconstants,EnergySourceConst_classname,&PField);
    DataFieldGetEntries(PField,(void**)&matconstants_h_const);

    DataBucketGetDataFieldByName(materialconstants,EnergySourceAdiabaticAdvection_classname,&PField);
    DataFieldGetEntries(PField,(void**)&matconstants_source_adi_adv);    

    conductivity_type = ENERGYCONDUCTIVITY_CONSTANT;
    density_type      = ENERGYDENSITY_CONSTANT;
    conductivity_type_asthen = ENERGYCONDUCTIVITY_TEMP_DEP_THRESHOLD;

    //ASTHENOSPHERE
    alpha   = 2.0e-5;
    beta    = 0.0;
    rho_ref = data->rhoa;
    Cp      = 1000.0;
    k0       = 1.0e-6 * rho_ref * Cp;
    k1       = 48.75; //artificial conductivity to conserve the heat flow
    T_threshold = 1350.0;
    dT       = 50.0;
    dTdy     = 0.4e-3;
    MaterialConstantsSetValues_EnergyMaterialConstants(0,matconstants_e,alpha,beta,rho_ref,Cp,density_type,conductivity_type_asthen,NULL);
    //EnergyConductivityConstSetField_k0(&matconstants_k_const[0],k);
    MaterialConstantsSetValues_ConductivityThreshold(0,matconstants_cond, k0, k1, T_threshold, dT);
    EnergySourceConstSetField_HeatSource(&matconstants_h_const[0],0.0);
    MaterialConstantsSetValues_SourceAdiabaticAdv(0, matconstants_source_adi_adv, dTdy);    


    //LITHOSPHERIC MANTLE
    alpha   = 2.0e-5;
    beta    = 0.0;
    rho_ref = data->rhom;
    Cp      = 1000.0;
    k       = 1.0e-6 * rho_ref * Cp;
    MaterialConstantsSetValues_EnergyMaterialConstants(1,matconstants_e,alpha,beta,rho_ref,Cp,density_type,conductivity_type,NULL);
    EnergyConductivityConstSetField_k0(&matconstants_k_const[1],k);
    EnergySourceConstSetField_HeatSource(&matconstants_h_const[1],0.0);
    
    //LOWER CRUST
    alpha   = 2.0e-5;
    beta    = 0.0;
    rho_ref = data->rhoc;
    Cp      = 1000.0;
    k       = 1.0e-6 * rho_ref * Cp;
    MaterialConstantsSetValues_EnergyMaterialConstants(2,matconstants_e,alpha,beta,rho_ref,Cp,density_type,conductivity_type,NULL);
    EnergyConductivityConstSetField_k0(&matconstants_k_const[2],k);
    EnergySourceConstSetField_HeatSource(&matconstants_h_const[2],0.9e-6);

    //LOWER CRUST MARKERS
    MaterialConstantsSetValues_EnergyMaterialConstants(3,matconstants_e,alpha,beta,rho_ref,Cp,density_type,conductivity_type,NULL);
    EnergyConductivityConstSetField_k0(&matconstants_k_const[3],k);
    EnergySourceConstSetField_HeatSource(&matconstants_h_const[3],0.9e-6);

    //UPPER CRUST
    MaterialConstantsSetValues_EnergyMaterialConstants(4,matconstants_e,alpha,beta,rho_ref,Cp,density_type,conductivity_type,NULL);
    EnergyConductivityConstSetField_k0(&matconstants_k_const[4],k);
    EnergySourceConstSetField_HeatSource(&matconstants_h_const[4],0.9e-6);

    //UPPER CRUST MARKERS
    MaterialConstantsSetValues_EnergyMaterialConstants(5,matconstants_e,alpha,beta,rho_ref,Cp,density_type,conductivity_type,NULL);
    EnergyConductivityConstSetField_k0(&matconstants_k_const[5],k);
    EnergySourceConstSetField_HeatSource(&matconstants_h_const[5],0.9e-6);


    //phase = 2,3;
    //kappa = 1.0e-6/data->length_bar/data->length_bar*data->time_bar;
    //H     = 0.9e-6/data->pressure_bar*data->time_bar;
    
    //phase = 0, 1;
    //kappa = 1.0e-6/data->length_bar/data->length_bar*data->time_bar;
    //H     = 0.0;
    
    EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[0],ENERGYSOURCE_NONE);
    EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[1],ENERGYSOURCE_NONE);
    EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[2],ENERGYSOURCE_NONE);
    EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[3],ENERGYSOURCE_NONE);
    EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[4],ENERGYSOURCE_NONE);
    EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[5],ENERGYSOURCE_NONE);
   
    EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[0],0,ENERGYSOURCE_CONSTANT);
    EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[1],0,ENERGYSOURCE_CONSTANT);
    EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[2],0,ENERGYSOURCE_CONSTANT);
    EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[3],0,ENERGYSOURCE_CONSTANT);
    EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[4],0,ENERGYSOURCE_CONSTANT);
    EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[5],0,ENERGYSOURCE_CONSTANT);
    //pseudo adiabat in the asthenosphere
    EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[0],1,ENERGYSOURCE_ADIABATIC_ADVECTION);
    
  }
  
  
  
	/* Read the options */
	/* cutoff */
	ierr = PetscOptionsGetBool(NULL,"-model_rift_fastscape_3D_apply_viscosity_cutoff_global",&rheology->apply_viscosity_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_eta_lower_cutoff_global",&rheology->eta_lower_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_eta_upper_cutoff_global",&rheology->eta_upper_cutoff_global,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-model_rift_fastscape_3D_runwithmises",&data->runmises,NULL);CHKERRQ(ierr);
	/* scaling */
	nondim = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-model_rift_fastscape_3D_nondimensional",&nondim,NULL);CHKERRQ(ierr);
	if (nondim) {
		data->dimensional = PETSC_FALSE;
	} else {
		ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_vis_bar",&data->viscosity_bar,NULL);CHKERRQ(ierr);
		ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_vel_bar",&data->velocity_bar,NULL);CHKERRQ(ierr);
		ierr = PetscOptionsGetReal(NULL,"-model_rift_fastscape_3D_length_bar",&data->length_bar,NULL);CHKERRQ(ierr);
	}
	
	/* compute vxdown */
	//dh_vx = data->hvbx1 - data->hvbx2;
	//data->vx_down = -data->vx_up * (data->hc + data->hm +dh_vx/2.0) / (data->ha - dh_vx/2.0);
	
	/* reports before scaling */
	PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  input: -model_rift_fastscape_3D_Lx : %+1.4e [SI]\n", data->Lx );
	PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  input: -model_rift_fastscape_3D_Ly : %+1.4e [SI]\n", data->Ly );
	PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  input: -model_rift_fastscape_3D_Lz : %+1.4e [SI]\n", data->Lz );
	PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  -model_rift_fastscape_3D_vx_up [m/s]:  %+1.4e \n", data->vx_up);
	PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  -model_rift_fastscape_3D_rhoc [kg/m^3] :%+1.4e \n", data->rhoc );
	PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  -model_rift_fastscape_3D_rhom [kg/m^3] :%+1.4e \n", data->rhom );
	PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  -model_rift_fastscape_3D_rhoa [kg/m^3] :%+1.4e \n", data->rhoa );
	
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		MaterialConstantsPrintAll(materialconstants,regionidx);
	}
	
	if (data->dimensional) {
		/* Compute additional scaling parameters */
		data->time_bar      = data->length_bar / data->velocity_bar;
		data->pressure_bar  = data->viscosity_bar/data->time_bar;
		data->density_bar   = data->pressure_bar / data->length_bar;
		
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]:  during the solve scaling will be done using \n");
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  L*    : %1.4e [m]\n", data->length_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  U*    : %1.4e [m.s^-1]\n", data->velocity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  t*    : %1.4e [s]\n", data->time_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  eta*  : %1.4e [Pa.s]\n", data->viscosity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  rho*  : %1.4e [kg.m^-3]\n", data->density_bar );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D]  P*    : %1.4e [Pa]\n", data->pressure_bar );
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
		//scale velocity
		data->vx_up    = data->vx_up  /data->velocity_bar;
		data->vx_down  = data->vx_down/data->velocity_bar;
		//data->vybottom = data->vybottom/data->velocity_bar;
		data->hvbx1    = data->hvbx1/data->length_bar;
		data->hvbx2    = data->hvbx2/data->length_bar;
		//scale rho
		data->rhoc = data->rhoc/data->density_bar;
		data->rhom = data->rhom/data->density_bar;
		data->rhoa = data->rhoa/data->density_bar;
		
    
  //  printf("[crust phase 2] kappa %1.9e H %1.9e\n",
   //        1.0e-6/data->length_bar/data->length_bar*data->time_bar,
   //        0.9e-6/data->pressure_bar*data->time_bar);
   // 0.9e-6*data->time_bar/data->pressure_bar);
    
		// scale material properties
		for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
			MaterialConstantsScaleAll(materialconstants,regionidx,data->length_bar,data->velocity_bar,data->time_bar,data->viscosity_bar,data->density_bar,data->pressure_bar);
      MaterialConstantsEnergyScaleAll(materialconstants,regionidx,data->length_bar,data->time_bar,data->pressure_bar);
		}
		
		/* Reports scaled values */		
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] scaled value    -model_rift_fastscape_3D_Lx   :  %+1.4e \n", data->Lx );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] scaled value    -model_rift_fastscape_3D_Ly   :  %+1.4e \n", data->Ly );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] scaled value    -model_rift_fastscape_3D_Lz   :  %+1.4e \n", data->Lz );
		
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] scaled value   -model_rift_fastscape_3D_vx_up:%+1.4e    -model_rift_fastscape_3D_vx_down:%+1.4e \n", data->vx_up ,data->vx_down);
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] scaled value   -model_rift_fastscape_3D_rhoc:%+1.4e \n", data->rhoc );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] scaled value   -model_rift_fastscape_3D_rhom:%+1.4e \n", data->rhom );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] scaled value   -model_rift_fastscape_3D_rhoa:%+1.4e \n", data->rhoa );
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] scaled value for material parameters\n");
		for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
			MaterialConstantsPrintAll(materialconstants,regionidx);
		}
	}
	
	data->output_markers = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-model_rift_fastscape_3D_output_markers",&data->output_markers,NULL);CHKERRQ(ierr);
	
	/* USE ENERGY EQUATION */
	if (use_energy) {
		ierr = PetscOptionsInsertString("-activate_energy");CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

/* SET AN INITIAL BACK GROUND STRAIN RATE, TEMPERATURE, PRESSURE */
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_rift_fastscape_3D"
PetscErrorCode ModelApplyInitialSolution_rift_fastscape_3D(pTatinCtx c,Vec X,void *ctx)
{
	Modelrift_fastscape_3DCtx *data = (Modelrift_fastscape_3DCtx*)ctx;
	PetscBool       use_energy;
	PhysCompStokes  stokes;
	DM              stokes_pack,dau,dap;
	DMDAVecTraverse3d_InterpCtx IntpCtx;
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	Vec            velocity,pressure;
	PetscReal      MeshMax[3],MeshMin[3],height,length,vxl,vxr,vybottom;
	PetscBool      use_initial_up_field = PETSC_FALSE;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscOptionsGetBool(NULL,"-model_rift_fastscape_3D_use_initial_up_field",&use_initial_up_field,NULL);
	
	if (use_initial_up_field) {
		PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] Using velocity from boundary condition and mantle hydrostatic pressure as initial condition\n");
		
		ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
                ierr = PhysCompStokesGetDMComposite(stokes,&stokes_pack);CHKERRQ(ierr);
		ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
		ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		
		/* velocity intial condition - background strain */
		ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
		vxl = -data->vx_up;
		vxr =  data->vx_up;
		ierr = DMDAGetBoundingBox(dau,MeshMin,MeshMax);CHKERRQ(ierr);
		height = MeshMax[1] - MeshMin[1];
		length = MeshMax[0] - MeshMin[0];
		vybottom = 2.0*data->vx_up * height / length;
		
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

	/* set swarm velocity and pressure vector once*/
	ierr = PSwarmAttachStateVecVelocityPressure(data->pswarm,X);CHKERRQ(ierr);


	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Modelrift_fastscape_3D_DefineBCList"
PetscErrorCode Modelrift_fastscape_3D_DefineBCList(BCList bclist,DM dav,pTatinCtx user,Modelrift_fastscape_3DCtx *data)
{
	PhysCompStokes stokes;
	DM             stokes_pack,dau,dap;
	PetscScalar    vxl,vxr,vybottom,zero,height,length;
	PetscReal      MeshMin[3],MeshMax[3];
	PetscInt       vbc_type;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	ierr = PhysCompStokesGetDMComposite(stokes,&stokes_pack);CHKERRQ(ierr);
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	
	zero=0.0;
	//    vbc_type = 1; /* in / out flow condition on the sides */
	vbc_type = 2; /* outflow condition on the sides, inflow condition on the base */
	
	/*if (vbc_type == 1) {
		
		// infilling free slip base
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		// free surface top
		
		//extension along face of normal x
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_rift_fastscape_3Dl,(void*)data);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_rift_fastscape_3Dr,(void*)data);CHKERRQ(ierr);
		
		//free slip base
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		// no flow in z
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}*/
	
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

/*#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_rift_fastscape_3Dl"
PetscBool BCListEvaluator_rift_fastscape_3Dl( PetscScalar position[], PetscScalar *value, void *data )
{
	Modelrift_fastscape_3DCtx *datal = (Modelrift_fastscape_3DCtx*)data;
	PetscBool   impose_dirichlet = PETSC_TRUE;
	PetscScalar vx,dh_vx;
	
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
#define __FUNCT__ "BCListEvaluator_rift_fastscape_3Dr"
PetscBool BCListEvaluator_rift_fastscape_3Dr( PetscScalar position[], PetscScalar *value, void *data )
{
	Modelrift_fastscape_3DCtx *datal = (Modelrift_fastscape_3DCtx*)data;
	PetscBool   impose_dirichlet = PETSC_TRUE;
	PetscScalar vx,dh_vx;
	
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
*/

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_rift_fastscape_3D"
PetscErrorCode ModelApplyBoundaryCondition_rift_fastscape_3D(pTatinCtx c,void *ctx)
{
	Modelrift_fastscape_3DCtx *data = (Modelrift_fastscape_3DCtx*)ctx;
	PhysCompStokes stokes;
	DM             stokes_pack,dav,dap;
    BCList         u_bclist;
	PetscBool      use_energy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	ierr = PhysCompStokesGetDMComposite(stokes,&stokes_pack);CHKERRQ(ierr);
    ierr = PhysCompStokesGetBCList(stokes,&u_bclist,NULL);CHKERRQ(ierr);
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

	ierr = Modelrift_fastscape_3D_DefineBCList(u_bclist,dav,c,data);CHKERRQ(ierr);
	
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
#define __FUNCT__ "ModelApplyBoundaryConditionMG_rift_fastscape_3D"
PetscErrorCode ModelApplyBoundaryConditionMG_rift_fastscape_3D(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	Modelrift_fastscape_3DCtx *data = (Modelrift_fastscape_3DCtx*)ctx;
	PetscInt       n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = Modelrift_fastscape_3D_DefineBCList(bclist[n],dav[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_rift_fastscape_3D"
PetscErrorCode ModelApplyInitialMeshGeometry_rift_fastscape_3D(pTatinCtx c,void *ctx)
{
	Modelrift_fastscape_3DCtx *data = (Modelrift_fastscape_3DCtx*)ctx;
	PhysCompStokes         stokes;
	DM                     stokes_pack,dau,dap;
	PetscErrorCode         ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	ierr = PhysCompStokesGetDMComposite(stokes,&stokes_pack);CHKERRQ(ierr);
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
    
	ierr = DMDASetUniformCoordinates(dau,0.0,data->Lx,0.0,data->Ly,0.0,data->Lz);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"[rift_fastscape_3D] Lx = %1.4e \n", data->Lx );

        // {
	//    PetscReal gvec[] = { 0.0, -10.0, 0.0 };
	//    ierr = PhysCompStokesSetGravityVector(c->stokes_ctx,gvec);CHKERRQ(ierr);
	//  }



	 //remesh vertically and preserve topography
	{
		PetscInt npoints,dir;
		PetscReal xref[10],xnat[10];

		npoints = 3;
		xref[0] = 0.0;
		xref[1] = 0.3;
		xref[2] = 1.0;

		xnat[0] = 0.0;
		xnat[1] = 0.6;
		xnat[2] = 1.0;


		dir = 1;
		ierr = DMDACoordinateRefinementTransferFunction(dau,dir,PETSC_TRUE,npoints,xref,xnat);CHKERRQ(ierr);
		ierr = DMDABilinearizeQ2Elements(dau);CHKERRQ(ierr);
	}


	
	{
	    PetscReal gvec[] = { 0.0, -10.0, 0.0 };
	    ierr = PhysCompStokesSetGravityVector(c->stokes_ctx,gvec);CHKERRQ(ierr);
	  }

	ierr = PSwarmSetUp(data->pswarm);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_rift_fastscape_3D"
PetscErrorCode ModelApplyInitialMaterialGeometry_rift_fastscape_3D(pTatinCtx c,void *ctx)
{
	Modelrift_fastscape_3DCtx *data = (Modelrift_fastscape_3DCtx*)ctx;
	PetscInt               p,n_mp_points;
	PetscInt               notch_type;
	DataBucket             db;
	DataField              PField_std,PField_stokes,PField_pls;
	PetscScalar            ha_dimensional,hm_dimensional,notch_height,notch_width,notch_base,x_center,z_center;
	PetscScalar            xp_dimensional,yp_dimensional,zp_dimensional;
	PetscErrorCode 		   ierr;
	PetscBool              use_energy;
	
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* define properties on material points */
	ierr = pTatinGetMaterialPoints(c,&db,NULL);CHKERRQ(ierr);

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
	x_center       = 0.5 * data->Lx * data->length_bar;
	z_center       = 0.5 * data->Lz * data->length_bar;
	notch_width    = data->notch_width;//20.0e3;
	notch_height   = data->notch_height;//10.0e3;
	notch_base     = data->notch_base;
	
	/* marker loop */
	DataBucketGetSizes(db,&n_mp_points,0,0);

	notch_type = data->notch_type;

	
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
		if ((position[1] > data->ha) && (position[1] < (data->ha+data->hm))) {
			phase = 1;
			eta=data->etam;
			rho=data->rhom;
		}
		if (position[1] > (data->ha+data->hm)) {
			phase = 2;
			eta=data->etac;
			rho=data->rhoc;
		}
/*		
		//make 25 km (0.25) stripes in the lower crust  in the x and z direction
		if (position[1] > data->ha + data->hm) {
			int j=0;
			while(j<100){
				if ((position[0] > (0.25*j + 0.25)) && (position[0] < (0.25*j + 0.4))){
					phase = 3;
					eta=data->etac;
					rho=data->rhoc;
				}
				if ((position[2] > (0.25*j + 0.25)) && (position[2] < (0.25*j + 0.4))){
					phase = 3;
					eta=data->etac;
					rho=data->rhoc;
				}
				j++;
				j++;
			}
		}

		//make 2 flat layers  in the upper crust
		if (position[1] > data->ha + data->hm + 0.2) {
			phase = 4;
			eta=data->etac;
			rho=data->rhoc;
			if (position[1] > data->ha + data->hm + 0.275){
				phase = 5;
				eta=data->etac;
				rho=data->rhoc;
			}
		}
*/


		xp_dimensional = position[0] * data->length_bar;
		yp_dimensional = position[1] * data->length_bar;
		zp_dimensional = position[2] * data->length_bar;
		plastic_strain = 0.0;


		if (notch_type == 1 ) {
			
			double beta,xref,fault_width;
			double convert=M_PI/180.;

			beta = data->beta; //rotation angle in x-z plane, in degree
			beta =  convert * beta; // convert deg to rad
			fault_width = data->fault_width; //Fault width in km
			
			if (yp_dimensional >= notch_base) {
				if (zp_dimensional >= 100e3 && zp_dimensional < 200e3){
					xref = 200e3 + (yp_dimensional-notch_base)/2.;
					if ( fabs(xp_dimensional - xref) <= fault_width) {
						plastic_strain = 1.0;
						phase = 3;
					}
					xref = 400e3 + (yp_dimensional-notch_base)/2.;
					if ( fabs(xp_dimensional - xref) <= fault_width) {
						plastic_strain = 1.0;
						phase = 3;
					}				

				}
				if (zp_dimensional >= 250e3 && zp_dimensional < 350e3){
					xref = 300e3 + (yp_dimensional-notch_base)/2.;
					if ( fabs(xp_dimensional - xref) <= 5e3) {
						plastic_strain = 1.0;
						phase = 3;
					}	
				}
				if (zp_dimensional >= 400e3 && zp_dimensional < 500e3){
					xref = 200e3 + (yp_dimensional-notch_base)/2.;
					if ( fabs(xp_dimensional - xref) <= fault_width) {
						plastic_strain = 1.0;
						phase = 3;
					}
					xref = 400e3 + (yp_dimensional-notch_base)/2.;
					if ( fabs(xp_dimensional - xref) <= fault_width) {
						plastic_strain = 1.0;
						phase = 3;
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
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialStokesVariableMarkers_rift_fastscape_3D"
/* ASSIGN INITIAL VISCOSITY BASED ON INITIAL STRAIN RATE PRESSURE */
PetscErrorCode ModelApplyInitialStokesVariableMarkers_rift_fastscape_3D(pTatinCtx user,Vec X,void *ctx)
{
	DM                           stokes_pack,dau,dap;
	PhysCompStokes               stokes;
	Vec                          Uloc,Ploc;
	PetscScalar                  *LA_Uloc,*LA_Ploc;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	ierr = PhysCompStokesGetDMComposite(stokes,&stokes_pack);CHKERRQ(ierr);
	
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

// adding particles on the lower boundary to accommodate inflow
// adding particles on the left and right boundary to accommodate inflow
#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_rift_fastscape_3D_semi_eulerian"
PetscErrorCode ModelApplyMaterialBoundaryCondition_rift_fastscape_3D_semi_eulerian(pTatinCtx c,void *ctx)
{
#if 0
	Modelrift_fastscape_3DCtx     *data = (Modelrift_fastscape_3DCtx*)ctx;
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
	ierr = PhysCompStokesGetDMComposite(stokes,&stokes_pack);CHKERRQ(ierr);
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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_rift_fastscape_3D_semi_eulerian"
/* DEFINE ALE */
PetscErrorCode ModelApplyUpdateMeshGeometry_rift_fastscape_3D_semi_eulerian(pTatinCtx c,Vec X,void *ctx)
{
	Modelrift_fastscape_3DCtx  *data = (Modelrift_fastscape_3DCtx*)ctx;
	PetscReal        step;
	//PetscReal        MeshMin[3],MeshMax[3],avg[3],height,length;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	PetscErrorCode   ierr;

	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	//if (c->step%(50) == 0) {	
	/* fully lagrangian update */
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	ierr = PhysCompStokesGetDMComposite(stokes,&stokes_pack);CHKERRQ(ierr);
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	


	/* ONLY VERTICAL REMESHING */
	//ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);
	
	/* SURFACE REMESHING */
	ierr = UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX(dav,velocity,NULL,step);
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
        /* UPDATE mesh refinement scheme */
	//remesh vertically and preserve topography
	{
		PetscInt npoints,dir;
		PetscReal xref[10],xnat[10];

		npoints = 3;
		xref[0] = 0.0;
		xref[1] = 0.3;
		xref[2] = 1.0;

		xnat[0] = 0.0;
		xnat[1] = 0.6;
		xnat[2] = 1.0;


		dir = 1;
		ierr = DMDACoordinateRefinementTransferFunction(dav,dir,PETSC_TRUE,npoints,xref,xnat);CHKERRQ(ierr);
		ierr = DMDABilinearizeQ2Elements(dav);CHKERRQ(ierr);
	}	
	
	/* PSwarm update */
	ierr = PSwarmFieldUpdateAll(data->pswarm);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelOutput_rift_fastscape_3D"
PetscErrorCode ModelOutput_rift_fastscape_3D(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	Modelrift_fastscape_3DCtx  *data = (Modelrift_fastscape_3DCtx*)ctx;
	PetscBool        active_energy;
	DataBucket       materialpoint_db;
	PetscInt         Step,outFre;
	PetscErrorCode   ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

	if (data->output_markers) { 
		ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	}
	
	if (data->output_markers) { 
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
		//  Write out just the stokes variable?
		//  const int nf = 1;
		//  const MaterialPointField mp_prop_list[] = { MPField_Stokes };
		//
		//  Write out just std, stokes and plastic variables
		const int nf = 3;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
		char mp_file_prefix[256];
		
		sprintf(mp_file_prefix,"%s_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
		//sprintf(mp_file_prefix,"%s_all_mp_data",prefix);
		//ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}

	/* MPOINTS_CELL OUTPUT */
	Step = c->step;
	outFre = c->output_frequency;
	if (Step%(outFre*2) == 0) {
		/*
			const int  nf = 4;
			const MaterialPointVariable mp_prop_list[] = { MPV_region, MPV_viscosity, MPV_density, MPV_plastic_strain };
			ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
		*/
		
		MaterialPointVariable vars[] = { MPV_region, MPV_viscosity, MPV_density, MPV_plastic_strain, MPV_diffusivity, MPV_heat_source };
		ierr = pTatin3dModelOutput_MarkerCellFieldsP0_PetscVec(c,PETSC_FALSE,sizeof(vars)/sizeof(MaterialPointVariable),vars,prefix);CHKERRQ(ierr);
		
		{
			const int                   nf = 4;
			const MaterialPointVariable mp_prop_list[] = { MPV_region, MPV_viscosity, MPV_density, MPV_plastic_strain };
			
			ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
		}

		/* Dump the PSwarm */

		ierr = PSwarmView(data->pswarm,PSW_VT_SINGLETON);CHKERRQ(ierr);
		ierr = PSwarmViewInfo(data->pswarm);CHKERRQ(ierr);
		
		ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
		if (active_energy) {
			PhysCompEnergy energy;
			Vec            temperature;
			
			ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
			ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);
			
			ierr = pTatin3d_ModelOutput_Temperature_Energy(c,temperature,prefix);CHKERRQ(ierr);
		}


		ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
		ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
		
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_rift_fastscape_3D"
PetscErrorCode ModelDestroy_rift_fastscape_3D(pTatinCtx c,void *ctx)
{
	Modelrift_fastscape_3DCtx *data = (Modelrift_fastscape_3DCtx*)ctx;
	PetscErrorCode         ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	/* destroy PSwarm*/
	ierr = PSwarmDestroy(&data->pswarm);CHKERRQ(ierr);
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_rift_fastscape_3D"
/* ASSIGN ALL FUNCTIONS DEFINED ABOVE */
PetscErrorCode pTatinModelRegister_rift_fastscape_3D(void)
{
	Modelrift_fastscape_3DCtx *data;
	pTatinModel            m;
	PetscErrorCode         ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(Modelrift_fastscape_3DCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(Modelrift_fastscape_3DCtx));CHKERRQ(ierr);
	
	/* set initial values for model parameters */
	//	data->Lx = 0.0;
	//	data->param2 = 0;
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"rift_fastscape_3D");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_rift_fastscape_3D);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_rift_fastscape_3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_rift_fastscape_3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_rift_fastscape_3D_semi_eulerian);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_rift_fastscape_3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_rift_fastscape_3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_rift_fastscape_3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,   (void (*)(void))ModelApplyInitialStokesVariableMarkers_rift_fastscape_3D);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_rift_fastscape_3D_semi_eulerian);CHKERRQ(ierr);
	
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_rift_fastscape_3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_rift_fastscape_3D);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

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
 **    filename:   model_rift3D_I_ops.c
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


#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_utils.h"
#include "dmda_bcs.h"
#include "data_bucket.h"
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
#include "geometry_object.h"
#include <material_constants_energy.h>
#include "dmda_remesh.h"
#include "pswarm.h"

#include "rift3D_I_ctx.h"

#define REMOVE_FACE_INJECTION

PetscErrorCode GeometryObjectSetFromOptions_Box(GeometryObject go);
PetscErrorCode GeometryObjectSetFromOptions_InfLayer(GeometryObject go);
PetscErrorCode GeometryObjectSetFromOptions_EllipticCylinder(GeometryObject go);

PetscErrorCode ModelApplyUpdateMeshGeometry_Rift3D_I_semi_eulerian(pTatinCtx c,Vec X,void *ctx);
PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D_I_semi_eulerian(pTatinCtx c,void *ctx);

PetscErrorCode ModelApplyInitialMaterialGeometry_Notchtest_Rift3D_I(pTatinCtx c,void *ctx);

PetscErrorCode ModelInitialize_Rift3D_I(pTatinCtx c,void *ctx)
{
  ModelRift3D_ICtx *data = (ModelRift3D_ICtx*)ctx;
  RheologyConstants       *rheology;
  DataBucket              materialconstants;
  PetscBool               nondim;
  PetscBool               use_passive_tracers = PETSC_FALSE;
  PetscScalar             vx,vy,vz,vx_I,vz_I,Sx,Sy,Sz;
  PetscInt                regionidx,region_idx;
  PetscInt                p,n_mp_points;
  PetscReal               cm_per_yer2m_per_sec = 1.0e-2 / ( 365.0 * 24.0 * 60.0 * 60.0 ) ;
  PetscReal               rho_ref,Cp;
  PetscReal               preexpA_uc,Ascale_uc,entalpy_uc,Vmol_uc,nexp_uc,Tref;
  PetscReal               preexpA_lc,Ascale_lc,entalpy_lc,Vmol_lc,nexp_lc;
  PetscReal               preexpA_ml,Ascale_ml,entalpy_ml,Vmol_ml,nexp_ml;
  PetscReal               preexpA_ma,Ascale_ma,entalpy_ma,Vmol_ma,nexp_ma;
  PetscReal               alpha,beta,rho_uc,rho_lc,rho_ml,rho_ma;
  PetscReal               phi,phi_inf,phi_rad,phi_inf_rad,Co,Co_inf,Hst_cutoff,Tens_cutoff,eps_min,eps_max;
  DataField               PField;
  EnergySourceConst       *data_Q;
  EnergyConductivityConst *data_k;
  DataField               PField_Q,PField_std,PField_energy,PField_k;
  PSwarm                  pswarm;
  DataBucket              db,materialpoint_db;
  EnergyMaterialConstants *matconstants_e;
  int                     source_type[7] = {0, 0, 0, 0, 0, 0, 0};
  PetscErrorCode          ierr;

  PetscFunctionBegin;

  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  PetscPrintf(PETSC_COMM_WORLD,"Rift model expects the following dimensions for input\n");
  PetscPrintf(PETSC_COMM_WORLD," Box geometry: [m] \n");
  PetscPrintf(PETSC_COMM_WORLD," Viscosity:    [Pa.s] \n");
  PetscPrintf(PETSC_COMM_WORLD," Velocity:     [m/sec] \n");
  PetscPrintf(PETSC_COMM_WORLD," Density:      [kg/m^3] \n");

  PetscPrintf(PETSC_COMM_WORLD,"if you wish to use non dimensional input you must add -model_rift3D_I_dimensional \n");
  ierr = pTatinGetRheology(c,&rheology);CHKERRQ(ierr);

  /* Particle Swarm */
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_use_passive_tracers",&use_passive_tracers,NULL);CHKERRQ(ierr);

  if (use_passive_tracers){ 
  // PSwarm create
  ierr = PSwarmCreate(PETSC_COMM_WORLD,&pswarm);CHKERRQ(ierr);
  ierr = PSwarmSetOptionsPrefix(pswarm,"passive_");CHKERRQ(ierr);
  ierr = PSwarmSetPtatinCtx(pswarm,c);CHKERRQ(ierr);
  ierr = PSwarmSetTransportModeType(pswarm,PSWARM_TM_LAGRANGIAN);CHKERRQ(ierr);
 
  ierr = PSwarmSetFromOptions(pswarm);CHKERRQ(ierr);

   /* Copy reference into model data for later use in different functions */
   data->pswarm = pswarm;
   }

  rheology->rheology_type = RHEOLOGY_VP_STD;
  /* force energy equation to be introduced */
  ierr = PetscOptionsInsertString(NULL,"-activate_energy true");CHKERRQ(ierr);
  /* I REALLY DONT LIKE THE FOLLOWING ONE, SHOULD BE  in model data */
  rheology->nphases_active = 4;
  rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
  rheology->eta_upper_cutoff_global = 1.e+25;
  rheology->eta_lower_cutoff_global = 1.e+19;
  data->runmises = PETSC_FALSE;
  /* set the deffault values of the material constant for this particular model */
  /*scaling */
  data->length_bar    = 100.0 * 1.0e3;
  data->viscosity_bar = 1e22;
  data->velocity_bar  = 1.0e-10;
  data->dimensional   = PETSC_TRUE;

  /* box geometry, m */
  data->Lx =  12.0e5;
  data->Ly =  0.0e5;
  data->Lz =  6.0e5;
  //data->Ox =  -6.0e5;
  data->Ox =  0.0e5;
  data->Oy =  -2.5e5;
  data->Oz =  0.0e5;
  /* velocity cm/y */
  vx = 1.0*cm_per_yer2m_per_sec;
  vz = 0.25*cm_per_yer2m_per_sec;
  vx_I = 0.5*cm_per_yer2m_per_sec;
  vz_I = 0.25*cm_per_yer2m_per_sec;
  /* rho0 for initial pressure*/
  data->rho0 = 3140.0;
  /*Temperature */
  data->Tbottom = 1364.0;
  data->Ttop    = 0.0;
  data->k_ref   = 3.3;
  data->h_prod  = 1.2e-6;
  data->h_prod_rcp = 4.4e-13;
  data->y_prod  = -40.0e3;
  data->qm      = 20.0e-3;
  data->ylab    = -122.0e3;
  data->Tlitho  = 1300.0;
 
  /* Material constant */
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  ierr = MaterialConstantsSetDefaults(materialconstants);CHKERRQ(ierr);

  DataBucketGetDataFieldByName(materialconstants,EnergyMaterialConstants_classname,&PField);
  DataFieldGetEntries(PField,(void**)&matconstants_e);

  DataBucketGetDataFieldByName(materialconstants,EnergySourceConst_classname,&PField_Q);
  DataFieldGetEntries(PField_Q,(void**)&data_Q);
  
  DataBucketGetDataFieldByName(materialconstants,EnergyConductivityConst_classname,&PField_k);
  DataFieldGetEntries(PField_k,(void**)&data_k);

  rho_ref = 1.0;
  Cp  = 1000.0;
  /* Constant parameters for all phases */
  Tref        = 273.0;
  phi         = 30.0;
  phi_inf     = 15.0;
  Co          = 2.0e7;
  Co_inf      = 2.0e7;
  Hst_cutoff  = 3.0e8;
  Tens_cutoff = 2.0e7;
  eps_min     = 0.0;
  eps_max     = 0.3;
  
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_phi",&phi,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_phi_inf",&phi_inf,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Co",&Co,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Co_inf",&Co_inf,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Hst_cutoff",&Hst_cutoff,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Tens_cutoff",&Tens_cutoff,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_eps_min",&eps_min,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_eps_max",&eps_max,NULL);CHKERRQ(ierr);
  
  phi_rad     = M_PI * phi/180.0;
  phi_inf_rad = M_PI * phi_inf/180.0;
  
  // UPPER CRUST WITH STRIPES OF 4
  // -------------------- //
  /* UPPER CRUST PHASE 0 */
  // ------------------- //  
  region_idx = 0;
  alpha = 2.0e-5;
  beta = 0.0;
  rho_uc = 2700;
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_rho_uc",&rho_uc,NULL);CHKERRQ(ierr);
  MaterialConstantsSetValues_MaterialType(materialconstants,region_idx,VISCOUS_ARRHENIUS_2,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);//DENSITY_BOUSSINESQ);
  MaterialConstantsSetValues_DensityBoussinesq(materialconstants,region_idx,rho_uc,alpha,beta);
  MaterialConstantsSetValues_DensityConst(materialconstants,region_idx,rho_uc);
     
  /* VISCOUS PARAMETERS */
  preexpA_uc  = 6.3e-6;
  Ascale_uc   = 1.e+6;
  entalpy_uc  = 156.e3;
  Vmol_uc     = 0.0;
  nexp_uc     = 2.4;

  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_preexpA_uc",&preexpA_uc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Ascale_uc",&Ascale_uc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_entaply_uc",&entalpy_uc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Vmol_uc",&Vmol_uc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_nexp_uc",&nexp_uc,NULL);CHKERRQ(ierr);
  MaterialConstantsSetValues_ViscosityArrh(materialconstants,region_idx,preexpA_uc,Ascale_uc,entalpy_uc,Vmol_uc,nexp_uc,Tref);
 
  /* PLASTIC PARAMETERS */
  MaterialConstantsSetValues_PlasticDP(materialconstants,region_idx,phi_rad,phi_inf_rad,Co,Co_inf,1.e7,2.e8);//Tens_cutoff,Hst_cutoff);
  MaterialConstantsSetValues_PlasticMises(materialconstants,region_idx,1.e8,1.e8);
  MaterialConstantsSetValues_SoftLin(materialconstants,region_idx,eps_min,eps_max);

  /* ENERGY */
  source_type[0] = ENERGYSOURCE_USE_MATERIALPOINT_VALUE;
  ierr = MaterialConstantsSetValues_EnergyMaterialConstants(region_idx,matconstants_e,alpha,beta,rho_uc,Cp,ENERGYDENSITY_CONSTANT,ENERGYCONDUCTIVITY_CONSTANT,source_type);CHKERRQ(ierr);
  EnergyConductivityConstSetField_k0(&data_k[region_idx],2.7);
//  EnergySourceConstSetField_HeatSource(&data_Q[region_idx],data->h_prod);

  // -------------------- //
  /* LOWER CRUST PHASE 1 */
  // ------------------- //  
  region_idx = 1;
  alpha      = 2.0e-5;
  beta       = 3.e-12;
  rho_lc     = 2800;
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_rho_lc",&rho_lc,NULL);CHKERRQ(ierr);
  MaterialConstantsSetValues_MaterialType(materialconstants,region_idx,VISCOUS_ARRHENIUS_2,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);//DENSITY_BOUSSINESQ); 
  MaterialConstantsSetValues_DensityBoussinesq(materialconstants,region_idx,rho_lc,alpha,beta);
  MaterialConstantsSetValues_DensityConst(materialconstants,region_idx,rho_lc);

  /* VISCOUS PARAMETERS */
  preexpA_lc  = 6.3e-6; //13.4637;
  Ascale_lc   = 1.e+6;
  entalpy_lc  = 156.e3; //345.e+3;
  Vmol_lc     = 0; //38.e-6;
  nexp_lc     = 2.4; //3;

  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_preexpA_lc",&preexpA_lc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Ascale_lc",&Ascale_lc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_entaply_lc",&entalpy_lc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Vmol_lc",&Vmol_lc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_nexp_lc",&nexp_lc,NULL);CHKERRQ(ierr);
  MaterialConstantsSetValues_ViscosityArrh(materialconstants,region_idx,preexpA_lc,Ascale_lc,entalpy_lc,Vmol_lc,nexp_lc,Tref);

  /* PLASTIC PARAMETERS */
  MaterialConstantsSetValues_PlasticDP(materialconstants,region_idx,phi_rad,phi_inf_rad,Co,Co_inf,1.e7,2.e8);//Tens_cutoff,Hst_cutoff);
  MaterialConstantsSetValues_PlasticMises(materialconstants,region_idx,1.e8,1.e8);
  MaterialConstantsSetValues_SoftLin(materialconstants,region_idx,eps_min,eps_max);

/* ENERGY */
  ierr = MaterialConstantsSetValues_EnergyMaterialConstants(region_idx,matconstants_e,alpha,beta,rho_lc,Cp,ENERGYDENSITY_CONSTANT,ENERGYCONDUCTIVITY_CONSTANT,source_type);CHKERRQ(ierr);
  EnergyConductivityConstSetField_k0(&data_k[region_idx],2.8);

  // -------------------------- //
  /* MANTLE LITHOSPHERE PHASE 2 */
  // -------------------------- //  
  region_idx = 2;
  alpha      = 2.0e-5;
  beta       = 3.e-12;
  rho_ml     = 3300;
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_rho_ml",&rho_ml,NULL);CHKERRQ(ierr);
  MaterialConstantsSetValues_MaterialType(materialconstants,region_idx,VISCOUS_ARRHENIUS_2,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);//DENSITY_BOUSSINESQ); 
  MaterialConstantsSetValues_DensityBoussinesq(materialconstants,region_idx,rho_ml,alpha,beta);
  MaterialConstantsSetValues_DensityConst(materialconstants,region_idx,rho_ml);

  /* VISCOUS PARAMETERS */
  preexpA_ml  = 1.1e+5;
  Ascale_ml   = 1.e+6;
  entalpy_ml  = 530.e3;
  Vmol_ml     = 18.e-6;
  nexp_ml     = 3.5;
 
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_preexpA_ml",&preexpA_ml,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Ascale_ml",&Ascale_ml,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_entaply_ml",&entalpy_ml,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Vmol_ml",&Vmol_ml,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_nexp_ml",&nexp_ml,NULL);CHKERRQ(ierr);
  MaterialConstantsSetValues_ViscosityArrh(materialconstants,region_idx,preexpA_ml,Ascale_ml,entalpy_ml,Vmol_ml,nexp_ml,Tref);

  /* PLASTIC PARAMETERS */
  MaterialConstantsSetValues_PlasticDP(materialconstants,region_idx,phi_rad,phi_inf_rad,Co,Co_inf,2.e7,4.e8);//Tens_cutoff,Hst_cutoff);
  MaterialConstantsSetValues_PlasticMises(materialconstants,region_idx,3.e8,3.e8);
  MaterialConstantsSetValues_SoftLin(materialconstants,region_idx,eps_min,eps_max);

/* ENERGY */
  ierr = MaterialConstantsSetValues_EnergyMaterialConstants(region_idx,matconstants_e,alpha,beta,rho_ml,Cp,ENERGYDENSITY_CONSTANT,ENERGYCONDUCTIVITY_CONSTANT,source_type);CHKERRQ(ierr);
  EnergyConductivityConstSetField_k0(&data_k[region_idx],3.3);

  // ---------------------------- //
  /* MANTLE ASTHENOSPHERE PHASE 3 */
  // ---------------------------- //  
  region_idx = 3;
  alpha      = 2.0e-5;
  beta       = 3.e-12;
  rho_ma     = 3300;
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_rho_ma",&rho_ma,NULL);CHKERRQ(ierr);
  MaterialConstantsSetValues_MaterialType(materialconstants,region_idx,VISCOUS_ARRHENIUS_2,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);//DENSITY_BOUSSINESQ); 
  MaterialConstantsSetValues_DensityBoussinesq(materialconstants,region_idx,rho_ma,alpha,beta);
  MaterialConstantsSetValues_DensityConst(materialconstants,region_idx,rho_ma);

  /* VISCOUS PARAMETERS */
  preexpA_ma  = 1.1e+5;
  Ascale_ma   = 1.e+6;
  entalpy_ma  = 530.e3;
  Vmol_ma     = 18.e-6;
  nexp_ma     = 3.5;
 
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_preexpA_ma",&preexpA_ma,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Ascale_ma",&Ascale_ma,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_entaply_ma",&entalpy_ma,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Vmol_ma",&Vmol_ma,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_nexp_ma",&nexp_ma,NULL);CHKERRQ(ierr);
  MaterialConstantsSetValues_ViscosityArrh(materialconstants,region_idx,preexpA_ma,Ascale_ma,entalpy_ma,Vmol_ma,nexp_ma,Tref);

  /* PLASTIC PARAMETERS */
  MaterialConstantsSetValues_PlasticDP(materialconstants,region_idx,phi_rad,phi_inf_rad,Co,Co_inf,2.e7,4.e8);//Tens_cutoff,Hst_cutoff);
  MaterialConstantsSetValues_PlasticMises(materialconstants,region_idx,3.e8,3.e8);
  MaterialConstantsSetValues_SoftLin(materialconstants,region_idx,eps_min,eps_max);

/* ENERGY */
  ierr = MaterialConstantsSetValues_EnergyMaterialConstants(region_idx,matconstants_e,alpha,beta,rho_ma,Cp,ENERGYDENSITY_CONSTANT,ENERGYCONDUCTIVITY_CONSTANT,source_type);CHKERRQ(ierr);
  EnergyConductivityConstSetField_k0(&data_k[region_idx],3.3);
  
  /* Read the options */
  /*cutoff */
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_apply_viscosity_cutoff_global",&rheology->apply_viscosity_cutoff_global,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_eta_lower_cutoff_global",&rheology->eta_lower_cutoff_global,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_eta_upper_cutoff_global",&rheology->eta_upper_cutoff_global,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_runwithmises",&data->runmises,NULL);CHKERRQ(ierr);
  /*scaling */
  nondim = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_nondimensional",&nondim,NULL);CHKERRQ(ierr);
  if (nondim){
    data->dimensional = PETSC_FALSE;
  } else {
    ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_vis_bar",&data->viscosity_bar,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_vel_bar",&data->velocity_bar,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_length_bar",&data->length_bar,NULL);CHKERRQ(ierr);
  }

  /* box geometry, m */
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Lx",&data->Lx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Ly",&data->Ly,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Lz",&data->Lz,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Ox",&data->Ox,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Oy",&data->Oy,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Oz",&data->Oz,NULL);CHKERRQ(ierr);

  /* velocity cm/y */
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_vx",&vx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_vz",&vz,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_vx_I",&vx_I,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_vz_I",&vz_I,NULL);CHKERRQ(ierr);

  /* rho0 for initial pressure */
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_rho0",&data->rho0,NULL);CHKERRQ(ierr);

  /* temperature initial condition */
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Tbot",&data->Tbottom,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_Ttop",&data->Ttop,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_age0",&data->thermal_age0,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_ageAnom",&data->thermal_age_anom,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_wx",&data->wx_anom,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_wz",&data->wz_anom,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_cx",&data->cx_anom,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_cz",&data->cz_anom,NULL);CHKERRQ(ierr);

  /* Material constant */
  for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
    PetscPrintf(PETSC_COMM_WORLD,"reading options");
    ierr= MaterialConstantsSetFromOptions(materialconstants,"model_rift3D_I",regionidx,PETSC_FALSE);CHKERRQ(ierr);
  }

  /*Compute velocity at bottom*/
  Sx = (data->Ly - data->Oy)*(data->Lz - data->Oz);
  Sz = (data->Ly - data->Oy)*(data->Lx - data->Ox);
  Sy = (data->Lx - data->Ox)*(data->Lz - data->Oz);
  vy = (2*vx*Sx-vz*Sz)/Sy;

  /* reports before scaling */
  PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_I_Ox %+1.4e [SI] -model_rift3D_I_Lx : %+1.4e [SI]\n", data->Ox ,data->Lx );
  PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_I_Oy %+1.4e [SI] -model_rift3D_I_Ly : %+1.4e [SI]\n", data->Oy ,data->Ly );
  PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_I_Oz %+1.4e [SI] -model_rift3D_I_Lz : %+1.4e [SI]\n", data->Oz ,data->Lz );
  PetscPrintf(PETSC_COMM_WORLD,"  -model_rift3D_I_vx [m/s]:  %+1.4e  -model_rift3D_I_vz [m/s]:  %+1.4e : computed vy [m/s]:  %+1.4e \n", vx,vz,vy);
  PetscPrintf(PETSC_COMM_WORLD,"-model_rift3D_I_rho0 [kg/m^3] :%+1.4e \n", data->rho0 );
  PetscPrintf(PETSC_COMM_WORLD,"-model_rift3D_I_Tbot:%+1.4e \t -model_rift3D_I_Ttop:%+1.4e \t -model_rift3D_I_age0:%+1.4e \n",data->Tbottom,data->Ttop,  data->thermal_age0);
  PetscPrintf(PETSC_COMM_WORLD,"ageAnom:%+1.4e \t wx:%+1.4e \t wz:%+1.4e cx:%+1.4e \t cz:%+1.4e \n",data->thermal_age_anom,data->wx_anom,data->wz_anom,data->cx_anom,data->cz_anom);

  for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
    MaterialConstantsPrintAll(materialconstants,regionidx);
    MaterialConstantsEnergyPrintAll(materialconstants,regionidx);
  }

  if (data->dimensional) {
    /*Compute additional scaling parameters*/
    data->time_bar      = data->length_bar / data->velocity_bar;
    data->pressure_bar  = data->viscosity_bar/data->time_bar;
    data->density_bar   = data->pressure_bar / data->length_bar;
    data->h_prod_bar    = data->pressure_bar / data->time_bar;
    data->k_bar         = data->pressure_bar*data->length_bar*data->length_bar/data->time_bar;

    PetscPrintf(PETSC_COMM_WORLD,"[rift3D_I]:  during the solve scaling will be done using \n");
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
    data->Lx     = data->Lx / data->length_bar;
    data->Ly     = data->Ly / data->length_bar;
    data->Lz     = data->Lz / data->length_bar;
    data->Ox     = data->Ox / data->length_bar;
    data->Oy     = data->Oy / data->length_bar;
    data->Oz     = data->Oz / data->length_bar;
    data->ylab   = data->ylab / data->length_bar;
    data->y_prod = data->y_prod / data->length_bar;
    //scale velocity
    data->vx   = vx/data->velocity_bar;
    data->vy   = vy/data->velocity_bar;
    data->vz   = vz/data->velocity_bar;
    data->vx_I = vx_I/data->velocity_bar;
    data->vz_I = vz_I/data->velocity_bar;
    //scale rho0
    data->rho0     = data->rho0 / data->density_bar;

    // scale material properties
    for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
      MaterialConstantsScaleAll(materialconstants,regionidx,data->length_bar,data->velocity_bar,data->time_bar,data->viscosity_bar,data->density_bar,data->pressure_bar);
      MaterialConstantsEnergyScaleAll(materialconstants,regionidx,data->length_bar,data->time_bar,
          data->pressure_bar);
    }
    // scale thermal paramatres 
    data->h_prod = data->h_prod / data->h_prod_bar;
    data->h_prod_rcp = data->h_prod_rcp / (data->h_prod_bar/(data->density_bar*(data->pressure_bar/data->density_bar)));
    data->k_ref  = data->k_ref / data->k_bar;
    data->qm     = data->qm / (data->pressure_bar*data->velocity_bar);

    /*Reports scaled values*/

    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_Ox   :  %+1.4e    -model_rift3D_I_Lx   :  %+1.4e  \n", data->Ox ,data->Lx );
    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_Oy   :  %+1.4e    -model_rift3D_I_Ly   :  %+1.4e \n", data->Oy, data->Ly );
    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_Oz   :  %+1.4e    -model_rift3D_I_Lz   :  %+1.4e\n", data->Oz , data->Lz );

    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_Vx:%+1.4e    -model_rift3D_I_vy:%+1.4e    -model_rift3D_I_vz:  %+1.4e \n", data->vx ,data->vy, data->vz);
    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_Vx_I:%+1.4e  -model_rift3D_I_vz_I:%+1.4e \n", data->vx_I , data->vz_I);
    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_rho0:%+1.4e \n", data->rho0 );
    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_h_prod:%+1.4e \n", data->h_prod );
    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_k_ref:%+1.4e \n", data->k_ref );
    PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_rift3D_I_qm:%+1.4e \n", data->qm );
    PetscPrintf(PETSC_COMM_WORLD,"scaled value for material parameters\n");

    for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
      MaterialConstantsPrintAll(materialconstants,regionidx);
      MaterialConstantsEnergyPrintAll(materialconstants,regionidx);
    }
  }
	
  data->use_semi_eulerian_mesh = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_use_semi_eulerian",&data->use_semi_eulerian_mesh,NULL);CHKERRQ(ierr);
  if (data->use_semi_eulerian_mesh) {
    pTatinModel model;

    PetscPrintf(PETSC_COMM_WORLD,"rift3D_I: activating semi Eulerian mesh advection\n");
    ierr = pTatinGetModel(c,&model);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Rift3D_I_semi_eulerian);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Rift3D_I_semi_eulerian);CHKERRQ(ierr);
  }

  data->output_markers = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_output_markers",&data->output_markers,NULL);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
   Returns the parameters and function need to define initial thermal field.
   The function returned can be used to define either the initial condition for T or the boundary condition for T.
 */
PetscBool DMDAVecTraverse3d_InitialContinentalGeothermFunction(PetscScalar pos[], PetscScalar *val,void *ctx)
{
  PetscScalar x,y,z;
  PetscReal   *coeffs;
  PetscBool   impose;
  PetscReal   Tbot,Ttop,Tlitho,h_prod,y_prod,ylab,k,ymin,qm;
  
  /*Get coordinates*/
  x = pos[0];
  y = pos[1];
  z = pos[2];

  coeffs = (PetscReal*)ctx;
    Tbot     = coeffs[0];
    Ttop     = coeffs[1];
    Tlitho   = coeffs[2];
    h_prod   = coeffs[3];
    y_prod   = coeffs[4];
    ylab     = coeffs[5];
    k        = coeffs[6];
    ymin     = coeffs[7];
    qm       = coeffs[8];

  *val = Ttop + qm*(-y)/k + h_prod*pow(y_prod,2)/k * (1-exp(-y/y_prod));

  if (*val >= Tlitho){
	*val = -((Tbot-Tlitho)/(ymin-ylab)) * (ylab-y) + Tlitho;
  }

  impose = PETSC_TRUE;
  return impose;
}

PetscErrorCode ModelRift3D_I_GetDescription_InitialThermalField(ModelRift3D_ICtx *data,PetscReal coeffs[],PetscBool (**func)(PetscScalar*,PetscScalar*,void*) )
{
  PetscFunctionBegin;
  /* assign params */
  coeffs[0] = data->Tbottom;
  coeffs[1] = data->Ttop;
  coeffs[2] = data->Tlitho;
  coeffs[3] = data->h_prod;
  coeffs[4] = data->y_prod;
  coeffs[5] = data->ylab;
  coeffs[6] = data->k_ref;
  coeffs[7] = data->Oy;
  coeffs[8] = data->qm;

  /* assign function to use */
  *func = DMDAVecTraverse3d_InitialContinentalGeothermFunction;

  PetscFunctionReturn(0);
}

PetscBool BCListEvaluator_SplitFace(PetscScalar position[],PetscScalar *value,void *ctx)
{
  PetscBool impose_dirichlet = PETSC_TRUE;
  BC_SplitFace data_ctx = (BC_SplitFace)ctx;

  PetscFunctionBegin;

  //data_ctx = (SplitFaceBCCtx*)ctx;

  if (position[data_ctx->dim]<data_ctx->x0) {
        *value = data_ctx->v0;
        impose_dirichlet = PETSC_TRUE;
  }else if (position[data_ctx->dim]>data_ctx->x1){
        *value = data_ctx->v1;
  }else{
        *value = data_ctx->v1+(position[data_ctx->dim]-data_ctx->x1)*(data_ctx->v0-data_ctx->v1)/(data_ctx->x0-data_ctx->x1);
  }

  return impose_dirichlet;
}

PetscBool BCListEvaluator_Free(PetscScalar position[],PetscScalar *value,void *ctx)
{
  PetscBool impose_dirichlet;

  PetscFunctionBegin;
  
  impose_dirichlet = PETSC_FALSE;
  
  return impose_dirichlet;
}
/*

   1/ Define boundary conditions in one place for this model.

   2/ Calling pattern should always be
   PetscErrorCode ModelRift3D_I_DefineBCList(BCList bclist,DM dav,pTatinCtx user,ModelRift3D_ICtx data)
   where ModelRift3D_ICtx data is a different type for each model.

   3/ Re-use this function in
   ModelApplyBoundaryCondition_Rift3D_I();
   ModelApplyBoundaryConditionMG_Rift3D_I();

*/
PetscErrorCode ModelRift3D_I_DefineBCList(BCList bclist,DM dav,pTatinCtx user,ModelRift3D_ICtx *data)
{
 // SplitFaceBCCtx *ctx = (SplitFaceBCCtx*)ctx;
  BC_SplitFace   bcdata;
  PetscScalar    vy;
  PetscReal      zero = 0.0;
  PetscScalar    Myr2sec = 3.14e13;
  PetscScalar    vxl,vxr,vzf,vzb,vx,vz,vx_I,vz_I,vx_P,vz_P,vx_init,vz_init;
  PetscScalar    time_pause_start,time_pause,time_invert_start,time_invert;
  PetscScalar    xcoord,xcentre;
  PetscScalar    Sx,Sy,Sz;
  PetscBool      strikeslip;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  strikeslip = PETSC_FALSE;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_strikeslip",&strikeslip,NULL);CHKERRQ(ierr);

  time_pause_start = 28;
  time_pause = 48;
  time_invert_start = 49;
  time_invert = 59;
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_time_pause_start",&time_pause_start,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_time_pause",&time_pause,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_time_invert_start",&time_invert_start,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_time_invert",&time_invert,NULL);CHKERRQ(ierr);
 /* Report before scaling */
  PetscPrintf(PETSC_COMM_WORLD,"  BC Time changes:  -model_rift3D_I_time_pause_start %+1.4e [Myr]\n", time_pause_start );
  PetscPrintf(PETSC_COMM_WORLD,"  BC Time changes:  -model_rift3D_I_time_pause %+1.4e [Myr]\n", time_pause );
  PetscPrintf(PETSC_COMM_WORLD,"  BC Time changes:  -model_rift3D_I_time_invert_start %+1.4e [Myr]\n", time_invert_start );
  PetscPrintf(PETSC_COMM_WORLD,"  BC Time changes:  -model_rift3D_I_time_invert %+1.4e [Myr]\n", time_invert );
 
  /* Scaling the time */
  time_pause_start  = time_pause_start * Myr2sec/data->time_bar;
  time_pause        = time_pause * Myr2sec/data->time_bar;
  time_invert_start = time_invert_start * Myr2sec/data->time_bar;
  time_invert       = time_invert * Myr2sec/data->time_bar;
  
  vx_init = data->vx;
  vz_init = data->vz;
  vx_P    = 0.0;
  vz_P    = 0.0;
  Sx = (data->Ly - data->Oy)*(data->Lz - data->Oz);
  Sz = (data->Ly - data->Oy)*(data->Lx - data->Ox);
  Sy = (data->Lx - data->Ox)*(data->Lz - data->Oz);
  
  if(user->time < time_pause_start){
	
	vx = data->vx;
	vz = data->vz;
       	vy = (2*vx*Sx-vz*Sz)/Sy;
        PetscPrintf(PETSC_COMM_WORLD,"Extension: scaled value   -model_rift3D_I_Vx: %+1.4e    -model_rift3D_I_Vz: %+1.4e    -model_rift3D_I_Vy: %+1.4e \n", vx ,vz,vy );
	
  }else if(user->time >= time_pause_start && user->time < time_pause){
	
	vx = ((vx_P-vx_init)/(time_pause-time_pause_start)) * (user->time-time_pause);
	vz = ((vz_P-vz_init)/(time_pause-time_pause_start)) * (user->time-time_pause);
        vy = (2*vx*Sx-vz*Sz)/Sy;
        PetscPrintf(PETSC_COMM_WORLD,"Velocity decreasing: scaled value   -model_rift3D_I_Vx: %+1.4e    -model_rift3D_I_Vz: %+1.4e    -model_rift3D_I_Vy: %+1.4e \n", vx ,vz, vy );

  }else if(user->time >= time_pause && user->time < time_invert_start){
	
	vx = vx_P;
	vz = vz_P;
        vy = (2*vx*Sx-vz*Sz)/Sy;
        PetscPrintf(PETSC_COMM_WORLD,"Relaxing: scaled value   -model_rift3D_I_Vx: %+1.4e    -model_rift3D_I_Vz: %+1.4e    -model_rift3D_I_Vy: %+1.4e \n", vx ,vz, vy );

  }else if(user->time >= time_invert_start && user->time < time_invert){

	vx_I = -data->vx_I;
	vz_I = data->vz_I;
	vx = ((vx_P-vx_I)/(time_invert_start-time_invert)) * (user->time-time_invert_start);
	vz = ((vz_P-vz_I)/(time_invert_start-time_invert)) * (user->time-time_invert_start);
 	vy = (2*vx*Sx-vz*Sz)/Sy;
        PetscPrintf(PETSC_COMM_WORLD,"Velocity increasing: scaled value   -model_rift3D_I_Vx: %+1.4e    -model_rift3D_I_Vz: %+1.4e    -model_rift3D_I_Vy: %+1.4e \n", vx ,vz, vy );

  }else if(user->time >= time_invert){

	vx = -data->vx_I;
	vz = data->vz_I;
 	vy = (2*vx*Sx-vz*Sz)/Sy;
        PetscPrintf(PETSC_COMM_WORLD,"Compression: scaled value   -model_rift3D_I_Vx: %+1.4e    -model_rift3D_I_Vz: %+1.4e    -model_rift3D_I_Vy: %+1.4e \n", vx ,vz, vy );
  }

  if (strikeslip && user->time < time_pause){
        ierr = PetscMalloc(sizeof(struct _p_BC_SplitFace),&bcdata);CHKERRQ(ierr);
	PetscScalar xc_front,xc_back,xc;
	PetscInt    direction;
  	PetscScalar notch_w2,notchspace;
//	PetscReal   x0_bc_back,x0_bc_front,x1_bc_back,x1_bc_front,v0_bc_back,v0_bc_front,v1_bc_back,v1_bc_front;

	/* Define the points which will split the face (based on the notch position and width) */
  	notch_w2   = 50.e3;
  	notchspace = 200.e3;
 	ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_notchspace",&notchspace,NULL);CHKERRQ(ierr);
  	ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_notchwidth",&notch_w2,NULL);CHKERRQ(ierr);
	/* scaling values */
	notch_w2   = notch_w2/data->length_bar;
  	notchspace = notchspace/data->length_bar; 	
	/* Compute the centre of the notch */	
	xc       = (data->Lx + data->Ox)/2.0;
	xc_back  = (data->Lx + data->Ox)/2.0 - notchspace/2.0;
      	xc_front = (data->Lx + data->Ox)/2.0 + notchspace/2.0;
        /* Direction in  which BC will vary */	
	bcdata->dim = 0;
	
	/* Before compression compute vy on bottom face as: */
	vy = vz*(xc_front - xc_back)*(data->Ly - data->Oy)/(data->Lx*data->Lz);
	
	PetscPrintf(PETSC_COMM_WORLD,"Applying BC SplitFace \n" );	
	/* Back Face */
	/* Define the points between which the velocity is lineraly interpolated (avoiding sharp velocity changes) */
	bcdata->x0 = xc_back - notch_w2/2;
        bcdata->x1 = xc_back + notch_w2/2;
	/* v0 is the velocity for position < x0 */
	bcdata->v0 = vz;
	/* v1 is the velocity for position > x1 */
	bcdata->v1 = 0.0;
	/*Vz*/
    	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_SplitFace,(void*)bcdata);CHKERRQ(ierr);
    	/*Vx*/
//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    	/*Vy*/
//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

	/* Front Face */
	bcdata->x0 = xc_front - notch_w2/2;
	bcdata->x1 = xc_front + notch_w2/2;
	/*Vz*/
    	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_SplitFace,(void*)bcdata);CHKERRQ(ierr);
    	/*Vx*/
//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    	/*Vy*/
//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
	/* Vy = 0 along faces of normal x */
//  	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
//  	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	/* Vz = vz along left face of normal x */
  	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&vz);CHKERRQ(ierr);
	/* Vz = 0 along right face of normal x */
  	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
  } else {
  	vzf = -vz;
  	vzb =  0.0;//data->vz;
    	/*compression along face of normal z */
    	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&(vzb));CHKERRQ(ierr);
    	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&(vzf));CHKERRQ(ierr);
//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_Free,(void*)&vz);CHKERRQ(ierr);
//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_Free,(void*)&vz);CHKERRQ(ierr);
	// Free the condition on the faces
// 	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_Free,(void*)&vx);CHKERRQ(ierr);
//  	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_Free,(void*)&vx);CHKERRQ(ierr);
//  	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_Free,(void*)&vx);CHKERRQ(ierr);
//  	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_Free,(void*)&vx);CHKERRQ(ierr);
	
 }
  vxl = -vx;
  vxr =  vx;
  /*extension or compression along face of normal x */
  ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&(vxl));CHKERRQ(ierr);
  ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&(vxr));CHKERRQ(ierr);
  

  /* infilling free slip base */
  ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&(vy));CHKERRQ(ierr);
  /* free surface top*/

    PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyBoundaryCondition_Rift3D_I(pTatinCtx user,void *ctx)
{
  ModelRift3D_ICtx *data = (ModelRift3D_ICtx*)ctx;
  PetscBool        active_energy;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = ModelRift3D_I_DefineBCList(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,user,data);CHKERRQ(ierr);

  /* set boundary conditions for temperature */
  ierr = pTatinContextValid_Energy(user,&active_energy);CHKERRQ(ierr);

  if (active_energy) {
    PetscReal      val_T;
    PhysCompEnergy energy;
    BCList         bclist;
    DM             daT;
    PetscBool      (*iterator_initial_thermal_field)(PetscScalar*,PetscScalar*,void*);
    PetscReal      coeffs[9];

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

PetscErrorCode ModelApplyBoundaryConditionMG_Rift3D_I(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
  ModelRift3D_ICtx *data = (ModelRift3D_ICtx*)ctx;
  PetscInt         n;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  for (n=0; n<nl; n++) {
    ierr = ModelRift3D_I_DefineBCList(bclist[n],dav[n],user,data);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D_I(pTatinCtx c,void *ctx)
{
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]] - Not implemented \n", PETSC_FUNCTION_NAME);
  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D_I_semi_eulerian(pTatinCtx c,void *ctx)
{
  PhysCompStokes  stokes;
  DM              stokes_pack,dav,dap;
  PetscInt        Nxp[2];
  PetscReal       perturb;
  DataBucket      material_point_db,material_point_face_db;
  PetscInt        f, n_face_list=2, face_list[] = { 3, 4 }; // ymin, zmax //
  //  PetscInt        f, n_face_list=1, face_list[] = { 3 }; /* base */
  int             p,n_mp_points;
  MPAccess        mpX;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

#ifndef REMOVE_FACE_INJECTION
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]] FACE_INJECTION IS BEING IGNORED - POTENTIAL BUG DETECTED \n", PETSC_FUNCTION_NAME);
#endif

#ifdef REMOVE_FACE_INJECTION

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
    ierr = MaterialPointSet_viscous_strain(mpX,p,0.0);CHKERRQ(ierr);
    ierr = MaterialPointSet_plastic_strain(mpX,p,0.0);CHKERRQ(ierr);
    ierr = MaterialPointSet_yield_indicator(mpX,p,0);CHKERRQ(ierr);
  }
  ierr = MaterialPointRestoreAccess(material_point_face_db,&mpX);CHKERRQ(ierr);

  /* delete */
  DataBucketDestroy(&material_point_face_db);

#endif

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D_I_semi_eulerian_v2(pTatinCtx c,void *ctx)
{
  PhysCompStokes  stokes;
  DM              stokes_pack,dav,dap;
  PetscInt        Nxp[2];
  PetscReal       perturb;
  DataBucket      material_point_db,material_point_face_db;
  PetscInt        f, n_face_list=2, face_list[] = { 3, 4 }; // ymin, zmax //
  //  PetscInt        f, n_face_list=1, face_list[] = { 3 }; /* base */
  int             p,n_mp_points;
  MPAccess        mpX;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

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

  /* Copy ONLY PHASE from nearest markers to newly inserted markers expect (xi,xip,pid) */
  ierr = MaterialPointRegionAssignment_v2(material_point_db,dav);CHKERRQ(ierr);

  /* delete */
  DataBucketDestroy(&material_point_face_db);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMeshGeometry_Rift3D_I(pTatinCtx c,void *ctx)
{
  ModelRift3D_ICtx *data = (ModelRift3D_ICtx*)ctx;
  PetscErrorCode   ierr;
  PetscBool        use_passive_tracers = PETSC_FALSE;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,data->Ox,data->Lx,data->Oy,data->Ly,data->Oz,data->Lz);
  CHKERRQ(ierr);

  /* note - Don't access the energy mesh here, its not yet created */
  /* note - The initial velocity mesh geometry will be copied into the energy mesh */

  PetscReal gvec[] = { 0.0, -10.0, 0.0 };
  ierr = PhysCompStokesSetGravityVector(c->stokes_ctx,gvec);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_use_passive_tracers",&use_passive_tracers,NULL);CHKERRQ(ierr);
  if (use_passive_tracers){
  ierr = PSwarmSetUp(data->pswarm);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_Rift3D_I(pTatinCtx c,void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = ModelApplyInitialMaterialGeometry_Notchtest_Rift3D_I(c,ctx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyUpdateMeshGeometry_Rift3D_I(pTatinCtx c,Vec X,void *ctx)
{
  PetscReal       step;
  PhysCompStokes  stokes;
  DM              stokes_pack,dav,dap;
  Vec             velocity,pressure;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* fully lagrangian update */
  ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);

  stokes_pack = stokes->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  //ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  //ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,N);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyUpdateMeshGeometry_Rift3D_I_semi_eulerian(pTatinCtx c,Vec X,void *ctx)
{
  ModelRift3D_ICtx *data = (ModelRift3D_ICtx*)ctx;
  PetscReal       step;
  PetscBool       use_passive_tracers = PETSC_FALSE;
  PhysCompStokes  stokes;
  DM              stokes_pack,dav,dap;
  Vec             velocity,pressure;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* fully lagrangian update */
  ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);

  stokes_pack = stokes->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
/* ONLY VERTICAL REMESHING */
 // ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);

 /* SURFACE REMESHING */
  ierr = UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX(dav,velocity,NULL,step);
  ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_use_passive_tracers",&use_passive_tracers,NULL);CHKERRQ(ierr);
  if (use_passive_tracers){
  /* PSwarm update */
  ierr = PSwarmFieldUpdateAll(data->pswarm);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelOutput_Rift3D_I(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
  ModelRift3D_ICtx *data = (ModelRift3D_ICtx*)ctx;
  PetscBool        active_energy;
  PetscBool        use_passive_tracers = PETSC_FALSE;
  DataBucket       materialpoint_db;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  //ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
  // just plot the velocity field (coords and vel stored in file as floats)
  ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(c,X,prefix);CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

  if (data->output_markers)
  {
    ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
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
    const int                   nf = 5;
    const MaterialPointVariable mp_prop_list[] = { MPV_region, MPV_viscosity, MPV_density, MPV_viscous_strain, MPV_plastic_strain };
	
    ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
  }

  /* standard viewer */
  ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
  if (active_energy) {
    PhysCompEnergy energy;
    Vec            temperature;

    ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
    ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);

    ierr = pTatin3d_ModelOutput_Temperature_Energy(c,temperature,prefix);CHKERRQ(ierr);
  }

  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_use_passive_tracers",&use_passive_tracers,NULL);CHKERRQ(ierr);
  if (use_passive_tracers){
  /* Dump the PSwarm */
  ierr = PSwarmView(data->pswarm,PSW_VT_SINGLETON);CHKERRQ(ierr);
  ierr = PSwarmViewInfo(data->pswarm);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelDestroy_Rift3D_I(pTatinCtx c,void *ctx)
{
  ModelRift3D_ICtx *data = (ModelRift3D_ICtx*)ctx;
  PetscBool        use_passive_tracers = PETSC_FALSE;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* Free contents of structure */
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_use_passive_tracers",&use_passive_tracers,NULL);CHKERRQ(ierr);
  if(use_passive_tracers){
  /* destroy PSwarm*/
  ierr = PSwarmDestroy(&data->pswarm);CHKERRQ(ierr);
  }

  /* Free structure */
  ierr = PetscFree(data);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialStokesVariableMarkers_Rift3D_I(pTatinCtx user,Vec X,void *ctx)
{
  DM                         stokes_pack,dau,dap;
  PhysCompStokes             stokes;
  Vec                        Uloc,Ploc;
  PetscScalar                *LA_Uloc,*LA_Ploc;
  ModelRift3D_ICtx           *data = (ModelRift3D_ICtx*)ctx;
  DataField                  PField;
  MaterialConst_MaterialType *truc;
  PetscInt                   regionidx;
  PetscErrorCode             ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

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

PetscErrorCode ModelApplyInitialCondition_Rift3D_I(pTatinCtx c,Vec X,void *ctx)
{
  ModelRift3D_ICtx                             *data = (ModelRift3D_ICtx*)ctx;
  DM                                           stokes_pack,dau,dap;
  Vec                                          velocity,pressure;
  PetscReal                                    vxl,vxr,vzb,vzf,vy;
  DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
  DMDAVecTraverse3d_InterpCtx                  IntpCtx;
  PetscReal                                    MeshMin[3],MeshMax[3],domain_height;
  PetscBool                                    active_energy;
  PetscBool                                    use_passive_tracers = PETSC_FALSE;
  PetscErrorCode                               ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  stokes_pack = c->stokes_ctx->stokes_pack;

  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  vxl = -data->vx;
  vxr =  data->vx;

  vzf = -data->vz;
  vzb =  0.0;//data->vz;

  ierr = VecZeroEntries(velocity);CHKERRQ(ierr);

  ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vxr-vxl)/(data->Lx-data->Ox),vxl,0.0);CHKERRQ(ierr);
  ierr = DMDAVecTraverse3d(dau,velocity,0,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
  ierr = DMDAVecTraverse3d_InterpCtxSetUp_Z(&IntpCtx,(vzf-vzb)/(data->Lz-data->Oz),vzb,0.0);CHKERRQ(ierr);
  ierr = DMDAVecTraverse3d(dau,velocity,2,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);

  vy= data->vy;
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
    PetscBool      (*iterator_initial_thermal_field)(PetscScalar*,PetscScalar*,void*);
    PetscReal      coeffs[9];

    ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
    ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);
    daT  = energy->daT;

    ierr = ModelRift3D_I_GetDescription_InitialThermalField(data,coeffs,&iterator_initial_thermal_field);CHKERRQ(ierr);
    ierr = DMDAVecTraverse3d(daT,temperature,0,iterator_initial_thermal_field,(void*)coeffs);CHKERRQ(ierr);
  }

  /* set swarm velocity and pressure vector once*/
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_use_passive_tracers",&use_passive_tracers,NULL);CHKERRQ(ierr);
  if(use_passive_tracers){
  ierr = PSwarmAttachStateVecVelocityPressure(data->pswarm,X);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinModelRegister_Rift3D_I(void)
{
  ModelRift3D_ICtx *data;
  pTatinModel      m;
  PetscErrorCode   ierr;

  PetscFunctionBegin;

  /* Allocate memory for the data structure for this model */
  ierr = PetscMalloc(sizeof(ModelRift3D_ICtx),&data);CHKERRQ(ierr);
  ierr = PetscMemzero(data,sizeof(ModelRift3D_ICtx));CHKERRQ(ierr);

  /* register user model */
  ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

  /* Set name, model select via -ptatin_model NAME */
  ierr = pTatinModelSetName(m,"rift3D_I");CHKERRQ(ierr);

  /* Set model data */
  ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);

  /* Set function pointers */
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Rift3D_I);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialCondition_Rift3D_I);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,   (void (*)(void))ModelApplyInitialStokesVariableMarkers_Rift3D_I);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Rift3D_I);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Rift3D_I);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Rift3D_I);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Rift3D_I);CHKERRQ(ierr);

  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Rift3D_I);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Rift3D_I);CHKERRQ(ierr);

  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Rift3D_I);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Rift3D_I);CHKERRQ(ierr);

  /* Insert model into list */
  ierr = pTatinModelRegister(m);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_Notchtest_Rift3D_I(pTatinCtx c,void *ctx)
{
  ModelRift3D_ICtx *data = (ModelRift3D_ICtx*)ctx;
  int              p,n_mp_points;
  PetscScalar      y_lab,y_moho,y_midcrust,notch_l,notch_w2,xc,notchspace;
  DataBucket       db;
  DataField        PField_std,PField_pls,PField_energy;
  int              phase;
  PetscBool        norandomiseplastic,double_notch;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* define properties on material points */
  db = c->materialpoint_db;
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  DataFieldVerifyAccess(PField_std,sizeof(MPntStd));

  DataBucketGetDataFieldByName(db,MPntPStokesPl_classname,&PField_pls);
  DataFieldGetAccess(PField_pls);
  DataFieldVerifyAccess(PField_pls,sizeof(MPntPStokesPl));
  
  /* Get bucket for material points energy */
  DataBucketGetDataFieldByName(db,MPntPEnergy_classname,&PField_energy);
  DataFieldGetAccess(PField_energy);
  DataFieldVerifyAccess(PField_energy,sizeof(MPntPEnergy));


  /* m */
  y_lab      = -120.0e3;
  y_moho     = -40.0e3;
  y_midcrust = -20.0e3;
  notch_w2   = 50.e3;
  notch_l    = 120.e3;
  xc         = (data->Lx + data->Ox)/2.0* data->length_bar;
  notchspace = 200.e3;
  //xc         = 0.0;
  DataBucketGetSizes(db,&n_mp_points,0,0);

  ptatin_RandomNumberSetSeedRank(PETSC_COMM_WORLD);

  double_notch = PETSC_FALSE;
  norandomiseplastic = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_norandom",&norandomiseplastic,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_rift3D_I_DoubleNotch",&double_notch,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_notchspace",&notchspace,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_notchwidth",&notch_w2,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_rift3D_I_notchlength",&notch_l,NULL);CHKERRQ(ierr);

  for (p=0; p<n_mp_points; p++) {
    MPntStd       *material_point;
    MPntPStokesPl *mpprop_pls;
    MPntPEnergy   *mpp_energy;
    PetscScalar   prod;
    double        *position,ycoord,xcoord,zcoord;
    float         pls;
    char          yield;

    DataFieldAccessPoint(PField_std,   p,(void**)&material_point);
    DataFieldAccessPoint(PField_pls,   p,(void**)&mpprop_pls);
    DataFieldAccessPoint(PField_energy,p,(void**)&mpp_energy);

    /* Access using the getter function provided for you (recommeneded for beginner user) */
    MPntStdGetField_global_coord(material_point,&position);
    
    /* Heat source production decreasing exponentially with depth */
    prod = data->h_prod*exp(-position[1]/data->y_prod);

    /* convert to scaled units */
    xcoord = position[0] * data->length_bar;
    ycoord = position[1] * data->length_bar;
    zcoord = position[2] * data->length_bar;

    if (ycoord < y_lab) {
      phase = 3;
    } else if (ycoord < y_moho) {
      phase = 2;
    } else if (ycoord < y_midcrust) {
      phase = 1;
    } else {
      phase = 0;
    }

    if (!double_notch){
      if (norandomiseplastic) {
        pls   = 0.0;
        if ( (fabs(xcoord - xc) < notch_w2) && (zcoord < notch_l) && (ycoord > y_lab) ) {
          pls = 0.05;
        }
      } else {
        pls = ptatin_RandomNumberGetDouble(0.0,0.03);
        if ( (fabs(xcoord - xc) < notch_w2) && (zcoord < notch_l) && (ycoord > y_lab) ) {
          pls = ptatin_RandomNumberGetDouble(0.0,0.3);
        }
      }
    }else{
      PetscScalar xc1,xc2,Lz;
      xc1 = (data->Lx + data->Ox)/2.0* data->length_bar - notchspace/2.0;
      xc2 = (data->Lx + data->Ox)/2.0* data->length_bar + notchspace/2.0;
      Lz  = data->Lz*data->length_bar;

      if (norandomiseplastic) {
        pls   = 0.0;
        if ( (fabs(xcoord - xc1) < notch_w2) && (zcoord < notch_l) && (ycoord > y_lab) )  {
          pls = 0.05;
        }
        if ( (fabs(xcoord - xc2) < notch_w2) && (zcoord > Lz-notch_l) && (ycoord > y_lab) )  {
          pls = 0.05;
        }
      } else {

        pls = ptatin_RandomNumberGetDouble(0.0,0.03);
        if ( (fabs(xcoord - xc1) < notch_w2) && (zcoord < notch_l) && (ycoord > y_lab) )  {
          pls = ptatin_RandomNumberGetDouble(0.0,0.3);
        }
        if ( (fabs(xcoord - xc2) < notch_w2) && (zcoord > Lz-notch_l) && (ycoord > y_lab) ) {

          pls = ptatin_RandomNumberGetDouble(0.0,0.3);
        }
      }
    }

    yield = 0;
    /* user the setters provided for you */
    MPntStdSetField_phase_index(material_point,phase);
    MPntPStokesPlSetField_yield_indicator(mpprop_pls,yield);
    MPntPStokesPlSetField_plastic_strain(mpprop_pls,pls);
    MPntPEnergySetField_heat_source(mpp_energy,prod);
  }

  DataFieldRestoreAccess(PField_std);
  DataFieldRestoreAccess(PField_pls);
  DataFieldRestoreAccess(PField_energy);
  
  PetscFunctionReturn(0);
}

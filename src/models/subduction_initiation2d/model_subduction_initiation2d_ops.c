
/*
 Developed by 
 Xiaolin Mao [xlmao@gps.caltech.edu]
 Mike Gurnis [gurnis@gps.caltech.edu]
 */


#define _GNU_SOURCE
#include "petsc.h"
#include "petscmath.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_utils.h"
#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_iterator.h"
#include "mesh_update.h"
#include "output_material_points.h"
#include "material_point_std_utils.h"
#include "energy_output.h"
#include "geometry_object.h"
#include "ptatin3d_energy.h"
#include "model_subduction_initiation2d_ctx.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Subduction_Initiation2d"
PetscErrorCode ModelInitialize_Subduction_Initiation2d(pTatinCtx c,void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	PetscReal cm_yr2m_s;
	RheologyConstants *rheology;
	DataBucket materialconstants = c->material_constants;
	PetscInt regionidx;
	
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	cm_yr2m_s = 1.0e-2 / ( 365.25 * 24.0 * 60.0 * 60.0 ) ;
	
	
	/*box geometry*/
	data->Lx = 12.5e5;
	data->Ly = 0.0e5;
	data->Lz = 0.01e5;
	data->Ox = 0.0;
	data->Oy = -3.0e5;
	data->Oz = 0.0e5;
	
	/*material constants*/	
	data->diffusivity[0] = 1.0e-6;
	data->alpha[0] = 3.0e-2;
	data->H[0] = 0.0;/*heat generation*/
	
	/*density*/
	data->rho[0] = 3300.0;
	
	/*plastic constants*/
	data->C0[0] = 1.0e8;
	data->mu[0] = 0.6;
	//	data->mu[0]=atan(data->mu[0]);/*internal angle of firction*/
	data->C0_inf[0] = 5.0e7;
	data->mu_inf[0] = 0.2;
	//	data->mu_inf[0]=atan(data->mu_inf[0]);
	
	
	
	/*velocity boundary*/
	data->velocity = 2.0;
	
	
	/*temperature*/
	data->Ttop = 0.0;
	data->Tbot = 1.4000; //	
	data->Thermal_age = 20.0; //Ma
	
	/*viscousity*/
	data->eta[0] = 5.3357e26;
	
	/*scaling bar*/
	data->length_bar    = 1000.0 * 1.0e3;
	data->viscosity_bar = 1.0e25;
	data->velocity_bar  = 1.0e-9;
	
	
	rheology                 = &c->rheology_constants;
	rheology->rheology_type  = RHEOLOGY_VP_STD;
	rheology->nphases_active = 4;
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global = 1.e+26;
	rheology->eta_lower_cutoff_global = 1.e+20;
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Lx",&data->Lx,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Ox",&data->Ox,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Ly",&data->Ly,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Oy",&data->Oy,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Lz",&data->Lz,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Oz",&data->Oz,PETSC_NULL);CHKERRQ(ierr);
	
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_velocity",&data->velocity,PETSC_NULL);CHKERRQ(ierr); 	
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_density_bar",&data->density_bar,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_length_bar",&data->length_bar,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_viscosity_bar",&data->viscosity_bar,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_velocity_bar",&data->velocity_bar,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Ttop",&data->Ttop,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Tbot",&data->Tbot,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_Thermal_age",&data->Thermal_age,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_eta0",&data->eta[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_rho0",&data->rho[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_diffusivity0",&data->diffusivity[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_alpha0",&data->alpha[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_H0",&data->H[0],PETSC_NULL);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_eta_upper_cutoff",&rheology->eta_upper_cutoff_global,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_SI_eta_lower_cutoff",&rheology->eta_lower_cutoff_global,PETSC_NULL);CHKERRQ(ierr);
	
	
	data->velocity = data->velocity*cm_yr2m_s;
	
	
	MaterialConstantsSetDefaults(materialconstants);
	
	MaterialConstantsSetValues_MaterialType(materialconstants,0,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	MaterialConstantsSetValues_ViscosityFK(materialconstants,0,data->eta[0],12.0);
	MaterialConstantsSetValues_ViscosityConst(materialconstants,0,data->eta[0]);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,0,data->rho[0],data->alpha[0],0.0);
	MaterialConstantsSetValues_DensityConst(materialconstants,0,data->rho[0]);
#if 1       
	MaterialConstantsSetValues_PlasticDP(materialconstants,0,data->mu[0],data->mu_inf[0],data->C0[0],data->C0_inf[0],data->C0_inf[0],3.5*data->C0[0]);
	//       MaterialConstantsSetValues_PlasticMises(materialconstants,0,10.0*data->C0[0],10.0*data->C0[0]);
	MaterialConstantsSetValues_SoftLin(materialconstants,0,0.1,0.6);
#endif
	MaterialConstantsSetValues_MaterialType(materialconstants,1,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	MaterialConstantsSetValues_ViscosityFK(materialconstants,1,data->eta[0],12.0);
	MaterialConstantsSetValues_ViscosityConst(materialconstants,1,data->eta[0]);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,1,data->rho[0],data->alpha[0],0.0);
	MaterialConstantsSetValues_DensityConst(materialconstants,1,data->rho[0]); 
	
#if 1       
	MaterialConstantsSetValues_PlasticDP(materialconstants,1,data->mu[0],data->mu_inf[0],0.1*data->C0[0],0.1*data->C0_inf[0],0.1*data->C0_inf[0],0.1*data->C0[0]);
	MaterialConstantsSetValues_SoftLin(materialconstants,1,0.0,0.2);
#endif      
	
	MaterialConstantsSetValues_MaterialType(materialconstants,2,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	MaterialConstantsSetValues_ViscosityFK(materialconstants,2,data->eta[0],12.0);
	MaterialConstantsSetValues_ViscosityConst(materialconstants,2,data->eta[0]);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,2,data->rho[0],data->alpha[0],0.0);
	MaterialConstantsSetValues_DensityConst(materialconstants,2,data->rho[0]);
#if 1       
	MaterialConstantsSetValues_PlasticDP(materialconstants,2,data->mu[0],data->mu_inf[0],0.1*data->C0[0],0.1*data->C0_inf[0],0.1*data->C0_inf[0],0.1*data->C0[0]);
	MaterialConstantsSetValues_SoftLin(materialconstants,2,0.0,0.2);
#endif
	
	MaterialConstantsSetValues_MaterialType(materialconstants,3,VISCOUS_FRANKK,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_BOUSSINESQ);
	MaterialConstantsSetValues_ViscosityFK(materialconstants,3,data->eta[0],12.0);
	MaterialConstantsSetValues_ViscosityConst(materialconstants,3,data->eta[0]);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,3,data->rho[0],data->alpha[0],0.0);
	MaterialConstantsSetValues_DensityConst(materialconstants,3,data->rho[0]);
#if 1       
	MaterialConstantsSetValues_PlasticDP(materialconstants,3,data->mu[0],data->mu_inf[0],0.1*data->C0[0],0.1*data->C0_inf[0],0.1*data->C0_inf[0],0.1*data->C0[0]);
	MaterialConstantsSetValues_SoftLin(materialconstants,3,0.0,0.2);
#endif
	
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) { 
		MaterialConstantsPrintAll(materialconstants,regionidx);
	} 
	
	/*Compute additional scaling parameters*/
	data->time_bar      = data->length_bar / data->velocity_bar;
	data->pressure_bar  = data->viscosity_bar/data->time_bar;
	data->density_bar   = data->pressure_bar / data->length_bar;
	
	
	
	PetscPrintf(PETSC_COMM_WORLD,"[subduction_initiation2d]:  during the solve scaling will be done using \n");
	PetscPrintf(PETSC_COMM_WORLD,"  L*    : %1.4e [m]\n", data->length_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  U*    : %1.4e [m.s^-1]\n", data->velocity_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  t*    : %1.4e [s]\n", data->time_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  eta*  : %1.4e [Pa.s]\n", data->viscosity_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  rho*  : %1.4e [kg.m^-3]\n", data->density_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  P*    : %1.4e [Pa]\n", data->pressure_bar );
	
	PetscPrintf(PETSC_COMM_WORLD,"[subduction_initiation2d]:  parameters \n");
	PetscPrintf(PETSC_COMM_WORLD,"  alpha    : %1.4e [K^-1]\n", data->alpha[0] );
	PetscPrintf(PETSC_COMM_WORLD,"  Thermal_age   : %1.4e [Ma]\n", data->Thermal_age );
	PetscPrintf(PETSC_COMM_WORLD,"  eta_upper_cutoff   : %1.4e [Pas]\n", rheology->eta_upper_cutoff_global );
	PetscPrintf(PETSC_COMM_WORLD,"  eta_lower_cutoff   : %1.4e [Pas]\n", rheology->eta_lower_cutoff_global );
	
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
	data->velocity = data->velocity/data->velocity_bar;
	
	//scale rho0
	data->rho[0] = data->rho[0]/data->density_bar;
	
	data->diffusivity[0]=data->diffusivity[0]/data->length_bar/data->length_bar*data->time_bar;
#if 1 	
	// scale material properties
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
		MaterialConstantsScaleAll(materialconstants,regionidx,data->length_bar,data->velocity_bar,data->time_bar,data->viscosity_bar,data->density_bar,data->pressure_bar);
	}
	
	ierr = PetscOptionsInsertString("-activate_energy");CHKERRQ(ierr);
#endif	
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_X_Subduction_Initiation2d"
PetscBool BCListEvaluator_X_Subduction_Initiation2d( PetscScalar position[], PetscScalar *value, void *ctx ) 
{
	PetscBool impose_dirichlet = PETSC_TRUE;
	pTatinCtx user = (pTatinCtx)ctx;
	PetscScalar vx,Vx_max;
	PetscReal Dprime;
	ModelSubduction_Initiation2dCtx *model_data_ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	ierr = pTatinModelGetUserData(user->model,(void**)&model_data_ctx);CHKERRQ(ierr);
	
	Dprime = (model_data_ctx->Ly-model_data_ctx->Oy)/2.0;
	//	Dprime = model_data_ctx->Ly/2.0;
	Vx_max = model_data_ctx->velocity;
	
	vx=-Vx_max*tanh(2.0*3.1415926*(position[1]-model_data_ctx->Oy-Dprime)/(2.0*Dprime));
#if 0
	if(position[1]>=-0.07){
		vx=-Vx_max*tanh(2.0*3.1415926*(position[1]+0.07)/(2.0*0.07));
	}else if(position[1]<=-0.18){
		vx=-Vx_max*tanh(2.0*3.1415926*(position[1]+0.18)/(2.0*0.07));
	}else{
		vx=0.0;
	}
	
#endif
	
	//        vx=-Vx_max*(tanh(2.0*3.1415926*(position[1]-model_data_ctx->Oy-Dprime)/(2.0*Dprime))+(1.0+position[1]/Dprime))/2.0; 
	
	//        vx=-Vx_max*(1.0+position[1]/Dprime);
	
	//	vx=-Vx_max*tanh(2.0*3.1415926*(position[1]-Dprime)/(2.0*Dprime));
	*value = vx;
	
	return impose_dirichlet;
}

#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_Y_Subduction_Initiation2d"
PetscBool BCListEvaluator_Y_Subduction_Initiation2d( PetscScalar position[], PetscScalar *value, void *ctx ) 
{
	PetscBool impose_dirichlet = PETSC_TRUE;
	pTatinCtx user = (pTatinCtx)ctx;
	PetscScalar vx,Vx_max;
	PetscReal Dprime;
	ModelSubduction_Initiation2dCtx *model_data_ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	ierr = pTatinModelGetUserData(user->model,(void**)&model_data_ctx);CHKERRQ(ierr);
	
	Dprime = (model_data_ctx->Lx-model_data_ctx->Ox);
	Vx_max = model_data_ctx->velocity;
	
	vx=0.0+(position[0]/Dprime)*Vx_max;
	//	vx=-Vx_max*tanh(2.0*3.1415926*(position[1]-Dprime)/(2.0*Dprime));
	*value = vx;
	
	return impose_dirichlet;
}

#undef __FUNCT__
#define __FUNCT__ "Subduction_Initiation2d_VelocityBC"
PetscErrorCode Subduction_Initiation2d_VelocityBC(BCList bclist,DM dav,pTatinCtx c,ModelSubduction_Initiation2dCtx *data)
{
	PetscScalar zero=0.0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
#if 1	
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_X_Subduction_Initiation2d,(void*)c);CHKERRQ(ierr);
	//	PetscScalar vxtemp;
	//	vxtemp=-data->velocity;
	//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
	//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
	
	//	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_Y_Subduction_Initiation2d,(void*)c);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
	
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
#endif	
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Subduction_Initiation2d"
PetscErrorCode ModelApplyBoundaryCondition_Subduction_Initiation2d(pTatinCtx c,void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	PetscScalar zero=0.0,velocity;
	PetscErrorCode    ierr;
	PetscBool active_energy;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = Subduction_Initiation2d_VelocityBC(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,c,data);CHKERRQ(ierr);
	
	
	
#if 1	
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		
		PetscReal      val_T;
		PhysCompEnergy energy;
		BCList         bclist;
		DM             daT;
		
		ierr   = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		daT    = energy->daT;
		bclist = energy->T_bclist;
		
		
		val_T = data->Tbot;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);
		
		val_T = data->Ttop;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);		
		
	}
#endif
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Subduction_Initiation2d"
PetscErrorCode ModelApplyBoundaryConditionMG_Subduction_Initiation2d(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	PetscInt n;
	PetscScalar zero = 0.0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = Subduction_Initiation2d_VelocityBC(bclist[n],dav[n],c,data);CHKERRQ(ierr);
		
	}	
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Subduction_Initiation2d"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Subduction_Initiation2d(pTatinCtx c,void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Subduction_Initiation2d"
PetscErrorCode ModelApplyInitialMeshGeometry_Subduction_Initiation2d(pTatinCtx c,void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	PetscReal    Lx,Ly,Lz;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, data->Ox,data->Lx, data->Oy,data->Ly, data->Oz, data->Lz);CHKERRQ(ierr);
	//	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, 0.0,data->Lx, 0.0,data->Ly, 0.0, data->Lz);CHKERRQ(ierr);	
	
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Subduction_Initiation2d"
PetscErrorCode ModelApplyInitialMaterialGeometry_Subduction_Initiation2d(pTatinCtx c,void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_pls,PField_stokes, PField_energy;
	int                    phase;
	PetscBool              active_energy;
	MPAccess               mpX;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
#if 1	
	
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
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd       *material_point;
		MPntPStokes   *mpprop_stokes;
		MPntPStokesPl *mpprop_pls;
		double        *position,ycoord,xcoord,zcoord;
		float         pls;
		char          yield;
		
		DataFieldAccessPoint(PField_std,p,(void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		DataFieldAccessPoint(PField_pls,p,(void**)&mpprop_pls);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		xcoord = position[0] * data->length_bar;
		ycoord = position[1] * data->length_bar;
		zcoord = position[2] * data->length_bar;
		
		phase = 0;
		
		pls   = 0.0;		
		yield = 0; 
#if 0  
		if(ycoord<=0.0&&ycoord>=-0.15e5){
			if(xcoord>=ycoord+6.25e5&&xcoord<=ycoord+6.35e5){
				phase = 1;
				
			}
		}  
#endif	
#if 0  
		if(xcoord>=6.25e5&&xcoord<=6.45e5){
			if(ycoord<=0.0&&ycoord>=-0.20e5){
				phase = 1;
				
				if(xcoord>=6.27e5&&xcoord<=6.43e5){
					if(ycoord<=0.0&&ycoord>=-0.18e5){
						phase = 2;
						
						if(xcoord>=6.29e5&&xcoord<=6.41e5){
							if(ycoord<=0.0&&ycoord>=-0.16e5){
								phase = 3;
							}
						}
						
					}
				}
				pls = ptatin_RandomNumberGetDouble(0.0,0.1);
			}
		}  
		
		
#endif		
		
#if 1       
		if(ycoord<=0.0&&ycoord>=-0.70e5){
			if(xcoord>=ycoord+6.25e5&&xcoord<=ycoord+6.45e5){
				phase = 1;
				
				
				pls = ptatin_RandomNumberGetDouble(0.0,0.1);
			}
		}  
#endif		
#if 0        
		if(ycoord<=0.0&&ycoord>=-0.1e5){
			if(xcoord>=ycoord+6.20e5&&xcoord<=ycoord+6.30e5){
				phase = 1;
				
				if(ycoord<=0.0&&ycoord>=-0.08e5){
					if(xcoord>=ycoord+6.22e5&&xcoord<=ycoord+6.38e5){                   
						phase = 2;
						
						if(ycoord<=0.0&&ycoord>=-0.06e5){
							if(xcoord>=ycoord+6.24e5&&xcoord<=ycoord+6.26e5){
								phase = 3;
							}
						}
						
					}
				}
				pls = ptatin_RandomNumberGetDouble(0.0,0.1);
			}
		}  
#endif		
		
#if 0        
		if(ycoord<=0.0&&ycoord>=-0.7e5){
			if(xcoord>=1.0*ycoord+6.15e5&&xcoord<=1.0*ycoord+6.35e5){
				phase = 1;
				
				//                   if(ycoord<=0.0&&ycoord>=-0.15e5){
				if(xcoord>=1.0*ycoord+6.18e5&&xcoord<=1.0*ycoord+6.32e5){                   
					phase = 2;
					
					//                       if(ycoord<=0.0&&ycoord>=-0.10e5){
					if(xcoord>=1.0*ycoord+6.21e5&&xcoord<=1.0*ycoord+6.29e5){
						phase = 3;
					}
					//                       }
					
				}
				//                    }
				
			}
		}  
#endif		
#if 0        
		if(ycoord<=0.0&&ycoord>=-0.55e5){
			if(xcoord>=2.0*ycoord+5.75e5&&xcoord<=5.35e5-0.35e5*log(1.04979-(ycoord+0.35e5)/0.35e5)){
				phase = 1;
				
				//                   if(ycoord<=0.0&&ycoord>=-0.15e5){
				if(xcoord>=2.0*ycoord+5.77e5&&xcoord<=5.33e5-0.35e5*log(1.04979-(ycoord+0.35e5)/0.35e5)){                   
					phase = 2;
					
					//                       if(ycoord<=0.0&&ycoord>=-0.10e5){
					if(xcoord>=2.0*ycoord+5.79e5&&xcoord<=5.31e5-0.35e5*log(1.04979-(ycoord+0.35e5)/0.35e5)){
						phase = 3;
					}
					//                       }
					
				}
				//                    }
				
			}
		}  
#endif		
#if 0
		double  eta_mp;
		PetscScalar temperature;
		
#if 1     
		PetscReal agesin[100],a;
		PetscInt i,i0;
		PetscScalar x,y,z;
		
		x = position[0];
		y = position[1];
		z = position[2];
		
		a=0.878132429842795;
		i0=0;
		for (i=0; i<100; i++) {
			agesin[i]=1.0/3.1536e7*pow(data->length_bar*(0.05+0.005*sin(6.0*3.1415926/1.25*(i*1.25/99)))/(2.0*a),2);
			if (x>=i*1.25/99.0&&x<=(i+1)*1.25/99.0){
				i0=i;
			}
		}
		temperature = data->Ttop+(data->Tbot-data->Ttop)*(1.0-erfc((-y)*data->length_bar/2.0/sqrt(1.e-6*agesin[i0]*3.1536e13)));
#endif		
		eta_mp=data->eta[0]/data->viscosity_bar*exp(-6.0*temperature);
		MPntPStokesSetField_eta_effective(mpprop_stokes,eta_mp);                               
#endif
		
		
		
		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,phase);
		MPntPStokesPlSetField_yield_indicator(mpprop_pls,yield);
		MPntPStokesPlSetField_plastic_strain(mpprop_pls,pls);
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	DataFieldRestoreAccess(PField_pls);
	
	ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mp_points; p++) {
		double kappa,H;
		
		ierr = MaterialPointGet_phase_index(mpX,p,&phase);CHKERRQ(ierr);
		
		kappa = data->diffusivity[0];
		H     = data->H[0];
		ierr = MaterialPointSet_diffusivity(mpX,p,kappa);CHKERRQ(ierr);
		ierr = MaterialPointSet_heat_source(mpX,p,H);CHKERRQ(ierr);
	}
	ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);
	
#endif
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_ERFC3DFunctionSubduction_Initiation2d"
PetscBool DMDAVecTraverse3d_ERFC3DFunctionSubduction_Initiation2d( PetscScalar position[], PetscScalar *val, void *ctx ) 
{
	PetscScalar x,y,z;
	PetscReal  *coeffs;
	PetscBool  impose;
	PetscReal  L_bar,Tbot,age,Ttop,Ly,Oy,length_bar,diffusivity;
	PetscReal cm_yr2m_s;
	
	cm_yr2m_s=1.0e-2 / ( 365.25 * 24.0 * 60.0 * 60.0 ) ;
	
	
	/* get coordinates */
	x = position[0];
	y = position[1];
	z = position[2];
	/* fetch user data */
	coeffs = (PetscReal*)ctx;
	Ttop       = coeffs[0]; 
	Tbot       = coeffs[1]; 
	age        = coeffs[2];
	Ly         = coeffs[3];
	length_bar = coeffs[4];
	diffusivity= coeffs[5];
	
	/*ridge*/   
	//       if(x>=0.61&&x<=0.64){
	//        *val = Ttop+(Tbot-Ttop)*(1.0-erfc((-y)*length_bar/2.0/sqrt(1.e-6*length_bar*(cabs(0.64-0.625)+1.e-16)/(2.0*cm_yr2m_s))));
	//       }else{
	//        *val = Ttop+(Tbot-Ttop)*(1.0-erfc((-y)*length_bar/2.0/sqrt(1.e-6*length_bar*(cabs(x-0.625)+1.e-16)/(2.0*cm_yr2m_s))));
	//        *val = Ttop+(Tbot-Ttop)*(1.0-erfc((-y)*length_bar/2.0/sqrt(1.e-6*length_bar*(cabs(x-0.625)+0.02)/(2.0*cm_yr2m_s))));
	//       }
	
	
#if 1
	if(x<=0.625){
		*val = Ttop+(Tbot-Ttop)*(1.0-erfc((-y)*length_bar/2.0/sqrt(1.e-6*30.0*3.1536e13)));
	}
	else{
		*val = Ttop+(Tbot-Ttop)*(1.0-erfc((-y)*length_bar/2.0/sqrt(1.e-6*30.0*3.1536e13)));
	}
#endif
	/*homogenous*/
	//	*val = Ttop+(Tbot-Ttop)*(1.0-erfc((-y)*length_bar/2.0/sqrt(1.e-6*age*3.1536e13)));//
	//        *val = Ttop+(Tbot-Ttop)*(1.0-erfc((Ly-y)*length_bar/2.0/sqrt(1.e-6*age*3.1536e13)));// 
	//	*val = Ttop+(Tbot-Ttop)*(1.0-erfc((Ly-y)/2.0/sqrt(diffusivity*age*3.1536e13)));// 
	
#if 0
	/*read age from file*/
	PetscReal agesin[100],a;
	PetscInt i,ii,i0;
	
	a=0.878132429842795;
	i0=0;
	for (i=0; i<100; i++) {
		agesin[i]=1.0/3.1536e7*pow(length_bar*(0.05+0.005*sin(6.0*3.1415926/1.25*(i*1.25/99)))/(2.0*a),2);
		if (x>=i*1.25/99.0&&x<=(i+1)*1.25/99.0){
			i0=i;
		}
	}
	*val = Ttop+(Tbot-Ttop)*(1.0-erfc((-y)*length_bar/2.0/sqrt(1.e-6*agesin[i0]*3.1536e13)));
#endif    
	/* indicate you want to set this value into the vector */
	impose = PETSC_TRUE;
	
	return impose;
}


#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_Interp_X_Subduction_Initiation2d"
PetscErrorCode DMDAVecTraverse3d_Interp_X_Subduction_Initiation2d(DM da,Vec X,PetscInt dof_idx,ModelSubduction_Initiation2dCtx *data)
{
#if 1
	PetscInt i,j,k,si,sj,sk,m,n,p,M,N,P,ndof;
	DM cda;
	Vec coords;
	DMDACoor3d ***LA_coords;	
	PetscScalar pos[3];
	PetscScalar val;
	PetscBool impose_value;
	PetscScalar ****LA_X;
	
	PetscScalar vx,Vx_max;
	PetscReal Dprime;
	
	PetscErrorCode ierr;
	
	PetscFunctionBegin;	
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = DMDAGetInfo(da,0, &M,&N,&P, 0,0,0, &ndof,0, 0,0,0, 0);CHKERRQ(ierr);
	if (dof_idx >= ndof) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"dof_index >= dm->blocksize"); }
	
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coords);CHKERRQ(ierr);
	if (!coords) { 
		PetscPrintf(PETSC_COMM_WORLD,"WARNING: coordinates not set\n"); 
	} else {
		ierr = DMDAVecGetArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	}
	ierr = DMDAVecGetArrayDOF(da,X,&LA_X);CHKERRQ(ierr);
	
	
	
	
	for (k=sk; k<sk+p; k++) {
		for (j=sj; j<sj+n; j++) {
			for (i=si; i<si+m; i++) {
				
				if (coords) {
					pos[0] = LA_coords[k][j][i].x;
					pos[1] = LA_coords[k][j][i].y;
					pos[2] = LA_coords[k][j][i].z;
				} else {
					pos[0] = pos[1] = pos[2] = 0.0;
				}
				
				Dprime = (data->Ly-data->Oy)/2.0;
				Vx_max = data->velocity;
				
				
				//                vx=-Vx_max*tanh(2.0*3.1415926*(pos[1]-data->Oy-Dprime)/(2.0*Dprime));
				vx=-Vx_max*(tanh(2.0*3.1415926*(pos[1]-data->Oy-Dprime)/(2.0*Dprime))+(1.0+pos[1]/Dprime))/2.0;
				
				
				LA_X[k][j][i][dof_idx] = 0.0+pos[0]/data->Lx*vx;
				
				
			}
		}
	}
	ierr = DMDAVecGetArrayDOF(da,X,&LA_X);CHKERRQ(ierr);
	if (coords) {
		ierr = DMDAVecRestoreArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
#endif
}




#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_Subduction_Initiation2d"
PetscErrorCode ModelApplyInitialSolution_Subduction_Initiation2d(pTatinCtx c,Vec X, void *ctx)
{
	
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	DM stokes_pack,dau,dap,daT;
	Vec velocity,pressure,temperature;
	PetscReal MeshMin[3],MeshMax[3],domain_height;
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	DMDAVecTraverse3d_InterpCtx IntpCtx;
	PetscReal vals[6];
	PetscScalar zero=0.0;
	PhysCompEnergy energy;
	PetscBool    active_energy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
#if 1	
	stokes_pack = c->stokes_ctx->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
#if 1	
	//	ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,-data->velocity/(data->Lx-data->Ox),0.0,0.0);CHKERRQ(ierr);
	//	ierr = DMDAVecTraverse3d(dau,velocity,0,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
	
	//	ierr = DMDAVecTraverse3d_Interp_X_Subduction_Initiation2d(dau,velocity,0,data);CHKERRQ(ierr);
	
	//	ierr = DMDAVecTraverse3d_InterpCtxSetUp_Y(&IntpCtx,0.0/(data->Ly-data->Oy),0.0,0.0);CHKERRQ(ierr);
	//	ierr = DMDAVecTraverse3d(dau,velocity,1,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
	
	//	ierr = DMDAVecTraverse3d_InterpCtxSetUp_Z(&IntpCtx,0.0/(data->Lz),0.0,0.0);CHKERRQ(ierr);
	//	ierr = DMDAVecTraverse3d(dau,velocity,2,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
#endif	
	
	ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
	
	ierr = DMDAGetBoundingBox(dau,MeshMin,MeshMax);CHKERRQ(ierr);
	domain_height = MeshMax[1] - MeshMin[1];
	
	HPctx.surface_pressure = 0.0;
	HPctx.ref_height = domain_height;
	HPctx.ref_N      = c->stokes_ctx->my-1;
	HPctx.grav       = 10.0;
	HPctx.rho        = data->rho[0];
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%e %d %e]]\n", HPctx.ref_height,HPctx.ref_N,HPctx.rho);
	
	ierr = DMDAVecTraverseIJK(dap,pressure,0,DMDAVecTraverseIJK_HydroStaticPressure_v2,     (void*)&HPctx);CHKERRQ(ierr); /* P = P0 + a.x + b.y + c.z, modify P0 (idx=0) */
	ierr = DMDAVecTraverseIJK(dap,pressure,2,DMDAVecTraverseIJK_HydroStaticPressure_dpdy_v2,(void*)&HPctx);CHKERRQ(ierr); /* P = P0 + a.x + b.y + c.z, modify b  (idx=2) */
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
#endif	
	
#if 1	
	
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	
	if (active_energy) {
		ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(c,&temperature,PETSC_NULL);CHKERRQ(ierr);
		daT  = energy->daT;
		
		vals[0]=data->Ttop; //top temperature K
		vals[1]=data->Tbot; //bottom temperature K
		vals[2]=data->Thermal_age;  //age Ma
		vals[3]=data->Ly;
		vals[4]=data->length_bar;
		vals[5]=data->diffusivity[0];
		
		ierr = DMDAVecTraverse3d(daT,temperature,0, DMDAVecTraverse3d_ERFC3DFunctionSubduction_Initiation2d, (void*)vals);CHKERRQ(ierr);
	}
#endif	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialStokesVariableMarkers_Subduction_Initiation2d"
PetscErrorCode ModelApplyInitialStokesVariableMarkers_Subduction_Initiation2d(pTatinCtx c,Vec X,void *ctx)
{
#if 1	
	DM                stokes_pack,dau,dap;
	PhysCompStokes    stokes;
	Vec               Uloc,Ploc;
	PetscScalar       *LA_Uloc,*LA_Ploc;
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	DataField                    PField;
	MaterialConst_MaterialType   *truc;
	PetscErrorCode    ierr;
	PetscInt regionidx;	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	ierr = pTatin_EvaluateRheologyNonlinearities(c,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	
	
#endif
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Subduction_Initiation2d"
PetscErrorCode ModelApplyUpdateMeshGeometry_Subduction_Initiation2d(pTatinCtx c,Vec X,void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	PetscReal        step;
	PhysCompStokes   stokes;
	Vec velocity,pressure;
	DM stokes_pack,dav,dap;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	//	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);
	ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Subduction_Initiation2d"
PetscErrorCode ModelOutput_Subduction_Initiation2d(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	PetscBool         active_energy;
	DataBucket        materialpoint_db;
	Vec               velocity,pressure;
	PetscBool         output_markers = PETSC_FALSE;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* output velocity/pressure */
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

	/*  Write out material points: standard, stokes, plastic variables, energy */
	if(output_markers){
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
		const int nf = 4;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
		char mp_file_prefix[256];
		
		sprintf(mp_file_prefix,"%s_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
	
	/*  Write out material point props on mesh: viscosity, density, inv2stress, plastic strain */
	/*		
	 {
	 const int                   nf = 4;
	 const MaterialPointVariable mp_prop_list[] = { MPV_viscosity, MPV_density, MPV_inv2stress, MPV_plastic_strain }; 
	 
	 ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
	 }	
	 */	
	
	/* output temperature */
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PhysCompEnergy energy;
		Vec            temperature;
		
		ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(c,&temperature,PETSC_NULL);CHKERRQ(ierr);
		
		ierr = pTatin3d_ModelOutput_Temperature_Energy(c,temperature,prefix);CHKERRQ(ierr);
		
		Quadrature     volQ;
		DM             daT;	
		
		daT  = energy->daT;
		volQ = energy->volQ;
	}
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Subduction_Initiation2d"
PetscErrorCode ModelDestroy_Subduction_Initiation2d(pTatinCtx c,void *ctx)
{
	ModelSubduction_Initiation2dCtx *data = (ModelSubduction_Initiation2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Subduction_Initiation2d"
PetscErrorCode pTatinModelRegister_Subduction_Initiation2d(void)
{
	ModelSubduction_Initiation2dCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelSubduction_Initiation2dCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelSubduction_Initiation2dCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"subduction_initiation2d");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,   (void (*)(void))ModelApplyInitialStokesVariableMarkers_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Subduction_Initiation2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Subduction_Initiation2d);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

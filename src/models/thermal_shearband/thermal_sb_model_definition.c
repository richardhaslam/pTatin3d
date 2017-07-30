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
 **    filename:   thermal_sb_model_definition.c
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


#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "data_bucket.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_iterator.h"
#include "mesh_deformation.h"
#include "mesh_update.h"
#include "dmda_remesh.h"
#include "output_material_points.h"
#include "material_constants.h"
#include "energy_output.h"

typedef struct {
    PetscReal L_bar,V_bar,E_bar,T_bar,t_bar,eta_bar,P_bar;
    PetscReal rhs_scale;
    PetscReal srate_xx;
    PetscReal srate_yy;
    PetscReal srate_zz;
    PetscReal inclusion_radius;
    PetscBool output_si;
} ThermalSBData;


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_ThermalSB"
PetscErrorCode ModelInitialize_ThermalSB(pTatinCtx ptatinctx,void *modelctx)
{
	ThermalSBData      *modeldata = (ThermalSBData*)modelctx;
	RheologyConstants  *rheology;
	DataBucket         materialconstants;
    PetscInt           regionidx;
    PetscReal          n_exp,F,eta_0,preexp_A;
    PetscInt           matrix_idx,inclusion_idx;
	PetscErrorCode     ierr;

    /* scaling params */
    modeldata->L_bar   = 35.0e3;  /* m */
    modeldata->E_bar   = 5.0e-14; /* strain rate  - 1/s */
    modeldata->eta_bar = 1.0e22; /* eta ~ 1e22 Pa s */
    modeldata->T_bar   = 673.0;   /* K */
    
    modeldata->V_bar = modeldata->E_bar * modeldata->L_bar; /* velocity */

    modeldata->P_bar = modeldata->eta_bar * modeldata->E_bar; /* stress */
    
    modeldata->t_bar = 1.0/modeldata->E_bar; /* time (sec) */
    
    modeldata->rhs_scale = modeldata->eta_bar * modeldata->V_bar / (modeldata->L_bar * modeldata->L_bar); /* [ eta.u/(L.L) ]^-1 */
    modeldata->rhs_scale = 1.0 / modeldata->rhs_scale;

    /* default strain rate for bcs */
    modeldata->srate_xx = -5.0e-14;
    modeldata->srate_yy = 0.0;
    modeldata->srate_zz = 0.0;

    
    /* scale strain rate */
    modeldata->srate_xx /= modeldata->E_bar;
    modeldata->srate_yy /= modeldata->E_bar;
    modeldata->srate_zz /= modeldata->E_bar;
    
    
    /* inclusion geometry */
    modeldata->inclusion_radius = 3.0e3;
    modeldata->inclusion_radius /= modeldata->L_bar;
    
    
    modeldata->output_si = PETSC_FALSE;
    PetscOptionsGetBool(NULL,NULL,"-model_thermal_sb_output_si",&modeldata->output_si,0);
    
    
    /* force energy equation to be introduced */
	ierr = PetscOptionsInsertString(NULL,"-activate_energy");CHKERRQ(ierr);

    /*
        Current Value -preexpA_0   :  3.2000e-20
        Current Value -Ascale_0    :  3.0285e-01
        Current Value -entalpy_0   :  2.7600e+05
        Current Value -Vmol_0      :  0.0000e+00
        Current Value -nexp_0      :  3.0000e+00
        Current Value -Tref_0      :  0.0000e+00
        Current Value -Eta_scale_0 :  1.0000e+22
        Current Value -P_scale_0   :  5.0000e+08
        Current Value -density_0   :  2.7000e+03

        Current Value -preexpA_1   :  3.1600e-26
        Current Value -Ascale_1    :  3.0154e-01
        Current Value -entalpy_1   :  1.9000e+05
        Current Value -Vmol_1      :  0.0000e+00
        Current Value -nexp_1      :  3.3000e+00
        Current Value -Tref_1      :  0.0000e+00
        Current Value -Eta_scale_1 :  1.0000e+22
        Current Value -P_scale_1   :  5.0000e+08
        Current Value -density_1   :  2.7000e+03  
    */
	/* Rheology prescription */
	ierr = pTatinGetRheology(ptatinctx,&rheology);CHKERRQ(ierr);
	rheology->rheology_type = RHEOLOGY_VP_STD;
    rheology->nphases_active = 2;
    
    /* Material constant */
	ierr = pTatinGetMaterialConstants(ptatinctx,&materialconstants);CHKERRQ(ierr);
    MaterialConstantsSetDefaults(materialconstants);

    matrix_idx    = 0;
    inclusion_idx = 1;
    
    /* matrix */
	MaterialConstantsSetValues_MaterialType(materialconstants,matrix_idx,VISCOUS_ARRHENIUS_2,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);

    n_exp = 3.0;
    F = pow( 2.0, (1.0-n_exp)/n_exp );
    F = F / pow( 3.0, (1.0+n_exp)/(2.0*n_exp) );
    eta_0 = 3.2e-20;
    preexp_A = eta_0;
    MaterialConstantsSetValues_DensityConst(materialconstants,matrix_idx,modeldata->rhs_scale * 2700.0);
    /* eta = F . pow(preexp_A,-1/n) . pow(e,1/n-1) . exp(E/nRT) */
    MaterialConstantsSetValues_ViscosityArrh(materialconstants,matrix_idx,preexp_A, F, 276.0e3, 0.0, n_exp, 0.0);
    MaterialConstantsScaleValues_ViscosityArrh(materialconstants,matrix_idx, modeldata->eta_bar, modeldata->P_bar );

    /* inclusion */
	MaterialConstantsSetValues_MaterialType(materialconstants,inclusion_idx,VISCOUS_ARRHENIUS_2,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);
    
    n_exp = 3.3;
    F = pow( 2.0, (1.0-n_exp)/n_exp );
    F = F / pow( 3.0, (1.0+n_exp)/(2.0*n_exp) );
    eta_0 = 3.16e-26;
    preexp_A = eta_0;
    MaterialConstantsSetValues_DensityConst(materialconstants,inclusion_idx,modeldata->rhs_scale * 2700.0);
    /* eta = F . pow(preexp_A,-1/n) . pow(e,1/n-1) . exp(E/nRT) */
    MaterialConstantsSetValues_ViscosityArrh(materialconstants,inclusion_idx,preexp_A, F, 190.0e3, 0.0, n_exp, 0.0);
    MaterialConstantsScaleValues_ViscosityArrh(materialconstants,inclusion_idx, modeldata->eta_bar, modeldata->P_bar );
    
    /* Material constant */
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		ierr= MaterialConstantsSetFromOptions(materialconstants,"model_thermal_sb",regionidx,PETSC_FALSE);CHKERRQ(ierr);
	}
    
    /* output */
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		MaterialConstantsPrintAll(materialconstants,regionidx);
	}
    
    /* scale material properties */
	//for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
    //    MaterialConstantsScaleAll(materialconstants,regionidx,data->length_bar,data->velocity_bar,data->time_bar,data->viscosity_bar,data->density_bar,data->pressure_bar);
    //}

    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_ThermalSB"
PetscErrorCode ModelApplyInitialMeshGeometry_ThermalSB(pTatinCtx ptatinctx,void *modelctx)
{
	ThermalSBData    *modeldata = (ThermalSBData*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    PetscReal        Lx,Ly,Lz;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* set initial velocity field */
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    
    /* dav -> velocity mesh */
    /* dap -> pressure mesh */
    /* sets cuboid mesh */
    Lx = 35.0 * 1.0e3; /* km */
    Ly = 30.0 * 1.0e3; /* km */
    Lz = 35.0 * 1.0e3; /* km */
	ierr = DMDASetUniformCoordinates(dav,0.0,Lx/modeldata->L_bar, 0.0,Ly/modeldata->L_bar, 0.0,Lz/modeldata->L_bar);CHKERRQ(ierr);
  {
    PetscReal gvec[] = { 0.0, -9.81, 0.0 };
    ierr = PhysCompStokesSetGravityVector(stokes,gvec);CHKERRQ(ierr);
  }
	
	PetscFunctionReturn(0);
}

/* velocity bcs */
#undef __FUNCT__
#define __FUNCT__ "ThermalSB_VelocityBC"
PetscErrorCode ThermalSB_VelocityBC(BCList bclist,DM dav,pTatinCtx ptatinctx,ThermalSBData *modelctx)
{
    ThermalSBData  *modeldata = (ThermalSBData*)modelctx;
    PetscReal      exx;
    PetscReal      Lx,Ly,vxE,vyN,gmin[3],gmax[3];
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
    exx = modeldata->srate_xx;
	//ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,exx,0);CHKERRQ(ierr);
    
    /*
     div(u) = exx + eyy + ezz
     
     2(vx)/Lx + 2(vy)/Ly + 0 = 0
     
     exx = 2.vx/Lx

     vE = exx.Lx/2
     vW = -vE
     vy = - 0.5 * exx * Ly
     
     vy = vx Ly/Lx on the top face
    */

    
    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    Lx = gmax[0] - gmin[0];
    Ly = gmax[1] - gmin[1];
    
    /* compute vx on the east face from exx */
    vxE = modeldata->srate_xx * Lx * 0.5;

	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&vxE);CHKERRQ(ierr);
    
    /* on the upper/north face */
    vyN = - 0.5 * exx * Ly;
    
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,NORTH_FACE,vyN);CHKERRQ(ierr);
    //ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,-vyN);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,0);CHKERRQ(ierr);

    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,FRONT_FACE,0);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,BACK_FACE,0);CHKERRQ(ierr);
    
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_ThermalSB"
PetscErrorCode ModelApplyBoundaryCondition_ThermalSB(pTatinCtx ptatinctx,void *modelctx)
{
	ThermalSBData         *modeldata = (ThermalSBData*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;
    
	
	PetscFunctionBegin;
    
	/* Define velocity boundary conditions */
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = ThermalSB_VelocityBC(stokes->u_bclist,dav,ptatinctx,modeldata);CHKERRQ(ierr);
    
	/* Define boundary conditions for any other physics */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_ThermalSB"
PetscErrorCode ModelApplyBoundaryConditionMG_ThermalSB(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx ptatinctx,void *modelctx)
{
	ThermalSBData    *modeldata = (ThermalSBData*)modelctx;
	PetscInt         n;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	
	/* Define velocity boundary conditions on each level within the MG hierarchy */
	for (n=0; n<nl; n++) {
		ierr = ThermalSB_VelocityBC(bclist[n],dav[n],ptatinctx,modeldata);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_ThermalSB"
PetscErrorCode ModelOutput_ThermalSB(pTatinCtx ptatinctx,Vec X,const char prefix[],void *modelctx)
{
	ThermalSBData    *modeldata = (ThermalSBData*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    Vec              coords,velocity,pressure;
    PetscBool        active_energy;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    if (modeldata->output_si) {
        /* get the velocity mesh */
        ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
        stokes_pack = stokes->stokes_pack;
        ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
        
        /* get the coordinates of the velocity mesh and scale into SI units <note, local and ghosted coordinates should be scaled> */
        ierr = DMGetCoordinates(dav,&coords);CHKERRQ(ierr);
        ierr = VecScale(coords,modeldata->L_bar);CHKERRQ(ierr);
        ierr = DMGetCoordinatesLocal(dav,&coords);CHKERRQ(ierr);
        ierr = VecScale(coords,modeldata->L_bar);CHKERRQ(ierr);
        
        /* unscale vel, p */
        ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
        ierr = VecScale(velocity,modeldata->V_bar);CHKERRQ(ierr);
        ierr = VecScale(pressure,modeldata->P_bar);CHKERRQ(ierr);
    }
    
	/* ---- Velocity-Pressure Mesh Output ---- */
	/* [1] Standard viewer: v,p written out as binary in double */
    ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);
	
    /* [2] Light weight viewer: Only v is written out. v and coords are expressed as floats */
    ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);

	/* [3] Write out v,p into PETSc Vec. These can be used to restart pTatin */
	/*
     ierr = pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
     */
    
    if (modeldata->output_si) {
        /* undo the coordinate scaling of velocity mesh <note, local and ghosted coordinates should be scaled> */
        ierr = DMGetCoordinates(dav,&coords);CHKERRQ(ierr);
        ierr = VecScale(coords,1.0/modeldata->L_bar);CHKERRQ(ierr);
        ierr = DMGetCoordinatesLocal(dav,&coords);CHKERRQ(ierr);
        ierr = VecScale(coords,1.0/modeldata->L_bar);CHKERRQ(ierr);

        /* unscale vel, p */
        ierr = VecScale(pressure,1.0/modeldata->P_bar);CHKERRQ(ierr);
        ierr = VecScale(velocity,1.0/modeldata->V_bar);CHKERRQ(ierr);
        ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    }
    
	/* ---- Material Point Output ---- */
	/* [1] Basic viewer: Only reports coords, regionid and other internal data */
	/*
     ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
     */
	
	/* [2] Customized viewer: User defines specific fields they want to view - NOTE not .pvd file will be created */
    {
        DataBucket                materialpoint_db;
        const int                 nf = 2;
        const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes };
        //const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
        char                      mp_file_prefix[256];
        
        ierr = pTatinGetMaterialPoints(ptatinctx,&materialpoint_db,NULL);CHKERRQ(ierr);
        sprintf(mp_file_prefix,"%s_mpoints",prefix);
        ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,ptatinctx->outputpath,mp_file_prefix);CHKERRQ(ierr);
    }

	/* [3] Customized marker->cell viewer: Marker data is projected onto the velocity mesh. User defines specific fields */
    {
        const int                    nf = 2;
        const MaterialPointVariable  mp_prop_list[] = { MPV_viscosity, MPV_density };
        
        ierr = pTatin3d_ModelOutput_MarkerCellFields(ptatinctx,nf,mp_prop_list,prefix);CHKERRQ(ierr);
    }	
	
    /* standard viewer for temperature */
	ierr = pTatinContextValid_Energy(ptatinctx,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PhysCompEnergy energy;
		Vec            temperature;
		
		ierr = pTatinGetContext_Energy(ptatinctx,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(ptatinctx,&temperature,NULL);CHKERRQ(ierr);
        
		ierr = pTatin3d_ModelOutput_Temperature_Energy(ptatinctx,temperature,prefix);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_ThermalSB"
PetscErrorCode ModelApplyInitialSolution_ThermalSB(pTatinCtx ptatinctx,Vec X,void *modelctx)
{
	ThermalSBData    *modeldata = (ThermalSBData*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    Vec              velocity,pressure;
    PetscReal        Lx,Ly,vxR,vxL,vyT,vyB,gmin[3],gmax[3];
	DMDAVecTraverse3d_InterpCtx IntpCtx;
    PetscBool        active_energy;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
    ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
    stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    
    /* interpolate vx and vy in the interior */
    /* velocity intial condition - background strain */
    ierr = VecZeroEntries(velocity);CHKERRQ(ierr);

    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    Lx = gmax[0] - gmin[0];
    Ly = gmax[1] - gmin[1];

    vxL = -modeldata->srate_xx * Lx * 0.5;
    vxR = -vxL;
    vxL = 0.0;
    
    /* on the upper/north face */
    vyT = -0.5 * modeldata->srate_xx * Ly;
    //vyB = -vyT;
    vyB = 0.0;
    
    ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vxR-vxL)/(Lx),vxL,0.0);CHKERRQ(ierr);
    ierr = DMDAVecTraverse3d(dav,velocity,0,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
    
    ierr = DMDAVecTraverse3d_InterpCtxSetUp_Y(&IntpCtx,(vyT-vyB)/(Ly),vyB,0.0);CHKERRQ(ierr);
    ierr = DMDAVecTraverse3d(dav,velocity,1,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
    
    /* set initial temperature field */
	/* initial condition for temperature */
	ierr = pTatinContextValid_Energy(ptatinctx,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PhysCompEnergy energy;
		Vec            temperature;
		/* DM             daT; */
		
		ierr = pTatinGetContext_Energy(ptatinctx,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(ptatinctx,&temperature,NULL);CHKERRQ(ierr);
		/* daT  = energy->daT; */
        
        ierr = VecSet(temperature,400.0+273.0);CHKERRQ(ierr);
	}
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_ThermalSB"
PetscErrorCode ModelApplyInitialMaterialGeometry_ThermalSB(pTatinCtx c,void *ctx)
{
	ThermalSBData    *data = (ThermalSBData*)ctx;
	MPAccess         mpX;
	int              p,n_mpoints;
	DataBucket       materialpoint_db;
	DataBucket       materialconstants;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    PetscReal        Ox[3],gmin[3],gmax[3],inc_rad2;
	PetscErrorCode   ierr;
	MaterialConst_DensityConst      *DensityConst_data;
	/* MaterialConst_ViscosityArrh     *ViscArrh_data; */
	DataField                       PField_DensityConst,PField_ViscArrh;
  PetscReal                       kappa;
  /* PetscReal                       H; */
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
	
    DataBucketGetDataFieldByName(materialconstants,MaterialConst_DensityConst_classname,&PField_DensityConst);
	DensityConst_data      = (MaterialConst_DensityConst*)PField_DensityConst->data;

	DataBucketGetDataFieldByName(materialconstants,MaterialConst_ViscosityArrh_classname,&PField_ViscArrh);
	/* ViscArrh_data          = (MaterialConst_ViscosityArrh*)PField_ViscArrh->data; */
	
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
    stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
    /* background */
    kappa = 1.0/(2700.0 * 1050.0) * 2.5 * data->t_bar / ( data->L_bar * data->L_bar );
    /* H = 1.0/(2700.0 * 1050.0) * data->eta_bar * data->E_bar * data->E_bar; */
    
    PetscPrintf(PETSC_COMM_WORLD,"thermal_sb: BACKGROUND kappa = %1.4e \n",kappa );
    for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointSet_phase_index(mpX,p, 0);CHKERRQ(ierr);
		ierr = MaterialPointSet_viscosity(mpX,  p, 1.0);CHKERRQ(ierr);
		ierr = MaterialPointSet_density(mpX,    p, DensityConst_data[0].density);CHKERRQ(ierr);

        ierr = MaterialPointSet_diffusivity(mpX,p,kappa);CHKERRQ(ierr);
        ierr = MaterialPointSet_heat_source(mpX,p,0.0);CHKERRQ(ierr);
    }

    inc_rad2 = data->inclusion_radius * data->inclusion_radius;

    /* inclusion */
    kappa = 1.0/(2700.0 * 1050.0 ) * 2.5 * data->t_bar / ( data->L_bar * data->L_bar );
    PetscPrintf(PETSC_COMM_WORLD,"thermal_sb: INCLUSION kappa = %1.4e \n",kappa );
	for (p=0; p<n_mpoints; p++) {
        double *position_p,r2;
        
		ierr = MaterialPointGet_global_coord(mpX,p,&position_p);CHKERRQ(ierr);
        
        //Ox[0] = (double)(gmin[0] + 0.5 * (gmax[0] - gmin[0]));
        //Ox[1] = (double)(gmin[1] + 0.5 * (gmax[1] - gmin[1]));
        //Ox[2] = (double)(gmin[2] + 0.5 * (gmax[2] - gmin[2]));
        
        Ox[0] = Ox[1] = Ox[2] = 0.0;
        
        r2 = (position_p[0]-Ox[0])*(position_p[0]-Ox[0]);
        r2 += (position_p[1]-Ox[1])*(position_p[1]-Ox[1]);
        r2 += (position_p[2]-Ox[2])*(position_p[2]-Ox[2]);
        
        if (r2 < inc_rad2) {
            ierr = MaterialPointSet_phase_index(mpX,p, 1);CHKERRQ(ierr);
            ierr = MaterialPointSet_viscosity(mpX,  p, 1.0);CHKERRQ(ierr);
            ierr = MaterialPointSet_density(mpX,    p, DensityConst_data[1].density);CHKERRQ(ierr);
            
            ierr = MaterialPointSet_diffusivity(mpX,p,kappa);CHKERRQ(ierr);
            ierr = MaterialPointSet_heat_source(mpX,p,0.0);CHKERRQ(ierr);
        }
	}
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_ThermalSB"
PetscErrorCode ModelApplyUpdateMeshGeometry_ThermalSB(pTatinCtx c,Vec X,void *ctx)
{
	PetscReal        step,gmin[3],gmax[3];
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	PetscErrorCode   ierr;
    
	/* fully lagrangian update */
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

    ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);
    
    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates1D(dav,0,gmin[0],gmax[0]);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates1D(dav,1,gmin[1],gmax[1]);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates1D(dav,2,gmin[2],gmax[2]);CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_ThermalSB"
PetscErrorCode ModelDestroy_ThermalSB(pTatinCtx c,void *ctx)
{
	ThermalSBData *data = (ThermalSBData*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* Free contents of structure */
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_ThermalSB"
PetscErrorCode pTatinModelRegister_ThermalSB(void)
{
	ThermalSBData   *data;
	pTatinModel     m;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ThermalSBData),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ThermalSBData));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
    
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"thermal_sb");CHKERRQ(ierr);
    
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_ThermalSB);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_ThermalSB);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

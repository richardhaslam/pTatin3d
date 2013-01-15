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
 **    Filename:      model_rift3D_ops.c
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
 **    $Id: model_rift3D_ops.c 3636 2013-01-14 10:36:00Z llepourhiet $
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"

#include "rift3D_ctx.h"




#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Rift3D"
PetscErrorCode ModelInitialize_Rift3D(pTatinCtx c,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
    RheologyConstants      *rheology;
	PetscBool flg;
	PetscReal max_eta;
    PetscScalar vx,vy,vz,Sx_bar,Sy_bar,Sz_bar;
	PetscInt n;
    PetscErrorCode ierr;
	
	PetscFunctionBegin;


	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
    rheology                = &c->rheology_constants;
	rheology->rheology_type = RHEOLOGY_VP_STD;
        
	data->dimensional   = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_rift3D_dimensional",&data->dimensional,&flg);CHKERRQ(ierr);

	/*
	data->density_bar   = 1000.0;
	data->length_bar    = 100.0 * 1.0e3;
	data->viscosity_bar = 1.0e22;
	data->velocity_bar  = 1.0e-10;
	data->time_bar      = data->length_bar / data->velocity_bar;
	data->pressure_bar  = data->length_bar * data->density_bar;
 */

	data->length_bar    = 100.0 * 1.0e3;
	data->viscosity_bar = 1e25;
	data->velocity_bar  = 1.0e-10;
	data->time_bar      = data->length_bar / data->velocity_bar;
	data->pressure_bar  = data->viscosity_bar*data->velocity_bar / data->length_bar;
	data->density_bar   = data->viscosity_bar*data->velocity_bar / ( data->length_bar * data->length_bar );
	
	
	/* box geometry */
	data->Lx =  6.0e5;
	data->Ly =  0.0e5;
	data->Lz =  6.0e5;
	data->Ox =  0.0e5;
	data->Oy =  -1.5e5;
	data->Oz =  0.0e5;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_Lx",&data->Lx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_Ly",&data->Ly,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_Lz",&data->Lz,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_Ox",&data->Ox,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_Oy",&data->Oy,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_Oz",&data->Oz,&flg);CHKERRQ(ierr);
    
    /* report */
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_Ox %1.4e [SI] -model_rift3D_Lx : %1.4e [SI]\n", data->Ox ,data->Lx );
    PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_Oy %1.4e [SI] -model_rift3D_Ly : %1.4e [SI]\n", data->Oy ,data->Ly );
    PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_Oz %1.4e [SI] -model_rift3D_Lz : %1.4e [SI]\n", data->Oz ,data->Lz );
    
    
    vx = 0.5;
    vz = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_vx",&vx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_vz",&vz,&flg);CHKERRQ(ierr);
    vx = vx*3.15e-10;
    vz = vz*3.15e-10;
    Sx_bar = (data->Ly - data->Oy)*(data->Lz - data->Oz);
    Sz_bar = (data->Ly - data->Oy)*(data->Lx - data->Ox);
    Sy_bar = (data->Lx - data->Ox)*(data->Lz - data->Oz);
    vy = (vx*Sx_bar+vz*Sz_bar)*2/Sy_bar;
    
    PetscPrintf(PETSC_COMM_WORLD,"  -model_rift3D_vx [m/s]:%e  -model_rift3D_vz [m/s]:%e computed vy [m/s]:%e \n", vx,vz,vy);

    /* parse from command line */
	rheology->nphases_active = 4;

	/* viscosity */
	rheology->const_eta0[0] = 2.0*1.0e23; /* crust */
	rheology->const_eta0[1] = 2.5*1.0e19; /* crust - lower */
	rheology->const_eta0[2] = 5.0*1.0e22; /* mantle */
	rheology->const_eta0[3] = 5.0*1.0e20; /* mantle - lower */

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_eta0",&rheology->const_eta0[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_eta1",&rheology->const_eta0[1],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_eta2",&rheology->const_eta0[2],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_eta3",&rheology->const_eta0[3],&flg);CHKERRQ(ierr);
	
    for (n=0; n<rheology->nphases_active; n++) {
		PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_eta%d [Pa.s]    : current value %1.4e [Pa.s]\n", n,rheology->const_eta0[n] );
	}

	/* density */
	rheology->const_rho0[0] = 2700.0; /* crust */
	rheology->const_rho0[1] = 2700.0; /* crust */
	rheology->const_rho0[2] = 3300.0; /* mantle */
	rheology->const_rho0[3] = 3280.0; /* mantle */

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_rho0",&rheology->const_rho0[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_rho1",&rheology->const_rho0[1],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_rho2",&rheology->const_rho0[2],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_rift3D_rho3",&rheology->const_rho0[3],&flg);CHKERRQ(ierr);
	
	for (n=0; n<rheology->nphases_active; n++) {
		PetscPrintf(PETSC_COMM_WORLD,"  input: -model_rift3D_rho%d [kg.m^-3] : current value %1.4e [kg.m^-3]\n", n,rheology->const_rho0[n] );
	}

	
	if (data->dimensional==PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"[rift3D]: using non-dimensional units\n");
		PetscPrintf(PETSC_COMM_WORLD,"  L*    : %1.4e [m]\n", data->length_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  U*    : %1.4e [m.s^-1]\n", data->velocity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  t*    : %1.4e [s]\n", data->time_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  eta*  : %1.4e [Pa.s]\n", data->viscosity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  rho*  : %1.4e [kg.m^-3]\n", data->density_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  P*    : %1.4e [Pa]\n", data->pressure_bar );

		//scale length 
        data->Lx = data->Lx / data->length_bar;
		data->Ly = data->Ly / data->length_bar;
		data->Lz = data->Lz / data->length_bar;
		data->Ox = data->Ox / data->length_bar;
		data->Oy = data->Oy / data->length_bar;
		data->Oz = data->Oz / data->length_bar;        
        PetscPrintf(PETSC_COMM_WORLD,"scaled values  -model_rift3D_Ox   :  %1.4e    -model_rift3D_Lx   :  %1.4e  \n", data->Ox ,data->Lx );
		PetscPrintf(PETSC_COMM_WORLD,"scaled values  -model_rift3D_Oy   :  %1.4e    -model_rift3D_Ly   :  %1.4e \n", data->Oy, data->Ly );
		PetscPrintf(PETSC_COMM_WORLD,"scaled values  -model_rift3D_Oz   :  %1.4e    -model_rift3D_Lz   :  %1.4e\n", data->Oz , data->Lz );
        //scale velocity
        data->vx = vx/data->velocity_bar;
        data->vy = vy/data->velocity_bar;
        data->vz = vz/data->velocity_bar;
        PetscPrintf(PETSC_COMM_WORLD,"scaled values  -model_rift3D_Vx:%1.4e    -model_rift3D_vy:%1.4e    -model_rift3D_vz:  %1.4e \n", data->vx ,data->vy, data->vz);
        
        // scale viscosity and density
		for (n=0; n<rheology->nphases_active; n++) {
			rheology->const_eta0[n] = rheology->const_eta0[n] / data->viscosity_bar;
			rheology->const_rho0[n] = rheology->const_rho0[n] / data->density_bar;
		}
		
		
		for (n=0; n<rheology->nphases_active; n++) {
			PetscPrintf(PETSC_COMM_WORLD,"  -model_rift3D_eta%d : scaled value %1.4e \n", n,rheology->const_eta0[n] );
		}
		
		for (n=0; n<rheology->nphases_active; n++) {
			PetscPrintf(PETSC_COMM_WORLD,"  -model_rift3D_rho%d : scaled value %1.4e \n", n,rheology->const_rho0[n] );
		}
		
	}

	/* set initial values for model parameters */
	/* material properties */
	data->nmaterials = rheology->nphases_active;
	data->eta[0] = rheology->const_eta0[0];
	data->eta[1] = rheology->const_eta0[1];
	data->eta[2] = rheology->const_eta0[2];
	data->eta[3] = rheology->const_eta0[3];
	
	data->rho[0] = rheology->const_rho0[0];
	data->rho[1] = rheology->const_rho0[1];
	data->rho[2] = rheology->const_rho0[2];
	data->rho[3] = rheology->const_rho0[3];
	
	
	
	PetscFunctionReturn(0);
}


/*
 
 ---------------
 |             |
 |             |
 |             |
 ---------------
 
 */

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Rift3D"
PetscErrorCode ModelApplyBoundaryCondition_Rift3D(pTatinCtx user,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscScalar vxl,vxr,vzf,vzb,vy;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    vxl  = - data->vx;
    vxr =  data->vx;
    vy = data->vy;
    vzf = - data->vz;
    vzb = data->vz;
        
	/* infilling free slip base */
	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&vy);CHKERRQ(ierr);
	/* free surface top*/
	
    
    /*extension along face of normal x */ 

        ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&(vxl));CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&(vxr));CHKERRQ(ierr);
    
    /*compression along face of normal z */ 

	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&(vzb));CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&(vzf));CHKERRQ(ierr);
    
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Rift3D"
PetscErrorCode ModelApplyBoundaryConditionMG_Rift3D(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscScalar vxl,vxr,vzf,vzb,vy;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    
    vxl  = - data->vx;
    vxr =  data->vx;
    vy = data->vy;
    vzf = - data->vz;
    vzb = data->vz;
 
	for (n=0; n<nl; n++) {
		
        
        /* infilling free slip base */
        ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&vy);CHKERRQ(ierr);
        /* free surface top*/
        
        
        /*extension along face of normal x */ 
        
        ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&(vxl));CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&(vxr));CHKERRQ(ierr);
        
        /*compression along face of normal z */ 
        
        ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&(vzb));CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&(vzf));CHKERRQ(ierr);
    }	
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Rift3D"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Rift3D(pTatinCtx c,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscErrorCode ierr;
	
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
	PetscInt               e,p,n_mp_points;
    PetscScalar            y_lab,y_moho,y_midcrust;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	int                    phase;
	PetscErrorCode ierr;
	
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
	
	
    y_lab      = -120.0e3; 
    y_moho     = -35.0e3;
    y_midcrust = -20.0e3;
    
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		double      *position,ycoord;
		double      eta,rho;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		if (data->dimensional == PETSC_FALSE) {
			ycoord = position[1] * data->length_bar;
		} else {
			ycoord = position[1];
		}
		
		
		if (ycoord<y_lab) {
			phase = 3;
			eta =  data->eta[3];
			rho =  data->rho[3];
		} else if (ycoord<y_moho) {
			phase = 2;
			eta =  data->eta[2];
			rho =  data->rho[2];
		} else if (ycoord<y_midcrust) {
			phase = 1;
			eta =  data->eta[1];
			rho =  data->rho[1];
		} else {
			phase = 0;
			eta =  data->eta[0];
			rho =  data->rho[0];
		}
		rho = -rho * 10.0;

		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,phase);
		
		MPntPStokesSetField_eta_effective(mpprop_stokes,eta);
		MPntPStokesSetField_density(mpprop_stokes,rho);
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Rift3D"
PetscErrorCode ModelApplyUpdateMeshGeometry_Rift3D(pTatinCtx c,Vec X,void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscErrorCode ierr;
	
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
#define __FUNCT__ "ModelOutput_Rift3D"
PetscErrorCode ModelOutput_Rift3D(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelRift3DCtx *data = (ModelRift3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = ModelOutput_Rift3D_CheckScales(c,X);CHKERRQ(ierr);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	
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
#define __FUNCT__ "pTatinModelRegister_Rift3D"
PetscErrorCode pTatinModelRegister_Rift3D(void)
{
	ModelRift3DCtx *data;
	pTatinModel m,model;
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

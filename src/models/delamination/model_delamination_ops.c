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
 **    filename:   model_delamination_ops.c
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
#include "dmda_iterator.h"
#include "mesh_update.h"

#include "delamination_ctx.h"

PetscReal cm_per_year2m_per_sec = 1.0e-2 / ( 365.0 * 24.0 * 60.0 * 60.0 ) ;

PetscErrorCode FormFunction_Stokes(SNES snes,Vec X,Vec F,void *ctx);


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Delamination"
PetscErrorCode ModelInitialize_Delamination(pTatinCtx c,void *ctx)
{
	ModelDelaminationCtx *data = (ModelDelaminationCtx*)ctx;
	RheologyConstants      *rheology;
	PetscBool flg;
	PetscScalar vx,vy,vz,Sx_bar,Sy_bar,Sz_bar;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	rheology                = &c->rheology_constants;
	rheology->rheology_type = RHEOLOGY_VISCOUS;
	//	rheology->rheology_type = RHEOLOGY_VP_STD;
	
	data->dimensional   = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,NULL,"-model_delamination_dimensional",&data->dimensional,&flg);CHKERRQ(ierr);
	if (data->dimensional) {
		PetscPrintf(PETSC_COMM_WORLD,"Delamination model expects the following dimensions for input\n");
		PetscPrintf(PETSC_COMM_WORLD," Box geometry: [m] \n");
		PetscPrintf(PETSC_COMM_WORLD," Viscosity:    [Pa.s] \n");
		PetscPrintf(PETSC_COMM_WORLD," Velocity:     [cm/yr] \n");
		PetscPrintf(PETSC_COMM_WORLD," Density:      [kg/m^3] \n");
	}
	
	/*
	 data->density_bar   = 1000.0;
	 data->length_bar    = 100.0 * 1.0e3;
	 data->viscosity_bar = 1.0e22;
	 data->velocity_bar  = 1.0e-10;
	 data->time_bar      = data->length_bar / data->velocity_bar;
	 data->pressure_bar  = data->length_bar * data->density_bar;
	 */
	
	data->length_bar    = 100.0 * 1.0e3;
	data->viscosity_bar = 1e23;
	data->velocity_bar  = 1.0e-10;
	
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_vis_bar",&data->viscosity_bar,&flg);CHKERRQ(ierr);
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_vel_bar",&data->velocity_bar,&flg);CHKERRQ(ierr);
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_length_bar",&data->length_bar,&flg);CHKERRQ(ierr);
	data->time_bar      = data->length_bar / data->velocity_bar;
	data->pressure_bar  = data->viscosity_bar*data->velocity_bar / data->length_bar;
	data->density_bar   = data->viscosity_bar*data->velocity_bar / ( data->length_bar * data->length_bar );
	/*
	 if (!data->dimensional) {
	 data->length_bar    = 1.0;
	 data->viscosity_bar = 1.0;
	 data->velocity_bar  = 1.0;
	 data->time_bar      = 1.0;
	 data->pressure_bar  = 1.0;
	 data->density_bar   = 1.0;
	 }
	 */
	
	/* box geometry, m */
	data->Lx =  8.0e5;
	data->Ly =  0.0e5;
	data->Lz =  8.0e5;
	data->Ox =  0.0e5;
	data->Oy = -4.0e5;
	data->Oz =  0.0e5;
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_Lx",&data->Lx,&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { data->Lx *= data->length_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_Ly",&data->Ly,&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { data->Ly *= data->length_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_Lz",&data->Lz,&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { data->Lz *= data->length_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_Ox",&data->Ox,&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { data->Ox *= data->length_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_Oy",&data->Oy,&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { data->Oy *= data->length_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_Oz",&data->Oz,&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { data->Oz *= data->length_bar; }
	
	/* report */
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_delamination_Ox %+1.4e [SI] -model_delamiantion_Lx : %+1.4e [SI]\n", data->Ox ,data->Lx );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_delamination_Oy %+1.4e [SI] -model_delamination_Ly : %+1.4e [SI]\n", data->Oy ,data->Ly );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_delamination_Oz %+1.4e [SI] -model_delamination_Lz : %+1.4e [SI]\n", data->Oz ,data->Lz );
	
	/* velocity cm/y */
	vx = 0.0;
	vz = 0.0;
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_vx",&vx,&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { vx *= data->velocity_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_vz",&vz,&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { vz *= data->velocity_bar; }
	/* convert to m/s */
	vx = vx * cm_per_year2m_per_sec;
	vz = vz * cm_per_year2m_per_sec;
	Sx_bar = (data->Ly - data->Oy)*(data->Lz - data->Oz);
	Sz_bar = (data->Ly - data->Oy)*(data->Lx - data->Ox);
	Sy_bar = (data->Lx - data->Ox)*(data->Lz - data->Oz);
	vy = (vx*Sx_bar+vz*Sz_bar)*2/Sy_bar;
	
	PetscPrintf(PETSC_COMM_WORLD,"  -model_delamination_vx [m/s]:  %+1.4e  -model_delamination_vz [m/s]:  %+1.4e : computed vy [m/s]:  %+1.4e \n", vx,vz,vy);
	PetscPrintf(PETSC_COMM_WORLD,"  -model_delamination_vx [cm/yr]:%+1.4e  -model_delamination_vz [cm/yr]:%+1.4e : computed vy [cm/yr]:%+1.4e \n", vx/cm_per_year2m_per_sec,vz/cm_per_year2m_per_sec,vy/cm_per_year2m_per_sec);
	
	/* parse from command line */
	rheology->nphases_active = 5;
	
	/* viscosity */
	rheology->const_eta0[0] = 2.0*1.0e23; /* crust */
	rheology->const_eta0[1] = 1.0*1.0e19; /* crust - lower */
	rheology->const_eta0[2] = 2.0*1.0e22; /* mantle - litho */
	rheology->const_eta0[3] = 1.0*1.0e19; /* mantle - astheno */
	rheology->const_eta0[4] = 1.0*1.0e21; /* mantle - lower */
	
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_eta0",&rheology->const_eta0[0],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_eta0[0] *= data->viscosity_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_eta1",&rheology->const_eta0[1],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_eta0[1] *= data->viscosity_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_eta2",&rheology->const_eta0[2],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_eta0[2] *= data->viscosity_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_eta3",&rheology->const_eta0[3],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_eta0[3] *= data->viscosity_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_eta4",&rheology->const_eta0[4],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_eta0[4] *= data->viscosity_bar; }
	
	
	for (n=0; n<rheology->nphases_active; n++) {
		PetscPrintf(PETSC_COMM_WORLD,"  input: -model_delamination_eta%d [Pa.s]    : current value %+1.4e [Pa.s]\n", n,rheology->const_eta0[n] );
	}
	
	/* density */
	data->rho0=3200.0;
	rheology->const_rho0[0] = 2700.0; /* crust */
	rheology->const_rho0[1] = 2700.0; /* crust */
	rheology->const_rho0[2] = 3300.0; /* mantle */
	rheology->const_rho0[3] = 3280.0; /* mantle */
	rheology->const_rho0[4] = 3300.0; /* mantle */
	
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_rho0",&rheology->const_rho0[0],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_rho0[0] *= data->density_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_rho1",&rheology->const_rho0[1],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_rho0[1] *= data->density_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_rho2",&rheology->const_rho0[2],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_rho0[2] *= data->density_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_rho3",&rheology->const_rho0[3],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_rho0[3] *= data->density_bar; }
	flg = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,"-model_delamination_rho4",&rheology->const_rho0[4],&flg);CHKERRQ(ierr); if (flg && !data->dimensional) { rheology->const_rho0[4] *= data->density_bar; }
	
	for (n=0; n<rheology->nphases_active; n++) {
		PetscPrintf(PETSC_COMM_WORLD,"  input: -model_delamination_rho%d [kg.m^-3] : current value %+1.4e [kg.m^-3]\n", n,rheology->const_rho0[n] );
	}
	
	PetscPrintf(PETSC_COMM_WORLD,"[delamination]: using non-dimensional units\n");
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
	//scale velocity
	data->vx = vx/data->velocity_bar;
	data->vy = vy/data->velocity_bar;
	data->vz = vz/data->velocity_bar;
	// scale viscosity and density
	for (n=0; n<rheology->nphases_active; n++) {
		rheology->const_eta0[n] = rheology->const_eta0[n] / data->viscosity_bar;
		rheology->const_rho0[n] = rheology->const_rho0[n] / data->density_bar;
	}
	data->rho0 = data->rho0/ data->density_bar;
	
	PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_delamination_Ox   :  %+1.4e    -model_delamination_Lx   :  %+1.4e  \n", data->Ox ,data->Lx );
	PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_delamination_Oy   :  %+1.4e    -model_delamination_Ly   :  %+1.4e \n", data->Oy, data->Ly );
	PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_delamination_Oz   :  %+1.4e    -model_delamination_Lz   :  %+1.4e\n", data->Oz , data->Lz );
	
	PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_delamination_Vx:%+1.4e    -model_delamination_vy:%+1.4e    -model_delamination_vz:  %+1.4e \n", data->vx ,data->vy, data->vz);
	for (n=0; n<rheology->nphases_active; n++) {
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_delamination_eta%d : %1.4e \n", n,rheology->const_eta0[n] );
	}
	for (n=0; n<rheology->nphases_active; n++) {
		PetscPrintf(PETSC_COMM_WORLD,"scaled value   -model_delamination_rho%d : scaled value %+1.4e \n", n,rheology->const_rho0[n] );
	}
	
	/* set initial values for model parameters */
	/* material properties */
	data->nmaterials = rheology->nphases_active;
	data->eta[0] = rheology->const_eta0[0];
	data->eta[1] = rheology->const_eta0[1];
	data->eta[2] = rheology->const_eta0[2];
	data->eta[3] = rheology->const_eta0[3];
	data->eta[4] = rheology->const_eta0[4];
	
	data->rho[0] = rheology->const_rho0[0];
	data->rho[1] = rheology->const_rho0[1];
	data->rho[2] = rheology->const_rho0[2];
	data->rho[3] = rheology->const_rho0[3];
	data->rho[4] = rheology->const_rho0[4];
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDelamination_DefineBCList"
PetscErrorCode ModelDelamination_DefineBCList(BCList bclist,DM dav,pTatinCtx user,ModelDelaminationCtx *data)
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
#define __FUNCT__ "ModelApplyBoundaryCondition_Delamination"
PetscErrorCode ModelApplyBoundaryCondition_Delamination(pTatinCtx user,void *ctx)
{
	ModelDelaminationCtx *data = (ModelDelaminationCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = ModelDelamination_DefineBCList(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,user,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Delamination"
PetscErrorCode ModelApplyBoundaryConditionMG_Delamination(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelDelaminationCtx *data = (ModelDelaminationCtx*)ctx;
	PetscInt       n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = ModelDelamination_DefineBCList(bclist[n],dav[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Delamination"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Delamination(pTatinCtx c,void *ctx)
{
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Delamination"
PetscErrorCode ModelApplyInitialMeshGeometry_Delamination(pTatinCtx c,void *ctx)
{
	ModelDelaminationCtx *data = (ModelDelaminationCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,data->Ox,data->Lx,data->Oy,data->Ly,data->Oz,data->Lz);CHKERRQ(ierr);
  {
    PetscReal gvec[] = { 0.0, -10.0, 0.0 };
    ierr = PhysCompStokesSetGravityVector(c->stokes_ctx,gvec);CHKERRQ(ierr);
  }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Delamination"
PetscErrorCode ModelApplyInitialMaterialGeometry_Delamination(pTatinCtx c,void *ctx)
{
	ModelDelaminationCtx *data = (ModelDelaminationCtx*)ctx;
	int                   p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	int                    phase;
	PetscScalar         y_midmantle = -300.0e3;
	PetscScalar     	y_lab2      = -150.0e3;
	PetscScalar         y_lab1      = -120.0e3; 
	PetscScalar     	y_moho      = -40.0e3;
	PetscScalar     	y_midcrust  = -20.0e3;
	
	PetscScalar         zG          = 500.0e3; 
	PetscScalar         zY          = 200.0e3;
	PetscScalar         xSN         = 200.0e3;
	PetscScalar         xWL         = 400.0e3;
	PetscScalar         xGV         = 300.0e3;
	PetscScalar         xDV         = 550.0e3;
	
	
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

	/* m */
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		double      *position,ycoord,xcoord,zcoord;
		double      eta,rho;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		/* convert to scaled units */
		xcoord = position[0] * data->length_bar;
		ycoord = position[1] * data->length_bar;
		zcoord = position[2] * data->length_bar;
		if (ycoord < y_midmantle) {
			phase = 4;
			eta =  data->eta[4];
			rho =  data->rho[4];
			
		} else if (ycoord < y_lab2) {
			phase = 3;
			eta =  data->eta[3];
			rho =  data->rho[3];
			
		} else if (ycoord < y_lab1) {
			
			if ((zcoord <= zG) && (zcoord >= zY) && (xcoord > xGV) && (xcoord <= xWL)){
				phase = 2;
				eta =  data->eta[2];
				rho =  data->rho[2];
			} else {
				phase = 3;
				eta =  data->eta[3];
				rho =  data->rho[3];    
			}    
			
		} else if (ycoord < y_moho) {
			
			if (zcoord < zG ){
				if ((xcoord > xWL) && (xcoord <= xDV) && (zcoord > zY)){
					phase = 3;
					eta =  data->eta[3];
					rho =  data->rho[3];    
				} else {
					phase = 2;
					eta =  data->eta[2];
					rho =  data->rho[2];
				}
			} else {
				
				if ((xcoord > xSN) && (xcoord <= xDV)){
					phase = 3;
					eta =  data->eta[3];
					rho =  data->rho[3];
				} else {
					phase = 2;
					eta =  data->eta[2];
					rho =  data->rho[2];
				}                
			}
			
		} else if (ycoord < y_midcrust) {
			if ((xcoord > xSN) && (xcoord <= xDV)) {
				phase = 1;
				eta =  data->eta[1];
				rho =  data->rho[1];
			} else {
				phase = 0;
				eta =  data->eta[0];
				rho =  data->rho[0];
			}
		} else {
			phase = 0;
			eta =  data->eta[0];
			rho =  data->rho[0];
		}
		
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
#define __FUNCT__ "ModelApplyInitialCondition_Delamination"
PetscErrorCode ModelApplyInitialCondition_Delamination(pTatinCtx c,Vec X,void *ctx)
{
	ModelDelaminationCtx *data = (ModelDelaminationCtx*)ctx;
	DM stokes_pack,dau,dap;
	Vec velocity,pressure;
	PetscReal vxl,vxr,vzb,vzf,vy;
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
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
	vzb =  0.0;//data->vz;
	
	//	ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
	/* apply -5 < vx 5 across the domain x \in [0,1] */
	
	/*	ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vxr-vxl)/(data->Lx-data->Ox),vxl,0.0);CHKERRQ(ierr);
	 ierr = DMDAVecTraverse3d(dau,velocity,0,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
	 ierr = DMDAVecTraverse3d_InterpCtxSetUp_Z(&IntpCtx,(vzf-vzb)/(data->Lz-data->Oz),vzb,0.0);CHKERRQ(ierr);
	 ierr = DMDAVecTraverse3d(dau,velocity,2,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
	 
	 ierr = DMDAVecTraverse3d_InterpCtxSetUp_Y(&IntpCtx,-vy/(data->Ly-data->Oy),0.0,0.0);CHKERRQ(ierr);
	 ierr = DMDAVecTraverse3d(dau,velocity,1,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
	 
	 */	
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
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Delamination"
PetscErrorCode ModelApplyUpdateMeshGeometry_Delamination(pTatinCtx c,Vec X,void *ctx)
{
	PetscReal        step;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	
	PetscErrorCode ierr;
	
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
	
	//ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	//ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,N);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Delamination_CheckScales"
PetscErrorCode ModelOutput_Delamination_CheckScales(pTatinCtx c,Vec X)
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
	
	PetscPrintf(PETSC_COMM_WORLD,"[Delamination]: check scales \n");
	
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
#define __FUNCT__ "ModelOutput_Delamination"
PetscErrorCode ModelOutput_Delamination(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = ModelOutput_Delamination_CheckScales(c,X);CHKERRQ(ierr);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Delamination"
PetscErrorCode ModelDestroy_Delamination(pTatinCtx c,void *ctx)
{
	ModelDelaminationCtx *data = (ModelDelaminationCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Delamination"
PetscErrorCode pTatinModelRegister_Delamination(void)
{
	ModelDelaminationCtx *data;
	pTatinModel m;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelDelaminationCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelDelaminationCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"delamination");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Delamination);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialCondition_Delamination);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

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
 **    Filename:      model_indentor_ops.c
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

#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"

#include "indentor_ctx.h"




#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Indentor"
PetscErrorCode ModelInitialize_Indentor(pTatinCtx c,void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
  RheologyConstants      *rheology;
	PetscBool flg;
	PetscReal max_eta;
	PetscInt n;
	PetscReal km2m,Ma2sec,cm_per_yer2m_per_sec;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;


	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
  rheology                = &c->rheology_constants;
	rheology->rheology_type = RHEOLOGY_VISCOUS;

	/* bc type */
	data->boundary_conditon_type = VSBC_FreeSlip;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-model_indentor_bc_type",(PetscInt*)&data->boundary_conditon_type,&flg);CHKERRQ(ierr);
	
	data->dimensional   = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_indentor_dimensional",&data->dimensional,&flg);CHKERRQ(ierr);

	/* default scales for non-dimensionalisation */
	km2m                 = 1.0e3;
	Ma2sec               = 1.0e6 * (365.0 * 24.0 * 60.0 * 60.0 );
	cm_per_yer2m_per_sec = 1.0e-2 / ( 365.0 * 24.0 * 60.0 * 60.0 ) ;
	
	/*
	data->density_bar   = 1000.0;
	data->length_bar    = 100.0 * 1.0e3;
	data->viscosity_bar = 1.0e22;
	data->velocity_bar  = 1.0e-10;
	data->time_bar      = data->length_bar / data->velocity_bar;
	data->pressure_bar  = data->length_bar * data->density_bar;
 */

	data->length_bar    = 120.0 * 1.0e3;
	data->viscosity_bar = 2.5e21;
	data->velocity_bar  = 1.0e-10;
	data->time_bar      = data->length_bar / data->velocity_bar;
	data->pressure_bar  = data->viscosity_bar*data->velocity_bar / data->length_bar;
	data->density_bar   = data->viscosity_bar*data->velocity_bar / ( data->length_bar * data->length_bar );
	
	
	/* box geometry */
	data->Lx = 500.0;
	data->Ly = 120.0;
	data->Lz = 500.0;

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_Lx",&data->Lx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_Ly",&data->Ly,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_Lz",&data->Lz,&flg);CHKERRQ(ierr);

	/* convert input to meters */
	data->Lx = data->Lx * 1.0e3;
	data->Ly = data->Ly * 1.0e3;
	data->Lz = data->Lz * 1.0e3;
	
	
	/* parse from command line */
	rheology->nphases_active = 4;

	/* viscosity */
	rheology->const_eta0[0] = 2.0*1.0e23; /* crust */
	rheology->const_eta0[1] = 2.5*1.0e21; /* crust - lower */
	rheology->const_eta0[2] = 5.0*1.0e23; /* mantle */
	rheology->const_eta0[3] = 5.0*1.0e21; /* mantle - lower */

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_eta0",&rheology->const_eta0[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_eta1",&rheology->const_eta0[1],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_eta2",&rheology->const_eta0[2],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_eta3",&rheology->const_eta0[3],&flg);CHKERRQ(ierr);
	
	/* density */
	rheology->const_rho0[0] = 2700.0; /* crust */
	rheology->const_rho0[1] = 2700.0; /* crust */
	rheology->const_rho0[2] = 3200.0; /* mantle */
	rheology->const_rho0[3] = 3200.0; /* mantle */

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_rho0",&rheology->const_rho0[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_rho1",&rheology->const_rho0[1],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_rho2",&rheology->const_rho0[2],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_rho3",&rheology->const_rho0[3],&flg);CHKERRQ(ierr);
	
	data->cutoff_time = 1.0;
	data->indentation_velocity = 1.5;
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_cutofftime",&data->cutoff_time,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_indentation_velocity",&data->indentation_velocity,&flg);CHKERRQ(ierr);

	/* convert input time Ma => sec */
	data->cutoff_time = data->cutoff_time * Ma2sec;
	
	/* convert input velocity cm/yr => m/s */
	data->indentation_velocity = data->indentation_velocity * cm_per_yer2m_per_sec;
	
	
	/* report */
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_indentor_Lx [km]        : current value %1.4e [km]\n", data->Lx*(1.0/km2m) );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_indentor_Ly [km]        : current value %1.4e [km]\n", data->Ly*(1.0/km2m) );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_indentor_Lz [km]        : current value %1.4e [km]\n", data->Lz*(1.0/km2m) );
  
	
	for (n=0; n<rheology->nphases_active; n++) {
		PetscPrintf(PETSC_COMM_WORLD,"  input: -model_indentor_eta%d [Pa.s]    : current value %1.4e [Pa.s]\n", n,rheology->const_eta0[n] );
	}

	for (n=0; n<rheology->nphases_active; n++) {
		PetscPrintf(PETSC_COMM_WORLD,"  input: -model_indentor_rho%d [kg.m^-3] : current value %1.4e [kg.m^-3]\n", n,rheology->const_rho0[n] );
	}

	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_indentor_cutofftime [Ma]                 : current value %1.4e [Ma]\n", data->cutoff_time*(1.0/Ma2sec) );
	PetscPrintf(PETSC_COMM_WORLD,"  input: -model_indentor_indentation_velocity [cm.yr^-1] : current value %1.4e [cm.yr^-1]\n", data->indentation_velocity*(1.0/cm_per_yer2m_per_sec) );
	
	if (data->dimensional==PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD,"[indentor]: using non-dimensional units\n");

		PetscPrintf(PETSC_COMM_WORLD,"  L*    : %1.4e [m]\n", data->length_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  U*    : %1.4e [m.s^-1]\n", data->velocity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  t*    : %1.4e [s]\n", data->time_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  eta*  : %1.4e [Pa.s]\n", data->viscosity_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  rho*  : %1.4e [kg.m^-3]\n", data->density_bar );
		PetscPrintf(PETSC_COMM_WORLD,"  P*    : %1.4e [Pa]\n", data->pressure_bar );

		PetscPrintf(PETSC_COMM_WORLD,"  Paraview: scaling factor for velocity [mm/yr] = %1.4e \n", 10.0 * data->velocity_bar/cm_per_yer2m_per_sec );
		PetscPrintf(PETSC_COMM_WORLD,"  Paraview: scaling factor for pressure [Mpa]   = %1.4e \n", 1.0e-6 * data->pressure_bar );
		
		data->Lx = data->Lx / data->length_bar;
		data->Ly = data->Ly / data->length_bar;
		data->Lz = data->Lz / data->length_bar;
		
		for (n=0; n<rheology->nphases_active; n++) {
			rheology->const_eta0[n] = rheology->const_eta0[n] / data->viscosity_bar;
			rheology->const_rho0[n] = rheology->const_rho0[n] / data->density_bar;
		}
		
		data->indentation_velocity = data->indentation_velocity / data->velocity_bar;
		data->cutoff_time = data->cutoff_time / data->time_bar;
		
		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_Lx   : scaled value %1.4e \n", data->Lx );
		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_Ly   : scaled value %1.4e \n", data->Ly );
		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_Lz   : scaled value %1.4e \n", data->Lz );

		for (n=0; n<rheology->nphases_active; n++) {
			PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_eta%d : scaled value %1.4e \n", n,rheology->const_eta0[n] );
		}
		
		for (n=0; n<rheology->nphases_active; n++) {
			PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_rho%d : scaled value %1.4e \n", n,rheology->const_rho0[n] );
		}
		
		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_cutofftime           : scaled value %1.4e \n", data->cutoff_time );
		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_indentation_velocity : scaled value %1.4e \n", data->indentation_velocity );
		
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"[indentor]: using dimensional units\n");

		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_Lx [m]         : current value %1.4e [m]\n", data->Lx );
		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_Ly [m]         : current value %1.4e [m]\n", data->Ly );
		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_Lz [m]         : current value %1.4e [m]\n", data->Lz );

		for (n=0; n<rheology->nphases_active; n++) {
			PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_eta%d [Pa.s]    : scaled value %1.4e [Pa.s]\n", n,rheology->const_eta0[n] );
		}
		
		for (n=0; n<rheology->nphases_active; n++) {
			PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_rho%d [kg.m^-3] : scaled value %1.4e [kg.m^-3]\n", n,rheology->const_rho0[n] );
		}

		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_cutofftime [s]               : scaled value %1.4e [s] \n", data->cutoff_time );
		PetscPrintf(PETSC_COMM_WORLD,"  -model_indentor_indentation_velocity [m.s^1] : scaled value %1.4e [m.s^1]\n", data->indentation_velocity );
		
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
#define __FUNCT__ "BCListEvaluator_indentor"
PetscBool BCListEvaluator_indentor( PetscScalar position[], PetscScalar *value, void *ctx ) 
{
	PetscBool impose_dirichlet = PETSC_TRUE;
	pTatinCtx user = (pTatinCtx)ctx;
	PetscScalar vx,Vx_max;
	PetscReal time,Dprime,cc;
	ModelIndentorCtx *model_data_ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = pTatinModelGetUserData(user->model,(void**)&model_data_ctx);CHKERRQ(ierr);
	time = user->time;

	Dprime = model_data_ctx->Lz/4.0;
	Vx_max = model_data_ctx->indentation_velocity;
	
	if ( (position[2] >= 3.0*Dprime) && (position[2] <= 4.0*Dprime) ) {
		vx = -Vx_max;
	} else if ( (position[2] >= 2.0*Dprime) && (position[2] <= 3.0*Dprime) ) {
		cc = cos( (M_PI/(2.0*Dprime) )*(position[2] - Dprime ) );
		vx = -Vx_max * cc * cc;
	} else {
		vx = 0.0;
	}
	*value = vx;
	return impose_dirichlet;
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Indentor"
PetscErrorCode ModelApplyBoundaryCondition_Indentor(pTatinCtx user,void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscScalar zero = 0.0;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* free slip base */
	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	/* free surface top */

	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

	/* indentor face */
	/* x is prescribed via function */
	/* y is unconstrained */
	/* z is prescribed as zero */
	// X
	//ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_indentor,(void*)user);CHKERRQ(ierr);
	// Y
	//ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	// Z
	ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Indentor"
PetscErrorCode ModelApplyBoundaryConditionMG_Indentor(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscScalar zero = 0.0;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	for (n=0; n<nl; n++) {
		
		/* free slip base */
		ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		/* free surface top */
		
		ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		/* indentor face */
		/* x is prescribed via function */
		/* y is unconstrained */
		/* z is prescribed as zero */
		// X
		//ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_indentor,(void*)user);CHKERRQ(ierr);
		// Y
		//ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		// Z
		ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Indentor"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Indentor(pTatinCtx c,void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Indentor"
PetscErrorCode ModelApplyInitialMeshGeometry_Indentor(pTatinCtx c,void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,0.0,data->Lx, 0.0,data->Ly, 0.0,data->Lz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Indentor"
PetscErrorCode ModelApplyInitialMaterialGeometry_Indentor(pTatinCtx c,void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscInt               e,p,n_mp_points;
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
		ycoord = ycoord * 1.0e-3; /* m => km */
		
		if (ycoord<72.0) {
			phase = 3;
			eta =  data->eta[3];
			rho =  data->rho[3];
		} else if (ycoord<(72.0+12.0)) {
			phase = 2;
			eta =  data->eta[2];
			rho =  data->rho[2];
		} else if (ycoord<(72.0+12.0+18.0)) {
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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Indentor"
PetscErrorCode ModelApplyUpdateMeshGeometry_Indentor(pTatinCtx c,Vec X,void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Indentor_CheckScales"
PetscErrorCode ModelOutput_Indentor_CheckScales(pTatinCtx c,Vec X)
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

	PetscPrintf(PETSC_COMM_WORLD,"[indentor]: check scales \n");

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
#define __FUNCT__ "ModelOutput_Indentor"
PetscErrorCode ModelOutput_Indentor(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = ModelOutput_Indentor_CheckScales(c,X);CHKERRQ(ierr);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Indentor"
PetscErrorCode ModelDestroy_Indentor(pTatinCtx c,void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Indentor"
PetscErrorCode pTatinModelRegister_Indentor(void)
{
	ModelIndentorCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelIndentorCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelIndentorCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"indentor");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Indentor);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Indentor);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Indentor);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Indentor);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Indentor);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Indentor);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Indentor);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Indentor);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Indentor);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}



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
	PetscErrorCode ierr;
	
	PetscFunctionBegin;


	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
  rheology                = &c->rheology_constants;
	rheology->rheology_type = RHEOLOGY_VISCOUS;
	
	/* box geometry */
	data->Lx = 5.0;
	data->Ly = 1.2;
	data->Lz = 5.0;
	
	/* bc type */
	data->boundary_conditon_type = VSBC_FreeSlip;
		
	/* parse from command line */
//	max_eta = 6.0*1.0e23;
	max_eta = 1.0*1.0e22;
	rheology->const_eta0[0] = 2.0*1.0e23/max_eta; /* crust */
	rheology->const_eta0[1] = 2.0*1.0e21/max_eta; /* crust - lower */
	rheology->const_eta0[2] = 6.0*1.0e23/max_eta; /* mantle */
	rheology->const_eta0[3] = 4.0*1.0e21/max_eta; /* mantle - lower */

	rheology->const_rho0[0] = 2700.0; /* crust */
	rheology->const_rho0[1] = 2700.0; /* crust */
	rheology->const_rho0[2] = 3000.0; /* mantle */
	rheology->const_rho0[3] = 3000.0; /* mantle */
	
	rheology->nphases_active = 2;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_eta0",&rheology->const_eta0[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_eta1",&rheology->const_eta0[1],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_eta2",&rheology->const_eta0[2],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_eta3",&rheology->const_eta0[3],&flg);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_rho0",&rheology->const_rho0[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_rho1",&rheology->const_rho0[1],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_rho2",&rheology->const_rho0[2],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_rho3",&rheology->const_rho0[3],&flg);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_Lx",&data->Lx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_Ly",&data->Ly,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_Lz",&data->Lz,&flg);CHKERRQ(ierr);
	
	
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
	
	ierr = PetscOptionsGetInt(PETSC_NULL,"-model_indentor_bc_type",(PetscInt*)&data->boundary_conditon_type,&flg);CHKERRQ(ierr);
	
	data->cutoff_time = 1e8;
	data->indentation_velocity = 1.0;
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_cutofftime",&data->cutoff_time,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_indentor_indentation_velocity",&data->indentation_velocity,&flg);CHKERRQ(ierr);
	
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
		double      *position;
		double      eta,rho;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		if (position[1]<0.75) {
			phase = 3;
			eta =  data->eta[3];
			rho = -data->rho[3]*9.8;
		} else if (position[1]<0.85) {
			phase = 2;
			eta =  data->eta[2];
			rho = -data->rho[2]*9.8;
		} else if (position[1]<1.05) {
			phase = 1;
			eta =  data->eta[1];
			rho = -data->rho[1]*9.8;
		} else {
			phase = 0;
			eta =  data->eta[0];
			rho = -data->rho[0]*9.8;
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
#define __FUNCT__ "ModelOutput_Indentor"
PetscErrorCode ModelOutput_Indentor(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelIndentorCtx *data = (ModelIndentorCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

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

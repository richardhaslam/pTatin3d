

#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"

#include "viscous_sinker_ctx.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_ViscousSinker"
PetscErrorCode ModelInitialize_ViscousSinker(pTatinCtx c,void *ctx)
{
	ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  RheologyConstants      *rheology;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;


	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
  rheology                = &c->rheology_constants;
	rheology->rheology_type = RHEOLOGY_VISCOUS;
	
	/* box geometry */
	data->Lx = 2.0;
	data->Ly = 3.0;
	data->Lz = 4.0;
	
	/* inclusion geometry */
	data->is_sphere = PETSC_TRUE;

	data->origin[0] = 0.5 * data->Lx;
	data->origin[1] = 0.5 * data->Ly;
	data->origin[2] = 0.5 * data->Lz;
	
	/* spheriod diameter or box length */
	data->length[0] = 0.5 * data->Lx;
	data->length[1] = 0.5 * data->Ly;
	data->length[2] = 0.5 * data->Lz;
	
	/* bc type */
	data->boundary_conditon_type = VSBC_FreeSlip;
		
	/* parse from command line */
	rheology->nphases_active = 2;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_eta0",&rheology->const_eta0[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_eta1",&rheology->const_eta0[1],&flg);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_rho0",&rheology->const_rho0[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_rho1",&rheology->const_rho0[1],&flg);CHKERRQ(ierr);

	/* set initial values for model parameters */
	/* material properties */
	data->nmaterials = rheology->nphases_active;
	data->eta0 = rheology->const_eta0[0];
	data->eta1 = rheology->const_eta0[1];
	data->rho0 = rheology->const_rho0[0];
	data->rho1 = rheology->const_rho0[1];
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_Lx",&data->Lx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_Ly",&data->Ly,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_Lz",&data->Lz,&flg);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_viscous_sinker_cube",&data->is_sphere,&flg);CHKERRQ(ierr);
	if (flg==PETSC_TRUE) { data->is_sphere = PETSC_FALSE; }
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_Ox",&data->origin[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_Oy",&data->origin[1],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_Oz",&data->origin[2],&flg);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_lx",&data->length[0],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_ly",&data->length[1],&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_viscous_sinker_lz",&data->length[2],&flg);CHKERRQ(ierr);

	ierr = PetscOptionsGetInt(PETSC_NULL,"-model_viscous_sinker_bc_type",(PetscInt*)&data->boundary_conditon_type,&flg);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_ViscousSinker"
PetscErrorCode ModelApplyBoundaryCondition_ViscousSinker(pTatinCtx user,void *ctx)
{
	ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
	PetscScalar zero = 0.0;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	switch (data->boundary_conditon_type) {
		case VSBC_FreeSlip:
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case VSBC_NoSlip:
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case VSBC_FreeSlipFreeSurface:
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case VSBC_NoSlipFreeSurface:
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		default:
			break;
	}
	
/*	
	{
		BCList flat;
		
		ierr = BCListFlattenedCreate(user->stokes_ctx->u_bclist,&flat);CHKERRQ(ierr);
		ierr = BCListDestroy(&user->stokes_ctx->u_bclist);CHKERRQ(ierr);
		user->stokes_ctx->u_bclist = flat;
	}
*/	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_ViscousSinker"
PetscErrorCode ModelApplyMaterialBoundaryCondition_ViscousSinker(pTatinCtx c,void *ctx)
{
	ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_ViscousSinker"
PetscErrorCode ModelApplyInitialMeshGeometry_ViscousSinker(pTatinCtx c,void *ctx)
{
	ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,0.0,data->Lx, 0.0,data->Ly, 0.0,data->Lz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_ViscousSinker"
PetscErrorCode ModelApplyInitialMaterialGeometry_ViscousSinker(pTatinCtx c,void *ctx)
{
	ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
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
		
		phase = 0;
		eta =  data->eta0;
		rho = -data->rho0;
		
		if (data->is_sphere) {
			double rx = (position[0]-data->origin[0])/(0.5*data->length[0]);
			double ry = (position[1]-data->origin[1])/(0.5*data->length[1]);
			double rz = (position[2]-data->origin[2])/(0.5*data->length[2]);
			
			if ( rx*rx + ry*ry + rz*rz < 1.0 ) {
				phase = 1;
				eta =  data->eta1;
				rho = -data->rho1;
			}
			
		} else { /* box */
			if ( (position[0]>data->origin[0] - 0.5*data->length[0]) && (position[0]<data->origin[0] + 0.5*data->length[0]) ) {
				if ( (position[1]>data->origin[1] - 0.5*data->length[1]) && (position[1]<data->origin[1] + 0.5*data->length[1]) ) {
					if ( (position[2]>data->origin[2] - 0.5*data->length[2]) && (position[2]<data->origin[1] + 0.5*data->length[2]) ) {
						phase = 1;
						eta =  data->eta1;
						rho = -data->rho1;
					}
				}
			}
			
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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_ViscousSinker"
PetscErrorCode ModelApplyUpdateMeshGeometry_ViscousSinker(pTatinCtx c,Vec X,void *ctx)
{
	ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_ViscousSinker"
PetscErrorCode ModelOutput_ViscousSinker(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_ViscousSinker"
PetscErrorCode ModelDestroy_ViscousSinker(pTatinCtx c,void *ctx)
{
	ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_ViscousSinker"
PetscErrorCode pTatinModelRegister_ViscousSinker(void)
{
	ModelViscousSinkerCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelViscousSinkerCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelViscousSinkerCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"viscous_sinker");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_ViscousSinker);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_ViscousSinker);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_ViscousSinker);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_ViscousSinker);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_ViscousSinker);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_ViscousSinker);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_ViscousSinker);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_ViscousSinker);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

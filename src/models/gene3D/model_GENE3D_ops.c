

#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "Phase_map.h"

#include "GENE3D_ctx.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_GENE3D"
PetscErrorCode ModelInitialize_GENE3D(pTatinCtx c,void *ctx)
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
  	RheologyConstants      *rheology;
	PetscBool flg;
    char   				   *option_name;
    char   				   model = "model_GENE3D_" ;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;


	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
        rheology                = &c->rheology_constants;


	/* model geometry */
	PetscPrintf(PETSC_COMM_WORLD,"reading model initial geometry from options\n");
 	ierr       = PetscOptionsGetInt(model,"-initial_geom",&data->initial_geom,&found);CHKERRQ(ierr);
   	if (found == PETSC_FALSE){
   	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide a type of material index initialisation \n");
  	 }

	/* box geometry */
        PetscPrintf(PETSC_COMM_WORLD,"reading box geometry from options\n");

	ierr = PetscOptionsGetReal(model,"-Lx",&data->Lx,&flg);CHKERRQ(ierr);
        if (found == PETSC_FALSE); {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide model length Lx \n";}

	ierr = PetscOptionsGetReal(model,"-Ly",&data->Ly,&flg);CHKERRQ(ierr);
        if (found == PETSC_FALSE); {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide model length Ly \n";}

	ierr = PetscOptionsGetReal(model,"-Lz",&data->Lz,&flg);CHKERRQ(ierr);
        if (found == PETSC_FALSE); {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide model length Lz \n";}


	ierr = PetscOptionsGetReal(model,"-Ox",&data->Ox,&flg);CHKERRQ(ierr);
        if (found == PETSC_FALSE); {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide model Origin Ox \n";}

	ierr = PetscOptionsGetReal(model,"-Oy",&data->Oy,&flg);CHKERRQ(ierr);
        if (found == PETSC_FALSE); {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide model Origin Oy \n";}

	ierr = PetscOptionsGetReal(model,"-Oz",&data->Oz,&flg);CHKERRQ(ierr);
        if (found == PETSC_FALSE); {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide model Origin Oz \n";}
	


	/* rheology type */
  	PetscPrintf(PETSC_COMM_WORLD,"reading rheology type from options\n");
 	ierr       = PetscOptionsGetInt(model,"-rheol",&rheology->rheology_type,&found);CHKERRQ(ierr);
   	if (found == PETSC_FALSE){
   	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value for rheology type \n");
  	 }

	/* material properties */
	PetscPrintf(PETSC_COMM_WORLD,"reading material properties from options\n");
	ierr       = PetscOptionsGetInt(model,"-nphase",&nphase,&found);CHKERRQ(ierr);
   	if (found == PETSC_FALSE){
   	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value for number of materials \n");
   	}
	for (i=0; i< nphase ;i++) {
   
		switch  (rheology->rheology_type) {
		case 0 :{
			asprintf(&option_name, "-eta_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->const_eta0[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

			asprintf(&option_name, "-rho_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->const_rho0[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);
			}
		break; 
		case 1: {
			asprintf(&option_name, "-eta_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->const_eta0[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

			asprintf(&option_name, "-rho_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->const_rho0[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);
        			
			asprintf(&option_name, "-Co_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->mises_tau_yield[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

 			asprintf(&option_name, "-Phi_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->dp_pressure_dependance[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

 			asprintf(&option_name, "-Tens_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->tens_cutoff[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

 			asprintf(&option_name, "-Hs_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->Hst_cutoff[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);
			}
		break; 
		case 2 : {
			asprintf(&option_name, "-Co_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->mises_tau_yield[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

 			asprintf(&option_name, "-Phi_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->dp_pressure_dependance[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

 			asprintf(&option_name, "-Tens_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->tens_cutoff[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

 			asprintf(&option_name, "-Hs_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->Hst_cutoff[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

			asprintf(&option_name, "-Co_inf_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->soft_Co_inf[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

 			asprintf(&option_name, "-Phi_inf_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->soft_phi_inf[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
    		free(option_name);

 			asprintf(&option_name, "-eps_min_%d",i);
   			ierr = PetscOptionsGetReal(model,option_name,&(rheology->soft_min_strain_cutoff[i]),&found);CHKERRQ(ierr);
   			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
			free(option_name);

 			asprintf(&option_name, "-eps_max_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&(rheology->soft_max_strain_cutoff[i]),&found);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
			free(option_name);
			}
		break;
		}
	}

	/* bc type */
	data->boundary_conditon_type = VSBC_FreeSlip;
		

	/* set initial values for model parameters */


	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_GENE3D"
PetscErrorCode ModelApplyBoundaryCondition_GENE3D(pTatinCtx user,void *ctx)
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
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
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_GENE3D"
PetscErrorCode ModelApplyMaterialBoundaryCondition_GENE3D(pTatinCtx c,void *ctx)
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_GENE3D"
PetscErrorCode ModelApplyInitialMeshGeometry_GENE3D(pTatinCtx c,void *ctx)
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,data->Ox,data->Lx, data->Oy,data->Ly, data->Oz,data->Lz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_GENE3D"
PetscErrorCode ModelApplyInitialMaterialGeometry_GENE3D(pTatinCtx c,void *ctx)
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
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

//=====================================================================================================================================

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_ModelSetMarkerIndexLayerCake_GENE3D"
PetscErrorCode pTatin3d_ModelSetMarkerIndexFromMap_GENE3D(pTatinCtx c,void *ctx)
/* define phase index on material points from a map file extruded in z direction */
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
	PetscInt               p,n_mp_points,i,nlayer;
	DataBucket             db;
	DataField              PField_std;
	int                    phase_init, phase, phase_index, is_valid;
	Petsc_int              phaseLayer[LAYER_MAX];
	PetscScalar            Ylayer[LAYER_MAX+1];
	char                   *option_name;

	PetscFunctionBegin;
  	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* define properties on material points */
	db = ctx->db;
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	DataBucketGetSizes(db,&n_mp_points,0,0);

	/* read layers from options */
	nlayer = 1; 
	ierr = PetscOptionsGetInt("Model_GENE3D_","-nlayer",&nlayer,&flg);CHKERRQ(ierr);
    Ylayer[0] = Oy ; 
	for (i=1,i<nlayer,i++){	
			asprintf(&option_name, "-Ylayer_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&YLayer[i],&flg);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
			free(option_name);
	if (YLayer[i] > Ylayer[i-1]){
			asprintf(&option_name, "-phaseLayer_%d",i);
			ierr = PetscOptionsGetReal(model,option_name,&phaseLayer[i],&flg);CHKERRQ(ierr);
			if (found == PETSC_FALSE) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
			free(option_name);
	} else { 	
    	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Layers must be entered so that ylayer_[i-1] is larger than ylayer_[i]\n");
	}


	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		double  Ypos;
		DataFieldAccessPoint(PField_std,p,(void**)&material_point);
		Ypos = material_point->coor[1];
		i=0
		while (Ypos >= Ylayer[i] & i < nlayer){ 		
		 i++
		}
       phase = phaseLayer[i]; 
  
		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,phase);

	}
	
	DataFieldRestoreAccess(PField_std);
	PetscFunctionReturn(0);
}
//===============================================================================================================================
#undef __FUNCT__  
#define __FUNCT__ "ModelSetMarkerIndexFromMap_GENE3D"
PetscErrorCode ModelSetMarkerIndexFromMap_GENE3D(pTatinCtx c,void *ctx)
/* define phase index on material points from a map file extruded in z direction */
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
	PhaseMap               phasemap;
	PetscInt               p,n_mp_points;
	DataBucket             db;
	DataField              PField_std;
	int                    phase_init, phase, phase_index, is_valid;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
		 

	ierr = PetscOptionsGetString(Model_GENE3D_,"-map_file",map_file,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  	if (flg == PETSC_FALSE) {
    	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide a map file \n");
  	}

  	asprintf(&name,"./inputdata/%s.pmap",map_file);
	PhaseMapLoadFromFile(name,&phasemap);
  	asprintf(&name,"./inputdata/%s_phase_map.gp",map_file);
	PhaseMapViewGnuplot(name,phasemap);
	free(name);


	/* define properties on material points */
	db = ctx->db;
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		double  XYposition[2];

		DataFieldAccessPoint(PField_std,p,(void**)&material_point);

		XYposition[0] = material_point->coor[0];
		XYposition[1] = material_point->coor[1];
		
		MPntStdGetField_phase_index(material_point,&phase_init);

		PhaseMapGetPhaseIndex(phasemap,XYposition,&phase_index);

		PhaseMapCheckValidity(phasemap,phase_index,&is_valid);
    //PetscPrintf(PETSC_COMM_WORLD,"Phase index : %d  is_valid %d \n", phase_index,is_valid);

		if (is_valid==1) { /* point located in the phase map */
			phase = phase_index;
		} else {
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"marker outside the domain\n your phasemap is smaller than the domain \n please check your parameters and retry");
		}	
		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,phase);

	}
	
	DataFieldRestoreAccess(PField_std);
	
	
	
	PetscFunctionReturn(0);
}


//======================================================================================================================================

#undef __FUNCT__  
#define __FUNCT__ "ModelSetInitialStokesVariableOnMarker_GENE3D"
PetscErrorCode ModelSetInitialStokesVariableOnMarker_GENE3D(pTatinCtx c,void *ctx)
/* define properties on material points */
{
	ModelGENE3DCtx *data = 	(ModelGENE3DCtx*)ctx;
	PhaseMap               	phasemap;
	PetscInt               	e,p,n_mp_points;
	DataBucket            	db;
	DataField             	PField_std,PField_stokes;
	int                    	phase_index, i;
	PetscErrorCode         	ierr;
	RheologyConstants      	*rheology;
    char 					*name;
	      
	PetscFunctionBegin;
  	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);



	rheology   = &ctx->rheology_constants;
  	db = ctx->db;

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
		DataFieldAccessPoint(PField_std,p,(void**)&material_point);
  		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
    	MPntStdGetField_phase_index(material_point,&phase_index);
    	MPntPStokesSetField_eta_effective(mpprop_stokes,rheology->const_eta0[phase_index]);
		MPntPStokesSetField_density(mpprop_stokes,rheology->const_rho0[phase_index]);	
	}

	DataFieldRestoreAccess(PField_std);
  	DataFieldRestoreAccess(PField_stokes);
	
	PetscFunctionReturn(0);
}

//======================================================================================================================================

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_GENE3D"
PetscErrorCode ModelApplyInitialMaterialGeometry_GENE3D(pTatinCtx c,void *ctx)
{
PetscErrorCode         ierr;
ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
PetscFunctionBegin;

	switch(data->initial_geom){
		/*Layered cake*/
		case 0 :{
		ierr = ModelSetMarkerIndexLayeredCake_GENE3D(c,ctx); CHKERRQ(ierr);
		}
		break  
 		/*Extrude from Map along Z*/
		case 1 :{
		ierr = ModelSetMarkerIndexFromMap_GENE3D(c,ctx); CHKERRQ(ierr);	 		
		}
		break; 
		/*Read from CAD file*/
		case 2 :{
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Reading from CAD is not implemented yet \n");
		}
		break;
	}
		ierr = ModelSetInitialStokesVariableOnMarker_GENE3D(c,ctx); CHKERRQ(ierr);

PetscFunctionReturn(0);
}



//======================================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_GENE3D"
PetscErrorCode ModelApplyUpdateMeshGeometry_GENE3D(pTatinCtx c,Vec X,void *ctx)
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_GENE3D"
PetscErrorCode ModelOutput_GENE3D(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_GENE3D"
PetscErrorCode ModelDestroy_GENE3D(pTatinCtx c,void *ctx)
{
	ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_GENE3D"
PetscErrorCode pTatinModelRegister_GENE3D(void)
{
	ModelGENE3DCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelGENE3DCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelGENE3DCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"GENE3D");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_GENE3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_GENE3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_GENE3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_GENE3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_GENE3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_GENE3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_GENE3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_GENE3D);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

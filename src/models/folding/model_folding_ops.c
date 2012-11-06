

#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "swarm_fields.h"

#include "dmda_bcs.h"
#include "dmda_remesh.h"
#include "dmda_update_coords.h"
#include "mesh_update.h"

#include "folding_ctx.h"




#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Folding"
PetscErrorCode ModelInitialize_Folding(pTatinCtx c,void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
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
	data->model_type = FOLD_MODEL_A;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-model_type",(PetscInt*)&data->model_type,&flg);CHKERRQ(ierr);
	
	data->dimensional   = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_dimensional",&data->dimensional,&flg);CHKERRQ(ierr);

	/* box geometry */
	data->Lx = 45.0;
	data->Ly = 45.0;
	data->Lz = 45.0;

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_Lx",&data->Lx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_Ly",&data->Ly,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_Lz",&data->Lz,&flg);CHKERRQ(ierr);

	data->cutoff_time = 1.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_cutofftime",&data->cutoff_time,&flg);CHKERRQ(ierr);


	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_Folding"
PetscBool BCListEvaluator_Folding( PetscScalar position[], PetscScalar *value, void *ctx ) 
{
	PetscBool impose_dirichlet = PETSC_TRUE;
	pTatinCtx user = (pTatinCtx)ctx;
	PetscScalar vx,Vx_max;
	PetscReal time,Dprime,cc;
	ModelFoldingCtx *model_data_ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = pTatinModelGetUserData(user->model,(void**)&model_data_ctx);CHKERRQ(ierr);
	time = user->time;


	*value = vx;
	return impose_dirichlet;
}


#undef __FUNCT__
#define __FUNCT__ "ApplyDirichletBoundaryConditions_Folding"
PetscErrorCode ApplyDirichletBoundaryConditions_Folding(ModelFoldingCtx *data,pTatinCtx user,BCList u_bclist,DM dav)
{
	PetscErrorCode ierr;

	PetscFunctionBegin;

	
	switch (data->model_type) {

		case FOLD_MODEL_A:
		{
			PetscScalar zero = 0.0;
			PetscScalar Vx   = 1.0;
			PetscScalar Vz   = 2.0;

			/* free slip base */
			ierr = DMDABCListTraverse3d(u_bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			/* free surface top */
			//ierr = DMDABCListTraverse3d(u_bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

			
			ierr = DMDABCListTraverse3d(u_bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(u_bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			ierr = DMDABCListTraverse3d(u_bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&Vx);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(u_bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&Vz);CHKERRQ(ierr);
						
			break;
		}
			
			
	}
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Folding"
PetscErrorCode ModelApplyBoundaryCondition_Folding(pTatinCtx user,void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = ApplyDirichletBoundaryConditions_Folding(data,user,user->stokes_ctx->u_bclist,user->stokes_ctx->dav);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Folding"
PetscErrorCode ModelApplyBoundaryConditionMG_Folding(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	for (n=0; n<nl; n++) {
		ierr = ApplyDirichletBoundaryConditions_Folding(data,user,bclist[n],dav[n]);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Folding"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Folding(pTatinCtx c,void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

/*
 
y = (y0 + j*dy) + G(x,z) * alpha
 
alpha = 1.0 at y = y_ref
 and 
alpha = 0.0 at y = y_min, y_max
 
*/
#undef __FUNCT__
#define __FUNCT__ "MeshDeformation_setup1"
PetscErrorCode MeshDeformation_setup1(DM da)
{
	PetscErrorCode ierr;
	PetscInt si,sj,sk,nx,ny,nz,i,j,k,MY;
	DM cda;
	Vec coord;
	DMDACoor3d ***_coord;
	PetscReal g_xz, alpha, y_ref, dy;
	PetscReal Gmin[3],Gmax[3];
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(da,Gmin,Gmax);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(da,0,0,&MY,0, 0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	/* constant spacing */
	dy = (Gmax[1]-Gmin[1])/(PetscReal)(MY-1);

	ierr = DMDAGetCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&_coord);CHKERRQ(ierr);
	
	/* apply bump */
	for( j=sj; j<sj+ny; j++ ) {
		for( k=sk; k<sk+nz; k++ ) {
			for( i=si; i<si+nx; i++ ) {
				PetscReal xn,yn,zn;
				
				xn = _coord[k][j][i].x;
				yn = _coord[k][j][i].y;
				zn = _coord[k][j][i].z;

				y_ref = 0.5 * 45.0;
				g_xz = 4.1 * sin((1.0/20.0) * xn * M_PI) * cos( (1.0/10.0)*zn * M_PI);
				

				if (yn < y_ref) {
					alpha = (1.0/(y_ref-Gmin[1])) * yn;
				} else {
					alpha = (-1.0/(Gmax[1]-y_ref)) * yn + 1.0 + (y_ref/(Gmax[1]-y_ref));
				}
				
				_coord[k][j][i].y = (Gmin[1] + j * dy) + g_xz * alpha;
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,coord,&_coord);CHKERRQ(ierr);
	/* update */
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

void interface_layer_a(PetscReal xn[],PetscReal *yn)
{
	PetscReal g_xz, y_ref;

	y_ref = 0.5 * 45.0;
	g_xz = 4.1 * sin((1.0/20.0) * xn[0] * M_PI) * cos( (1.0/10.0)*xn[2] * M_PI);
	
	*yn = y_ref + g_xz;
}

void interface_layer_b(PetscReal xn[],PetscReal *yn)
{
	PetscReal g_xz, y_ref;
	
	y_ref = 0.5 * 45.0;
	g_xz = 4.1 * sin((1.0/20.0) * xn[0] * M_PI) * cos( (1.0/10.0)*xn[2] * M_PI) + 2.0;
	
	*yn = y_ref + g_xz;
}

#undef __FUNCT__
#define __FUNCT__ "MeshDeformation_insert_layer"
PetscErrorCode MeshDeformation_insert_layer(DM da,const PetscInt j,void (*interface_func)(PetscReal*,PetscReal*))
{
	PetscErrorCode ierr;
	PetscInt si,sj,sk,nx,ny,nz,i,k,MY;
	DM cda;
	Vec coord;
	DMDACoor3d ***_coord;
	PetscReal Gmin[3],Gmax[3];
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(da,Gmin,Gmax);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(da,0,0,&MY,0, 0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&_coord);CHKERRQ(ierr);
	
	/* apply bump */
	if( (j>=sj) && (j<sj+ny) ) {
		
		
		for( k=sk; k<sk+nz; k++ ) {
			for( i=si; i<si+nx; i++ ) {
				PetscReal xn[3],yn;
				
				xn[0] = _coord[k][j][i].x;
				xn[1] = _coord[k][j][i].y;
				xn[2] = _coord[k][j][i].z;

				interface_func(xn,&yn);
				_coord[k][j][i].y = yn;
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,coord,&_coord);CHKERRQ(ierr);
	/* update */
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Folding"
PetscErrorCode ModelApplyInitialMeshGeometry_Folding(pTatinCtx c,void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
	DM dm;
	PetscInt MX,MY,MZ;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	dm = c->stokes_ctx->dav;
	
	ierr = DMDASetUniformCoordinates(dm,0.0,data->Lx, 0.0,data->Ly, 0.0,data->Lz);CHKERRQ(ierr);
//	ierr = MeshDeformation_setup1(c->stokes_ctx->dav);CHKERRQ(ierr);

	
	ierr = DMDAGetInfo(dm,0,&MX,&MY,&MZ, 0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
/*	
	ierr = MeshDeformation_insert_layer(dm,5,interface_layer_a);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dm,0,6);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dm,5,MY);CHKERRQ(ierr);
*/

	ierr = MeshDeformation_insert_layer(dm,20,interface_layer_a);CHKERRQ(ierr);	
	ierr = MeshDeformation_insert_layer(dm,42,interface_layer_b);CHKERRQ(ierr);

	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dm,0,21);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dm,20,43);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dm,42,MY);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Folding"
PetscErrorCode ModelApplyInitialMaterialGeometry_Folding(pTatinCtx c,void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
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
		ycoord = position[1];
		
		phase = 0;
		eta = 1.0e-2;
		if ( (ycoord >= 20.0) && (ycoord <= 25.0) ) {
			phase = 1;
			eta =  1.0;
		}
		rho = 0.0;

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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Folding"
PetscErrorCode ModelApplyUpdateMeshGeometry_Folding(pTatinCtx c,Vec X,void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
	Vec velocity,pressure;
	DM dav,dap;
	PetscReal step;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	step = c->dt;

	ierr = DMCompositeGetEntries(c->stokes_ctx->stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(c->stokes_ctx->stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(c->stokes_ctx->stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_CheckScales"
PetscErrorCode ModelOutput_CheckScales(pTatinCtx c,Vec X)
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

	PetscPrintf(PETSC_COMM_WORLD,"[folding]: check scales \n");

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
#define __FUNCT__ "ModelOutput_Folding"
PetscErrorCode ModelOutput_Folding(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	//ierr = ModelOutput_CheckScales(c,X);CHKERRQ(ierr);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Folding"
PetscErrorCode ModelDestroy_Folding(pTatinCtx c,void *ctx)
{
	ModelFoldingCtx *data = (ModelFoldingCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Folding"
PetscErrorCode pTatinModelRegister_Folding(void)
{
	ModelFoldingCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelFoldingCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelFoldingCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"folding");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Folding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Folding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Folding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Folding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Folding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Folding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Folding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Folding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Folding);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_models.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_update_coords.h"
#include "dmda_element_q2p1.h"
#include "mesh_update.h"
#include "dmda_remesh.h"
#include "output_material_points.h"
#include "mesh_quality_metrics.h"

#include "model_wrench_fold_ctx.h"
#include "model_utils.h"
#include "math.h"

const PetscReal L_star = 1.0e2;


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_WrenchFold"
PetscErrorCode ModelInitialize_WrenchFold(pTatinCtx c,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	PetscInt n_int,n;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* assign defaults */
	data->max_layers = 100;
	
	data->n_interfaces = 2;
	PetscOptionsGetInt(NULL,"-model_wrench_fold_n_interfaces",&data->n_interfaces,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the number of interfaces including the top and bottom boundaries (-model_wrench_fold_n_interfaces)");
	}
	
	pTatinModelGetOptionReal("-model_wrench_fold_Lx", &data->Lx, "User must provide the length along the x direction", NULL,PETSC_TRUE);
	/*PetscOptionsGetReal(NULL,"-model_wrench_fold_Lx",&data->Lx,&flg);
	 if (!flg) {
	 SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the x direction (-model_wrench_fold_Lx)");
	 }*/
	
	PetscOptionsGetReal(NULL,"-model_wrench_fold_Ly",&data->Ly,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the y direction (-model_wrench_fold_Ly)");
	}

	/* scale length */
	data->Lx = data->Lx / L_star;
	data->Ly = data->Ly / L_star;
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_wrench_fold_interface_heights",data->interface_heights,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide interface heights relative from the base of the model including the top and bottom boundaries (-model_wrench_fold_interface_heights)");
	}
	if (n_int != data->n_interfaces) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d interface heights relative from the base of the model including the top and bottom boundaries (-model_wrench_fold_interface_heights)",data->n_interfaces);
	}
	/* scale length */
	for (n=0; n<data->n_interfaces; n++) {
		data->interface_heights[n] = data->interface_heights[n] / L_star;
	}
	
	n_int = data->max_layers;
	PetscOptionsGetIntArray(NULL,"-model_wrench_fold_layer_res_k",data->layer_res_k,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer resolution list (-model_wrench_fold_layer_res_k)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_wrench_fold_layer_res_k)",data->n_interfaces-1);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_wrench_fold_layer_eta",data->eta,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list (-model_wrench_fold_layer_eta)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity (-model_wrench_fold_layer_eta)",data->n_interfaces-1);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_wrench_fold_layer_rho",data->rho,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer density list (-model_wrench_fold_layer_rho)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer density (-model_wrench_fold_layer_rho)",data->n_interfaces-1);
	}
	
	/* define the mesh size the z-direction for the global problem */
	c->mz = 0;
	for (n=0; n<data->n_interfaces-1; n++) {
		c->mz += data->layer_res_k[n];
	}
	
	/* define the domain size in the z-direction for the global problem */
	data->Lz = data->interface_heights[data->n_interfaces-1];
	
	
	data->bc_type = 0; /* 0 for basal shear ; 1 for lateral shear */
	data->exy             = -1.0e-3;
	
	/* parse from command line or input file */
	ierr = PetscOptionsGetInt(NULL,"-model_wrench_fold_bc_type",&data->bc_type,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_wrench_fold_exy",&data->exy,&flg);CHKERRQ(ierr);
	
	/* define the coefficient for the BC.*/
	if(data->bc_type == 0){
		data->A = (data->exy * data->Ly/(atan(data->Ly/2.0) - atan(-data->Ly/2.0)))/(M_PI*0.5);
	}else if(data->bc_type == 1){
		data->A = data->exy * data->Ly;
	}
	PetscPrintf(PETSC_COMM_WORLD,"ModelReport: \"Wrench Fold\"\n");
	PetscPrintf(PETSC_COMM_WORLD," Domain: [0 , %1.4e] x [0 , %1.4e] x [0 , %1.4e] ", data->Lx,data->Ly,data->Lz );
	PetscPrintf(PETSC_COMM_WORLD," Mesh:   %.4D x %.4D x %.4D \n", c->mx,c->my,c->mz ); 
	for (n=data->n_interfaces-1; n>=1; n--) {
		PetscPrintf(PETSC_COMM_WORLD," ---------------------------- z = %1.4e ----------------------------\n",data->interface_heights[n]);
		PetscPrintf(PETSC_COMM_WORLD,"|\n"); 
		PetscPrintf(PETSC_COMM_WORLD,"|      eta = %1.4e , rho = %1.4e , mz = %.4D \n",data->eta[n-1],data->rho[n-1],data->layer_res_k[n-1]);
		PetscPrintf(PETSC_COMM_WORLD,"|\n");
	}
	//PetscPrintf(PETSC_COMM_WORLD,"|\n");
	PetscPrintf(PETSC_COMM_WORLD," ---------------------------- z = %1.4e ----------------------------\n",data->interface_heights[0],data->layer_res_k[0]);
	
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_WrenchFold"
PetscBool BCListEvaluator_WrenchFold( PetscScalar position[], PetscScalar *value, void *ctx ) 
{
	PetscBool impose_dirichlet = PETSC_TRUE;
	pTatinCtx user = (pTatinCtx)ctx;
	PetscScalar exy;
	PetscReal A, Ly, vx;
	ModelWrenchFoldCtx *model_data_ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = pTatinModelGetUserData(user->model,(void**)&model_data_ctx);CHKERRQ(ierr);
	//time = user->time;
	
	A = model_data_ctx->A;
	Ly = model_data_ctx->Ly;
	
	if (model_data_ctx->bc_type == 0) {
		vx = A*atan(position[1] - 0.5*Ly);
	} else if (model_data_ctx->bc_type == 1){
		vx = (position[1] - 0.5*Ly>0)?A:-A;
	}else{
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknonwn boundary condition type: 0 for basal shear, 1 for lateral shear (-model_wrench_fold_bc_type)");
		
	}
	*value = vx;
	return impose_dirichlet;
}

#undef __FUNCT__
#define __FUNCT__ "BoundaryCondition_WrenchFold"
PetscErrorCode BoundaryCondition_WrenchFold(DM dav,BCList bclist,pTatinCtx user,ModelWrenchFoldCtx *data)
{
	PetscScalar zero = 0.0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	if(data->bc_type == 0){
        
		/* inflow/outflow faces constrained such that flow has zero y,z-compoenent */
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
        
		/* free slip lateral */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		/* free surface top */
		
		/* arctan bottom */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_WrenchFold,(void*)user);CHKERRQ(ierr);


	}else if(data->bc_type == 1){
		/* lateral shear*/
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_WrenchFold,(void*)user);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_WrenchFold,(void*)user);CHKERRQ(ierr);
		
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		/* free surface top */

		/* free slip bottom */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr); 

		/* inflow/outflow faces constrained such that flow has zero y,z-compoenent */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_WrenchFold"
PetscErrorCode ModelApplyBoundaryCondition_WrenchFold(pTatinCtx user,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	BCList            bclist;
	DM                dav;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	bclist = user->stokes_ctx->u_bclist;
	dav    = user->stokes_ctx->dav;
	ierr = BoundaryCondition_WrenchFold(dav,bclist,user,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_WrenchFold"
PetscErrorCode ModelApplyBoundaryConditionMG_WrenchFold(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		ierr = BoundaryCondition_WrenchFold(dav[n],bclist[n],user,data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_WrenchFold"
PetscErrorCode ModelApplyMaterialBoundaryCondition_WrenchFold(pTatinCtx c,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WrenchFoldSetPerturbedInterfaces"
PetscErrorCode WrenchFoldSetPerturbedInterfaces(DM dav, PetscScalar interface_heights[], PetscInt layer_res_k[], PetscInt n_interfaces,PetscReal amp)
{
	PetscErrorCode ierr;
	PetscInt i,j,si,sj,sk,nx,ny,nz, interf, kinter, proc, M,N,P;
	PetscScalar *random, dz;
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	
	kinter = 0;
	ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dav),&proc);CHKERRQ(ierr);
	srand(proc+2);
	for(interf = 1; interf < n_interfaces-1; interf++){
		kinter += 2*layer_res_k[interf-1];
		PetscPrintf(PETSC_COMM_WORLD,"kinter = %d (max=%d)\n", kinter,P-1 );
		if ( (kinter>=sk) && (kinter<sk+nz) ) {
			
			dz = 0.5*((interface_heights[interf+1] - interface_heights[interf])/(PetscScalar)(layer_res_k[interf]) + (interface_heights[interf] - interface_heights[interf-1])/(PetscScalar)(layer_res_k[interf-1]) );
			PetscPrintf(PETSC_COMM_SELF," interface %d: using dz computed from avg %1.4e->%1.4e / mz=%d :: %1.4e->%1.4e / mz=%d \n", interf,
									interface_heights[interf+1],interface_heights[interf],layer_res_k[interf],
									interface_heights[interf],interface_heights[interf-1],layer_res_k[interf-1] );
			for(i = si; i<si+nx; i++) {
				for(j = sj; j<sj+ny; j++) {
					PetscScalar random = 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
					LA_coord[kinter][j][i].z += amp * dz * random;
				}
			}
			
		}
	}
	
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryMaterialPoints_WrenchFold"
PetscErrorCode InitialMaterialGeometryMaterialPoints_WrenchFold(pTatinCtx c,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
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
		//double      *position;
		PetscReal      eta,rho;
		PetscInt    phase;
		PetscInt    layer, kmaxlayer, kminlayer;
		PetscInt    I, J, K;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		//MPntStdGetField_global_coord(material_point,&position)
		
    MPntGetField_global_element_IJKindex(c->stokes_ctx->dav,material_point, &I, &J, &K);
		phase = -1;
		eta =  0.0;
		rho = 0.0;
		kmaxlayer = kminlayer = 0;
		layer = 0;
		// gets the global element index (i,j,k)
		//....
		
		//Set the properties
		while( (phase == -1) && (layer < data->n_interfaces-1) ){
			kmaxlayer += data->layer_res_k[layer];
			
			if( (K<kmaxlayer) && (K>=kminlayer) ){
				phase = layer + 1;
				eta = data->eta[layer];
				rho = data->rho[layer];
				
				rho = - rho * GRAVITY;
			}
			kminlayer += data->layer_res_k[layer];
			layer++;
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
#define __FUNCT__ "InitialMaterialGeometryQuadraturePoints_WrenchFold"
PetscErrorCode InitialMaterialGeometryQuadraturePoints_WrenchFold(pTatinCtx c,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	PhysCompStokes         user;
	QPntVolCoefStokes      *all_gausspoints,*cell_gausspoints;
	PetscInt               nqp,qp;
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
	
	
	/* get the quadrature points */
	user = c->stokes_ctx;
	ierr = VolumeQuadratureGetAllCellData_Stokes(user->volQ,&all_gausspoints);CHKERRQ(ierr);
	nqp = user->volQ->npoints;
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		//double      *position;
		PetscReal      eta,rho;
		PetscInt    phase;
		PetscInt    layer, kmaxlayer, kminlayer, localeid_p;
		PetscInt    I, J, K;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		
    MPntGetField_global_element_IJKindex(c->stokes_ctx->dav,material_point, &I, &J, &K);
		
		//Set the properties
		phase = -1;
		eta =  0.0;
		rho = 0.0;
		kmaxlayer = kminlayer = 0;
		layer = 0;
		while( (phase == -1) && (layer < data->n_interfaces-1) ){
			kmaxlayer += data->layer_res_k[layer];
			
			if( (K<kmaxlayer) && (K>=kminlayer) ){
				phase = layer + 1;
				eta = data->eta[layer];
				rho = data->rho[layer];
			}
			kminlayer += data->layer_res_k[layer];
			layer++;
		}
		
		
		MPntStdGetField_local_element_index(material_point,&localeid_p);
		ierr = VolumeQuadratureGetCellData_Stokes(user->volQ,all_gausspoints,localeid_p,&cell_gausspoints);CHKERRQ(ierr);
		
		for (qp=0; qp<nqp; qp++) {
			cell_gausspoints[qp].eta  = eta;
			cell_gausspoints[qp].rho  = rho;
			
			cell_gausspoints[qp].Fu[0] = 0.0;
			cell_gausspoints[qp].Fu[1] = -rho * GRAVITY;
			cell_gausspoints[qp].Fu[2] = 0.0;
			
			cell_gausspoints[qp].Fp = 0.0;
		}		
		
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_WrenchFold"
PetscErrorCode ModelApplyInitialMaterialGeometry_WrenchFold(pTatinCtx c,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = InitialMaterialGeometryMaterialPoints_WrenchFold(c,ctx);CHKERRQ(ierr);
	ierr = InitialMaterialGeometryQuadraturePoints_WrenchFold(c,ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_WrenchFold"
PetscErrorCode ModelApplyInitialMeshGeometry_WrenchFold(pTatinCtx c,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	PetscReal         Lx,Ly,dx,dy,dz,Lz;
	PetscInt          mx,my,mz, itf;
	PetscReal         amp,factor;
  char              mesh_outputfile[PETSC_MAX_PATH_LEN];
	Vec               coord;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* step 1 - create structured grid */
	Lx = data->Lx;
	Ly = data->Ly;
	Lz = data->Lz;
	/*
	 The length of the model in z-direction is determined by the grid spacing in x and y.
	 We do this so that the elements do not have a large aspect ratio. This would occur
	 if we hard coded Lz to a constant number which is independnet of the grid resolution in x and y.
	 We choose Lz to be mz * min(dx,dy).
	 */
	
	
	mx = c->mx; 
	my = c->my; 
	mz = c->mz; 
	
	mz=0;
	for(itf = 0; itf<data->n_interfaces -1; ++itf){
		mz += data->layer_res_k[itf];
	}
	
	
	dx = Lx / ((PetscReal)mx);
	dy = Ly / ((PetscReal)my);
	dz = Lz / ((PetscReal)mz);
	
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, 0.0,Lx,0.0 ,Ly, data->interface_heights[0],Lz);CHKERRQ(ierr);
	factor = 0.1;
	ierr = PetscOptionsGetReal(NULL,"-model_WrenchFold_amp_factor",&factor,NULL);CHKERRQ(ierr);
	amp = factor * 1.0; /* this is internal scaled by dy inside WrenchFoldSetPerturbedInterfaces() */
	if ( (amp < 0.0) || (amp >1.0) ) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-model_WrenchFold_amp_factor must be 0 < amp < 1");
	}
	
	/* step 2 - define two interfaces and perturb coords along the interface */
	ierr = WrenchFoldSetPerturbedInterfaces(c->stokes_ctx->dav, data->interface_heights, data->layer_res_k, data->n_interfaces,amp);CHKERRQ(ierr);
	
	ierr = DMDABilinearizeQ2Elements(c->stokes_ctx->dav);CHKERRQ(ierr);

	/*
	ierr = DMGetCoordinates(c->stokes_ctx->dav,&coord);CHKERRQ(ierr);
	ierr = VecScale(coord,1.0e-2);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(c->stokes_ctx->dav,&coord);CHKERRQ(ierr);
	ierr = VecScale(coord,1.0e-2);CHKERRQ(ierr);
*/
	
	PetscFunctionReturn(0);
}

/*
 
 0/ Full lagrangian update
 1/ Check mesh quality metrics
 2/ If mesh quality metrics are not satisfied (on the first failure only)
 a) set projection type to Q1
 
 3/ set advection vel = 0
 4/ remesh
 5/ Check mesh quality metrics
 
 
 */
#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_WrenchFold"
PetscErrorCode ModelApplyUpdateMeshGeometry_WrenchFold(pTatinCtx c,Vec X,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	PetscReal      step;
	PhysCompStokes stokes;
	DM             stokes_pack,dav,dap;
	Vec            velocity,pressure;
	PetscInt       M,N,P;
	PetscInt           metric_L = 5; 
	MeshQualityMeasure metric_list[] = { MESH_QUALITY_ASPECT_RATIO, MESH_QUALITY_DISTORTION, MESH_QUALITY_DIAGONAL_RATIO, MESH_QUALITY_VERTEX_ANGLE, MESH_QUALITY_FACE_AREA_RATIO };
	PetscReal          value[100];
	PetscBool          remesh;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	/* check mesh quality */
	ierr = DMDAComputeMeshQualityMetricList(dav,metric_L,metric_list,value);CHKERRQ(ierr);
	remesh = PETSC_FALSE;
	if (value[0] > 2.0) {
		remesh = PETSC_TRUE;
	}
	if ( (value[1] < 0.7) || (value[1] > 1.0)) {
		remesh = PETSC_TRUE;
	}
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_ASPECT_RATIO\"    %1.4e \n", value[0]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_DISTORTION\"      %1.4e \n", value[1]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_DIAGONAL_RATIO\"  %1.4e \n", value[2]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_VERTEX_ANGLE\"    %1.4e \n", value[3]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_FACE_AREA_RATIO\" %1.4e \n", value[4]);
	
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]] Remeshing currently deactivated \n", __FUNCT__);
	
#if 0
	/* activate marker interpolation */	
	if (remesh) {
		c->coefficient_projection_type = 1;
		
		ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
		ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,N);CHKERRQ(ierr);
	}
#endif	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialCondition_WrenchFold"
PetscErrorCode ModelInitialCondition_WrenchFold(pTatinCtx c,Vec X,void *ctx)
{
	/*
	 ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	 DM stokes_pack,dau,dap;
	 Vec velocity,pressure;
	 PetscReal rho0;
	 DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	 DMDAVecTraverse3d_InterpCtx IntpCtx;*/
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	/*
	 
	 stokes_pack = c->stokes_ctx->stokes_pack;
	 
	 ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	 ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	 
	 ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
	 ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
	 ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	 */
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_WrenchFold"
PetscErrorCode ModelOutput_WrenchFold(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	//char           name[256];
	DataBucket     materialpoint_db;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	
	{
		const int                   nf = 2;
		const MaterialPointVariable mp_prop_list[] = { MPV_viscosity, MPV_density }; 
		
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
		//sprintf(name,"%s_mpoints_cell",prefix);
		//ierr = pTatinOutputParaViewMarkerFields(c->stokes_ctx->stokes_pack,materialpoint_db,nf,mp_prop_list,c->outputpath,name);CHKERRQ(ierr);
		ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_WrenchFold"
PetscErrorCode ModelDestroy_WrenchFold(pTatinCtx c,void *ctx)
{
	ModelWrenchFoldCtx *data = (ModelWrenchFoldCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_WrenchFold"
PetscErrorCode pTatinModelRegister_WrenchFold(void)
{
	ModelWrenchFoldCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelWrenchFoldCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelWrenchFoldCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"wrench_fold");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_WrenchFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_WrenchFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_WrenchFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_WrenchFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_WrenchFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_WrenchFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_WrenchFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_WrenchFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_WrenchFold);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

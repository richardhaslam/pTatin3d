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

#include "model_fault_fold_plastic_ctx.h"
#include "model_utils.h"
#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_FaultFoldPlastic"
PetscErrorCode ModelInitialize_FaultFoldPlastic(pTatinCtx c,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	PetscInt n_int,n;
    RheologyConstants   *rheology;
	DataBucket          materialconstants = c->material_constants;
	PetscBool flg;
    const PetscReal one_third = 0.333333333333333;

	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* assign defaults */
	data->max_layers = 100;
    
	data->max_fault_layers = 1;
    
	data->n_interfaces = 2;
    
    rheology = &c->rheology_constants;
	rheology->rheology_type = RHEOLOGY_VP_STD;
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global = 1.0e+25;
	rheology->eta_lower_cutoff_global = 1.0e-20;
    
    
	PetscOptionsGetInt(PETSC_NULL,"-model_fault_fold_plastic_n_interfaces",&data->n_interfaces,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the number of interfaces including the top and bottom boundaries (-model_fault_fold_plastic_n_interfaces)");
	}

	pTatinModelGetOptionReal("-model_fault_fold_plastic_Lx", &data->Lx, "User must provide the length along the x direction", PETSC_NULL,PETSC_TRUE);
	/*PetscOptionsGetReal(PETSC_NULL,"-model_fault_fold_plastic_Lx",&data->Lx,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the x direction (-model_fault_fold_plastic_Lx)");
	}*/

	pTatinModelGetOptionReal("-model_fault_fold_plastic_Ly", &data->Ly, "User must provide the length along the y direction", PETSC_NULL,PETSC_TRUE);

        PetscOptionsGetReal(PETSC_NULL,"-model_fault_fold_plastic_sigma",&data->sigma,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the spreading of the weak zone (-model_fault_fold_plastic_sigma)");
	}

	n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_interface_heights",data->interface_heights,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide interface heights relative from the base of the model including the top and bottom boundaries (-model_fault_fold_plastic_interface_heights)");
	}
	if (n_int != data->n_interfaces) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d interface heights relative from the base of the model including the top and bottom boundaries (-model_fault_fold_plastic_interface_heights)",data->n_interfaces);
	}
	/*Set the resolutions: */
        /*In z*/
	n_int = data->max_layers;
	PetscOptionsGetIntArray(PETSC_NULL,"-model_fault_fold_plastic_layer_res_k",data->layer_res_k,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer resolution list (-model_fault_fold_plastic_layer_res_k)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_fault_fold_plastic_layer_res_k)",data->n_interfaces-1);
	}
        /*In y (with the fault res.)*/
    data->domain_res_j;
    n_int = 3;
	PetscOptionsGetIntArray(PETSC_NULL,"-model_fault_fold_plastic_domain_res_j",data->domain_res_j,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide domain resolution list in the j direction: sinistral/fault/dextral resolution (-model_fault_fold_plastic_domain_res_j)");
	}
	if (n_int != 3) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide domain resolution list in the j direction: sinistral/fault/dextral resolution  (-model_fault_fold_plastic_domain_res_j)",3);
	}
    
	/* ---------------------------Layer parameters---------------------------------*/
	n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_layer_eta",data->eta,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list (-model_fault_fold_plastic_layer_eta)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity (-model_fault_fold_plastic_layer_eta)",data->n_interfaces-1);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_layer_C0",data->C0,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer cohesion list (-model_fault_fold_plastic_layer_C0)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer cohesion list (-model_fault_fold_plastic_layer_C0)",data->n_interfaces-1);
	}

    n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_layer_mu",data->mu,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer friction list (-model_fault_fold_plastic_layer_mu)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer friction list (-model_fault_fold_plastic_layer_mu)",data->n_interfaces-1);
	}
    
    n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_layer_C0_inf",data->C0_inf,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer cohesion inf value list (-model_fault_fold_plastic_layer_C0_inf)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer cohesion inf value list (-model_fault_fold_plastic_layer_C0_inf)",data->n_interfaces-1);
	} 
    
    n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_layer_mu_inf",data->mu_inf,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer friction inf value list (-model_fault_fold_plastic_layer_mu_inf)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer friction inf value list (-model_fault_fold_plastic_layer_mu_inf)",data->n_interfaces-1);
	}
    
	n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_layer_rho",data->rho,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer density list (-model_fault_fold_plastic_layer_rho)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer density (-model_fault_fold_plastic_layer_rho)",data->n_interfaces-1);
	}
	/*--------------------Central fault parameters-------------------------*/
    n_int = data->max_fault_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_fault_eta",data->fault_eta,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault viscosity list (-model_fault_fold_plastic_fault_eta)");
	}
	if (n_int != data->max_fault_layers) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d fault viscosity (-model_fault_fold_plastic_fault_eta)",data->max_fault_layers);
	}
	
	n_int = data->max_fault_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_fault_C0",data->fault_C0,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault cohesion list (-model_fault_fold_plastic_fault_C0)");
	}
	if (n_int != data->max_fault_layers) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault cohesion list (-model_fault_fold_plastic_fault_C0)",data->max_fault_layers);
	}
    
    n_int = data->max_fault_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_fault_mu",data->fault_mu,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault friction list (-model_fault_fold_plastic_fault_mu)");
	}
	if (n_int != data->max_fault_layers) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault friction list (-model_fault_fold_plastic_fault_mu)",data->max_fault_layers);
	}
    
    n_int = data->max_fault_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_fault_C0_inf",data->fault_C0_inf,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault cohesion inf value list (-model_fault_fold_plastic_fault_C0_inf)");
	}
	if (n_int != data->max_fault_layers) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault cohesion inf value list (-model_fault_fold_plastic_fault_C0_inf)",data->max_fault_layers);
	} 
    
    n_int = data->max_fault_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_fault_mu_inf",data->fault_mu_inf,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault friction inf value list (-model_fault_fold_plastic_fault_mu_inf)");
	}
	if (n_int != data->max_fault_layers) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault friction inf value list (-model_fault_fold_plastic_fault_mu_inf)",data->max_fault_layers);
	}
    
	n_int = data->max_fault_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_fault_rho",data->fault_rho,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide fault density list (-model_fault_fold_plastic_fault_rho)");
	}
	if (n_int != data->max_fault_layers) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d fault density (-model_fault_fold_plastic_fault_rho)",data->max_fault_layers);
	}
	/*-------------------------Set the rheologic parameters------------------------*/
    /* Material constant */
	MaterialConstantsSetDefaults(materialconstants);
    /*For the layer*/
    for(n=0; n<data->n_interfaces-1; n++){
    	MaterialConstantsSetValues_MaterialType(materialconstants,  n, VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
        MaterialConstantsSetValues_ViscosityConst(materialconstants,n, data->eta[n]);
        MaterialConstantsSetValues_DensityConst(materialconstants,  n, data->rho[n]);
        MaterialConstantsSetValues_PlasticDP(materialconstants,     n, data->mu[n],data->mu_inf[n], data->C0[n],data->C0_inf[n], one_third*data->C0[n],1.0e20);
        MaterialConstantsSetValues_SoftLin(materialconstants,       n, 0.5,1.0);
    
    }
    /*For the fault*/
    n=0;
    MaterialConstantsSetValues_MaterialType(materialconstants,  data->n_interfaces, VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
    MaterialConstantsSetValues_ViscosityConst(materialconstants,data->n_interfaces, data->fault_eta[n]);
    MaterialConstantsSetValues_DensityConst(materialconstants,  data->n_interfaces, data->fault_rho[n]);
    MaterialConstantsSetValues_PlasticDP(materialconstants,     data->n_interfaces, data->fault_mu[n],data->fault_mu_inf[n], data->fault_C0[n],data->fault_C0_inf[n], one_third*data->fault_C0[n],1.0e20);
    MaterialConstantsSetValues_SoftLin(materialconstants,       data->n_interfaces, 0.5,1.0);
    
    rheology->nphases_active = data->n_interfaces+1;
    
    /*-------------------Perturbation parameters--------------------*/
    n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_fault_fold_plastic_fold_centers",data->fold_centers,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the initial fold separation a,b,c,d (-model_fault_fold_plastic_fold_centers). The perturbation is at  the position (a/b).Ly at the front and (c/d).Ly in the back. ");
	}
	if (n_int != 4) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_fault_fold_plastic_fold_centers)",data->n_interfaces-1);
	}
    
	/* define the mesh size the z-direction for the global problem */
	c->mz = 0;
	for (n=0; n<data->n_interfaces-1; n++) {
		c->mz += data->layer_res_k[n];
	}
	c->my = 0;
	for (n=0; n<3; n++) {
		c->my += data->domain_res_j[n];
	}
	/* define the domain size in the z-direction for the global problem */
	data->Lz = data->interface_heights[data->n_interfaces-1];
		
    data->bc_type = 0; /* 0 use vx compression ; 1 use exx compression */
	data->exx             = -1.0e-3;
	data->vx_commpression = 1.0;
	
	/* parse from command line or input file */
	ierr = PetscOptionsGetInt(PETSC_NULL,"-model_fault_fold_plastic_bc_type",&data->bc_type,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_fault_fold_plastic_exx",&data->exx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_fault_fold_plastic_vx_commpression",&data->exx,&flg);CHKERRQ(ierr);
    
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
#define __FUNCT__ "BoundaryCondition_FaultFoldPlastic"
PetscErrorCode BoundaryCondition_FaultFoldPlastic(DM dav,BCList bclist,pTatinCtx c,ModelFaultFoldPlasticCtx *data)
{
	PetscReal         exx, zero = 0.0, vx_E=-data->vx_commpression, vx_W = data->vx_commpression;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	exx = data->exx;
    
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    
	if (data->bc_type == 0) {
		/* compression east/west in the x-direction (0) [east-west] using constant velocity */
        
		/*ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,EAST_FACE,-data->vx_commpression);CHKERRQ(ierr);
		ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,WEST_FACE, data->vx_commpression);CHKERRQ(ierr);*/

        
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&vx_W);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&vx_E);CHKERRQ(ierr);
	} else if (data->bc_type == 1) {
		/* compression east/west in the x-direction (0) [east-west] using constant strain rate */
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"BC not Implemented yet!");

		//ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,exx,0);CHKERRQ(ierr);
	} else {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unknonwn boundary condition type");
	}
	
	/* free slip south (base) */
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr); 
	//ierr = DirichletBC_FreeSlip(bclist,dav,SOUTH_FACE);CHKERRQ(ierr);
	
	/* free surface north */
	/* do nothing! */
	
	/* free slip front/back to mimic 2d behaviour */
	//ierr = DirichletBC_FreeSlip(bclist,dav,FRONT_FACE);CHKERRQ(ierr);
	//ierr = DirichletBC_FreeSlip(bclist,dav,BACK_FACE);CHKERRQ(ierr);
    
    /* free slip lateral */
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_FaultFoldPlastic"
PetscErrorCode ModelApplyBoundaryCondition_FaultFoldPlastic(pTatinCtx c,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	PetscReal         exx;
	BCList            bclist;
	DM                dav;
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	
	exx = data->exx;

	bclist = c->stokes_ctx->u_bclist;
	dav    = c->stokes_ctx->dav;
	ierr = BoundaryCondition_FaultFoldPlastic(dav,bclist,c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_FaultFoldPlastic"
PetscErrorCode ModelApplyBoundaryConditionMG_FaultFoldPlastic(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		ierr = BoundaryCondition_FaultFoldPlastic(dav[n],bclist[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_FaultFoldPlastic"
PetscErrorCode ModelApplyMaterialBoundaryCondition_FaultFoldPlastic(pTatinCtx c,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FaultFoldPlasticSetMeshGeometry"
PetscErrorCode FaultFoldPlasticSetMeshGeometry(DM dav, void *ctx)
{
	PetscErrorCode ierr;
    ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz,M,N,P, kinter_max, kinter_min, interf, n_interfaces,*layer_res_k;
    PetscScalar dz, *interface_heights;
    
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	n_interfaces = data->n_interfaces;
    layer_res_k = data->layer_res_k;
    interface_heights = data->interface_heights;
    
    
    kinter_max = 0;
	for(interf = 0; interf < n_interfaces-1; interf++){ 
        kinter_min = kinter_max;
        kinter_max += 2*layer_res_k[interf];
            dz = (interface_heights[interf+1] - interface_heights[interf])/(PetscReal)(2.0*layer_res_k[interf]);
            for(i=si; i<si+nx; i++){
                for(j=sj; j<sj+ny; j++){
                    PetscScalar h;
                    h = data->interface_heights[interf];
                    for(k=sk;k<sk+nz;k++){
                        if((k <= kinter_max)&&(k >= kinter_min)){
                            LA_coord[k][j][i].z = h + (PetscReal)dz*(k-kinter_min); 
                            
                        }   
                    }
                }
            }
        
    }
    
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FaultFoldPlasticSetPerturbedInterfaces"
PetscErrorCode FaultFoldPlasticSetPerturbedInterfaces(DM dav,void *ctx)// PetscScalar interface_heights[], PetscInt layer_res_k[], PetscInt n_interfaces,PetscReal amp, PetscReal Lx)
{
	PetscErrorCode ierr;
    ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
    PetscReal *interface_heights, Lx, amp, *fold_centers; 
    PetscInt *layer_res_k, n_interfaces;
    
	PetscInt i,j,si,sj,sk,nx,ny,nz,M,N,P, interf, kinter;
	PetscScalar dz, dy;
    PetscReal fold_center_front, fold_center_back;
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	
	PetscFunctionBegin;
    
    interface_heights = data->interface_heights;
    layer_res_k = data->layer_res_k;
    n_interfaces = data->n_interfaces;
    amp = data->amp;
    Lx = data->Lx;
    fold_centers = data->fold_centers;
    fold_center_front = Lx*fold_centers[0]/fold_centers[1];
    fold_center_back = Lx*fold_centers[2]/fold_centers[3];
    
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);

	dy = data->Ly/(PetscScalar)(N/2-1);
    
	/*Perturbes the interface for cylindrical folding*/
    kinter = 0;
	for(interf = 1; interf < n_interfaces-1; interf++){
		kinter += 2*layer_res_k[interf-1];
		PetscPrintf(PETSC_COMM_WORLD,"jinter = %d (max=%d)\n", kinter,N-1 );

		if ( (kinter>=sk) && (kinter<sk+nz) ) {
			
			dz = 0.5*((interface_heights[interf+1] - interface_heights[interf])/(PetscScalar)(layer_res_k[interf]) + (interface_heights[interf] - interface_heights[interf-1])/(PetscScalar)(layer_res_k[interf-1]) );
			PetscPrintf(PETSC_COMM_SELF," interface %d: using dy computed from avg %1.4e->%1.4e / my=%d :: %1.4e->%1.4e / my=%d \n", interf,
									interface_heights[interf+1],interface_heights[interf],layer_res_k[interf],
									interface_heights[interf],interface_heights[interf-1],layer_res_k[interf-1] );
			for(i = si; i<si+nx; i++) {

				if((sj+ny == N) && (sj == 0)){
                    PetscScalar center = 0;
                    j=sj+ny-1;
                    center = LA_coord[kinter][j][i].x-fold_center_back;
					LA_coord[kinter][j][i].z += amp * dz * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j-1][i].z += 0.75*amp * dz * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j-2][i].z += 0.5*amp * dz * exp(-center*center/(6.*dy*dy));
                    j=0;
                    center = LA_coord[kinter][j][i].x-fold_center_front;
                    LA_coord[kinter][j][i].z += amp * dz * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j+1][i].z += 0.75*amp * dz * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j+2][i].z += 0.5*amp * dz * exp(-center*center/(6.*dy*dy));
				}else if ((sj+ny == N) || (sj == 0)){
                    PetscScalar center = 0;
                    PetscInt sgn = 1;
                    j= (sj == 0)?0:(sj+ny-1);
                    sgn = (sj == 0)?1:-1;
                    
                    center = (sj == 0)?(LA_coord[kinter][j][i].x-fold_center_front):(LA_coord[kinter][j][i].x-fold_center_back);
					LA_coord[kinter][j][i].z += amp * dz * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j+sgn*1][i].z += 0.75*amp * dz * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j+sgn*2][i].z += 0.5*amp * dz * exp(-center*center/(6.*dy*dy));
				}
			}
			
		}
	}
	
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryMaterialPoints_FaultFoldPlastic"
PetscErrorCode InitialMaterialGeometryMaterialPoints_FaultFoldPlastic(pTatinCtx c,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	int                    p,n_mp_points;
    MPAccess       mpX;
	DataBucket             db, material_points;
	DataField              PField_std,PField_stokes;
	PetscErrorCode ierr;
			
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
			
	
    /**/
    
    ierr = pTatinGetMaterialPoints(c,&material_points,PETSC_NULL);CHKERRQ(ierr);
	DataBucketGetSizes(material_points,&n_mp_points,0,0);
	
	ierr = MaterialPointGetAccess(material_points,&mpX);CHKERRQ(ierr);
    
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		double      *position;
		PetscInt    phase;
		PetscInt    layer, kmaxlayer, kminlayer;
		PetscInt    I, J, K;

		DataFieldAccessPoint(mpX->PField[MPField_Std],p,   (void**)&material_point);
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		//MPntStdGetField_global_coord(material_point,&position);
        MPntGetField_global_element_IJKindex(c->stokes_ctx->dav,material_point, &I, &J, &K);
		phase = -1;
		kmaxlayer = kminlayer = 0;
		layer = 0;
		// gets the global element index (i,j,k)
		//....
		
		/*Set the phases*/
		while( (phase == -1) && (layer < data->n_interfaces-1) ){
			kmaxlayer += data->layer_res_k[layer];
			
			if( (K<kmaxlayer) && (K>=kminlayer) ){
				phase = layer + 1;
			}
			kminlayer += data->layer_res_k[layer];
			layer++;
		}
        /*Check if we are in the fault:*/
        if ((data->domain_res_j[0]<=J) && (data->domain_res_j[0]+data->domain_res_j[1]>J)){
            phase = data->n_interfaces;
        } 
		/* user the setters provided for you */
        ierr = MaterialPointSet_phase_index(mpX,p,phase);CHKERRQ(ierr);
	}
	ierr = MaterialPointRestoreAccess(material_points,&mpX);CHKERRQ(ierr);			

			
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryQuadraturePoints_FaultFoldPlastic"
PetscErrorCode InitialMaterialGeometryQuadraturePoints_FaultFoldPlastic(pTatinCtx c,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
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
		double      *position;
		PetscReal      eta,rho, center2, sigma2;
		PetscInt    phase;
		PetscInt    layer, kmaxlayer, kminlayer, localeid_p;
		PetscInt    I, J, K;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		MPntStdGetField_global_coord(material_point,&position);
        
		sigma2  = data->sigma*data->sigma;
        center2 = (position[1]-0.5*data->Ly)*(position[1]-0.5*data->Ly);        
        
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
				//eta = data->eta[layer]- (data->eta[layer]-data->etaweak[layer])*exp(-center2/(2.0*sigma2));
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
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_FaultFoldPlastic"
PetscErrorCode ModelApplyInitialMaterialGeometry_FaultFoldPlastic(pTatinCtx c,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	//int                    p,n_mp_points;
	//DataBucket             db;
	//DataField              PField_std,PField_stokes;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = InitialMaterialGeometryMaterialPoints_FaultFoldPlastic(c,data);CHKERRQ(ierr);
	ierr = InitialMaterialGeometryQuadraturePoints_FaultFoldPlastic(c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

		
		
		
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_FaultFoldPlastic"
PetscErrorCode ModelApplyInitialMeshGeometry_FaultFoldPlastic(pTatinCtx c,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	PetscReal         Lx,Ly,Lz;
	PetscInt          itf;
	PetscReal         amp,factor;
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* step 1 - create structured grid */
	Lx = data->Lx;
	Ly = data->Ly;
	Lz = data->Lz;


	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, 0.0,Lx, 0.0,Ly, data->interface_heights[0], Lz);CHKERRQ(ierr);
	factor = 0.1;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_FaultFoldPlastic_amp_factor",&factor,PETSC_NULL);CHKERRQ(ierr);
	amp = factor * 1.0; /* this is internal scaled by dy inside FaultFoldPlasticSetPerturbedInterfaces() */
    data->amp = amp;
	if ( (amp < 0.0) || (amp >1.0) ) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-model_FaultFoldPlastic_amp_factor must be 0 < amp < 1");
	}
	ierr = FaultFoldPlasticSetMeshGeometry(c->stokes_ctx->dav, data);CHKERRQ(ierr);
	/* step 2 - define two interfaces and perturb coords along the interface */
	ierr = FaultFoldPlasticSetPerturbedInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);//data->interface_heights, data->layer_res_k, data->n_interfaces,amp, Lx);CHKERRQ(ierr);
	
	ierr = DMDABilinearizeQ2Elements(c->stokes_ctx->dav);CHKERRQ(ierr);
    

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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_FaultFoldPlastic"
PetscErrorCode ModelApplyUpdateMeshGeometry_FaultFoldPlastic(pTatinCtx c,Vec X,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	PetscReal      step;
	PhysCompStokes stokes;
	DM             stokes_pack,dav,dap;
	Vec            velocity,pressure;
	//PetscInt       M,N,P;
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
#define __FUNCT__ "ModelInitialCondition_FaultFoldPlastic"
PetscErrorCode ModelInitialCondition_FaultFoldPlastic(pTatinCtx c,Vec X,void *ctx)
{
    /*
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
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
#define __FUNCT__ "ModelOutput_FaultFoldPlastic"
PetscErrorCode ModelOutput_FaultFoldPlastic(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	//char           name[256];
	DataBucket     materialpoint_db;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	
	{
		const int                   nf = 2;
		const MaterialPointVariable mp_prop_list[] = { MPV_viscosity, MPV_density }; 
		
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
		//sprintf(name,"%s_mpoints_cell",prefix);
		//ierr = pTatinOutputParaViewMarkerFields(c->stokes_ctx->stokes_pack,materialpoint_db,nf,mp_prop_list,c->outputpath,name);CHKERRQ(ierr);
		ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_FaultFoldPlastic"
PetscErrorCode ModelDestroy_FaultFoldPlastic(pTatinCtx c,void *ctx)
{
	ModelFaultFoldPlasticCtx *data = (ModelFaultFoldPlasticCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_FaultFoldPlastic"
PetscErrorCode pTatinModelRegister_FaultFoldPlastic(void)
{
	ModelFaultFoldPlasticCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelFaultFoldPlasticCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelFaultFoldPlasticCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"fault_fold_plastic");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_FaultFoldPlastic);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_FaultFoldPlastic);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_FaultFoldPlastic);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_FaultFoldPlastic);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_FaultFoldPlastic);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_FaultFoldPlastic);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_FaultFoldPlastic);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_FaultFoldPlastic);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_FaultFoldPlastic);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

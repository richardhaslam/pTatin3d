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
#include "dmda_redundant.h"
#include "model_fault_fold_ctx.h"
#include "model_utils.h"
#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_FaultFold"
PetscErrorCode ModelInitialize_FaultFold(pTatinCtx c,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	PetscInt n_int,n;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* assign defaults */
	data->max_layers = 100;
	
	data->n_interfaces = 2;
	PetscOptionsGetInt(NULL,"-model_fault_fold_n_interfaces",&data->n_interfaces,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the number of interfaces including the top and bottom boundaries (-model_fault_fold_n_interfaces)");
	}

	pTatinModelGetOptionReal("-model_fault_fold_Lx", &data->Lx, "User must provide the length along the x direction", NULL,PETSC_TRUE);
	/*PetscOptionsGetReal(NULL,"-model_fault_fold_Lx",&data->Lx,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the x direction (-model_fault_fold_Lx)");
	}*/

	pTatinModelGetOptionReal("-model_fault_fold_Ly", &data->Ly, "User must provide the length along the y direction", NULL,PETSC_TRUE);

        PetscOptionsGetReal(NULL,"-model_fault_fold_sigma",&data->sigma,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the spreading of the weak zone (-model_fault_fold_sigma)");
	}

	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_fault_fold_interface_heights",data->interface_heights,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide interface heights relative from the base of the model including the top and bottom boundaries (-model_fault_fold_interface_heights)");
	}
	if (n_int != data->n_interfaces) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d interface heights relative from the base of the model including the top and bottom boundaries (-model_fault_fold_interface_heights)",data->n_interfaces);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetIntArray(NULL,"-model_fault_fold_layer_res_k",data->layer_res_k,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer resolution list (-model_fault_fold_layer_res_k)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_fault_fold_layer_res_k)",data->n_interfaces-1);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_fault_fold_layer_eta",data->eta,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list (-model_fault_fold_layer_eta)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity (-model_fault_fold_layer_eta)",data->n_interfaces-1);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_fault_fold_layer_etaweak",data->etaweak,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list for the weak zone(-model_fault_fold_layer_etaweak)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity for the weak zone (-model_fault_fold_layer_etaweak)",data->n_interfaces-1);
	}

	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_fault_fold_layer_rho",data->rho,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer density list (-model_fault_fold_layer_rho)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer density (-model_fault_fold_layer_rho)",data->n_interfaces-1);
	}
	
    n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_fault_fold_fold_centers",data->fold_centers,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the initial fold separation a,b,c,d (-model_fault_fold_fold_centers). The perturbation is at  the position (a/b).Ly at the front and (c/d).Ly in the back. ");
	}
	if (n_int != 4) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_fault_fold_fold_centers)",data->n_interfaces-1);
	}
    data->amp = 0.05;
	ierr = PetscOptionsGetReal(NULL,"-model_fault_fold_amp",&data->amp,NULL);CHKERRQ(ierr);


	/* define the mesh size the z-direction for the global problem */
	c->mz = 0;
	for (n=0; n<data->n_interfaces-1; n++) {
		c->mz += data->layer_res_k[n];
	}

	/* define the domain size in the z-direction for the global problem */
	data->Lz = data->interface_heights[data->n_interfaces-1];
		
    data->bc_type = 0; /* 0 use vx compression ; 1 use exx compression */
	data->exx             = -1.0e-3;
	data->vx_commpression = 1.0;
	data->perturbation_type = 0;
	data->perturbation_width = (PetscInt)((2*c->my+1)/4);

	/* parse from command line or input file */
    ierr = PetscOptionsGetInt(NULL,"-model_fault_fold_perturbation_type",&data->perturbation_type,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,"-model_fault_fold_perturbation_width",&data->perturbation_width,&flg);CHKERRQ(ierr);
    if (data->perturbation_width > (PetscInt)((2*c->my+1)/2) ){
        data->perturbation_width = (PetscInt)((2*c->my+1)/2);
    }
    data->phi = 0;
    ierr = PetscOptionsGetReal(NULL,"-model_fault_fold_phi",&data->phi,&flg);CHKERRQ(ierr);
    data->phi *= 2.0*M_PI;
    data->offset = 1.0;
    ierr = PetscOptionsGetReal(NULL,"-model_fault_fold_offset",&data->offset,&flg);CHKERRQ(ierr);
    data->psig1 = 0;
    ierr = PetscOptionsGetReal(NULL,"-model_fault_fold_psig1",&data->psig1,&flg);CHKERRQ(ierr);
    data->psig2 = 0;
    ierr = PetscOptionsGetReal(NULL,"-model_fault_fold_psig2",&data->psig2,&flg);CHKERRQ(ierr);
	/* parse from command line or input file */
	ierr = PetscOptionsGetInt(NULL,"-model_fault_fold_bc_type",&data->bc_type,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_fault_fold_exx",&data->exx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_fault_fold_vx_commpression",&data->exx,&flg);CHKERRQ(ierr);
    
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
#define __FUNCT__ "BoundaryCondition_FaultFold"
PetscErrorCode BoundaryCondition_FaultFold(DM dav,BCList bclist,pTatinCtx c,ModelFaultFoldCtx *data)
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
#define __FUNCT__ "ModelApplyBoundaryCondition_FaultFold"
PetscErrorCode ModelApplyBoundaryCondition_FaultFold(pTatinCtx c,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	PetscReal         exx;
	BCList            bclist;
	DM                dav;
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	
	exx = data->exx;

	bclist = c->stokes_ctx->u_bclist;
	dav    = c->stokes_ctx->dav;
	ierr = BoundaryCondition_FaultFold(dav,bclist,c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_FaultFold"
PetscErrorCode ModelApplyBoundaryConditionMG_FaultFold(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		ierr = BoundaryCondition_FaultFold(dav[n],bclist[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_FaultFold"
PetscErrorCode ModelApplyMaterialBoundaryCondition_FaultFold(pTatinCtx c,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FaultFoldSetMeshGeometry"
PetscErrorCode FaultFoldSetMeshGeometry(DM dav, void *ctx)
{
	PetscErrorCode ierr;
    ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz,M,N,P, kinter_max, kinter_min, interf, n_interfaces,*layer_res_k;
    PetscScalar dz, *interface_heights;
    
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
#define __FUNCT__ "FaultFoldRemeshingAccordingToTheInterfaces"
PetscErrorCode FaultFoldRemeshingAccordingToTheInterfaces(DM dav, void *ctx)
{   
    
    ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz, rank, kinter_max, kinter_min, interf;
	PetscScalar dz;
    PetscInt *layer_res_k, n_interfaces;
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
    PetscErrorCode ierr;
    
 	PetscFunctionBegin;

    layer_res_k = data->layer_res_k;
    n_interfaces = data->n_interfaces;
    
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
    
    /* Loop again through the layers  to set the perturbation around the layer.*/    
    kinter_max = 0;
	for(interf = 0; interf < n_interfaces-1; interf++){ 
        DM botinterface_da, topinterface_da;
        Vec botinterface_coords, topinterface_coords;
        PetscScalar *botinterface_nodes, *topinterface_nodes;
        
        kinter_min = kinter_max;
        kinter_max += 2*layer_res_k[interf];
        ierr = DMDACreate3dRedundant( dav, si,si+nx, sj,sj+ny, kinter_min, kinter_min + 1, 1, &botinterface_da );CHKERRQ(ierr);
        ierr = DMDACreate3dRedundant( dav, si,si+nx, sj,sj+ny, kinter_max, kinter_max+1, 1, &topinterface_da );CHKERRQ(ierr);
        
        ierr = DMGetCoordinates( botinterface_da,&botinterface_coords );CHKERRQ(ierr);
        ierr = VecGetArray(botinterface_coords,&botinterface_nodes);CHKERRQ(ierr);
        
        ierr = DMGetCoordinates( topinterface_da,&topinterface_coords );CHKERRQ(ierr);
        ierr = VecGetArray(topinterface_coords,&topinterface_nodes);CHKERRQ(ierr);
        
        for(i=0; i<nx; i++){
            for(j=0; j<ny; j++){
                for(k=sk;k<sk+nz;k++){
                    if((k <= kinter_max)&&(k >= kinter_min)){
                        dz = (  topinterface_nodes[3*(0+nx*j+i)+2]  -  botinterface_nodes[3*(0+nx*j+i)+2]  )/((PetscReal)(2.0*layer_res_k[interf]));
                        LA_coord[k][j+sj][i+si].z = botinterface_nodes[3*(0+nx*j+i)+2] + (PetscReal)dz*(k-kinter_min); 
                    }   
                }
            }
        }
        
        
        ierr = VecRestoreArray(botinterface_coords,&botinterface_nodes);CHKERRQ(ierr);
        ierr = VecRestoreArray(topinterface_coords,&topinterface_nodes);CHKERRQ(ierr);
        
        ierr = DMDestroy(&botinterface_da);CHKERRQ(ierr);
        ierr = DMDestroy(&topinterface_da);CHKERRQ(ierr);
        
    }
    
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);

    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FaultFoldSetPerturbedInterfaces"
PetscErrorCode FaultFoldSetPerturbedInterfaces(DM dav, void *ctx)
{
	PetscErrorCode ierr;
    ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz,M,N,P, interf, kinter, rank, kinter_max, kinter_min;
	PetscScalar dz, dy, pertu, attenuation, H;
    PetscReal *interface_heights, lamb, phi, offset, psig1, psig2;
    PetscInt *layer_res_k, n_interfaces, pwidth;
    PetscReal amp;
    PetscReal fold_center_front, fold_center_back, Lx, Ly, *fold_centers;
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    interface_heights = data->interface_heights;

    layer_res_k = data->layer_res_k;
    n_interfaces = data->n_interfaces;
    fold_centers = data->fold_centers;
    amp = data->amp;
    phi = data->phi;
    psig1 = data->psig1;
    psig2 = data->psig2;
    Lx = data->Lx;
    Ly = data->Ly;
    offset = data->offset;

	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	


    

    kinter = 0;
    MPI_Comm_rank(PetscObjectComm((PetscObject)dav),&rank);
	for(interf = 1; interf < n_interfaces-1; interf++){
		kinter += 2*layer_res_k[interf-1];
		PetscPrintf(PETSC_COMM_WORLD,"jinter = %d (max=%d)\n", kinter,P-1 );
        srand(rank*interf+2);//The seed changes with the interface and the process.
        
		if ( (kinter>=sk) && (kinter<sk+nz) ) {
            /*Take the dominant wavelength of the viscous layer*/
                lamb = 0.0;
                PetscBool flg;
                
                if(data->eta[interf-1] < data->eta[interf]){
                    H = interface_heights[interf+1] - interface_heights[interf];
                }else{
                    H = (interface_heights[interf] - interface_heights[interf-1]);
                }
                
                PetscOptionsGetReal(NULL,"-model_fault_fold_dominantlength",&lamb,&flg);
                if (!flg) {
                    /*Take the dominant wavelength of the viscous layer in the thin plate approximation*/
                    if(data->eta[interf-1] < data->eta[interf]){
                        lamb = 2.0*M_PI*H*pow(data->eta[interf]/(6.0*data->eta[interf-1]), 1.0/3.0);
                        
                    }else{
                        lamb = 2.0*M_PI*H*pow(data->eta[interf-1]/(6.0*data->eta[interf]), 1.0/3.0);
                    }
                }


        switch(data->perturbation_type){ 
            case 0:{
                fold_centers = data->fold_centers;
                fold_center_front = 0.25*Lx + 0.5*Lx*fold_centers[0]/fold_centers[1];
                fold_center_back = 0.25*Lx + 0.5*Lx*fold_centers[2]/fold_centers[3];
                dy = data->Ly/(PetscScalar)(N/2-1);
            for(i = si; i<si+nx; i++) {
                if((sj+ny == N) && (sj == 0)){
                    PetscScalar center = 0;
                    j=sj+ny-1;
                    center = LA_coord[kinter][j][i].x-fold_center_back;
					LA_coord[kinter][j][i].z += amp * H * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j-1][i].z += 0.75*amp * H * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j-2][i].z += 0.5*amp * H * exp(-center*center/(6.*dy*dy));
                    j=0;
                    center = LA_coord[kinter][j][i].x-fold_center_front;
                    LA_coord[kinter][j][i].z += amp * H * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j+1][i].z += 0.75*amp * H * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j+2][i].z += 0.5*amp * H * exp(-center*center/(6.*dy*dy));
				}else if ((sj+ny == N) || (sj == 0)){
                    PetscScalar center = 0;
                    PetscInt sgn = 1;
                    j= (sj == 0)?0:(sj+ny-1);
                    sgn = (sj == 0)?1:-1;
                    
                    center = (sj == 0)?(LA_coord[kinter][j][i].x-fold_center_front):(LA_coord[kinter][j][i].x-fold_center_back);
					LA_coord[kinter][j][i].z += amp * dz * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j+sgn*1][i].z += 0.75*amp * H * exp(-center*center/(6.*dy*dy));
                    LA_coord[kinter][j+sgn*2][i].z += 0.5*amp * H * exp(-center*center/(6.*dy*dy));
				}
			}
            break;
            }
            case 1:{
            attenuation = -log(0.01)/(PetscScalar)(data->perturbation_width);
            PetscReal sigp = 2.0*lamb;    
            for(i = si; i<si+nx; i++) {
                PetscScalar center = 0;
				if((sj+ny == N) && (sj == 0)){
                    j=sj+ny-1;
                    center = (LA_coord[kinter][j][i].x-0.5*Lx);
                    pertu =  amp*H*cos(center/(offset*lamb)+phi)*exp(-0.5*center*center/(2.*sigp*sigp));
                    for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                        LA_coord[kinter][j-pwidth][i].z += pertu*exp(-pwidth*attenuation);
                    }
                    j=0;
                    center = (LA_coord[kinter][j][i].x-0.5*Lx);
                    pertu =  amp*H*cos(center/lamb)*exp(-0.5*center*center/(sigp*sigp));
                    for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                        LA_coord[kinter][j+pwidth][i].z += pertu*exp(-pwidth*attenuation);
                    }

				}else if ((sj+ny == N) || (sj == 0)){
                    PetscReal dz = 0.0;
                    PetscInt sgn = 1;
                    sgn = (sj == 0)?1:-1;
                    j = (sj == 0)?0:(sj+ny-1);
                    center = (LA_coord[kinter][j][i].x-0.5*Lx);
                    pertu = (sj == 0)?amp*H*cos(center/lamb)*exp(-0.5*center*center/(sigp*sigp)):amp*H*cos(center/(offset*lamb)+phi)*exp(-0.5*center*center/(2.*sigp*sigp));
                    for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                        LA_coord[kinter][j+sgn*pwidth][i].z += pertu *exp(-pwidth*attenuation);
                    }

				}
			}
            break;
            }
            case 2:
            { 
                attenuation = -log(0.01)/(PetscScalar)(data->perturbation_width);
                fold_centers = data->fold_centers;
                fold_center_front = 0.25*Lx + 0.5*Lx*fold_centers[0]/fold_centers[1];
                fold_center_back = 0.25*Lx + 0.5*Lx*fold_centers[2]/fold_centers[3];
                dy = data->Ly/(PetscScalar)(N/2-1);
                for(i = si; i<si+nx; i++) {
                    if((sj+ny == N) && (sj == 0)){
                        PetscScalar center = 0;
                        j=sj+ny-1;
                        center = LA_coord[kinter][j][i].x-fold_center_back;
                        pertu = amp * H * exp(-0.5*center*center/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j-pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                        j=0;
                        pertu = amp * H *(exp(-0.5*(center-offset*psig2)*(center-offset*psig2)/(psig2*psig2)) + exp(-0.5*(center+offset*psig2)*(center+offset*psig2)/(psig2*psig2)));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j+pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }else if ((sj+ny == N) || (sj == 0)){
                        PetscScalar center = 0;
                        PetscInt sgn = 1;
                        j= (sj == 0)?0:(sj+ny-1);
                        sgn = (sj == 0)?1:-1;
                        
                        center = (sj == 0)?(LA_coord[kinter][j][i].x-fold_center_front):(LA_coord[kinter][j][i].x-fold_center_back);
                        pertu = (sj == 0)?(amp * H *(exp(-0.5*(center-offset*psig2)*(center-offset*psig2)/(psig2*psig2)) + exp(-0.5*(center+offset*psig2)*(center+offset*psig2)/(psig2*psig2)))):amp * H * exp(-0.5*center*center/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j+sgn*pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }
                }
                break;    
                }
                case 3:
                { //Mainly to see the absorbed embryonic folds, if it exists?
                attenuation = -log(0.01)/(PetscScalar)(data->perturbation_width);
                fold_centers = data->fold_centers;
                fold_center_front = 0.25*Lx + 0.5*Lx*fold_centers[0]/fold_centers[1];
                fold_center_back = 0.25*Lx + 0.5*Lx*fold_centers[2]/fold_centers[3];
                dy = data->Ly/(PetscScalar)(N/2-1);
                for(i = si; i<si+nx; i++) {
                    if(sj+ny == N){
                        PetscScalar center = 0;
                        j=sj+ny-1;
                        center = LA_coord[kinter][j][i].x-fold_center_back;
                        pertu = amp * H * exp(-0.5*center*center/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j-pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }

                    for(j = sj; j<sj+ny; j++){
                        
                        PetscScalar centerr, centerl;
                        centerr = ((LA_coord[kinter][j][i].x-fold_center_front) - offset*psig2)*((LA_coord[kinter][j][i].x-fold_center_front) - offset*psig2) + (0.5*Ly - LA_coord[kinter][j][i].y)*(0.5*Ly - LA_coord[kinter][j][i].y);
                        centerl = ((LA_coord[kinter][j][i].x-fold_center_front) + offset*psig2)*((LA_coord[kinter][j][i].x-fold_center_front) + offset*psig2) + (0.5*Ly - LA_coord[kinter][j][i].y)*(0.5*Ly - LA_coord[kinter][j][i].y);
                        if ( (sj<= (int)(0.5*N)) && (sj+ny >= (int)(0.5*N)) ){
                            LA_coord[kinter][j][i].z += 1.0*amp * H *(exp(-0.5*(centerr)*(centerr)/(psig2*psig2)) + exp(-0.5*(centerl)*(centerl)/(psig2*psig2)))*exp(-absolute(j-(int)(0.5*N))*attenuation);
                        }else{
                            LA_coord[kinter][j][i].z += 1.0*amp * H *(exp(-0.5*(centerr)*(centerr)/(psig2*psig2)) + exp(-0.5*(centerl)*(centerl)/(psig2*psig2)));
                        }
                        
                    }
                        
                        
                }
                break;
                }
                case 4:
                { 
                attenuation = -log(0.01)/(PetscScalar)(data->perturbation_width);
                fold_centers = data->fold_centers;
                fold_center_front = 0.25*Lx + 0.5*Lx*fold_centers[0]/fold_centers[1];
                fold_center_back = 0.25*Lx + 0.5*Lx*fold_centers[2]/fold_centers[3];
                dy = data->Ly/(PetscScalar)(N/2-1);
                psig1 = 0.25*lamb;
                psig2 = psig1;
                
                for(i = si; i<si+nx; i++) {
                    if((sj+ny == N) && (sj == 0)){
                        PetscScalar center = 0;
                        j=sj+ny-1;
                        center = LA_coord[kinter][j][i].x-fold_center_back;
                        pertu = amp * H * exp(-0.5*center*center/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j-pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                        j=0;
                        pertu = amp * H *(exp(-0.5*(center-offset*lamb)*(center-offset*lamb)/(psig2*psig2)) + exp(-0.5*(center+offset*lamb)*(center+offset*lamb)/(psig2*psig2)));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j+pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }else if ((sj+ny == N) || (sj == 0)){
                        PetscScalar center = 0;
                        PetscInt sgn = 1;
                        j= (sj == 0)?0:(sj+ny-1);
                        sgn = (sj == 0)?1:-1;
                        
                        center = (sj == 0)?(LA_coord[kinter][j][i].x-fold_center_front):(LA_coord[kinter][j][i].x-fold_center_back);
                        pertu = (sj == 0)?(amp * H *(exp(-0.5*(center-offset*lamb)*(center-offset*lamb)/(psig2*psig2)) + exp(-0.5*(center+offset*lamb)*(center+offset*lamb)/(psig2*psig2)) )):amp * H * exp(-0.5*center*center/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j+sgn*pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }
                }
                break;
                }
                case 5:
                { 
                attenuation = -log(0.01)/(PetscScalar)(data->perturbation_width);
                fold_centers = data->fold_centers;
                fold_center_front = 0.25*Lx + 0.5*Lx*fold_centers[0]/fold_centers[1];
                fold_center_back = 0.25*Lx + 0.5*Lx*fold_centers[2]/fold_centers[3];
                dy = data->Ly/(PetscScalar)(N/2-1);
                psig1 = 0.25*lamb;
                psig2 = psig1;
                
                for(i = si; i<si+nx; i++) {
                    if((sj+ny == N) && (sj == 0)){
                        PetscScalar center = 0;
                        j=sj+ny-1;
                        center = LA_coord[kinter][j][i].x-fold_center_back;
                        pertu = amp * H * exp(-0.5*(center-offset*0.5*lamb)*(center-offset*0.5*lamb)/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j-pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                        j=0;
                        pertu = amp * H *exp(-0.5*(center+offset*0.5*lamb)*(center+offset*0.5*lamb)/(psig2*psig2));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j+pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }else if ((sj+ny == N) || (sj == 0)){
                        PetscScalar center = 0;
                        PetscInt sgn = 1;
                        j= (sj == 0)?0:(sj+ny-1);
                        sgn = (sj == 0)?1:-1;
                        
                        center = (sj == 0)?(LA_coord[kinter][j][i].x-fold_center_front):(LA_coord[kinter][j][i].x-fold_center_back);
                        pertu = (sj == 0)?(amp * H *exp(-0.5*(center+offset*0.5*lamb)*(center+offset*0.5*lamb)/(psig2*psig2))):amp * H * exp(-0.5*(center-offset*0.5*lamb)*(center-offset*0.5*lamb)/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j+sgn*pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }
                }
                break;
                }
                case 6://not centered perturbation
                { 
                attenuation = -log(0.01)/(PetscScalar)(data->perturbation_width);
                fold_centers = data->fold_centers;
                fold_center_front = 0.25*Lx + 0.5*Lx*fold_centers[0]/fold_centers[1];
                fold_center_back = 0.25*Lx + 0.5*Lx*fold_centers[2]/fold_centers[3];
                dy = data->Ly/(PetscScalar)(N/2-1);
                psig1 = 0.25*lamb;
                psig2 = psig1;
                
                for(i = si; i<si+nx; i++) {
                    if((sj+ny == N) && (sj == 0)){
                        PetscScalar center = 0;
                        j=sj+ny-1;
                        center = LA_coord[kinter][j][i].x-fold_center_back;
                        pertu = amp * H * exp(-0.5*(center)*(center)/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j-pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                        j=0;
                        pertu = amp * H *exp(-0.5*(center+offset*lamb)*(center+offset*lamb)/(psig2*psig2));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j+pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }else if ((sj+ny == N) || (sj == 0)){
                        PetscScalar center = 0;
                        PetscInt sgn = 1;
                        j= (sj == 0)?0:(sj+ny-1);
                        sgn = (sj == 0)?1:-1;
                        
                        center = (sj == 0)?(LA_coord[kinter][j][i].x-fold_center_front):(LA_coord[kinter][j][i].x-fold_center_back);
                        pertu = (sj == 0)?(amp * H *exp(-0.5*(center+offset*lamb)*(center+offset*lamb)/(psig2*psig2))):amp * H * exp(-0.5*(center)*(center)/(psig1*psig1));
                        for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                            LA_coord[kinter][j+sgn*pwidth][i].z += pertu *exp(-pwidth*attenuation);
                        }
                    }
                }
                break;
            }                
            }
                
        
            

	}
    }

	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);

    
	PetscFunctionReturn(0);


}

#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryMaterialPoints_FaultFold"
PetscErrorCode InitialMaterialGeometryMaterialPoints_FaultFold(pTatinCtx c,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
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
		double      *position;
		PetscReal      eta,rho, sigma2, center2;
		PetscInt    phase;
		PetscInt    layer, kmaxlayer, kminlayer;
		PetscInt    I, J, K;

		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
        
		sigma2  = data->sigma*data->sigma;
        center2 = (position[1]-0.5*data->Ly)*(position[1]-0.5*data->Ly);
        
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
				eta = data->eta[layer]- (data->eta[layer]-data->etaweak[layer])*exp(-center2/(2.0*sigma2));
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
#define __FUNCT__ "InitialMaterialGeometryQuadraturePoints_FaultFold"
PetscErrorCode InitialMaterialGeometryQuadraturePoints_FaultFold(pTatinCtx c,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
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
				eta = data->eta[layer]- (data->eta[layer]-data->etaweak[layer])*exp(-center2/(2.0*sigma2));
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
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_FaultFold"
PetscErrorCode ModelApplyInitialMaterialGeometry_FaultFold(pTatinCtx c,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	//int                    p,n_mp_points;
	//DataBucket             db;
	//DataField              PField_std,PField_stokes;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = InitialMaterialGeometryMaterialPoints_FaultFold(c,data);CHKERRQ(ierr);
	ierr = InitialMaterialGeometryQuadraturePoints_FaultFold(c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

		
		
		
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_FaultFold"
PetscErrorCode ModelApplyInitialMeshGeometry_FaultFold(pTatinCtx c,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
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

	ierr = FaultFoldSetMeshGeometry(c->stokes_ctx->dav, data);CHKERRQ(ierr);
	/* step 2 - define two interfaces and perturb coords along the interface */
	ierr = FaultFoldSetPerturbedInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);//data->interface_heights, data->layer_res_k, data->n_interfaces,amp, Lx);CHKERRQ(ierr);
    ierr =  FaultFoldRemeshingAccordingToTheInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);
	
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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_FaultFold"
PetscErrorCode ModelApplyUpdateMeshGeometry_FaultFold(pTatinCtx c,Vec X,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
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
#define __FUNCT__ "ModelInitialCondition_FaultFold"
PetscErrorCode ModelInitialCondition_FaultFold(pTatinCtx c,Vec X,void *ctx)
{
    /*
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
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
#define __FUNCT__ "ModelOutput_FaultFold"
PetscErrorCode ModelOutput_FaultFold(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
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
#define __FUNCT__ "ModelDestroy_FaultFold"
PetscErrorCode ModelDestroy_FaultFold(pTatinCtx c,void *ctx)
{
	ModelFaultFoldCtx *data = (ModelFaultFoldCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_FaultFold"
PetscErrorCode pTatinModelRegister_FaultFold(void)
{
	ModelFaultFoldCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelFaultFoldCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelFaultFoldCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"fault_fold");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_FaultFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_FaultFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_FaultFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_FaultFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_FaultFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_FaultFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_FaultFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_FaultFold);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_FaultFold);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

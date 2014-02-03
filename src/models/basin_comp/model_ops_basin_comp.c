
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
#include "element_utils_q2.h"
#include "ptatin3d_defs.h"
#include "math.h"
#include "dmda_redundant.h"
#include <time.h>

#include "model_basin_comp_ctx.h"
#include "model_utils.h"




#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_BasinComp"
PetscErrorCode ModelInitialize_BasinComp(pTatinCtx c,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscInt n_int,n;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);


	/* assign defaults */
	data->max_layers = 100;
	
	PetscOptionsGetInt(NULL,"-model_basin_comp_n_interfaces",&data->n_interfaces,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the number of interfaces including the top and bottom boundaries (-model_basin_comp_n_interfaces)");
	}
    
	PetscOptionsGetReal(NULL,"-model_basin_comp_Lx",&data->Lx,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the x direction (-model_basin_comp_Lx)");
	}
	
	PetscOptionsGetReal(NULL,"-model_basin_comp_Ly",&data->Ly,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the y direction (-model_basin_comp_Ly)");
	}
    
	n_int = data->n_interfaces;
	PetscOptionsGetRealArray(NULL,"-model_basin_comp_interface_heights_f",data->interface_heights_f,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide front interface heights relative from the base of the model including the top and bottom boundaries. Interface heights taken from the top of the slope (-model_basin_comp_interface_heights_f)");
	}
	if (n_int != data->n_interfaces) {
        	    //printf("------>%d %f   %f    %f\n",n_int, data->interface_heights[0], data->interface_heights[1], data->interface_heights[2]);
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d front interface heights relative from the base of the model including the top and bottom boundaries (-model_basin_comp_interface_heights_f)",data->n_interfaces);

    }
    
	n_int = data->n_interfaces;
	PetscOptionsGetRealArray(NULL,"-model_basin_comp_interface_heights_b",data->interface_heights_b,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide back interface heights relative from the base of the model including the top and bottom boundaries. Interface heights taken from the top of the slope (-model_basin_comp_interface_heights_b)");
	}
	if (n_int != data->n_interfaces) {
        //printf("------>%d %f   %f    %f\n",n_int, data->interface_heights[0], data->interface_heights[1], data->interface_heights[2]);
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d back interface heights relative from the base of the model including the top and bottom boundaries (-model_basin_comp_interface_heights_b)",data->n_interfaces);
        
    }    
	data->Lz = data->interface_heights_f[data->n_interfaces-1];

    PetscOptionsGetIntArray(NULL,"-model_basin_comp_layer_res_k",data->layer_res_k,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer resolution list (-model_basin_comp_layer_res_k)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_basin_comp_layer_res_k)",data->n_interfaces-1);
	}
    
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_basin_comp_layer_eta_f",data->eta_f,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list at the front (-model_basin_comp_layer_eta_f)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity at the front (-model_basin_comp_layer_eta_f)",data->n_interfaces-1);
	}
    
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_basin_comp_layer_eta_b",data->eta_b,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list at the back (-model_basin_comp_layer_eta_b)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity at the back (-model_basin_comp_layer_eta_b)",data->n_interfaces-1);
	}	
    
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_basin_comp_layer_rho_f",data->rho_f,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer density list at the front (-model_basin_comp_layer_rho_f)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer density at the front (-model_basin_comp_layer_rho_f)",data->n_interfaces-1);
	}

	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_basin_comp_layer_rho_b",data->rho_b,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer density list at the back (-model_basin_comp_layer_rho_b)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer density at the back (-model_basin_comp_layer_rho_b)",data->n_interfaces-1);
	}
	
    
    /* define the mesh size the z-direction for the global problem */
	c->mz = 0;
	for (n=0; n<data->n_interfaces-1; n++) {
		c->mz += data->layer_res_k[n];
	}
    
	data->bc_type = 0; /* 0 use vx compression ; 1 use exx compression */
	data->exx             = -1.0e-3;
	data->vx_commpression = 1.0;
	data->perturbation_type = 0;
	data->perturbation_width = (PetscInt)((2*c->my+1)/4);
    data->layering_type = 0;
	/* parse from command line or input file */
	ierr = PetscOptionsGetInt(NULL,"-model_basin_comp_bc_type",&data->bc_type,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,"-model_basin_comp_perturbation_type",&data->perturbation_type,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,"-model_basin_comp_perturbation_width",&data->perturbation_width,&flg);CHKERRQ(ierr);
    if (data->perturbation_width > (PetscInt)((2*c->my+1)/2) ){
        data->perturbation_width = (PetscInt)((2*c->my+1)/2);
    }
	ierr = PetscOptionsGetReal(NULL,"-model_basin_comp_exx",&data->exx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_basin_comp_vx",&data->vx_commpression,&flg);CHKERRQ(ierr);
    
    ierr = PetscOptionsGetInt(NULL,"-model_basin_comp_layering_type",&data->layering_type,&flg);CHKERRQ(ierr);
    /* Data for layering type 1*/
    data->depth = 1.0;
    data->midleslope = data->Ly/4.0;
    data->slope = data->depth/data->midleslope;
    
    ierr = PetscOptionsGetReal(NULL,"-model_basin_comp_depth",&data->depth,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,"-model_basin_comp_midleslope",&data->midleslope,&flg);CHKERRQ(ierr);
    data->slope = 1.0/data->midleslope;
    ierr = PetscOptionsGetReal(NULL,"-model_basin_comp_slope",&data->slope,&flg);CHKERRQ(ierr);

    
	PetscPrintf(PETSC_COMM_WORLD,"ModelReport: \"Basin Compression\"\n");
	PetscPrintf(PETSC_COMM_WORLD," Domain: [0 , %1.4e] x [0 , %1.4e] x [0 , %1.4e]\n", data->Lx,data->Ly,data->Lz );
	PetscPrintf(PETSC_COMM_WORLD," Mesh:   %.4D x %.4D x %.4D \n", c->mx,c->my,c->mz ); 
    
    n=data->n_interfaces-1;
    	/*
    PetscPrintf(PETSC_COMM_WORLD," ---------------------------- z = %1.4e ----------------------------\n",data->interface_heights[n]);
    PetscPrintf(PETSC_COMM_WORLD,"|\n"); 
    PetscPrintf(PETSC_COMM_WORLD,"|      eta = %1.4e , rho = %1.4e , my = %.4D \n",data->eta[n-1],data->rho[n-1],data->layer_res_k[n-1]);
    PetscPrintf(PETSC_COMM_WORLD,"|\n");
    //PetscPrintf(PETSC_COMM_WORLD,"-\n -\n  -\n   -\n    -\n");

    for (n=data->n_interfaces-1; n>=1; n--) {
		PetscPrintf(PETSC_COMM_WORLD," ---------------------------- z = %1.4e ----------------------------\n",data->interface_heights[n]);
		PetscPrintf(PETSC_COMM_WORLD,"|\n"); 
		PetscPrintf(PETSC_COMM_WORLD,"|      eta = %1.4e , rho = %1.4e , my = %.4D \n",data->eta[n-1],data->rho[n-1],data->layer_res_k[n-1]);
		PetscPrintf(PETSC_COMM_WORLD,"|\n");
	}
	//PetscPrintf(PETSC_COMM_WORLD,"|\n");
	PetscPrintf(PETSC_COMM_WORLD," ---------------------------- z = %1.4e ----------------------------\n",data->interface_heights[0],data->layer_res_k[0]);
	*/
	
	PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "BoundaryCondition_BasinComp"
PetscErrorCode BoundaryCondition_BasinComp(DM dav,BCList bclist,pTatinCtx c,ModelBasinCompCtx *data)
{
	PetscReal         exx, zero = 0.0, vx_E=0.0, vx_W = 0.0;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	exx = data->exx;
    vx_E =-data->vx_commpression;
    vx_W = data->vx_commpression;
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    
    //ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    //ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    
	if (data->bc_type == 0) {
		/* compression east/west in the x-direction (0) [east-west] using constant velocity */

        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&vx_W);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&vx_E);CHKERRQ(ierr);
	} else if (data->bc_type == 1) {
		/* compression east/west in the x-direction (0) [east-west] using constant strain rate */
        ierr = DirichletBC_ApplyStrainRateExx(bclist,dav,exx);CHKERRQ(ierr);
        
		//ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,exx,0);CHKERRQ(ierr);
	} else {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unknonwn boundary condition type");
	}
    

    
	/* free slip south (base) */
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr); 

	
	/* free surface north */
	/* do nothing! */
	;
    
    /* free slip lateral */
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_BasinComp"
PetscErrorCode ModelApplyBoundaryCondition_BasinComp(pTatinCtx c,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscReal         exx;
	BCList            bclist;
	DM                dav;
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	
	exx = data->exx;

	bclist = c->stokes_ctx->u_bclist;
	dav    = c->stokes_ctx->dav;
	ierr = BoundaryCondition_BasinComp(dav,bclist,c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_BasinComp"
PetscErrorCode ModelApplyBoundaryConditionMG_BasinComp(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		ierr = BoundaryCondition_BasinComp(dav[n],bclist[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_BasinComp"
PetscErrorCode ModelApplyMaterialBoundaryCondition_BasinComp(pTatinCtx c,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "BasinCompSetMeshGeometry"
PetscErrorCode BasinCompSetMeshGeometry(DM dav, void *ctx)
{
	PetscErrorCode ierr;
    ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz,M,N,P, kinter_max, kinter_min, interf;
    PetscScalar *dzs, a_b, a_t;
    PetscReal *interface_heights_f, *interface_heights_b, Ly;
    PetscInt *layer_res_k, n_interfaces;

	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    interface_heights_f = data->interface_heights_f;//front heights
    interface_heights_b = data->interface_heights_b;//back heights
    layer_res_k = data->layer_res_k;
    n_interfaces = data->n_interfaces;
    Ly = data->Ly;
    
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
    
    /*Layering type 0: the remeshing is already included*/
    if (data->layering_type == 0 ){
        ierr = PetscMalloc1(ny, &dzs);CHKERRQ(ierr);
        
        kinter_max = 0;
        for(interf = 0; interf < n_interfaces-1; interf++){ 
            kinter_min = kinter_max;
            kinter_max += 2*layer_res_k[interf];
            for(i=si; i<si+nx; i++){
                a_b = (interf == 0)?0.0:(interface_heights_b[interf] - interface_heights_f[interf])/Ly;
                a_t = (interface_heights_b[interf+1] - interface_heights_f[interf+1])/Ly;
                for(j = sj; j<ny+sj; j++){
                    
                    dzs[j-sj] = ((a_t*LA_coord[sk][j][i].y + interface_heights_f[interf+1]) - (a_b*LA_coord[sk][j][i].y + interface_heights_f[interf]))/(PetscReal)(2.0*layer_res_k[interf]);
                }
                for(j=sj; j<sj+ny; j++){
                    PetscScalar h;
                    h = (a_b*LA_coord[sk][j][i].y + interface_heights_f[interf]);
                    for(k=sk;k<sk+nz;k++){
                        if((k <= kinter_max)&&(k >= kinter_min)){
                            LA_coord[k][j][i].z = h + (PetscReal)dzs[j-sj]*(k-kinter_min); 
                            
                        }   
                    }
                }
            }
        
        }
        ierr = PetscFree(dzs);CHKERRQ(ierr);
    }else if (data->layering_type == 1){
            /*Layering type 1: the remeshing between interfaces is not done. Use the function BasinCompRemeshingAccordingToTheInterfaces*/
        kinter_max = 0;
        for(interf = 0; interf < n_interfaces-1; interf++){ 
            kinter_min = kinter_max;
            kinter_max += 2*layer_res_k[interf];
            for(i=si; i<si+nx; i++){
                for(j=sj; j<sj+ny; j++){
                    for(k=sk;k<sk+nz;k++){
                        if((k == kinter_max)&&((interf+1)%2 == 0)){
                            LA_coord[k][j][i].z = interface_heights_f[interf+1];
                        }else if( (k == kinter_min) && (interf%2 != 0)){
                            PetscReal a,b,c,d;
                            a=data->midleslope;
                            b=data->slope*(2.0*a);
                            c=data->Ly/2.0;
                            d=data->depth;
                            LA_coord[k][j][i].z = interface_heights_f[interf] - d/(1.0 + pow(absolute((LA_coord[k][j][i].y - c)/a), 2.0*b));
                                                                          
                        }
                            
                           
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
#define __FUNCT__ "BasinCompRemeshingAccordingToTheInterfaces"
PetscErrorCode BasinCompRemeshingAccordingToTheInterfaces(DM dav, void *ctx)
{   
    
    ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
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
#define __FUNCT__ "BasinCompSetPerturbedInterfaces"
PetscErrorCode BasinCompSetPerturbedInterfaces(DM dav, void *ctx)
{
	PetscErrorCode ierr;
    ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz,M,N,P, interf, kinter, rank, kinter_max, kinter_min;
	PetscScalar pertu, dz_f, dz_b, pertu_f, pertu_b, dz, attenuation, H_f, H_b;
    PetscReal *interface_heights_f, *interface_heights_b, lamb_b, lamb_f;
    PetscInt *layer_res_k, n_interfaces, pwidth;
    PetscReal amp;
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    interface_heights_f = data->interface_heights_f;
    interface_heights_b = data->interface_heights_b;
    layer_res_k = data->layer_res_k;
    n_interfaces = data->n_interfaces;
    amp = data->amp;
    
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	


    
    if (data->perturbation_type <= 1){
    kinter = 0;
    attenuation = -log(0.01)/(PetscScalar)(data->perturbation_width);
    
    MPI_Comm_rank(PetscObjectComm((PetscObject)dav),&rank);
	for(interf = 1; interf < n_interfaces-1; interf++){
		kinter += 2*layer_res_k[interf-1];
		PetscPrintf(PETSC_COMM_WORLD,"jinter = %d (max=%d)\n", kinter,N-1 );
        srand(rank*interf+2);//The seed changes with the interface and the process.
        
		if ( (kinter>=sk) && (kinter<sk+nz) ) {
            /*Take the dominant wavelength of the viscous layer*/
            if(data->eta_b[interf-1] < data->eta_b[interf]){
                H_b = interface_heights_b[interf+1] - interface_heights_b[interf];
                lamb_b = H_b*pow(data->eta_b[interf]/(6.0*data->eta_b[interf-1]), 1.0/3.0);
                
            }else{
                H_b = (interface_heights_b[interf] - interface_heights_b[interf-1]);
                lamb_b = H_b*pow(data->eta_b[interf-1]/(6.0*data->eta_b[interf]), 1.0/3.0);
            }
            if(data->eta_f[interf-1] < data->eta_f[interf]){
                H_f = (interface_heights_f[interf+1] - interface_heights_f[interf]);
                lamb_f = H_f*pow(data->eta_f[interf]/(6.0*data->eta_f[interf-1]), 1.0/3.0);
            }else{
                H_f = (interface_heights_f[interf] - interface_heights_f[interf-1]);
                lamb_f = H_f*pow(data->eta_f[interf-1]/(6.0*data->eta_f[interf]), 1.0/3.0);
            }
            
               
            for(i = si; i<si+nx; i++) {
                
				if((sj+ny == N) && (sj == 0)){
                    j=sj+ny-1;
                    pertu_b =  (data->perturbation_type != 1)?H_b*(2.0 * rand()/(RAND_MAX+1.0) - 1.0):H_b*cos(LA_coord[kinter][j][i].x/lamb_b);
                    for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                        LA_coord[kinter][j-pwidth][i].z += amp * pertu_b*exp(-pwidth*attenuation);
                    }
                    j=0;
                    pertu_f =  (data->perturbation_type != 1)?H_f*(2.0 * rand()/(RAND_MAX+1.0) - 1.0):H_f*cos(LA_coord[kinter][j][i].x/lamb_f);
                    for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                        LA_coord[kinter][j+pwidth][i].z += amp * pertu_f*exp(-pwidth*attenuation);
                    }

				}else if ((sj+ny == N) || (sj == 0)){
                    PetscReal dz = 0.0;
                    PetscInt sgn = 1;
                    sgn = (sj == 0)?1:-1;
                    j = (sj == 0)?0:(sj+ny-1);
                    dz = (sj == 0)?dz_f:dz_b;
                    if (data->perturbation_type == 1){
                        pertu = (sj == 0)?H_f*cos(LA_coord[kinter][j][i].x/lamb_f):H_b*cos(LA_coord[kinter][j][i].x/lamb_b);
                    }else{
                        pertu = (sj == 0)?H_f*(2.0 * rand()/(RAND_MAX+1.0) - 1.0):H_b*(2.0 * rand()/(RAND_MAX+1.0) - 1.0);

                    }
                    for(pwidth = 0; pwidth<data->perturbation_width; pwidth++){
                        LA_coord[kinter][j+sgn*pwidth][i].z += amp * pertu*exp(-pwidth*attenuation);
                    }

				}
			}
			
		}
	}
    }else if (data->perturbation_type == 2){/*perturbe all the interfaces with a red noise (not sure that works properly)*/
        PetscInt kinter_min, kinter_max;
        PetscReal *rnoise;
        
        kinter_min = 0;
        kinter_max = 0;
        MPI_Comm_rank(PetscObjectComm((PetscObject)dav),&rank);
        for(interf = 1; interf < n_interfaces-1; interf++){
            DM botinterface_da, topinterface_da;
            Vec botinterface_coords, topinterface_coords;
            PetscScalar *botinterface_nodes, *topinterface_nodes;
            
            kinter_min = kinter_max;
            kinter_max += 2*layer_res_k[interf-1];
            
            ierr = PetscMalloc1(N*M,&rnoise);CHKERRQ(ierr);
            ierr = PetscMemzero(rnoise,N*M*sizeof(PetscReal));CHKERRQ(ierr);
            
            ierr = rednoise(rnoise, N*M,time(NULL)*(interf+2));CHKERRQ(ierr);


            /*Find the viscous layer*/
            #if 0
            if(data->eta_b[interf-1] < data->eta_b[interf]){
                H_b = interface_heights_b[interf+1] - interface_heights_b[interf];
            }else{
                H_b = (interface_heights_b[interf] - interface_heights_b[interf-1]);           
            }
            if(data->eta_f[interf-1] < data->eta_f[interf]){
                H_f = (interface_heights_f[interf+1] - interface_heights_f[interf]);
            }else{
                H_f = (interface_heights_f[interf] - interface_heights_f[interf-1]);
            }
            #endif
            
            
            if(data->eta_f[interf-1] < data->eta_f[interf]){
                ierr = DMDACreate3dRedundant( dav, si,si+nx, sj,sj+ny, kinter_max, kinter_max + 1, 1, &botinterface_da );CHKERRQ(ierr);
                ierr = DMDACreate3dRedundant( dav, si,si+nx, sj,sj+ny, kinter_max + 2*layer_res_k[interf], kinter_max + 2*layer_res_k[interf]+1, 1, &topinterface_da );CHKERRQ(ierr);
        
            }else{
                ierr = DMDACreate3dRedundant( dav, si,si+nx, sj,sj+ny, kinter_min, kinter_min + 1, 1, &botinterface_da );CHKERRQ(ierr);
                ierr = DMDACreate3dRedundant( dav, si,si+nx, sj,sj+ny, kinter_max, kinter_max + 1, 1, &topinterface_da );CHKERRQ(ierr);   
          
                
            }        

            ierr = DMGetCoordinates( botinterface_da,&botinterface_coords );CHKERRQ(ierr);
            ierr = VecGetArray(botinterface_coords,&botinterface_nodes);CHKERRQ(ierr);
            
            ierr = DMGetCoordinates( topinterface_da,&topinterface_coords );CHKERRQ(ierr);
            ierr = VecGetArray(topinterface_coords,&topinterface_nodes);CHKERRQ(ierr);
            
            if ( (kinter_max>=sk) && (kinter_max<sk+nz) ) {
                for(i = si; i<si+nx; i++) {
                    for(j = sj; j<sj+ny; j++){
                        LA_coord[kinter_max][j][i].z += amp *(topinterface_nodes[3*(0+nx*(j-sj)+(i-si))+2] - botinterface_nodes[3*(0+nx*(j-sj)+(i-si))+2])*rnoise[(i)*ny+(j)];//Perturbe with an amplitude prop to the local thickness.
                        
                    }
                        
                    }
                }
            ierr = VecRestoreArray(botinterface_coords,&botinterface_nodes);CHKERRQ(ierr);
            ierr = VecRestoreArray(topinterface_coords,&topinterface_nodes);CHKERRQ(ierr);
            
            ierr = DMDestroy(&botinterface_da);CHKERRQ(ierr);
            ierr = DMDestroy(&topinterface_da);CHKERRQ(ierr);
            ierr = PetscFree(rnoise);CHKERRQ(ierr);
            
            }
            
    }else if (data->perturbation_type == 3){/*White noise*/
        kinter = 0;
        
        MPI_Comm_rank(PetscObjectComm((PetscObject)dav),&rank);
        
        for(interf = 1; interf < n_interfaces-1; interf++){
            kinter += 2*layer_res_k[interf-1];
            PetscPrintf(PETSC_COMM_WORLD,"jinter = %d (max=%d)\n", kinter,N-1 );
            srand((rank+1)*(interf+1)+1);//The seed changes with the interface and the process.
            
            if ( (kinter>=sk) && (kinter<sk+nz) ) {
                PetscReal H;
                /*Take the dominant wavelength of the viscous layer*/
                if(data->eta_b[interf-1] < data->eta_b[interf]){
                    H_b = interface_heights_b[interf+1] - interface_heights_b[interf];
                }else{
                    H_b = (interface_heights_b[interf] - interface_heights_b[interf-1]);           
                }
                if(data->eta_f[interf-1] < data->eta_f[interf]){
                    H_f = (interface_heights_f[interf+1] - interface_heights_f[interf]);
                }else{
                    H_f = (interface_heights_f[interf] - interface_heights_f[interf-1]);
                }
                H = 0.5*(H_b + H_f);
                for(i = si; i<si+nx; i++) {
                    for(j = sj; j<sj+ny; j++){
                        pertu = 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
                        LA_coord[kinter][j][i].z += amp * H * pertu;
                        
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
#define __FUNCT__ "InitialMaterialGeometryMaterialPoints_BasinComp"
PetscErrorCode InitialMaterialGeometryMaterialPoints_BasinComp(pTatinCtx c,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
    PetscReal              Ly;
	PetscErrorCode ierr;
			
	PetscFunctionBegin;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);CHKERRQ(ierr);
			
			
	/* define properties on material points */
	db = c->materialpoint_db;
    Ly = data->Ly;
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
		PetscReal      eta,rho;
		PetscInt    phase;
		PetscInt    layer, kmaxlayer, kminlayer;
		PetscInt    I, J, K;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);

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
				eta = data->eta_f[layer] + position[1]*(data->eta_b[layer]-data->eta_f[layer])/Ly;
				rho = data->rho_f[layer] + position[1]*(data->rho_b[layer]-data->rho_f[layer])/Ly;

				rho = -rho * GRAVITY;
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
#define __FUNCT__ "InitialMaterialGeometryQuadraturePoints_BasinComp"
PetscErrorCode InitialMaterialGeometryQuadraturePoints_BasinComp(pTatinCtx c,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	int                     p,n_mp_points;
	DataBucket              db;
	DataField               PField_std,PField_stokes;
	PhysCompStokes          user;
	QPntVolCoefStokes       *all_gausspoints,*cell_gausspoints;
	PetscInt                nqp,qp;
    DM                      dav, cda;
    Vec                     gcoords;
    PetscScalar             *LA_gcoords;
    PetscInt                nel,nen_v;
    const PetscInt          *elnidx_v;   
	PetscErrorCode          ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	/* define properties on material points */
	db = c->materialpoint_db;
    dav = c->stokes_ctx->dav;
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));
	
	
	DataBucketGetSizes(db,&n_mp_points,0,0);

    ierr = DMGetCoordinateDM(dav, &cda);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dav,&gcoords );CHKERRQ(ierr);
    ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
    
    ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_v,&elnidx_v);CHKERRQ(ierr);

	/* get the quadrature points */
	user = c->stokes_ctx;
	ierr = VolumeQuadratureGetAllCellData_Stokes(user->volQ,&all_gausspoints);CHKERRQ(ierr);
	nqp = user->volQ->npoints;
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		double      *position;
		PetscReal      eta,rho, Ly;
		PetscInt    phase;
		PetscInt    layer, kmaxlayer, kminlayer, localeid_p;
		PetscInt    I, J, K;
        PetscScalar     elcoords[Q2_NODES_PER_EL_3D*NSD];
        PetscScalar     Ni_p[Q2_NODES_PER_EL_3D], coord_qp[NSD];

		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
        MPntStdGetField_global_coord(material_point,&position);
        
        MPntGetField_global_element_IJKindex(c->stokes_ctx->dav,material_point, &I, &J, &K);
        Ly = data->Ly;
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
				eta = data->eta_f[layer] + position[1]*(data->eta_b[layer]-data->eta_f[layer])/Ly;
				rho = data->rho_f[layer] + position[1]*(data->rho_b[layer]-data->rho_f[layer])/Ly;
			}
			kminlayer += data->layer_res_k[layer];
			layer++;
		}

		
		MPntStdGetField_local_element_index(material_point,&localeid_p);

		ierr = VolumeQuadratureGetCellData_Stokes(user->volQ,all_gausspoints,localeid_p,&cell_gausspoints);CHKERRQ(ierr);

        ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_v[nen_v*localeid_p],LA_gcoords);CHKERRQ(ierr);

            
		for (qp=0; qp<nqp; qp++) {
            PetscReal *xi_qp;
            int i;
            // get local coords for point p
            xi_qp = &user->volQ->q_xi_coor[ NSD * qp ];
            
            P3D_ConstructNi_Q2_3D( xi_qp,Ni_p);
            
            // interpolate nodal coords to quadrature point
            coord_qp[0] = coord_qp[1] = coord_qp[2] = 0.0;
            for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
                coord_qp[0] += Ni_p[i] * elcoords[NSD*i+0];
                coord_qp[1] += Ni_p[i] * elcoords [NSD*i+1];
                coord_qp[2] += Ni_p[i] * elcoords[NSD*i+2];
            }
            
			cell_gausspoints[qp].eta  = data->eta_f[phase-1] + coord_qp[1]*(data->eta_b[phase-1]-data->eta_f[phase-1])/Ly;
			cell_gausspoints[qp].rho  = data->rho_f[phase-1] + coord_qp[1]*(data->rho_b[phase-1]-data->rho_f[phase-1])/Ly;

			cell_gausspoints[qp].Fu[0] = 0.0;
			cell_gausspoints[qp].Fu[1] = -rho * GRAVITY;
			cell_gausspoints[qp].Fu[2] = 0.0;

			cell_gausspoints[qp].Fp = 0.0;
		}		
		
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
    ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_BasinComp"
PetscErrorCode ModelApplyInitialMaterialGeometry_BasinComp(pTatinCtx c,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = InitialMaterialGeometryMaterialPoints_BasinComp(c,ctx);CHKERRQ(ierr);
	ierr = InitialMaterialGeometryQuadraturePoints_BasinComp(c,ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

		
		
		
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_BasinComp"
PetscErrorCode ModelApplyInitialMeshGeometry_BasinComp(pTatinCtx c,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscReal         Lx,Ly,dx,dy,dz,Lz;
	PetscInt          mx,my,mz;
	PetscReal         amp,factor;
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* step 1 - create structured grid */
	Lx = data->Lx;
	Ly = data->Ly;
    Lz = data->Lz;
	
	mx = c->mx; 
	my = c->my; 
	mz = c->mz; 
	
	dx = Lx / ((PetscReal)mx);
	dy = Ly / ((PetscReal)my);
	dz = Lz / ((PetscReal)mz);

	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, 0.0,Lx, 0.0,Ly, data->interface_heights_f[0],Lz);CHKERRQ(ierr);
	factor = 0.1;
	ierr = PetscOptionsGetReal(NULL,"-model_basin_comp_amp_factor",&factor,NULL);CHKERRQ(ierr);
	amp = factor * 1.0; /* this is internal scaled by dy inside BasinCompSetPerturbedInterfaces() */
	/*if ( (amp < 0.0) || (amp >1.0) ) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-model_basin_comp_amp_factor must be 0 < amp < 1");
	}*/
	data->amp = amp;
	/* step 2 - define two interfaces and perturb coords along the interface */
	ierr = BasinCompSetMeshGeometry(c->stokes_ctx->dav, data);CHKERRQ(ierr);
	ierr = BasinCompSetPerturbedInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);
    ierr = BasinCompRemeshingAccordingToTheInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);
    
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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_BasinComp"
PetscErrorCode ModelApplyUpdateMeshGeometry_BasinComp(pTatinCtx c,Vec X,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscReal      step;
	PhysCompStokes stokes;
	DM             stokes_pack,dav,dap;
	Vec            velocity,pressure;
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
#define __FUNCT__ "ModelInitialCondition_BasinComp"
PetscErrorCode ModelInitialCondition_BasinComp(pTatinCtx c,Vec X,void *ctx)
{
    /*
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
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
#define __FUNCT__ "ModelOutput_BasinComp"
PetscErrorCode ModelOutput_BasinComp(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
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
#define __FUNCT__ "ModelDestroy_BasinComp"
PetscErrorCode ModelDestroy_BasinComp(pTatinCtx c,void *ctx)
{
	ModelBasinCompCtx *data = (ModelBasinCompCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_BasinComp"
PetscErrorCode pTatinModelRegister_BasinComp(void)
{
	ModelBasinCompCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelBasinCompCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelBasinCompCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"basin_comp");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_BasinComp);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_BasinComp);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_BasinComp);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_BasinComp);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_BasinComp);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_BasinComp);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_BasinComp);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_BasinComp);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_BasinComp);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

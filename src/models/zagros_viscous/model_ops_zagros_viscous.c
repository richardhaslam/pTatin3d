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

#include "model_zagros_viscous_ctx.h"
#include "model_utils.h"




#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_ZagrosViscous"
PetscErrorCode ModelInitialize_ZagrosViscous(pTatinCtx c,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
	PetscInt n_int,n;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);


	/* assign defaults */
	data->max_layers = 100;
	
	PetscOptionsGetInt(PETSC_NULL,"-model_zagros_viscous_n_interfaces",&data->n_interfaces,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the number of interfaces including the top and bottom boundaries (-model_zagros_viscous_n_interfaces)");
	}
    
	PetscOptionsGetReal(PETSC_NULL,"-model_zagros_viscous_Lx",&data->Lx,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the x direction (-model_zagros_viscous_Lx)");
	}
	
	PetscOptionsGetReal(PETSC_NULL,"-model_zagros_viscous_Ly",&data->Ly,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the y direction (-model_zagros_viscous_Ly)");
	}
    
	n_int = data->n_interfaces;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_zagros_viscous_interface_heights",data->interface_heights,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide back interface heights relative from the base of the model including the top and bottom boundaries.");
	}
	if (n_int != data->n_interfaces) {
        //printf("------>%d %f   %f    %f\n",n_int, data->interface_heights[0], data->interface_heights[1], data->interface_heights[2]);
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d back interface heights relative from the base of the model including the top and bottom boundaries (-model_zagros_viscous_interface_heights)",data->n_interfaces);
        
    }    
	data->Lz = data->interface_heights[data->n_interfaces-1];

    PetscOptionsGetIntArray(PETSC_NULL,"-model_zagros_viscous_layer_res_k",data->layer_res_k,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer resolution list (-model_zagros_viscous_layer_res_k)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_zagros_viscous_layer_res_k)",data->n_interfaces-1);
	}
    
	n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_zagros_viscous_layer_eta",data->eta,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list.(-model_zagros_viscous_layer_eta)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity. (-model_zagros_viscous_layer_eta)",data->n_interfaces-1);
	}
    
    n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_zagros_viscous_layer_rho",data->rho,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer density list. (-model_zagros_viscous_layer_rho)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer density. (-model_zagros_viscous_layer_rho)",data->n_interfaces-1);
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
	/* parse from command line or input file */
	ierr = PetscOptionsGetInt(PETSC_NULL,"-model_zagros_viscous_bc_type",&data->bc_type,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-model_zagros_viscous_perturbation_type",&data->perturbation_type,&flg);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_zagros_viscous_exx",&data->exx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_zagros_viscous_vx",&data->vx_commpression,&flg);CHKERRQ(ierr);
    
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
#define __FUNCT__ "BoundaryCondition_ZagrosViscous"
PetscErrorCode BoundaryCondition_ZagrosViscous(DM dav,BCList bclist,pTatinCtx c,ModelZagrosViscousCtx *data)
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
#define __FUNCT__ "ModelApplyBoundaryCondition_ZagrosViscous"
PetscErrorCode ModelApplyBoundaryCondition_ZagrosViscous(pTatinCtx c,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
	PetscReal         exx;
	BCList            bclist;
	DM                dav;
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	
	exx = data->exx;

	bclist = c->stokes_ctx->u_bclist;
	dav    = c->stokes_ctx->dav;
	ierr = BoundaryCondition_ZagrosViscous(dav,bclist,c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_ZagrosViscous"
PetscErrorCode ModelApplyBoundaryConditionMG_ZagrosViscous(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		ierr = BoundaryCondition_ZagrosViscous(dav[n],bclist[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_ZagrosViscous"
PetscErrorCode ModelApplyMaterialBoundaryCondition_ZagrosViscous(pTatinCtx c,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ZagrosViscousSetMeshGeometry"
PetscErrorCode ZagrosViscousSetMeshGeometry(DM dav, void *ctx)
{
	PetscErrorCode ierr;
    ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
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
#define __FUNCT__ "ZagrosViscousRemeshingAccordingToTheInterfaces"
PetscErrorCode ZagrosViscousRemeshingAccordingToTheInterfaces(DM dav, void *ctx)
{   
    
    ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
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
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(dav,&coord);CHKERRQ(ierr);
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
        
        ierr = DMDAGetCoordinates( botinterface_da,&botinterface_coords );CHKERRQ(ierr);
        ierr = VecGetArray(botinterface_coords,&botinterface_nodes);CHKERRQ(ierr);
        
        ierr = DMDAGetCoordinates( topinterface_da,&topinterface_coords );CHKERRQ(ierr);
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
#define __FUNCT__ "ZagrosViscousSetPerturbedInterfaces"
//Perturbes the interfaces with a whitenoise

PetscErrorCode ZagrosViscousSetPerturbedInterfaces(DM dav, void *ctx)
{
	PetscErrorCode ierr;
    ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz,M,N,P, interf, kinter, rank;
	PetscScalar pertu, dz, H;
    PetscReal *interface_heights;
    PetscInt *layer_res_k, n_interfaces;
    PetscReal amp;
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    interface_heights = data->interface_heights;
    layer_res_k = data->layer_res_k;
    n_interfaces = data->n_interfaces;
    amp = data->amp;
    
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);

    kinter = 0;

    MPI_Comm_rank(((PetscObject)dav)->comm,&rank);
    
	for(interf = 1; interf < n_interfaces-1; interf++){
		kinter += 2*layer_res_k[interf-1];
		//PetscPrintf(PETSC_COMM_WORLD,"kinter = %d (max=%d)\n", kinter,P-1 );
        srand((rank+1)*(interf+1)+1);//The seed changes with the interface and the process.
        
		if ( (kinter>=sk) && (kinter<sk+nz) ) {
            /*Take the dominant wavelength of the viscous layer*/
            if(data->eta[interf-1] < data->eta[interf]){
                H = interface_heights[interf+1] - interface_heights[interf];
                
            }else{
                H = (interface_heights[interf] - interface_heights[interf-1]);
            }
              //PetscPrintf(PETSC_COMM_WORLD,"kinter = %d H = %f\n", kinter,H );
            for(i = si; i<si+nx; i++) {
                for(j = sj; j<sj+ny; j++){
                    
                    pertu = 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
                    LA_coord[kinter][j][i].z += amp * H * pertu;
                    
                }
            
            }
        }
    }
        
   	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);

    
	PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryMaterialPoints_ZagrosViscous"
PetscErrorCode InitialMaterialGeometryMaterialPoints_ZagrosViscous(pTatinCtx c,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
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
				eta = data->eta[layer];
				rho = data->rho[layer];

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
#define __FUNCT__ "InitialMaterialGeometryQuadraturePoints_ZagrosViscous"
PetscErrorCode InitialMaterialGeometryQuadraturePoints_ZagrosViscous(pTatinCtx c,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
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

    ierr = DMDAGetCoordinateDA(dav, &cda);CHKERRQ(ierr);
    ierr = DMDAGetGhostedCoordinates(dav,&gcoords );CHKERRQ(ierr);
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
				eta = data->eta[layer];
				rho = data->rho[layer];
                rho = -rho * GRAVITY;
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
            
			cell_gausspoints[qp].eta  = data->eta[phase-1];
			cell_gausspoints[qp].rho  = data->rho[phase-1];

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
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_ZagrosViscous"
PetscErrorCode ModelApplyInitialMaterialGeometry_ZagrosViscous(pTatinCtx c,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = InitialMaterialGeometryMaterialPoints_ZagrosViscous(c,ctx);CHKERRQ(ierr);
	ierr = InitialMaterialGeometryQuadraturePoints_ZagrosViscous(c,ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

		
		
		
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_ZagrosViscous"
PetscErrorCode ModelApplyInitialMeshGeometry_ZagrosViscous(pTatinCtx c,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
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
	/*
	dx = Lx / ((PetscReal)mx);
	dy = Ly / ((PetscReal)my);
	dz = Lz / ((PetscReal)mz);
    */
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, 0.0,Lx, 0.0,Ly, data->interface_heights[0],Lz);CHKERRQ(ierr);
	factor = 0.1;
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_zagros_viscous_amp_factor",&factor,PETSC_NULL);CHKERRQ(ierr);
	amp = factor * 1.0; /* this is internal scaled by dy inside ZagrosViscousSetPerturbedInterfaces() */
	/*if ( (amp < 0.0) || (amp >1.0) ) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-model_zagros_viscous_amp_factor must be 0 < amp < 1");
	}*/
	data->amp = amp;
	/* step 2 - define two interfaces and perturb coords along the interface */
	ierr = ZagrosViscousSetMeshGeometry(c->stokes_ctx->dav, data);CHKERRQ(ierr);
	ierr = ZagrosViscousSetPerturbedInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);
    ierr = ZagrosViscousRemeshingAccordingToTheInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);
    
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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_ZagrosViscous"
PetscErrorCode ModelApplyUpdateMeshGeometry_ZagrosViscous(pTatinCtx c,Vec X,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
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
#define __FUNCT__ "ModelInitialCondition_ZagrosViscous"
PetscErrorCode ModelInitialCondition_ZagrosViscous(pTatinCtx c,Vec X,void *ctx)
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
#define __FUNCT__ "ModelOutput_ZagrosViscous"
PetscErrorCode ModelOutput_ZagrosViscous(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
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
#define __FUNCT__ "ModelDestroy_ZagrosViscous"
PetscErrorCode ModelDestroy_ZagrosViscous(pTatinCtx c,void *ctx)
{
	ModelZagrosViscousCtx *data = (ModelZagrosViscousCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_ZagrosViscous"
PetscErrorCode pTatinModelRegister_ZagrosViscous(void)
{
	ModelZagrosViscousCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelZagrosViscousCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelZagrosViscousCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"zagros_viscous");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_ZagrosViscous);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_ZagrosViscous);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_ZagrosViscous);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_ZagrosViscous);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_ZagrosViscous);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_ZagrosViscous);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_ZagrosViscous);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_ZagrosViscous);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_ZagrosViscous);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

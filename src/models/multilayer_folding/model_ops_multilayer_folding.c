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
#include "dmda_iterator.h"
#include "material_point_std_utils.h"
#include <time.h>

#include "model_multilayer_folding_ctx.h"
#include "model_utils.h"
#include "ptatin_utils.h"


#undef __FUNCT__ 
#define __FUNCT__ "ModelInitialize_MultilayerFolding"
PetscErrorCode ModelInitialize_MultilayerFolding(pTatinCtx c,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	PetscInt       n_int,n;
	PetscBool      flg;
	RheologyConstants   *rheology;
	DataBucket          materialconstants;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* assign defaults */
	data->max_layers = 100;
	
	data->quasi2d = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_multilayer_folding_quasi2d",&data->quasi2d,&flg);CHKERRQ(ierr);
	
	PetscOptionsGetInt(NULL,"-model_multilayer_folding_n_interfaces",&data->n_interfaces,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the number of interfaces including the top and bottom boundaries (-model_multilayer_folding_n_interfaces)");
	}
	
	PetscOptionsGetReal(NULL,"-model_multilayer_folding_Lx",&data->Lx,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the x direction (-model_multilayer_folding_Lx)");
	}
	
	PetscOptionsGetReal(NULL,"-model_multilayer_folding_Lz",&data->Lz,&flg);
	if (!data->quasi2d) {
		if (!flg) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the z direction (-model_multilayer_folding_Lz)");
		}
	} else {
		data->Lz = 1.0;		
	}
	
	n_int = data->n_interfaces;
	PetscOptionsGetRealArray(NULL,"-model_multilayer_folding_interface_heights",data->interface_heights,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide back interface heights relative from the base of the model including the top and bottom boundaries.");
	}
	if (n_int != data->n_interfaces) {
		//printf("------>%d %f   %f    %f\n",n_int, data->interface_heights[0], data->interface_heights[1], data->interface_heights[2]);
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d back interface heights relative from the base of the model including the top and bottom boundaries (-model_multilayer_folding_interface_heights)",data->n_interfaces);
	}
	
	PetscOptionsGetIntArray(NULL,"-model_multilayer_folding_layer_res_j",data->layer_res_j,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer resolution list (-model_multilayer_folding_layer_res_j)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_multilayer_folding_layer_res_j)",data->n_interfaces-1);
	}
    n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_multilayer_folding_layer_eta",data->eta,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list.(-model_multilayer_folding_layer_eta)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity. (-model_multilayer_folding_layer_eta)",data->n_interfaces-1);
	}
    
    n_int=data->n_interfaces-1;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_multilayer_seed_eta",&data->eta[n_int],&flg); CHKERRQ(ierr);

	data->visco_plastic = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_multilayer_folding_visco_plastic",&data->visco_plastic,&flg);CHKERRQ(ierr);
	if (data->visco_plastic) {
		n_int = data->max_layers;
		PetscOptionsGetRealArray(PETSC_NULL,"-model_multilayer_folding_layer_cohesion",data->cohesion,&n_int,&flg);
		if (!flg) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list.(-model_multilayer_folding_layer_cohesion)");
		}
		if (n_int != data->n_interfaces-1) {
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer cohesion. (-model_multilayer_folding_layer_cohesion)",data->n_interfaces-1);
		}
	}
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(NULL,"-model_multilayer_folding_layer_rho",data->rho,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer density list. (-model_multilayer_folding_layer_rho)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer density. (-model_multilayer_folding_layer_rho)",data->n_interfaces-1);
	}
	
    n_int=data->n_interfaces-1;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-model_multilayer_seed_rho",&data->rho[n_int],&flg); CHKERRQ(ierr);
    
    
	/* define the mesh size the y-direction for the global problem  (gravity in y direction) */
	c->my = 0;
	for (n=0; n<data->n_interfaces-1; n++) {
		c->my += data->layer_res_j[n];
	}
	data->Ly = data->interface_heights[ data->n_interfaces - 1 ];
	
	/* initialize values */
	data->bc_type           = 2; /* 0 use vx/vz compression ; 1 use exx/ezz compression ; 2 use exx/ezz with moving base which conserves volume */
	data->vx_compression    = 1.0;
	data->vz_compression    = 1.0;
	data->exx               = -1.0e-3;
	data->ezz               = -1.0e-3;
	
	data->perturbation_type = 0;
	data->kx                = 0.2;
	data->kz                = 0.2;
	data->A0                = 1.0e-2;
	
	/* parse from command line or input file */
	ierr = PetscOptionsGetInt(NULL,"-model_multilayer_folding_bc_type",&data->bc_type,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,"-model_multilayer_folding_perturbation_type",&data->perturbation_type,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_kx",&data->kx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_kz",&data->kz,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_A0",&data->A0,&flg);CHKERRQ(ierr);
	data->L_char = 1.0;
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_L_char",&data->L_char,&flg); CHKERRQ(ierr); 
	
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_vx",&data->vx_compression,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_vz",&data->vz_compression,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_exx",&data->exx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_ezz",&data->ezz,&flg); CHKERRQ(ierr); 
  data->seed_layer_1 = -1;
  ierr = PetscOptionsGetInt(NULL,"-model_multilayer_seed_layer_1",&data->seed_layer_1,&flg);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"ModelReport: \"Multilayer Folding\"\n");
	PetscPrintf(PETSC_COMM_WORLD," Domain: [0 , %1.4e] x [0 , %1.4e] x [0 , %1.4e]\n", data->Lx,data->Ly,data->Lz );
	PetscPrintf(PETSC_COMM_WORLD," Mesh:   %.4D x %.4D x %.4D \n", c->mx,c->my,c->mz ); 
	
	n = data->n_interfaces-1;
	
	PetscPrintf(PETSC_COMM_WORLD," Layer Depth Profile:\n");
	for (n=data->n_interfaces-1; n>=1; n--) {
		PetscPrintf(PETSC_COMM_WORLD," ---------------------------- y = %1.4e ----------------------------\n",data->interface_heights[n]);
		PetscPrintf(PETSC_COMM_WORLD,"|\n"); 
		if (!data->visco_plastic) {
			PetscPrintf(PETSC_COMM_WORLD,"|      eta = %1.4e , rho = %1.4e , my = %.4D \n",data->eta[n-1],data->rho[n-1],data->layer_res_j[n-1]);
		} else {
			PetscPrintf(PETSC_COMM_WORLD,"|      eta = %1.4e , cohesion = %1.4e , rho = %1.4e , my = %.4D \n",data->eta[n-1],data->cohesion[n-1],data->rho[n-1],data->layer_res_j[n-1]);
		}
		PetscPrintf(PETSC_COMM_WORLD,"|\n");
	}
	PetscPrintf(PETSC_COMM_WORLD," ---------------------------- y = %1.4e ----------------------------\n",data->interface_heights[0],data->layer_res_j[0]);

	/* Rheology prescription */
	ierr = pTatinGetRheology(c,&rheology);CHKERRQ(ierr);
	rheology->rheology_type = RHEOLOGY_VP_STD;
	
	/* Material constant */
	ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
	MaterialConstantsSetDefaults(materialconstants);
	for (n=0; n<data->n_interfaces-1; n++) {
		if (data->visco_plastic) {
			ierr = MaterialConstantsSetValues_MaterialType(materialconstants,   n ,VISCOUS_CONSTANT,PLASTIC_MISES,SOFTENING_NONE,DENSITY_CONSTANT);CHKERRQ(ierr);
			ierr = MaterialConstantsSetValues_PlasticMises(materialconstants,   n, data->cohesion[n],data->cohesion[n]);CHKERRQ(ierr);
		} else {
			ierr = MaterialConstantsSetValues_MaterialType(materialconstants,   n ,VISCOUS_CONSTANT,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);CHKERRQ(ierr);
		}
		ierr = MaterialConstantsSetValues_ViscosityConst(materialconstants, n, data->eta[n]);CHKERRQ(ierr);
		ierr = MaterialConstantsSetValues_DensityConst(materialconstants,   n, data->rho[n]);CHKERRQ(ierr);
	}
	
    if (data->seed_layer_1 == 1){
        ierr = MaterialConstantsSetValues_ViscosityConst(materialconstants, data->n_interfaces-1, data->eta[data->n_interfaces-1]);CHKERRQ(ierr);
        ierr = MaterialConstantsSetValues_DensityConst(materialconstants, data->n_interfaces-1, data->rho[data->n_interfaces-1]);CHKERRQ(ierr);
    }
	
	for (n=0; n<data->n_interfaces-1; n++) {
        PetscPrintf(PETSC_COMM_WORLD,"============== Region index [%D] ============================ \n",n);
        ierr = MaterialConstantsSetFromOptions(materialconstants,"mc_",n,PETSC_FALSE);CHKERRQ(ierr);
        ierr = MaterialConstantsPrintAll(materialconstants,n);CHKERRQ(ierr);
    }
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BoundaryCondition_MultilayerFolding"
PetscErrorCode BoundaryCondition_MultilayerFolding(DM dav,BCList bclist,pTatinCtx user,ModelMultilayerFoldingCtx *data)
{
	PetscReal         zero;
	PetscErrorCode    ierr;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/*          Boundary conditions 
	 
	 base = south face defined by the plan ZOX at Y = 0  
	 top  = north face             // to ZOX at Y = N-1
	 East  Face                    // to ZOY at X = M-1
	 West  Face                    // to ZOY at X = 0
	 lateral face                  // to YOX at Z = 0 
	 lateral face                  // to YOX at Z = K-1
	 */
	switch(data->bc_type) {
			
		case 0:
		{
			PetscReal vx_E=0.0,vx_W = 0.0;
			
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			
			/* compression east/west in the x-direction (0) [east-west] using constant velocity */
			vx_E = -data->vx_compression;
			vx_W =  data->vx_compression;
			//vz_B =  data->vz_compression;
			//vz_F = -data->vz_compression;
			
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&vx_W);CHKERRQ(ierr);
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&vx_E);CHKERRQ(ierr);
			
			/* free slip south (base) */
			zero = 0.0;
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr); 
		}
			break;
			
		case 1:
		{
			PetscReal ezz,exx;
			
			/* compression east/west in the x-direction (0) [east-west] using constant strain rate */
			exx = data->exx;
			ezz = data->ezz;
			
			ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,exx,0);CHKERRQ(ierr);
			ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,ezz,2);CHKERRQ(ierr);
			
			/* free slip south (base) */
			zero = 0.0;
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr); 
		}
			break;
			
		case 2:
		{
			PetscReal ezz,exx;
			PetscReal Max[3],Min[3],Ly,vyB;
			
			/* compression east/west in the x-direction (0) [east-west] using constant strain rate */
			exx = data->exx;
			ezz = data->ezz;
			
			ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,exx,0);CHKERRQ(ierr);
			ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,ezz,2);CHKERRQ(ierr);
			
			/* move base down to accomodate pure thickening */
			ierr = DMDAGetBoundingBox(dav,Min,Max);CHKERRQ(ierr);
			Ly = Max[1]-Min[1];
			
			vyB = 0.5 * (data->exx + data->ezz)*Ly;
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&vyB);CHKERRQ(ierr); 
		}
			break;
			
		default:		
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unknonwn boundary condition type. Valid range is: -model_multilayer_folding_bc_type {0,1,2}");
	}
	
	/* free surface north */
	/* do nothing! */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_MultilayerFolding"
PetscErrorCode ModelApplyBoundaryCondition_MultilayerFolding(pTatinCtx c,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode  ierr;
	BCList          bclist;
	DM              dav;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	bclist = c->stokes_ctx->u_bclist;
	dav    = c->stokes_ctx->dav;
	ierr = BoundaryCondition_MultilayerFolding(dav,bclist,c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_MultilayerFolding"
PetscErrorCode ModelApplyBoundaryConditionMG_MultilayerFolding(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	PetscInt       n;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		ierr = BoundaryCondition_MultilayerFolding(dav[n],bclist[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_MultilayerFolding"
PetscErrorCode ModelApplyMaterialBoundaryCondition_MultilayerFolding(pTatinCtx c,void *ctx)
{
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFoldingSetMeshGeometry"
PetscErrorCode MultilayerFoldingSetMeshGeometry(DM dav, void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	PetscInt       i,j,k,si,sj,sk,nx,ny,nz,M,N,P, jinter_max, jinter_min, interf, n_interfaces,*layer_res_j;
	PetscScalar    dy, *interface_heights;
	DM             cda;
	Vec            coord;
	DMDACoor3d     ***LA_coord;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	n_interfaces = data->n_interfaces;
	layer_res_j = data->layer_res_j;
	interface_heights = data->interface_heights;
	
	jinter_max = 0;
	for (interf=0; interf<n_interfaces-1; interf++) { 
		jinter_min = jinter_max;
		jinter_max += 2*layer_res_j[interf];
		dy = (interface_heights[interf+1] - interface_heights[interf])/(PetscReal)(2.0*layer_res_j[interf]);
		for (i=si; i<si+nx; i++) {
			for (k=sk; k<sk+nz; k++) {
				PetscScalar h;
				
				h = data->interface_heights[interf];
				for (j=sj; j<sj+ny; j++) {
					if ((j <= jinter_max) && (j >= jinter_min)) {
						LA_coord[k][j][i].y = h + (PetscReal)dy*(j-jinter_min); 
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
#define __FUNCT__ "MultilayerFoldingRemeshingAccordingToTheInterfaces"
PetscErrorCode MultilayerFoldingRemeshingAccordingToTheInterfaces(DM dav, void *ctx)
{   
	
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	PetscInt       i,j,k,si,sj,sk,nx,ny,nz, jinter_max, jinter_min, interf;
	PetscScalar    dy;
	PetscInt       *layer_res_j, n_interfaces;
	DM             cda;
	Vec            coord;
	DMDACoor3d     ***LA_coord;
	
    
 	PetscFunctionBegin;
	
	layer_res_j = data->layer_res_j;
	n_interfaces = data->n_interfaces;
	
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	/* Loop again through the layers  to set the perturbation around the layer.*/    
	jinter_max = 0;
	for (interf=0; interf<n_interfaces-1; interf++) { 
		DM botinterface_da, topinterface_da;
		Vec botinterface_coords, topinterface_coords;
		PetscScalar *botinterface_nodes, *topinterface_nodes;
		
		jinter_min = jinter_max;
		jinter_max += 2*layer_res_j[interf];
		
		//		ierr = DMDACreate3dRedundant( dav, si,si+nx, sj,sj+ny, kinter_min, kinter_min + 1, 1, &botinterface_da );CHKERRQ(ierr);
		//		ierr = DMDACreate3dRedundant( dav, si,si+nx, sj,sj+ny, kinter_max, kinter_max+1, 1, &topinterface_da );CHKERRQ(ierr);
		
		ierr = DMDACreate3dRedundant( dav, si,si+nx, jinter_min, jinter_min+1, sk, sk+nz, 1, &botinterface_da );CHKERRQ(ierr);
		ierr = DMDACreate3dRedundant( dav, si,si+nx, jinter_max, jinter_max+1, sk, sk+nz, 1, &topinterface_da );CHKERRQ(ierr);
		
		ierr = DMGetCoordinates( botinterface_da,&botinterface_coords );CHKERRQ(ierr);
		ierr = VecGetArray(botinterface_coords,&botinterface_nodes);CHKERRQ(ierr);
		
		ierr = DMGetCoordinates( topinterface_da,&topinterface_coords );CHKERRQ(ierr);
		ierr = VecGetArray(topinterface_coords,&topinterface_nodes);CHKERRQ(ierr);
		
		for (i=0; i<nx; i++) {
			for (k=0; k<nz; k++) {
				for (j=sj; j<sj+ny; j++) {
					PetscInt idx = i + 0*nx + k*nx*1;
					
					if ((j <= jinter_max)&&(j >= jinter_min)){
						dy = (  topinterface_nodes[3*(idx)+1]  -  botinterface_nodes[3*(idx)+1]  )/((PetscReal)(2.0*layer_res_j[interf]));
						LA_coord[k+sk][j][i+si].y = botinterface_nodes[3*(idx)+1] + (PetscReal)dy*(j-jinter_min); 
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
#define __FUNCT__ "MultilayerFoldingSetPerturbedInterfaces"
PetscErrorCode MultilayerFoldingSetPerturbedInterfaces(DM dav, void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	PetscInt       i,k,si,sj,sk,nx,ny,nz,M,N,P, interf, jinter;
	PetscScalar    pertu, H;
	PetscReal      *interface_heights;
	PetscInt       *layer_res_j, n_interfaces;
	PetscReal      amp;
	DM             cda;
	Vec            coord;
	DMDACoor3d     ***LA_coord;
	PetscMPIInt    rank;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	interface_heights = data->interface_heights;
	layer_res_j = data->layer_res_j;
	n_interfaces = data->n_interfaces;
	amp = data->amp;
	
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	jinter = 0;
	
	ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dav),&rank);CHKERRQ(ierr);
	
	for (interf=1; interf<n_interfaces-1; interf++) {
		jinter += 2*layer_res_j[interf-1];
		//PetscPrintf(PETSC_COMM_WORLD,"kinter = %d (max=%d)\n", kinter,P-1 );
				
		switch (data->perturbation_type) {

			case 0:

				srand((rank+1)*(interf+1)+1); // The seed changes with the interface and the process.		
				if ( (jinter>=sj) && (jinter<sj+ny) ) {
					/* Take the dominant wavelength of the viscous layer */
					if (data->eta[interf-1] < data->eta[interf]) {
						H = interface_heights[interf+1] - interface_heights[interf];
					} else {
						H = (interface_heights[interf] - interface_heights[interf-1]);
					}
					for (i=si; i<si+nx; i++) {
						for (k=sk; k<sk+nz; k++) {
							pertu = 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
							LA_coord[k][jinter][i].y += amp * H * pertu;
						}
					}
				}
				
				break;
			
			case 1:

				if ( (jinter>=sj) && (jinter<sj+ny) ) {
					for (i=si; i<si+nx; i++) {
						for (k=sk; k<sk+nz; k++) {
							PetscReal kx, kz, A0, L_char, x_dimensional, z_dimensional;
							
							kx = data->kx;
							kz = data->kz;
							A0 = data->A0;
							L_char = data->L_char;
							x_dimensional = LA_coord[k][jinter][i].x * L_char;
							z_dimensional = LA_coord[k][jinter][i].z * L_char;
							
							pertu = cos( kx * x_dimensional ) * cos( kz * z_dimensional );
							//						pertu_ = pertu * 1.0/L_char; 
							LA_coord[k][jinter][i].y += A0 * pertu;
						}
					}
				}
				
				break;
            
            case 2: 
            
                if ( (jinter>=sj) && (jinter<sj+ny) ) {
					for (i=si; i<si+nx; i++) {
						for (k=sk; k<sk+nz; k++) {
                
                            if (( (i>0.5*M-2) && (i<0.5*M+2)) || ( (i>0.25*M-2) && (i<0.25*M+2)) || ( (i>0.75*M-2) && (i<0.75*M+2)) ){
                                pertu = 1; 
                            } else{
                                pertu = 0; 
                            }
                         LA_coord[k][jinter][i].y += amp * 0.1 * pertu;
                        
                        }
                    }
                }

                break;

			default:
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown perturbation type selected. Should be {0,1,2}");
				break;
		}
	}
	
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryMaterialPoints_MultilayerFolding"
PetscErrorCode InitialMaterialGeometryMaterialPoints_MultilayerFolding(pTatinCtx c,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	int            p,n_mp_points;
	DataBucket     db;
	DataField      PField_std,PField_stokes;
	
    
	PetscFunctionBegin;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);CHKERRQ(ierr);
	
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
		PetscReal   eta,rho;
		PetscInt    region_idx,phase;
		PetscInt    layer, jmaxlayer, jminlayer;
		PetscInt    nI, nJ, nK, nJmid;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		MPntGetField_global_element_nInJnKindex(c->stokes_ctx->dav,material_point, &nI, &nJ, &nK);
		phase = -1;
		eta =  0.0;
		rho = 0.0;
		jmaxlayer = jminlayer = 0;
		layer = 0;
		region_idx = -1;
		// gets the global element index (i,j,k)
		//....
        
        nJmid = data->layer_res_j[0] + data->layer_res_j[1]/2.0;     
		
		//Set the properties
		while( (phase == -1) && (layer < data->n_interfaces-1) ){
			jmaxlayer += data->layer_res_j[layer];
			
			if( (nJ<jmaxlayer) && (nJ>=jminlayer) ){
				phase = layer + 1;
                
                if( ((data->seed_layer_1 ==1) && (nJ<nJmid+2) && (nJ>nJmid-2) && (nI<c->mx*0.25+2) && (nI>c->mx*0.25-2))
                    ||
                    ((data->seed_layer_1 ==1) && (nJ<nJmid+2) && (nJ>nJmid-2) && (nI<c->mx*0.75+2) && (nI>c->mx*0.75-2)) ){
                eta = data->eta[data->n_interfaces-1]; 
                rho = data->rho[data->n_interfaces-1];
                region_idx = data->n_interfaces-1;   
                    
                }
                else{
				eta = data->eta[layer];
                rho = data->rho[layer];  
                region_idx = layer;
				}
            
				rho = -rho * GRAVITY;
		
            }
			jminlayer += data->layer_res_j[layer];
			layer++;
		}
		if (region_idx == -1) {
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Failed to located material point within any layer");
		}
		/* user the setters provided for you */
		
		MPntStdSetField_phase_index(material_point,region_idx);
		MPntPStokesSetField_eta_effective(mpprop_stokes,eta);
		MPntPStokesSetField_density(mpprop_stokes,rho);
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFolding_InitialMaterialGeometry_DamageMP"
PetscErrorCode MultilayerFolding_InitialMaterialGeometry_DamageMP(pTatinCtx c,ModelMultilayerFoldingCtx *data)
{
	PetscErrorCode ierr;
	MPAccess       mpX;
	int            p,n_mp_points;
	DataBucket     db;
    DM             stokes_pack,dav,dap;
	PhysCompStokes stokes;
	PetscReal dmin,dmax,drange;
    PetscInt  damagelist[100],maxlayers=100;
    PetscBool layer2damge[100];
    
    
	PetscFunctionBegin;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);CHKERRQ(ierr);
    
    dmin = 0.0;
    dmax = 0.0;
    PetscOptionsGetReal(PETSC_NULL,"-model_multilayer_folding_damage_init_min",&dmin,0);
    PetscOptionsGetReal(PETSC_NULL,"-model_multilayer_folding_damage_init_max",&dmax,0);
    drange = dmax - dmin;

    for (p=0; p<100; p++) {
        damagelist[p] = -1;
        layer2damge[p] = PETSC_FALSE;
    }
    PetscOptionsGetIntArray(PETSC_NULL,"-model_multilayer_folding_damage_layer",damagelist,&maxlayers,0);
    for (p=0; p<maxlayers; p++) {
        /* check indices are valid */
        if (damagelist[p] < 0) {
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"InitialDamage: Cannot set a layer index to be less than 0");
        }
        if (damagelist[p] >= data->n_interfaces-1) {
            SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"InitialDamage: Cannot set a layer index greater than max layers (%D)",data->n_interfaces-2);
        }
        
        layer2damge[ damagelist[p] ] = PETSC_TRUE;
    }
	
    for (p=0; p<data->n_interfaces-1; p++) {
        if (!layer2damge[p]) {
            PetscPrintf(PETSC_COMM_WORLD,"  layer[%D] : undamaged \n",p);
        } else {
            PetscPrintf(PETSC_COMM_WORLD,"  layer[%D] : initial damage [+%1.4e,+%1.4e] \n",p,dmin,dmax);
        }
    }
    
    /* get velocity dm */
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    
	/* get material material points */
	ierr = pTatinGetMaterialPoints(c,&db,PETSC_NULL);CHKERRQ(ierr);
	ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
	DataBucketGetSizes(db,&n_mp_points,0,0);

    ptatin_RandomNumberSetSeedRank(PETSC_COMM_WORLD);
    
	for (p=0; p<n_mp_points; p++) {
		double    *position_p;
		double    eta,rho;
        float     pl_strain;
		int       region_idx,phase;
		PetscInt  layer, jmaxlayer, jminlayer;
		PetscInt  nI,nJ,nK;
        int       wil_p;
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		ierr = MaterialPointGet_global_coord(mpX,p,&position_p);CHKERRQ(ierr);
		
        ierr = MaterialPointGet_local_element_index(mpX,p,&wil_p);CHKERRQ(ierr);

		//MPntGetField_global_element_nInJnKindex(c->stokes_ctx->dav,material_point, &nI, &nJ, &nK);
        ierr = DMDAConvertLocalElementIndex2GlobalnInJnK(dav,wil_p,&nI,&nJ,&nK);CHKERRQ(ierr);
		phase      = -1;
		eta        =  0.0;
		rho        = 0.0;
		jmaxlayer  = jminlayer = 0;
		layer      = 0;
		region_idx = -1;
		
		//Set the properties
		while ((phase == -1) && (layer < data->n_interfaces-1)) {
			jmaxlayer += data->layer_res_j[layer];
			
			if ((nJ < jmaxlayer) && (nJ >= jminlayer)) {
				phase      = layer + 1;
				eta        = (double)data->eta[layer];
				rho        = (double)data->rho[layer];
				region_idx = layer;
				
				rho = -rho * GRAVITY;
			}
			jminlayer += data->layer_res_j[layer];
			layer++;
		}
		if (region_idx == -1) {
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Failed to located material point within any layer");
		}
        
		/* user the setters provided for you */
        ierr = MaterialPointSet_phase_index(mpX,p, region_idx);CHKERRQ(ierr);
        ierr = MaterialPointSet_viscosity(mpX,  p, eta);CHKERRQ(ierr);
        ierr = MaterialPointSet_density(mpX,    p, rho);CHKERRQ(ierr);

        /* set initial plastic strain */
        pl_strain = 0.0;
        if (layer2damge[region_idx]) {
            double tmp;
            
            tmp = ptatin_RandomNumberGetDouble(dmin,dmax);
            pl_strain = (float)tmp;
        }
        ierr = MaterialPointSet_plastic_strain(mpX,p,pl_strain);CHKERRQ(ierr);
    }
    ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* DAML This function is totally fucking inefficient */
#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryQuadraturePoints_MultilayerFolding"
PetscErrorCode InitialMaterialGeometryQuadraturePoints_MultilayerFolding(pTatinCtx c,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode          ierr;
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
		PetscInt    phase;
		PetscInt    layer, jmaxlayer, jminlayer, localeid_p;
		PetscInt    nI, nJ, nK, nJmid;
		PetscScalar elcoords[Q2_NODES_PER_EL_3D*NSD];
		PetscScalar Ni_p[Q2_NODES_PER_EL_3D], coord_qp[NSD];
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		MPntStdGetField_global_coord(material_point,&position);
		
		MPntGetField_global_element_nInJnKindex(c->stokes_ctx->dav,material_point, &nI, &nJ, &nK);
		//Set the properties
		phase     = -1;
		jmaxlayer = jminlayer = 0;
		layer     = 0;
        
        nJmid = data->layer_res_j[0] + data->layer_res_j[1]/2.0;
        
		while ( (phase == -1) && (layer < data->n_interfaces-1) ) {
			jmaxlayer += data->layer_res_j[layer];
			
			if ( (nJ<jmaxlayer) && (nJ>=jminlayer) ) {
				phase = layer + 1;
                
                if( ((data->seed_layer_1 ==1) && (nJ<nJmid+2) && (nJ>nJmid-2) && (nI<c->mx*0.25+2) && (nI>c->mx*0.25-2))
                   ||
                   ((data->seed_layer_1 ==1) && (nJ<nJmid+2) && (nJ>nJmid-2) && (nI<c->mx*0.75+2) && (nI>c->mx*0.75-2)) ){
                    phase = data->n_interfaces;
                }
			}
			jminlayer += data->layer_res_j[layer];
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
				coord_qp[1] += Ni_p[i] * elcoords[NSD*i+1];
				coord_qp[2] += Ni_p[i] * elcoords[NSD*i+2];
			}
			
			cell_gausspoints[qp].eta  = data->eta[phase-1];
			cell_gausspoints[qp].rho  = data->rho[phase-1];
			
			cell_gausspoints[qp].Fu[0] = 0.0;
			cell_gausspoints[qp].Fu[1] = -data->rho[phase-1] * GRAVITY;
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
#define __FUNCT__ "_InitialMaterialGeometryQuadraturePoints_MultilayerFolding"
PetscErrorCode _InitialMaterialGeometryQuadraturePoints_MultilayerFolding(pTatinCtx c,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PhysCompStokes            stokes;
	QPntVolCoefStokes         *all_gausspoints,*cell_gausspoints;
	DM                        stokes_pack,dav,dap;
	PetscInt                  e,nel,nen_v,nqp,qp;
	const PetscInt            *elnidx_v;   
	PetscInt                  phase;
	PetscInt                  layer,jmaxlayer,jminlayer;
	PetscInt                  nI,nJ,nK,nJmid;
	PetscErrorCode            ierr;
	

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
		
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_v,&elnidx_v);CHKERRQ(ierr);
	
	/* get the quadrature points */
	ierr = VolumeQuadratureGetAllCellData_Stokes(stokes->volQ,&all_gausspoints);CHKERRQ(ierr);
	nqp = stokes->volQ->npoints;
	
     nJmid = data->layer_res_j[0] + data->layer_res_j[1]/2.0;  
    
	for (e=0; e<nel; e++) {
		ierr = DMDAConvertLocalElementIndex2GlobalnInJnK(dav,e, &nI,&nJ,&nK);CHKERRQ(ierr);

		// Determine phase from layer
		phase = -1;
		jmaxlayer = jminlayer = 0;
		layer = 0;
		while ( (phase == -1) && (layer < data->n_interfaces-1) ) {
			jmaxlayer += data->layer_res_j[layer];
			
			if ( (nJ < jmaxlayer) && (nJ >= jminlayer) ) {
                if( ((data->seed_layer_1 ==1) && (nJ<nJmid+2) && (nJ>nJmid-2) && (nI<c->mx*0.25+2) && (nI>c->mx*0.25-2))
                    ||
                    ((data->seed_layer_1 ==1) && (nJ<nJmid+2) && (nJ>nJmid-2) && (nI<c->mx*0.75+2) && (nI>c->mx*0.75-2)) ){
                    phase = data->n_interfaces;    
                }
                else{
				phase = layer + 1;
                }
			}
			jminlayer += data->layer_res_j[layer];
			layer++;
		}
		if (phase == -1) {
			SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot locate the phase element %D is associated with",e);
		}
		
		ierr = VolumeQuadratureGetCellData_Stokes(stokes->volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		for (qp=0; qp<nqp; qp++) {
			cell_gausspoints[qp].eta  = data->eta[phase-1];
			cell_gausspoints[qp].rho  = data->rho[phase-1];
			
			cell_gausspoints[qp].Fu[0] = 0.0;
			cell_gausspoints[qp].Fu[1] = -data->rho[phase-1] * GRAVITY;
			cell_gausspoints[qp].Fu[2] = 0.0;
			
			cell_gausspoints[qp].Fp = 0.0;
        }
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFolding_SetMaterialPointPropertiesFromLayer"
PetscErrorCode MultilayerFolding_SetMaterialPointPropertiesFromLayer(pTatinCtx c,ModelMultilayerFoldingCtx *data)
{
	PetscErrorCode   ierr;
	int              p,n_mpoints;
	DataBucket       materialpoint_db;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	MPAccess         mpX;
	
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	/* define properties on material points */
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	for (p=0; p<n_mpoints; p++) {
		PetscInt    phase;
		PetscInt    layer,ei,localeid_p;
		PetscInt    nI,nJ,nK,nJmid;
		
		ierr = MaterialPointGet_local_element_index(mpX,p,&localeid_p);CHKERRQ(ierr);
		ierr = DMDAConvertLocalElementIndex2GlobalnInJnK(dav,localeid_p,&nI,&nJ,&nK);CHKERRQ(ierr);

		/* Set the properties based on the J index of the element containing the marker */
		phase = -1;
		ei = 0;
        nJmid = data->layer_res_j[0] + data->layer_res_j[1]/2.0; 
        
		for (layer=0; layer<data->n_interfaces-1; layer++) {
		
			if ((nJ >= ei) && (nJ <ei + data->layer_res_j[layer])) {
                if( ((data->seed_layer_1 ==1) && (nJ<nJmid+2) && (nJ>nJmid-2) && (nI<c->mx*0.25+2) && (nI>c->mx*0.25-2))
                    ||
                    ((data->seed_layer_1 ==1) && (nJ<nJmid+2) && (nJ>nJmid-2) && (nI<c->mx*0.75+2) && (nI>c->mx*0.75-2)) ){
                    phase = data->n_interfaces-1;    
                }
                else{
				phase = layer;
                }    
                break;
			}
			ei = ei + data->layer_res_j[layer];
		}
		if (phase == -1) {
			SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Couldn't identify marker %D in any layer",p);
		}
		
		ierr = MaterialPointSet_viscosity(  mpX,p, data->eta[phase]);CHKERRQ(ierr);
		ierr = MaterialPointSet_density(    mpX,p,-data->rho[phase] * GRAVITY);CHKERRQ(ierr);
		ierr = MaterialPointSet_phase_index(mpX,p,phase);CHKERRQ(ierr);
        
    }
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_MultilayerFolding"
PetscErrorCode ModelApplyInitialMaterialGeometry_MultilayerFolding(pTatinCtx c,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode         ierr;
    PetscInt               mp_geom;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
    mp_geom = 0;
    PetscOptionsGetInt(PETSC_NULL,"-model_multilayer_folding_mp_geom",&mp_geom,0);
    
    if (mp_geom == 0) {
        ierr = InitialMaterialGeometryMaterialPoints_MultilayerFolding(c,ctx);CHKERRQ(ierr);
        ierr = InitialMaterialGeometryQuadraturePoints_MultilayerFolding(c,ctx);CHKERRQ(ierr);
    } else if(mp_geom == 1) {
        ierr = MultilayerFolding_InitialMaterialGeometry_DamageMP(c,data);CHKERRQ(ierr);
    }
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_MultilayerFolding"
PetscErrorCode ModelApplyInitialMeshGeometry_MultilayerFolding(pTatinCtx c,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode    ierr;
	PetscReal         Lx,Ly,Lz;
	PetscReal         amp,factor;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* step 1 - create structured grid */
	Lx = data->Lx;
	Ly = data->Ly;
	Lz = data->Lz;
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, 0.0,Lx, data->interface_heights[0],Ly, 0.0,Lz);CHKERRQ(ierr);
	factor = 0.1;
	ierr = PetscOptionsGetReal(NULL,"-model_multilayer_folding_amp_factor",&factor,NULL);CHKERRQ(ierr);
	amp = factor * 1.0; /* this is internal scaled by dy inside MultilayerFoldingSetPerturbedInterfaces() */
	/*if ( (amp < 0.0) || (amp >1.0) ) {
	 SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-model_multilayer_folding_amp_factor must be 0 < amp < 1");
	 }*/
	data->amp = amp;
	/* step 2 - define two interfaces and perturb coords along the interface */
	ierr = MultilayerFoldingSetMeshGeometry(c->stokes_ctx->dav, data);CHKERRQ(ierr);
	ierr = MultilayerFoldingSetPerturbedInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);
	ierr = MultilayerFoldingRemeshingAccordingToTheInterfaces(c->stokes_ctx->dav, data);CHKERRQ(ierr);
	
	ierr = DMDABilinearizeQ2Elements(c->stokes_ctx->dav);CHKERRQ(ierr);

	if (data->quasi2d) {
		ierr = pTatin3d_DefineVelocityMeshGeometryQuasi2D(c);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

/* 
DAM: I'm really not convinced this function is
 i) a good idea, and 
 ii) even needed

 For instance, this function assumes that there are 2x2x2 material points per element.
 It's not need as 
   ModelApplyInitialMaterialGeometry_MultilayerFolding()
 calls both
   InitialMaterialGeometryMaterialPoints_MultilayerFolding();
   InitialMaterialGeometryQuadraturePoints_MultilayerFolding()
 
 The only useful thing I see that this function does is switch the projection type.
*/
#undef __FUNCT__
#define __FUNCT__ "MultilayerFolding_Mesh2MarkerRemesh"
PetscErrorCode MultilayerFolding_Mesh2MarkerRemesh(pTatinCtx c,ModelMultilayerFoldingCtx *data)
{
	PetscErrorCode ierr;
	PetscInt       Nxp[] = { 2, 2, 2 };
	PetscReal      perturb = 0.1;
	DataBucket     materialpoint_db;
	PhysCompStokes stokes;
	DM             stokes_pack,dav,dap;
	
    
	PetscPrintf(PETSC_COMM_WORLD,"[[MultilayerFoldingMesh2MarkerRemesh]]\n");
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);

	/* Force new coordinate layout. This won't allocate any memory, it will just reset points into the current mesh configuration */
	ierr = SwarmMPntStd_CoordAssignment_LatticeLayout3d(dav,Nxp,perturb,materialpoint_db);CHKERRQ(ierr);
	
	/* Re-assign material properties based on element index */
	ierr = MultilayerFolding_SetMaterialPointPropertiesFromLayer(c,data);CHKERRQ(ierr);
	
	/* switch coefficient projection type to use Q1 rather than null projection */
	c->coefficient_projection_type = 1;
	
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
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_MultilayerFolding"
PetscErrorCode ModelApplyUpdateMeshGeometry_MultilayerFolding(pTatinCtx c,Vec X,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	PetscReal      step;
	PhysCompStokes stokes;
	DM             stokes_pack,dav,dap;
	Vec            velocity,pressure;
	PetscInt           metric_L = 5; 
	MeshQualityMeasure metric_list[] = { MESH_QUALITY_ASPECT_RATIO, MESH_QUALITY_DISTORTION, MESH_QUALITY_DIAGONAL_RATIO, MESH_QUALITY_VERTEX_ANGLE, MESH_QUALITY_FACE_AREA_RATIO };
	PetscReal          value[100];
	PetscBool          basal_remesh = PETSC_FALSE;
	PetscBool          marker_remesh = PETSC_FALSE;
	static PetscBool   tracking_layer_phase = PETSC_TRUE;
	
    
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
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_ASPECT_RATIO\"    %1.4e \n", value[0]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_DISTORTION\"      %1.4e \n", value[1]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_DIAGONAL_RATIO\"  %1.4e \n", value[2]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_VERTEX_ANGLE\"    %1.4e \n", value[3]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_FACE_AREA_RATIO\" %1.4e \n", value[4]);
	

	ierr = PetscOptionsGetBool(NULL,"-model_multilayer_basal_remesh",&basal_remesh,NULL);CHKERRQ(ierr);
	if (basal_remesh && (value[0] > 20.0)) {
		/* Remesh only the first layer if the aspect ratio is > 20 */
		PetscInt   res_j;
		DMDACoor3d span_xz[4];
		PetscReal  gmin[3],gmax[3];
		
		PetscPrintf(PETSC_COMM_WORLD,"[[%s]] Activating basal layer remeshing\n", __FUNCT__);
		//ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,"before");CHKERRQ(ierr);

		/* clean up base */
		ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
		span_xz[0].x = gmin[0];    span_xz[0].y = gmin[1];    span_xz[0].z = gmin[2];
		span_xz[1].x = gmin[0];    span_xz[1].y = gmin[1];    span_xz[1].z = gmax[2];
		span_xz[2].x = gmax[0];    span_xz[2].y = gmin[1];    span_xz[2].z = gmax[2];
		span_xz[3].x = gmax[0];    span_xz[3].y = gmin[1];    span_xz[3].z = gmin[2];
		ierr = DMDARemeshSetUniformCoordinatesInPlane_IK(dav,0,span_xz);CHKERRQ(ierr);
		//ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,"baseremesh");CHKERRQ(ierr);
		
		/* interpolate between base and top of first layer */
		res_j = 2 * data->layer_res_j[0] + 1; /* use nodal indices */
		ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,res_j);CHKERRQ(ierr);
		//ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,"after");CHKERRQ(ierr);
	}
	
	ierr = PetscOptionsGetBool(NULL,"-model_multilayer_marker_remesh",&marker_remesh,NULL);CHKERRQ(ierr);
	if (marker_remesh) {
		PetscInt   JMAX;
		DMDACoor3d span_xz[4];
		PetscReal  gmin[3],gmax[3];
		
		ierr = DMDAGetInfo(dav,0,0,&JMAX,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

		ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
		span_xz[0].x = gmin[0];    span_xz[0].y = gmin[1];    span_xz[0].z = gmin[2];
		span_xz[1].x = gmin[0];    span_xz[1].y = gmin[1];    span_xz[1].z = gmax[2];
		span_xz[2].x = gmax[0];    span_xz[2].y = gmin[1];    span_xz[2].z = gmax[2];
		span_xz[3].x = gmax[0];    span_xz[3].y = gmin[1];    span_xz[3].z = gmin[2];
		
		if ((tracking_layer_phase == PETSC_TRUE) && (value[0] > 20.0)) {
		//if ((tracking_layer_phase == PETSC_TRUE) && (c->step >= 30)) {
			PetscPrintf(PETSC_COMM_WORLD,"[[%s]] Activating marker remeshing\n", __FUNCT__);
		
			/* DAM: I think this function farts higher than its arse and should be removed */
			/* reset material point coordinates and set eta/rho */
			//ierr = MultilayerFolding_Mesh2MarkerRemesh(c,data);CHKERRQ(ierr);

			/* 0 switch projection type */
			c->coefficient_projection_type = 1;

			/* 1 - re-initialize the basement layer element spacing */
			ierr = DMDARemeshSetUniformCoordinatesInPlane_IK(dav,0,span_xz);CHKERRQ(ierr);

			/* 2 - clean up the interior */
			ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,JMAX);CHKERRQ(ierr);
			tracking_layer_phase = PETSC_FALSE;
			
		} else if (tracking_layer_phase == PETSC_FALSE) {
			PetscPrintf(PETSC_COMM_WORLD,"[[%s]] Subsequent marker remeshing\n", __FUNCT__);
			/* 1 - re-initialize the basement layer element spacing */
			ierr = DMDARemeshSetUniformCoordinatesInPlane_IK(dav,0,span_xz);CHKERRQ(ierr);
			
			/* 2 - clean up the interior */
			ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,JMAX);CHKERRQ(ierr);
		}
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFoldingUpdate_RemeshBasalLayer"
PetscErrorCode MultilayerFoldingUpdate_RemeshBasalLayer(PetscReal AR,DM dav,ModelMultilayerFoldingCtx *data)
{
	PetscInt       res_j;
	DMDACoor3d     span_xz[4];
	PetscReal      gmin[3],gmax[3];
	PetscReal      AR_max = 20.0;
	PetscErrorCode ierr;
	
    
	/* Remesh only the first layer if the aspect ratio is > AR */
	if (AR > AR_max) {
		PetscPrintf(PETSC_COMM_WORLD,"[[%s]] Activating basal layer remeshing\n", __FUNCT__);
		
		/* clean up base */
		ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
		span_xz[0].x = gmin[0];    span_xz[0].y = gmin[1];    span_xz[0].z = gmin[2];
		span_xz[1].x = gmin[0];    span_xz[1].y = gmin[1];    span_xz[1].z = gmax[2];
		span_xz[2].x = gmax[0];    span_xz[2].y = gmin[1];    span_xz[2].z = gmax[2];
		span_xz[3].x = gmax[0];    span_xz[3].y = gmin[1];    span_xz[3].z = gmin[2];
		ierr = DMDARemeshSetUniformCoordinatesInPlane_IK(dav,0,span_xz);CHKERRQ(ierr);
		
		/* interpolate between base and top of first layer */
		res_j = 2 * data->layer_res_j[0] + 1; /* use nodal indices */
		ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,res_j);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);	
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFoldingUpdate_RemeshMarkerProjection_v1"
PetscErrorCode MultilayerFoldingUpdate_RemeshMarkerProjection_v1(PetscReal AR,DM dav,pTatinCtx c)
{
	PetscInt         JMAX;
	DMDACoor3d       span_xz[4];
	PetscReal        gmin[3],gmax[3];
	PetscReal        AR_max = 20.0;
	static PetscBool tracking_layer_phase = PETSC_TRUE;
	PetscErrorCode   ierr;
	
    
	ierr = DMDAGetInfo(dav,0,0,&JMAX,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	
	ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
	span_xz[0].x = gmin[0];    span_xz[0].y = gmin[1];    span_xz[0].z = gmin[2];
	span_xz[1].x = gmin[0];    span_xz[1].y = gmin[1];    span_xz[1].z = gmax[2];
	span_xz[2].x = gmax[0];    span_xz[2].y = gmin[1];    span_xz[2].z = gmax[2];
	span_xz[3].x = gmax[0];    span_xz[3].y = gmin[1];    span_xz[3].z = gmin[2];
	
	if (AR > AR_max) {
		PetscPrintf(PETSC_COMM_WORLD,"[[%s]] Activating marker remeshing [projection]\n", __FUNCT__);

		if (tracking_layer_phase ) {
			PetscPrintf(PETSC_COMM_WORLD,"[[%s]]   Switching to material point projection\n", __FUNCT__);
			
			/* 0 switch projection type */
			c->coefficient_projection_type = 1;
			tracking_layer_phase = PETSC_FALSE;
		} 

		PetscPrintf(PETSC_COMM_WORLD,"[[%s]]   Performing remeshing\n", __FUNCT__);
		/* 1 - re-initialize the basement layer element spacing */
		ierr = DMDARemeshSetUniformCoordinatesInPlane_IK(dav,0,span_xz);CHKERRQ(ierr);
		/* 2 - clean up the interior */
		ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,JMAX);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFoldingUpdate_RemeshMarkerProjection_v2"
PetscErrorCode MultilayerFoldingUpdate_RemeshMarkerProjection_v2(PetscReal AR,DM dav,Vec velocity,PetscReal dt,ModelMultilayerFoldingCtx *data,pTatinCtx c)
{
	Vec              mesh_velocity;
	PetscReal        AR_max = 20.0;
	static PetscBool tracking_layer_phase = PETSC_TRUE;
	DMDAVecTraverse3d_InterpCtx IntpCtx;
	PetscReal        Lx,Lz;
	PetscReal        gmin[3],gmax[3];
	PetscErrorCode   ierr;

	
	if (AR > AR_max) {
		PetscPrintf(PETSC_COMM_WORLD,"[[%s]] Activating marker remeshing [projection]\n", __FUNCT__);
		
		/*
		 Reset position of mesh 
		 This is required as (i) AR is estimated on the deformed mesh
		 (ii) UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX() advects the free surface 
		 */
		ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,-dt);CHKERRQ(ierr);
		
		if (tracking_layer_phase ) {
			PetscPrintf(PETSC_COMM_WORLD,"[[%s]]   Switching to material point projection\n", __FUNCT__);
			
			/* 0 switch projection type */
			c->coefficient_projection_type = 1;
			tracking_layer_phase = PETSC_FALSE;
		} 
		
		PetscPrintf(PETSC_COMM_WORLD,"[[%s]]   Performing remeshing\n", __FUNCT__);
		
		/* [A] create mesh advection velocity field in x-z */
		ierr = DMCreateGlobalVector(dav,&mesh_velocity);CHKERRQ(ierr);
		
		ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
		Lx = gmax[0] - gmin[0];
		Lz = gmax[2] - gmin[2];
		
		if (data->bc_type == 0) {
			PetscReal vx_E,vx_W;
			/* compression east/west in the x-direction (0) [east-west] using constant velocity */

			vx_E = -data->vx_compression;
			vx_W =  data->vx_compression;

			ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vx_E - vx_W)/(Lx),vx_W,0.0);CHKERRQ(ierr);
			ierr = DMDAVecTraverse3d(dav,mesh_velocity, 0, DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);

			ierr = VecStrideSet(mesh_velocity,1,0.0);CHKERRQ(ierr); // zero y component //
			ierr = VecStrideSet(mesh_velocity,2,0.0);CHKERRQ(ierr); // zero z component //
		} else if (data->bc_type == 1) {
			PetscReal vx_E,vx_W;
			PetscReal vz_F,vz_B;
			PetscReal cx[3];
			
			/* center of domain */
			cx[0] = 0.5 * (gmax[0] + gmin[0]);
			cx[1] = 0.0;
			cx[2] = 0.5 * (gmax[2] + gmin[2]);
			
			/* x component */
			vx_E = data->exx * ( gmax[0] - cx[0] );
			vx_W = data->exx * ( gmin[0] - cx[0] );
			PetscPrintf(PETSC_COMM_WORLD,"  x: pos(min/max) = %+1.4e / %+1.4e  vx = %+1.4e / %+1.4e \n",gmin[0],gmax[0],vx_W,vx_E);
			
			ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vx_E - vx_W)/(Lx),vx_W,0.0);CHKERRQ(ierr);
			ierr = DMDAVecTraverse3d(dav,mesh_velocity, 0, DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);

			/* y component */
			ierr = VecStrideSet(mesh_velocity,1,0.0);CHKERRQ(ierr); // zero y component //

			/* z component */
			vz_F = data->ezz * ( gmax[2] - cx[2] );
			vz_B = data->ezz * ( gmin[2] - cx[2] );
			PetscPrintf(PETSC_COMM_WORLD,"  z: pos(min/max) = %+1.4e / %+1.4e  vz = %+1.4e / %+1.4e \n",gmin[2],gmax[2],vz_B,vz_F);

			ierr = DMDAVecTraverse3d_InterpCtxSetUp_Z(&IntpCtx,(vz_F - vz_B)/(Lz),vz_B,0.0);CHKERRQ(ierr);
			ierr = DMDAVecTraverse3d(dav,mesh_velocity, 2, DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
		} else if (data->bc_type == 2) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unsupported boundary condition type <move base down to accomodate pure thickening : case 2>");			
		} else {
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unsupported boundary condition type %D",data->bc_type);
		}
		
		/* [B] Advect surface with fluid velocity, internal mesh geometry in x-z is advected with background strain-rate */
		ierr = UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX(dav,velocity,mesh_velocity,dt);CHKERRQ(ierr);
		
		ierr = VecDestroy(&mesh_velocity);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFoldingUpdate_RemeshResampleSurface"
PetscErrorCode MultilayerFoldingUpdate_RemeshResampleSurface(PetscReal AR,DM dav,Vec velocity,PetscReal dt,ModelMultilayerFoldingCtx *data,pTatinCtx c)
{
	Vec              mesh_velocity;
	PetscReal        AR_max = 20.0;
	DMDAVecTraverse3d_InterpCtx IntpCtx;
	PetscReal        Lx,Lz,Ox,Oz;
	PetscReal        gmin[3],gmax[3];
	static PetscBool remesh = PETSC_FALSE;
	PetscErrorCode   ierr;
	
    
	if (AR > AR_max) {
		PetscPrintf(PETSC_COMM_WORLD,"[[%s]] Activating remeshing \n", __FUNCT__);
		remesh = PETSC_TRUE;
	}
	
	if (remesh) {
		PetscPrintf(PETSC_COMM_WORLD,"Performing surface remeshing \n", __FUNCT__);
		/*
		 Reset position of mesh 
		 This is required as (i) AR is estimated on the deformed mesh
		 (ii) UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX() advects the free surface 
		 */
		ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,-dt);CHKERRQ(ierr);
		
		/* [A] create mesh advection velocity field in x-z */
		ierr = DMCreateGlobalVector(dav,&mesh_velocity);CHKERRQ(ierr);
		
		ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
		Lx = gmax[0] - gmin[0];		Ox = gmin[0];
		Lz = gmax[2] - gmin[2];   Oz = gmin[2];
		
		if (data->bc_type == 0) {
			PetscReal vx_E,vx_W;
			/* compression east/west in the x-direction (0) [east-west] using constant velocity */
			
			vx_E = -data->vx_compression;
			vx_W =  data->vx_compression;
			
			ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vx_E - vx_W)/(Lx),vx_W,0.0);CHKERRQ(ierr);
			ierr = DMDAVecTraverse3d(dav,mesh_velocity, 0, DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
			
			ierr = VecStrideSet(mesh_velocity,1,0.0);CHKERRQ(ierr); // zero y component //
			ierr = VecStrideSet(mesh_velocity,2,0.0);CHKERRQ(ierr); // zero z component //
		} else if (data->bc_type == 1) {
			PetscReal vx_E,vx_W;
			PetscReal vz_F,vz_B;
						
			/* x component */
			//vx_E = data->exx * ( gmax[0] - cx[0] );
			//vx_W = data->exx * ( gmin[0] - cx[0] );
			vx_E = 0.5 * data->exx * Lx;
			vx_W = -vx_E;
			PetscPrintf(PETSC_COMM_WORLD,"  x: pos(min/max) = %+1.4e / %+1.4e  vx = %+1.4e / %+1.4e \n",gmin[0],gmax[0],vx_W,vx_E);
			
			ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vx_E - vx_W)/(Lx),vx_W,Ox);CHKERRQ(ierr);
			ierr = DMDAVecTraverse3d(dav,mesh_velocity, 0, DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
			
			/* y component */
			ierr = VecStrideSet(mesh_velocity,1,0.0);CHKERRQ(ierr); // zero y component //
			
			/* z component */
			//vz_F = data->ezz * ( gmax[2] - cx[2] );
			//vz_B = data->ezz * ( gmin[2] - cx[2] );
			vz_F = 0.5 * data->ezz * Lz;
			vz_B = -vz_F;
			PetscPrintf(PETSC_COMM_WORLD,"  z: pos(min/max) = %+1.4e / %+1.4e  vz = %+1.4e / %+1.4e \n",gmin[2],gmax[2],vz_B,vz_F);
			
			ierr = DMDAVecTraverse3d_InterpCtxSetUp_Z(&IntpCtx,(vz_F - vz_B)/(Lz),vz_B,Oz);CHKERRQ(ierr);
			ierr = DMDAVecTraverse3d(dav,mesh_velocity, 2, DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
		} else if (data->bc_type == 2) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unsupported boundary condition type <move base down to accomodate pure thickening : case 2>");			
		} else {
			SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unsupported boundary condition type %D",data->bc_type);
		}
		
		/* [B] Advect surface with fluid velocity, internal mesh geometry in x-z is advected with background strain-rate */
		ierr = UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX(dav,velocity,mesh_velocity,dt);CHKERRQ(ierr);

		
		ierr = VecDestroy(&mesh_velocity);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_ModelApplyUpdateMeshGeometry_MultilayerFolding"
PetscErrorCode _ModelApplyUpdateMeshGeometry_MultilayerFolding(pTatinCtx c,Vec X,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	PetscReal      step;
	PhysCompStokes stokes;
	DM             stokes_pack,dav,dap;
	Vec            velocity,pressure;
	PetscInt           metric_L = 5; 
	MeshQualityMeasure metric_list[] = { MESH_QUALITY_ASPECT_RATIO, MESH_QUALITY_DISTORTION, MESH_QUALITY_DIAGONAL_RATIO, MESH_QUALITY_VERTEX_ANGLE, MESH_QUALITY_FACE_AREA_RATIO };
	PetscReal          value[100];
	PetscBool          basal_remesh    = PETSC_FALSE;
	PetscBool          marker_remesh_v1 = PETSC_FALSE;
	PetscBool          marker_remesh_v2 = PETSC_FALSE;
	PetscBool          surface_remesh = PETSC_FALSE;
	
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);

	/* advect the mesh with the full velocity field */
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	/* check mesh quality */
	ierr = DMDAComputeMeshQualityMetricList(dav,metric_L,metric_list,value);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_ASPECT_RATIO\"    %1.4e \n", value[0]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_DISTORTION\"      %1.4e \n", value[1]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_DIAGONAL_RATIO\"  %1.4e \n", value[2]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_VERTEX_ANGLE\"    %1.4e \n", value[3]);
	PetscPrintf(PETSC_COMM_WORLD,"  Mesh metrics \"MESH_QUALITY_FACE_AREA_RATIO\" %1.4e \n", value[4]);
	
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_multilayer_basal_remesh",&basal_remesh,PETSC_NULL);CHKERRQ(ierr);
	if (basal_remesh) {
		ierr = MultilayerFoldingUpdate_RemeshBasalLayer(value[0],dav,data);CHKERRQ(ierr);
	}
	
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_multilayer_marker_remesh",&marker_remesh_v1,PETSC_NULL);CHKERRQ(ierr);
	if (marker_remesh_v1) {
		ierr = MultilayerFoldingUpdate_RemeshMarkerProjection_v1(value[0],dav,c);CHKERRQ(ierr);
	}

	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_multilayer_marker_remesh_v2",&marker_remesh_v2,PETSC_NULL);CHKERRQ(ierr);
	if (marker_remesh_v2) {
		/* 
		 Reset position of mesh 
		 This is required as (i) AR is estimated on the deformed mesh
		                     (ii) UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX() advects the free surface 
		*/
		ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		
		ierr = MultilayerFoldingUpdate_RemeshMarkerProjection_v2(value[0],dav,velocity,step,data,c);CHKERRQ(ierr);

		ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	}
	
	ierr = PetscOptionsGetBool(PETSC_NULL,"-model_multilayer_surface_remesh",&surface_remesh,PETSC_NULL);CHKERRQ(ierr);
	if (surface_remesh) {
		/* 
		 Reset position of mesh 
		 This is required as (i) AR is estimated on the deformed mesh
		 (ii) UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX() advects the free surface 
		 */
		ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		
		ierr = MultilayerFoldingUpdate_RemeshResampleSurface(value[0],dav,velocity,step,data,c);CHKERRQ(ierr);
		
		ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialCondition_MultilayerFolding"
PetscErrorCode ModelInitialCondition_MultilayerFolding(pTatinCtx c,Vec X,void *ctx)
{
	/*
	 ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	 DM stokes_pack,dau,dap;
	 Vec velocity,pressure;
	 PetscReal rho0;
	 DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	 DMDAVecTraverse3d_InterpCtx IntpCtx;*/
	
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

/* 
 Fletcher, Three-dimensional folding of an embedded viscous layer in pure shear, 1991, eq 29 a
 */
double compute_q_3d(double H,double eta_l,double eta_m,double l,double m)
{
	double R,k,alpha,beta,beta1,q;
	double lambda;
	
	R = eta_m/eta_l;
	lambda = sqrt(l*l + m*m);
	k = lambda * H;
	
	alpha = -2.0 * (1.0 - R);
	beta  = 1.0 - R*R;
	beta1 = -1.0 * ( (1.0+R*R)*(exp(k)-exp(-k)) + 2.0*R*(exp(k)+exp(-k))  )/(2.0*k);
	q = alpha / ( beta + beta1 );
	
	return q;
}

double compute_q_2d(double H, double eta_l,double eta_m,double kx)
{
	double R,k,alpha,beta,beta1,q;
	
	R = eta_m/eta_l;
	k = kx*H;
	
	alpha = -2.0 * (1.0 - R);
	beta  = 1.0 - R*R;
	beta1 = -1.0 * ( (1.0+R)*(1.0+R)*exp(k) - (1.0-R)*(1.0-R)*exp(-k)  )/(2.0*k);
	q = alpha / ( beta + beta1 );
	
	return q;
}

double compute_amplitude(double A0,double H,double time,double exx,double ezz,double eta_l,double eta_m,double l,double m)
{
	double eyy,lambda,lambda2,arg,q,Anew;
	
	q = compute_q_3d(H,eta_l,eta_m,l,m);
	
	eyy = -(exx + ezz);
	lambda = sqrt(l*l + m*m);
	lambda2 = lambda*lambda;
	// solution of equation 29 from Fletcher, 1991 //
	arg = eyy - 0.5 * q * ( (l*l/lambda2)*exx + (m*m/lambda2)*ezz - eyy );
	Anew = A0 * exp( arg * time );
	
	return Anew;
}

double extract_q_cylindrical(double amp, double exx,double time, double A0)
{
	// eq 28 Fletcher 1977 dA_dt = (1+q)(-exx)*A where dA_dt can be approx by (A-A0)/time if time at t0 = 0 
	double q, alpha;
	
	alpha = (1.0-A0/amp)/time; 
	q = -(alpha/exx)-1.0; 
	
	return q; 
}

double extract_q_3d(double amp, double exx, double ezz, double l, double m, double time, double A0)
{
	// eq 29 Fletcher 1991 dA_dt = eyy*A-(q/2)[(l*l/lambda2)*exx + (m*m/lambda2)*ezz - eyy]*A     where dA_dt can be approx by (A-A0)/time if time at t0 = 0 
	double eyy, q, lambda, lambda2, arg, alpha, beta; 
	
	arg     =  (1.0-A0/amp)/time; 
	eyy     =  -(exx + ezz);
	lambda  =  sqrt(l*l + m*m); 
	lambda2 =  lambda*lambda;
	
	alpha = -2.0 * (arg-eyy);
	beta  = (l*l/lambda2)*exx + (m*m/lambda2)*ezz - eyy; 
	
	q = alpha/beta; 
	
	return q; 
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFoldingOutputAmplitudeMax"
PetscErrorCode MultilayerFoldingOutputAmplitudeMax(pTatinCtx c,ModelMultilayerFoldingCtx *data)
{
	PetscInt jinter,i,k,si,sj,sk,nx,ny,nz;
	DM dav,cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	PetscErrorCode ierr;
	PetscReal peak,xz[2],H;
	static int been_here = 0;
	
	
	if (been_here == 0) {
		PetscPrintf(PETSC_COMM_WORLD,"# kx %+1.4e \n",data->kx);
		PetscPrintf(PETSC_COMM_WORLD,"# kz %+1.4e \n",data->kz);
		
		if (data->bc_type == 0) {
			PetscPrintf(PETSC_COMM_WORLD,"# vx %+1.4e \n",data->vx_compression);
			PetscPrintf(PETSC_COMM_WORLD,"# vz %+1.4e \n",data->vz_compression);
		} else if (data->bc_type == 1) {
			PetscPrintf(PETSC_COMM_WORLD,"# exx %+1.4e \n",data->exx);
			PetscPrintf(PETSC_COMM_WORLD,"# ezz %+1.4e \n",data->ezz);
		} else if (data->bc_type == 2) {
			PetscPrintf(PETSC_COMM_WORLD,"# exx %+1.4e \n",data->exx);
			PetscPrintf(PETSC_COMM_WORLD,"# ezz %+1.4e \n",data->ezz);
		}
		
		PetscPrintf(PETSC_COMM_WORLD,"# H  %+1.4e \n", data->interface_heights[2]-data->interface_heights[1]);
		PetscPrintf(PETSC_COMM_WORLD,"# A0 %+1.4e \n", data->A0);
		
		PetscPrintf(PETSC_COMM_WORLD,"# interface3 %d [node index]\n", 2*(data->layer_res_j[0]+data->layer_res_j[1]+data->layer_res_j[2]));
		PetscPrintf(PETSC_COMM_WORLD,"# interface2 %d \n", 2*(data->layer_res_j[0]+data->layer_res_j[1]));
		PetscPrintf(PETSC_COMM_WORLD,"# interface1 %d \n", 2*data->layer_res_j[0]);
		PetscPrintf(PETSC_COMM_WORLD,"# interface0 %d \n", 0);
		been_here = 1;
	}
	
	H = data->interface_heights[2]-data->interface_heights[1];
	
	jinter  = 2 * data->layer_res_j[0];
	jinter += 2 * data->layer_res_j[1];
	
	dav    = c->stokes_ctx->dav;
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	peak = -1.0e32;
	xz[0] = xz[1] = -1.0e32;
	for (i=si; i<si+nx; i++) {
		for (k=sk; k<sk+nz; k++) {
			
			if (LA_coord[k][jinter][i].y > peak) {
				peak = LA_coord[k][jinter][i].y;
				xz[0] = LA_coord[k][jinter][i].x;
				xz[1] = LA_coord[k][jinter][i].z;
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"# A_max %1.8e [x,z %1.4e,%1.4e] time %1.4e \n", peak,xz[0],xz[1],c->time);
	{
		double A_init,A_anl;
		
		A_init = 0.5 * H + data->A0 * cos(data->kx * xz[0]) * cos(data->kz * xz[1]);
		A_anl = compute_amplitude(A_init,H,c->time,data->exx,data->ezz,data->eta[1],data->eta[0],data->kx,data->kz);
		PetscPrintf(PETSC_COMM_WORLD,"# ** A_anl ** %1.8e : variation %1.4f%% \n",A_anl,100.0*(peak-A_anl)/A_anl);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFoldingOutput_q"
PetscErrorCode MultilayerFoldingOutput_q(pTatinCtx c,ModelMultilayerFoldingCtx *data)
{
	PetscInt jinter,i,k,si,sj,sk,nx,ny,nz;
	DM dav,cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	PetscErrorCode ierr;
	PetscReal peak_max,peak_min,g_peak_max,g_peak_min,xz[2],yz[2],H;
	static int been_here = 0;
	PetscReal arg,arg_3d,q_cylindrical,q_3d,q_th_3D,q_th_2D,err,amp;
	
    
	if (been_here == 0) {
		PetscPrintf(PETSC_COMM_WORLD,"# kx %+1.4e \n",data->kx);
		PetscPrintf(PETSC_COMM_WORLD,"# kz %+1.4e \n",data->kz);
		
		if (data->bc_type == 0) {
			PetscPrintf(PETSC_COMM_WORLD,"# vx %+1.4e \n",data->vx_compression);
			PetscPrintf(PETSC_COMM_WORLD,"# vz %+1.4e \n",data->vz_compression);
		} else if (data->bc_type == 1) {
			PetscPrintf(PETSC_COMM_WORLD,"# exx %+1.4e \n",data->exx);
			PetscPrintf(PETSC_COMM_WORLD,"# ezz %+1.4e \n",data->ezz);
		} else if (data->bc_type == 2) {
			PetscPrintf(PETSC_COMM_WORLD,"# exx %+1.4e \n",data->exx);
			PetscPrintf(PETSC_COMM_WORLD,"# ezz %+1.4e \n",data->ezz);
		}
		
		PetscPrintf(PETSC_COMM_WORLD,"# H  %+1.4e \n", data->interface_heights[2]-data->interface_heights[1]);
		PetscPrintf(PETSC_COMM_WORLD,"# A0 %+1.4e \n", data->A0);
		
		PetscPrintf(PETSC_COMM_WORLD,"# interface3 %d [node index]\n", 2*(data->layer_res_j[0]+data->layer_res_j[1]+data->layer_res_j[2]));
		PetscPrintf(PETSC_COMM_WORLD,"# interface2 %d \n", 2*(data->layer_res_j[0]+data->layer_res_j[1]));
		PetscPrintf(PETSC_COMM_WORLD,"# interface1 %d \n", 2*data->layer_res_j[0]);
		PetscPrintf(PETSC_COMM_WORLD,"# interface0 %d \n", 0);
		been_here = 1;
	}
	
	H = data->interface_heights[2]-data->interface_heights[1];
	
	jinter  = 2 * data->layer_res_j[0];
	jinter += 2 * data->layer_res_j[1];
	
	dav    = c->stokes_ctx->dav;
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	peak_max = -1.0e32;
	peak_min = 1.0e32;
	xz[0] = xz[1] = -1.0e32;
	yz[0] = yz[1] = -1.0e32;
	for (i=si; i<si+nx; i++) {
		for (k=sk; k<sk+nz; k++) {
			
			if (LA_coord[k][jinter][i].y > peak_max) {
				peak_max = LA_coord[k][jinter][i].y;
				xz[0] = LA_coord[k][jinter][i].x;
				xz[1] = LA_coord[k][jinter][i].z;
			}
			
			if (LA_coord[k][jinter][i].y < peak_min) {
				peak_min = LA_coord[k][jinter][i].y;
				yz[0] = LA_coord[k][jinter][i].x;
				yz[1] = LA_coord[k][jinter][i].z;
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&peak_min,&g_peak_min,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&peak_max,&g_peak_max,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);

	
	PetscPrintf(PETSC_COMM_WORLD,"# A_max %1.8e [x,z %1.4e,%1.4e] time %1.4e \n", g_peak_max,xz[0],xz[1],c->time);
	PetscPrintf(PETSC_COMM_WORLD,"# A_min %1.8e [x,z %1.4e,%1.4e] time %1.4e \n", g_peak_min,yz[0],yz[1],c->time);
	
	amp = 0.5 * ( g_peak_max - g_peak_min ); 
	
	/* for the case of the 2d cylindrical problem dA_dt = (1+q)(-exx)*A eq(28) */
	q_cylindrical = extract_q_cylindrical(amp, data->exx, c->time, data->A0);
	
    arg = arg_3d = -1.0e32;
	q_3d = extract_q_3d(amp, data->exx, data->ezz, data->kx, data->kz, c->time, data->A0);
	PetscPrintf(PETSC_COMM_WORLD,"amp %1.8e, arg %1.8e, arg_3d %1.8e \n",amp, arg, arg_3d); 
		
	q_th_3D = compute_q_3d(H, data->eta[1], data->eta[0], data->kx, data->kz);
	q_th_2D = compute_q_2d(H, data->eta[1], data->eta[0],data->kx); 
	PetscPrintf(PETSC_COMM_WORLD,"# q_cylindrical %1.8e, q_3d %1.8e, q_th_3D %1.8e, q_th_2D %1.8e \n", q_cylindrical, q_3d, q_th_3D, q_th_2D);

	err = fabs( 100.0*(q_3d - q_th_3D)/q_th_3D );
	PetscPrintf(PETSC_COMM_WORLD,"# ** q_th: error %1.4f%% \n",err);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultilayerFoldingOutputAmplitude"
PetscErrorCode MultilayerFoldingOutputAmplitude(pTatinCtx c,ModelMultilayerFoldingCtx *data)
{
	PetscInt jinter,i,k,si,sj,sk,nx,ny,nz;
	DM dav,cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	PetscErrorCode ierr;
	PetscReal peak,xz[2],H;
	static int been_here = 0;
	double error[2],gerror[2];
	
	
	if (been_here == 0) {
		PetscPrintf(PETSC_COMM_WORLD,"# kx %+1.4e \n",data->kx);
		PetscPrintf(PETSC_COMM_WORLD,"# kz %+1.4e \n",data->kz);
		
		if (data->bc_type == 0) {
			PetscPrintf(PETSC_COMM_WORLD,"# vx %+1.4e \n",data->vx_compression);
			PetscPrintf(PETSC_COMM_WORLD,"# vz %+1.4e \n",data->vz_compression);
		} else if (data->bc_type == 1) {
			PetscPrintf(PETSC_COMM_WORLD,"# exx %+1.4e \n",data->exx);
			PetscPrintf(PETSC_COMM_WORLD,"# ezz %+1.4e \n",data->ezz);
		} else if (data->bc_type == 2) {
			PetscPrintf(PETSC_COMM_WORLD,"# exx %+1.4e \n",data->exx);
			PetscPrintf(PETSC_COMM_WORLD,"# ezz %+1.4e \n",data->ezz);
		}
		
		PetscPrintf(PETSC_COMM_WORLD,"# H  %+1.4e \n", data->interface_heights[2]-data->interface_heights[1]);
		PetscPrintf(PETSC_COMM_WORLD,"# A0 %+1.4e \n", data->A0);
		
		PetscPrintf(PETSC_COMM_WORLD,"# interface3 %d [node index]\n", 2*(data->layer_res_j[0]+data->layer_res_j[1]+data->layer_res_j[2]));
		PetscPrintf(PETSC_COMM_WORLD,"# interface2 %d \n", 2*(data->layer_res_j[0]+data->layer_res_j[1]));
		PetscPrintf(PETSC_COMM_WORLD,"# interface1 %d \n", 2*data->layer_res_j[0]);
		PetscPrintf(PETSC_COMM_WORLD,"# interface0 %d \n", 0);
		been_here = 1;
	}
	
	H = data->interface_heights[2]-data->interface_heights[1];
	
	jinter  = 2 * data->layer_res_j[0];
	jinter += 2 * data->layer_res_j[1];
	
	dav    = c->stokes_ctx->dav;
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	error[0] = 1.0e32;
	error[1] = 0.0;
	for (i=si; i<si+nx; i++) {
		for (k=sk; k<sk+nz; k++) {
			double A_init,A_anl,err;
			
			peak = LA_coord[k][jinter][i].y;
			xz[0] = LA_coord[k][jinter][i].x;
			xz[1] = LA_coord[k][jinter][i].z;
			
			A_init = 0.5 * H + data->A0 * cos(data->kx * xz[0]) * cos(data->kz * xz[1]);
			A_anl = compute_amplitude(A_init,H,c->time,data->exx,data->ezz,data->eta[1],data->eta[0],data->kx,data->kz);
			err = fabs( 100.0*(peak-A_anl)/A_anl );
			
			error[0] = fmin(error[0],err);
			error[1] = fmax(error[1],err);
		}
	}
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&error[0],&gerror[0],1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&error[1],&gerror[1],1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"# ** A_anl: errors [min/max] %1.4f%% %1.4f%% \n",gerror[0],gerror[1]);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_MultilayerFolding"
PetscErrorCode ModelOutput_MultilayerFolding(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	DataBucket     materialpoint_db;
	PetscBool      verify_with_analytics = PETSC_FALSE;
	PetscBool      output_markers = PETSC_FALSE;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);

	PetscOptionsGetBool(NULL,"-model_multilayer_folding_output_markers",&output_markers,0);
	/* output raw marker fields */
	if (output_markers) {
		const int                 nf = 3;
		const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
		char                      name[PETSC_MAX_PATH_LEN];

		sprintf(name,"%s_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,name);CHKERRQ(ierr);
	}
	
	/* output marker->cell projected fields */
	{
		const int                   nf = 2;
		const MaterialPointVariable mp_prop_list[] = { MPV_viscosity, MPV_density }; 
		
		ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
	}
	
	PetscOptionsGetBool(NULL,"-verify_with_analytics",&verify_with_analytics,NULL);
	if (verify_with_analytics) {
		//ierr = MultilayerFoldingOutputAmplitudeMax(c,data);CHKERRQ(ierr);
		//ierr = MultilayerFoldingOutputAmplitude(c,data);CHKERRQ(ierr);
		ierr = MultilayerFoldingOutput_q(c,data); CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_MultilayerFolding"
PetscErrorCode ModelDestroy_MultilayerFolding(pTatinCtx c,void *ctx)
{
	ModelMultilayerFoldingCtx *data = (ModelMultilayerFoldingCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_MultilayerFolding"
PetscErrorCode pTatinModelRegister_MultilayerFolding(void)
{
	ModelMultilayerFoldingCtx *data;
	PetscErrorCode ierr;
	pTatinModel    m;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelMultilayerFoldingCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelMultilayerFoldingCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"multilayer_folding");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_MultilayerFolding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_MultilayerFolding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_MultilayerFolding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_MultilayerFolding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_MultilayerFolding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_MultilayerFolding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))_ModelApplyUpdateMeshGeometry_MultilayerFolding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_MultilayerFolding);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_MultilayerFolding);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

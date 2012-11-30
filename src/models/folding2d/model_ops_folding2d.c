
#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_models.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"

#include "model_folding2d_ctx.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Folding2d"
PetscErrorCode ModelInitialize_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* assign defaults */
	data->eta1 = 1.0;
	data->eta2 = 400.0;
	data->eta3 = 1.0;

	data->rho1 = data->rho2 = data->rho3 = 1.0;
	data->exx = 1.0e-3;
	
	/* parse from command line or input file */
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_eta1",&data->eta1,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_eta2",&data->eta2,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_eta3",&data->eta3,&flg);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_rho1",&data->rho1,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_rho2",&data->rho2,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_rho3",&data->rho3,&flg);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_exx",&data->exx,&flg);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Folding2d"
PetscErrorCode ModelApplyBoundaryCondition_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscReal         exx;
	BCList            bclist;
	DM                dav;
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	
	exx = data->exx;

	bclist = c->stokes_ctx->u_bclist;
	dav    = c->stokes_ctx->dav;
	
	/* extension/compression east/west in the x-direction (0) */
	ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,exx,0);CHKERRQ(ierr);
	
	/* free slip south */
	ierr = DirichletBC_FreeSlip(bclist,dav,SOUTH_FACE);CHKERRQ(ierr);
	
	/* free surface north */
	/* do nothing! */
	
	/* free slip front/back to mimic 2d behaviour */
	ierr = DirichletBC_FreeSlip(bclist,dav,FRONT_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,BACK_FACE);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Folding2d"
PetscErrorCode ModelApplyBoundaryConditionMG_Folding2d(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Folding2d"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Folding2dPerturbInterfaces"
PetscErrorCode Folding2dPerturbInterfaces(DM dav,PetscReal dy,PetscInt *jb,PetscInt *jt)
{
	PetscErrorCode ierr;
	PetscReal Ly1,Ly2;
	PetscInt j_max_layer_top, j_max_layer_bottom;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz,M,N,P;
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	PetscReal MeshMin[3],MeshMax[3];
	PetscReal min_sep_bottom,min_sep_top,sep,amp;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	Ly1 = 1.0; /* height of layer (bottom) */
	Ly2 = 2.0; /* height of layer (top) */


	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	j_max_layer_bottom = -1;
	j_max_layer_top    = -1;
	min_sep_bottom = 1.0e32;
	min_sep_top = 1.0e32;
	
	ierr = DMDAGetLocalBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);	
	
	if ( (Ly1 >= MeshMin[1]) && (Ly1 < MeshMax[1]) ) {
		k = sk;
		i = si;
		for (j=sj; j<sj+ny; j++) {
			PetscReal x,y,z;
			
			x = LA_coord[k][j][i].x;
			y = LA_coord[k][j][i].y;
			z = LA_coord[k][j][i].z;
			
			sep = sqrt( (Ly1-y) * (Ly1-y) );
			if (sep < min_sep_bottom) {
				min_sep_bottom = sep;
				j_max_layer_bottom = j;
			}
		}
	}

	if ( (Ly2 >= MeshMin[1]) && (Ly2 < MeshMax[1]) ) {
		k = sk;
		i = si;
		for (j=sj; j<sj+ny; j++) {
			PetscReal x,y,z;
			
			x = LA_coord[k][j][i].x;
			y = LA_coord[k][j][i].y;
			z = LA_coord[k][j][i].z;
			
			sep = sqrt( (Ly2-y) * (Ly2-y) );
			if (sep < min_sep_top) {
				min_sep_top = sep;
				j_max_layer_top= j;
			}
		}
	}
	

	/* NOTES:
	 Reset set on each z level.
	 Use a different seed for upper and lower layers
	 */
	amp = 1.0e-1 * 1; /* make this 5 to see the perturb visible for testing */
	if (j_max_layer_bottom != -1) {
		j = j_max_layer_bottom;
		
		for (k=sk; k<sk+nz; k++) {
			srand(0);
			for (i=0; i<M; i++) {	
				PetscReal y;
				double random;

				random = 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
				if ( (i>=si) && (i<si+nx) ) {
					y = LA_coord[k][j][i].y;
					y = y + amp * dy * random;
					LA_coord[k][j][i].y = y;
				}
			}
		}
	}
	

	if (j_max_layer_top != -1) {
		j = j_max_layer_top;
		
		for (k=sk; k<sk+nz; k++) {
			srand(1);
			for (i=0; i<M; i++) {
				PetscReal y;
				double random;
				
				random = 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
				if ( (i>=si) && (i<si+nx) ) {
					y = LA_coord[k][j][i].y;
					y = y + amp * dy * random;
					LA_coord[k][j][i].y = y;
				}
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);

	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	
	*jb = j_max_layer_bottom;
	*jt = j_max_layer_top;

	printf("Ly1 = %1.4e => j = %d \n", Ly1,*jb);
	printf("Ly2 = %1.4e => j = %d \n", Ly2,*jt);

	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Folding2d"
PetscErrorCode ModelApplyInitialMeshGeometry_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscReal Lx,Ly,dx,dy,dz,Lz;
	PetscInt mx,my,mz;
	PetscInt j_int1,j_int2;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* step 1 - create structured grid */
	Lx = 50.0;
	Ly = 3.0;
	
	/*
	 The length of the model in z-direction is determined by the grid spacing in x and y.
	 We do this so that the elements do not have a large aspect ratio. This would occur
	 if we hard coded Lz to a constant number which is independnet of the grid resolution in x and y.
	 We choose Lz to be mz * min(dx,dy).
	 */
	mx = c->mx; 
	my = c->my; 
	mz = c->mz; 
	
	dx = Lx / ((PetscReal)mx);
	dy = Ly / ((PetscReal)my);
	dz = dx;
	if (dz < dy) {
		dz = dy;
	}
	Lz = mz * dz;
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, 0.0,Lx, 0.0,Ly, 0.0,Lz);CHKERRQ(ierr);

	
	/* step 2 - define two interfaces and perturb coords along the interface */
	ierr = Folding2dPerturbInterfaces(c->stokes_ctx->dav,dy, &j_int1,&j_int2);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Folding2d"
PetscErrorCode ModelApplyInitialMaterialGeometry_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Folding2d"
PetscErrorCode ModelApplyUpdateMeshGeometry_Folding2d(pTatinCtx c,Vec X,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Folding2d"
PetscErrorCode ModelOutput_Folding2d(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Folding2d"
PetscErrorCode ModelDestroy_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Folding2d"
PetscErrorCode pTatinModelRegister_Folding2d(void)
{
	ModelFolding2dCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelFolding2dCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelFolding2dCtx));CHKERRQ(ierr);
	
	/* set initial values for model parameters */
	data->eta1 = data->eta2 = data->eta3 = 0.0;
	data->rho1 = data->rho2 = data->rho3 = 0.0;
	data->exx = 0.0;
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"folding2d");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Folding2d);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

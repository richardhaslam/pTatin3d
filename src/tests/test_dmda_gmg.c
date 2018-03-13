/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   test_dmda_gmg.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

static const char help[] = "Stokes solver using Q2-Pm1 mixed finite elements.\n"
"3D prototype of the (p)ragmatic version of Tatin. (pTatin3d_v0.0)\n\n";

#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscdm.h"
#include "petsc/private/dmdaimpl.h" 

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "stokes_form_function.h"
#include "dmda_project_coords.h"
#include "dmda_element_q2p1.h"
#include "stokes_operators.h"


#undef __FUNCT__  
#define __FUNCT__ "FormJacobian_Stokes"
PetscErrorCode FormJacobian_Stokes(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
	
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}
	

#undef __FUNCT__  
#define __FUNCT__ "test_pTatin3d_gmg_galerkin"
PetscErrorCode test_pTatin3d_gmg_galerkin(int argc,char **argv)
{
	DM             multipys_pack,dav;
	DM             *dav_hierarchy;
	Mat            *interpolatation;
	pTatinCtx      user;
	PetscInt       nlevels,k;
	PetscMPIInt    rank;
	Vec            X,F;
	Mat            A,B;
	SNES           snes;
	KSP            ksp;
	PC             pc;
	IS             *isg;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);

	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);

	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	user->pack = user->stokes_ctx->stokes_pack;

	/* fetch some local variables */
	multipys_pack = user->pack;
	dav           = user->stokes_ctx->dav;
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	/* interpolate point coordinates (needed if mesh was modified) */
	//ierr = QuadratureStokesCoordinateSetUp(user->stokes_ctx->Q,dav);CHKERRQ(ierr);
	//for (e=0; e<QUAD_EDGES; e++) {
	//	ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
	//}
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);


	/* set up mg */
	user->stokes_ctx->dav->ops->coarsenhierarchy = DMCoarsenHierarchy2_DA;
	
	nlevels = 1;
	PetscOptionsGetInt(NULL,NULL,"-dau_nlevels",&nlevels,0);
	ierr = PetscMalloc(sizeof(DM)*nlevels,&dav_hierarchy);CHKERRQ(ierr);
	dav_hierarchy[ nlevels-1 ] = dav;
	ierr = PetscObjectReference((PetscObject)dav);CHKERRQ(ierr);
	
	/* option 1 - simply coarsen nlevels - 1 times */
	{
		DM *coarsened_list;
		ierr = PetscMalloc(sizeof(DM)*(nlevels-1),&coarsened_list);CHKERRQ(ierr);
		ierr = DMCoarsenHierarchy(dav,nlevels-1,coarsened_list);CHKERRQ(ierr);
		for (k=0; k<nlevels-1; k++) {
			dav_hierarchy[ nlevels-2-k ] = coarsened_list[k];
		}
		PetscFree(coarsened_list);
/*		
		for (k=0; k<nlevels; k++) {
			PetscPrintf(PETSC_COMM_WORLD,"level[%d] :\n", k );
			ierr = DMView(dav_hierarchy[k],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		}
*/ 
	}

	/* set to be q2 */
	for (k=0; k<nlevels-1; k++) {
		ierr = DMDASetElementType_Q2(dav_hierarchy[k]);CHKERRQ(ierr);
	}
	
	/* test */
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	for (k=0; k<nlevels; k++) {
		PetscInt nels,nen;
		const PetscInt *els;
		PetscInt lmx,lmy,lmz,si,sj,sk;
		
		ierr = DMDAGetElements_pTatinQ2P1(dav_hierarchy[k],&nels,&nen,&els);CHKERRQ(ierr);
		

		ierr = DMDAGetSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"         level [%2D]: global Q2 elements (%D x %D x %D) \n", k,lmx,lmy,lmz );

		ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: local Q2 elements  (%D x %D x %D) \n", rank, k,lmx,lmy,lmz );

		ierr = DMDAGetCornersElementQ2(dav_hierarchy[k],&si,&sj,&sk,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		si = si/2;
		sj = sj/2;
		sk = sk/2;
		PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: element range [%D - %D] x [%D - %D] x [%D - %D] \n", rank, k,si,si+lmx-1,sj,sj+lmy-1,sk,sk+lmz-1 );
		
	}
	
	/* inject coordinates */
	ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);
	
	/* define interpolation operators */
	ierr = PetscMalloc(sizeof(Mat)*nlevels,&interpolatation);CHKERRQ(ierr);
	interpolatation[0] = NULL;
	for (k=0; k<nlevels-1; k++) {
		ierr = DMCreateInterpolation(dav_hierarchy[k],dav_hierarchy[k+1],&interpolatation[k+1],NULL);CHKERRQ(ierr);
	}
	
#if 0
	/* define boundary conditions */
	for (k=0; k<nlevels-1; k++) {
		ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);
	}
	
	
	/* define material properties on gauss points */
	for (k=0; k<nlevels-1; k++) {
		PetscInt ncells,lmx,lmy,lmz;
		PetscInt np_per_dim;
		
		np_per_dim = 3;
		ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		ncells = lmx * lmy * lmz;
		ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&volQ[k]);CHKERRQ(ierr);
	}
#endif	
	
	
	/* define operators */
	ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
	/* I cheat here - I assemble A11 - OUCH */
	ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,0,1,1,&B);CHKERRQ(ierr);
	
	/* Basic solver configuration using SNES - FIELDSPLIT */
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
	
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	ierr = SNESSetDM(snes,multipys_pack);CHKERRQ(ierr);
	ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  

//	if (user->use_mf_stokes) {
//		ierr = SNESSetJacobian(snes,B,B,FormJacobian_Stokes,user);CHKERRQ(ierr);
//	} else {
		ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);
//	}
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

	/* configure for fieldsplit */
	ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	
	ierr = DMCompositeGetGlobalISs(multipys_pack,&isg);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc,"u",isg[0]);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc,"p",isg[1]);CHKERRQ(ierr);
	
	/* configure uu split for galerkin multi-grid */
	{
		PetscInt nsplits;
		PC       pc_i;
		KSP      *sub_ksp;
		
		ierr = KSPSetUp(ksp);CHKERRQ(ierr);
		ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
		
		ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
		ierr = PCSetType(pc_i,PCMG);CHKERRQ(ierr);
		ierr = PCMGSetLevels(pc_i,nlevels,NULL);CHKERRQ(ierr);
		ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
		ierr = PCMGSetGalerkin(pc_i,PETSC_TRUE);CHKERRQ(ierr); /* OUCH - GALERKIN */
		
		for( k=1; k<nlevels; k++ ){
			ierr = PCMGSetInterpolation(pc_i,k,interpolatation[k]);CHKERRQ(ierr);
		}
		
	}
	
	ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
	
	
	
	
	
	ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	ierr = ISDestroy(&isg[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&isg[1]);CHKERRQ(ierr);
	ierr = PetscFree(isg);CHKERRQ(ierr);
	
	
	
	for (k=1; k<nlevels; k++) {
		ierr = MatDestroy(&interpolatation[k]);CHKERRQ(ierr);
	}
	ierr = PetscFree(interpolatation);CHKERRQ(ierr);
	
	for (k=0; k<nlevels; k++) {
		ierr = DMDestroy(&dav_hierarchy[k]);CHKERRQ(ierr);
	}
	ierr = PetscFree(dav_hierarchy);CHKERRQ(ierr);
	
	
	
	ierr = pTatin3dDestroyContext(&user);

	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "test_pTatin3d_gmg_mf"
PetscErrorCode test_pTatin3d_gmg_mf(int argc,char **argv)
{
	DM             multipys_pack,dav;
	DM             *dav_hierarchy;
	Mat            *interpolatation;
	pTatinCtx      user;
	PetscInt       nlevels,k;
	Quadrature     volQ[10];
	BCList         u_bclist[10];
	PetscMPIInt    rank;
	Vec            X,F;
	Mat            A,B,A11MF[10];
	SNES           snes;
	KSP            ksp;
	PC             pc;
	IS             *isg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
	
	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
	
	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	user->pack = user->stokes_ctx->stokes_pack;
	
	/* fetch some local variables */
	multipys_pack = user->pack;
	dav           = user->stokes_ctx->dav;
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	/* interpolate point coordinates (needed if mesh was modified) */
	//ierr = QuadratureStokesCoordinateSetUp(user->stokes_ctx->Q,dav);CHKERRQ(ierr);
	//for (e=0; e<QUAD_EDGES; e++) {
	//	ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
	//}
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);
	
	
	/* set up mg */
	user->stokes_ctx->dav->ops->coarsenhierarchy = DMCoarsenHierarchy2_DA;
	
	nlevels = 1;
	PetscOptionsGetInt(NULL,NULL,"-dau_nlevels",&nlevels,0);
	ierr = PetscMalloc(sizeof(DM)*nlevels,&dav_hierarchy);CHKERRQ(ierr);
	dav_hierarchy[ nlevels-1 ] = dav;
	ierr = PetscObjectReference((PetscObject)dav);CHKERRQ(ierr);
	
	/* option 1 - simply coarsen nlevels - 1 times */
	{
		DM *coarsened_list;
		ierr = PetscMalloc(sizeof(DM)*(nlevels-1),&coarsened_list);CHKERRQ(ierr);
		ierr = DMCoarsenHierarchy(dav,nlevels-1,coarsened_list);CHKERRQ(ierr);
		for (k=0; k<nlevels-1; k++) {
			dav_hierarchy[ nlevels-2-k ] = coarsened_list[k];
		}
		PetscFree(coarsened_list);
		/*		
		 for (k=0; k<nlevels; k++) {
		 PetscPrintf(PETSC_COMM_WORLD,"level[%d] :\n", k );
		 ierr = DMView(dav_hierarchy[k],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		 }
		 */ 
	}
	
	/* set to be q2 */
	for (k=0; k<nlevels-1; k++) {
		ierr = DMDASetElementType_Q2(dav_hierarchy[k]);CHKERRQ(ierr);
	}
	
	/* test */
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	for (k=0; k<nlevels; k++) {
		PetscInt nels,nen;
		const PetscInt *els;
		PetscInt lmx,lmy,lmz,si,sj,sk;
		
		ierr = DMDAGetElements_pTatinQ2P1(dav_hierarchy[k],&nels,&nen,&els);CHKERRQ(ierr);
		
		
		ierr = DMDAGetSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"         level [%2D]: global Q2 elements (%D x %D x %D) \n", k,lmx,lmy,lmz );
		
		ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: local Q2 elements  (%D x %D x %D) \n", rank, k,lmx,lmy,lmz );
		
		ierr = DMDAGetCornersElementQ2(dav_hierarchy[k],&si,&sj,&sk,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		si = si/2;
		sj = sj/2;
		sk = sk/2;
		PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: element range [%D - %D] x [%D - %D] x [%D - %D] \n", rank, k,si,si+lmx-1,sj,sj+lmy-1,sk,sk+lmz-1 );
		
	}
	
	/* inject coordinates */
	ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);
	
	/* define interpolation operators */
	ierr = PetscMalloc(sizeof(Mat)*nlevels,&interpolatation);CHKERRQ(ierr);
	interpolatation[0] = NULL;
	for (k=0; k<nlevels-1; k++) {
		ierr = DMCreateInterpolation(dav_hierarchy[k],dav_hierarchy[k+1],&interpolatation[k+1],NULL);CHKERRQ(ierr);
	}
	
	/* define boundary conditions */
	for (k=0; k<nlevels-1; k++) {
		ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);
	}
	u_bclist[nlevels-1] = user->stokes_ctx->u_bclist;
	
	
	/* define material properties on gauss points */
	for (k=0; k<nlevels-1; k++) {
		PetscInt ncells,lmx,lmy,lmz;
		PetscInt np_per_dim;
		
		np_per_dim = 3;
		ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		ncells = lmx * lmy * lmz;
		ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&volQ[k]);CHKERRQ(ierr);
	}
	volQ[nlevels-1] = user->stokes_ctx->volQ;
	
	
	/* define operators */
	ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
	/* FINE GRID */
	ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,1,1,1,&B);CHKERRQ(ierr);

	for (k=0; k<nlevels-1; k++) {
		MatA11MF ctx;
		
		ierr = MatA11MFCreate(&ctx);CHKERRQ(ierr);
		ierr = MatA11MFSetup(ctx,dav_hierarchy[k],volQ[k],u_bclist[k]);CHKERRQ(ierr);
		
		ierr = StokesQ2P1CreateMatrix_MFOperator_A11(ctx,&A11MF[k]);CHKERRQ(ierr);
		ierr = MatA11MFDestroy(&ctx);CHKERRQ(ierr);
	}
	A11MF[nlevels-1] = NULL;
	
	ierr = DMCompositeGetGlobalISs(multipys_pack,&isg);CHKERRQ(ierr);
	/* Fetch from the nest */
	{
		Mat sub_A11;
		
		ierr = MatCreateSubMatrix(B,isg[0],isg[0],MAT_INITIAL_MATRIX,&sub_A11);CHKERRQ(ierr);
		A11MF[nlevels-1] = sub_A11;
		ierr = MatDestroy(&sub_A11);CHKERRQ(ierr);

		ierr = PetscObjectReference((PetscObject)A11MF[nlevels-1]);CHKERRQ(ierr);
	}
	
	
	
	/* Basic solver configuration using SNES - FIELDSPLIT */
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

	ierr = VecSetRandom(F,NULL);CHKERRQ(ierr);
	
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	ierr = SNESSetDM(snes,multipys_pack);CHKERRQ(ierr);
	ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
	
	//	if (user->use_mf_stokes) {
	//		ierr = SNESSetJacobian(snes,B,B,FormJacobian_Stokes,user);CHKERRQ(ierr);
	//	} else {
	ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);
	//	}
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
	
	/* configure for fieldsplit */
	ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

	ierr = PCFieldSplitSetIS(pc,"u",isg[0]);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc,"p",isg[1]);CHKERRQ(ierr);
	
	
	/* configure uu split for galerkin multi-grid */
	{
		PetscInt nsplits;
		PC       pc_i;
		KSP      *sub_ksp,ksp_smoother;
		
		ierr = KSPSetUp(ksp);CHKERRQ(ierr);
		ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
		
		ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
		ierr = PCSetType(pc_i,PCMG);CHKERRQ(ierr);
		ierr = PCMGSetLevels(pc_i,nlevels,NULL);CHKERRQ(ierr);
		ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
		ierr = PCMGSetGalerkin(pc_i,PETSC_FALSE);CHKERRQ(ierr);
		ierr = PCSetDM(pc_i,NULL);CHKERRQ(ierr);
		
		for( k=1; k<nlevels; k++ ){
			ierr = PCMGSetInterpolation(pc_i,k,interpolatation[k]);CHKERRQ(ierr);
		}

		for( k=1; k<nlevels; k++ ){
			ierr = PCMGSetInterpolation(pc_i,k,interpolatation[k]);CHKERRQ(ierr);
		}
		
		ierr = PCMGGetCoarseSolve(pc_i,&ksp_smoother);CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp_smoother,A11MF[0],A11MF[0]);CHKERRQ(ierr);
		for( k=1; k<nlevels; k++ ){
			ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);
			ierr = KSPSetOperators(ksp_smoother,A11MF[k],A11MF[k]);CHKERRQ(ierr);
		}
		
	}
	
	ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
	
	
	
	
	
	ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	ierr = ISDestroy(&isg[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&isg[1]);CHKERRQ(ierr);
	ierr = PetscFree(isg);CHKERRQ(ierr);

	for (k=0; k<nlevels-1; k++) {
		ierr = BCListDestroy(&u_bclist[k]);CHKERRQ(ierr);
	}
	for (k=0; k<nlevels-1; k++) {
		ierr = QuadratureDestroy(&volQ[k]);CHKERRQ(ierr);
	}
	
	
	for (k=0; k<nlevels; k++) {
		ierr = MatDestroy(&A11MF[k]);CHKERRQ(ierr);
	}
	
	for (k=1; k<nlevels; k++) {
		ierr = MatDestroy(&interpolatation[k]);CHKERRQ(ierr);
	}
	ierr = PetscFree(interpolatation);CHKERRQ(ierr);
	
	for (k=0; k<nlevels; k++) {
		ierr = DMDestroy(&dav_hierarchy[k]);CHKERRQ(ierr);
	}
	ierr = PetscFree(dav_hierarchy);CHKERRQ(ierr);
	
	
	
	ierr = pTatin3dDestroyContext(&user);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "test_putatin"
PetscErrorCode test_putatin(int argc,char **argv)
{
	DM             dav;
	pTatinCtx      user;
	Mat            A,B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
	
	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
	
	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	user->pack = user->stokes_ctx->stokes_pack;
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	
	/* fetch some local variables */
	dav           = user->stokes_ctx->dav;
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	/* interpolate point coordinates (needed if mesh was modified) */
	//ierr = QuadratureStokesCoordinateSetUp(user->stokes_ctx->Q,dav);CHKERRQ(ierr);
	//for (e=0; e<QUAD_EDGES; e++) {
	//	ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
	//}
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);
	
	
	/* define operators A */
	ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
	/* define operators B */
	ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,1,1,1,&B);CHKERRQ(ierr);
	
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	ierr = pTatin3dDestroyContext(&user);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
//	ierr = test_pTatin3d_gmg_galerkin(argc,argv);CHKERRQ(ierr);
	ierr = test_pTatin3d_gmg_mf(argc,argv);CHKERRQ(ierr);
//	ierr = test_putatin(argc,argv);CHKERRQ(ierr);
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

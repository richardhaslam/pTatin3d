/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      ptatin_driver_asmsolve2.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

static const char help[] = "Stokes solver using Q2-Pm1 mixed finite elements.\n"
"3D prototype of the (p)ragmatic version of Tatin. (pTatin3d_v0.0)\n\n";

#include "private/daimpl.h" 

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "stokes_form_function.h"
#include "stokes_operators.h"
#include "stokes_operators_mf.h"
#include "stokes_assembly.h"
#include "dmda_element_q2p1.h"
#include "dmda_duplicate.h"
#include "dmda_project_coords.h"
#include "monitors.h"

#include "ptatin3d_energy.h"
#include "energy_assembly.h"


typedef enum { OP_TYPE_REDISC_ASM=0, OP_TYPE_REDISC_MF, OP_TYPE_GALERKIN } OperatorType;

extern PetscErrorCode SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,Mat R[],DM da[],Quadrature Q[]);

PetscErrorCode MatMultTransposeAdd_generic(Mat mat,Vec v1,Vec v2,Vec v3)
{
	Vec vt;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = VecDuplicate(v1,&vt);CHKERRQ(ierr);
	ierr = VecCopy(v1,vt);CHKERRQ(ierr);
	
	ierr = MatMult(mat,v1,vt);CHKERRQ(ierr);
	ierr = VecAXPY(vt,1.0,v2);
	ierr = VecCopy(vt,v3);CHKERRQ(ierr);
	ierr = VecDestroy(&vt);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

typedef struct {
	PetscInt     nlevels;
	OperatorType *level_type;
	Mat          *operatorA11;
	Mat          *operatorB11;
	DM           *dav_hierarchy;
	Mat          *interpolatation_v;
	Mat          *interpolatation_eta;
	Quadrature   *volQ;
	BCList       *u_bclist;
	IS           *is_stokes_field;
} AuuMultiLevelCtx;

#undef __FUNCT__  
#define __FUNCT__ "FormJacobian_Stokes"
PetscErrorCode FormJacobian_Stokes(SNES snes,Vec X,Mat *A,Mat *B,MatStructure *mstr,void *ctx)
{
  pTatinCtx         user;
  DM                stokes_pack,dau,dap;
	IS                *is;
	PhysCompStokes    stokes;
    Vec               Uloc,Ploc;
    PetscScalar       *LA_Uloc,*LA_Ploc;
	PetscBool         is_mffd = PETSC_FALSE;
	PetscBool         is_nest = PETSC_FALSE;
	PetscBool         is_shell = PETSC_FALSE;
  PetscErrorCode    ierr;
	
  PetscFunctionBegin;
	
	user = (pTatinCtx)ctx;

	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;

  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	ierr = DMCompositeGetGlobalISs(stokes_pack,&is);CHKERRQ(ierr);
	
	/* Jacobian */
	ierr = pTatin_EvaluateRheologyNonlinearities(user,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
	
	ierr = PetscTypeCompare((PetscObject)(*A),MATMFFD, &is_mffd);CHKERRQ(ierr);
	ierr = PetscTypeCompare((PetscObject)(*A),MATNEST, &is_nest);CHKERRQ(ierr);
	ierr = PetscTypeCompare((PetscObject)(*A),MATSHELL,&is_shell);CHKERRQ(ierr);

	if (is_nest) {
		Mat Auu;
		
		ierr = MatGetSubMatrix(*A,is[0],is[0],MAT_INITIAL_MATRIX,&Auu);CHKERRQ(ierr);

		is_shell = PETSC_FALSE;
		ierr = PetscTypeCompare((PetscObject)Auu,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
			ierr = MatAssemble_StokesA_AUU(Auu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		ierr = MatDestroy(&Auu);CHKERRQ(ierr);
	}
	/* If shell, do nothing */
	/* If mffd,  do nothing */
	
	ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	/* preconditioner for Jacobian */
	{
		Mat Buu,Bpp;
		
		ierr = MatGetSubMatrix(*B,is[0],is[0],MAT_INITIAL_MATRIX,&Buu);CHKERRQ(ierr);
		ierr = MatGetSubMatrix(*B,is[1],is[1],MAT_INITIAL_MATRIX,&Bpp);CHKERRQ(ierr);
		
//		ierr = Assemble_Stokes_A11_Q2(user,dau,u,dap,p,Buu);CHKERRQ(ierr);
//		ierr = Assemble_Stokes_B22_P1(user,dau,u,dap,p,Bpp);CHKERRQ(ierr);
	
		is_shell = PETSC_FALSE;
		ierr = PetscTypeCompare((PetscObject)Buu,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Buu);CHKERRQ(ierr);
			ierr = MatAssemble_StokesA_AUU(Buu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		is_shell = PETSC_FALSE;
		ierr = PetscTypeCompare((PetscObject)Bpp,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Bpp);CHKERRQ(ierr);
			ierr = MatAssemble_StokesPC_ScaledMassMatrix(Bpp,dau,dap,user->stokes_ctx->p_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		ierr = MatDestroy(&Buu);CHKERRQ(ierr);
		ierr = MatDestroy(&Bpp);CHKERRQ(ierr);		
  }
  ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	*mstr = DIFFERENT_NONZERO_PATTERN;
	
	/* clean up */
	ierr = ISDestroy(&is[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is[1]);CHKERRQ(ierr);
	ierr = PetscFree(is);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "FormJacobian_StokesMGAuu"
PetscErrorCode FormJacobian_StokesMGAuu(SNES snes,Vec X,Mat *A,Mat *B,MatStructure *mstr,void *ctx)
{
  pTatinCtx         user;
	AuuMultiLevelCtx  *mlctx;
  DM                stokes_pack,dau,dap;
	PhysCompStokes    stokes;
  Vec               Uloc,Ploc;
  PetscScalar       *LA_Uloc,*LA_Ploc;
	PetscBool         is_mffd = PETSC_FALSE;
	PetscBool         is_nest = PETSC_FALSE;
	PetscBool         is_shell = PETSC_FALSE;
	PetscContainer    container;
	PetscInt          k;
  PetscErrorCode    ierr;
	
  PetscFunctionBegin;
	
	user = (pTatinCtx)ctx;
	
	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	
	ierr = PetscObjectQuery((PetscObject)snes,"AuuMultiLevelCtx",(PetscObject*)&container);CHKERRQ(ierr);
	if (!container) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"No data with name \"AuuMultiLevelCtx\" was composed with SNES");
	ierr = PetscContainerGetPointer(container,(void*)&mlctx);CHKERRQ(ierr);
	
  
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
	/* nonlinearitiers: markers => quad points */
	ierr = pTatin_EvaluateRheologyNonlinearities(user,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);

#if 1
	/* interpolate coefficients */
	{
		int               npoints;
		DataField         PField_std;
		DataField         PField_stokes;
		MPntStd           *mp_std;
		MPntPStokes       *mp_stokes;
		
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
		
		DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
		mp_std    = PField_std->data; /* should write a function to do this */
		mp_stokes = PField_stokes->data; /* should write a function to do this */
		
		ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(npoints,mp_std,mp_stokes,mlctx->nlevels,mlctx->interpolatation_eta,mlctx->dav_hierarchy,mlctx->volQ);CHKERRQ(ierr);
	}
#endif	
	
	/* operator */
	ierr = PetscTypeCompare((PetscObject)(*A),MATMFFD, &is_mffd);CHKERRQ(ierr);
	ierr = PetscTypeCompare((PetscObject)(*A),MATNEST, &is_nest);CHKERRQ(ierr);
	ierr = PetscTypeCompare((PetscObject)(*A),MATSHELL,&is_shell);CHKERRQ(ierr);

	if (is_nest) {
		Mat Auu;
		
		ierr = MatGetSubMatrix(*A,mlctx->is_stokes_field[0],mlctx->is_stokes_field[0],MAT_INITIAL_MATRIX,&Auu);CHKERRQ(ierr);
		
		is_shell = PETSC_FALSE;
		ierr = PetscTypeCompare((PetscObject)Auu,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
			ierr = MatAssemble_StokesA_AUU(Auu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		ierr = MatDestroy(&Auu);CHKERRQ(ierr);
	}
	/* If shell, do nothing */
	/* If mffd,  do nothing */
	
	ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	/* preconditioner operator for Jacobian */
	{
		Mat Buu,Bpp;
		
		ierr = MatGetSubMatrix(*B,mlctx->is_stokes_field[0],mlctx->is_stokes_field[0],MAT_INITIAL_MATRIX,&Buu);CHKERRQ(ierr);
		ierr = MatGetSubMatrix(*B,mlctx->is_stokes_field[1],mlctx->is_stokes_field[1],MAT_INITIAL_MATRIX,&Bpp);CHKERRQ(ierr);
		
		is_shell = PETSC_FALSE;
		ierr = PetscTypeCompare((PetscObject)Buu,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Buu);CHKERRQ(ierr);
			ierr = MatAssemble_StokesA_AUU(Buu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		is_shell = PETSC_FALSE;
		ierr = PetscTypeCompare((PetscObject)Bpp,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Bpp);CHKERRQ(ierr);
			ierr = MatAssemble_StokesPC_ScaledMassMatrix(Bpp,dau,dap,user->stokes_ctx->p_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		ierr = MatDestroy(&Buu);CHKERRQ(ierr);
		ierr = MatDestroy(&Bpp);CHKERRQ(ierr);		
  }
  ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	*mstr = DIFFERENT_NONZERO_PATTERN;
#if 1	
	/* Buu preconditioner for all other levels in the hierarchy */
	for (k=mlctx->nlevels-2; k>=0; k--) {
		
		switch (mlctx->level_type[k]) {
				
			case OP_TYPE_REDISC_ASM:
			{
				ierr = MatZeroEntries(mlctx->operatorB11[k]);CHKERRQ(ierr);
				ierr = MatAssemble_StokesA_AUU(mlctx->operatorB11[k],mlctx->dav_hierarchy[k],mlctx->u_bclist[k],mlctx->volQ[k]);CHKERRQ(ierr);
			}
				break;
				
			case OP_TYPE_REDISC_MF:
			{
			}
				break;
				
			case OP_TYPE_GALERKIN:
			{
				Mat Auu_k;
				
				if (k==mlctx->nlevels-1) {
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
				}	
				PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Galerkin :: assembled operator \n", k);
				
				/* should move coarse grid assembly into jacobian */
				ierr = MatPtAP(mlctx->operatorA11[k+1],mlctx->interpolatation_v[k+1],MAT_INITIAL_MATRIX,1.0,&Auu_k);CHKERRQ(ierr);
				
				mlctx->operatorA11[k] = Auu_k;
				mlctx->operatorB11[k] = Auu_k;
			}
				break;
				
			default:
				break;
		}
	}	
#endif	
	
	
	/* clean up */
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_gmg2_material_points"
PetscErrorCode pTatin3d_gmg2_material_points(int argc,char **argv)
{
	DM        multipys_pack,dav,dap;
	pTatinCtx user;
	Mat       A,B,operatorA11[10],operatorB11[10];
	Vec       X,F;
	IS        *is_stokes_field;
	SNES      snes;
	KSP       ksp;
	PC        pc;
	PetscBool flg;
	OperatorType  level_type[10];
	
	PetscMPIInt    rank;
	DM             dav_hierarchy[10];
	Mat            interpolatation_v[10],interpolatation_eta[10];
	PetscInt       nlevels,k,max;
	Quadrature     volQ[10];
	BCList         u_bclist[10];
	PetscInt       kk;
	AuuMultiLevelCtx mlctx;
	PetscInt newton_its,picard_its;
	PhysCompEnergy  energy;
	PetscBool active_energy;
	Vec T;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dParseOptions(user);CHKERRQ(ierr);
	
	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	/* Check if model is being restarted from a checkpointed file */
	ierr = pTatin3dRestart(user);CHKERRQ(ierr);
	
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
	dap           = user->stokes_ctx->dap;
	
	/* IF I DON'T DO THIS, THE IS's OBTAINED FROM DMCompositeGetGlobalISs() are wrong !! */
	{
		Vec X;

		ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = DMRestoreGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	}
	 ierr = DMCompositeGetGlobalISs(multipys_pack,&is_stokes_field);CHKERRQ(ierr);	
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	/* generate energy solver */
	/* NOTE - Generating the thermal solver here will ensure that the initial geometry on the mechanical model is copied */
	/* NOTE - Calling pTatinPhysCompActivate_Energy() after pTatin3dCreateMaterialPoints() is essential */
	{
		PetscBool load_energy = PETSC_FALSE;
		
		PetscOptionsGetBool(PETSC_NULL,"-activate_energy",&load_energy,PETSC_NULL);
		ierr = pTatinPhysCompActivate_Energy(user,load_energy);CHKERRQ(ierr);
		ierr = pTatinContextValid_Energy(user,&active_energy);CHKERRQ(ierr);
	}
	if (active_energy) {
		ierr = pTatinGetContext_Energy(user,&energy);CHKERRQ(ierr);
		
		ierr = DMCreateGlobalVector(energy->daT,&T);CHKERRQ(ierr);
		ierr = pTatinPhysCompAttachData_Energy(user,T,PETSC_NULL);CHKERRQ(ierr);
	}
	
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);
	
	/* set up mg */
	user->stokes_ctx->dav->ops->coarsenhierarchy = DMCoarsenHierarchy2_DA;
	
	nlevels = 1;
	PetscOptionsGetInt(PETSC_NULL,"-dau_nlevels",&nlevels,0);
	PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%d x %d x %d) : MG levels %d  \n", user->mx,user->my,user->mz,nlevels );
	dav_hierarchy[ nlevels-1 ] = dav;
	ierr = PetscObjectReference((PetscObject)dav);CHKERRQ(ierr);
	
	/* Coarsen nlevels - 1 times, and add levels into list so that level 0 is the coarsest */
	{
		DM *coarsened_list;
		
		ierr = PetscMalloc(sizeof(DM)*(nlevels-1),&coarsened_list);CHKERRQ(ierr);
		ierr = DMCoarsenHierarchy(dav,nlevels-1,coarsened_list);CHKERRQ(ierr);
		for (k=0; k<nlevels-1; k++) {
			dav_hierarchy[ nlevels-2-k ] = coarsened_list[k];
		}
		
		PetscFree(coarsened_list);
	}
	/* Set all dav's to be of type Q2 */
	for (k=0; k<nlevels-1; k++) {
		ierr = PetscObjectSetOptionsPrefix((PetscObject)dav_hierarchy[k],"stk_velocity_");CHKERRQ(ierr);
		ierr = DMDASetElementType_Q2(dav_hierarchy[k]);CHKERRQ(ierr);
	}
	
	/* Report mesh sizes */
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	for (k=0; k<nlevels; k++) {
		PetscInt lmx,lmy,lmz;
		ierr = DMDAGetSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"         level [%2D]: global Q2 elements (%D x %D x %D) \n", k,lmx,lmy,lmz );
	}		
	for (k=0; k<nlevels; k++) {
		PetscInt nels,nen;
		const PetscInt *els;
		PetscInt lmx,lmy,lmz,si,sj,sk;
		
		ierr = DMDAGetElements_pTatinQ2P1(dav_hierarchy[k],&nels,&nen,&els);CHKERRQ(ierr);
		ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		if (rank<10) {
			PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: local Q2 elements  (%D x %D x %D) \n", rank, k,lmx,lmy,lmz );
		}
		
		ierr = DMDAGetCornersElementQ2(dav_hierarchy[k],&si,&sj,&sk,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		si = si/2;
		sj = sj/2;
		sk = sk/2;
		if (rank<10) {
			PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: element range [%D - %D] x [%D - %D] x [%D - %D] \n", rank, k,si,si+lmx-1,sj,sj+lmy-1,sk,sk+lmz-1 );
		}
	}
	
	mlctx.is_stokes_field = is_stokes_field;
	mlctx.nlevels = nlevels;
	mlctx.dav_hierarchy = dav_hierarchy;
	
	/* inject coordinates */
	ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);
	
	/* Define interpolation operators for velocity space */
	interpolatation_v[0] = PETSC_NULL;
	for (k=0; k<nlevels-1; k++) {
		ierr = DMGetInterpolation(dav_hierarchy[k],dav_hierarchy[k+1],&interpolatation_v[k+1],PETSC_NULL);CHKERRQ(ierr);
	}
	
	mlctx.interpolatation_v = interpolatation_v;

	/* Define interpolation operators for scalar space */
	interpolatation_eta[0] = PETSC_NULL;
	for (k=1; k<nlevels; k++) {
		ierr = MatMAIJRedimension(interpolatation_v[k],1,&interpolatation_eta[k]);CHKERRQ(ierr);
	}
	
	mlctx.interpolatation_eta = interpolatation_eta;
	
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

	mlctx.volQ = volQ;
	
	/* update markers = >> gauss points */
	{
		int               npoints;
		DataField         PField_std;
		DataField         PField_stokes;
		MPntStd           *mp_std;
		MPntPStokes       *mp_stokes;
		
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
		
		DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
		mp_std    = PField_std->data; /* should write a function to do this */
		mp_stokes = PField_stokes->data; /* should write a function to do this */
		
		//ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(npoints,mp_std,mp_stokes,nlevels,interpolatation_eta,dav_hierarchy,volQ);CHKERRQ(ierr);
	}
	
	/* define boundary conditions */
	for (k=0; k<nlevels-1; k++) {
		ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);
	}
	u_bclist[nlevels-1] = user->stokes_ctx->u_bclist;

	ierr = pTatinModel_ApplyBoundaryConditionMG(nlevels,u_bclist,dav_hierarchy,user->model,user);CHKERRQ(ierr);

	mlctx.u_bclist = u_bclist;
	
	
	/* ============================================== */
	/* CONFIGURE OPERATORS */

	/* A operator */
	ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
	/* memory saving - only need daU IF you want to split A11 into A11uu,A11vv,A11ww */
	{
		MatStokesMF mf;
		
		ierr = MatShellGetMatStokesMF(A,&mf);CHKERRQ(ierr);
		ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
		mf->daU = PETSC_NULL;
	}

	/* B operator */
//	ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,PETSC_FALSE,PETSC_TRUE,PETSC_TRUE,&B);CHKERRQ(ierr);
	{
		Mat         Aup,Apu,Spp,bA[2][2];
		MatStokesMF StkCtx;
		PetscInt    i,j;
	
		ierr = MatShellGetMatStokesMF(A,&StkCtx);CHKERRQ(ierr);

		/* Schur complement */
		ierr = DMGetMatrix(dap,MATSBAIJ,&Spp);CHKERRQ(ierr);
		ierr = MatSetOptionsPrefix(Spp,"S*_");CHKERRQ(ierr);
		ierr = MatSetOption(Spp,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
		ierr = MatSetFromOptions(Spp);CHKERRQ(ierr);
	
		/* A12 */
		ierr = StokesQ2P1CreateMatrix_MFOperator_A12(StkCtx,&Aup);CHKERRQ(ierr);
		ierr = MatSetOptionsPrefix(Aup,"Bup_");CHKERRQ(ierr);
		ierr = MatSetFromOptions(Aup);CHKERRQ(ierr);

		/* A21 */
		ierr = StokesQ2P1CreateMatrix_MFOperator_A21(StkCtx,&Apu);CHKERRQ(ierr);
		ierr = MatSetOptionsPrefix(Apu,"Bpu_");CHKERRQ(ierr);
		ierr = MatSetFromOptions(Apu);CHKERRQ(ierr);
	
		/* nest */
		bA[0][0] = PETSC_NULL; bA[0][1] = Aup;
		bA[1][0] = Apu;        bA[1][1] = Spp;
	
		/* Set fine A11 into nest */
		//bA[0][0] = operatorA11[nlevels-1];		
		
		ierr = MatCreateNest(PETSC_COMM_WORLD,2,is_stokes_field,2,is_stokes_field,&bA[0][0],&B);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		//ierr = MatNestSetSubMat(B,0,0,operatorA11[nlevels-1]);CHKERRQ(ierr);
		
		/* tidy up - hand back destruction to B */
		ierr = MatDestroy(&Aup);CHKERRQ(ierr);
		ierr = MatDestroy(&Apu);CHKERRQ(ierr);
		ierr = MatDestroy(&Spp);CHKERRQ(ierr);
	}
	//ierr = MatShellSetOperation(B,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))MatMultTransposeAdd_generic);CHKERRQ(ierr);
	
	
	/* A11 operator */
	
	/* defaults */
	level_type[0] = OP_TYPE_REDISC_ASM;
	for (k=1; k<nlevels; k++) {
		level_type[k] = OP_TYPE_REDISC_MF;
	}
	
	max = nlevels;
	ierr = PetscOptionsGetIntArray(PETSC_NULL,"-A11_operator_type",(PetscInt*)level_type,&max,&flg);CHKERRQ(ierr);
	for (k=nlevels-1; k>=0; k--) {

		switch (level_type[k]) {

			case 0:
			{
				Mat Auu;
				PetscBool same1 = PETSC_FALSE,same2 = PETSC_FALSE,same3 = PETSC_FALSE;
				
				/* use -stk_velocity_da_mat_type sbaij or -Buu_da_mat_type sbaij */
				PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Re-discretisation :: assembled operator \n", k);
				ierr = DMGetMatrix(dav_hierarchy[k],MATSBAIJ,&Auu);CHKERRQ(ierr);
				ierr = MatSetOptionsPrefix(Auu,"Buu_");CHKERRQ(ierr);
				ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);
				ierr = PetscTypeCompare((PetscObject)Auu,MATSBAIJ,&same1);CHKERRQ(ierr);
				ierr = PetscTypeCompare((PetscObject)Auu,MATSEQSBAIJ,&same2);CHKERRQ(ierr);
				ierr = PetscTypeCompare((PetscObject)Auu,MATMPISBAIJ,&same3);CHKERRQ(ierr);
				if (same1||same2||same3) {
					ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
				}
				/* should move assembly into jacobian */
				ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
				//ierr = MatAssemble_StokesA_AUU(Auu,dav_hierarchy[k],u_bclist[k],volQ[k]);CHKERRQ(ierr);
				
				operatorA11[k] = Auu;
				operatorB11[k] = Auu;
				ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
				
			}
				break;

			case 1:
			{
				Mat Auu;
				MatA11MF mf,A11Ctx;

				PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Re-discretisation :: matrix free operator \n", k);
				ierr = MatA11MFCreate(&A11Ctx);CHKERRQ(ierr);
				ierr = MatA11MFSetup(A11Ctx,dav_hierarchy[k],volQ[k],u_bclist[k]);CHKERRQ(ierr);
				
				ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
				ierr = MatShellGetMatA11MF(Auu,&mf);CHKERRQ(ierr);
				ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
				mf->daU = PETSC_NULL;				
				operatorA11[k] = Auu;
				ierr = MatShellSetOperation(Auu,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))MatMultTransposeAdd_generic);CHKERRQ(ierr);
				
				{
					PetscBool use_low_order_geometry = PETSC_FALSE;
					
					ierr = PetscOptionsGetBool(PETSC_NULL,"-use_low_order_geometry",&use_low_order_geometry,PETSC_NULL);CHKERRQ(ierr);
					if (use_low_order_geometry==PETSC_TRUE) {
						Mat Buu;
						
						PetscPrintf(PETSC_COMM_WORLD,"Activiting low order A11 operator \n");
						ierr = StokesQ2P1CreateMatrix_MFOperator_A11LowOrder(A11Ctx,&Buu);CHKERRQ(ierr);
						ierr = MatShellGetMatA11MF(Buu,&mf);CHKERRQ(ierr);
						ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
						mf->daU = PETSC_NULL;				
						operatorB11[k] = Buu;
						ierr = MatShellSetOperation(Buu,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))MatMultTransposeAdd_generic);CHKERRQ(ierr);
					} else {
						operatorB11[k] = Auu;
						ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
					}
				}
				
				
				ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
			}
				break;
				
			case 2:
			{
				Mat Auu;

				if (k==nlevels-1) {
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
				}	
				
				PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Galerkin :: assembled operator \n", k);
				
				/* should move coarse grid assembly into jacobian */
				ierr = MatPtAP(operatorA11[k+1],interpolatation_v[k+1],MAT_INITIAL_MATRIX,1.0,&Auu);CHKERRQ(ierr);
				
				operatorA11[k] = Auu;
				operatorB11[k] = Auu;
				ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
			}
				break;
				
			default:
				break;
		}
	}	
	
	mlctx.level_type  = level_type;
	mlctx.operatorA11 = operatorA11;
	mlctx.operatorB11 = operatorB11;
	
	
	
	/* Set fine A11 into nest */
	ierr = MatNestSetSubMat(B,0,0,operatorA11[nlevels-1]);CHKERRQ(ierr);
	
	/* ============================================== */
	
	
	

	/* work vector for solution and residual */
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
    /* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	if (active_energy) {
		ierr = pTatinPhysCompEnergy_MPProjectionQ1(user);CHKERRQ(ierr);
	}
	
	/* initial condition */
	ierr = pTatinModel_ApplyInitialSolution(user->model,user,X);CHKERRQ(ierr);
    
    /* initial viscosity  */
	ierr = pTatinModel_ApplyInitialStokesVariableMarkers(user->model,user,X);CHKERRQ(ierr);
    
	/* boundary condition */
	{
		Vec velocity,pressure;
		
		ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		ierr = BCListInsert(user->stokes_ctx->u_bclist,velocity);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);

		if (active_energy) {
			ierr = BCListInsert(energy->T_bclist,T);CHKERRQ(ierr);
		}
	}
	
	
	
	
	
	
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	
	ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
	
	// activate mffd via -snes_mf_operator
	ierr = SNESSetJacobian(snes,A,B,FormJacobian_StokesMGAuu,user);CHKERRQ(ierr);
	// Force mffd
	//ierr = SNESSetJacobian(snes,B,B,FormJacobian_StokesMGAuu,user);CHKERRQ(ierr);
	
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
	

	{
		PetscContainer container;
		
		ierr = PetscContainerCreate(PETSC_COMM_WORLD,&container);CHKERRQ(ierr);
		ierr = PetscContainerSetPointer(container,(void*)&mlctx);CHKERRQ(ierr);
		ierr = PetscObjectCompose((PetscObject)snes,"AuuMultiLevelCtx",(PetscObject)container);CHKERRQ(ierr);

	}
	
	
	
	/* configure for fieldsplit */
	ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);

	ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_StdoutStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
//	ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_ParaviewStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
	ierr = SNESMonitorSet(snes,pTatin_SNESMonitor_ParaviewStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	
	ierr = PCFieldSplitSetIS(pc,"u",is_stokes_field[0]);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc,"p",is_stokes_field[1]);CHKERRQ(ierr);

	/* configure uu split for galerkin multi-grid */
	{
		PetscInt nsplits;
		PC       pc_i;
		KSP      *sub_ksp,ksp_coarse,ksp_smoother;
		
		ierr = KSPSetUp(ksp);CHKERRQ(ierr);
		ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
		
		ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
		ierr = PCSetType(pc_i,PCMG);CHKERRQ(ierr);
		ierr = PCMGSetLevels(pc_i,nlevels,PETSC_NULL);CHKERRQ(ierr);
		ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
		ierr = PCMGSetGalerkin(pc_i,PETSC_FALSE);CHKERRQ(ierr);
		
		for( k=1; k<nlevels; k++ ){
			ierr = PCMGSetInterpolation(pc_i,k,interpolatation_v[k]);CHKERRQ(ierr);
		}
		
		/* drop the operators in - i presume this will also need to be performed inside the jacobian each time the operators are modified */
		/* No - it looks like PCSetUp_MG will call set operators on all levels if the SetOperators was called on the finest, which should/is done by the SNES */
		ierr = PCMGGetCoarseSolve(pc_i,&ksp_coarse);CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp_coarse,operatorA11[0],operatorA11[0],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		for( k=1; k<nlevels; k++ ){
			PetscBool use_low_order_geometry = PETSC_FALSE;

			
			ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);

			// use A for smoother, B for residual
			ierr = PetscOptionsGetBool(PETSC_NULL,"-use_low_order_geometry",&use_low_order_geometry,PETSC_NULL);CHKERRQ(ierr);
			if (use_low_order_geometry==PETSC_TRUE) {
				ierr = KSPSetOperators(ksp_smoother,operatorB11[k],operatorB11[k],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				
				//ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorB11[k],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
			} else {
				// Use A for smoother, lo
				ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
			}
		
		}
	}

	ierr = pTatinModel_Output(user->model,user,X,"icbc");CHKERRQ(ierr);
	
	
#if 0
          
	SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);        
	ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
	ierr = pTatinModel_Output(user->model,user,X,"continuation");CHKERRQ(ierr);
        
#endif
	
#if 1
{
	PetscInt snes_its;

	SNESGetTolerances(snes,0,0,0,&snes_its,0);
    
	/* switch to linear rheology to use the viscosity set on marker by the initialStokeVariables */

	user->rheology_constants.rheology_type = RHEOLOGY_VISCOUS;

	/* do a linear solve */
	SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);
	PetscPrintf(PETSC_COMM_WORLD,"############## LINEAR STAGE ##############\n");
	ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
	ierr = pTatinModel_Output(user->model,user,X,"linear_stage");CHKERRQ(ierr);
	
	/* switch to non-linear rheology */
	user->rheology_constants.rheology_type = RHEOLOGY_VP_STD;

	/* do a picard solve */
	//ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
	
	//SNESGetTolerances(snes,0,0,0,&snes_its,0);
	picard_its = snes_its;
	PetscOptionsGetInt(PETSC_NULL,"-picard_its",&picard_its,0);
	SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,picard_its,PETSC_DEFAULT);
	PetscPrintf(PETSC_COMM_WORLD,"############## PICARD STAGE ##############\n");
	ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
	ierr = pTatinModel_Output(user->model,user,X,"picard_stage");CHKERRQ(ierr);

	
	/*
	newton_its = 0;
	PetscOptionsGetInt(PETSC_NULL,"-newton_its",&newton_its,0);
	if (newton_its>0) {
		// Force mffd
		ierr = SNESSetJacobian(snes,B,B,FormJacobian_StokesMGAuu,user);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"############## NEWTON STAGE ##############\n");
		ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
		ierr = pTatinModel_Output(user->model,user,X,"newton_stage");CHKERRQ(ierr);
	}
	*/
}
#endif

	newton_its = 0;
	PetscOptionsGetInt(PETSC_NULL,"-newton_its",&newton_its,0);
	if (newton_its>0) {
		SNES snes_newton;
		PetscContainer container;
		
		ierr = SNESCreate(PETSC_COMM_WORLD,&snes_newton);CHKERRQ(ierr);
		ierr = SNESSetOptionsPrefix(snes_newton,"n_");CHKERRQ(ierr);
		ierr = SNESSetFunction(snes_newton,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
		// Force mffd
		ierr = SNESSetJacobian(snes_newton,B,B,FormJacobian_StokesMGAuu,user);CHKERRQ(ierr);
		
		ierr = SNESSetFromOptions(snes_newton);CHKERRQ(ierr);
		
		ierr = PetscContainerCreate(PETSC_COMM_WORLD,&container);CHKERRQ(ierr);
		ierr = PetscContainerSetPointer(container,(void*)&mlctx);CHKERRQ(ierr);
		ierr = PetscObjectCompose((PetscObject)snes_newton,"AuuMultiLevelCtx",(PetscObject)container);CHKERRQ(ierr);
		
		/* configure for fieldsplit */
		ierr = SNESGetKSP(snes_newton,&ksp);CHKERRQ(ierr);
		ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_StdoutStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
		ierr = SNESMonitorSet(snes_newton,pTatin_SNESMonitor_ParaviewStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
		
		ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
		ierr = PCFieldSplitSetIS(pc,"u",is_stokes_field[0]);CHKERRQ(ierr);
		ierr = PCFieldSplitSetIS(pc,"p",is_stokes_field[1]);CHKERRQ(ierr);
		
		/* configure uu split for galerkin multi-grid */
		{
			PetscInt nsplits;
			PC       pc_i;
			KSP      *sub_ksp,ksp_coarse,ksp_smoother;
			
			ierr = KSPSetUp(ksp);CHKERRQ(ierr);
			ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
			
			ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
			ierr = PCSetType(pc_i,PCMG);CHKERRQ(ierr);
			ierr = PCMGSetLevels(pc_i,nlevels,PETSC_NULL);CHKERRQ(ierr);
			ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
			ierr = PCMGSetGalerkin(pc_i,PETSC_FALSE);CHKERRQ(ierr);
			
			for( k=1; k<nlevels; k++ ){
				ierr = PCMGSetInterpolation(pc_i,k,interpolatation_v[k]);CHKERRQ(ierr);
			}
			
			/* drop the operators in - i presume this will also need to be performed inside the jacobian each time the operators are modified */
			/* No - it looks like PCSetUp_MG will call set operators on all levels if the SetOperators was called on the finest, which should/is done by the SNES */
			ierr = PCMGGetCoarseSolve(pc_i,&ksp_coarse);CHKERRQ(ierr);
			ierr = KSPSetOperators(ksp_coarse,operatorA11[0],operatorA11[0],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
			for( k=1; k<nlevels; k++ ){
				PetscBool use_low_order_geometry = PETSC_FALSE;
				
				ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);
				// use A for smoother, B for residual
				ierr = PetscOptionsGetBool(PETSC_NULL,"-use_low_order_geometry",&use_low_order_geometry,PETSC_NULL);CHKERRQ(ierr);
				if (use_low_order_geometry==PETSC_TRUE) {
					ierr = KSPSetOperators(ksp_smoother,operatorB11[k],operatorB11[k],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				} else {
					// Use A for smoother, low
					ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				}
			}
		}
		
		SNESSetTolerances(snes_newton,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,newton_its,PETSC_DEFAULT);
		PetscPrintf(PETSC_COMM_WORLD,"############## NEWTON STAGE ##############\n");
		ierr = SNESSolve(snes_newton,PETSC_NULL,X);CHKERRQ(ierr);
		ierr = pTatinModel_Output(user->model,user,X,"newton_stage");CHKERRQ(ierr);

		ierr = SNESDestroy(&snes_newton);CHKERRQ(ierr);
	}
	
	
	
	PetscPrintf(PETSC_COMM_WORLD,"[Initial condition] Timestep[%d]: time %lf Myr \n", user->step, user->time );
	for (kk=0; kk<user->nsteps; kk++) {
		PetscInt tk = user->step+1;
		
		/* do solve */

		/* update markers */
		
		
		
		user->time += 0.12;
		user->step++;
		PetscPrintf(PETSC_COMM_WORLD,"Timestep[%d] : Cycle[%d/%d] : time %lf Myr \n", tk, kk, user->nsteps-1, user->time );

		
		if ((kk+1)%5==0) {
			char name[100];
			
			sprintf(name,"step%.6d",tk);
			ierr = pTatinModel_Output(user->model,user,X,name);CHKERRQ(ierr);
		}
		
		if ((kk+1)%10==0) {
			char name[100];
			
			PetscPrintf(PETSC_COMM_WORLD,"  checkpointing ptatin :: Model timestep %d : time %lf Myr : cycle[%d/%d] \n", tk,user->time,kk, user->nsteps-1 );
			/* check point test */
			//	ierr = pTatin3dContextSave(user,"checkpoint.file");CHKERRQ(ierr);
			//	ierr = pTatin3dContextLoad(user,"checkpoint.file");CHKERRQ(ierr);
			
			sprintf(name,"step%.6d",tk);
			ierr = pTatin3dCheckpoint(user,X,name);CHKERRQ(ierr);
		}

	}
	
	/* test viewer */
	//ierr = pTatinModel_Output(user->model,user,X,"step001");CHKERRQ(ierr);
	
	
	/* Clean up */
	for (k=0; k<nlevels-1; k++) {
		ierr = BCListDestroy(&u_bclist[k]);CHKERRQ(ierr);
	}
	for (k=0; k<nlevels-1; k++) {
		ierr = QuadratureDestroy(&volQ[k]);CHKERRQ(ierr);
	}
	for (k=0; k<nlevels; k++) {
		ierr = MatDestroy(&operatorA11[k]);CHKERRQ(ierr);
		ierr = MatDestroy(&operatorB11[k]);CHKERRQ(ierr);
	}
	
	for (k=1; k<nlevels; k++) {
		ierr = MatDestroy(&interpolatation_v[k]);CHKERRQ(ierr);
		ierr = MatDestroy(&interpolatation_eta[k]);CHKERRQ(ierr);
	}
	
	for (k=0; k<nlevels; k++) {
		ierr = DMDestroy(&dav_hierarchy[k]);CHKERRQ(ierr);
	}
	
	ierr = ISDestroy(&is_stokes_field[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is_stokes_field[1]);CHKERRQ(ierr);
	ierr = PetscFree(is_stokes_field);CHKERRQ(ierr);
	
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);

	{	
		PetscContainer    container;
		ierr = PetscObjectQuery((PetscObject)snes,"AuuMultiLevelCtx",(PetscObject*)&container);CHKERRQ(ierr);
		if (!container) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"No data with name \"AuuMultiLevelCtx\" was composed with SNES");
		ierr = PetscContainerDestroy(&container);CHKERRQ(ierr);
	}
	ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = pTatin3dDestroyContext(&user);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	
	ierr = PetscInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	ierr = pTatin3d_gmg2_material_points(argc,argv);CHKERRQ(ierr);

	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}

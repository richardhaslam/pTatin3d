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
 **    Filename:      ptatin_driver_nostokessolve.c
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


#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "stokes_form_function.h"
#include "stokes_operators.h"

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_material_points_gmg"
PetscErrorCode pTatin3d_material_points_gmg(int argc,char **argv)
{
	DM              multipys_pack,dav,dap;
	PetscErrorCode ierr;
	pTatinCtx user;

	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dParseOptions(user);CHKERRQ(ierr);

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
	dap           = user->stokes_ctx->dap;
	
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

/*	
	{	
		PetscScalar zero = 0.0;
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}
	{
		BCList flat;
		
		ierr = BCListFlattenedCreate(user->stokes_ctx->u_bclist,&flat);CHKERRQ(ierr);
		ierr = BCListDestroy(&user->stokes_ctx->u_bclist);CHKERRQ(ierr);
		user->stokes_ctx->u_bclist = flat;
	}
*/	
	
	/* test form function */
	{
		Vec X,F;
		SNES snes;
		Mat JMF;
		
		ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
		
		ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
		ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
		ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
		ierr = SNESDestroy(&snes);CHKERRQ(ierr);
		ierr = VecDestroy(&X);CHKERRQ(ierr);
		ierr = VecDestroy(&F);CHKERRQ(ierr);
	}
	
	{
		Vec X,Y;

		ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = VecDuplicate(X,&Y);CHKERRQ(ierr);

		ierr = MF_Stokes(X,Y,(void*)user);CHKERRQ(ierr);

		ierr = VecDestroy(&X);CHKERRQ(ierr);
		ierr = VecDestroy(&Y);CHKERRQ(ierr);
	}
	
#if 0
	{
		Vec X,F;
		SNES snes;
		Mat JMF;
		
		ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

		ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
		ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
		ierr = MatCreateSNESMF(snes,&JMF);CHKERRQ(ierr);
		ierr = SNESSetJacobian(snes,JMF,JMF,0,0);CHKERRQ(ierr);
		
//		ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
		ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

		ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
		
		ierr = SNESDestroy(&snes);CHKERRQ(ierr);
		ierr = VecDestroy(&X);CHKERRQ(ierr);
		ierr = VecDestroy(&F);CHKERRQ(ierr);
	}
#endif	
	
	
	{
		Mat B,subA;
		Vec X,F;
		IS *is;
		Vec u,p;
		PetscInt n,mu,mp,Mu,Mp,start,offset;
		IS isU,isV,isW;
		Mat Aii;
		
		ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
		
		ierr = DMCompositeGetGlobalISs(user->stokes_ctx->stokes_pack,&is);CHKERRQ(ierr);

		
		/* Sizes */
		ierr = DMCompositeGetAccess(multipys_pack,X,&u,&p);CHKERRQ(ierr);
		ierr = VecGetSize(u,&Mu);CHKERRQ(ierr);
		ierr = VecGetLocalSize(u,&mu);CHKERRQ(ierr);
		ierr = VecGetSize(p,&Mp);CHKERRQ(ierr);
		ierr = VecGetLocalSize(p,&mp);CHKERRQ(ierr);
		ierr = VecGetOwnershipRange(u,&start,PETSC_NULL);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&u,&p);CHKERRQ(ierr);
		
		n = (mu/3);
		offset = start + 0;
		ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&isU);CHKERRQ(ierr);
		offset = start + 1;
		ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&isV);CHKERRQ(ierr);
		offset = start + 2;
		ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&isW);CHKERRQ(ierr);
			
		
		/* test operator, A */
		PetscPrintf(PETSC_COMM_WORLD,"operatorA\n");
		ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&B);CHKERRQ(ierr);
		ierr = MatMult(B,X,F);CHKERRQ(ierr);
		
		// A11
		ierr = MatGetSubMatrix(B,is[0],is[0],MAT_INITIAL_MATRIX,&subA);CHKERRQ(ierr);
		ierr = MatDestroy(&subA);CHKERRQ(ierr);

		// A11vv
		ierr = MatGetSubMatrix(B,isV,isV,MAT_INITIAL_MATRIX,&Aii);CHKERRQ(ierr);
		ierr = MatDestroy(&Aii);CHKERRQ(ierr);

		
		// A12
		ierr = MatGetSubMatrix(B,is[0],is[1],MAT_INITIAL_MATRIX,&subA);CHKERRQ(ierr);
		ierr = MatDestroy(&subA);CHKERRQ(ierr);

		// A21
		ierr = MatGetSubMatrix(B,is[1],is[0],MAT_INITIAL_MATRIX,&subA);CHKERRQ(ierr);
		ierr = MatDestroy(&subA);CHKERRQ(ierr);

		// A22
		ierr = MatGetSubMatrix(B,is[1],is[1],MAT_INITIAL_MATRIX,&subA);CHKERRQ(ierr);
		ierr = MatDestroy(&subA);CHKERRQ(ierr);

		ierr = MatDestroy(&B);CHKERRQ(ierr);

		/* ========================================= */
			
		/* test pc operator, B */
		PetscPrintf(PETSC_COMM_WORLD,"operatorB\n");
		ierr = StokesQ2P1CreateMatrixNest_Operator(user->stokes_ctx,1,1,1,&B);CHKERRQ(ierr);
		ierr = MatMult(B,X,F);CHKERRQ(ierr);

		/* Get A11 */
		ierr = MatGetSubMatrix(B,is[0],is[0],MAT_INITIAL_MATRIX,&subA);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
		ierr = MatView(subA,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		
		/* Get Auu,Avv,Aww */	
		ierr = MatGetSubMatrix(subA,isU,isU,MAT_INITIAL_MATRIX,&Aii);CHKERRQ(ierr);
		ierr = MatDestroy(&Aii);CHKERRQ(ierr);

		ierr = MatGetSubMatrix(subA,isV,isV,MAT_INITIAL_MATRIX,&Aii);CHKERRQ(ierr);
		ierr = MatDestroy(&Aii);CHKERRQ(ierr);

		ierr = MatGetSubMatrix(subA,isW,isW,MAT_INITIAL_MATRIX,&Aii);CHKERRQ(ierr);
		ierr = MatDestroy(&Aii);CHKERRQ(ierr);
	
		ierr = MatDestroy(&subA);CHKERRQ(ierr);

		/* Get A12 */
		ierr = MatGetSubMatrix(B,is[0],is[1],MAT_INITIAL_MATRIX,&subA);CHKERRQ(ierr);
		ierr = MatDestroy(&subA);CHKERRQ(ierr);
		/* A21 */
		ierr = MatGetSubMatrix(B,is[1],is[0],MAT_INITIAL_MATRIX,&subA);CHKERRQ(ierr);
		ierr = MatDestroy(&subA);CHKERRQ(ierr);
		/* A22 */
		ierr = MatGetSubMatrix(B,is[1],is[1],MAT_INITIAL_MATRIX,&subA);CHKERRQ(ierr);
		ierr = MatDestroy(&subA);CHKERRQ(ierr);
		
		ierr = MatDestroy(&B);CHKERRQ(ierr);

		
		ierr = ISDestroy(&isU);CHKERRQ(ierr);
		ierr = ISDestroy(&isV);CHKERRQ(ierr);
		ierr = ISDestroy(&isW);CHKERRQ(ierr);

		ierr = VecDestroy(&X);CHKERRQ(ierr);
		ierr = VecDestroy(&F);CHKERRQ(ierr);

		ierr = ISDestroy(&is[0]);CHKERRQ(ierr);
		ierr = ISDestroy(&is[1]);CHKERRQ(ierr);
		ierr = PetscFree(is);CHKERRQ(ierr);
	}
	
	
	/* test viewer */
	DataBucketView(((PetscObject)multipys_pack)->comm, user->materialpoint_db,"materialpoint_stokes",DATABUCKET_VIEW_STDOUT);
	
	
	{
		Vec X;
		
		ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = pTatinModel_Output(user->model,user,X,"icbc");CHKERRQ(ierr);
		ierr = DMRestoreGlobalVector(multipys_pack,&X);CHKERRQ(ierr);

	}
	
	/* dummy outputting */
	{
		Vec X;
		
		ierr = pTatin3d_ModelOutput_MPntStd(user,"step0");CHKERRQ(ierr);

		ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(user,X,"step0");CHKERRQ(ierr);
		ierr = DMRestoreGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	}
	
	/* test generic viewer */
	{
//		const int nf=2;
//		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes }; 
		const int nf=1;
		const MaterialPointField mp_prop_list[] = { MPField_Stokes }; 
		ierr = SwarmViewGeneric_ParaView(user->materialpoint_db,nf,mp_prop_list,user->outputpath,"test0");CHKERRQ(ierr);
	}
	
	
	ierr = pTatin3dDestroyContext(&user);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	ierr = pTatin3d_material_points_gmg(argc,argv);CHKERRQ(ierr);
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

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
 **    filename:   ptatin_driver_basic.c
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

#include "petsc/private/dmdaimpl.h" 

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "ptatin_log.h"
#include "stokes_form_function.h"
#include "stokes_operators.h"
#include "stokes_operators_mf.h"
#include "stokes_assembly.h"
#include "dmda_element_q2p1.h"
#include "dmda_duplicate.h"
#include "dmda_project_coords.h"
#include "monitors.h"


#undef __FUNCT__  
#define __FUNCT__ "FormJacobian_Stokes"
PetscErrorCode FormJacobian_Stokes(SNES snes,Vec X,Mat A,Mat B,void *ctx)
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
	
	ierr = PetscObjectTypeCompare((PetscObject)A,MATMFFD, &is_mffd);CHKERRQ(ierr);
	ierr = PetscObjectTypeCompare((PetscObject)A,MATNEST, &is_nest);CHKERRQ(ierr);
	ierr = PetscObjectTypeCompare((PetscObject)A,MATSHELL,&is_shell);CHKERRQ(ierr);

	if (is_nest) {
		Mat Auu;
		
		ierr = MatCreateSubMatrix(A,is[0],is[0],MAT_INITIAL_MATRIX,&Auu);CHKERRQ(ierr);

		is_shell = PETSC_FALSE;
		ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
			ierr = MatAssemble_StokesA_AUU(Auu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		ierr = MatDestroy(&Auu);CHKERRQ(ierr);
	}
	/* If shell, do nothing */
	/* If mffd,  do nothing */
	
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	/* preconditioner for Jacobian */
	{
		Mat Buu,Bpp;
		
		ierr = MatCreateSubMatrix(B,is[0],is[0],MAT_INITIAL_MATRIX,&Buu);CHKERRQ(ierr);
		ierr = MatCreateSubMatrix(B,is[1],is[1],MAT_INITIAL_MATRIX,&Bpp);CHKERRQ(ierr);
		
//		ierr = Assemble_Stokes_A11_Q2(user,dau,u,dap,p,Buu);CHKERRQ(ierr);
//		ierr = Assemble_Stokes_B22_P1(user,dau,u,dap,p,Bpp);CHKERRQ(ierr);
	
		is_shell = PETSC_FALSE;
		ierr = PetscObjectTypeCompare((PetscObject)Buu,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Buu);CHKERRQ(ierr);
			ierr = MatAssemble_StokesA_AUU(Buu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		is_shell = PETSC_FALSE;
		ierr = PetscObjectTypeCompare((PetscObject)Bpp,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Bpp);CHKERRQ(ierr);
			ierr = MatAssemble_StokesPC_ScaledMassMatrix(Bpp,dau,dap,user->stokes_ctx->p_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		ierr = MatDestroy(&Buu);CHKERRQ(ierr);
		ierr = MatDestroy(&Bpp);CHKERRQ(ierr);		
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
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
#define __FUNCT__ "pTatin3d_material_points"
PetscErrorCode pTatin3d_material_points(int argc,char **argv)
{
	DM        multipys_pack,dav;
	pTatinCtx user;
	Mat       A,B;
	Vec       X,F;
	IS        *is;
	SNES      snes;
	KSP       ksp;
	PC        pc;
	PetscReal dt;
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
	ierr = pTatinLogBasicDMDA(user,"velDMDA",dav);CHKERRQ(ierr);
	
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

	
	/* work vector for solution and residual */
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
	
	{
		Vec velocity,pressure;
		
		ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		ierr = BCListInsert(user->stokes_ctx->u_bclist,velocity);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	}

	/* INITIAL CONDITION */
	ierr = pTatinLogBasic(user);
	user->step = 0;
	user->time = 0.0;
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
	/* operators */
	ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
	ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,PETSC_FALSE,PETSC_TRUE,PETSC_TRUE,&B);CHKERRQ(ierr);
	ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
	/* configure for fieldsplit */
	ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
	ierr = KSPMonitorSet(ksp,pTatinStokesKSPMonitorBlocks,(void*)user,NULL);CHKERRQ(ierr);

	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = DMCompositeGetGlobalISs(multipys_pack,&is);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc,"u",is[0]);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc,"p",is[1]);CHKERRQ(ierr);
	ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);

	ierr = pTatinLogBasicSNES(user,"ic_stk",snes);CHKERRQ(ierr);
	ierr = pTatinLogBasicStokesSolution(user,multipys_pack,X);CHKERRQ(ierr);
	ierr = pTatinLogBasicStokesSolutionResiduals(user,snes,multipys_pack,X);CHKERRQ(ierr);
	
	/* COMPUTE DT */
	dt = 0.1;
	user->dt = dt;
	
	/* OUTPUT */
	ierr = pTatinModel_Output(user->model,user,X,"step000000");CHKERRQ(ierr);
	/* CHECKPOINT */
	ierr = pTatin3dCheckpointManager(user,X);CHKERRQ(ierr);
	
	while (user->step < user->nsteps) {
		char prefix[PETSC_MAX_PATH_LEN];

				
		/* UPDATE */
		ierr = pTatinModel_UpdateMeshGeometry(user->model,user,X);CHKERRQ(ierr);
				
		/* SOLVE */
		ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);

		/* increment time step */
		user->time += user->dt;
		user->step++;
		
		PetscPrintf(PETSC_COMM_WORLD,"[Step %1.6D:] time = %1.4e dt = %1.4e \n",user->step,user->time,user->dt);
		ierr = pTatinLogBasic(user);
		
		
		/* COMPUTE DT */

		
		
		/* ------------------- */
		/* OUTPUT */
		if (user->step%user->output_frequency==0) {
			PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step%1.6D",user->step);
			ierr = pTatinModel_Output(user->model,user,X,prefix);CHKERRQ(ierr);
		}
		/* CHECKPOINT */
		ierr = pTatin3dCheckpointManager(user,X);CHKERRQ(ierr);

		/* Terminate time stepping */
		if (user->time >= user->time_max) {
			break;
		}
	}
	
	
	ierr = ISDestroy(&is[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is[1]);CHKERRQ(ierr);
	ierr = PetscFree(is);CHKERRQ(ierr);

	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	ierr = SNESDestroy(&snes);CHKERRQ(ierr);

	
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = pTatin3dDestroyContext(&user);

	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_material_points_restart"
PetscErrorCode pTatin3d_material_points_restart(int argc,char **argv)
{
	pTatinCtx      user;
	PetscErrorCode ierr;
	Vec X;
	DM        multipys_pack,dav;
	IS        *is_stokes_field;

	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
	
	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	
	/* call model init */
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
	user->stokes_ctx->dav->ops->coarsenhierarchy = DMCoarsenHierarchy2_DA;
	dav                                          = user->stokes_ctx->dav;
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);

	
	
	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	user->pack = user->stokes_ctx->stokes_pack;
	/* fetch some local variables */
	multipys_pack = user->pack;
	/* IF I DON'T DO THIS, THE IS's OBTAINED FROM DMCompositeGetGlobalISs() are wrong !! */
	ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = DMCompositeGetGlobalISs(multipys_pack,&is_stokes_field);CHKERRQ(ierr);	
	
	

	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	
	
	while (user->step < user->nsteps) {
		char prefix[PETSC_MAX_PATH_LEN];
		
		/* UPDATE */
		
		/* SOLVE */
		
		/* increment time step */
		user->time += user->dt;
		user->step++;
		
		PetscPrintf(PETSC_COMM_WORLD,"[Step %1.6D:] time = %1.4e dt = %1.4e \n",user->step,user->time,user->dt);
		
		
		/* COMPUTE DT */
		
		
		/* ------------------- */
		/* OUTPUT */
		if (user->step%user->output_frequency==0) {
			PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step%1.6D",user->step);
			ierr = pTatinModel_Output(user->model,user,X,prefix);CHKERRQ(ierr);
		}
		/* CHECKPOINT */
		ierr = pTatin3dCheckpointManager(user,X);CHKERRQ(ierr);
		
		/* Terminate time stepping */
		if (user->time >= user->time_max) {
			break;
		}
	}
	
	
	
	
	
	ierr = ISDestroy(&is_stokes_field[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is_stokes_field[1]);CHKERRQ(ierr);
	ierr = PetscFree(is_stokes_field);CHKERRQ(ierr);
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = pTatin3dDestroyContext(&user);
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscBool restart;
	PetscErrorCode ierr;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	restart = PETSC_FALSE;
	PetscOptionsGetBool(NULL,NULL,"-test_restart",&restart,0);
	if (restart) {
		ierr = pTatin3d_material_points_restart(argc,argv);CHKERRQ(ierr);
	} else {
		ierr = pTatin3d_material_points(argc,argv);CHKERRQ(ierr);
	}
	
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

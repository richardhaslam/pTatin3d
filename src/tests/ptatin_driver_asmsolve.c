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
 **    filename:   ptatin_driver_asmsolve.c
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
#include "stokes_form_function.h"
#include "stokes_operators.h"
#include "stokes_operators_mf.h"
#include "stokes_assembly.h"
#include "stokes_output.h"
#include "dmda_element_q2p1.h"
#include "dmda_duplicate.h"
#include "dmda_project_coords.h"
#include "monitors.h"

typedef enum { OP_TYPE_REDISC_ASM=0, OP_TYPE_REDISC_MF, OP_TYPE_GALERKIN } OperatorType;

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

//    ierr = Assemble_Stokes_A11_Q2(user,dau,u,dap,p,Buu);CHKERRQ(ierr);
//    ierr = Assemble_Stokes_B22_P1(user,dau,u,dap,p,Bpp);CHKERRQ(ierr);

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
  //  ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
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



  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);

  /* operators */
  ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
//  ierr = StokesQ2P1CreateMatrixNest_Operator(user->stokes_ctx,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE,&A);CHKERRQ(ierr);

  ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,PETSC_FALSE,PETSC_TRUE,PETSC_TRUE,&B);CHKERRQ(ierr);

//  ierr = SNESSetJacobian(snes,B,B,FormJacobian_Stokes,user);CHKERRQ(ierr); //-snes_mf_operator -pc_fieldsplit_real_diagonal
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


  /* test viewer */
  ierr = pTatinModel_Output(user->model,user,X,"model");CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(user,X,"up");CHKERRQ(ierr);
  /* test generic viewer */
  {
    //    const int nf=2;
    //    const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes };
    const int nf=1;
    const MaterialPointField mp_prop_list[] = { MPField_Stokes };
    ierr = SwarmViewGeneric_ParaView(user->materialpoint_db,nf,mp_prop_list,user->outputpath,"stokesmp");CHKERRQ(ierr);
  }

  DataBucketView(PetscObjectComm((PetscObject)multipys_pack), user->materialpoint_db,"materialpoint_stokes",DATABUCKET_VIEW_STDOUT);



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

PetscErrorCode pTatin3d_galerkin_mg_material_points(int argc,char **argv)
{
  DM        multipys_pack,dav;
  pTatinCtx user;
  Mat       A,B;
  Vec       X,F;
  IS        *is;
  SNES      snes;
  KSP       ksp;
  PC        pc;

  PetscMPIInt    rank;
  DM             *dav_hierarchy;
  Mat            *interpolatation;
  PetscInt       nlevels,k;

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
  //  ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
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
  PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%D x %D x %D) : MG levels %D  \n", user->mx,user->my,user->mz,nlevels );
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
  }
  /* set to be q2 */
  for (k=0; k<nlevels-1; k++) {
    ierr = DMDASetElementType_Q2(dav_hierarchy[k]);CHKERRQ(ierr);
  }

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

  /* work vector for solution and residual */
  ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);

  /* operators */
  ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
  //  ierr = StokesQ2P1CreateMatrixNest_Operator(user->stokes_ctx,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE,&A);CHKERRQ(ierr);

  ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,PETSC_FALSE,PETSC_TRUE,PETSC_TRUE,&B);CHKERRQ(ierr);

  //  ierr = SNESSetJacobian(snes,B,B,FormJacobian_Stokes,user);CHKERRQ(ierr); //-snes_mf_operator -pc_fieldsplit_real_diagonal
  ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* configure for fieldsplit */
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPMonitorSet(ksp,pTatinStokesKSPMonitorBlocks,(void*)user,NULL);CHKERRQ(ierr);


  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

  ierr = DMCompositeGetGlobalISs(multipys_pack,&is);CHKERRQ(ierr);
  ierr = PCFieldSplitSetIS(pc,"u",is[0]);CHKERRQ(ierr);
  ierr = PCFieldSplitSetIS(pc,"p",is[1]);CHKERRQ(ierr);

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
    ierr = PCMGSetGalerkin(pc_i,PC_MG_GALERKIN_BOTH);CHKERRQ(ierr); /* OUCH - GALERKIN */

    for( k=1; k<nlevels; k++ ){
      ierr = PCMGSetInterpolation(pc_i,k,interpolatation[k]);CHKERRQ(ierr);
    }

  }


  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);



  /* test viewer */
  ierr = pTatinModel_Output(user->model,user,X,"model");CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(user,X,"up");CHKERRQ(ierr);
  /* test generic viewer */
  {
    //    const int nf=2;
    //    const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes };
    const int nf=1;
    const MaterialPointField mp_prop_list[] = { MPField_Stokes };
    ierr = SwarmViewGeneric_ParaView(user->materialpoint_db,nf,mp_prop_list,user->outputpath,"stokesmp");CHKERRQ(ierr);
  }

  DataBucketView(PetscObjectComm((PetscObject)multipys_pack), user->materialpoint_db,"materialpoint_stokes",DATABUCKET_VIEW_STDOUT);



  for (k=1; k<nlevels; k++) {
    ierr = MatDestroy(&interpolatation[k]);CHKERRQ(ierr);
  }
  ierr = PetscFree(interpolatation);CHKERRQ(ierr);

  for (k=0; k<nlevels; k++) {
    ierr = DMDestroy(&dav_hierarchy[k]);CHKERRQ(ierr);
  }
  ierr = PetscFree(dav_hierarchy);CHKERRQ(ierr);

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

PetscErrorCode pTatin3d_gmg_material_points(int argc,char **argv)
{
  DM        multipys_pack,dav;
  pTatinCtx user;
  Mat       A,B,operatorA11[10];
  Vec       X,F;
  IS        *is_stokes_field;
  SNES      snes;
  KSP       ksp;
  PC        pc;

  PetscMPIInt    rank;
  DM             *dav_hierarchy;
  Mat            *interpolatation;
  PetscInt       nlevels,k;
  Quadrature     volQ[10];
  BCList         u_bclist[10];

  Mat  R1[10];
  DM   scalarDA[10];
  Vec  sca[10];

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
  //  ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
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
  PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%D x %D x %D) : MG levels %D  \n", user->mx,user->my,user->mz,nlevels );
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
  }
  /* set to be q2 */
  for (k=0; k<nlevels-1; k++) {
    ierr = DMDASetElementType_Q2(dav_hierarchy[k]);CHKERRQ(ierr);
  }

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

  { /* test */
    for (k=0; k<nlevels; k++) {
      ierr = DMDADuplicateLayout(dav_hierarchy[k],1,2,DMDA_STENCIL_BOX,&scalarDA[k]);CHKERRQ(ierr);
      ierr = DMCreateGlobalVector(scalarDA[k],&sca[k]);CHKERRQ(ierr);
    }

    ierr = VecSetRandom(sca[nlevels-1],NULL);CHKERRQ(ierr);
    for (k=nlevels-1; k>=1; k--) {
      ierr = MatMAIJRedimension(interpolatation[k],1,&R1[k]);CHKERRQ(ierr);
      ierr = MatRestrict(R1[k],sca[k],sca[k-1]);CHKERRQ(ierr);
    }

    /* update markers = >> gauss points */
    {
      int               npoints;
      DataField         PField_std;
      DataField         PField_stokes;
      MPntStd           *mp_std;
      MPntPStokes       *mp_stokes;

      DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
      DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);

      DataBucketGetSizes(user->materialpoint_db,&npoints,NULL,NULL);
      mp_std    = PField_std->data; /* should write a function to do this */
      mp_stokes = PField_stokes->data; /* should write a function to do this */

      ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(user->coefficient_projection_type,npoints,mp_std,mp_stokes,nlevels,R1,dav_hierarchy,volQ);CHKERRQ(ierr);
    }
  }

  /* define boundary conditions */
  for (k=0; k<nlevels-1; k++) {
    PetscScalar zero = 0.0;

    ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);

    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

  }
  u_bclist[nlevels-1] = user->stokes_ctx->u_bclist;



  /* work vector for solution and residual */
  ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);

  /* operators */
  ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
  //  ierr = StokesQ2P1CreateMatrixNest_Operator(user->stokes_ctx,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE,&A);CHKERRQ(ierr);

  ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,PETSC_FALSE,PETSC_TRUE,PETSC_TRUE,&B);CHKERRQ(ierr);

  //  ierr = SNESSetJacobian(snes,B,B,FormJacobian_Stokes,user);CHKERRQ(ierr); //-snes_mf_operator -pc_fieldsplit_real_diagonal
  ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);


  /* generate the operators */
  for (k=0; k<nlevels-1; k++) {
    Mat Auu;

    ierr = DMSetMatType(dav_hierarchy[k],MATSBAIJ);CHKERRQ(ierr);
    ierr = DMCreateMatrix(dav_hierarchy[k],&Auu);CHKERRQ(ierr);
    ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);

    ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
    ierr = MatAssemble_StokesA_AUU(Auu,dav_hierarchy[k],u_bclist[k],volQ[k]);CHKERRQ(ierr);

    operatorA11[k] = Auu;
  }
  operatorA11[nlevels-1] = NULL;

  /* Fetch from the nest */
  ierr = DMCompositeGetGlobalISs(multipys_pack,&is_stokes_field);CHKERRQ(ierr);
  {
    Mat sub_A11;

    ierr = MatCreateSubMatrix(B,is_stokes_field[0],is_stokes_field[0],MAT_INITIAL_MATRIX,&sub_A11);CHKERRQ(ierr);
    operatorA11[nlevels-1] = sub_A11;
    ierr = MatDestroy(&sub_A11);CHKERRQ(ierr);

    ierr = PetscObjectReference((PetscObject)operatorA11[nlevels-1]);CHKERRQ(ierr);
  }



  /* configure for fieldsplit */
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPMonitorSet(ksp,pTatinStokesKSPMonitorBlocks,(void*)user,NULL);CHKERRQ(ierr);


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
    ierr = PCMGSetLevels(pc_i,nlevels,NULL);CHKERRQ(ierr);
    ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
    ierr = PCMGSetGalerkin(pc_i,PC_MG_GALERKIN_NONE);CHKERRQ(ierr);
    ierr = PCSetDM(pc_i,NULL);CHKERRQ(ierr);

    for( k=1; k<nlevels; k++ ){
      ierr = PCMGSetInterpolation(pc_i,k,interpolatation[k]);CHKERRQ(ierr);
    }

    /* drop the operators in */
    ierr = PCMGGetCoarseSolve(pc_i,&ksp_coarse);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp_coarse,operatorA11[0],operatorA11[0]);CHKERRQ(ierr);
    for( k=1; k<nlevels; k++ ){
      ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);
      ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k]);CHKERRQ(ierr);
    }

  }

  ierr = pTatinModel_ApplyInitialSolution(user->model,user,X);CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);



  /* test viewer */
  ierr = pTatinModel_Output(user->model,user,X,"model");CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(user,X,"up");CHKERRQ(ierr);
  /* test generic viewer */
  {
    //    const int nf=2;
    //    const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes };
    const int nf=1;
    const MaterialPointField mp_prop_list[] = { MPField_Stokes };
    ierr = SwarmViewGeneric_ParaView(user->materialpoint_db,nf,mp_prop_list,user->outputpath,"stokesmp");CHKERRQ(ierr);
  }

  DataBucketView(PetscObjectComm((PetscObject)multipys_pack), user->materialpoint_db,"materialpoint_stokes",DATABUCKET_VIEW_STDOUT);

  for (k=1; k<nlevels; k++) {
    ierr = MatDestroy(&R1[k]);CHKERRQ(ierr);
  }
  for (k=0; k<nlevels; k++) {
    ierr = DMDestroy(&scalarDA[k]);CHKERRQ(ierr);
  }
  for (k=0; k<nlevels; k++) {
    ierr = VecDestroy(&sca[k]);CHKERRQ(ierr);
  }

  for (k=0; k<nlevels-1; k++) {
    ierr = BCListDestroy(&u_bclist[k]);CHKERRQ(ierr);
  }
  for (k=0; k<nlevels-1; k++) {
    ierr = QuadratureDestroy(&volQ[k]);CHKERRQ(ierr);
  }
  for (k=0; k<nlevels; k++) {
    ierr = MatDestroy(&operatorA11[k]);CHKERRQ(ierr);
  }

  for (k=1; k<nlevels; k++) {
    ierr = MatDestroy(&interpolatation[k]);CHKERRQ(ierr);
  }
  ierr = PetscFree(interpolatation);CHKERRQ(ierr);

  for (k=0; k<nlevels; k++) {
    ierr = DMDestroy(&dav_hierarchy[k]);CHKERRQ(ierr);
  }
  ierr = PetscFree(dav_hierarchy);CHKERRQ(ierr);

  ierr = ISDestroy(&is_stokes_field[0]);CHKERRQ(ierr);
  ierr = ISDestroy(&is_stokes_field[1]);CHKERRQ(ierr);
  ierr = PetscFree(is_stokes_field);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = pTatin3dDestroyContext(&user);

  PetscFunctionReturn(0);
}

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
  PetscInt  level_type[10];

  PetscMPIInt    rank;
  DM             dav_hierarchy[10];
  Mat            interpolatation_v[10],interpolatation_eta[10];
  PetscInt       nlevels,k,max;
  Quadrature     volQ[10];
  BCList         u_bclist[10];
  PetscInt       kk,e;

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

  /* interpolate point coordinates (needed if mesh was modified) */
  //ierr = QuadratureStokesCoordinateSetUp(user->stokes_ctx->Q,dav);CHKERRQ(ierr);
  for (e=0; e<HEX_EDGES; e++) {
    ierr = SurfaceQuadratureGeometrySetUpStokes(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
    //ierr = SurfaceQuadratureOrientationViewGnuplotStokes(user->stokes_ctx->surfQ[e],dav,"def");CHKERRQ(ierr);
    ierr = SurfaceQuadratureViewParaview_Stokes(user->stokes_ctx,user->outputpath,"def");CHKERRQ(ierr);
  }
  /* interpolate material point coordinates (needed if mesh was modified) */
  ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);

  /* material geometry */
  ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
  //ierr = MaterialPointPopulationControl_v1(user);CHKERRQ(ierr);

  /* boundary conditions */
  ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);

  /* set up mg */
  user->stokes_ctx->dav->ops->coarsenhierarchy = DMCoarsenHierarchy2_DA;

  nlevels = 1;
  PetscOptionsGetInt(NULL,NULL,"-dau_nlevels",&nlevels,0);
  PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%D x %D x %D) : MG levels %D  \n", user->mx,user->my,user->mz,nlevels );
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
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
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

  /* inject coordinates */
  ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);

  /* Define interpolation operators for velocity space */
  interpolatation_v[0] = NULL;
  for (k=0; k<nlevels-1; k++) {
    ierr = DMCreateInterpolation(dav_hierarchy[k],dav_hierarchy[k+1],&interpolatation_v[k+1],NULL);CHKERRQ(ierr);
  }

  /* Define interpolation operators for scalar space */
  interpolatation_eta[0] = NULL;
  for (k=1; k<nlevels; k++) {
    ierr = MatMAIJRedimension(interpolatation_v[k],1,&interpolatation_eta[k]);CHKERRQ(ierr);
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
  volQ[nlevels-1] = user->stokes_ctx->volQ;

  /* update markers = >> gauss points */
  {
    int               npoints;
    DataField         PField_std;
    DataField         PField_stokes;
    MPntStd           *mp_std;
    MPntPStokes       *mp_stokes;

    DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
    DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);

    DataBucketGetSizes(user->materialpoint_db,&npoints,NULL,NULL);
    mp_std    = PField_std->data; /* should write a function to do this */
    mp_stokes = PField_stokes->data; /* should write a function to do this */

    ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(user->coefficient_projection_type,npoints,mp_std,mp_stokes,nlevels,interpolatation_eta,dav_hierarchy,volQ);CHKERRQ(ierr);
    /* alternative hierachy construction */
    /*
    // a
    for (k=0; k<nlevels; k++) {
      ierr = MProjection_Q1Projection_onto_Q2_MPntPStokes_Level(npoints,mp_std,mp_stokes,nlevels,dav_hierarchy,k,volQ[k]);CHKERRQ(ierr);
    }
    */
    /*
    // b
    for (k=0; k<nlevels; k++) {
      ierr = MProjection_P0Projection_onto_Q2_MPntPStokes_Level(npoints,mp_std,mp_stokes,nlevels,dav_hierarchy,k,volQ[k]);CHKERRQ(ierr);
    }
    */
    /*
    //
    for (k=nlevels-1; k<nlevels; k++) {
      ierr = MProjection_P0Projection_onto_Q2_MPntPStokes_Level(npoints,mp_std,mp_stokes,nlevels,dav_hierarchy,k,volQ[k]);CHKERRQ(ierr);
    }
    */
  }

  /* define boundary conditions - HARDCODED */
  /*
  for (k=0; k<nlevels-1; k++) {
    PetscScalar zero = 0.0;

    ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);

    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(u_bclist[k],dav_hierarchy[k],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

  }
  */
  for (k=0; k<nlevels-1; k++) {
    ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);
  }
  u_bclist[nlevels-1] = user->stokes_ctx->u_bclist;

  ierr = pTatinModel_ApplyBoundaryConditionMG(nlevels,u_bclist,dav_hierarchy,user->model,user);CHKERRQ(ierr);


  /* ============================================== */
  /* CONFIGURE OPERATORS */

  /* A operator */
  ierr = StokesQ2P1CreateMatrix_Operator(user->stokes_ctx,&A);CHKERRQ(ierr);
  /* memory saving - only need daU IF you want to split A11 into A11uu,A11vv,A11ww */
  {
    MatStokesMF mf;

    ierr = MatShellGetMatStokesMF(A,&mf);CHKERRQ(ierr);
    ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
    mf->daU = NULL;
  }

  /* B operator */
//  ierr = StokesQ2P1CreateMatrixNest_PCOperator(user->stokes_ctx,PETSC_FALSE,PETSC_TRUE,PETSC_TRUE,&B);CHKERRQ(ierr);
  {
    Mat         Aup,Apu,Spp,bA[2][2];
    MatStokesMF StkCtx;

    ierr = MatShellGetMatStokesMF(A,&StkCtx);CHKERRQ(ierr);

    /* Schur complement */
    ierr = DMSetMatType(dap,MATSBAIJ);CHKERRQ(ierr);
    ierr = DMCreateMatrix(dap,&Spp);CHKERRQ(ierr);
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
    bA[0][0] = NULL; bA[0][1] = Aup;
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



  /* A11 operator */

  /* OPTION 1 - JUST ASSEMBLE EVERYTHING */
/*
  for (k=0; k<nlevels; k++) {
    Mat Auu;

    ierr = DMSetMatType(dav_hierarchy[k],MATSBAIJ);CHKERRQ(ierr);
    ierr = DMCreateMatrix(dav_hierarchy[k],&Auu);CHKERRQ(ierr);
    ierr = MatSetOptionsPrefix(Auu,"Buu_");CHKERRQ(ierr);
    ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);

    ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
    ierr = MatAssemble_StokesA_AUU(Auu,dav_hierarchy[k],u_bclist[k],volQ[k]);CHKERRQ(ierr);

    operatorA11[k] = Auu;

    operatorB11[k] = Auu;
    ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
  }
*/

  /* defaults */
  level_type[0] = (PetscInt)OP_TYPE_REDISC_ASM;
  for (k=1; k<nlevels; k++) {
    level_type[k] = (PetscInt)OP_TYPE_REDISC_MF;
  }

  max = nlevels;
  ierr = PetscOptionsGetIntArray(NULL,NULL,"-A11_operator_type",level_type,&max,&flg);CHKERRQ(ierr);
  for (k=nlevels-1; k>=0; k--) {

    switch ((OperatorType)level_type[k]) {

      case OP_TYPE_REDISC_ASM:
      {
        Mat Auu;
        PetscBool same1 = PETSC_FALSE,same2 = PETSC_FALSE,same3 = PETSC_FALSE;

        /* use -stk_velocity_da_mat_type sbaij or -Buu_da_mat_type sbaij */
        PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Re-discretisation :: assembled operator \n", k);
        ierr = DMSetMatType(dav_hierarchy[k],MATSBAIJ);CHKERRQ(ierr);
        ierr = DMCreateMatrix(dav_hierarchy[k],&Auu);CHKERRQ(ierr);
        ierr = MatSetOptionsPrefix(Auu,"Buu_");CHKERRQ(ierr);
        ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSBAIJ,&same1);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSEQSBAIJ,&same2);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATMPISBAIJ,&same3);CHKERRQ(ierr);
        if (same1||same2||same3) {
          ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
        }
        /* should move assembly into jacobian */
        ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
        ierr = MatAssemble_StokesA_AUU(Auu,dav_hierarchy[k],u_bclist[k],volQ[k]);CHKERRQ(ierr);

        operatorA11[k] = Auu;
        operatorB11[k] = Auu;
        ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);

      }
        break;

      case OP_TYPE_REDISC_MF:
      {
        Mat Auu;
        MatA11MF mf,A11Ctx;

        PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Re-discretisation :: matrix free operator \n", k);
        ierr = MatA11MFCreate(&A11Ctx);CHKERRQ(ierr);
        ierr = MatA11MFSetup(A11Ctx,dav_hierarchy[k],volQ[k],u_bclist[k]);CHKERRQ(ierr);

        ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
        ierr = MatShellGetMatA11MF(Auu,&mf);CHKERRQ(ierr);
        ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
        mf->daU = NULL;
        operatorA11[k] = Auu;

        {
          PetscBool use_low_order_geometry = PETSC_FALSE;

          ierr = PetscOptionsGetBool(NULL,NULL,"-use_low_order_geometry",&use_low_order_geometry,NULL);CHKERRQ(ierr);
          if (use_low_order_geometry==PETSC_TRUE) {
            Mat Buu;

            PetscPrintf(PETSC_COMM_WORLD,"Activiting low order A11 operator \n");
            ierr = StokesQ2P1CreateMatrix_MFOperator_A11LowOrder(A11Ctx,&Buu);CHKERRQ(ierr);
            ierr = MatShellGetMatA11MF(Buu,&mf);CHKERRQ(ierr);
            ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
            mf->daU = NULL;
            operatorB11[k] = Buu;

          } else {
            operatorB11[k] = Auu;
            ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
          }
        }


        ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
      }
        break;

      case OP_TYPE_GALERKIN:
      {
        Mat Auu;

        if (k==nlevels-1) {
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
        }

        PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Galerkin :: assembled operator \n", k);

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




  /* Set fine A11 into nest */
  ierr = MatNestSetSubMat(B,0,0,operatorA11[nlevels-1]);CHKERRQ(ierr);

  /* ============================================== */




  /* work vector for solution and residual */
  ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

  /* initial condition */
  ierr = pTatinModel_ApplyInitialSolution(user->model,user,X);CHKERRQ(ierr);

  /* boundary condition */
  {
    Vec velocity,pressure;

    ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    ierr = BCListInsert(user->stokes_ctx->u_bclist,velocity);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  }

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);

  ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* configure for fieldsplit */
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);

  ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_StdoutStokesResiduals3d,(void*)user,NULL);CHKERRQ(ierr);
//  ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_ParaviewStokesResiduals3d,(void*)user,NULL);CHKERRQ(ierr);

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
    ierr = PCMGSetLevels(pc_i,nlevels,NULL);CHKERRQ(ierr);
    ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
    ierr = PCMGSetGalerkin(pc_i,PC_MG_GALERKIN_NONE);CHKERRQ(ierr);
    ierr = PCSetDM(pc_i,NULL);CHKERRQ(ierr);

    for( k=1; k<nlevels; k++ ){
      ierr = PCMGSetInterpolation(pc_i,k,interpolatation_v[k]);CHKERRQ(ierr);
    }

    /* drop the operators in - i presume this will also need to be performed inside the jacobian each time the operators are modified */
    /* No - it looks like PCSetUp_MG will call set operators on all levels if the SetOperators was called on the finest, which should/is done by the SNES */
    ierr = PCMGGetCoarseSolve(pc_i,&ksp_coarse);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp_coarse,operatorA11[0],operatorA11[0]);CHKERRQ(ierr);
    for( k=1; k<nlevels; k++ ){
      PetscBool use_low_order_geometry = PETSC_FALSE;


      ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);

      // use A for smoother, B for residual
      ierr = PetscOptionsGetBool(NULL,NULL,"-use_low_order_geometry",&use_low_order_geometry,NULL);CHKERRQ(ierr);
      if (use_low_order_geometry==PETSC_TRUE) {
        ierr = KSPSetOperators(ksp_smoother,operatorB11[k],operatorB11[k]);CHKERRQ(ierr);

        //ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorB11[k]);CHKERRQ(ierr);
      } else {
        // Use A for smoother, lo
        ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k]);CHKERRQ(ierr);
      }

    }
  }

  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);





  {
    char name[PETSC_MAX_PATH_LEN];

    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"ic_step%.6D",user->step);
    ierr = pTatinModel_Output(user->model,user,X,name);CHKERRQ(ierr);
  }

  PetscPrintf(PETSC_COMM_WORLD,"[Initial condition] Timestep[%D]: time %lf Myr \n", user->step, user->time );
  for (kk=0; kk<user->nsteps; kk++) {
    PetscInt tk = user->step+1;

    /* do solve */

    /* update markers */



    user->time += 0.12;
    user->step++;
    PetscPrintf(PETSC_COMM_WORLD,"Timestep[%D] : Cycle[%D/%D] : time %lf Myr \n", tk, kk, user->nsteps-1, user->time );


    if ((kk+1)%5==0) {
      char name[PETSC_MAX_PATH_LEN];

      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"step%.6D",tk);
      ierr = pTatinModel_Output(user->model,user,X,name);CHKERRQ(ierr);
    }

    if ((kk+1)%10==0) {
      char name[PETSC_MAX_PATH_LEN];

      PetscPrintf(PETSC_COMM_WORLD,"  checkpointing ptatin :: Model timestep %D : time %lf Myr : cycle[%D/%D] \n", tk,user->time,kk, user->nsteps-1 );
      /* check point test */
      //  ierr = pTatin3dContextSave(user,"checkpoint.file");CHKERRQ(ierr);
      //  ierr = pTatin3dContextLoad(user,"checkpoint.file");CHKERRQ(ierr);

      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"step%.6D",tk);
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

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = pTatin3dDestroyContext(&user);

  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;

  ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);

//  ierr = pTatin3d_material_points(argc,argv);CHKERRQ(ierr);

/*
  ./ptatin_driver_asmsolve.app -lattice_layout_Nx 3 -lattice_layout_Ny 3 -lattice_layout_Nz 3 -ptatin_model viscous_sinker -output_path crap -mx 4 -my 4 -mz 4 -snes_view -pc_type fieldsplit -fieldsplit_u_ksp_type fgmres -fieldsplit_u_ksp_monitor -ksp_monitor_true_residual -ksp_type fgmres -fieldsplit_p_pc_type jacobi -mx 10 -my 10 -mz 10 -snes_view  -snes_type ksponly -snes_rtol 1.0e-6 -snes_atol 1.0e-6 -ksp_rtol 1.0e-6 -ksp_atol 1.0e-8   -model_viscous_sinker_eta1 1.0e4 -fieldsplit_u_ksp_rtol 1.0e-6 -fieldsplit_u_ksp_max_it 1000  -model_viscous_sinker_dx 0.25 -model_viscous_sinker_dy 0.25 -model_viscous_sinker_dz 0.25 -mx 30 -my 30 -mz 30 -pc_fieldsplit_type ADDITIVE -model_viscous_sinker_eta0 1.0e-4 -model_viscous_sinker_eta1 1.0 -mx 12 -my 12 -mz 12 -pc_fieldsplit_type schur -pc_fieldsplit_schur_factorization_type upper -fieldsplit_p_ksp_type preonly -dau_nlevels 2 -log_summary
*/
//  ierr = pTatin3d_galerkin_mg_material_points(argc,argv);CHKERRQ(ierr);

/*
  ./ptatin_driver_asmsolve.app -lattice_layout_Nx 3 -lattice_layout_Ny 3 -lattice_layout_Nz 3 -ptatin_model viscous_sinker -output_path crap -snes_view -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_factorization_type upper -fieldsplit_p_ksp_type preonly -fieldsplit_u_mg_coarse_pc_type cholesky -dau_nlevels 3 -mx 20 -my 20 -mz 20 -fieldsplit_u_ksp_type fgmres -fieldsplit_u_ksp_monitor -ksp_monitor_true_residual -snes_type ksponly -snes_rtol 1.0e-6 -snes_atol 1.0e-6 -ksp_rtol 1.0e-6 -ksp_atol 1.0e-8   -model_viscous_sinker_eta0 1.0e-8 -log_summary -fieldsplit_u_mg_levels_ksp_max_it 4 -ksp_type fgmres -fieldsplit_u_mg_levels_pc_type jacobi
*/
//  ierr = pTatin3d_gmg_material_points(argc,argv);CHKERRQ(ierr);

  ierr = pTatin3d_gmg2_material_points(argc,argv);CHKERRQ(ierr);


  ierr = pTatinFinalize();CHKERRQ(ierr);
  return 0;
}

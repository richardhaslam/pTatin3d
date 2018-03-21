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
 **    filename:   test_mp_advection.c
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

static const char help[] = "Advection perforamnce / profiling test.\n"
"  Options:\n"
"    -flow_field 0 : "
"    -flow_field 1 : "
"    -flow_field 2 : "
"    -use_model_vel_field : defines flow field for a given model - solves once, then advects without updating velocity field\n\n";


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
#include "mp_advection.h"
#include "material_point_popcontrol.h"
#include "output_paraview.h"
#include "dmda_iterator.h"
#include "pswarm.h"

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

PetscErrorCode test_mp_advection(int argc,char **argv)
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
  PetscInt       kk;

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
    PetscInt nels,nen;
    const PetscInt *els;
    PetscInt lmx,lmy,lmz,si,sj,sk;

    ierr = DMDAGetElements_pTatinQ2P1(dav_hierarchy[k],&nels,&nen,&els);CHKERRQ(ierr);


    ierr = DMDAGetSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"         level [%2D]: global Q2 elements (%D x %D x %D) \n", k,lmx,lmy,lmz );

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
  }

  /* define boundary conditions - HARDCODED */
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
          operatorB11[k] = Auu;
          ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);

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
      ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);
      ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorB11[k]);CHKERRQ(ierr);
    }
  }

  /* solve once, then advect */
  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"[Initial condition] Timestep[%D]: time %lf Myr \n", user->step, user->time );
  for (kk=0; kk<user->nsteps; kk++) {
    PetscInt tk = user->step+1;
    PetscReal timestep;

    /* do solve */

    /* update markers */
    /* 1 compute timestep */
    {
      Vec velocity,pressure;

      ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
      ierr = SwarmUpdatePosition_ComputeCourantStep(dav_hierarchy[nlevels-1],velocity,&timestep);CHKERRQ(ierr);

      ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"  timestep[%D] dt_courant = %1.4e \n", tk, timestep );

      user->dt = timestep;
      user->dt = 1.0e-1 * user->dt;
    }
    PetscPrintf(PETSC_COMM_WORLD,"  timestep[%D] dt = %1.4e \n", tk, user->dt );

    ierr = SwarmUpdateProperties_MPntStd(user->materialpoint_db,user,X);CHKERRQ(ierr);

    /* 2 Advect markers */
    {
      Vec velocity,pressure;

      ierr = DMCompositeGetAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
      ierr = MaterialPointStd_UpdateGlobalCoordinates(user->materialpoint_db,dav_hierarchy[nlevels-1],velocity,user->dt);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);

    }


    /* 3 Update local coordinates and communicate */
    ierr = MaterialPointStd_UpdateCoordinates(user->materialpoint_db,dav_hierarchy[nlevels-1],user->materialpoint_ex);CHKERRQ(ierr);

    /* add / remove points if cells are over populated or depleted of points */
    //    ierr = MaterialPointPopulationControl(user);CHKERRQ(ierr);


    user->time += user->dt;
    user->step++;
    PetscPrintf(PETSC_COMM_WORLD,"Timestep[%D] : Cycle[%D/%D] : time %lf Myr \n", tk, kk, user->nsteps-1, user->time );


    if ((kk+1)%user->output_frequency==0) {
      char name[PETSC_MAX_PATH_LEN];

      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"step%.6D",tk);
      ierr = pTatinModel_Output(user->model,user,X,name);CHKERRQ(ierr);
    }

    if ((kk+1)%user->output_frequency==0) {
      char name[PETSC_MAX_PATH_LEN];

      PetscPrintf(PETSC_COMM_WORLD,"  checkpointing ptatin :: Model timestep %D : time %lf Myr : cycle[%D/%D] \n", tk,user->time,kk, user->nsteps-1 );

      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"step%.6D",tk);
      ierr = pTatin3dCheckpoint(user,X,name);CHKERRQ(ierr);
    }

  }

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

PetscBool EvaluateTestVelocityField_0(PetscScalar coor[],PetscScalar *value,void *ctx)
{
  PetscInt dof;
  PetscReal fac,x,y,z;

  dof = *((PetscInt*)ctx);
  x = coor[0];
  y = coor[1];
  z = coor[2];

  switch (dof) {
    case 0:
      fac = 0.5*(cos(M_PI*x)+1.0);
      *value = fac * ( sin(y*3.3) + sin(z) );
      break;
    case 1:
      fac = 0.5*(cos(M_PI*y)+1.0);
      *value = fac * ( sin(y*1.3) * x + z );
      break;
    case 2:
      fac = 0.5*(cos(M_PI*z)+1.0);
      *value = fac * ( cos(x*1.1) * y*z + x );
      break;
  }

  return PETSC_TRUE;
}

#include <petsc/ptatin_petsc_ex43-solcx.h>
PetscBool EvaluateTestVelocityField_1(PetscScalar coor[],PetscScalar *value,void *ctx)
{
  PetscInt  dof;
  PetscReal x,y,z;
  PetscReal coorxy[2],coorxz[2],cooryz[2];
  PetscReal vxy[2],vxz[2],vyz[2];

  dof = *((PetscInt*)ctx);
  x = coor[0];
  y = coor[1];
  z = coor[2];

  coorxy[0] = x;
  coorxy[1] = y;

  coorxz[0] = x;
  coorxz[1] = x;

  cooryz[0] = y;
  cooryz[1] = z;

  evaluate_solCx(coorxy,
      1.0,50.0,   /* Input parameters: coord, viscosity_A, viscosity_B */
      0.5,1,       /* Input parameters: viscosity jump location, wave number in y */
      vxy,NULL,NULL,NULL);

  evaluate_solCx(coorxz,
      1.0,1.0,   /* Input parameters: coord, viscosity_A, viscosity_B */
      0.5,2,       /* Input parameters: viscosity jump location, wave number in y */
      vxz,NULL,NULL,NULL);

  evaluate_solCx(cooryz,
      1.0,1.0,   /* Input parameters: coord, viscosity_A, viscosity_B */
      0.5,1,       /* Input parameters: viscosity jump location, wave number in y */
      vyz,NULL,NULL,NULL);

  *value = 0.0;
  switch (dof) {
    case 0:
      *value = vxy[0];
      break;
    case 1:
      *value = vxy[1];// + 0.1*vyz[0];
      break;
    case 2:
      *value = 0.0;// + 0.1*vyz[1];
      break;
  }

  return PETSC_TRUE;
}

PetscErrorCode DefineTestVelocityField(PetscInt vfield_idx,DM dmv,Vec velocity)
{
  PetscErrorCode ierr;
  PetscInt dof_idx;

  switch (vfield_idx) {
    case 0:
      dof_idx = 0;
      ierr = DMDAVecTraverse3d(dmv,velocity,dof_idx,EvaluateTestVelocityField_0,(void*)&dof_idx);CHKERRQ(ierr);
      dof_idx = 1;
      ierr = DMDAVecTraverse3d(dmv,velocity,dof_idx,EvaluateTestVelocityField_0,(void*)&dof_idx);CHKERRQ(ierr);
      dof_idx = 2;
      ierr = DMDAVecTraverse3d(dmv,velocity,dof_idx,EvaluateTestVelocityField_0,(void*)&dof_idx);CHKERRQ(ierr);
      break;
    case 1:
      dof_idx = 0;
      ierr = DMDAVecTraverse3d(dmv,velocity,dof_idx,EvaluateTestVelocityField_1,(void*)&dof_idx);CHKERRQ(ierr);
      dof_idx = 1;
      ierr = DMDAVecTraverse3d(dmv,velocity,dof_idx,EvaluateTestVelocityField_1,(void*)&dof_idx);CHKERRQ(ierr);
      dof_idx = 2;
      ierr = DMDAVecTraverse3d(dmv,velocity,dof_idx,EvaluateTestVelocityField_1,(void*)&dof_idx);CHKERRQ(ierr);
      break;
    case 2:
      break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MaterialPointAdvectionTest2(void)
{
  pTatinCtx       user;
  PhysCompStokes  stokes;
  DM              multipys_pack,dmv,dmp,dmx;
  Vec             X;
  IS              *is_up_field;
  PetscInt        f,nfields;
  DataBucket      materialpoint_db;
  PetscInt        step;
  Vec             velocity,pressure;
  PetscReal       timestep;
  PetscReal       dt_factor = 1.0;
  PetscInt        vfield_idx = 0;
  PetscErrorCode  ierr;
  PSwarm          pswarm,*pswarm2,*psi;

  PetscFunctionBegin;

  ierr = PetscOptionsGetInt(NULL,NULL,"-flow_field",&vfield_idx,NULL);CHKERRQ(ierr);

  ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
  ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);

  ierr = pTatinModelRegisterAll();CHKERRQ(ierr); /* Register all models */

  /* Generate physics modules [Stokes] */
  ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr); /* Fetch local variables */
  ierr = PhysCompStokesGetDMs(stokes,&dmv,&dmp);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dmv,&dmx);CHKERRQ(ierr);

  /* Pack all physics together */
  user->pack = stokes->stokes_pack;
  PetscObjectReference((PetscObject)user->pack);
  multipys_pack = user->pack;

  /* -- NOTE -- DM{Get/Create}GlobalVector must be called before DMCompositeGetGlobalISs() */
  ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(multipys_pack,&X);CHKERRQ(ierr);

  ierr = DMCompositeGetNumberDM(multipys_pack,&nfields);CHKERRQ(ierr);
  ierr = DMCompositeGetGlobalISs(multipys_pack,&is_up_field);CHKERRQ(ierr);

  ierr = pTatin3dCreateMaterialPoints(user,dmv);CHKERRQ(ierr); /* Allocate data bucket for material points */
  ierr = pTatinGetMaterialPoints(user,&materialpoint_db,NULL);CHKERRQ(ierr);

  ierr = PSwarmCreate(PETSC_COMM_WORLD,&pswarm);CHKERRQ(ierr);
  ierr = PSwarmSetOptionsPrefix(pswarm,"passive_");CHKERRQ(ierr);
  ierr = PSwarmSetPtatinCtx(pswarm,user);CHKERRQ(ierr);
  ierr = PSwarmSetTransportModeType(pswarm,PSWARM_TM_LAGRANGIAN);CHKERRQ(ierr);

  ierr = PSwarmSetFromOptions(pswarm);CHKERRQ(ierr);

  ierr = PSwarmCreateMultipleInstances(PETSC_COMM_WORLD,&pswarm2);CHKERRQ(ierr);
  psi = &pswarm2[0];
  while (*psi && pswarm2) {
    ierr = PSwarmSetPtatinCtx(*psi,user);CHKERRQ(ierr);
    ierr = PSwarmSetTransportModeType(*psi,PSWARM_TM_LAGRANGIAN);CHKERRQ(ierr);
    ierr = PSwarmSetFromOptions(*psi);CHKERRQ(ierr);
    psi++;
  }

  //ierr = pTatinModel_ApplyInitialMeshGeometry(model,user);CHKERRQ(ierr); /* <<<< User call back >>>> */
  //ierr = DMDASetUniformCoordinates(dmv,-1.0,1.0,-1.0,1.0,-1.0,1.0);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(dmv,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);

  /* interpolate material point coordinates (needed if mesh was modified) */
  ierr = MaterialPointCoordinateSetUp(user,dmv);CHKERRQ(ierr); /* <<<< User call back >>>> */
  //ierr = pTatinModel_ApplyInitialMaterialGeometry(model,user);CHKERRQ(ierr);
  {
    DataField PField_std;
    int p,n_mp_points;

    DataBucketGetDataFieldByName(materialpoint_db,MPntStd_classname,&PField_std);
    DataFieldGetAccess(PField_std);
    DataBucketGetSizes(materialpoint_db,&n_mp_points,0,0);
    for (p=0; p<n_mp_points; p++) {
      MPntStd *material_point;
      int     region,regionIJK[3],II,JJ,KK;
      double  *position;
      double  ddx=0.5,ddy=0.5,ddz=0.5;

      DataFieldAccessPoint(PField_std,p,(void**)&material_point);
      MPntStdGetField_global_coord(material_point,&position);
      region = 0;
      II = (int)(position[0]/ddx);
      JJ = (int)(position[1]/ddy);
      KK = (int)(position[2]/ddz);

      regionIJK[0] = regionIJK[1] = regionIJK[2] = 0;
      if (II%2) {
        regionIJK[0] = 1;
      }
      if (JJ%2) {
        regionIJK[1] = 1;
      }
      if (KK%2) {
        regionIJK[2] = 1;
      }
      region = regionIJK[0] + regionIJK[1]*2 + regionIJK[2]*4;
      MPntStdSetField_phase_index(material_point,region);
    }
    DataFieldRestoreAccess(PField_std);
  }

  //ierr = pTatinModel_ApplyBoundaryCondition(model,user);CHKERRQ(ierr); /* <<<< User call back >>>> */

  /* work vector for solution and residual */
  ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);

  ierr = PSwarmAttachStateVecVelocityPressure(pswarm,X);CHKERRQ(ierr);
  psi = &pswarm2[0];
  while (*psi) {
    ierr = PSwarmAttachStateVecVelocityPressure(*psi,X);CHKERRQ(ierr);
    psi++;
  }

  /* insert boundary conditions into solution vector X */
  {
    BCList bclist_u;

    ierr = PhysCompStokesGetBCList(stokes,&bclist_u,NULL);CHKERRQ(ierr);
    ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    ierr = BCListInsert(bclist_u,velocity);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  }


  /* initial condition */
  //ierr = pTatinModel_ApplyInitialSolution(model,user,X);CHKERRQ(ierr); /* <<<< User call back >>>> */
  /* Define user velocity field */
  ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  ierr = DefineTestVelocityField(vfield_idx,dmv,velocity);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  /* output ic */
  ierr = pTatinOutputLiteParaViewMeshVelocity(multipys_pack,X,user->outputpath,"step000000_v");CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_MPntStd(user,"step000000");CHKERRQ(ierr);


  /* compute timestep */
  user->dt = 1.0e32;
  /* -- NOTE -- We can ignore last physics entry */
  ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  ierr = SwarmUpdatePosition_ComputeCourantStep(dmv,velocity,&timestep);CHKERRQ(ierr);
  ierr = pTatin_SetTimestep(user,"StkCourant",timestep);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"  timestep[stokes] dt_courant = %1.4e \n", user->dt );

  user->step = 1;
  user->time = user->time + user->dt;

  /* Begin time stepping */
  for (step=1; step <= user->nsteps; step++) {
    char  stepname[PETSC_MAX_PATH_LEN];

    PetscPrintf(PETSC_COMM_WORLD,"<<----------------------------------------------------------------------------------------------->>\n");
    PetscPrintf(PETSC_COMM_WORLD,"   [[ EXECUTING TIME STEP : %D ]]\n", step );
    PetscPrintf(PETSC_COMM_WORLD,"     dt    : %1.4e \n", user->dt );
    PetscPrintf(PETSC_COMM_WORLD,"     time  : %1.4e \n", user->time );

    ierr = pTatinLogBasic(user);CHKERRQ(ierr);

    /* update marker time dependent terms */
    /* e.g. e_plastic^1 = e_plastic^0 + dt * [ strain_rate_inv(u^0) ] */
    /*
       NOTE: for a consistent forward difference time integration we evaluate u^0 at x^0
       - thus this update is performed BEFORE we advect the markers
    */
    ierr = pTatin_UpdateCoefficientTemporalDependence_Stokes(user,X);CHKERRQ(ierr);

    /* update marker positions */
    ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    ierr = MaterialPointStd_UpdateGlobalCoordinates(materialpoint_db,dmv,velocity,user->dt);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);

    ierr = PSwarmFieldUpdateAll(pswarm);CHKERRQ(ierr);
    psi = &pswarm2[0];
    while (*psi) {
      ierr = PSwarmFieldUpdateAll(*psi);CHKERRQ(ierr);
      psi++;
    }

    /* update mesh */
    //ierr = pTatinModel_UpdateMeshGeometry(model,user,X);CHKERRQ(ierr);

    /* 3 Update local coordinates and communicate */
    ierr = MaterialPointStd_UpdateCoordinates(materialpoint_db,dmv,user->materialpoint_ex);CHKERRQ(ierr);

    /* 3a - Add material */
    //ierr = pTatinModel_ApplyMaterialBoundaryCondition(model,user);CHKERRQ(ierr);

    /* add / remove points if cells are over populated or depleted of points */
    ierr = MaterialPointPopulationControl_v1(user);CHKERRQ(ierr);

    /* output material points */
    if (step%user->output_frequency == 0) {
      PetscSNPrintf(stepname,PETSC_MAX_PATH_LEN-1,"step%1.6D",step);
      ierr = pTatin3d_ModelOutput_MPntStd(user,stepname);CHKERRQ(ierr);

      ierr = PSwarmView(pswarm,PSW_VT_SINGLETON);CHKERRQ(ierr);
      psi = &pswarm2[0];
      while (*psi) {
        ierr = PSwarmView(*psi,PSW_VT_SINGLETON);CHKERRQ(ierr);
        psi++;
      }
    }

    /* compute timestep */
    user->dt = 1.0e32;
    ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    ierr = SwarmUpdatePosition_ComputeCourantStep(dmv,velocity,&timestep);CHKERRQ(ierr);

    timestep = timestep/dt_factor;
    ierr = pTatin_SetTimestep(user,"StkCourant",timestep);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"  timestep_stokes[%D] dt_courant = %1.4e \n", step,user->dt );

    /* Terminate time stepping */
    if (user->time >= user->time_max) {
      break;
    }

    /* update time */
    user->step++;
    user->time = user->time + user->dt;
  }

  /* clean up */
  ierr = PSwarmViewInfo(pswarm);CHKERRQ(ierr);
  ierr = PSwarmDestroy(&pswarm);CHKERRQ(ierr);

  psi = &pswarm2[0];
  while (*psi) {
    PSwarm *next = psi + 1;

    ierr = PSwarmViewInfo(*psi);CHKERRQ(ierr);
    printf("PSwarmDestroy - freeing %p \n",*psi);
    ierr = PSwarmDestroy(psi);CHKERRQ(ierr);
    psi = next;
  }
  PetscFree(pswarm2);

  ierr = VecDestroy(&X);CHKERRQ(ierr);
  for (f=0; f<nfields; f++) {
    ierr = ISDestroy(&is_up_field[f]);CHKERRQ(ierr);
  }
  ierr = PetscFree(is_up_field);CHKERRQ(ierr);
  ierr = pTatin3dDestroyContext(&user);

  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscBool      model_vel_field = PETSC_FALSE;

  ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);


  ierr = PetscOptionsGetBool(NULL,NULL,"-use_model_vel_field",&model_vel_field,NULL);CHKERRQ(ierr);

  /*
     This requires a model to be specified.
     It performs a single solve, then continuously advects material points
  */
  if (model_vel_field) {
    ierr = test_mp_advection(argc,argv);CHKERRQ(ierr);
  } else {
    ierr = MaterialPointAdvectionTest2();CHKERRQ(ierr);
  }

  ierr = pTatinFinalize();CHKERRQ(ierr);
  return 0;
}

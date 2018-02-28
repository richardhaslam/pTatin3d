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
 **    filename:   ptatin_eigen_analysis.c
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
"3D prototype of the (p)ragmatic version of pTatin with eigen value analysis. (pTatin3d_v0.0)\n"
"Accepts -asm_dump and -block_dump flags to output PETSc binary matrices for external analysis\n\n";

#include "slepceps.h"

#include "petsc/private/dmdaimpl.h" 

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_log.h"

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
#include "mp_advection.h"
#include "eigen_operators.h"

static PetscErrorCode DumpA11(pTatinCtx user)
{
  Mat            A11;
  PhysCompStokes stk = user->stokes_ctx;
  DM             da = stk->dav;
  PetscViewer    viewer;
  PetscBool      same;
  const char     *filename="A11_out";
  PetscErrorCode ierr; 

  PetscFunctionBeginUser;
  ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
  ierr = DMCreateMatrix(da,&A11);CHKERRQ(ierr);

  ierr = PetscObjectTypeCompare((PetscObject)A11,MATSBAIJ,&same);CHKERRQ(ierr);
  if (same) {
    ierr = MatSetOption(A11,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
  }

  ierr = MatAssemble_StokesA_AUU(A11,da,stk->u_bclist,stk->volQ);CHKERRQ(ierr);
  {
    PetscInt m,n;
    ierr = MatGetSize(A11,&m,&n);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing %s of size %d x %d\n",filename,m,n);CHKERRQ(ierr);       
  }
  ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)A11),filename,FILE_MODE_WRITE,&viewer);
  ierr = MatView(A11,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = MatDestroy(&A11);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DumpA12(pTatinCtx user)
{
  PhysCompStokes stk=user->stokes_ctx;
  Mat            A12;
  PetscViewer    viewer;
  const char     *filename="A12_out";
  PetscErrorCode ierr; 

  PetscFunctionBeginUser;
  ierr = StokesQ2P1CreateMatrix_A12(stk,&A12);CHKERRQ(ierr);
  ierr = MatAssemble_StokesA_A12(A12,stk->dav,stk->dap,stk->u_bclist,stk->p_bclist,stk->volQ);CHKERRQ(ierr);
  {
    PetscInt m,n;
    ierr = MatGetSize(A12,&m,&n);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing %s of size %d x %d\n",filename,m,n);CHKERRQ(ierr);       
  }
  ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)A12),"A12_out",FILE_MODE_WRITE,&viewer);
  ierr = MatView(A12,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = MatDestroy(&A12);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode BlockDump(pTatinCtx user)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DumpA11(user);CHKERRQ(ierr);
  ierr = DumpA12(user);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* PS : quick and dirty viewing of subdomain mats */
static PetscErrorCode ASMDump(PC pc)
{
  PetscBool isfs;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscObjectTypeCompare((PetscObject)pc,PCFIELDSPLIT,&isfs);
  if(isfs){
    PetscInt n;
    KSP *subksp;
    PC subpc;
    /* If the PC type is fieldsplit, call this function on the (1,1) block */
    ierr = PCFieldSplitGetSubKSP(pc,&n,&subksp);CHKERRQ(ierr);
    /* Choose the first split, dangerously */
    ierr = KSPGetPC(subksp[0],&subpc);CHKERRQ(ierr);
    ierr = ASMDump(subpc);CHKERRQ(ierr);
  }else{ 
    /* Else if the PC type is MG, call this function on all the level PCs*/
    PetscBool ismg;
    KSP subksp;
    PC subpc;
    ierr = PetscObjectTypeCompare((PetscObject)pc,PCMG,&ismg);
    if(ismg){
      PetscInt nlevels,k;
      ierr = PCMGGetLevels(pc,&nlevels);CHKERRQ(ierr);
      for(k=0;k<nlevels;++k){
        ierr = PCMGGetSmoother(pc,k,&subksp);CHKERRQ(ierr);
        ierr = KSPGetPC(subksp,&subpc);CHKERRQ(ierr);
        ierr = ASMDump(subpc);CHKERRQ(ierr); 
      }
      if(nlevels > 1){
        ierr = PCMGGetCoarseSolve(pc,&subksp);CHKERRQ(ierr);
        ierr = KSPGetPC(subksp,&subpc);CHKERRQ(ierr);
        ierr = ASMDump(subpc);CHKERRQ(ierr); 
      }
    }else{
      /* Else, confirm that this is PCASM */
      PetscBool isasm;
      ierr = PetscObjectTypeCompare((PetscObject)pc,PCASM,&isasm);
      if(!isasm){
        PCType type;
        ierr = PCGetType(pc,&type);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"[Skipping level PC of type %s]\n",type);
      }else{
        /* Dump the local Blocks */
        KSP *subksp;
        PetscInt n_local,first_local,i;
        ierr = PCASMGetSubKSP(pc,&n_local,&first_local,&subksp);CHKERRQ(ierr);
        for(i=0;i<n_local;++i){ 
          Mat         Asub;
          PetscMPIInt rank;
          const char* prefix;
          char        filename[PETSC_MAX_PATH_LEN];
          PetscViewer viewer;

          MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
          ierr = KSPGetOperators(subksp[i],&Asub,NULL);CHKERRQ(ierr);
          ierr = PCGetOptionsPrefix(pc,&prefix);CHKERRQ(ierr);
          if(prefix) {
            PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%ssub%d_%d_out",prefix,rank,i);
          }else{ 
            PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"sub%d_%d_out",rank,i);
          }
          ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)Asub),filename,FILE_MODE_WRITE,&viewer);
          ierr = MatView(Asub,viewer);CHKERRQ(ierr);
          {
            PetscInt m,n;
            ierr = MatGetSize(Asub,&m,&n);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF,"Writing %s of size %d x %d\n",filename,m,n);CHKERRQ(ierr);       
          }
          ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
        }
      }
    }
  }
  PetscFunctionReturn(0); 
}

typedef enum { OP_TYPE_REDISC_ASM=0, OP_TYPE_REDISC_MF, OP_TYPE_GALERKIN } OperatorType;

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
  
  ierr = PetscObjectTypeCompare((PetscObject)(A),MATMFFD, &is_mffd);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)(A),MATNEST, &is_nest);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)(A),MATSHELL,&is_shell);CHKERRQ(ierr);
  
  if (is_nest) {
    Mat Auu;
    
    ierr = MatGetSubMatrix(A,is[0],is[0],MAT_INITIAL_MATRIX,&Auu);CHKERRQ(ierr);
    
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
    
    ierr = MatGetSubMatrix(B,is[0],is[0],MAT_INITIAL_MATRIX,&Buu);CHKERRQ(ierr);
    ierr = MatGetSubMatrix(B,is[1],is[1],MAT_INITIAL_MATRIX,&Bpp);CHKERRQ(ierr);
    
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
#define __FUNCT__ "pTatin3dStokesBuildMeshHierarchy"
PetscErrorCode pTatin3dStokesBuildMeshHierarchy(DM dav,PetscInt nlevels,DM dav_hierarchy[])
{
  PetscErrorCode ierr;
  DM *coarsened_list;
  PetscInt k;
  
  PetscFunctionBegin;
  
  /* set up mg */
  dav->ops->coarsenhierarchy = DMCoarsenHierarchy2_DA;
  
  dav_hierarchy[ nlevels-1 ] = dav;
  ierr = PetscObjectReference((PetscObject)dav);CHKERRQ(ierr);
  
  /* Coarsen nlevels - 1 times, and add levels into list so that level 0 is the coarsest */
  ierr = PetscMalloc(sizeof(DM)*(nlevels-1),&coarsened_list);CHKERRQ(ierr);
  ierr = DMCoarsenHierarchy(dav,nlevels-1,coarsened_list);CHKERRQ(ierr);
  for (k=0; k<nlevels-1; k++) {
    dav_hierarchy[ nlevels-2-k ] = coarsened_list[k];
  }
  PetscFree(coarsened_list);
  
  /* Set all dav's to be of type Q2 */
  for (k=0; k<nlevels-1; k++) {
    ierr = PetscObjectSetOptionsPrefix((PetscObject)dav_hierarchy[k],"stk_velocity_");CHKERRQ(ierr);
    ierr = DMDASetElementType_Q2(dav_hierarchy[k]);CHKERRQ(ierr);
  }
  
  /* inject coordinates */
  ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dStokesReportMeshHierarchy"
PetscErrorCode pTatin3dStokesReportMeshHierarchy(PetscInt nlevels,DM dav_hierarchy[])
{
  PetscErrorCode ierr;
  PetscInt       k,lmx,lmy,lmz,si,sj,sk;
  PetscInt       nels,nen;
  const PetscInt *els;
  PetscMPIInt    rank;
  
  PetscFunctionBegin;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  /* Report mesh sizes */
  for (k=0; k<nlevels; k++) {
    ierr = DMDAGetSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"         level [%2D]: global Q2 elements (%D x %D x %D) \n", k,lmx,lmy,lmz );
  }   
  
  for (k=0; k<nlevels; k++) {
    
    ierr = DMDAGetElements_pTatinQ2P1(dav_hierarchy[k],&nels,&nen,&els);CHKERRQ(ierr);
    ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
    if (rank<10) {
      PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: local Q2 elements  (%D x %D x %D) \n", rank, k,lmx,lmy,lmz );
    }
  }
  
  for (k=0; k<nlevels; k++) {
    ierr = DMDAGetElements_pTatinQ2P1(dav_hierarchy[k],&nels,&nen,&els);CHKERRQ(ierr);
    
    ierr = DMDAGetCornersElementQ2(dav_hierarchy[k],&si,&sj,&sk,&lmx,&lmy,&lmz);CHKERRQ(ierr);
    si = si/2;
    sj = sj/2;
    sk = sk/2;
    if (rank<10) {
      PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: element range [%D - %D] x [%D - %D] x [%D - %D] \n", rank, k,si,si+lmx-1,sj,sj+lmy-1,sk,sk+lmz-1 );
    }
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dCreateStokesOperators"
PetscErrorCode pTatin3dCreateStokesOperators(PhysCompStokes stokes_ctx,IS is_stokes_field[],
                                             PetscInt nlevels,DM dav_hierarchy[],Mat interpolation_v[],
                                             BCList u_bclist[],Quadrature volQ[],
                                             Mat *_A,Mat operatorA11[],Mat *_B,Mat operatorB11[])
{
  Mat            A,B;
  OperatorType   level_type[10];
  DM             dap;
  PetscInt       k,max;
  PetscBool      flg;
  static int     been_here = 0;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  dap = stokes_ctx->dap;
  
  /* A operator */
  ierr = StokesQ2P1CreateMatrix_Operator(stokes_ctx,&A);CHKERRQ(ierr);
  /* memory saving - only need daU IF you want to split A11 into A11uu,A11vv,A11ww */
  {
    MatStokesMF mf;
    
    ierr = MatShellGetMatStokesMF(A,&mf);CHKERRQ(ierr);
    ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
    mf->daU = NULL;
  }
  
  /* B operator */
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
    
    ierr = MatCreateNest(PETSC_COMM_WORLD,2,is_stokes_field,2,is_stokes_field,&bA[0][0],&B);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    /* tidy up - hand back destruction to B */
    ierr = MatDestroy(&Aup);CHKERRQ(ierr);
    ierr = MatDestroy(&Apu);CHKERRQ(ierr);
    ierr = MatDestroy(&Spp);CHKERRQ(ierr);
  }
  
  /* A11 operator */  
  /* defaults */
  level_type[0] = OP_TYPE_REDISC_ASM;
  for (k=1; k<nlevels; k++) {
    level_type[k] = OP_TYPE_REDISC_MF;
  }
  
  max = nlevels;
  ierr = PetscOptionsGetIntArray(NULL,NULL,"-A11_operator_type",(PetscInt*)level_type,&max,&flg);CHKERRQ(ierr);
  for (k=nlevels-1; k>=0; k--) {
    
    switch (level_type[k]) {
        
      case 0:
      {
        Mat Auu;
        PetscBool same1 = PETSC_FALSE,same2 = PETSC_FALSE,same3 = PETSC_FALSE;
        
        /* use -stk_velocity_da_mat_type sbaij or -Buu_da_mat_type sbaij */
        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Re-discretisation :: assembled operator \n", k);
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
        
      case 1:
      {
        Mat Auu;
        MatA11MF mf,A11Ctx;
        
        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Re-discretisation :: matrix free operator \n", k);
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
            
            if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Activating low order A11 operator \n");
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
        
      case 2:
      {
        Mat Auu;
        
        if (k==nlevels-1) {
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
        } 
        
        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Galerkin :: assembled operator \n", k);
        
        /* should move coarse grid assembly into jacobian */
        ierr = MatPtAP(operatorA11[k+1],interpolation_v[k+1],MAT_INITIAL_MATRIX,1.0,&Auu);CHKERRQ(ierr);
        
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
  
  *_A = A;
  *_B = B;
  
  been_here = 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dStokesKSPConfigureFSGMG"
PetscErrorCode pTatin3dStokesKSPConfigureFSGMG(KSP ksp,PetscInt nlevels,Mat operatorA11[],Mat operatorB11[],Mat interpolation_v[])
{
  PetscInt k,nsplits;
  PC       pc,pc_i;
  KSP      *sub_ksp,ksp_coarse,ksp_smoother;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
  
  ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
  ierr = PCSetType(pc_i,PCMG);CHKERRQ(ierr);
  ierr = PCMGSetLevels(pc_i,nlevels,NULL);CHKERRQ(ierr);
  ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
  ierr = PCMGSetGalerkin(pc_i,PETSC_FALSE);CHKERRQ(ierr);
  
  for( k=1; k<nlevels; k++ ){
    ierr = PCMGSetInterpolation(pc_i,k,interpolation_v[k]);CHKERRQ(ierr);
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
  PetscFunctionReturn(0);
}

/* analysers */
#undef __FUNCT__  
#define __FUNCT__ "_slepc_eigen"
PetscErrorCode _slepc_eigen(Mat A,EPSProblemType prob_type,const char description[])
{
  EPS            eps;         /* eigenproblem solver context */
  EPSType  type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  Vec            xr,xi;
  PetscInt       i,nev,maxit,its,nconv;
  PetscBool      compute_errors = PETSC_FALSE;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;

	ierr = PetscOptionsGetBool(NULL,NULL,"-slepc_compute_errors",&compute_errors,NULL);CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"===== EigenAnalysis: \"%s\" ===== \n",description);
  ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);
	
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps,prob_type);CHKERRQ(ierr);
	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
	ierr = EPSSolve(eps);CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"===== EigenAnalysisReport: \"%s\" ===== \n",description);
	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
	ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);
	
	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);
	
	if (nconv>0) {
		/* Display eigenvalues and relative errors */
		if (compute_errors) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,
												 "               k                   ||Ax-kx||/||kx||\n"
												 "   --------------------------    --------------------\n");CHKERRQ(ierr);
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,
												 "               k              \n"
												 "   -------------------------- \n");CHKERRQ(ierr);
		}
		
		for (i=0;i<nconv;i++) {
			/* 
			 Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and ki (imaginary part)
			 */
			ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
			re = kr;
			im = ki;

      /*
       Compute the relative error associated to each eigenpair
       */
      if (compute_errors) {
        ierr = EPSComputeRelativeError(eps,i,&error);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," %1.6e %+1.6e j       %1.6e\n",re,im,error);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %1.6e %+1.6e j \n",re,im);CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"## GNUPLOT ## \n");CHKERRQ(ierr);
    for (i=0;i<nconv;i++) {
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
      re = kr;
      im = ki;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"%1.6e %1.6e\n",re,im);CHKERRQ(ierr);
    }
    
  }
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = VecDestroy(&xr);CHKERRQ(ierr);
  ierr = VecDestroy(&xi);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ptatinEigenAnalyser_Stokes"
PetscErrorCode ptatinEigenAnalyser_Stokes(SNES snes,PetscBool view)
{
  KSP ksp_stokes;
  PC pc_stokes;
  Mat As,Bs;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;

  ierr = SNESGetKSP(snes,&ksp_stokes);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_stokes,&As,&Bs);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_stokes,&pc_stokes);CHKERRQ(ierr);
  
  ierr = _slepc_eigen(As,EPS_NHEP,"StokesA");CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ptatinEigenAnalyser_StokesPC"
PetscErrorCode ptatinEigenAnalyser_StokesPC(SNES snes,PetscBool view)
{
  KSP ksp_stokes;
  PC pc_stokes;
  Mat As,Bs,Ae;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;

  ierr = SNESGetKSP(snes,&ksp_stokes);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_stokes,&As,&Bs);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_stokes,&pc_stokes);CHKERRQ(ierr);

  ierr = MatCreateEigenOperatorFromKSPOperators(ksp_stokes,&Ae);CHKERRQ(ierr);
  if (view) {
    PetscPrintf(PETSC_COMM_WORLD,"ptatinEigenAnalyser(StokesAPC): Using operators\n");
    ierr = KSPView(ksp_stokes,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  ierr = _slepc_eigen(Ae,EPS_NHEP,"StokesAPC");CHKERRQ(ierr);
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ptatinEigenAnalyser_A11PC"
PetscErrorCode ptatinEigenAnalyser_A11PC(SNES snes,PetscBool view)
{
  KSP ksp_stokes,ksp_A11,*sub_ksp;
  PC pc_stokes;
  Mat As,Bs,Ae;
  PetscInt nsplits;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = SNESGetKSP(snes,&ksp_stokes);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_stokes,&As,&Bs);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_stokes,&pc_stokes);CHKERRQ(ierr);

  ierr = PCFieldSplitGetSubKSP(pc_stokes,&nsplits,&sub_ksp);CHKERRQ(ierr);
  ksp_A11 = sub_ksp[0];
  
  ierr = MatCreateEigenOperatorFromKSPOperators(ksp_A11,&Ae);CHKERRQ(ierr);
  if (view) {
    PetscPrintf(PETSC_COMM_WORLD,"ptatinEigenAnalyser(StokesA11PC): Using operators\n");
    ierr = KSPView(ksp_A11,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }   
  
  ierr = _slepc_eigen(Ae,EPS_HEP,"StokesA11PC");CHKERRQ(ierr);
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ptatinEigenAnalyser_A11PCToMatlab"
PetscErrorCode ptatinEigenAnalyser_A11PCToMatlab(SNES snes,PetscBool view)
{
  KSP ksp_stokes,ksp_A11,*sub_ksp;
  PC pc_stokes,pc_A11;
  Mat As,Bs,Ae;
  PetscInt nsplits;
  PetscLogDouble t0,t1;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = SNESGetKSP(snes,&ksp_stokes);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_stokes,&As,&Bs);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_stokes,&pc_stokes);CHKERRQ(ierr);
  
  ierr = PCFieldSplitGetSubKSP(pc_stokes,&nsplits,&sub_ksp);CHKERRQ(ierr);
  ksp_A11 = sub_ksp[0];
  ierr = KSPGetPC(ksp_A11,&pc_A11);CHKERRQ(ierr);
  
  if (view) {
    PetscPrintf(PETSC_COMM_WORLD,"ptatinEigenAnalyser(StokesA11PC): Using operators\n");
    ierr = KSPView(ksp_A11,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }   
  
  PetscPrintf(PETSC_COMM_WORLD,"ptatinEigenAnalyser(A11PCToMatlab): Assembling explicit operator... this may take a while\n");
  
  PetscTime(&t0);
  ierr = PCComputeExplicitOperator(pc_A11,&Ae);CHKERRQ(ierr);
  PetscTime(&t1);
  PetscPrintf(PETSC_COMM_WORLD,"PCComputeExplicitOperator(time) %1.4e (sec)\n",t1-t0);
  
  ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
  ierr = MatView(Ae,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  {
    PC mm;
    Mat A11,Aii;
    
    ierr = KSPGetOperators(ksp_A11,&A11,0);CHKERRQ(ierr);
    ierr = PCCreate(PETSC_COMM_WORLD,&mm);CHKERRQ(ierr);
    ierr = PCSetOperators(mm,A11,A11);CHKERRQ(ierr);
    ierr = PCSetType(mm,PCMAT);CHKERRQ(ierr);

    PetscTime(&t0);
    ierr = PCComputeExplicitOperator(mm,&Aii);CHKERRQ(ierr);
    PetscTime(&t1);
    
    PetscPrintf(PETSC_COMM_WORLD,"A11PCComputeExplicitOperator(time) %1.4e (sec)\n",t1-t0);
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
    ierr = MatView(Aii,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    
    ierr = MatDestroy(&Aii);CHKERRQ(ierr);
}
  
  
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ptatinEigenAnalyser_A11PCSmoother"
PetscErrorCode ptatinEigenAnalyser_A11PCSmoother(SNES snes,PetscBool view)
{
  KSP ksp_stokes,ksp_A11,*sub_ksp,ksp_level;
  PC pc_stokes,pc_A11;
  Mat As,Bs,Ae;
  PetscInt nsplits,k,nlevels;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = SNESGetKSP(snes,&ksp_stokes);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_stokes,&As,&Bs);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_stokes,&pc_stokes);CHKERRQ(ierr);
  
  ierr = PCFieldSplitGetSubKSP(pc_stokes,&nsplits,&sub_ksp);CHKERRQ(ierr);
  ksp_A11 = sub_ksp[0];
  ierr = KSPGetPC(ksp_A11,&pc_A11);CHKERRQ(ierr);

  ierr = PCMGGetLevels(pc_A11,&nlevels);CHKERRQ(ierr);
  
  ierr = PCMGGetCoarseSolve(pc_A11,&ksp_level);CHKERRQ(ierr);
  //ierr = KSPGetOperators(ksp_level,&A_lv,&B_lv);CHKERRQ(ierr);

  ierr = MatCreateEigenOperatorFromKSPOperators(ksp_level,&Ae);CHKERRQ(ierr);
  if (view) {
    PetscPrintf(PETSC_COMM_WORLD,"ptatinEigenAnalyser(StokesA11_Lv[0]): Using operators\n");
    ierr = KSPView(ksp_level,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  
  ierr = _slepc_eigen(Ae,EPS_HEP,"StokesPCA11_Lv[0]");CHKERRQ(ierr);
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  
  for( k=1; k<nlevels; k++ ){
    char name[256];
    
    ierr = PCMGGetSmoother(pc_A11,k,&ksp_level);CHKERRQ(ierr);
  
    sprintf(name,"StokesPCA11_Lv[%d]",k);

    ierr = MatCreateEigenOperatorFromKSPOperators(ksp_level,&Ae);CHKERRQ(ierr);
    if (view) {
      PetscPrintf(PETSC_COMM_WORLD,"ptatinEigenAnalyser(%s): Using operators\n",name);
      ierr = KSPView(ksp_level,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    
    ierr = _slepc_eigen(Ae,EPS_HEP,name);CHKERRQ(ierr);
    ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ptatinEigenAnalyser_A11SmootherComputeExplicitOperator"
PetscErrorCode ptatinEigenAnalyser_A11SmootherComputeExplicitOperator(SNES snes,PetscBool view)
{
  KSP            ksp_stokes,ksp_A11,*sub_ksp,ksp_level;
  PC             pc_stokes,pc_A11;
  Mat            As,Bs,Ae,Ai;
  PetscInt       nsplits,k,nlevels;
  PetscViewer    viewer;
  PetscBool      ascii_view = PETSC_FALSE;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = SNESGetKSP(snes,&ksp_stokes);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_stokes,&As,&Bs);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_stokes,&pc_stokes);CHKERRQ(ierr);
  
  ierr = PCFieldSplitGetSubKSP(pc_stokes,&nsplits,&sub_ksp);CHKERRQ(ierr);
  ksp_A11 = sub_ksp[0];
  ierr = KSPGetPC(ksp_A11,&pc_A11);CHKERRQ(ierr);
  
  /* matlab agrees */
  /*
  {
    ierr = KSPGetOperators(ksp_A11,&Ae,NULL);
    ierr = _slepc_eigen(Ae,EPS_HEP,"StokesA11");CHKERRQ(ierr);
  }
  */
  ierr = PCMGGetLevels(pc_A11,&nlevels);CHKERRQ(ierr);
  
  ierr = PCMGGetCoarseSolve(pc_A11,&ksp_level);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_level,&Ai,0);CHKERRQ(ierr);
  ierr = MatComputeExplicitOperator(Ai,&Ae);CHKERRQ(ierr);

  if (ascii_view) {
    ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)Ae),"A_coarse.mat",&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)Ae),"A_coarse.mat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  }
  ierr = MatView(Ae,viewer);CHKERRQ(ierr);
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  if (ascii_view) { ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr); }
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  for( k=1; k<nlevels; k++ ){
    char     name[256];
    
    sprintf(name,"A_level_%d.mat",k);

    ierr = PCMGGetSmoother(pc_A11,k,&ksp_level);CHKERRQ(ierr);
    ierr = KSPGetOperators(ksp_level,&Ai,0);CHKERRQ(ierr);
    ierr = MatComputeExplicitOperator(Ai,&Ae);CHKERRQ(ierr);
 
    if (ascii_view) {
      ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)Ae),name,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)Ae),name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    }
    ierr = MatView(Ae,viewer);CHKERRQ(ierr);
    ierr = MatDestroy(&Ae);CHKERRQ(ierr);
    if (ascii_view) { ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr); }
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ptatinEigenAnalyser_A11PCSmootherComputeExplicitOperator"
PetscErrorCode ptatinEigenAnalyser_A11PCSmootherComputeExplicitOperator(SNES snes,PetscBool view)
{
  KSP ksp_stokes,ksp_A11,*sub_ksp,ksp_level;
  PC             pc_stokes,pc_A11;
  Mat            As,Bs,Ae;
  PetscInt       nsplits,k,nlevels;
  PetscViewer    viewer;
  PetscBool      ascii_view = PETSC_FALSE;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = SNESGetKSP(snes,&ksp_stokes);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_stokes,&As,&Bs);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_stokes,&pc_stokes);CHKERRQ(ierr);
  
  ierr = PCFieldSplitGetSubKSP(pc_stokes,&nsplits,&sub_ksp);CHKERRQ(ierr);
  ksp_A11 = sub_ksp[0];
  ierr = KSPGetPC(ksp_A11,&pc_A11);CHKERRQ(ierr);
  
  ierr = PCMGGetLevels(pc_A11,&nlevels);CHKERRQ(ierr);
  
  ierr = PCMGGetCoarseSolve(pc_A11,&ksp_level);CHKERRQ(ierr);
  ierr = KSPComputeExplicitOperator(ksp_level,&Ae);CHKERRQ(ierr);

  if (ascii_view) {
    ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)Ae),"A_ksp_coarse.mat",&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)Ae),"A_ksp_coarse.mat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  }
  
  ierr = MatView(Ae,viewer);CHKERRQ(ierr);
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  if (ascii_view) { ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr); }
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  for( k=1; k<nlevels; k++ ){
    char name[256];
    
    sprintf(name,"A_ksp_level_%d.mat",k);
    
    ierr = PCMGGetSmoother(pc_A11,k,&ksp_level);CHKERRQ(ierr);
    ierr = KSPComputeExplicitOperator(ksp_level,&Ae);CHKERRQ(ierr);
    if (ascii_view) {
      ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)Ae),name,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)Ae),name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    }
    ierr = MatView(Ae,viewer);CHKERRQ(ierr);
    ierr = MatDestroy(&Ae);CHKERRQ(ierr);
    if (ascii_view) { ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr); }
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* complete fine grid */
  ierr = KSPComputeExplicitOperator(ksp_A11,&Ae);CHKERRQ(ierr);
  if (ascii_view) {
    ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)Ae),"A_ksp_mg.mat",&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)Ae),"A_ksp_mg.mat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  }
  ierr = MatView(Ae,viewer);CHKERRQ(ierr);

  /* matlab agrees */
  /*
  {
    ierr = _slepc_eigen(Ae,EPS_NHEP,"StokesA.inv(A')");CHKERRQ(ierr);
  }
  */
   
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  if (ascii_view) { ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr); }
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = PCComputeExplicitOperator(pc_A11,&Ae);CHKERRQ(ierr);
  if (ascii_view) {
    ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)Ae),"A_pc_mg.mat",&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerBinaryOpen(PetscObjectComm((PetscObject)Ae),"A_pc_mg.mat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  }
  ierr = MatView(Ae,viewer);CHKERRQ(ierr);
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  if (ascii_view) { ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr); }
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ptatinEigenAnalyser_A11KSPSmoother"
PetscErrorCode ptatinEigenAnalyser_A11KSPSmoother(SNES snes,PetscBool view)
{
  KSP ksp_stokes,ksp_A11,*sub_ksp,ksp_level;
  PC pc_stokes,pc_A11;
  Mat As,Bs,Ae;
  PetscInt nsplits,k,nlevels;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = SNESGetKSP(snes,&ksp_stokes);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp_stokes,&As,&Bs);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_stokes,&pc_stokes);CHKERRQ(ierr);
  
  ierr = PCFieldSplitGetSubKSP(pc_stokes,&nsplits,&sub_ksp);CHKERRQ(ierr);
  ksp_A11 = sub_ksp[0];
  ierr = KSPGetPC(ksp_A11,&pc_A11);CHKERRQ(ierr);
  
  ierr = PCMGGetLevels(pc_A11,&nlevels);CHKERRQ(ierr);
  
  ierr = PCMGGetCoarseSolve(pc_A11,&ksp_level);CHKERRQ(ierr);
  
  ierr = MatCreateEigenOperatorFromKSP(ksp_level,&Ae);CHKERRQ(ierr);
  if (view) {
    PetscPrintf(PETSC_COMM_WORLD,"ptatinEigenAnalyser(StokesA11_Lv[0]): Using operators\n");
    ierr = KSPView(ksp_level,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  
  ierr = _slepc_eigen(Ae,EPS_HEP,"StokesKSPA11_Lv[0]");CHKERRQ(ierr);
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  
  for( k=1; k<nlevels; k++ ){
    char name[256];
    
    ierr = PCMGGetSmoother(pc_A11,k,&ksp_level);CHKERRQ(ierr);
    sprintf(name,"StokesKSPA11_Lv[%d]",k);
    
    ierr = MatCreateEigenOperatorFromKSP(ksp_level,&Ae);CHKERRQ(ierr);
    if (view) {
      PetscPrintf(PETSC_COMM_WORLD,"ptatinEigenAnalyser(%s): Using operators\n",name);
      ierr = KSPView(ksp_level,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    
    ierr = _slepc_eigen(Ae,EPS_HEP,name);CHKERRQ(ierr);
    ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}




#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_linear_viscous_forward_model_driver"
PetscErrorCode pTatin3d_linear_viscous_forward_model_driver(int argc,char **argv)
{
  DM             multipys_pack,dav,dap;
  pTatinCtx      user;
  Mat            A,B;
  Vec            X,F;
  IS             *is_stokes_field;
  SNES           snes;
  KSP            ksp;
  PC             pc;
  PetscInt       nlevels,k;
  Mat            operatorA11[10],operatorB11[10];
  DM             dav_hierarchy[10];
  Mat            interpolation_v[10],interpolation_eta[10];
  Quadrature     volQ[10];
  BCList         u_bclist[10];
  pTatinModel    model;
  PetscErrorCode ierr;
  PetscBool      asm_dump=PETSC_FALSE,block_dump=PETSC_FALSE;
  
  PetscFunctionBegin;
  
  /* PS: collect new flags for custom dumps */
  ierr = PetscOptionsGetBool(NULL,NULL,"-asm_dump",&asm_dump,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-block_dump",&block_dump,NULL);CHKERRQ(ierr);
  
  ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
  ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
  ierr = pTatinLogNote(user,"[ptatin_driver_linear_ts] -> new simulation");CHKERRQ(ierr);
  
  /* Register all models */
  ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
  /* Load model, call an initialization routines */
  ierr = pTatinModelLoad(user);CHKERRQ(ierr);
  ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);
  ierr = pTatinLogNote2(user,"[ptatin_driver_linear_ts] -> model loaded:",model->model_name);CHKERRQ(ierr);
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
  
  ierr = pTatinLogBasicDMDA(user,"Velocity",dav);CHKERRQ(ierr);
  ierr = pTatinLogBasicDMDA(user,"Pressure",dap);CHKERRQ(ierr);
  
  /* IF I DON'T DO THIS, THE IS's OBTAINED FROM DMCompositeGetGlobalISs() are wrong !! */
  {
    ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
  }
  ierr = DMCompositeGetGlobalISs(multipys_pack,&is_stokes_field);CHKERRQ(ierr); 
  
  ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
  
  /* mesh geometry */
  ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
  
  /* interpolate point coordinates (needed if mesh was modified) */
  //for (e=0; e<QUAD_EDGES; e++) {
  //  ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
  //}
  
  /* interpolate material point coordinates (needed if mesh was modified) */
  ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
  
  /* material geometry */
  ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
  
  /* boundary conditions */
  ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);
  
  
  /* setup mg */
  nlevels = 1;
  PetscOptionsGetInt(NULL,NULL,"-dau_nlevels",&nlevels,0);
  PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%d x %d x %d) : MG levels %d  \n", user->mx,user->my,user->mz,nlevels );
  ierr = pTatin3dStokesBuildMeshHierarchy(dav,nlevels,dav_hierarchy);CHKERRQ(ierr);
  ierr = pTatin3dStokesReportMeshHierarchy(nlevels,dav_hierarchy);CHKERRQ(ierr);
  for (k=nlevels-1; k>=0; k--) {
    char name[128];
    sprintf(name,"vel_dmda_Lv%d",k);
    ierr = pTatinLogBasicDMDA(user,name,dav_hierarchy[k]);CHKERRQ(ierr);
  }
  
  /* Define interpolation operators for velocity space */
  interpolation_v[0] = NULL;
  for (k=0; k<nlevels-1; k++) {
    ierr = DMCreateInterpolation(dav_hierarchy[k],dav_hierarchy[k+1],&interpolation_v[k+1],NULL);CHKERRQ(ierr);
  }
  
  /* Define interpolation operators for scalar space */
  interpolation_eta[0] = NULL;
  for (k=1; k<nlevels; k++) {
    ierr = MatMAIJRedimension(interpolation_v[k],1,&interpolation_eta[k]);CHKERRQ(ierr);
  }
  
  /* define material properties on gauss points on coarse grids */
  for (k=0; k<nlevels-1; k++) {
    PetscInt ncells,lmx,lmy,lmz;
    PetscInt np_per_dim;
    
    np_per_dim = 3;
    ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
    ncells = lmx * lmy * lmz;
    ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&volQ[k]);CHKERRQ(ierr);
  }
  volQ[nlevels-1] = user->stokes_ctx->volQ;
  
  
  /* define bounary list on coarse grids */
  for (k=0; k<nlevels-1; k++) {
    ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);
  }
  u_bclist[nlevels-1] = user->stokes_ctx->u_bclist;
  
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
    
    ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(user->coefficient_projection_type,npoints,mp_std,mp_stokes,nlevels,interpolation_eta,dav_hierarchy,volQ);CHKERRQ(ierr);
  }
  
  
  /* define bc's for hiearchy */
  ierr = pTatinModel_ApplyBoundaryConditionMG(nlevels,u_bclist,dav_hierarchy,user->model,user);CHKERRQ(ierr);
  
  
  /* configure stokes opertors */
  ierr = pTatin3dCreateStokesOperators(user->stokes_ctx,is_stokes_field,
                                       nlevels,dav_hierarchy,interpolation_v,u_bclist,volQ,
                                       &A,operatorA11,&B,operatorB11);CHKERRQ(ierr);
  
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
  
  /* write the initial fields */
  ierr = pTatinModel_Output(user->model,user,X,"icbc");CHKERRQ(ierr);
  
  
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
  ierr = PCSetType(pc,PCFIELDSPLIT);CHKERRQ(ierr);
  ierr = PCFieldSplitSetIS(pc,"u",is_stokes_field[0]);CHKERRQ(ierr);
  ierr = PCFieldSplitSetIS(pc,"p",is_stokes_field[1]);CHKERRQ(ierr);
  
  /* configure uu split for galerkin multi-grid */
  ierr = pTatin3dStokesKSPConfigureFSGMG(ksp,nlevels,operatorA11,operatorB11,interpolation_v);CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"   [[ COMPUTING FLOW FIELD FOR STEP : %D ]]\n", 0 );
  ierr = pTatinLogBasic(user);CHKERRQ(ierr);

  {
    KSP ksp;
    
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1);CHKERRQ(ierr);
  }
  
  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
  ierr = pTatinLogBasicSNES(user,"Stokes",snes);CHKERRQ(ierr);
  
  /* dump */
  ierr = pTatinModel_Output(user->model,user,X,"step000000");CHKERRQ(ierr);
  
  /* === do eigen analysis === */
  {
    PetscInt eigen_anal = 0;
    PetscBool view = PETSC_FALSE;
    
    ierr = PetscOptionsGetInt(NULL,NULL,"-eigen_anal",&eigen_anal,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL,NULL,"-eigen_anal_view_operators",&view,NULL);CHKERRQ(ierr);

    switch (eigen_anal) {
      case 0:
        ierr = ptatinEigenAnalyser_Stokes(snes,view);CHKERRQ(ierr);
        break;
      case 1:
        ierr = ptatinEigenAnalyser_StokesPC(snes,view);CHKERRQ(ierr);
        break;
      case 2:
        ierr = ptatinEigenAnalyser_A11PC(snes,view);CHKERRQ(ierr);
        break;
      case 3:
        ierr = ptatinEigenAnalyser_A11PCSmoother(snes,view);CHKERRQ(ierr);
        break;
      case 4:
        ierr = ptatinEigenAnalyser_A11KSPSmoother(snes,view);CHKERRQ(ierr);
        break;
        
      case 5:
        ierr = ptatinEigenAnalyser_A11PCToMatlab(snes,view);CHKERRQ(ierr);
        break;
        
      case 6:
        ierr = ptatinEigenAnalyser_A11SmootherComputeExplicitOperator(snes,view);CHKERRQ(ierr);
        ierr = ptatinEigenAnalyser_A11PCSmootherComputeExplicitOperator(snes,view);CHKERRQ(ierr);
        break;
    }
  }
  /* ========================= */

  /* Dump ASM submatrices if requested */
  if(asm_dump){
    PC pc;
    KSP ksp;
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = ASMDump(pc);CHKERRQ(ierr);
  }

  /* Dump (assembled) A11 and A12 blocks if requested */
  if(block_dump){
    ierr = BlockDump(user);CHKERRQ(ierr);
  }
  
  /* tidy up */
  for (k=0; k<nlevels; k++) {
    ierr = MatDestroy(&operatorA11[k]);CHKERRQ(ierr);
    ierr = MatDestroy(&operatorB11[k]);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  
  /* Clean up */
  for (k=0; k<nlevels-1; k++) {
    ierr = BCListDestroy(&u_bclist[k]);CHKERRQ(ierr);
  }
  for (k=0; k<nlevels-1; k++) {
    ierr = QuadratureDestroy(&volQ[k]);CHKERRQ(ierr);
  }
  for (k=1; k<nlevels; k++) {
    ierr = MatDestroy(&interpolation_v[k]);CHKERRQ(ierr);
    ierr = MatDestroy(&interpolation_eta[k]);CHKERRQ(ierr);
  }
  for (k=0; k<nlevels; k++) {
    ierr = DMDestroy(&dav_hierarchy[k]);CHKERRQ(ierr);
  }
  
  ierr = ISDestroy(&is_stokes_field[0]);CHKERRQ(ierr);
  ierr = ISDestroy(&is_stokes_field[1]);CHKERRQ(ierr);
  ierr = PetscFree(is_stokes_field);CHKERRQ(ierr);
  
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
  
  ierr = SlepcInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);

  ierr = pTatin3d_linear_viscous_forward_model_driver(argc,argv);CHKERRQ(ierr);
  
  ierr = SlepcFinalize();CHKERRQ(ierr);
  return 0;
}

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
 **    filename:   pc_dmdarepart.c
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

#include <petsc/private/pcimpl.h> /*I "petscpc.h" I*/
#include <petscksp.h>             /*I "petscksp.h" I*/
#include <petscdm.h>
#include <petscdmda.h>
#include <sub_comm.h>
#include <dmda_duplicate.h>

PetscLogStage _PCDMDARepart_PCApplyStage = -1;
PetscLogStage _PCDMDARepart_PCSetUpMatrixStage = -1;

/*

   Vec x;    = [ a b c ] [ d e f ] [ g h ]          <input>
   Vec xred; = [ a b c d e ][ 0 ][ f g h ][ 0 ]     <after scatter>
   Vec xsub; = [ a b c d e ] [ f g h ]              <defined on subcomm>

*/

typedef struct {
  PetscInt        nsubcomm_factor;
  PetscMPIInt     nsubcomm_size;
  PetscMPISubComm subcomm;
  KSP             ksp;
  IS              isin;
  VecScatter      scatter;
  Vec             xred,xsub,ysub;
  Mat             permutation;
  Mat             Bsub;
  DM              dmrepart;
  PetscInt        Mp_re,Np_re,Pp_re;
  PetscInt        *range_i_re,*range_j_re,*range_k_re;
  PetscInt        *start_i_re,*start_j_re,*start_k_re;
  PetscBool       log;
} PC_DMDARepart;

PetscErrorCode _DMDARepartitionDetermineRankFromGlobalIJK(PetscInt i,PetscInt j,PetscInt k,
    PetscInt Mp,PetscInt Np,PetscInt Pp,
    PetscInt start_i[],PetscInt start_j[],PetscInt start_k[],
    PetscInt span_i[],PetscInt span_j[],PetscInt span_k[],
    PetscMPIInt *_pi,PetscMPIInt *_pj,PetscMPIInt *_pk,PetscMPIInt *rank_re)
{
  PetscInt pi,pj,pk,n;

  pi = pj = pk = -1;
  for (n=0; n<Mp; n++) {
    if ( (i >= start_i[n]) && (i < start_i[n]+span_i[n]) ) {
      pi = n;
      break;
    }
  }
  for (n=0; n<Np; n++) {
    if ( (j >= start_j[n]) && (j < start_j[n]+span_j[n]) ) {
      pj = n;
      break;
    }
  }
  for (n=0; n<Pp; n++) {
    if ( (k >= start_k[n]) && (k < start_k[n]+span_k[n]) ) {
      pk = n;
      break;
    }
  }
  if (pi == -1) { SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"  [dmdarepart][pi] cannot be determined : range %D, val %D",Mp,i); }
  if (pj == -1) { SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"  [dmdarepart][pj] cannot be determined : range %D, val %D",Np,j); }
  if (pk == -1) { SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"  [dmdarepart][pk] cannot be determined : range %D, val %D",Pp,k); }

  *_pi = pi;
  *_pj = pj;
  *_pk = pk;
  *rank_re = pi + pj * Mp + pk * (Mp*Np);

  PetscFunctionReturn(0);
}

PetscErrorCode _DMDARepartitionDetermineGlobalS0(PetscMPIInt rank_re,PetscInt Mp_re,PetscInt Np_re,PetscInt Pp_re,
    PetscInt range_i_re[],PetscInt range_j_re[],PetscInt range_k_re[],PetscInt *s0)
{
  PetscInt i,j,k,start_IJK;

  start_IJK = 0;
  for (k=0; k<Pp_re; k++) {
    for (j=0; j<Np_re; j++) {
      for (i=0; i<Mp_re; i++) {
        PetscInt rank_ijk;

        rank_ijk = i + j*Mp_re + k*Mp_re*Np_re;

        if (rank_ijk < rank_re) {
          start_IJK += range_i_re[i]*range_j_re[j]*range_k_re[k];
        }

      }
    }
  }
  *s0 = start_IJK;
  PetscFunctionReturn(0);
}

/*
   Vec x;    = [ a b c ] [ d e f ] [ g h ]          <input>
   Vec xred; = [ a b c d e ][ 0 ][ f g h ][ 0 ]     <after scatter>
   Vec xsub; = [ a b c d e ] [ f g h ]              <defined on subcomm>
*/
PetscErrorCode _DMDARepart_SetupScatters(PC pc,PC_DMDARepart *red)
{
  PetscErrorCode ierr;
  PetscBool      active;
  Vec            xsub,xred,x;
  MPI_Comm       comm;
  PetscInt       st,ed,n,N,m;
  IS             isin;
  VecScatter     scatter;
  Mat            B;

  ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
  ierr = PCGetOperators(pc,NULL,&B);CHKERRQ(ierr);
  ierr = MatCreateVecs(B,&x,NULL);CHKERRQ(ierr);

  xsub = NULL;
  ierr = PetscMPISubCommGetActive(red->subcomm,&active);CHKERRQ(ierr);
  if (active) {
    ierr = DMCreateGlobalVector(red->dmrepart,&xsub);CHKERRQ(ierr);
  }
  //ierr = _VecCreateFromSubcomm(xsub,comm,&xred);CHKERRQ(ierr);
  //VecView(xred,PETSC_VIEWER_STDOUT_(comm));CHKERRQ(ierr);

  m = -1;
  if (xsub) {
    ierr = VecGetLocalSize(xsub,&m);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(xsub,&st,&ed);CHKERRQ(ierr);
    ierr = ISCreateStride(comm,ed-st,st,1,&isin);CHKERRQ(ierr);
  } else {
    /* fetch some local owned data - just to deal with avoiding zero length ownership on range */
    ierr = VecGetOwnershipRange(x,&st,&ed);CHKERRQ(ierr);
    ierr = ISCreateStride(comm,1,st,1,&isin);CHKERRQ(ierr);
  }
  if (active) {
    ierr = VecDestroy(&xsub);CHKERRQ(ierr);
  }

  ierr = ISGetLocalSize(isin,&n);CHKERRQ(ierr);
  ierr = ISGetSize(isin,&N);CHKERRQ(ierr);
  ierr = VecCreate(PetscObjectComm((PetscObject)isin),&xred);CHKERRQ(ierr);
  ierr = VecSetSizes(xred,n,N);CHKERRQ(ierr);
  ierr = VecSetType(xred,((PetscObject)x)->type_name);CHKERRQ(ierr);

  /* VecView(xred,PETSC_VIEWER_STDOUT_(comm));CHKERRQ(ierr); */
  ierr = VecScatterCreate(x,isin,xred,NULL,&scatter);CHKERRQ(ierr);

  xsub = NULL;
  if (active) {
    ierr = VecCreateMPIWithArray(red->subcomm->sub_comm,1,m,PETSC_DECIDE,NULL,&xsub);CHKERRQ(ierr);
  }

  if (active) {
    ierr = DMCreateGlobalVector(red->dmrepart,&red->ysub);CHKERRQ(ierr);
  }
  red->xred    = xred;
  red->isin    = isin;
  red->scatter = scatter;
  red->xsub    = xsub;

  /* testing */
  /*
     ierr = VecGetOwnershipRange(x,&st,&ed);CHKERRQ(ierr);
     for (i=st; i<ed; i++) {
     VecSetValue(x,i,(PetscScalar)i,INSERT_VALUES);
     }
     VecAssemblyBegin(x);
     VecAssemblyEnd(x);
     VecView(x,PETSC_VIEWER_STDOUT_(comm));CHKERRQ(ierr);

     ierr = VecScatterBegin(scatter,x,xred,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
     ierr = VecScatterEnd(scatter,x,xred,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

     if (active) {
     PetscScalar *array;

     ierr = VecCreateMPIWithArray(red->subcomm->sub_comm,1,m,PETSC_DECIDE,NULL,&xsub);CHKERRQ(ierr);

  // define xred //
  ierr = VecGetArray(xred,&array);CHKERRQ(ierr);
  // we created xred with empty local arrays, now we fill it in //
  ierr = VecPlaceArray(xsub,(const PetscScalar*)array);CHKERRQ(ierr);

  VecView(xsub,PETSC_VIEWER_STDOUT_(red->subcomm->sub_comm));

  ierr = VecResetArray(xsub);CHKERRQ(ierr);
  ierr = VecRestoreArray(xred,&array);CHKERRQ(ierr);
  }
  */

  ierr = VecDestroy(&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode _DMDARepart_SetupPMatrix(PC pc,PC_DMDARepart *red)
{
  PetscErrorCode ierr;
  DM             dm,dmsc;
  MPI_Comm       comm;
  Mat            Pscalar,P;
  PetscInt       ndof;
  PetscInt       i,j,k,location,startI[3],endI[3],lenI[3],nx,ny,nz;
  PetscInt       sr,er,Mr;
  PetscMPIInt    rank;
  Vec            V;

  ierr = PCGetDM(pc,&dm);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  PetscInfo(pc,"Setting up the permutation matrix\n");

  ierr = DMDADuplicateLayout(dm,1,0,DMDA_STENCIL_BOX,&dmsc);CHKERRQ(ierr); /* stencil type (box/star) ignored as stencil width = 0 */

  ierr = DMCreateGlobalVector(dmsc,&V);CHKERRQ(ierr);
  ierr = VecGetSize(V,&Mr);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(V,&sr,&er);CHKERRQ(ierr);
  ierr = VecDestroy(&V);CHKERRQ(ierr);

  ierr = MatCreate(comm,&Pscalar);CHKERRQ(ierr);
  ierr = MatSetSizes(Pscalar,(er-sr),(er-sr),Mr,Mr);CHKERRQ(ierr);
  ierr = MatSetType(Pscalar,MATAIJ);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(Pscalar,2,NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(Pscalar,2,NULL,2,NULL);CHKERRQ(ierr);

  ierr = DMDAGetCorners(dm,NULL,NULL,NULL,&lenI[0],&lenI[1],&lenI[2]);CHKERRQ(ierr);
  ierr = DMDAGetCorners(dm,&startI[0],&startI[1],&startI[2],&endI[0],&endI[1],&endI[2]);CHKERRQ(ierr);
  endI[0] += startI[0];
  endI[1] += startI[1];
  endI[2] += startI[2];

  ierr = DMDAGetInfo(dm,NULL, &nx,&ny,&nz, NULL,NULL,NULL, &ndof,NULL,NULL,NULL,NULL, NULL);CHKERRQ(ierr);

  for (k=startI[2]; k<endI[2]; k++) {
    for (j=startI[1]; j<endI[1]; j++) {
      for (i=startI[0]; i<endI[0]; i++) {
        PetscMPIInt rank_ijk_re,rank_reI[3];
        PetscInt    s0_re;
        PetscInt    ii,jj,kk,local_ijk_re,mapped_ijk;
        /* PetscInt    natural_ijk; */
        PetscInt    lenI_re[3];

        location = (i - startI[0]) + (j - startI[1])*lenI[0] + (k - startI[2])*lenI[0]*lenI[1];

        ierr = _DMDARepartitionDetermineRankFromGlobalIJK(i,j,k,  red->Mp_re,red->Np_re,red->Pp_re,
                                                          red->start_i_re,red->start_j_re,red->start_k_re,
                                                          red->range_i_re,red->range_j_re,red->range_k_re,
                                                          &rank_reI[0],&rank_reI[1],&rank_reI[2],&rank_ijk_re);

        ierr = _DMDARepartitionDetermineGlobalS0(rank_ijk_re, red->Mp_re,red->Np_re,red->Pp_re, red->range_i_re,red->range_j_re,red->range_k_re, &s0_re);CHKERRQ(ierr);

        ii = i - red->start_i_re[ rank_reI[0] ];
        if (ii < 0) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"  [dmdarepart] index error ii \n"); }

        jj = j - red->start_j_re[ rank_reI[1] ];
        if (jj < 0) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"  [dmdarepart] index error jj \n"); }

        kk = k - red->start_k_re[ rank_reI[2] ];
        if (kk < 0) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"  [dmdarepart] index error kk \n"); }

        lenI_re[0] = red->range_i_re[ rank_reI[0] ];
        lenI_re[1] = red->range_j_re[ rank_reI[1] ];
        lenI_re[2] = red->range_k_re[ rank_reI[2] ];

        local_ijk_re = ii + jj * lenI_re[0] + kk * lenI_re[0] * lenI_re[1];
        mapped_ijk = s0_re + local_ijk_re;
        /* natural_ijk = i + j*nx + k*nx*ny; */

        ierr = MatSetValue(Pscalar,sr+location,mapped_ijk,1.0,INSERT_VALUES);CHKERRQ(ierr);
        //PetscPrintf(PETSC_COMM_SELF,"[%d] (%D,%D,%D) --> local %D [g=%D] [natural %D] --> repart %D \n",(int)rank,i,j,k,location,sr+location,natural_ijk,mapped_ijk);
      }
    }
  }
  ierr = MatAssemblyBegin(Pscalar,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Pscalar,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatCreateMAIJ(Pscalar,ndof,&P);CHKERRQ(ierr);
  ierr = MatDestroy(&Pscalar);CHKERRQ(ierr);

  red->permutation = P;

  ierr = DMDestroy(&dmsc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode _DMDARepart_UpdateOperator(PC pc,PC_DMDARepart *red)
{
  PetscErrorCode ierr;
  DM             dm;
  PetscBool      active,new_Bsub;
  MPI_Comm       comm;
  Mat            P,B,Bperm;
  PetscInt       ndof,i,j,sr,er,Mc,m_repart;
  PetscMPIInt    rank;
  IS             isrow,iscol;
  Mat            Blocal,*_Blocal,Bsub;

  ierr = PCGetDM(pc,&dm);CHKERRQ(ierr);
  ierr = PCGetOperators(pc,NULL,&B);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)B,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  PetscInfo(pc,"Updating the redundant preconditioned operator\n");

  P = red->permutation;
  ierr = MatPtAP(B,P,MAT_INITIAL_MATRIX,1.1,&Bperm);CHKERRQ(ierr);

  /* Get submatrices */
  ierr = MatGetSize(B,NULL,&Mc);CHKERRQ(ierr);
  ierr = ISCreateStride(comm,Mc,0,1,&iscol);CHKERRQ(ierr);

  ierr = PetscMPISubCommGetActive(red->subcomm,&active);CHKERRQ(ierr);

  if (active) {
    ierr = VecGetOwnershipRange(red->xsub,&sr,&er);CHKERRQ(ierr);
    ierr = ISCreateStride(comm,(er-sr),sr,1,&isrow);CHKERRQ(ierr);
  } else {
    PetscInt start,end;

    /* fetch some local owned data - just to deal with avoiding zero length ownership on range */
    ierr = MatGetOwnershipRange(B,&start,&end);CHKERRQ(ierr);
    ierr = ISCreateStride(comm,1,start,1,&isrow);CHKERRQ(ierr);
  }

  ierr = MatCreateSubMatrices(Bperm,1,&isrow,&iscol,MAT_INITIAL_MATRIX,&_Blocal);CHKERRQ(ierr);
  Blocal = *_Blocal;
  ierr = ISDestroy(&isrow);CHKERRQ(ierr);
  ierr = ISDestroy(&iscol);CHKERRQ(ierr);

  m_repart = 0;
  if (active) {
    ierr = VecGetLocalSize(red->xsub,&m_repart);CHKERRQ(ierr);
  }

  new_Bsub = PETSC_FALSE;
  if (active) {
    /* if it exists, fetch it */
    if (red->Bsub) {
      PetscInfo(pc,"Fetching existing redundant matrix\n");
      Bsub = red->Bsub;
      ierr = MatZeroEntries(Bsub);CHKERRQ(ierr);
    } else {
      PetscInfo(pc,"Creating new redundant matrix\n");
      ierr = DMDAGetInfo(dm,NULL, NULL,NULL,NULL, NULL,NULL,NULL, &ndof,NULL,NULL,NULL,NULL, NULL);CHKERRQ(ierr);

      ierr = MatCreate(red->subcomm->sub_comm,&Bsub);CHKERRQ(ierr);
      ierr = MatSetSizes(Bsub,m_repart,m_repart,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
      ierr = MatSetBlockSize(Bsub,ndof);CHKERRQ(ierr);
      ierr = MatSetFromOptions(Bsub);CHKERRQ(ierr);
      ierr = MatSetUp(Bsub);CHKERRQ(ierr);
      /* if I just created it, set the pointer */
      red->Bsub = Bsub;
      new_Bsub = PETSC_TRUE;
    }
  } else {
    Bsub = NULL;
  }

  /* insert entries */
  if (active) {
    PetscInt          start,end;
    PetscInt          ncols,startc,endc,rowidx;
    const PetscInt    *cols;
    const PetscScalar *vals;

    if (new_Bsub) {
      PetscInt *nnz,*onnz;

      PetscInfo(pc,"Defining preallocation for redundant matrix\n");
      /* preallocation */
      ierr = VecGetOwnershipRange(red->xsub,&start,&end);CHKERRQ(ierr);
      ierr = VecGetOwnershipRange(red->xsub,&startc,&endc);CHKERRQ(ierr);

      PetscMalloc(sizeof(PetscInt)*(end-start),&nnz);
      PetscMalloc(sizeof(PetscInt)*(end-start),&onnz);
      PetscMemzero(nnz,sizeof(PetscInt)*(end-start));
      PetscMemzero(onnz,sizeof(PetscInt)*(end-start));

      for (i=0; i<(end-start); i++) {
        ierr = MatGetRow(Blocal,i,&ncols,&cols,NULL);CHKERRQ(ierr);
        for (j=0; j<ncols; j++) {
          if ( (cols[j] >= startc) && (cols[j] < endc) ) {
            nnz[i]++;
          } else {
            onnz[i]++;
          }
        }
        ierr = MatRestoreRow(Blocal,i,&ncols,&cols,NULL);CHKERRQ(ierr);
      }
      ierr = MatSeqAIJSetPreallocation(Bsub,PETSC_DEFAULT,nnz);CHKERRQ(ierr);
      ierr = MatMPIAIJSetPreallocation(Bsub,PETSC_DEFAULT,nnz,PETSC_DEFAULT,onnz);CHKERRQ(ierr);

      PetscFree(nnz);
      PetscFree(onnz);
    }

    /* insert */
    ierr = MatGetOwnershipRange(Bsub,&start,&end);CHKERRQ(ierr);
    for (i=0; i<(end-start); i++) {
      ierr = MatGetRow(Blocal,i,&ncols,&cols,&vals);CHKERRQ(ierr);

      rowidx = i + start;
      ierr = MatSetValues(Bsub,1,&rowidx,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);

      ierr = MatRestoreRow(Blocal,i,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(Bsub,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Bsub,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  ierr = MatDestroy(&Bperm);CHKERRQ(ierr);
  ierr = MatDestroyMatrices(1,&_Blocal);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* implementations for DMDARepart */
static PetscErrorCode PCSetUp_DMDARepart(PC pc)
{
  PetscErrorCode ierr;
  PC_DMDARepart  *red = (PC_DMDARepart*)pc->data;

  PetscFunctionBegin;
  /* construction phase */
  if (!pc->setupcalled) {
    DM                    dm;
    PetscBool             isdmda;
    MPI_Comm              comm;
    PetscMPISubComm       subcomm;
    PetscMPIInt           rank,comm_size,nsubcomm_size;
    PetscInt              sum,k,nx,ny,nz,ndof,nsw,dim;
    const PetscInt        *_range_i_re;
    const PetscInt        *_range_j_re;
    const PetscInt        *_range_k_re;
    DMDAStencilType       stencil;
    DMBoundaryType        bx,by,bz;
    DMDAInterpolationType itype;
    PetscInt              refine_x,refine_y,refine_z;

    ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
    ierr = PCGetDM(pc,&dm);CHKERRQ(ierr);
    if (!dm) SETERRQ(comm,PETSC_ERR_SUP,"You must attach a DM to the PC via KSPSetDM, or PCSetDM");
    ierr = PetscObjectTypeCompare((PetscObject)dm,DMDA,&isdmda);CHKERRQ(ierr);
    if (!isdmda) SETERRQ(comm,PETSC_ERR_SUP,"DM attached to the PC must be of type DMDA");

    /* set up sub communicator */
    ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

    ierr = PetscMPISubCommCreate(comm,red->nsubcomm_factor,&subcomm);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&comm_size);CHKERRQ(ierr);
    ierr = MPI_Comm_size(subcomm->sub_comm,&nsubcomm_size);CHKERRQ(ierr);

    PetscInfo2(pc,"Creating PetscMPISubComm. Restricting from %d to %d ranks\n",(int)comm_size,(int)nsubcomm_size);
    red->nsubcomm_size = nsubcomm_size;
    red->subcomm       = subcomm;

    /* setup krylov solver */
    red->ksp = NULL;
    if (red->subcomm->parent_rank_active_in_subcomm) {
      const char *prefix;
      PetscInt   tablevel;

      ierr = PetscObjectGetTabLevel((PetscObject)pc,&tablevel);CHKERRQ(ierr);
      tablevel += 4;

      ierr = KSPCreate(red->subcomm->sub_comm,&red->ksp);CHKERRQ(ierr);
      ierr = PetscObjectIncrementTabLevel((PetscObject)red->ksp,(PetscObject)pc,tablevel);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)pc,(PetscObject)red->ksp);CHKERRQ(ierr);

      ierr = PCGetOptionsPrefix(pc,&prefix);CHKERRQ(ierr);
      ierr = KSPSetOptionsPrefix(red->ksp,prefix);CHKERRQ(ierr);
      ierr = KSPAppendOptionsPrefix(red->ksp,"dmdarepart_");CHKERRQ(ierr);

      ierr = KSPSetType(red->ksp,KSPPREONLY);CHKERRQ(ierr);
    }

    /* setup repartitioned dm */
    ierr = DMDAGetInfo(dm,&dim,&nx,&ny,&nz,0,0,0, &ndof,&nsw,&bx,&by,&bz,&stencil);CHKERRQ(ierr);
    if (dim != 3) SETERRQ(comm,PETSC_ERR_SUP,"DM attached to the PC must be a DMDA defined in three dimensions");
    ierr = DMDAGetInterpolationType(dm,&itype);CHKERRQ(ierr);
    ierr = DMDAGetRefinementFactor(dm,&refine_x,&refine_y,&refine_z);CHKERRQ(ierr);

    red->dmrepart = NULL;
    _range_i_re = _range_j_re = _range_k_re = NULL;
    if (red->subcomm->parent_rank_active_in_subcomm) {

      //ierr = DMDACreate3d(subcomm->sub_comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,stencil,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, ndof,nsw, NULL,NULL,NULL,&red->dmrepart);CHKERRQ(ierr);
      // ierr = DMSetUp(red->dmrepart);CHKERRQ(ierr);
      /* Note - I just use stencil_width = 1 here - this allows me to tunnel down deep with gmg without getting errors about stencil width > overlap */
      ierr = DMDACreate3d(subcomm->sub_comm,bx,by,bz,stencil,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, ndof,1, NULL,NULL,NULL,&red->dmrepart);CHKERRQ(ierr);
      ierr = DMSetUp(red->dmrepart);CHKERRQ(ierr);
      ierr = DMSetOptionsPrefix(red->dmrepart,"repart_");CHKERRQ(ierr);

      ierr = DMDASetRefinementFactor(red->dmrepart,refine_x,refine_y,refine_z);CHKERRQ(ierr);
      ierr = DMDASetInterpolationType(red->dmrepart,itype);CHKERRQ(ierr);

      //ierr = DMView(red->dmrepart, PETSC_VIEWER_STDOUT_(subcomm->sub_comm));CHKERRQ(ierr);

      ierr = DMDAGetInfo(red->dmrepart,0,0,0,0,&red->Mp_re,&red->Np_re,&red->Pp_re, 0,0,0,0,0,0);CHKERRQ(ierr);
      ierr = DMDAGetOwnershipRanges(red->dmrepart,&_range_i_re,&_range_j_re,&_range_k_re);CHKERRQ(ierr);
    }

    /* generate ranges for repartitioned dm */
    /* note - assume rank 0 always participates */
    ierr = MPI_Bcast(&red->Mp_re,1,MPIU_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(&red->Np_re,1,MPIU_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(&red->Pp_re,1,MPIU_INT,0,comm);CHKERRQ(ierr);

    PetscMalloc(sizeof(PetscInt)*red->Mp_re,&red->range_i_re);
    PetscMalloc(sizeof(PetscInt)*red->Np_re,&red->range_j_re);
    PetscMalloc(sizeof(PetscInt)*red->Pp_re,&red->range_k_re);

    if (_range_i_re != NULL) { PetscMemcpy(red->range_i_re,_range_i_re,sizeof(PetscInt)*red->Mp_re); }
    if (_range_j_re != NULL) { PetscMemcpy(red->range_j_re,_range_j_re,sizeof(PetscInt)*red->Np_re); }
    if (_range_k_re != NULL) { PetscMemcpy(red->range_k_re,_range_k_re,sizeof(PetscInt)*red->Pp_re); }

    ierr = MPI_Bcast(red->range_i_re,red->Mp_re,MPIU_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(red->range_j_re,red->Np_re,MPIU_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(red->range_k_re,red->Pp_re,MPIU_INT,0,comm);CHKERRQ(ierr);

    PetscMalloc(sizeof(PetscInt)*red->Mp_re,&red->start_i_re);
    PetscMalloc(sizeof(PetscInt)*red->Np_re,&red->start_j_re);
    PetscMalloc(sizeof(PetscInt)*red->Pp_re,&red->start_k_re);

    sum = 0;
    for (k=0; k<red->Mp_re; k++) {
      red->start_i_re[k] = sum;
      sum += red->range_i_re[k];
    }

    sum = 0;
    for (k=0; k<red->Np_re; k++) {
      red->start_j_re[k] = sum;
      sum += red->range_j_re[k];
    }

    sum = 0;
    for (k=0; k<red->Pp_re; k++) {
      red->start_k_re[k] = sum;
      sum += red->range_k_re[k];
    }

    /* view result */
    /*
       if (red->subcomm->parent_rank_active_in_subcomm) {
       PetscMPIInt rank_re;
       PetscInt    start_IJK;

       start_IJK = -1;
       ierr = MPI_Comm_rank(subcomm->sub_comm,&rank_re);CHKERRQ(ierr);
       ierr = _DMDARepartitionDetermineGlobalS0(rank_re,red->Mp_re,red->Np_re,red->Pp_re,red->range_i_re,red->range_j_re,red->range_k_re,&start_IJK);CHKERRQ(ierr);
       PetscPrintf(PETSC_COMM_SELF,"  [dmdarepart] rank[%d]: subrank[%d]: start idx = %D \n",(int)rank,(int)rank_re,start_IJK);
       } else {
       PetscPrintf(PETSC_COMM_SELF,"  [dmdarepart] rank[%d]: dmrepart doesn't live on this rank \n",(int)rank);
       }
    */
    /* Determine (i,j,k) value of subcomm ranks */
    /*
       if (red->subcomm->parent_rank_active_in_subcomm) {
       PetscMPIInt rank_re;
       PetscInt rankIJ_re;
       PetscInt pI_re,pJ_re,pK_re;

       pI_re = pJ_re = pK_re = -1;
       ierr = MPI_Comm_rank(subcomm->sub_comm,&rank_re);CHKERRQ(ierr);
       pK_re = rank_re/(red->Mp_re*red->Np_re);
       rankIJ_re = rank_re - pK_re * (red->Mp_re*red->Np_re);
       pJ_re = rankIJ_re/red->Mp_re;
       pI_re = rankIJ_re - pJ_re*red->Mp_re;

       PetscPrintf(PETSC_COMM_SELF,"  [dmdarepart] rank[%d]: subrank[%d] (%D %D %D)\n",(int)rank,(int)rank_re,pI_re,pJ_re,pK_re);
       }
    */
    /* attach dm to ksp on sub communicator */
    if (subcomm->parent_rank_active_in_subcomm) {
      ierr = KSPSetDM(red->ksp,red->dmrepart);CHKERRQ(ierr);
      ierr = KSPSetDMActive(red->ksp,PETSC_FALSE);CHKERRQ(ierr);
    }
  }

  /* setup scatters and permutation matrix */
  if (!pc->setupcalled) {
    ierr = _DMDARepart_SetupScatters(pc,red);CHKERRQ(ierr);
    ierr = _DMDARepart_SetupPMatrix(pc,red);CHKERRQ(ierr);
  }

  /* fetch redundant matrix - PCSetUp_DMDARepart is called we need to update the entires in the matrix */
  if (red->log) { ierr = PetscLogStagePush(_PCDMDARepart_PCSetUpMatrixStage);CHKERRQ(ierr); }
  ierr = _DMDARepart_UpdateOperator(pc,red);CHKERRQ(ierr);
  if (red->log) { ierr = PetscLogStagePop();CHKERRQ(ierr); }

  /* common - no construction */
  if (red->Bsub) {
    ierr = KSPSetOperators(red->ksp,red->Bsub,red->Bsub);CHKERRQ(ierr);
    if (pc->setfromoptionscalled && !pc->setupcalled){
      ierr = KSPSetFromOptions(red->ksp);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PCApply_DMDARepart(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;
  PC_DMDARepart  *red = (PC_DMDARepart*)pc->data;
  PetscBool      active;
  PetscScalar    *LA_red,*LA_sub;
  Vec            xtmp;

  PetscFunctionBegin;
  if (red->log) { ierr = PetscLogStagePush(_PCDMDARepart_PCApplyStage);CHKERRQ(ierr); }
  ierr = VecDuplicate(x,&xtmp);CHKERRQ(ierr);
  ierr = MatMultTranspose(red->permutation,x,xtmp);CHKERRQ(ierr);
  //ierr = MatMult(red->permutation,x,xtmp);CHKERRQ(ierr);

  /* pull in vector */
  ierr = VecScatterBegin(red->scatter,xtmp,red->xred,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(red->scatter,xtmp,red->xred,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

  ierr = PetscMPISubCommGetActive(red->subcomm,&active);CHKERRQ(ierr);
  if (active) {
    PetscScalar *array;

    /* define xsub */
    ierr = VecGetArray(red->xred,&array);CHKERRQ(ierr);
    /* we created xsub with empty local arrays, now we fill it in */
    ierr = VecPlaceArray(red->xsub,(const PetscScalar*)array);CHKERRQ(ierr);

    ierr = KSPSolve(red->ksp,red->xsub,red->ysub);CHKERRQ(ierr); /* Note: We could override rhs, however VecGetArray() with objects created via VecCreateWithArray() causes segv... */

    ierr = VecResetArray(red->xsub);CHKERRQ(ierr);
    ierr = VecRestoreArray(red->xred,&array);CHKERRQ(ierr);
  }

  /* return vector */
  ierr = VecGetArray(red->xred,&LA_red);CHKERRQ(ierr);
  if (active) {
    PetscInt i,start,end;

    ierr = VecGetOwnershipRange(red->ysub,&start,&end);CHKERRQ(ierr);
    ierr = VecGetArray(red->ysub,&LA_sub);CHKERRQ(ierr);
    for (i=0; i<end-start; i++) {
      LA_red[i] = LA_sub[i];
    }
    ierr = VecRestoreArray(red->ysub,&LA_sub);CHKERRQ(ierr);
  } else {
    /* fill in the dummy entry */
    LA_red[0] = 0.0; /* insert zero as we will ADD_VALUES in the scatter below */
  }
  ierr = VecRestoreArray(red->xred,&LA_red);CHKERRQ(ierr);

  /* scatter back from xred -> y */
  ierr = VecZeroEntries(xtmp);CHKERRQ(ierr);
  ierr = VecScatterBegin(red->scatter,red->xred,xtmp,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  ierr = VecScatterEnd(red->scatter,red->xred,xtmp,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);

  ierr = MatMult(red->permutation,xtmp,y);CHKERRQ(ierr);
  //ierr = MatMultTranspose(red->permutation,xtmp,y);CHKERRQ(ierr);

  ierr = VecDestroy(&xtmp);CHKERRQ(ierr);
  if (red->log) { ierr = PetscLogStagePop();CHKERRQ(ierr); }

  PetscFunctionReturn(0);
}

static PetscErrorCode PCReset_DMDARepart(PC pc)
{
  PetscErrorCode ierr;
  PC_DMDARepart  *red = (PC_DMDARepart*)pc->data;

  PetscFunctionBegin;
  if (red->xred) { ierr = VecDestroy(&red->xred);CHKERRQ(ierr);}
  if (red->xsub) { ierr = VecDestroy(&red->xsub);CHKERRQ(ierr);}
  if (red->ysub) { ierr = VecDestroy(&red->ysub);CHKERRQ(ierr);}
  if (red->Bsub) { ierr = MatDestroy(&red->Bsub);CHKERRQ(ierr);}
  if (red->permutation) { ierr = MatDestroy(&red->permutation);CHKERRQ(ierr);}

  if (red->dmrepart) { ierr = DMDestroy(&red->dmrepart);CHKERRQ(ierr); }

  ierr = ISDestroy(&red->isin);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&red->scatter);CHKERRQ(ierr);

  PetscFree(red->range_i_re);
  PetscFree(red->range_j_re);
  PetscFree(red->range_k_re);

  PetscFree(red->start_i_re);
  PetscFree(red->start_j_re);
  PetscFree(red->start_k_re);

  if (red->ksp) {  ierr = KSPReset(red->ksp);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

static PetscErrorCode PCDestroy_DMDARepart(PC pc)
{
  PetscErrorCode ierr;
  PC_DMDARepart  *red = (PC_DMDARepart*)pc->data;

  PetscFunctionBegin;
  ierr = PCReset_DMDARepart(pc);CHKERRQ(ierr);
  if (red->ksp) {
    ierr = KSPDestroy(&red->ksp);CHKERRQ(ierr);
  }
  ierr = PetscMPISubCommDestroy(&red->subcomm);CHKERRQ(ierr);
  ierr = PetscFree(pc->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PCSetFromOptions_DMDARepart(PetscOptionItems *PetscOptionsObject,PC pc)
{
  PetscErrorCode ierr;
  PC_DMDARepart  *red = (PC_DMDARepart*)pc->data;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject, "DMDARepart options");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-pc_dmdarepart_factor","Factor to reduce parent communication size by","PCDMDARepartSetFactor",red->nsubcomm_factor,&red->nsubcomm_factor,0);CHKERRQ(ierr);
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PCView_DMDARepart(PC pc,PetscViewer viewer)
{
  PetscErrorCode  ierr;
  PC_DMDARepart   *red = (PC_DMDARepart*)pc->data;
  PetscBool       iascii,isstring;
  PetscViewer     subviewer;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERSTRING,&isstring);CHKERRQ(ierr);
  if (iascii) {
    if (!red->subcomm) {
      ierr = PetscViewerASCIIPrintf(viewer,"  DMDARepart: preconditioner not yet setup\n");CHKERRQ(ierr);
    } else {
      /* ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant preconditioner:\n");CHKERRQ(ierr); */
      ierr = PetscViewerASCIIPrintf(viewer,"  DMDARepart: parent comm size reduction factor = %D\n",red->nsubcomm_factor);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"  DMDARepart: subcomm_size = %D\n",red->nsubcomm_size);CHKERRQ(ierr);
      ierr = PetscViewerGetSubViewer(viewer,red->subcomm->sub_comm,&subviewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      if (red->subcomm->parent_rank_active_in_subcomm) {
        /*ierr = DMView(red->dmrepart,subviewer);CHKERRQ(ierr);*/
        ierr = KSPView(red->ksp,subviewer);CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
      ierr = PetscViewerRestoreSubViewer(viewer,red->subcomm->sub_comm,&subviewer);CHKERRQ(ierr);
    }
  } else if (isstring) {
    ierr = PetscViewerStringSPrintf(viewer," DMDARepart preconditioner");CHKERRQ(ierr);
  } else {
    SETERRQ1(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"Viewer type %s not supported for PC DMDARepart",((PetscObject)viewer)->type_name);
  }
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
PetscErrorCode PCCreate_DMDARepart(PC pc)
{
  PetscErrorCode   ierr;
  PC_DMDARepart    *red;
  PetscMPIInt      size;
  PetscBool        log;

  PetscFunctionBegin;

  /*
     We cannot provide -XXX_pc_dmdarepart_log as an arg to enable logging of separate pcdmdarepart objects.
     This is a limitation enforced by petsc logging infrastructure which doesn't accept a communicator as an arguement.
     As a result all logging objects must be valid on PETSC_COMM_WORLD.
   */
  log = PETSC_FALSE;
  PetscOptionsGetBool(NULL,NULL,"-pc_dmdarepart_log",&log,NULL);
  if (log) {
    if (_PCDMDARepart_PCSetUpMatrixStage == -1) { ierr = PetscLogStageRegister("PCRprt_SetUpMat",&_PCDMDARepart_PCSetUpMatrixStage);CHKERRQ(ierr); }
    if (_PCDMDARepart_PCApplyStage == -1) { ierr = PetscLogStageRegister("PCRprt_Apply",&_PCDMDARepart_PCApplyStage);CHKERRQ(ierr); }
  }

  ierr = PetscNewLog(pc,&red);CHKERRQ(ierr);
  pc->data            = (void*)red;

  red->nsubcomm_factor = 2;
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pc),&size);CHKERRQ(ierr);
  red->nsubcomm_size   = size;
  if (size == 1) {
    red->nsubcomm_factor = 1;
  }

  red->log         = log;
  red->isin        = NULL;
  red->scatter     = NULL;
  red->xred        = NULL;
  red->xsub        = NULL;
  red->ysub        = NULL;
  red->permutation = NULL;
  red->Bsub        = NULL;
  red->dmrepart    = NULL;
  red->ksp         = NULL;

  pc->ops->apply           = PCApply_DMDARepart;
  pc->ops->applytranspose  = 0;
  pc->ops->setup           = PCSetUp_DMDARepart;
  pc->ops->destroy         = PCDestroy_DMDARepart;
  pc->ops->reset           = PCReset_DMDARepart;
  pc->ops->setfromoptions  = PCSetFromOptions_DMDARepart;
  pc->ops->view            = PCView_DMDARepart;

  PetscFunctionReturn(0);
}
EXTERN_C_END

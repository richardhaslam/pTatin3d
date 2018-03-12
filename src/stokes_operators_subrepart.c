/*
   A prototype of a matrix-free method designed to operate on heterogeneous
   compute nodes.

   It is assumed that moving data is cheap (or even free if we are on a
   shared-memory/UMA domain with the appropriate tools), but that we would like
   to process arbitrary subsets of elements on each rank, utilizing different
   kernels (in particular, we're interested in using AVX on some ranks and CUDA
   on others).

   The use case in mind is a heterogeneous compute node with a single
   coprocessor and several CPU cores, when we would like to use a flat MPI
   paradigm.

   This prototype does not implement general redistribution of work. Instead,
   it assumes that the MPI communicator on which the velocity DA lives (comm_u)
   can be decomposed into sub-communicators representing shared-memory domains,
   each with exactly one CUDA-enabled GPU available.

   This implementation assumes working
    * MPI-3 Shared Memory features,
    * CUDA, and
    * AVX.
   */

#include <petscfe.h>
#include <petsc/ptatin_petsc_hash.h> /* ptatin supplied header copied from PETSc source tree
                                        In PETSc 3.8, this header is moved to $PETSC_DIR/include/petsc/private */
#include <ptatin3d.h>
#include <ptatin3d_stokes.h>
#include <dmda_element_q2p1.h>
#include <stokes_operators.h>
#include <element_utils_q2.h>
#include <element_utils_q1.h>
#include <immintrin.h>

extern PetscLogEvent MAT_MultMFA11_SUP;
extern PetscLogEvent MAT_MultMFA11_stp;
extern PetscLogEvent MAT_MultMFA11_sub;
extern PetscLogEvent MAT_MultMFA11_rto;
extern PetscLogEvent MAT_MultMFA11_rfr;
extern PetscLogEvent MAT_MultMFA11_cto;
extern PetscLogEvent MAT_MultMFA11_ker;
extern PetscLogEvent MAT_MultMFA11_cfr;

#ifndef __FMA__
#  define _mm256_fmadd_pd(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b),c)
#endif

#define ALIGN32 __attribute__((aligned(32))) /* AVX packed instructions need 32-byte alignment */

typedef struct _p_MFA11SubRepart *MFA11SubRepart;

struct _p_MFA11SubRepart {
  /*
   <no postfix> : refers to quantities used in the "normal" case without repartitioning
   _repart      : refers to quantities which are used after repartitioning (these can
                  be the same as the unrepartitioned values, say if we use the same array)
   _remote      : refers to quantities which are used to communicate with rank 0

   Note that not all fields are used on all ranks. For example, elnidx_u_repart is only
   used when rank_sub == 0.
   */
  PetscObjectState state;
  const PetscInt   *elnidx_u;
  MPI_Comm         comm_sub;
  PetscMPIInt      size_sub,rank_sub;
  PetscInt         nel_sub,nen_u,nel,nel_remote,nel_repart,el_offset,nnodes,
                   nnodes_remote,nnodes_repart, nodes_offset;
  PetscInt         *nel_remote_in,*nnodes_remote_in,*nodes_remote,*elnidx_u_repart;
  PetscScalar      *mem_ufield_repart,*ufield_repart_base,*mem_Yu_repart,*Yu_repart_base;
  PetscBool        win_ufield_repart_allocated,win_Yu_repart_allocated;
  MPI_Win          win_ufield_repart,win_Yu_repart;
  MFA11CUDA        cudactx;
};

#undef __FUNCT__
#define __FUNCT__ "TransferQPData_A11_SubRepart"
static PetscErrorCode TransferQPData_A11_SubRepart(MFA11SubRepart ctx,PetscReal (*wp)[NQP],Quadrature volQ,QPntVolCoefStokes *all_gausspoints,PetscReal *gaussdata_w_repart_base,MPI_Win win_gaussdata_w_repart)
{
  PetscErrorCode          ierr;
  PetscInt                i,e;
  const QPntVolCoefStokes *cell_gausspoints;

  PetscFunctionBeginUser;

  /* Create and populate data in gaussdata_w_repart_base. This is the quadrature-point-wise
     information needed to apply the operator. On rank_sub 0, this
     has information for all the elements which will be processed.  */
  if (ctx->rank_sub) {
    PetscReal * const gaussdata_w_remote = &gaussdata_w_repart_base[NQP*ctx->el_offset];
    for (e=ctx->nel_repart; e<ctx->nel; e++) {
      PetscInt e_remote = e - ctx->nel_repart; /* element number in set to send */
      ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,(QPntVolCoefStokes**)&cell_gausspoints);CHKERRQ(ierr);
      for (i=0; i<NQP; i++) gaussdata_w_remote[e_remote*NQP + i] = cell_gausspoints[i].eta * (*wp)[i];
    }
  } else {
    for (e=0; e<ctx->nel; e++) {
      ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,(QPntVolCoefStokes**)&cell_gausspoints);CHKERRQ(ierr);
      for (i=0; i<NQP; i++) gaussdata_w_repart_base[e*NQP + i] = cell_gausspoints[i].eta * (*wp)[i];
    }
  }

  /* Note: a possible optimization is to delay this until we actually need this data */
  /* Synchronize (not sure if this is the optimal set of commands) */
  ierr = MPI_Win_sync(win_gaussdata_w_repart);CHKERRQ(ierr);
  ierr = MPI_Barrier(ctx->comm_sub);CHKERRQ(ierr);
  ierr = MPI_Win_sync(win_gaussdata_w_repart);CHKERRQ(ierr); /* apparently required on some systems */

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TransferCoordinates_A11_SubRepart"
static PetscErrorCode TransferCoordinates_A11_SubRepart(MFA11SubRepart ctx,const PetscReal *LA_gcoords,PetscReal *LA_gcoords_repart_base,MPI_Win win_LA_gcoords_repart)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;

  if (ctx->rank_sub) {
    PetscReal * const LA_gcoords_remote = &LA_gcoords_repart_base[NSD*ctx->nodes_offset];
    for (i=0;i<ctx->nnodes_remote;++i) {
      PetscInt d;
      for(d=0;d<NSD;++d){
        LA_gcoords_remote[NSD*i+d] = LA_gcoords[NSD*ctx->nodes_remote[i]+d];
      }
    }
  } else {
    ierr = PetscMemcpy(LA_gcoords_repart_base,LA_gcoords,NSD*ctx->nnodes*sizeof(PetscScalar));CHKERRQ(ierr);
  }

  /* Note: a possible optimization is to delay this until we actually need this data */
  /* Synchronize (not sure if this is the optimal set of commands) */
  ierr = MPI_Win_sync(win_LA_gcoords_repart);CHKERRQ(ierr);
  ierr = MPI_Barrier(ctx->comm_sub);CHKERRQ(ierr);
  ierr = MPI_Win_sync(win_LA_gcoords_repart);CHKERRQ(ierr); /* apparently required on some systems */

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TransferUfield_A11_SubRepart"
static PetscErrorCode TransferUfield_A11_SubRepart(MFA11SubRepart ctx,PetscScalar *ufield)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;

  /* Gather velocity field entries from all ranks in shared mem array
     starting at ufield_repart_base. */

  /* Rank_sub 1,2,.. poke data directly into the shared array, and rank_sub 0
     simply copies */
  if (ctx->rank_sub) {
    PetscScalar * const ufield_remote = &ctx->ufield_repart_base[NSD*ctx->nodes_offset];
    for (i=0;i<ctx->nnodes_remote;++i) {
      PetscInt d;
      for(d=0;d<NSD;++d){
        ufield_remote[NSD*i+d] = ufield[NSD*ctx->nodes_remote[i]+d];
      }
    }
  } else {
    ierr = PetscMemcpy(ctx->ufield_repart_base,ufield,NSD*ctx->nnodes*sizeof(PetscScalar));CHKERRQ(ierr);
  }

  /* Synchronize (not sure if this is the optimal set of commands) */
  ierr = MPI_Win_sync(ctx->win_ufield_repart);CHKERRQ(ierr);
  ierr = MPI_Barrier(ctx->comm_sub);CHKERRQ(ierr);
  ierr = MPI_Win_sync(ctx->win_ufield_repart);CHKERRQ(ierr); /* apparently required on some systems */

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TransferYu_A11_SubRepart"
static PetscErrorCode TransferYu_A11_SubRepart(MFA11SubRepart ctx,PetscScalar *Yu)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;

  /* Synchronize (not sure if this is the optimal set of commands) */
  ierr = MPI_Win_sync(ctx->win_ufield_repart);CHKERRQ(ierr);
  ierr = MPI_Barrier(ctx->comm_sub);CHKERRQ(ierr);
  ierr = MPI_Win_sync(ctx->win_ufield_repart);CHKERRQ(ierr); /* apparently required on some systems */

  /* Accumulate into Yu on rank_sub 1,2,.. just copy on rank_sub 0*/
  if (ctx->rank_sub) {
    PetscScalar * const Yu_remote = &ctx->Yu_repart_base[NSD*ctx->nodes_offset];
    for (i=0; i<ctx->nnodes_remote; ++i) {
      PetscInt d;
      for(d=0;d<NSD;++d){
        Yu[NSD*ctx->nodes_remote[i]+d] += Yu_remote[NSD*i+d];
      }
    }
  } else {
#if 1
    for (i=0;i<NSD*ctx->nnodes;++i){
      Yu[i] += ctx->Yu_repart_base[i];
    }
#else
    ierr = PetscMemcpy(Yu,ctx->Yu_repart_base,NSD*ctx->nnodes*sizeof(PetscScalar));CHKERRQ(ierr);
#endif
  }

  /* Synchronize (not sure if this is the optimal set of commands) */
  ierr = MPI_Win_sync(ctx->win_ufield_repart);CHKERRQ(ierr);
  ierr = MPI_Barrier(ctx->comm_sub);CHKERRQ(ierr);
  ierr = MPI_Win_sync(ctx->win_ufield_repart);CHKERRQ(ierr); /* apparently required on some systems */

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFA11SetUp_SubRepart"
PetscErrorCode MFA11SetUp_SubRepart(MatA11MF mf)
{
  PetscErrorCode  ierr;
  DM              dau = mf->daUVW;
  PetscInt        i;
  MFA11SubRepart  ctx;
  PetscHashI      nodes_remote_inv;

  // TODO: this whole function wantonly uses MPI calls. These should be collected
  //       to minimize communication, and the heavy transfer for the elements should be done with shared memory.
  //       It might be best just to do one MPI_Allgather to get the number of remote elements and nodes from each rank,
  //       on all ranks. Then each rank can independently compute its offsets.

  PetscFunctionBeginUser;
  if (mf->ctx) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(MAT_MultMFA11_SUP,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);

  ctx->state = 0;

  ctx->win_ufield_repart_allocated = PETSC_FALSE;
  ctx->ufield_repart_base          = NULL;
  ctx->mem_ufield_repart           = NULL;
  ctx->win_Yu_repart_allocated     = PETSC_FALSE;
  ctx->Yu_repart_base              = NULL;
  ctx->mem_Yu_repart               = NULL;

  /* Define a subcomm relative to the local shared-memory domain.
     This implementation assumes that you have exactly one CUDA-enabled
     GPU per shared-memory domain. */
  ierr = MPI_Comm_split_type(PetscObjectComm((PetscObject)dau),MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&ctx->comm_sub);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(ctx->comm_sub,&ctx->rank_sub);CHKERRQ(ierr);
  ierr = MPI_Comm_size(ctx->comm_sub,&ctx->size_sub);CHKERRQ(ierr);

  ierr = DMDAGetElements_pTatinQ2P1(dau,&ctx->nel,&ctx->nen_u,&ctx->elnidx_u);CHKERRQ(ierr);

  /* Decide on a repartitioning of the elements
     We assume that we send a contiguous chunk of elements from
     each non-zero rank on comm_sub to rank 0 on comm_sub, and simply process
     fewer elements from the original list on ranks 1,... size_comm-1. On rank 0,
     we process a larger array.
     */
  ierr = MPI_Allreduce(&ctx->nel,&ctx->nel_sub,1,MPIU_INT,MPI_SUM,ctx->comm_sub);CHKERRQ(ierr);

  /* To compute how many elements to send to rank_sub 0, we compute an equal partition
     of the elements among the non-zero rank_subs, compute a fraction of this to
     send to rank_sub 0, and round down. */
  if (ctx->size_sub == 1) {
    ctx->nel_repart = ctx->nel;
    ierr = PetscPrintf(PETSC_COMM_SELF,"WARNING: sub-communicator of size 1. Unless you are testing, something is wrong with your use of the repartitioned operator. Perhaps you are in an environment in which MPI_Comm_split_type() is not behaving as expected.\n");CHKERRQ(ierr);
  } else {
    PetscReal remotefrac = 0.78; /* How much of the work on each rank to offload to rank_sub 0 */
    ierr = PetscOptionsGetReal(NULL,NULL,"-subrepart_frac",&remotefrac,NULL);CHKERRQ(ierr);
    if (remotefrac < 0.0 || remotefrac > 1.0) {
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-subrepart_frac must be in [0,1]");
    }
    const PetscInt nel_repart_small = (PetscInt) ((1.0-remotefrac) * ctx->nel_sub / (ctx->size_sub-1));
    const PetscInt nel_repart_big   = ctx->nel_sub - ((ctx->size_sub-1) * nel_repart_small);
    ctx->nel_repart                 = ctx->rank_sub ?  nel_repart_small : nel_repart_big;
  }

  /* (non-zero ranks) Determine the unique nodes corresponding to the elements
     which will be sent to rank_sub 0. Obtain a count. */
  ctx->nel_remote = ctx->rank_sub ? ctx->nel - ctx->nel_repart: 0;

  if (ctx->nel_remote < 0) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"SubRepart's current partitioning scheme does not allow for the current partition. Perhaps you have too few ranks in a sub communicator.");
  }

  /* Overestimate the number of remote nodes, then get the actual number
     of unique nodes by sorting. This could, with the help of a hash table
     be done in linear time, and the wasted memory here could be recovered */
  if (ctx->rank_sub) {
    ierr = PetscMalloc1(ctx->nel_remote*ctx->nen_u,&ctx->nodes_remote);CHKERRQ(ierr);
    ierr = PetscMemcpy(ctx->nodes_remote,&ctx->elnidx_u[ctx->nel_repart*ctx->nen_u],
        ctx->nel_remote*ctx->nen_u*sizeof(PetscInt));CHKERRQ(ierr);
    ctx->nnodes_remote = ctx->nen_u*ctx->nel_remote; /* will be updated */
    ierr = PetscSortRemoveDupsInt(&ctx->nnodes_remote,ctx->nodes_remote);CHKERRQ(ierr);
  } else {
    ctx->nnodes_remote = 0;
  }

  /* Create a map from the entries of the local vector here
     to just the entries which will be needed on rank_sub 0, without the offset.
     That is, the inverse of nodes_remote[].  This is a hash table.
     This lets us do a rank-local "scatter" to an array
     which we send to the right place on rank_sub 0. This table is only used
     here, and is then destroyed */
  PetscHashICreate(nodes_remote_inv);
  for(i=0;i<ctx->nnodes_remote;++i){
    PetscHashIAdd(nodes_remote_inv,ctx->nodes_remote[i],i);
  }

  /* On rank_sub 0, obtain the number of nodes being sent from
     each other rank in comm_sub */
  if (!ctx->rank_sub) {
    ierr = PetscMalloc1(ctx->size_sub,&ctx->nnodes_remote_in);CHKERRQ(ierr);
  }
  ierr = MPI_Gather(&ctx->nnodes_remote,1,MPIU_INT,ctx->nnodes_remote_in,1,MPIU_INT,0,ctx->comm_sub);CHKERRQ(ierr);

  /* Compute the number of nodes on each rank
     This is only actually used on sub_rank 0.  */
  {
    PetscInt X,Y,Z; /* local sizes */
    ierr = DMDAGetGhostCorners(dau,NULL,NULL,NULL, &X,&Y,&Z);CHKERRQ(ierr);
    ctx->nnodes = X * Y * Z;
  }

  /* Calculate the number of nodes on each rank, after repartitioning.
     This only changes on rank_sub 0 */
  if (ctx->rank_sub) {
    ctx->nnodes_repart = ctx->nnodes; /* We use the existing array */
  } else {
    ctx->nnodes_repart = ctx->nnodes;
    for (i=1;i<ctx->size_sub;++i) ctx->nnodes_repart += ctx->nnodes_remote_in[i];
  }

  /* Compute an offset for each rank_sub > 0 to process its element
     list to correspond to the numbering on rank_sub 0 */
  {
    PetscInt *nodes_offsets = NULL;
    if (!ctx->rank_sub) {
      ierr = PetscMalloc1(ctx->size_sub,&nodes_offsets);
      nodes_offsets[0] = 0;
      nodes_offsets[1] = ctx->nnodes;
      for (i=2;i<ctx->size_sub;++i) nodes_offsets[i] = nodes_offsets[i-1] + ctx->nnodes_remote_in[i-1];
    }
    ierr = MPI_Scatter(nodes_offsets,1,MPIU_INT,&ctx->nodes_offset,1,MPIU_INT,0,ctx->comm_sub);CHKERRQ(ierr);
    ierr = PetscFree(nodes_offsets);
  }

  /* On rank 0, we need to know how many elements from each other rank, in order
     to receive data in the correct parts of the concatenated element, nodes,
     and coordinate arrays */
  ierr = PetscMalloc1(ctx->size_sub,&ctx->nel_remote_in);CHKERRQ(ierr);
  ierr = MPI_Gather(&ctx->nel_remote,1,MPIU_INT,ctx->nel_remote_in,1,MPIU_INT,0,ctx->comm_sub);CHKERRQ(ierr);

  /* Compute an element offset for each rank */
  {
    PetscInt *el_offsets = NULL;
    if (!ctx->rank_sub) {
      ierr = PetscMalloc1(ctx->size_sub,&el_offsets);
      el_offsets[0] = 0;
      el_offsets[1] = ctx->nel;
      for (i=2;i<ctx->size_sub;++i) el_offsets[i] = el_offsets[i-1] + ctx->nel_remote_in[i-1];
    }
    ierr = MPI_Scatter(el_offsets,1,MPIU_INT,&ctx->el_offset,1,MPIU_INT,0,ctx->comm_sub);CHKERRQ(ierr);
    ierr = PetscFree(el_offsets);
  }

  /* Send and receive elements being offloaded to rank 0, defining element-->node
     maps for the repartition.  On rank_sub > 0, this involves computing the
     element indices to send to rank 0 by mapping them into the local ordering
     for the nodes which are used there, and then offsetting.  */
  if (ctx->rank_sub) {
    PetscInt *elnidx_u_remote;
    PetscInt i;
    ierr = PetscMalloc1(ctx->nel_remote*ctx->nen_u,&elnidx_u_remote);CHKERRQ(ierr);
    for (i=0;i<ctx->nen_u*ctx->nel_remote;++i){
      PetscInt ind = ctx->elnidx_u[ctx->nen_u*ctx->nel_repart + i];
      {
        PetscInt val ;
        PetscHashIMap(nodes_remote_inv,ind,val); /* local index --> index in ctx->nodes_remote */
        elnidx_u_remote[i] = val + ctx->nodes_offset; /* index in ctx->nodes_remote --> index for ctx->elnidx_u_repart on rank 0 */
      }
    }

    PetscHashIDestroy(nodes_remote_inv);
    ierr = MPI_Send(elnidx_u_remote,ctx->nen_u*ctx->nel_remote,MPIU_INT,0,0,ctx->comm_sub);CHKERRQ(ierr);
    ierr = PetscFree(elnidx_u_remote);CHKERRQ(ierr);
  } else {
    PetscMPIInt i;
    PetscInt    el_offset=ctx->nel;
    ierr = PetscMalloc1(ctx->nel_repart*ctx->nen_u,&ctx->elnidx_u_repart);CHKERRQ(ierr);
    ierr = PetscMemcpy(ctx->elnidx_u_repart,ctx->elnidx_u,ctx->nel*ctx->nen_u*sizeof(PetscInt));CHKERRQ(ierr);

    /* Receive a chunk of elements (sets of 27 nodes for Q2 elements) from each other rank */
    // TODO: replace with shared memory operation (just as we have for other arrays in this impl)
    for(i=1;i<ctx->size_sub;++i)  {
      ierr = MPI_Recv(&ctx->elnidx_u_repart[ctx->nen_u*el_offset],ctx->nen_u*ctx->nel_remote_in[i],MPIU_INT,i,0,ctx->comm_sub,MPI_STATUS_IGNORE);CHKERRQ(ierr);
      el_offset += ctx->nel_remote_in[i];
    }
  }

#if 0
  { /* Print stats on nodes and el for normal, remote, and repart */
    PetscErrorCode ierr;
    PetscMPIInt    rank,size,r;
    ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\033[32m= SubRepart A11 SetUp Info ===\033[0m\n");CHKERRQ(ierr);
    for(r=0;r<size;++r) {
      if(r==rank) {
        ierr = PetscPrintf(PETSC_COMM_SELF,"\033[32m[%d/%d]\033[0m ",rank,size);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF,"\033[32m[n%d/%d]\033[0m\n",ctx->rank_sub,ctx->size_sub);CHKERRQ(ierr);
        if (!ctx->rank_sub) {
          ierr = PetscPrintf(PETSC_COMM_SELF, " nel_sub:%d \n  nnodes_remote_in[]: ",ctx->nel_sub);CHKERRQ(ierr);
          for (i=0;i<ctx->size_sub;++i) {
            ierr = PetscPrintf(PETSC_COMM_SELF, "%d ",ctx->nnodes_remote_in[i]);CHKERRQ(ierr);
          }
          ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_SELF, " nel_remote_in[]: ",ctx->nel_sub);CHKERRQ(ierr);
          for (i=0;i<ctx->size_sub;++i) {
            ierr = PetscPrintf(PETSC_COMM_SELF, "%d ",ctx->nel_remote_in[i]);CHKERRQ(ierr);
          }
          ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, " [Unpartitioned] nel: %4d  nnodes: %4d \n",ctx->nel,ctx->nnodes);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " [Remote]        nel: %4d  nnodes: %4d\n",ctx->nel_remote,ctx->nnodes_remote);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " [Repartitioned] nel: %4d  nnodes: %4d \n",ctx->nel_repart,ctx->nnodes_repart);CHKERRQ(ierr);
#       if 1
        /* Compute the range of indices used after repartitioning, and the percentage of indices in this range that are used by
           the repartitioned set of elements */
        {
          PetscInt    imin=PETSC_MAX_INT, imax=0, icount, icountUnique;
          PetscInt    *indices;
          PetscScalar efficiency;
          icount = ctx->nel_repart*ctx->nen_u;
          ierr = PetscMalloc1(icount,&indices);CHKERRQ(ierr);
          for (i=0; i<icount; ++i) {
            PetscInt curr;
            if (ctx->rank_sub == 0) {
              curr = ctx->elnidx_u_repart[i];
            } else {
              curr = ctx->elnidx_u[i];
            }
            imin = PetscMin(imin,curr);
            imax = PetscMax(imax,curr);
            indices[i] = curr;
          }
          icountUnique = icount; /* will be overwritten */
          ierr = PetscSortRemoveDupsInt(&icountUnique,indices);CHKERRQ(ierr);
          efficiency = icountUnique/((PetscScalar)(imax-imin));
          ierr = PetscPrintf(PETSC_COMM_SELF, " [Repartitioned] element dof range: %d-%d #unique: %d, usage: %0.2f pct.\n",imin,imax,icountUnique,efficiency*100.0);CHKERRQ(ierr);
          ierr = PetscFree(indices);CHKERRQ(ierr);
        }
#       endif
      }
      ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\033[32m========================\033[0m\n");CHKERRQ(ierr);
  }
#endif

  /* Set up CUDA context */
  if (!ctx->rank_sub) {
    ierr = PetscMalloc1(1,&ctx->cudactx);CHKERRQ(ierr);
    ierr = MFA11CUDA_SetUp(ctx->cudactx);CHKERRQ(ierr);
  }

  mf->ctx=ctx;
  ierr = PetscLogEventEnd(MAT_MultMFA11_SUP,0,0,0,0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFA11Destroy_SubRepart"
PetscErrorCode MFA11Destroy_SubRepart(MatA11MF mf)
{
  PetscErrorCode  ierr;
  MFA11SubRepart  ctx;

  PetscFunctionBeginUser;
  ctx = (MFA11SubRepart)mf->ctx;

  /* Free shared memory windows */
  if (ctx->win_ufield_repart_allocated) {
    ierr = MPI_Win_unlock_all(ctx->win_ufield_repart);CHKERRQ(ierr);
    ierr = MPI_Win_free(&ctx->win_ufield_repart);CHKERRQ(ierr);
    ctx->win_ufield_repart_allocated = PETSC_FALSE;
    ctx->ufield_repart_base = NULL;
    ctx->mem_ufield_repart = NULL;
  }
  if (ctx->win_Yu_repart_allocated) {
    ierr = MPI_Win_unlock_all(ctx->win_Yu_repart);CHKERRQ(ierr);
    ierr = MPI_Win_free(&ctx->win_Yu_repart);CHKERRQ(ierr);
    ctx->win_Yu_repart_allocated = PETSC_FALSE;
    ctx->Yu_repart_base = NULL;
    ctx->mem_Yu_repart = NULL;
  }

  if (ctx->rank_sub) {
    ierr = PetscFree(ctx->nodes_remote);CHKERRQ(ierr);
  } else {
    ierr = PetscFree(ctx->elnidx_u_repart);CHKERRQ(ierr);
    ierr = PetscFree(ctx->nnodes_remote_in);CHKERRQ(ierr);
    ierr = PetscFree(ctx->nel_remote_in);CHKERRQ(ierr);
  }

  /* Destroy CUDA ctx */
  if (!ctx->rank_sub) {
    ierr = MFA11CUDA_CleanUp(ctx->cudactx);CHKERRQ(ierr);
    ierr = PetscFree(ctx->cudactx);CHKERRQ(ierr);
  }

  ierr = PetscFree(mf->ctx);CHKERRQ(ierr);
  mf->ctx = NULL;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFStokesWrapper_A11_SubRepart"
PetscErrorCode MFStokesWrapper_A11_SubRepart(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{
  PetscErrorCode          ierr;
  MFA11SubRepart          ctx;
  DM                      cda;
  Vec                     gcoords;
  const PetscReal         *LA_gcoords;
  PetscReal               x1[3],w1[3],B[3][3],D[3][3],w[NQP];
  PetscInt                i,j,k,e;
  QPntVolCoefStokes       *all_gausspoints;
  PetscReal               *mem_gaussdata_w_repart,*gaussdata_w_repart_base,*mem_LA_gcoords_repart,*LA_gcoords_repart_base;
  PetscBool               win_LA_gcoords_repart_allocated=PETSC_FALSE,win_gaussdata_w_repart_allocated=PETSC_FALSE ;
  MPI_Win                 win_LA_gcoords_repart,win_gaussdata_w_repart;

  PetscFunctionBeginUser;
  ctx = (MFA11SubRepart)mf->ctx;

  ierr = PetscLogEventBegin(MAT_MultMFA11_stp,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscDTGaussQuadrature(3,-1,1,x1,w1);CHKERRQ(ierr);
  for (i=0; i<3; i++) {
    B[i][0] = .5*(PetscSqr(x1[i]) - x1[i]);
    B[i][1] = 1 - PetscSqr(x1[i]);
    B[i][2] = .5*(PetscSqr(x1[i]) + x1[i]);
    D[i][0] = x1[i] - .5;
    D[i][1] = -2*x1[i];
    D[i][2] = x1[i] + .5;
  }
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        w[(i*3+j)*3+k] = w1[i] * w1[j] * w1[k];}}}

  /* setup for coords */
  ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
  ierr = VecGetArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);

  ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);

  ierr = PetscLogEventEnd(MAT_MultMFA11_stp,0,0,0,0);CHKERRQ(ierr);

  if(ctx->state != mf->state) {
    ierr = PetscLogEventBegin(MAT_MultMFA11_SUP,0,0,0,0);CHKERRQ(ierr);
    /* Allocate space for repartitioned node-wise fields, and repartitioned elementwise data */

    /* (re)Allocate shared memory windows on rank_sub 0*/
    if (ctx->win_ufield_repart_allocated) {
      ctx->ufield_repart_base = NULL;
      ctx->mem_ufield_repart = NULL;
      ierr = MPI_Win_unlock_all(ctx->win_ufield_repart);CHKERRQ(ierr);
      ierr = MPI_Win_free(&ctx->win_ufield_repart);CHKERRQ(ierr);
      ctx->win_ufield_repart_allocated = PETSC_FALSE;
    }
    {
      MPI_Aint          sz;
      PetscMPIInt       du;
      const PetscMPIInt nnodes_win = ctx->rank_sub ? 0 : ctx->nnodes_repart;
      ierr = MPI_Win_allocate_shared(NSD*nnodes_win*sizeof(PetscScalar),sizeof(PetscScalar),MPI_INFO_NULL,ctx->comm_sub,&ctx->mem_ufield_repart,&ctx->win_ufield_repart);CHKERRQ(ierr);
      ctx->win_ufield_repart_allocated = PETSC_TRUE;
      ierr = MPI_Win_shared_query(ctx->win_ufield_repart,MPI_PROC_NULL,&sz,&du,&ctx->ufield_repart_base);CHKERRQ(ierr);
      ierr = MPI_Win_lock_all(MPI_MODE_NOCHECK,ctx->win_ufield_repart);CHKERRQ(ierr);
    }
    if (ctx->win_Yu_repart_allocated) {
      ctx->Yu_repart_base = NULL;
      ctx->mem_Yu_repart = NULL;
      ierr = MPI_Win_unlock_all(ctx->win_Yu_repart);CHKERRQ(ierr);
      ierr = MPI_Win_free(&ctx->win_Yu_repart);CHKERRQ(ierr);
      ctx->win_Yu_repart_allocated = PETSC_FALSE;
    }
    {
      MPI_Aint          sz;
      PetscMPIInt       du;
      const PetscMPIInt nnodes_win = ctx->rank_sub ? 0 : ctx->nnodes_repart;
      ierr = MPI_Win_allocate_shared(NSD*nnodes_win*sizeof(PetscScalar),sizeof(PetscScalar),MPI_INFO_NULL,ctx->comm_sub,&ctx->mem_Yu_repart,&ctx->win_Yu_repart);CHKERRQ(ierr);
      ctx->win_Yu_repart_allocated = PETSC_TRUE;
      ierr = MPI_Win_shared_query(ctx->win_Yu_repart,MPI_PROC_NULL,&sz,&du,&ctx->Yu_repart_base);CHKERRQ(ierr);
      ierr = MPI_Win_lock_all(MPI_MODE_NOCHECK,ctx->win_Yu_repart);CHKERRQ(ierr);
    }
    {
      MPI_Aint          sz;
      PetscMPIInt       du;
      const PetscMPIInt nnodes_win = ctx->rank_sub ? 0 : ctx->nnodes_repart;
      ierr = MPI_Win_allocate_shared(NSD*nnodes_win*sizeof(PetscReal),sizeof(PetscReal),MPI_INFO_NULL,ctx->comm_sub,&mem_LA_gcoords_repart,&win_LA_gcoords_repart);CHKERRQ(ierr);
      win_LA_gcoords_repart_allocated = PETSC_TRUE;
      ierr = MPI_Win_shared_query(win_LA_gcoords_repart,MPI_PROC_NULL,&sz,&du,&LA_gcoords_repart_base);CHKERRQ(ierr);
      ierr = MPI_Win_lock_all(MPI_MODE_NOCHECK,win_LA_gcoords_repart);CHKERRQ(ierr);
    }
    {
      MPI_Aint          sz;
      PetscMPIInt       du;
      const PetscMPIInt nel_win = ctx->rank_sub ? 0 : ctx->nel_repart;
      ierr = MPI_Win_allocate_shared(NQP*nel_win*sizeof(PetscReal),sizeof(PetscReal),MPI_INFO_NULL,ctx->comm_sub,&mem_gaussdata_w_repart,&win_gaussdata_w_repart);CHKERRQ(ierr);
      win_gaussdata_w_repart_allocated = PETSC_TRUE;
      ierr = MPI_Win_shared_query(win_gaussdata_w_repart,MPI_PROC_NULL,&sz,&du,&gaussdata_w_repart_base);CHKERRQ(ierr);
      ierr = MPI_Win_lock_all(MPI_MODE_NOCHECK,win_gaussdata_w_repart);CHKERRQ(ierr);
    }

    /* Send required quadrature-pointwise data to rank_sub 0 */
    ierr = TransferQPData_A11_SubRepart(ctx,&w,volQ,all_gausspoints,gaussdata_w_repart_base,win_gaussdata_w_repart);CHKERRQ(ierr);

    /* Send required coordinate data to rank_sub 0 */
    ierr = TransferCoordinates_A11_SubRepart(ctx,LA_gcoords,LA_gcoords_repart_base,win_LA_gcoords_repart);CHKERRQ(ierr);

    ctx->state = mf->state;
    ierr = PetscLogEventEnd(MAT_MultMFA11_SUP,0,0,0,0);CHKERRQ(ierr);
  }

  /* Send required ufield data to rank_sub 0 */
  ierr = PetscLogEventBegin(MAT_MultMFA11_rto,0,0,0,0);CHKERRQ(ierr);
  ierr = TransferUfield_A11_SubRepart(ctx,ufield);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(MAT_MultMFA11_rto,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscLogEventBegin(MAT_MultMFA11_sub,0,0,0,0);CHKERRQ(ierr);

#if defined(_OPENMP)
#define OPENMP_CHKERRQ(x)
#else
#define OPENMP_CHKERRQ(x)   CHKERRQ(x)
#endif
  if (ctx->rank_sub) {
    /* Rank_sub > 0 implementation. This is the same as the AVX implementation,
       over a smaller set of elements (ctx->nel_repart) */
#if defined(_OPENMP)
#pragma omp parallel for private(i)
#endif
    for (e=0;e<ctx->nel_repart;e+=NEV) {
      PetscScalar elu[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32,elx[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32,elv[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32;
      PetscScalar dx[3][3][NQP][NEV] ALIGN32,dxdet[NQP][NEV],du[3][3][NQP][NEV] ALIGN32,dv[3][3][NQP][NEV] ALIGN32;
      const QPntVolCoefStokes *cell_gausspoints[NEV];
      PetscInt ee,l;

      for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
        for (ee=0; ee<NEV; ee++) {
          PetscInt E = ctx->elnidx_u[ctx->nen_u*PetscMin(e+ee,ctx->nel_repart-1)+i]; /* Pad up to length NEV by duplicating last element */
          for (l=0; l<3; l++) {
            elx[l][i][ee] = LA_gcoords[3*E+l];
            elu[l][i][ee] = ufield[3*E+l];
          }
        }
      }
      for (ee=0; ee<NEV; ee++) {
        ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,PetscMin(e+ee,ctx->nel_repart-1),(QPntVolCoefStokes**)&cell_gausspoints[ee]);OPENMP_CHKERRQ(ierr);
      }

      ierr = PetscMemzero(dx,sizeof dx);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(D,B,B,GRAD,elx,dx[0]);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(B,D,B,GRAD,elx,dx[1]);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(B,B,D,GRAD,elx,dx[2]);OPENMP_CHKERRQ(ierr);

      ierr = JacobianInvertNEV_AVX(dx,dxdet);OPENMP_CHKERRQ(ierr);

      ierr = PetscMemzero(du,sizeof du);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(D,B,B,GRAD,elu,du[0]);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(B,D,B,GRAD,elu,du[1]);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(B,B,D,GRAD,elu,du[2]);OPENMP_CHKERRQ(ierr);

      ierr = QuadratureAction_A11_AVX(cell_gausspoints,dx,dxdet,w,du,dv);OPENMP_CHKERRQ(ierr);

      ierr = PetscMemzero(elv,sizeof elv);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(D,B,B,GRAD_TRANSPOSE,dv[0],elv);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(B,D,B,GRAD_TRANSPOSE,dv[1],elv);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(B,B,D,GRAD_TRANSPOSE,dv[2],elv);OPENMP_CHKERRQ(ierr);

#if defined(_OPENMP)
#pragma omp critical
#endif
      for (ee=0; ee<PetscMin(NEV,ctx->nel_repart-e); ee++) {
        for (i=0; i<NQP; i++) {
          PetscInt E = ctx->elnidx_u[ctx->nen_u*(e+ee)+i];
          for (l=0; l<3; l++) {
            Yu[3*E+l] += elv[l][i][ee];
          }
        }
      }
    }
  } else {
    /* Rank_sub 0 CUDA Implementation */
    ierr = PetscLogEventBegin(MAT_MultMFA11_cto,0,0,0,0);CHKERRQ(ierr);
    ierr = CopyTo_A11_CUDA(mf,ctx->cudactx,ctx->ufield_repart_base,LA_gcoords_repart_base,gaussdata_w_repart_base,ctx->nel_repart,ctx->nen_u,ctx->elnidx_u_repart,ctx->nnodes_repart);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(MAT_MultMFA11_cto,0,0,0,0);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(MAT_MultMFA11_ker,0,0,0,0);CHKERRQ(ierr);
    ierr = ProcessElements_A11_CUDA(ctx->cudactx,ctx->nen_u,NSD*ctx->nnodes_repart);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(MAT_MultMFA11_ker,0,0,0,0);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(MAT_MultMFA11_cfr,0,0,0,0);CHKERRQ(ierr);
    ierr = CopyFrom_A11_CUDA(ctx->cudactx,ctx->Yu_repart_base,NSD*ctx->nnodes_repart);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(MAT_MultMFA11_cfr,0,0,0,0);CHKERRQ(ierr);
  }
#undef OPENMP_CHKERRQ

  if (win_LA_gcoords_repart_allocated) {
    LA_gcoords_repart_base = NULL;
    mem_LA_gcoords_repart = NULL;
    ierr = MPI_Win_unlock_all(win_LA_gcoords_repart);CHKERRQ(ierr);
    ierr = MPI_Win_free(&win_LA_gcoords_repart);CHKERRQ(ierr);
    win_LA_gcoords_repart_allocated = PETSC_FALSE;
  }
  if (win_gaussdata_w_repart_allocated) {
    gaussdata_w_repart_base = NULL;
    mem_gaussdata_w_repart = NULL;
    ierr = MPI_Win_unlock_all(win_gaussdata_w_repart);CHKERRQ(ierr);
    ierr = MPI_Win_free(&win_gaussdata_w_repart);CHKERRQ(ierr);
    win_gaussdata_w_repart_allocated = PETSC_FALSE;
  }

  ierr = VecRestoreArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(MAT_MultMFA11_sub,0,0,0,0);CHKERRQ(ierr);

  /* Transfer and accumulate contributions to Yu from rank_sub 0 */
  ierr = PetscLogEventBegin(MAT_MultMFA11_rfr,0,0,0,0);CHKERRQ(ierr);
  ierr = TransferYu_A11_SubRepart(ctx,Yu);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(MAT_MultMFA11_rfr,0,0,0,0);CHKERRQ(ierr);

  /* Note that this is identical to the AVX impl, under the assumption
     that we are logging "useful flops" and that identical work is being performed.
     This has not been scrutinized */
  PetscLogFlops((ctx->nel * 9) * 3*NQP*(6+6+6));           /* 9 tensor contractions per element */
  PetscLogFlops(ctx->nel*NQP*(14 + 1/* division */ + 27)); /* 1 Jacobi inversion per element */
  PetscLogFlops(ctx->nel*NQP*(5*9+6+6+6*9));               /* 1 quadrature action per element */

  PetscFunctionReturn(0);
}

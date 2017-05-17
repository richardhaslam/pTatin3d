/*
A prototype of a matrix-free method
designed to operate on heterogeneous compute nodes.

That is, it is assumed that moving data is cheap (or even free if
we are on a shared-memory/UMA domain with the appropriate tools),
but that we would like to process arbitrary subsets of elements
on each rank, utilizing different kernels (in particular,
we're interested in using AVX on some ranks and CUDA on others).

This prototype does not implement general redistribution of work.
Instead, it assumes that the MPI communicator on which the
velocity DA lives (comm_u) can be decomposed into sub-communicators,
where different computational resources are available on
rank_sub 0 than on the remaining ranks. The use case in mind
is a heterogeneous compute node with a single coprocessor
and several CPU cores, when we would like to use a flat MPI
paradigm.

This current prototype is intended to allow you to choose
to use a modified AVX implementation, or the CUDA implementation,
on rank_sub 0, and the usual AVX implementation on other ranks.

Farther into the future, it would be natural to investigate
shared-memory abstractions (such as those supported by MPI-3)
which would allow for a more natural set of operations on a single
shared set elements to be processed per shared memory domain.
*/

#include <petscfe.h>
#include <../src/sys/utils/hash.h> /* not portable to prefix installs
                                     In master (PETSc 3.8), this header is
                                     moved to $PETSC_DIR/include/petsc/private) */
#include <ptatin3d.h>
#include <ptatin3d_stokes.h>
#include <dmda_element_q2p1.h>
#include <stokes_operators.h>
#include <element_utils_q2.h>
#include <element_utils_q1.h>
#include <immintrin.h>

extern PetscLogEvent MAT_MultMFA11_setup;
extern PetscLogEvent MAT_MultMFA11_sub;
extern PetscLogEvent MAT_MultMFA11_rto;
extern PetscLogEvent MAT_MultMFA11_rfr;

#ifndef __FMA__
#  define _mm256_fmadd_pd(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b),c)
#endif

#define ALIGN32 __attribute__((aligned(32))) /* AVX packed instructions need 32-byte alignment */

__attribute__((noinline))
static PetscErrorCode QuadratureAction_A11AVX_mod(const PetscReal *gaussdata_w,
					   PetscScalar dx[3][3][Q2_NODES_PER_EL_3D][NEV],
					   PetscScalar dxdet[Q2_NODES_PER_EL_3D][NEV],
					   PetscScalar du[3][3][Q2_NODES_PER_EL_3D][NEV],
					   PetscScalar dv[3][3][Q2_NODES_PER_EL_3D][NEV])
{
	PetscInt i,l,k,e;

	for (i=0; i<NQP; i++) {
		PetscScalar Du[6][NEV] ALIGN32,Dv[6][NEV] ALIGN32; /* Symmetric gradient with respect to physical coordinates, xx, yy, zz, xy+yx, xz+zx, yz+zy */
		__m256d dux[3][3],mhalf = _mm256_set1_pd(0.5),dvx[3][3];
		__m256d mdxdet = _mm256_load_pd(dxdet[i]);

		for (k=0; k<3; k++) { /* directions */
			__m256d dxk[3] = {_mm256_load_pd(dx[k][0][i]),_mm256_load_pd(dx[k][1][i]),_mm256_load_pd(dx[k][2][i])};
			for (l=0; l<3; l++) { /* fields */
				dux[k][l] = _mm256_mul_pd(_mm256_load_pd(du[0][l][i]),dxk[0]);
				dux[k][l] = _mm256_fmadd_pd(_mm256_load_pd(du[1][l][i]),dxk[1],dux[k][l]);
				dux[k][l] = _mm256_fmadd_pd(_mm256_load_pd(du[2][l][i]),dxk[2],dux[k][l]);
			}
		}
		_mm256_store_pd(Du[0],dux[0][0]);
		_mm256_store_pd(Du[1],dux[1][1]);
		_mm256_store_pd(Du[2],dux[2][2]);
		_mm256_store_pd(Du[3],_mm256_mul_pd(mhalf,_mm256_add_pd(dux[0][1],dux[1][0])));
		_mm256_store_pd(Du[4],_mm256_mul_pd(mhalf,_mm256_add_pd(dux[0][2],dux[2][0])));
		_mm256_store_pd(Du[5],_mm256_mul_pd(mhalf,_mm256_add_pd(dux[1][2],dux[2][1])));

		for (e=0; e<NEV; e++) {
			for (k=0; k<6; k++) { /* Stress is coefficient of test function */
				Dv[k][e] = 2 * gaussdata_w[NQP*e+i]* Du[k][e]; 
			}
		}

		dvx[0][0] = _mm256_load_pd(Dv[0]);
		dvx[0][1] = _mm256_load_pd(Dv[3]);
		dvx[0][2] = _mm256_load_pd(Dv[4]);
		dvx[1][0] = _mm256_load_pd(Dv[3]);
		dvx[1][1] = _mm256_load_pd(Dv[1]);
		dvx[1][2] = _mm256_load_pd(Dv[5]);
		dvx[2][0] = _mm256_load_pd(Dv[4]);
		dvx[2][1] = _mm256_load_pd(Dv[5]);
		dvx[2][2] = _mm256_load_pd(Dv[2]);

		for (l=0; l<3; l++) { /* fields  */
			for (k=0; k<3; k++) { /* directions */
				__m256d sum = _mm256_mul_pd(dvx[0][l],_mm256_load_pd(dx[0][k][i]));
				sum = _mm256_fmadd_pd(dvx[1][l],_mm256_load_pd(dx[1][k][i]),sum);
				sum = _mm256_fmadd_pd(dvx[2][l],_mm256_load_pd(dx[2][k][i]),sum);
				_mm256_store_pd(dv[k][l][i],_mm256_mul_pd(mdxdet,sum));
			}
		}
	}
	return 0;
}

typedef struct _p_MFA11SubRepart *MFA11SubRepart;

struct _p_MFA11SubRepart {
  PetscObjectState state;

  PetscScalar       *ufield;
  PetscReal         *LA_gcoords;
  const PetscInt    *elnidx_u;
  PetscScalar       *Yu;

  /* 
   <no postfix> : refers to quantities used in the "normal" case without repartitioning
   _repart      : refers to quantities which are used after repartitioning (these can 
                    be the same as the unrepartitioned values, say if we use the same array)
   _remote      : refers to quantities which are used to communicate with rank 0
   */
  MPI_Comm       comm_sub;
  PetscMPIInt    size_sub;
  PetscInt       nel_sub,nen_u;
  PetscInt       nel,nel_remote,nel_repart;
  PetscInt       nnodes,nnodes_remote,nnodes_repart;
  PetscInt       *nnodes_remote_in,nodes_offset,*nodes_remote,*nel_remote_in;
  PetscInt       *elnidx_u_repart;

#ifdef TATIN_HAVE_CUDA
  MFA11CUDA     cudactx;
#endif
};

#undef __FUNCT__
#define __FUNCT__ "TransferQPData_A11_SubRepart"
static PetscErrorCode TransferQPData_A11_SubRepart(MFA11SubRepart ctx,PetscReal (*wp)[NQP],Quadrature volQ,QPntVolCoefStokes *all_gausspoints,PetscReal *gaussdata_w_remote, PetscReal *gaussdata_w_repart)
{
  PetscErrorCode          ierr;
  PetscInt                i,e;
  PetscMPIInt             rank_sub;
  const QPntVolCoefStokes *cell_gausspoints;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(ctx->comm_sub,&rank_sub);CHKERRQ(ierr);

  /* Create and populate gaussdata_w_remote and _repart. This is the quadrature-point-wise
     information needed to apply the operator. On rank_sub 0, this
     has information for all the elements which will be processed. On  other 
     ranks it has information on the elements which will be offloaded 
     to rank_sub 0. */
  if (rank_sub ) {
    for (e=ctx->nel_repart; e<ctx->nel; e++) {
      PetscInt e_remote = e - ctx->nel_repart; /* element number in set to send */
      ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,(QPntVolCoefStokes**)&cell_gausspoints);CHKERRQ(ierr);
      for (i=0; i<NQP; i++) gaussdata_w_remote[e_remote*NQP + i] = cell_gausspoints[i].eta * (*wp)[i];
    }
  } else {
    for (e=0;               e<ctx->nel; e++) {
      ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,(QPntVolCoefStokes**)&cell_gausspoints);CHKERRQ(ierr);
      for (i=0; i<NQP; i++) gaussdata_w_repart[e       *NQP + i] = cell_gausspoints[i].eta * (*wp)[i];
    }
  }

  /* Send required quadrature-pointwise data to rank_sub 0 */
  if (rank_sub) {
    ierr = MPI_Send(gaussdata_w_remote,NQP*ctx->nel_remote,MPIU_REAL,0,0,ctx->comm_sub);CHKERRQ(ierr);
  } else {
    PetscInt el_offset=ctx->nel;
    for(i=1; i<ctx->size_sub; ++i) {
      ierr = MPI_Recv(&gaussdata_w_repart[NQP * el_offset],NQP * ctx->nel_remote_in[i],MPIU_REAL,i,0,ctx->comm_sub,MPI_STATUS_IGNORE);CHKERRQ(ierr);
      el_offset += ctx->nel_remote_in[i];
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TransferCoordinates_A11_SubRepart"
static PetscErrorCode TransferCoordinates_A11_SubRepart(MFA11SubRepart ctx,const PetscReal *LA_gcoords,PetscReal *LA_gcoords_remote,PetscReal *LA_gcoords_repart)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscMPIInt    rank_sub;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(ctx->comm_sub,&rank_sub);CHKERRQ(ierr);

  /* Copy to arrays to send, or to receive */
  if (rank_sub) {
    for (i=0;i<ctx->nnodes_remote;++i) {
      PetscInt d;
      for(d=0;d<NSD;++d){
        LA_gcoords_remote[NSD*i+d] = LA_gcoords[NSD*ctx->nodes_remote[i]+d];
      }
    }
  } else {
    /* Rank 0  - populate local contributions to arrays we'll receive with */
    ierr = PetscMemcpy(LA_gcoords_repart ,LA_gcoords ,NSD*ctx->nnodes*sizeof(PetscScalar));CHKERRQ(ierr);
  }

  if (rank_sub) {
    ierr = MPI_Send(LA_gcoords_remote,NSD*ctx->nnodes_remote,MPIU_SCALAR,0,0,ctx->comm_sub);CHKERRQ(ierr);
  } else {
    PetscInt nodes_offset = ctx->nnodes;
    for(i=1;i<ctx->size_sub;++i){
      ierr = MPI_Recv(&LA_gcoords_repart[NSD*nodes_offset],NSD*ctx->nnodes_remote_in[i],MPIU_SCALAR,i,0,ctx->comm_sub,MPI_STATUS_IGNORE);CHKERRQ(ierr);
      nodes_offset += ctx->nnodes_remote_in[i];
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TransferUfield_A11_SubRepart"
static PetscErrorCode TransferUfield_A11_SubRepart(MFA11SubRepart ctx,PetscScalar *ufield,PetscScalar *ufield_remote,PetscScalar *ufield_repart)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscMPIInt    rank_sub;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(ctx->comm_sub,&rank_sub);CHKERRQ(ierr);

  /* Copy to arrays to send, or to receive */
  if (rank_sub) {
    /* Ranks 1,2,..  - move data to arrays to be sent */
    for (i=0;i<ctx->nnodes_remote;++i) {
      PetscInt d;
      for(d=0;d<NSD;++d){
        ufield_remote[NSD*i+d] = ufield[NSD*ctx->nodes_remote[i]+d];
      }
    }
  } else {
    /* Rank 0  - populate local contributions to arrays we'll receive with */
    ierr = PetscMemcpy(ufield_repart,ufield,NSD*ctx->nnodes*sizeof(PetscScalar));CHKERRQ(ierr);
  }

  /* Send data to rank_sub 0 */
    if (rank_sub) {
      ierr = MPI_Send(ufield_remote,NSD*ctx->nnodes_remote,MPIU_SCALAR,0,0,ctx->comm_sub);CHKERRQ(ierr);
    } else {
      PetscInt nodes_offset = ctx->nnodes;
      for(i=1;i<ctx->size_sub;++i){
        ierr = MPI_Recv(&ufield_repart[NSD*nodes_offset],NSD*ctx->nnodes_remote_in[i],MPIU_SCALAR,i,0,ctx->comm_sub,MPI_STATUS_IGNORE);CHKERRQ(ierr);
        nodes_offset += ctx->nnodes_remote_in[i];
      }
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TransferYu_A11_SubRepart"
static PetscErrorCode TransferYu_A11_SubRepart(MFA11SubRepart ctx,PetscScalar *Yu,PetscScalar *Yu_remote,PetscScalar *Yu_repart)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscMPIInt    rank_sub;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(ctx->comm_sub,&rank_sub);CHKERRQ(ierr);

  if (rank_sub) {
    ierr = MPI_Recv(Yu_remote,NSD*ctx->nnodes_remote,MPIU_SCALAR,0,0,ctx->comm_sub,MPI_STATUS_IGNORE);CHKERRQ(ierr);
  } else {
    PetscInt nodes_offset = ctx->nnodes;
    for (i=1; i<ctx->size_sub; ++i) {
      ierr = MPI_Send(&Yu_repart[NSD*nodes_offset],NSD*ctx->nnodes_remote_in[i],MPIU_SCALAR,i,0,ctx->comm_sub);CHKERRQ(ierr);
      nodes_offset += ctx->nnodes_remote_in[i];
    }
  }

  /* Accumulate into Yu */
  if (rank_sub) {
    for (i=0; i<ctx->nnodes_remote; ++i) {
      PetscInt d;
      for(d=0;d<NSD;++d){
        Yu[NSD*ctx->nodes_remote[i]+d] += Yu_remote[NSD*i+d];
      }
    }
  } else {
    for (i=0; i<ctx->nnodes; ++i) {
      PetscInt d;
      for(d=0;d<NSD;++d){
        Yu[NSD*i+d] += Yu_repart[NSD*i+d];
      }
    }
  }

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
  PetscMPIInt     rank,rank_sub;
  MPI_Comm        comm_u = PetscObjectComm((PetscObject)dau);
  PetscHashI      nodes_remote_inv;

  PetscFunctionBeginUser;
  if (mf->ctx) PetscFunctionReturn(0);
  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);

  // TODO: set up state

  // TODO: as with CUDA, put this in the apply function but guard it 
  //       with checks on the state. This should make profiling easier, as this stuff is not included in the "setup" log stage!

  /* Define a subcomm. We would hope that this would
     work with MPI_Comm_split_type to split by shared-memory
     domains, but for now we hard-code picking every Nth rank */
  {
    PetscMPIInt size_sub_nominal = 12;
    ierr = PetscOptionsGetInt(NULL,NULL,"-subrepart_size",&size_sub_nominal,NULL);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm_u,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_split(comm_u,rank/size_sub_nominal,0,&ctx->comm_sub);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(ctx->comm_sub,&rank_sub);CHKERRQ(ierr);
    ierr = MPI_Comm_size(ctx->comm_sub,&ctx->size_sub);CHKERRQ(ierr);
  }

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
  } else {
    PetscReal remotefrac = 0.5; /* How much of the work on each rank to offload to rank_sub 0 */
    ierr = PetscOptionsGetReal(NULL,NULL,"-subrepart_frac",&remotefrac,NULL);CHKERRQ(ierr);
    if (remotefrac < 0.0 || remotefrac > 1.0) {
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-subrepart_frac must be in [0,1]");
    }
    const PetscInt nel_repart_small = (PetscInt) ((1.0-remotefrac) * ctx->nel_sub / (ctx->size_sub-1));
    const PetscInt nel_repart_big   = ctx->nel_sub - ((ctx->size_sub-1) * nel_repart_small);
    ctx->nel_repart                 = rank_sub ?  nel_repart_small : nel_repart_big;
  }

  /* (non-zero ranks) Determine the unique nodes corresponding to the elements 
     which will be sent to rank_sub 0. Obtain a count. */
  ctx->nel_remote = rank_sub ? ctx->nel - ctx->nel_repart: 0;

  if (ctx->nel_remote < 0) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"SubRepart's current partitioning scheme does not allow for the current partition. Perhaps you have too few ranks in a sub communicator.");
  }

  /* Overestimate the number of remote nodes, then get the actual number
     of unique nodes by sorting. This could, with the help of a hash table
     be done in linear time, and the wasted memory here could be recovered */
  if (rank_sub) {  
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
    if (!rank_sub) {
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
  if (rank_sub) {
    ctx->nnodes_repart = ctx->nnodes; /* We use the existing array */
  } else {
    ctx->nnodes_repart = ctx->nnodes;
    for (i=1;i<ctx->size_sub;++i) ctx->nnodes_repart += ctx->nnodes_remote_in[i];
  }

  /* Compute an offset for each rank_sub > 0 to process its element
     list to correspond to the numbering on rank_sub 0 */
  {
    PetscInt *nodes_offsets = NULL;
    if (!rank_sub) {
      ierr = PetscMalloc1(ctx->size_sub,&nodes_offsets);
      ctx->nodes_offset=0;
      nodes_offsets[0] = ctx->nodes_offset;
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

  /* Send and receive elements being offloaded to rank 0, defining element-->node 
     maps for the repartition.  On rank_sub > 0, this involves computing the 
     element indices to send to rank 0 by mapping them into the local ordering 
     for the nodes which are used there, and then offsetting.  */
  if (rank_sub) {
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
    for(i=1;i<ctx->size_sub;++i)  {
      ierr = MPI_Recv(&ctx->elnidx_u_repart[ctx->nen_u*el_offset],ctx->nen_u*ctx->nel_remote_in[i],MPIU_INT,i,0,ctx->comm_sub,MPI_STATUS_IGNORE);CHKERRQ(ierr);
      el_offset += ctx->nel_remote_in[i];
    }
  }

#if 1
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
        ierr = PetscPrintf(PETSC_COMM_SELF,"\033[32m[n%d/%d]\033[0m\n",rank_sub,ctx->size_sub);CHKERRQ(ierr);
        if (!rank_sub) {
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
      }
      ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\033[32m========================\033[0m\n");CHKERRQ(ierr);
  }
#endif

#ifdef TATIN_HAVE_CUDA
  /* Set up CUDA context */
  if (!rank_sub) { 
    ierr = PetscMalloc1(1,&ctx->cudactx);CHKERRQ(ierr);
    ierr = MFA11CUDA_SetUp(ctx->cudactx);CHKERRQ(ierr);
  }
#endif

  mf->ctx=ctx;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFA11Destroy_SubRepart"
PetscErrorCode MFA11Destroy_SubRepart(MatA11MF mf)
{
  PetscErrorCode  ierr;
  MFA11SubRepart  ctx;
  PetscMPIInt     rank_sub;

  PetscFunctionBeginUser;
  ctx = (MFA11SubRepart)mf->ctx;
  ierr = MPI_Comm_rank(ctx->comm_sub,&rank_sub);CHKERRQ(ierr);

  if (rank_sub) {
    ierr = PetscFree(ctx->nodes_remote);CHKERRQ(ierr);
  }else {
    ierr = PetscFree(ctx->elnidx_u_repart);CHKERRQ(ierr);
    ierr = PetscFree(ctx->nnodes_remote_in);CHKERRQ(ierr); 
    ierr = PetscFree(ctx->nel_remote_in);CHKERRQ(ierr);
  }

#ifdef TATIN_HAVE_CUDA
  /* Destroy CUDA ctx */
  if (!rank_sub) {
    ierr = MFA11CUDA_CleanUp(ctx->cudactx);CHKERRQ(ierr);
    ierr = PetscFree(ctx->cudactx);CHKERRQ(ierr);
  }
#endif
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
  PetscMPIInt             rank_sub;
  PetscScalar             *ufield_remote,     *ufield_repart;
  PetscScalar             *Yu_remote,         *Yu_repart; 
  PetscReal               *gaussdata_w_remote,*gaussdata_w_repart;
  PetscReal               *LA_gcoords_remote, *LA_gcoords_repart;

  PetscFunctionBeginUser;
  ctx = (MFA11SubRepart)mf->ctx;
  ierr = MPI_Comm_rank(ctx->comm_sub,&rank_sub);

  ierr = PetscLogEventBegin(MAT_MultMFA11_setup,0,0,0,0);CHKERRQ(ierr);
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
  ierr = PetscLogEventEnd(MAT_MultMFA11_setup,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscLogEventBegin(MAT_MultMFA11_rto,0,0,0,0);CHKERRQ(ierr);
  // TODO: only do this for the coords and gauss data if the state has changed,
  //       or if we don't have CUDA (since aren't storing the repart arrays)
  /* Allocate space for repartitioned node-wise fields, and repartitioned elementwise data */
  if (rank_sub) {
    ierr = PetscMalloc1(NSD*ctx->nnodes_remote,&ufield_remote     );CHKERRQ(ierr);
    ierr = PetscMalloc1(NSD*ctx->nnodes_remote,&Yu_remote         );CHKERRQ(ierr); 
    ierr = PetscMalloc1(NSD*ctx->nnodes_remote,&LA_gcoords_remote );CHKERRQ(ierr);
    ierr = PetscMalloc1(NQP*ctx->nel_remote   ,&gaussdata_w_remote);CHKERRQ(ierr);
  } else {
    ierr = PetscMalloc1(NSD*ctx->nnodes_repart,&ufield_repart     );CHKERRQ(ierr);
    ierr = PetscCalloc1(NSD*ctx->nnodes_repart,&Yu_repart         );CHKERRQ(ierr); /* Note Calloc (zeroed) */
    ierr = PetscMalloc1(NSD*ctx->nnodes_repart,&LA_gcoords_repart );CHKERRQ(ierr);
    ierr = PetscMalloc1(NQP*ctx->nel_repart   ,&gaussdata_w_repart);CHKERRQ(ierr);
  }

  /* Send required quadrature-pointwise data to rank_sub 0 */
  ierr = TransferQPData_A11_SubRepart(ctx,&w,volQ,all_gausspoints,gaussdata_w_remote,gaussdata_w_repart);CHKERRQ(ierr);

  /* Send required coordinate data to rank_sub 0 */
  // TODO: only do this if state has changed
  ierr = TransferCoordinates_A11_SubRepart(ctx,LA_gcoords,LA_gcoords_remote,LA_gcoords_repart);CHKERRQ(ierr);

  /* Send required ufield data to rank_sub 0 */
  ierr = TransferUfield_A11_SubRepart(ctx,ufield,ufield_remote,ufield_repart);CHKERRQ(ierr);

  if (rank_sub) {
    ierr = PetscFree(ufield_remote);CHKERRQ(ierr);
    // TODO: only if state has changed (hence we used these)
    ierr = PetscFree(LA_gcoords_remote);CHKERRQ(ierr);
    ierr = PetscFree(gaussdata_w_remote);CHKERRQ(ierr);
  }

  ierr = PetscLogEventEnd(MAT_MultMFA11_rto,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(MAT_MultMFA11_sub,0,0,0,0);CHKERRQ(ierr);

#if defined(_OPENMP)
#define OPENMP_CHKERRQ(x)
#else
#define OPENMP_CHKERRQ(x)   CHKERRQ(x)
#endif
  if (rank_sub) {
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
#if TATIN_HAVE_CUDA
    /* Rank_sub 0 CUDA Implementation */

    ierr = CopyTo_A11_CUDA(mf,ctx->cudactx,ufield_repart,LA_gcoords_repart,gaussdata_w_repart,ctx->nel_repart,ctx->nen_u,ctx->elnidx_u_repart,NSD*ctx->nnodes_repart);CHKERRQ(ierr); 

    // TODO: at this point it'd be safe to free gaussdata_w_repart and LA_gcoords_repart, as they wouldn't be needed until the state changed

    ierr = ProcessElements_A11_CUDA(ctx->cudactx,ctx->nen_u,NSD*ctx->nnodes_repart);CHKERRQ(ierr);

    ierr = CopyFrom_A11_CUDA(ctx->cudactx,Yu_repart,NSD*ctx->nnodes_repart);CHKERRQ(ierr);

#else
    /* Rank_sub 0 AVX Implementation */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\033[31m!!!!!!!!\nWARNING: using AVX implementation on rank_sub 0, because CUDA isn't available: THIS WILL NEVER PERFORM BETTER THAN THE AVX IMPLEMENTATION\n!!!!!!!!\033[0m\n");CHKERRQ(ierr);
#if defined(_OPENMP)
#pragma omp parallel for private(i)
#endif
    for (e=0;e<ctx->nel_repart;e+=NEV) {
      PetscScalar elu[3][Q2_NODES_PER_EL_3D][NEV]={},elx[3][Q2_NODES_PER_EL_3D][NEV]={},elv[3][Q2_NODES_PER_EL_3D][NEV];
      PetscScalar dx[3][3][NQP][NEV],dxdet[NQP][NEV],du[3][3][NQP][NEV],dv[3][3][NQP][NEV];
      PetscInt    ee,l;
      PetscReal   gaussdata_w_local[NQP*NEV];

      for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
        for (ee=0; ee<NEV; ee++) {
          PetscInt E = ctx->elnidx_u_repart[ctx->nen_u*PetscMin(e+ee,ctx->nel_repart-1)+i]; /* Pad up to length NEV by duplicating last element */
          for (l=0; l<3; l++) {
            elx[l][i][ee] = LA_gcoords_repart[3*E+l];
            elu[l][i][ee] = ufield_repart[3*E+l];
          }
        }
      }

      for (ee=0; ee<NEV; ++ee){
        const PetscInt el = PetscMin(e + ee,ctx->nel_repart-1);
        for (i=0;i<NQP;++i){
          gaussdata_w_local[ee*NQP+i] = gaussdata_w_repart[el*NQP+i];
        }
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

      ierr = QuadratureAction_A11AVX_mod(gaussdata_w_local,dx,dxdet,du,dv);OPENMP_CHKERRQ(ierr);

      ierr = PetscMemzero(elv,sizeof elv);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(D,B,B,GRAD_TRANSPOSE,dv[0],elv);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(B,D,B,GRAD_TRANSPOSE,dv[1],elv);OPENMP_CHKERRQ(ierr);
      ierr = TensorContractNEV_AVX(B,B,D,GRAD_TRANSPOSE,dv[2],elv);OPENMP_CHKERRQ(ierr);

#if defined(_OPENMP)
#pragma omp critical
#endif
      for (ee=0; ee<PetscMin(NEV,ctx->nel_repart-e); ee++) {
        for (i=0; i<NQP; i++) {
          PetscInt E = ctx->elnidx_u_repart[ctx->nen_u*(e+ee)+i];
          for (l=0; l<3; l++) {
            Yu_repart[3*E+l] += elv[l][i][ee];
          }
        }
      }
    }
#endif
  }
#undef OPENMP_CHKERRQ
  ierr = PetscLogEventEnd(MAT_MultMFA11_sub,0,0,0,0);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr); 

  /* Transfer and accumulate contributions to Yu from rank_sub 0 */
  ierr = PetscLogEventBegin(MAT_MultMFA11_rfr,0,0,0,0);CHKERRQ(ierr);
  ierr = TransferYu_A11_SubRepart(ctx,Yu,Yu_remote,Yu_repart);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(MAT_MultMFA11_rfr,0,0,0,0);CHKERRQ(ierr);

  // Outdated. TODO update once AVX and CUDA guts are in
#if 0
  PetscLogFlops((ctx->nel * 9) * 3*NQP*(6+6+6));           /* 9 tensor contractions per element */
  PetscLogFlops(ctx->nel*NQP*(14 + 1/* division */ + 27)); /* 1 Jacobi inversion per element */
  PetscLogFlops(ctx->nel*NQP*(5*9+6+6+6*9));               /* 1 quadrature action per element */
#endif

  if (rank_sub) {
    ierr = PetscFree(Yu_remote);CHKERRQ(ierr);
  } else {
    ierr = PetscFree(ufield_repart);CHKERRQ(ierr);
    ierr = PetscFree(Yu_repart);CHKERRQ(ierr);
    ierr = PetscFree(LA_gcoords_repart);CHKERRQ(ierr); // TODO: only destroy if state changed
    ierr = PetscFree(gaussdata_w_repart);CHKERRQ(ierr); // TODO: only destroy if state changed
  }

  PetscFunctionReturn(0);
}

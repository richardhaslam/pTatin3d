/*
   A prototype of a matrix-free method designed to use available CUDA devices.

   Work is shared between CUDA and AVX implementations, assuming that
   one may launch CUDA data transfers and kernels from each MPI rank.
 */

#if defined(__AVX__) && defined(TATIN_HAVE_CUDA)

#include <petscfe.h>
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

#define PTATIN_AVXCUDA_CPU_FRAC_OPT "-avxcuda_cpu_frac"

typedef struct _p_MFA11AVXCUDA *MFA11AVXCUDA;

struct _p_MFA11AVXCUDA {
  PetscObjectState state;
  // TODO
  PetscInt  nel,nel_cpu,nel_gpu,nen_u;
  PetscInt  *elnidx_u;
  MFA11CUDA cudactx;
};

PetscErrorCode MFA11SetUp_AVXCUDA(MatA11MF mf)
{
  PetscErrorCode ierr;
  MFA11AVXCUDA   ctx;

  PetscFunctionBeginUser;
  if (mf->ctx) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(MAT_MultMFA11_SUP,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);

  ctx->state = 0;

  ierr = DMDAGetElements_pTatinQ2P1(dau,&ctx->nel,&ctx->nen_u,&ctx->elnidx_u);CHKERRQ(ierr);

  /* Determine element partition between the two kernels */
  {
    PetscBool set;
    PetscReal cpufrac = 0.25;
    ierr = PetscOptionsGetReal(NULL,NULL,PTATIN_AVXCUDA_CPU_FRAC_OPT,&cpufrac,&set);CHKERRQ(ierr);
    if (!set) {
      ierr = PetscPrintf(PetscObjectComm((PetscObject)mf->dmUVW),"Warning: %s not provided, so default value of %f used (likely to be unbalanced!)\n",PTATIN_AVXCUDA_CPU_FRAC_OPT,cpufrac);CHKERRQ(ierr);
    }
    if (cpufrac < 0.0 || cpufrac > 1.0) {
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"%s must be in [0,1]",PTATIN_AVXCUDA_CPU_FRAC_OPT);
    }
    // TODO consider forcing to a layer boundary
    ctx->nel_cpu = (PetscInt) (cpufrac * ctx->nel);
    ctx->nel_gpu = ctx->nel - ctx->nel_cpu;
  }
  // TODO compute an upper bound on the dof #s for cpu, and a lower bound for gpu

  /* Set up CUDA context */
  ierr = PetscMalloc1(1,&ctx->cudactx);CHKERRQ(ierr);
  ierr = MFA11CUDA_SetUp(ctx->cudactx);CHKERRQ(ierr);

  mf->ctx=ctx;
  ierr = PetscLogEventEnd(MAT_MultMFA11_SUP,0,0,0,0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode MFA11Destroy_AVXCUDA(MatA11MF mf)
{
  PetscErrorCode ierr;
  MFA11AVXCUDA   ctx;

  PetscFunctionBeginUser;
  ctx = (MFA11AVXCUDA)mf->ctx;

  /* Destroy CUDA ctx */
  ierr = MFA11CUDA_CleanUp(ctx->cudactx);CHKERRQ(ierr);
  ierr = PetscFree(ctx->cudactx);CHKERRQ(ierr);

  ierr = PetscFree(mf->ctx);CHKERRQ(ierr);
  mf->ctx = NULL;

  PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A11_AVXCUDA(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{
  PetscErrorCode        ierr;
  MFA11AVXCUDA          ctx;

  PetscFunctionBeginUser;
  ctx = (MFA11AVXCUDA)mf->ctx;

  // TODO

  // TODO placeholder: just 2 copies of the AVX kernel (test logic, set up passing tests)

  PetscFunctionReturn(0);
}

#endif /* defined(__AVX__) && defined(TATIN_HAVE_CUDA) */

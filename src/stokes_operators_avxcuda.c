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

typedef struct _p_MFA11AVXCUDA *MFA11AVXCUDA;

struct _p_MFA11AVXCUDA {
  PetscObjectState state;
  // TODO
  MFA11CUDA        cudactx;
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

  // TODO

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

  // TODO

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

  PetscFunctionReturn(0);
}

#endif /* defined(__AVX__) && defined(TATIN_HAVE_CUDA) */

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
  PetscInt        nel,nel_cpu,nel_gpu,nen_u;
  const PetscInt  *elnidx_u;
  MFA11CUDA       cudactx;
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

  ierr = DMDAGetElements_pTatinQ2P1(mf->daUVW,&ctx->nel,&ctx->nen_u,&ctx->elnidx_u);CHKERRQ(ierr);

  /* Determine element partition between the two kernels */
  {
    PetscBool set;
    PetscReal cpufrac = 0.25;
    ierr = PetscOptionsGetReal(NULL,NULL,PTATIN_AVXCUDA_CPU_FRAC_OPT,&cpufrac,&set);CHKERRQ(ierr);
    if (!set) {
      //ierr = PetscPrintf(PetscObjectComm((PetscObject)mf->daUVW),"Warning: %s not provided, so default value of %f used (likely to be unbalanced!)\n",PTATIN_AVXCUDA_CPU_FRAC_OPT,cpufrac);CHKERRQ(ierr);
      SETERRQ2(PetscObjectComm((PetscObject)mf->daUVW),PETSC_ERR_SUP,"%s not provided, so default value of %f used (likely to be unbalanced!)\n",PTATIN_AVXCUDA_CPU_FRAC_OPT,cpufrac);CHKERRQ(ierr); // TODO convert back to warning, but for now we don't want this printing to mess with profiling
    }
    if (cpufrac < 0.0 || cpufrac > 1.0) {
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"%s must be in [0,1]",PTATIN_AVXCUDA_CPU_FRAC_OPT);
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

  // TODO this is just the AVX kernel data
  DM cda;
  Vec gcoords;
  const PetscReal *LA_gcoords;
  PetscInt nel,nen_u;
  const PetscInt *elnidx_u;
  QPntVolCoefStokes *all_gausspoints;
  PetscReal x1[3],w1[3],B[3][3],D[3][3],w[NQP];
  PetscInt i,j,k,e;

  PetscFunctionBeginUser;
  ctx = (MFA11AVXCUDA)mf->ctx;

  // TODO

  // TODO placeholder: just 2 copies of the AVX kernel (test logic, set up passing tests)

  // TODO placeholder: this is just the AVX kernel
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

  ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);

  ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);

#if defined(_OPENMP)
  #define OPENMP_CHKERRQ(x)
#else
  #define OPENMP_CHKERRQ(x)   CHKERRQ(x)
#endif

#if defined(_OPENMP)
    #pragma omp parallel for private(i)
#endif
  for (e=0;e<nel;e+=NEV) {
    PetscScalar elu[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32,elx[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32,elv[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32;
    PetscScalar dx[3][3][NQP][NEV] ALIGN32,dxdet[NQP][NEV],du[3][3][NQP][NEV] ALIGN32,dv[3][3][NQP][NEV] ALIGN32;
    const QPntVolCoefStokes *cell_gausspoints[NEV];
    PetscInt ee,l;

    for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
      for (ee=0; ee<NEV; ee++) {
        PetscInt E = elnidx_u[nen_u*PetscMin(e+ee,nel-1)+i]; // Pad up to length NEV by duplicating last element
        for (l=0; l<3; l++) {
          elx[l][i][ee] = LA_gcoords[3*E+l];
          elu[l][i][ee] = ufield[3*E+l];
        }
      }
    }
    for (ee=0; ee<NEV; ee++) {
      ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,PetscMin(e+ee,nel-1),(QPntVolCoefStokes**)&cell_gausspoints[ee]);OPENMP_CHKERRQ(ierr);
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
    for (ee=0; ee<PetscMin(NEV,nel-e); ee++) {
      for (i=0; i<NQP; i++) {
        PetscInt E = elnidx_u[nen_u*(e+ee)+i];
        for (l=0; l<3; l++) {
          Yu[3*E+l] += elv[l][i][ee];
        }
      }
    }

  }

#undef OPENMP_CHKERRQ

  PetscLogFlops((nel * 9) * 3*NQP*(6+6+6));           /* 9 tensor contractions per element */
  PetscLogFlops(nel*NQP*(14 + 1/* division */ + 27)); /* 1 Jacobi inversion per element */
  PetscLogFlops(nel*NQP*(5*9+6+6+6*9));               /* 1 quadrature action per element */

  ierr = VecRestoreArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif /* defined(__AVX__) && defined(TATIN_HAVE_CUDA) */

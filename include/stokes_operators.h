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
 **    filename:   stokes_operators.h
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


#ifndef __ptatin_stokes_operators_h__
#define __ptatin_stokes_operators_h__

typedef struct _p_MatStokesMF *MatStokesMF;
typedef struct _p_MatA11MF *MatA11MF;

struct _p_MatStokesMF {
  PetscInt    mu,mp,Mu,Mp;
  DM          stokes_pack,daUVW,dap;
  BCList      u_bclist,p_bclist;
  Quadrature  volQ;
  DM          daU;
  IS          isUVW,isU,isV,isW,isP; /* Need the IS's for GetSubMatrix */
  PetscInt    refcnt;
};

struct _p_MatA11MF {
  PetscInt    mu,Mu;
  DM          daUVW;
  BCList      u_bclist;
  Quadrature  volQ;
  DM          daU; /* Optionally need this */
  IS          isUVW; /* Needed for full column space */
  IS          isU,isV,isW; /* Optionally: Need the IS's for GetSubMatrix */
  PetscInt    refcnt;
  PetscObjectState state;
  PetscBool   is_setup;
  void        *ctx; /* special context for cuda, opencl, openmp spmv implementations */
  PetscErrorCode (*SpMVOp_MatMult)(MatA11MF,Quadrature,DM,PetscScalar[],PetscScalar[]);
  PetscErrorCode (*SpMVOp_SetUp)(MatA11MF);
  PetscErrorCode (*SpMVOp_Destroy)(MatA11MF);
};

PetscErrorCode StokesQ2P1CreateMatrix_Operator(PhysCompStokes user,Mat *B);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_QuasiNewtonX(PhysCompStokes user,Mat *B);

PetscErrorCode StokesQ2P1CreateMatrixNest_Operator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B);

PetscErrorCode StokesQ2P1CreateMatrixNest_PCOperator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B);

PetscErrorCode MatStokesMFCreate(MatStokesMF *B);
PetscErrorCode MatA11MFCreate(MatA11MF *B);
PetscErrorCode MatStokesMFSetup(MatStokesMF StkCtx,PhysCompStokes user);
PetscErrorCode MatA11MFSetup(MatA11MF A11Ctx,DM dav,Quadrature volQ,BCList u_bclist);
PetscErrorCode MatStokesMFDestroy(MatStokesMF *B);
PetscErrorCode MatA11MFDestroy(MatA11MF *B);
PetscErrorCode MatStokesMFCopy(MatStokesMF A,MatStokesMF *B);
PetscErrorCode MatA11MFCopy(MatA11MF A,MatA11MF *B);
PetscErrorCode MatCopy_StokesMF_A11MF(MatStokesMF A,MatA11MF *B);

PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11(MatA11MF A11,Mat *A);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11LowOrder(MatA11MF A11,Mat *A);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A12(MatStokesMF Stk,Mat *A12);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A21(MatStokesMF Stk,Mat *A21);

PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A_QuasiNewtonX(MatStokesMF Stk,Mat *A11);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11_QuasiNewtonX(MatA11MF A11,Mat *A);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A12_QuasiNewtonX(MatStokesMF Stk,Mat *A12);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A21_QuasiNewtonX(MatStokesMF Stk,Mat *A21);

/* non-optimal preallocation routines for A12,A21 */
PetscErrorCode StokesQ2P1CreateMatrix_A12(PhysCompStokes user,Mat *mat);
PetscErrorCode StokesQ2P1CreateMatrix_A21(PhysCompStokes user,Mat *mat);

PetscErrorCode MatShellGetMatStokesMF(Mat A,MatStokesMF *mf);
PetscErrorCode MatShellGetMatA11MF(Mat A,MatA11MF *mf);

PetscErrorCode MatMultAdd_basic(Mat A,Vec v1,Vec v2,Vec v3);
PetscErrorCode MatMultTransposeAdd_generic(Mat mat,Vec v1,Vec v2,Vec v3);

/* Implementation-specific structures and routines shared amongst multiple implementations */
#define NQP 27  /* Number of quadrature points per element; must equal Q2_NODES_PER_EL_3D (27) */
#define NEV 4   /* Number of elements over which to vectorize */

typedef enum {
  GRAD,
  GRAD_TRANSPOSE
} GradMode;

PetscErrorCode TensorContractNEV_AVX(PetscReal Rf[][3],PetscReal Sf[][3],PetscReal Tf[][3],GradMode gmode,PetscReal x[][NQP][NEV],PetscReal y[][NQP][NEV]);
__attribute__((noinline))
PetscErrorCode JacobianInvertNEV_AVX(PetscScalar dx[3][3][NQP][NEV],PetscScalar dxdet[NQP][NEV]);
__attribute__((noinline))
PetscErrorCode QuadratureAction_A11_AVX(const QPntVolCoefStokes *gausspt[],
             PetscScalar dx[3][3][Q2_NODES_PER_EL_3D][NEV],
             PetscScalar dxdet[Q2_NODES_PER_EL_3D][NEV],
             PetscReal w[Q2_NODES_PER_EL_3D],
             PetscScalar du[3][3][Q2_NODES_PER_EL_3D][NEV],
             PetscScalar dv[3][3][Q2_NODES_PER_EL_3D][NEV]);


typedef struct _p_MFA11CUDA *MFA11CUDA;

struct _p_MFA11CUDA {
  PetscObjectState state;

  PetscScalar *ufield;
  PetscReal   *LA_gcoords;
  PetscReal   *gaussdata_w;  // Data at Gauss points multiplied by respective quadrature weight
  PetscInt    element_colors;
  PetscInt    *elements_per_color;
  PetscInt    **el_ids_colored;
  PetscInt    *elnidx_u;
  PetscScalar *Yu;
};

#ifdef __cplusplus
extern "C" {
#endif
PetscErrorCode MFA11CUDA_SetUp(MFA11CUDA);
PetscErrorCode MFA11CUDA_CleanUp(MFA11CUDA);
PetscErrorCode CopyTo_A11_CUDA(MatA11MF,MFA11CUDA,const PetscScalar*,const PetscReal*,const PetscReal*,PetscInt,PetscInt,const PetscInt*,PetscInt);
PetscErrorCode ProcessElements_A11_CUDA(MFA11CUDA,PetscInt,PetscInt);
PetscErrorCode CopyFrom_A11_CUDA(MFA11CUDA,PetscScalar*,PetscInt);
#ifdef __cplusplus
}
#endif

#endif


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
	PetscInt    level;
	PetscInt    ii;
	DM          stokes_pack,daUVW,dap;
	BCList      u_bclist,p_bclist;
	Quadrature  volQ;
	DM          daU;
	IS          isUVW,isU,isV,isW,isP; /* Need the IS's for GetSubMatrix */
	PetscInt    refcnt;
	Vec         nodal_viscosity;
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
	/* Not sure I need this at all */
	PetscInt    level;
	PetscInt    ii;
	Vec         nodal_viscosity;
	PetscErrorCode (*Wrapper_A11)(Quadrature,DM,PetscScalar[],PetscScalar[]);
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

#endif


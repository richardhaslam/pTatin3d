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
 **    filename:   stokes_operators_mf.h
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

#ifndef __ptatin3d_stokes_operators_mf_h__
#define __ptatin3d_stokes_operators_mf_h__

#include "stokes_operators.h"

/* MatMult wrappers */
PetscErrorCode MFStokesWrapper_A11(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A11_Tensor(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A11_AVX(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A11_OpenCL(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFA11SetUp_OpenCL(MatA11MF mf);
PetscErrorCode MFA11Destroy_OpenCL(MatA11MF mf);
PetscErrorCode MFStokesWrapper_A11PC(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A11PC_1x1x1(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A11PC_2x2x2(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A11_UPX(Quadrature volQ,DM dau,PetscScalar ufield[],DM dax,PetscScalar xfield[],PetscScalar Yu[]);

PetscErrorCode MFStokesWrapper_A(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],PetscScalar Yu[],PetscScalar Yp[]);
PetscErrorCode MFStokesWrapper_A_AVX(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],PetscScalar Yu[],PetscScalar Yp[]);
PetscErrorCode MFStokesWrapper_A_UPX(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],DM dax,PetscScalar xfield[],PetscScalar Yu[],PetscScalar Yp[]);

PetscErrorCode MFStokesWrapper_A12(Quadrature volQ,DM dau,DM dap,PetscScalar Xp[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A12_AVX(Quadrature volQ,DM dau,DM dap,PetscScalar pfield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A12_UPX(Quadrature volQ,DM dau,DM dap,PetscScalar pfield[],DM dax,PetscScalar xfield[],PetscScalar Yu[]);

PetscErrorCode MFStokesWrapper_A21(Quadrature volQ,DM dau,DM dap,PetscScalar Xu[],PetscScalar Yp[]);
PetscErrorCode MFStokesWrapper_A21_UPX(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,DM dax,PetscScalar xfield[],PetscScalar Yp[]);

#ifdef __cplusplus
extern "C" {
#endif
PetscErrorCode MFStokesWrapper_A11_CUDA(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFA11SetUp_CUDA(MatA11MF mf);
PetscErrorCode MFA11Destroy_CUDA(MatA11MF mf);
#ifdef __cplusplus
}
#endif

PetscErrorCode MFStokesWrapper_A11_SubRepart(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFA11SetUp_SubRepart(MatA11MF mf);
PetscErrorCode MFA11Destroy_SubRepart(MatA11MF mf);
PetscErrorCode MFStokesWrapper_A11_AVXCUDA(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFA11SetUp_AVXCUDA(MatA11MF mf);
PetscErrorCode MFA11Destroy_AVXCUDA(MatA11MF mf);

/* MatGetDiagonal wrappers */
PetscErrorCode MFStokesWrapper_diagA11(Quadrature volQ,DM dau,PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_diagA11LowOrder(Quadrature volQ,DM dau,PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_diagA11_UPX(Quadrature volQ,DM dau,DM dax,PetscScalar xfield[],PetscScalar Yu[]);

#endif


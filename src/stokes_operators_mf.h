/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      stokes_operators_mf.h
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

#ifndef __ptatin3d_stokes_operators_mf_h__
#define __ptatin3d_stokes_operators_mf_h__

/* MatMult wrappers */
PetscErrorCode MFStokesWrapper_A11(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A11PC(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);

PetscErrorCode MFStokesWrapper_A(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],PetscScalar Yu[],PetscScalar Yp[]);

PetscErrorCode MFStokesWrapper_A12(Quadrature volQ,DM dau,DM dap,PetscScalar Xp[],PetscScalar Yu[]);

PetscErrorCode MFStokesWrapper_A21(Quadrature volQ,DM dau,DM dap,PetscScalar Xu[],PetscScalar Yp[]);

/* MatGetDiagonal wrappers */
PetscErrorCode MFStokesWrapper_diagA11(Quadrature volQ,DM dau,PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_diagA11LowOrder(Quadrature volQ,DM dau,PetscScalar Yu[]);

#endif

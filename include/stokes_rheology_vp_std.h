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
 **    filename:   stokes_rheology_vp_std.h
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

#ifndef __ptatin_stokes_rheology_vp_std_h__
#define __ptatin_stokes_rheology_vp_std_h__

PetscErrorCode EvaluateRheologyNonlinearitiesMarkers_VPSTD(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[]);
PetscErrorCode ApplyViscosityCutOffMarkers_VPSTD(pTatinCtx user);
PetscErrorCode StokesCoefficient_UpdateTimeDependentQuantities_VPSTD(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[]);

#endif


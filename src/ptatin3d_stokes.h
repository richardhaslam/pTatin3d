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
 **    Filename:      ptatin3d_stokes.h
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

#ifndef __ptatin3d_stokes_h__
#define __ptatin3d_stokes_h__

//#include "petsc.h"
//#include "petscdm.h"
//#include "dmda_bcs.h"
#include "QPntVolCoefStokes_def.h"

PetscErrorCode StokesVelocity_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);
PetscErrorCode StokesPressure_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);
PetscErrorCode StokesVelocityScalar_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);

PetscErrorCode PhysCompCreate_Stokes(PhysCompStokes *ctx);
PetscErrorCode PhysCompDestroy_Stokes(PhysCompStokes *ctx);

PetscErrorCode PhysCompCreateMesh_Stokes3d(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx);
PetscErrorCode PhysCompCreateBoundaryList_Stokes(PhysCompStokes ctx);
PetscErrorCode PhysCompCreateVolumeQuadrature_Stokes(PhysCompStokes ctx);

PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[]);
PetscErrorCode DMDASetValuesLocalStencil_InsertValues_Stokes_Velocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[]);
PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(PetscScalar *fields_F,PetscInt p_eqn[],PetscScalar Fe_p[]);
PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_ScalarVelocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[]);

PetscErrorCode PhysCompSaveMesh_Stokes3d(PhysCompStokes ctx,const char fname_vel[],const char fname_p[],const char fname_coors[]);
PetscErrorCode PhysCompLoadMesh_Stokes3d(PhysCompStokes ctx,const char fname_vel[],const char fname_p[]);

PetscErrorCode pTatinStokesKSPMonitorBlocks(KSP ksp,PetscInt n,PetscReal rnorm,void *data);

PetscErrorCode VolumeQuadratureCreate_GaussLegendreStokes(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature);
PetscErrorCode VolumeQuadratureGetAllCellData_Stokes(Quadrature Q,QPntVolCoefStokes *coeffs[]);
PetscErrorCode VolumeQuadratureGetCellData_Stokes(Quadrature Q,QPntVolCoefStokes coeffs[],PetscInt cidx,QPntVolCoefStokes *cell[]);

#endif


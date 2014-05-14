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
#include "QPntSurfCoefStokes_def.h"

PetscErrorCode StokesVelocity_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);
PetscErrorCode StokesPressure_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);
PetscErrorCode StokesVelocityScalar_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);

PetscErrorCode PhysCompCreate_Stokes(PhysCompStokes *ctx);
PetscErrorCode PhysCompDestroy_Stokes(PhysCompStokes *ctx);

PetscErrorCode PhysCompCreateMesh_Stokes3d(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx);
PetscErrorCode PhysCompCreateBoundaryList_Stokes(PhysCompStokes ctx);
PetscErrorCode PhysCompCreateVolumeQuadrature_Stokes(PhysCompStokes ctx);
PetscErrorCode PhysCompCreateSurfaceQuadrature_Stokes(PhysCompStokes ctx);

PetscErrorCode PhysCompStokesSetGravityUnitVector(PhysCompStokes ctx,PetscReal grav[]);
PetscErrorCode PhysCompStokesScaleGravityVector(PhysCompStokes ctx,PetscReal fac);
PetscErrorCode PhysCompStokesSetGravityVector(PhysCompStokes ctx,PetscReal grav[]);

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

PetscErrorCode SurfaceQuadratureCreate_GaussLegendreStokes(DM da,HexElementFace index,SurfaceQuadrature *quadrature);
PetscErrorCode SurfaceQuadratureGeometrySetUpStokes(SurfaceQuadrature Q,DM da);
PetscErrorCode SurfaceQuadratureOrientationSetUpStokes(SurfaceQuadrature Q,DM da);
PetscErrorCode SurfaceQuadratureGetAllCellData_Stokes(SurfaceQuadrature Q,QPntSurfCoefStokes *coeffs[]);
PetscErrorCode SurfaceQuadratureGetCellData_Stokes(SurfaceQuadrature Q,QPntSurfCoefStokes coeffs[],PetscInt cidx,QPntSurfCoefStokes *cell[]);

PetscErrorCode SurfaceQuadratureOrientationViewGnuplotStokes(SurfaceQuadrature Q,DM da,const char name[]);

PetscErrorCode SNESStokes_ConvergenceTest_UPstol(SNES snes,PetscInt it,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason,void *dummy);
PetscErrorCode SNESStokes_SetConvergenceTest_UPstol(SNES snes,pTatinCtx user);

PetscErrorCode PetscOptionsInsertPrefixString(const char prefix[],const char option[]);
PetscErrorCode SNESStokesPCSetOptions_A(SNES snes);
PetscErrorCode SNESStokesPCMGSetOptions(SNES snes,PetscInt maxits,PetscBool mglog);
PetscErrorCode SNESStokesPCMGCoarseSetOptions_SparseDirect(SNES snes);
PetscErrorCode SNESStokesPCMGCoarseSetOptions_IterativeASM(SNES snes,PetscReal rtol,PetscInt maxits,PetscInt overlap);
PetscErrorCode SNESStokesPCMGCoarseSetOptions_NestedIterativeASM(SNES snes,PetscReal rtol,PetscInt maxits,PetscInt maxitsnested,PetscInt overlap);

PetscErrorCode PhysCompStokesGetDMComposite(PhysCompStokes stokes,DM *dmc);
PetscErrorCode PhysCompStokesGetDMs(PhysCompStokes stokes,DM *dmv,DM *dmp);
PetscErrorCode PhysCompStokesGetBCList(PhysCompStokes stokes,BCList *ulist,BCList *plist);
PetscErrorCode PhysCompStokesGetVolumeQuadrature(PhysCompStokes stokes,Quadrature *q);
PetscErrorCode PhysCompStokesGetVolumeQuadratureAllCellData(PhysCompStokes stokes,QPntVolCoefStokes *coeffs[]);
PetscErrorCode PhysCompStokesGetSurfaceQuadrature(PhysCompStokes stokes,HexElementFace fid,SurfaceQuadrature *sq);
PetscErrorCode PhysCompStokesGetSurfaceQuadratureAllCellData(PhysCompStokes stokes,HexElementFace fid,QPntSurfCoefStokes *coeffs[]);

#endif


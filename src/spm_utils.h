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
 **    Filename:      spm_utils.h
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

#ifndef __spm_utils_h__
#define __spm_utils_h__

#include "petsc.h"
#include "petscdm.h"

/* MPI to SEQ */
PetscErrorCode DMDAGatherIKRedundantSurfaceDMDA(DM dm_mech,DM *dm_msurf0);
PetscErrorCode InterpolateMSurf0ToSPMSurfIKGrid(DM dm_msurf0,PetscInt spm_mi,PetscInt spm_mj,PetscReal *spm_coords,PetscReal *spm_H);
PetscErrorCode InterpolateSPMSurfIKGridToMSurf0(PetscInt spm_mi,PetscInt spm_mj,PetscReal *spm_coords,PetscReal *spm_H,DM dm_msurf0);
PetscErrorCode DMDAScatterIKRedundantSurfaceDMDA(DM dm_msurf0,DM dm_mech);

/* MPI to MPI (semi-redundant) */
/* TODO */

/* MPI to MPI */
PetscErrorCode DMDAGatherIKSurfaceDMDA(DM dm_mech,DM *_dm_msurf,Vec *_elevation);
PetscErrorCode DMDAScatterIKSurfaceDMDA(DM dm_msurf,Vec height,DM dm_mech);

#endif

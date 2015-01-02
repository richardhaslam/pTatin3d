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
 **    filename:   dmda_remesh.h
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
#ifndef __PTATIN3d_DMDA_REMESH_H__
#define __PTATIN3d_DMDA_REMESH_H__

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

PetscErrorCode DMDARemeshSetUniformCoordinatesInPlane_IJ(DM da,PetscInt K,DMDACoor3d coords[]);
PetscErrorCode DMDARemeshSetUniformCoordinatesInPlane_IK(DM da,PetscInt J,DMDACoor3d coords[]);

PetscErrorCode DMDAGetCornerCoordinatesInPlane_IJ(DM da,PetscInt K,DMDACoor3d coords[]);

PetscErrorCode DMDARemeshSetUniformCoordinatesBetweenJLayers3d( DM da, PetscInt startJ, PetscInt endJ );
PetscErrorCode DMDARemeshSetUniformCoordinatesBetweenKLayers3d( DM da, PetscInt startK, PetscInt endK );
PetscErrorCode DMDARemeshJMAX_UpdateHeightsFromInterior(DM da);
PetscErrorCode DMDASetCoordinatesColumnRefinement(DM da,PetscInt dir,PetscReal factor,PetscReal x1_frac,PetscReal x2_frac);

#endif


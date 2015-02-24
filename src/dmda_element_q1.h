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
 **    filename:   dmda_element_q1.h
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

#ifndef __ptatin_dmda_element_q1_h__
#define __ptatin_dmda_element_q1_h__

#include <petsc.h>
#include <petscdm.h>

PetscErrorCode DMDASetElementType_Q1(DM da);
PetscErrorCode DMDAGetElements_DA_Q1_3D(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx);

PetscErrorCode DMDACreateOverlappingQ1FromQ2(DM dmq2,PetscInt ndofs,DM *dmq1);
PetscErrorCode DMDACreateNestedQ1FromQ2(DM dmq2,PetscInt ndofs,DM *dmq1);

PetscErrorCode DMDAEQ1_GetElementCoordinates_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[]);
PetscErrorCode DMDAEQ1_GetScalarElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDAEQ1_GetVectorElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDAEQ1_SetValuesLocalStencil_AddValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[]);
PetscErrorCode DMDAEQ1_GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[]);

PetscErrorCode DMDAGetElementsQ1(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx);

PetscErrorCode DMDAProjectVectorQ2toOverlappingQ1_3d(DM daq2,Vec x2,DM daq1,Vec x1);
PetscErrorCode DMDAProjectVectorQ2toNestedQ1_3d(DM daq2,Vec x2,DM daq1,Vec x1);
PetscErrorCode DMDAProjectVectorQ2toQ1(DM daq2,Vec x2,DM daq1,Vec x1,PetscInt mesh_type);
PetscErrorCode DMDAProjectCoordinatesQ2toQ1(DM daq2,DM daq1,PetscInt mesh_type);

PetscErrorCode DMDAGetLocalSizeElementQ1(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz);
PetscErrorCode DMDAGetCornersElementQ1(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz);
PetscErrorCode DMDACreateQ1(MPI_Comm comm,PetscInt MX,PetscInt MY,PetscInt MZ,PetscInt ndofs,DM *dmq1);

#endif


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
 **    filename:   dmda_element_q1macrop1.h
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

#ifndef __ptatin_dmda_element_q1macrop1_h__
#define __ptatin_dmda_element_q1macrop1_h__

#include <petsc.h>
#include <petscdm.h>

PetscErrorCode DMDASetElementType_Q1Macro(DM da);
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ);
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetLocalSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ);
PetscErrorCode DMDAEQ1Macro_NaturalSpace_GetSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ);
PetscErrorCode DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ);
PetscErrorCode DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ);
PetscErrorCode DMDAEGetElements_Q1MacroMixedSpace(DM da,PetscInt *nel,PetscInt *nen,PetscInt *e[]);
PetscErrorCode DMDAEGetElements_Q1MacroNaturalSpace(DM da,PetscInt *nel,PetscInt *nen,PetscInt *e[]);
PetscErrorCode DMDAEGetElementMap_Q1MacroNaturalToMixedSpace(DM da,PetscInt *nen,PetscInt *e[]);
PetscErrorCode DMDAEGetElementMap_Q1MacroNaturalToMixedLocalSpace(DM da,PetscInt *nen,PetscInt *e[]);
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetOwnershipRangesElement(DM da,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt **si,PetscInt **sj,PetscInt **sk,PetscInt **_mx,PetscInt **_my,PetscInt **_mz);

PetscErrorCode DMDAESetType_Q1Macro(DM da);

PetscErrorCode DMDAEQ1Macro_GetElementCoordinates_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[]);
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetElementCoordinates_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[]);
PetscErrorCode DMDAEQ1Macro_GetScalarElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetScalarElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDAEQ1Macro_GetVectorElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetVectorElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDAEQ1Macro_SetValuesLocalStencil_AddValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[]);
PetscErrorCode DMDAEQ1Macro_MixedSpace_SetValuesLocalStencil_AddValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[]);
PetscErrorCode DMDAEQ1Macro_GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[]);
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[]);


#endif


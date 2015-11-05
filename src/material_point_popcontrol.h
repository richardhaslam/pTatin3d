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
 **    filename:   material_point_popcontrol.h
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

#ifndef __ptatin_material_point_popcontrol_h__
#define __ptatin_material_point_popcontrol_h__

typedef struct _p_PSortCtx {
    PetscInt point_index;
    PetscInt cell_index;
} PSortCtx;

PetscErrorCode MaterialPointPopulationControl_v1(pTatinCtx ctx);
PetscErrorCode MaterialPointRegionAssignment_v1(DataBucket db,DM da);
PetscErrorCode MaterialPointRegionAssignment_v2(DataBucket db,DM da);

PetscErrorCode MPPCCreateSortedCtx(DataBucket db,DM da,PetscInt *_np,PetscInt *_nc,PSortCtx **_plist,PetscInt **_pcell_list);
PetscErrorCode MPPCDestroySortedCtx(DataBucket db,DM da,PSortCtx **_plist,PetscInt **_pcell_list);
PetscErrorCode MPPCSortedCtxGetNumberOfPointsPerCell(DataBucket db,PetscInt cell_idx,PetscInt pcell_list[],PetscInt *np);
PetscErrorCode MPPCSortedCtxGetPointByCell(DataBucket db,PetscInt cell_idx,PetscInt pidx,PSortCtx plist[],PetscInt pcell_list[],MPntStd **point);

#endif

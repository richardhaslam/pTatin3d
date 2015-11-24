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
 **    filename:   output_material_points_p0.h
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

#ifndef __ptatin_output_material_points_p0_h__
#define __ptatin_output_material_points_p0_h__

PetscErrorCode MarkerCellFieldsP0Write_ParaView(DM pack,DataBucket material_points,const char basename[],
                                                const int nvars,const MaterialPointVariable vars[],
                                                PetscBool low_precision,
                                                const char path[],const char prefix[]);

PetscErrorCode pTatin3dModelOutput_MarkerCellFieldsP0_ParaView(pTatinCtx ctx,const int nvars,const MaterialPointVariable vars[],PetscBool low_precision,const char prefix[]);
PetscErrorCode pTatin3dModelOutput_MarkerCellFieldsP0_PetscVec(pTatinCtx ctx,PetscBool dm_velocity_data_required,const int nvars,const MaterialPointVariable vars[],PetscBool low_precision,const char prefix[]);

#endif

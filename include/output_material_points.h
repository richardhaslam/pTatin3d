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
 **    filename:   output_material_points.h
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

#ifndef __ptatin_output_material_points_h__
#define __ptatin_output_material_points_h__

typedef enum { 
	MPV_region=0, 
	MPV_viscosity, 
	MPV_density, 
	MPV_plastic_strain, 
	MPV_yield_indicator, 
	MPV_diffusivity, 
	MPV_heat_source 
} MaterialPointVariable;

PetscErrorCode pTatinOutputParaViewMarkerFields_VTS(DM dau,DataBucket material_points,const int nvars,const MaterialPointVariable vars[],const char name[]);
PetscErrorCode pTatinOutputParaViewMarkerFields_PVTS(DM dau,const int nvars,const MaterialPointVariable vars[],const char prefix[],const char name[]);
PetscErrorCode pTatinOutputParaViewMarkerFields(DM pack,DataBucket material_points,const int nvars,const MaterialPointVariable vars[],const char path[],const char prefix[]);
PetscErrorCode pTatin3d_ModelOutput_MarkerCellFields(pTatinCtx ctx,const int nvars,const MaterialPointVariable vars[],const char prefix[]);

#endif

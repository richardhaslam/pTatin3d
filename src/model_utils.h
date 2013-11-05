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
 **    Filename:      model_utils.h
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

#ifndef __ptatin3d_model_utils_h__
#define __ptatin3d_model_utils_h__

#define GRAVITY 9.8

PetscErrorCode MPntGetField_global_element_IJKindex(DM da, MPntStd *material_point, PetscInt *I, PetscInt *J, PetscInt *K);
PetscErrorCode pTatinModelGetOptionReal(const char option[],PetscReal *val,const char error[],const char default_opt[],PetscBool essential);
PetscReal absolute(PetscReal a);
PetscErrorCode detrend(PetscReal array[],PetscInt n);
PetscErrorCode rednoise(PetscReal rnoise[],PetscInt n,PetscInt seed);

typedef struct {
        PetscScalar nlayers;
        PetscReal lscale;
        PetscReal cond[20];
        PetscReal hp[20];
        PetscReal qbase[20];
        PetscReal ytop[20];
        PetscReal ttop[20];
        PetscReal thick[20];
} DMDA_thermalfield_init_params;

PetscBool DMDAVecTraverse_InitialThermalField3D(PetscScalar pos[],PetscScalar *val,void *ctx);
PetscErrorCode DMDAConvertGlobalElementIndex2IJK(DM da,PetscInt localeid,PetscInt *I,PetscInt *J,PetscInt *K);

#endif

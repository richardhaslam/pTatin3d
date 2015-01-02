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
 **    filename:   material_point_point_location.h
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


#ifndef __ptatin3d_MATERIAL_POINT_POINT_LOCATION_H__
#define __ptatin3d_MATERIAL_POINT_POINT_LOCATION_H__

void InverseMappingDomain_2dQ2(PetscReal tolerance,PetscInt max_its,
                               PetscBool use_nonzero_guess,
                               PetscBool monitor,
                               const PetscReal coords[],const PetscInt mx,const PetscInt my,const PetscInt element[],
                               PetscInt np,MPntStd marker[]);

void InverseMappingDomain_3dQ2(PetscReal tolerance, PetscInt max_its,
                               PetscBool use_nonzero_guess,
                               PetscBool monitor,
                               const PetscReal coords[],const PetscInt mx,const PetscInt my,const PetscInt mz,const PetscInt element[],
                               int np,MPntStd marker[]);

#endif


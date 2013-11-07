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
 **    Filename:      iplus_ctx.h
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
 **    $Id:$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#ifndef __ptatin3d_model_iplus_ctx_h__
#define __ptatin3d_model_iplus_ctx_h__

typedef enum { iPLUsModelA=0 } iPLUSModelType;
typedef enum { iPLUSMatMantle = 0, iPLUSMatPlume, iPLUSMatSlab } iPLUSMaterialType;

/* define user model */
typedef struct {
	iPLUSModelType modeltype;
	PetscReal slab_eta,slab_rho;
	PetscReal mantle_eta,mantle_rho;
	PetscReal plume_eta,plume_rho;
	PetscReal plume_pos[3];
	PetscReal plume_A0;
	PetscReal plume_radius;
	PetscInt  refinement_type;
	PetscInt  nplume_elements,*plume_element;
	PetscReal intial_domain_volume;
	PetscViewer logviewer;
	PetscInt    np_plume_x,np_plume_z;
} iPLUSCtx;

#endif

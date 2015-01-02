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
 **    filename:   iplus_ctx.h
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


#ifndef __ptatin3d_model_iplus_ctx_h__
#define __ptatin3d_model_iplus_ctx_h__

#include "geometry_object.h"

typedef enum { iPLUsModelSlab = 0, iPLUsModelPlume, iPLUsModelSlabPlume } iPLUSModelType;
typedef enum { iPLUSMatMantle = 0, iPLUSMatPlume, iPLUSMatSlab } iPLUSMaterialType;

/* define user model */
typedef struct {
	iPLUSModelType modeltype;
	PetscReal   slab_eta,slab_rho;
	PetscReal   mantle_eta,mantle_rho;
	PetscReal   plume_eta,plume_rho;
	PetscReal   plume_pos[3];
	PetscReal   plume_A0;
	PetscReal   plume_radius;
	PetscInt    refinement_type;
	PetscInt    nplume_elements,*plume_element;
	PetscReal   intial_domain_volume;
	PetscViewer logviewer;
	PetscInt    np_plume_x,np_plume_z;
	GeometryObject slab_geometry;
	PetscInt    iplus_output_frequency;
	PetscReal   eta_scale,vel_scale,time_scale,length_scale;
} iPLUSCtx;


/* slab prototypes */
PetscErrorCode iPLUS_DefineSlabMaterial(DM dav,DataBucket materialpoint_db,iPLUSCtx *data);
PetscErrorCode iPLUS_CreateSlabGeometry(iPLUSCtx *data);

/* plume prototypes */
PetscErrorCode iPLUS_DetermineElementsContainingPlumeInlet(DM dav,iPLUSCtx *data);
PetscErrorCode iPLUS_InsertPlumeMaterial(DM dav,DataBucket materialpoint_db,iPLUSCtx *data);
PetscBool iPLUS_Inflow_BCListEvaluator(PetscScalar position[],PetscScalar *value,void *ctx);
PetscErrorCode iPLUS_ApplyMaterialBoundaryCondition_Plume(pTatinCtx c,iPLUSCtx *data);

#endif

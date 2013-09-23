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
 **    Filename:      ptatin3d.h
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


#ifndef __ptatin3d_h__
#define __ptatin3d_h__

#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"



typedef struct _p_pTatinCtx *pTatinCtx;
typedef struct _p_pTatinModel *pTatinModel;
typedef struct _p_PhysCompStokes *PhysCompStokes;
typedef struct _p_PhysCompEnergy *PhysCompEnergy;
typedef struct _p_Quadrature *Quadrature;
typedef struct _p_SurfaceQuadrature *SurfaceQuadrature;

typedef enum { LINE_QUAD=0,SURFACE_QUAD,VOLUME_QUAD } QuadratureType;

typedef struct _p_RheologyConstants RheologyConstants;


/*
#include "dmda_bcs.h"
#include "dmda_checkpoint.h"
#include "dmda_compare.h"
#include "dmda_duplicate.h"
#include "dmda_element_q2p1.h"
#include "dmda_project_coords.h"
#include "dmda_redundant.h"
#include "dmda_remesh.h"
#include "dmda_update_coords.h"
#include "dmda_view_petscvtk.h"

#include "swarm_fields.h"
#include "data_exchanger.h"

#include "ptatin3d_defs.h"
#include "ptatin_utils.h"
#include "ptatin3d_stokes.h"
#include "ptatin_models.h"
#include "rheology.h"
*/

#include "swarm_fields.h"
#include "data_exchanger.h"
//#include "rheology.h"
#include "material_point_load.h"
#include "material_point_utils.h"


PetscErrorCode pTatin3d_PhysCompStokesCreate(pTatinCtx user);
PetscErrorCode pTatin3d_ModelOutput_VelocityPressure_Stokes(pTatinCtx ctx,Vec X,const char prefix[]);
PetscErrorCode pTatin3d_ModelOutputLite_Velocity_Stokes(pTatinCtx ctx,Vec X,const char prefix[]);
PetscErrorCode pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(pTatinCtx ctx,Vec X,const char prefix[]);

PetscErrorCode pTatin3dCreateMaterialPoints(pTatinCtx ctx,DM dav);
PetscErrorCode MaterialPointCoordinateSetUp(pTatinCtx ctx,DM da);
PetscErrorCode pTatin3d_ModelOutput_MPntStd(pTatinCtx ctx,const char prefix[]);

PetscErrorCode pTatin3dCreateContext(pTatinCtx *ctx);
PetscErrorCode pTatin3dDestroyContext(pTatinCtx *ctx);
PetscErrorCode pTatin3dParseOptions(pTatinCtx ctx);
PetscErrorCode pTatinModelLoad(pTatinCtx ctx);

PetscErrorCode pTatinGetTimestep(pTatinCtx ctx,PetscReal *dt);
PetscErrorCode pTatinGetMaterialPoints(pTatinCtx ctx,DataBucket *db,DataEx *de);
PetscErrorCode pTatinGetModel(pTatinCtx ctx,pTatinModel *m);
PetscErrorCode pTatinGetRheology(pTatinCtx ctx,RheologyConstants **r);
PetscErrorCode pTatinGetStokesContext(pTatinCtx ctx,PhysCompStokes *s);
PetscErrorCode pTatinGetMaterialConstants(pTatinCtx ctx,DataBucket *db);

PetscErrorCode pTatin3dContextLoad(pTatinCtx *ctx,const char filename[]);
PetscErrorCode pTatin3dContextSave(pTatinCtx ctx,const char filename[]);
PetscErrorCode pTatin3dCheckpoint(pTatinCtx ctx,Vec X,const char prefix[]);

PetscErrorCode pTatin3d_PhysCompStokesLoad(pTatinCtx user,const char vname[],const char pname[]);
PetscErrorCode pTatin3dRestart(pTatinCtx ctx);

PetscErrorCode pTatinCtxGetModelData(pTatinCtx ctx,const char name[],void **data);
PetscErrorCode pTatinCtxAttachModelData(pTatinCtx ctx,const char name[],void *data);

PetscErrorCode pTatin3dCheckpointManager(pTatinCtx ctx,Vec X);
PetscErrorCode DMCoarsenHierarchy2_DA(DM da,PetscInt nlevels,DM dac[]);

PetscErrorCode pTatin_SetTimestep(pTatinCtx ctx,const char timescale_name[],PetscReal dt_trial);

#endif

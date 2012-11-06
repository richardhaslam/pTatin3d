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
 **    Filename:      ptatin3d_ctx_impl.h
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


#ifndef __private_ptatin3d_ctx_impl_h__
#define __private_ptatin3d_ctx_impl_h__

#include "petsc.h"
#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "data_exchanger.h"
#include "ptatin3d_stokes.h"
#include "ptatin_models.h"
#include "rheology.h"
#include "material_constants.h"

struct _p_pTatinCtx {
	PhysCompStokes stokes_ctx;
	//	PhysCompEnergy energy_ctx;
	//	PhysCompCoords coords_ctx;

	PetscBool  restart_from_file;
	char       restart_prefix[PETSC_MAX_PATH_LEN];
	char			 restart_dir[PETSC_MAX_PATH_LEN];
	PetscInt   checkpoint_every, checkpoint_every_nsteps;
	PetscReal  checkpoint_every_ncpumins;
	PetscBool  use_mf_stokes;
	
	/* rheology */
  //RheologyConstants rheology_constants;
	
	/* Mesh size */
	PetscInt   mx,my,mz;
	
  DM  pack; /* all physics gets jammed in here */
	
	char       formatted_timestamp[PETSC_MAX_PATH_LEN];
	char       outputpath[PETSC_MAX_PATH_LEN];
	/* material points */
	PetscInt   coefficient_projection_type;
	DataBucket materialpoint_db;
	DataEx     materialpoint_ex;	
	/* options */
	PetscBool solverstatistics;
	/* snes continuation paramter */
	PetscInt continuation_m, continuation_M;
	/* time stepping */
	PetscInt  nsteps,step;
	PetscReal dt,dt_max,dt_min,dt_adv;
	PetscInt  output_frequency;
	PetscReal time_max,time;
	
	/* rheology */
  RheologyConstants rheology_constants;
	DataBucket material_constants;

	/* model function pointers */
	pTatinModel model;
	PetscContainer model_data;
	/* logger */
	PetscViewer log;
};

#endif



#ifndef __private_ptatin3d_ctx_impl_h__
#define __private_ptatin3d_ctx_impl_h__

#include "petsc.h"
#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "data_exchanger.h"
#include "ptatin3d_stokes.h"
#include "ptatin_models.h"
#include "rheology.h"

struct _p_pTatinCtx {
	PhysCompStokes stokes_ctx;
	//	PhysCompEnergy energy_ctx;
	//	PhysCompCoords coords_ctx;

	PetscBool  restart_from_file;
	char       restart_prefix[1256];
	PetscBool  use_mf_stokes;
	
	/* rheology */
  //RheologyConstants rheology_constants;
	
	/* Mesh size */
	PetscInt   mx,my,mz;
	
  DM  pack; /* all physics gets jammed in here */
	
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

	/* model function pointers */
	pTatinModel model;
	PetscContainer model_data;
};

#endif

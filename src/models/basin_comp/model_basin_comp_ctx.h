
#ifndef __ptatin3d_model_basin_comp_ctx_h__
#define __ptatin3d_model_basin_comp_ctx_h__

/* define user model */
typedef struct {
	PetscInt  max_layers;
	PetscInt  n_interfaces;
	PetscReal interface_heights[3];
    PetscInt  layer_res_j[2];
	PetscReal eta[2];
	PetscReal rho[2];
	PetscInt  bc_type;
	PetscReal exx;
	PetscReal vx_commpression;
	PetscReal Lx, Lz, Ly, alpha;
} ModelBasinCompCtx;

#endif

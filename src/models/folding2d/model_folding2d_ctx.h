
#ifndef __ptatin3d_model_folding2d_ctx_h__
#define __ptatin3d_model_folding2d_ctx_h__

/* define user model */
typedef struct {
	PetscInt  max_layers;
	PetscInt  n_interfaces;
	PetscReal interface_heights[100];
	PetscInt  layer_res_j[100];
	PetscReal eta[100];
	PetscReal rho[100];
	PetscInt  bc_type;
	PetscReal exx;
	PetscReal vx_commpression;
} ModelFolding2dCtx;

#endif

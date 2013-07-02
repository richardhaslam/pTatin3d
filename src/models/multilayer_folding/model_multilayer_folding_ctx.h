#ifndef __ptatin3d_model_multilayer_folding_ctx_h__
#define __ptatin3d_model_multilayer_folding_ctx_h__

/* define user model */
typedef struct {
	PetscInt  max_layers;
	PetscInt  n_interfaces;
	PetscReal interface_heights[100];
    PetscInt  layer_res_k[99];
	PetscReal eta[100];
	PetscReal rho[100];
	PetscInt  bc_type, perturbation_type;
	PetscReal exx;
	PetscReal vx_commpression;
	PetscReal Lx, Lz, Ly;
    PetscReal amp;
} ModelMultilayerFoldingCtx;

#endif


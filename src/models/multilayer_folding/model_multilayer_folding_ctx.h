#ifndef __ptatin3d_model_multilayer_folding_ctx_h__
#define __ptatin3d_model_multilayer_folding_ctx_h__

/* define user model */
typedef struct {
	PetscInt  max_layers;
	PetscInt  n_interfaces;
	PetscReal interface_heights[100];
    PetscInt  layer_res_j[99];
	PetscReal eta[100];
	PetscReal rho[100];
	PetscInt  bc_type, perturbation_type;
	PetscReal ezz; 
	PetscReal exx;
	PetscReal vx_compression;
	PetscReal vz_compression;
	PetscReal Lx, Lz, Ly, L_char;
    PetscReal amp;
    PetscReal kx,kz; 
    PetscReal A0;
} ModelMultilayerFoldingCtx;

#endif



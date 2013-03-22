
#ifndef __ptatin3d_model_wrench_fold_ctx_h__
#define __ptatin3d_model_wrench_fold_ctx_h__

/* define user model */
typedef struct {
	PetscInt  max_layers;
	PetscInt  n_interfaces;
	PetscReal interface_heights[101];
        PetscInt  layer_res_k[100];
	PetscReal eta[100];
	PetscReal rho[100];
	PetscInt  bc_type;
	PetscReal exy, A;
	PetscReal Lx, Lz, Ly;
} ModelWrenchFoldCtx;

#endif

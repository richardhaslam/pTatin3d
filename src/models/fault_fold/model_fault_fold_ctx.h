
#ifndef __ptatin3d_model_faul_fold_ctx_h__
#define __ptatin3d_model_faul_fold_ctx_h__

/* define user model */
typedef struct {
	PetscInt  max_layers;
	PetscInt  n_interfaces;
	PetscReal interface_heights[101];
    PetscInt  layer_res_k[100];
	PetscReal eta[100];
	PetscReal etaweak[100];
	PetscReal rho[100];
    PetscReal fold_centers[4];
    PetscReal sigma;
	PetscInt  bc_type;
	PetscReal exx, vx_commpression;
	PetscReal Lx, Lz, Ly;
    PetscReal amp;
    PetscReal phi;
    PetscReal psig1;
    PetscReal psig2;
    PetscReal offset;
    PetscInt perturbation_type, perturbation_width;
     
} ModelFaultFoldCtx;

#endif


#ifndef __ptatin3d_model_faul_fold_plastic_ctx_h__
#define __ptatin3d_model_faul_fold_plastic_ctx_h__

/* define user model */
typedef struct {
	PetscInt  max_layers;
    PetscInt  max_fault_layers;
	PetscInt  n_interfaces;
	PetscReal interface_heights[101];
    PetscInt  layer_res_k[100];
	PetscReal eta[100];
	PetscReal C0[100];
    PetscReal mu[100];
	PetscReal C0_inf[100];
    PetscReal mu_inf[100];
	PetscReal rho[100];
    
    PetscInt  domain_res_j[3];    
 	PetscReal fault_eta[1];
	PetscReal fault_C0[1];
    PetscReal fault_mu[1];
	PetscReal fault_C0_inf[1];
    PetscReal fault_mu_inf[1];
	PetscReal fault_rho[1]; 
    
    PetscReal fold_centers[4];
    PetscReal sigma;
	PetscInt  bc_type;
	PetscReal exx, vx_commpression;
	PetscReal Lx, Lz, Ly;
    PetscReal amp;
} ModelFaultFoldPlasticCtx;

#endif

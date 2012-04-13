

#ifndef __ptatinmodel_indendor_ctx_h__
#define __ptatinmodel_indentor_ctx_h__

typedef enum { VSBC_FreeSlip=0, VSBC_NoSlip, VSBC_FreeSlipFreeSurface, VSBC_NoSlipFreeSurface } ModelBC;

typedef struct {
	PetscInt  nmaterials;
	PetscReal Lx,Ly,Lz;
	PetscReal eta[20];
	PetscReal rho[20];
	ModelBC  boundary_conditon_type; /* [ 0 free slip | 1 no slip | 2 free surface + free slip | 3 free surface + no slip ] */
	PetscReal cutoff_time;
	PetscReal indentation_velocity;
	
	PetscBool dimensional;
	PetscReal density_bar;
	PetscReal length_bar;
	PetscReal viscosity_bar;
	PetscReal velocity_bar;
	PetscReal time_bar;
	PetscReal pressure_bar;
} ModelIndentorCtx;


#endif

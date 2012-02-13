

#ifndef __ptatinmodel_viscous_sinker_ctx_h__
#define __ptatinmodel_viscous_sinker_ctx_h__

typedef enum { VSBC_FreeSlip=0, VSBC_NoSlip, VSBC_FreeSlipFreeSurface, VSBC_NoSlipFreeSurface } ViscousSinkerBC;

typedef struct {
	PetscInt  nmaterials;
	PetscReal eta0,eta1;
	PetscReal rho0,rho1;
	PetscBool is_sphere;
	PetscReal Lx,Ly,Lz;
	PetscReal origin[3];
	PetscReal length[3];
	ViscousSinkerBC  boundary_conditon_type; /* [ 0 free slip | 1 no slip | 2 free surface + free slip | 3 free surface + no slip ] */
} ModelViscousSinkerCtx;


#endif

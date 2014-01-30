
#ifndef __ptatin3d_model_convection2d_ctx_h__
#define __ptatin3d_model_convection2d_ctx_h__

/* define user model */
typedef struct {
	PetscReal Lx,Ly,Lz,Ox,Oy,Oz;
	PetscReal velocity;
	PetscReal eta[1],rho[1];
	PetscReal C0[1],mu[1],C0_inf[1],mu_inf[1];
	PetscReal diffusivity[1],alpha[1],H[1];
	PetscReal Ttop,Tbot;
	PetscReal density_bar,length_bar,viscosity_bar,velocity_bar,time_bar,pressure_bar;
} ModelThermal_Convection2dCtx;

#endif

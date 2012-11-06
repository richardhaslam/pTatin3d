

#ifndef __ptatinmodel_folding_ctx_h__
#define __ptatinmodel_folding_ctx_h__

typedef enum { FOLD_MODEL_A=0 } FoldModel;
/* 
  FOLD_MODEL_A => Benchmark from ""

*/

typedef struct {
	PetscInt  nmaterials;
	PetscReal Lx,Ly,Lz;
	PetscReal eta[20];
	PetscReal rho[20];
	FoldModel  model_type; 
	PetscReal cutoff_time;
	
	PetscBool dimensional;
	PetscReal density_bar;
	PetscReal length_bar;
	PetscReal viscosity_bar;
	PetscReal velocity_bar;
	PetscReal time_bar;
	PetscReal pressure_bar;
} ModelFoldingCtx;


#endif


#ifndef __ptatin_external_packages_h__
#define __ptatin_external_packages_h__

PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_SPMA(pTatinCtx pctx,Vec X);


PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X,
																																	PetscInt refinement_factor,
																																	PetscReal Lstar,PetscReal Vstar,
																																	PetscReal dt_mechanical,
																																	int _law,double _m,double _kf,double _kd,int _bc);


#endif



#ifndef __ptatin_external_packages_h__
#define __ptatin_external_packages_h__

PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_SPMA(pTatinCtx pctx,Vec X);


PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X,
																																	 PetscInt refinement_factor,
																																	 PetscReal Lstar,PetscReal Vstar,
																																	 PetscReal dt_mechanical,PetscReal dt_spm,
																																	 PetscInt _law,PetscReal _m,PetscReal _kf,PetscReal _kd,PetscInt _bc);


#endif


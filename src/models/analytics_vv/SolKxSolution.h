
#ifndef __SolKxSolution_h__
#define __SolKxSolution_h__

#include <petsc.h>

typedef struct {
  PetscReal m;
  PetscInt n;
  PetscReal B;
} ParamsSolKx;

PetscErrorCode SolKxSolution(const PetscReal pos[], PetscReal m, PetscInt n, PetscReal B,
                             PetscScalar vel[], PetscScalar *p, PetscScalar s[], PetscScalar gamma[], PetscScalar *nu);

PetscErrorCode EvaluateV_SolKx(PetscReal pos[],PetscReal vel[],void *ctx);
PetscErrorCode EvaluateE_SolKx(PetscReal pos[],PetscReal E[],void *ctx);
PetscErrorCode EvaluateP_SolKx(PetscReal pos[],PetscReal p[],void *ctx);
PetscErrorCode SolKxParamsSetValues(ParamsSolKx *data,PetscReal m,PetscInt n,PetscReal B);

#endif

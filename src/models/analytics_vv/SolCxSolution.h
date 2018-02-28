
#ifndef __SolCxSolution_h__
#define __SolCxSolution_h__

#include <petsc.h>

typedef struct {
  PetscReal etaA,etaB,xc;
  PetscInt n;
} ParamsSolCx;

PetscErrorCode SolCxSolution(const PetscReal pos[], PetscReal etaA, PetscReal etaB, PetscReal xc,
                             PetscInt n,
                             PetscScalar vel[], PetscScalar *p, PetscScalar s[], PetscScalar gamma[]);

PetscErrorCode EvaluateV_SolCx(PetscReal pos[],PetscReal vel[],void *ctx);
PetscErrorCode EvaluateE_SolCx(PetscReal pos[],PetscReal E[],void *ctx);
PetscErrorCode EvaluateP_SolCx(PetscReal pos[],PetscReal p[],void *ctx);
PetscErrorCode SolCxParamsSetValues(ParamsSolCx *data,PetscReal etaA,PetscReal etaB,PetscReal xc,PetscInt n);

#endif

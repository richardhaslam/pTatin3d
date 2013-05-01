

#ifndef __ptatin_eigen_operators_h__
#define __ptatin_eigen_operators_h__

PetscErrorCode MatCreateEigenOperatorFromKSPOperators(KSP ksp,Mat *A);
PetscErrorCode MatCreateEigenOperatorFromKSP(KSP ksp,Mat *A);

#endif

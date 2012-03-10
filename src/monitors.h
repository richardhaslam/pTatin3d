
#ifndef __ptatin_monitors_h__
#define __ptatin_monitors_h__

PetscErrorCode pTatin_KSPMonitor_StdoutStokesResiduals3d(KSP ksp,PetscInt n,PetscReal rnorm,void *data);
PetscErrorCode pTatin_KSPMonitor_ParaviewStokesResiduals3d(KSP ksp,PetscInt n,PetscReal rnorm,void *data);

#endif

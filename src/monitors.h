
#ifndef __ptatin_monitors_h__
#define __ptatin_monitors_h__

PetscErrorCode pTatin_KSPMonitorStokesResiduals3d(KSP ksp,PetscInt n,PetscReal rnorm,void *data);
PetscErrorCode pTatin_KSPMonitor_ViewStokesResiduals3d(KSP ksp,PetscInt n,PetscReal rnorm,void *data);

#endif

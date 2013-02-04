
#ifndef __ptatin3d_energy_h__
#define __ptatin3d_energy_h__

PetscErrorCode pTatinGetContext_Energy(pTatinCtx ctx,PhysCompEnergy *e);
PetscErrorCode pTatinContextValid_Energy(pTatinCtx ctx,PetscBool *exists);
PetscErrorCode pTatinPhysCompCreate_Energy(pTatinCtx user);
PetscErrorCode pTatinPhysCompActivate_Energy(pTatinCtx user,PetscBool load);

#endif

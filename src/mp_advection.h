
#ifndef __ptatin_mp_advection_h__
#define __ptatin_mp_advection_h__

PetscErrorCode SwarmUpdatePosition_MPntStd_Euler(DM da,Vec velocity,PetscReal step,int npoints,MPntStd marker[]);
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep(DM da,Vec velocity,PetscReal *step);
PetscErrorCode SwarmUpdateProperties_MPntStd(DataBucket db,pTatinCtx ctx,Vec X);
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates(DataBucket materialpoints,DM dav,Vec velocity,PetscReal dt);
PetscErrorCode MaterialPointStd_UpdateLocalCoordinates(DataBucket materialpoints,DM dav);
PetscErrorCode MaterialPointStd_Removal(DataBucket materialpoints);
PetscErrorCode SwarmUpdatePosition_Communication_Generic(DataBucket db,DM da,DataEx de);
PetscErrorCode MaterialPointStd_UpdateCoordinates(DataBucket materialpoints,DM dav,DataEx de);


#endif

#ifndef __PTATIN3d_DMDA_CHECKPOINT_H__
#define __PTATIN3d_DMDA_CHECKPOINT_H__

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>


PetscErrorCode DMDAPackDataToFile(DM da,const char name[]);
PetscErrorCode DMDACreateFromPackDataToFile(MPI_Comm comm,const char name[],DM *da);
PetscErrorCode DMDALoadGlobalVectorFromFile(DM da,const char name[],Vec *da_x);
PetscErrorCode DMDALoadCoordinatesFromFile(DM da,const char name[]);


#endif

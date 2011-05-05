
#ifndef __PTATIN3d_DMDA_COMPARE_H__
#define __PTATIN3d_DMDA_COMPARE_H__

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

PetscErrorCode DMDACompareStructures(DM da1,DM da2,PetscBool *flg);

PetscErrorCode DMDA_CheckNodeIndex3d(DM da,PetscBool ghosted,PetscInt i,PetscInt j,PetscInt k );
PetscErrorCode DMDA_CheckNodeIndex2d(DM da,PetscBool ghosted,PetscInt i,PetscInt j,PetscInt k );

#endif



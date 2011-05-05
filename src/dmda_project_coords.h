
#ifndef __PTATIN3d_DMDA_PROJECT_COORDS_H__
#define __PTATIN3d_DMDA_PROJECT_COORDS_H__

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>


PetscErrorCode DMDARestrictCoordinateHierarchy(DM da[],PetscInt nlevels);
PetscErrorCode DMDARestrictCoordinates(DM daf,DM dac);

#endif

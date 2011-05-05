


#ifndef __PTATIN3d_DMAD_UPDATE_COORDS_H__
#define __PTATIN3d_DMAD_UPDATE_COORDS_H__


#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>


PetscErrorCode DMDAUpdateGhostedCoordinates(DM da);
PetscErrorCode DMDASetCoordinatesFromLocalVector(DM da,Vec local_coords);
PetscErrorCode DMDASetCoordinatesU(DM da,Vec coords);
PetscErrorCode DMDACloneCoordinates(DM da,DM da_clone);

#endif

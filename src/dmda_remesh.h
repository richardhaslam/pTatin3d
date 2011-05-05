#ifndef __PTATIN3d_DMDA_REMESH_H__
#define __PTATIN3d_DMDA_REMESH_H__

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

PetscErrorCode DMDARemeshSetUniformCoordinatesInPlane_IJ(DM da,PetscInt K,DMDACoor3d coords[]);

PetscErrorCode DMDAGetCornerCoordinatesInPlane_IJ(DM da,PetscInt K,DMDACoor3d coords[]);
PetscErrorCode DMDARemeshSetUniformCoordinatesBetweenKLayers3d( DM da, PetscInt startK, PetscInt endK );


#endif

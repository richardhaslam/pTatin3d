#ifndef __PTATIN3d_DMDA_VIEW_PETSCVTK_H__
#define __PTATIN3d_DMDA_VIEW_PETSCVTK_H__

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

PetscErrorCode DMDAViewPetscVTK(DM da,Vec field,const char name[]);

#endif

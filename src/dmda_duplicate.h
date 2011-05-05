#ifndef __PTATIN3d_DMDA_DUPLICATE_H__
#define __PTATIN3d_DMDA_DUPLICATE_H__

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

PetscErrorCode DMDADuplicateLayout(DM da1,PetscInt dof2,PetscInt sw2,DMDAStencilType st2,DM *da2);

#endif

#ifndef __PTATIN3d_DMDA_REDUNDANT_H__
#define __PTATIN3d_DMDA_REDUNDANT_H__

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

PetscErrorCode DMDACreate3dRedundant(DM da,PetscInt si, PetscInt ei, PetscInt sj, PetscInt ej, PetscInt sk, PetscInt ek, PetscInt n_dofs, DM *_seq_DM );


#endif

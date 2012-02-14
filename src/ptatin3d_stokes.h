
#ifndef __ptatin3d_stokes_h__
#define __ptatin3d_stokes_h__

//#include "petsc.h"
//#include "petscdm.h"
//#include "dmda_bcs.h"

PetscErrorCode StokesVelocity_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);
PetscErrorCode StokesPressure_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);

PetscErrorCode PhysCompCreate_Stokes(PhysCompStokes *ctx);
PetscErrorCode PhysCompDestroy_Stokes(PhysCompStokes *ctx);

PetscErrorCode PhysCompCreateMesh_Stokes3d(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx);
PetscErrorCode PhysCompCreateBoundaryList_Stokes(PhysCompStokes ctx);
PetscErrorCode PhysCompCreateVolumeQuadrature_Stokes(PhysCompStokes ctx);

PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[]);
PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(PetscScalar *fields_F,PetscInt p_eqn[],PetscScalar Fe_p[]);


#endif

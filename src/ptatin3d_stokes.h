
#ifndef __ptatin3d_stokes_h__
#define __ptatin3d_stokes_h__

#include "petsc.h"
#include "petscdm.h"
#include "dmda_bcs.h"


typedef struct _p_PhysCompStokes *PhysCompStokes;

struct _p_PhysCompStokes {
	PetscInt                mx,my,mz; /* global mesh size */
	DM                      dav,dap;
  DM                      stokes_pack;
	BCList                  u_bclist,p_bclist;
//	QuadratureStokes        Q;
//	SurfaceQuadratureStokes surfQ[QUAD_EDGES]; /* four edges */
	PetscBool               use_mf_stokes;
};

PetscErrorCode StokesVelocity_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);
PetscErrorCode StokesPressure_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[]);

PetscErrorCode PhysCompCreate_Stokes(PhysCompStokes *ctx);
PetscErrorCode PhysCompDestroy_Stokes(PhysCompStokes *ctx);

PetscErrorCode PhysCompCreateMesh_Stokes3d(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx);
PetscErrorCode PhysCompCreateBoundaryList_Stokes(PhysCompStokes ctx);


#endif
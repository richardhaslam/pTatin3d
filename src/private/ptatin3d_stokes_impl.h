
#ifndef __private_ptatin3d_stokes_impl_h__
#define __private_ptatin3d_stokes_impl_h__

#include "petsc.h"
#include "petscdm.h"
#include "dmda_bcs.h"

struct _p_PhysCompStokes {
	PetscInt                mx,my,mz; /* global mesh size */
	DM                      dav,dap;
  DM                      stokes_pack;
	BCList                  u_bclist,p_bclist;
	Quadrature              volQ;
	//	SurfaceQuadratureStokes surfQ[QUAD_EDGES]; /* four edges */
	PetscBool               use_mf_stokes;
};

#endif

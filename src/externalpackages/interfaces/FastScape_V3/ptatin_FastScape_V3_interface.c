/*
 Interface to call the surface process model FastScape_V3
 - In general, this interface could be used with all SPM's
 - FastScape specific calls should be embedded inside the
  #ifdef PTATIN_HAVE_FASTSCAPE_V3
 
  #endif
 
*/


#include "ptatin3d.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_redundant.h"
#include "dmda_remesh.h"
#include "dmda_view_petscvtk.h"
#include "mesh_update.h"
#include "spm_utils.h"

PetscErrorCode _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X);


#undef __FUNCT__
#define __FUNCT__ "ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3"
PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
#ifdef PTATIN_HAVE_FASTSCAPE_V3
	ierr = _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pctx,X);CHKERRQ(ierr);
#else
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatind3D must be compiled with external package <FastScape_V3: Landscape evolution model of Jean Braun>");
#endif
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3"
PetscErrorCode _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X)
{
	PetscBool       debug = PETSC_TRUE;
	PhysCompStokes  stokes;
	DM              stokes_pack,dav,dap;
	DM              dm_spmsurf0;
	PetscInt        JMAX;
	PetscErrorCode  ierr;
	
	
	PetscFunctionBegin;
	
	ierr = pTatinGetStokesContext(pctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

	/* scatter parallel surface onto rank 0 */
	ierr = DMDAGatherIKRedundantSurfaceDMDA(dav,&dm_spmsurf0);CHKERRQ(ierr);
	/* debug - write out extracted surface */
	if (debug) {
		if (dm_spmsurf0) {
			ierr = DMDAViewPetscVTK(dm_spmsurf0,PETSC_NULL,"surf_extraction_ic.vtk");CHKERRQ(ierr);
		}
	}

	if (dm_spmsurf0) {

#ifdef PTATIN_HAVE_FASTSCAPE_V3

#endif

	}

	/* scatter resulting sequential surface to mechanical model */
	ierr = DMDAScatterIKRedundantSurfaceDMDA(dm_spmsurf0,dav);CHKERRQ(ierr);
	if (dm_spmsurf0) {
		ierr = DMDestroy(&dm_spmsurf0);CHKERRQ(ierr);
	}

	/* clean up mechanical model mesh */
	ierr = DMDAGetInfo(dav,0,0,&JMAX,0,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(dav,0,JMAX);CHKERRQ(ierr);
	ierr = DMDABilinearizeQ2Elements(dav);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

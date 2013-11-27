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
#include "dmda_element_q2p1.h"
#include "dmda_redundant.h"
#include "dmda_remesh.h"
#include "dmda_view_petscvtk.h"
#include "mesh_update.h"
#include "spm_utils.h"

PetscErrorCode _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X,
																																	 PetscInt refinement_factor,
																																	 PetscReal Lstar,PetscReal Vstar,
																																	 PetscReal dt_mechanical,
																																	 int _law,double _m,double _kf,double _kd,int _bc);


#undef __FUNCT__
#define __FUNCT__ "ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3"
PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
#ifdef PTATIN_HAVE_FASTSCAPE_V3
	//ierr = _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pctx,X);CHKERRQ(ierr);
#else
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatind3D must be compiled with external package <FastScape_V3: Landscape evolution model of Jean Braun>");
#endif
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3"
PetscErrorCode _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X,
																																	 PetscInt refinement_factor,
																																	 PetscReal Lstar,PetscReal Vstar,
																																	 PetscReal dt_mechanical,
																																	 int _law,double _m,double _kf,double _kd,int _bc)
{
	PetscBool       debug = PETSC_TRUE;
	PhysCompStokes  stokes;
	DM              stokes_pack,dav,dap;
	DM              dm_spmsurf0;
	PetscInt        mx,my,mz,JMAX;
	PetscReal       gmin[3],gmax[3];
	PetscLogDouble  t0,t1;
	PetscErrorCode  ierr;
	
	
	PetscFunctionBegin;
	
	ierr = pTatinGetStokesContext(pctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

	ierr = DMDAGetSizeElementQ2(dav,&mx,&my,&mz);CHKERRQ(ierr);
	ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
	
	/* scatter parallel surface onto rank 0 */
	ierr = DMDAGatherIKRedundantSurfaceDMDA(dav,&dm_spmsurf0);CHKERRQ(ierr);
	/* debug - write out extracted surface */
	if (debug) {
		if (dm_spmsurf0) {
			ierr = DMDAViewPetscVTK(dm_spmsurf0,PETSC_NULL,"surf_extraction_ic.vtk");CHKERRQ(ierr);
		}
	}

	if (dm_spmsurf0) {
		double *sheight;
		int    nx,ny,law,bc;
		double dx,dy,dt,m,kf,kd;
		int    ii,jj,smx,smy;
		double *scoord;
		double Lx,Ly;
		
		/* generate regular 2d mesh */
		smx = mx * refinement_factor;
		smy = mx * refinement_factor;
		
		nx = (int)smx + 1;
		ny = (int)smy + 1;
		
		PetscMalloc(sizeof(PetscReal)*2*nx*ny,&scoord);
		PetscMalloc(sizeof(double)*nx*ny,&sheight);
		
		for (jj=0; jj<smy+1; jj++) {
			for (ii=0; ii<smx+1; ii++) {
				scoord[2*(ii+nx*jj)+0] = 0.0 + ii * 1.0/((double)smx);
				scoord[2*(ii+nx*jj)+1] = 0.0 + jj * 1.0/((double)smy);
				sheight[ii+nx*jj] = 0.0;
			}
		}
		ierr = InterpolateMSurf0ToSPMSurfIKGrid(dm_spmsurf0,(PetscInt)smx,(PetscInt)smy,scoord,sheight);CHKERRQ(ierr);
		
		/* interpolate topo: da_spmsurf0 -> 2d mesh */
		
		PetscGetTime(&t0);
#ifdef PTATIN_HAVE_FASTSCAPE_V3
		//fastscape_(height,&nx,&ny,&dx,&dy,&dt,&law,&m,&kf,&kd,&bc);
#endif
		PetscGetTime(&t1);

		/* interpolate topo: 2d mesh -> da_spmsurf0  */
		ierr = InterpolateSPMSurfIKGridToMSurf0((PetscInt)smx,(PetscInt)smy,scoord,sheight,dm_spmsurf0);CHKERRQ(ierr);
		if (debug) {
			ierr = DMDAViewPetscVTK(dm_spmsurf0,PETSC_NULL,"surf_extraction_interp.vtk");CHKERRQ(ierr);
		}
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

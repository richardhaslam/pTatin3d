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
PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X,
																																	PetscInt refinement_factor,
																																	PetscReal Lstar,PetscReal Vstar,
																																	PetscReal dt_mechanical,
																																	int _law,double _m,double _kf,double _kd,int _bc)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
#ifdef PTATIN_HAVE_FASTSCAPE_V3
	ierr = _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pctx,X,refinement_factor,Lstar,Vstar,dt_mechanical,_law,_m,_kf,_kd,_bc);CHKERRQ(ierr);
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

	PetscPrintf(PETSC_COMM_WORLD,"Mechanical mx: %D \n",mx);
	PetscPrintf(PETSC_COMM_WORLD,"Mechanical my: %D \n",my);
	PetscPrintf(PETSC_COMM_WORLD,"Mechanical mz: %D \n",mz);
	PetscPrintf(PETSC_COMM_WORLD,"Mechanical dt: %1.4e \n",dt_mechanical);
	
	if (dm_spmsurf0) {
		double *sheight;
		int    snx,sny,law,bc;
		double dx,dy,dt,m,kf,kd;
		int    ii,jj,smx,smy;
		double *scoord;
		double Lx,Ly;
		
		
		/* generate regular 2d mesh */
		smx = mx * refinement_factor;
		smy = mz * refinement_factor;
		
		snx = (int)smx + 1;
		sny = (int)smy + 1;
		/* FastScape requires nx be a multiple of 4.... Beware! */
		if ((snx/4)*4 != snx) {
			snx=(snx/4)*4;
		}
		if ((sny/4)*4 != sny) {
			sny=(sny/4)*4;
		}
		smx = snx - 1;
		smy = sny - 1;
		
		
		PetscPrintf(PETSC_COMM_WORLD,"SPM snx: %d \n",snx);
		PetscPrintf(PETSC_COMM_WORLD,"SPM sny: %d \n",sny);

		PetscPrintf(PETSC_COMM_WORLD,"SPM smx: %d \n",smx);
		PetscPrintf(PETSC_COMM_WORLD,"SPM smy: %d \n",smy);
		
		PetscMalloc(sizeof(PetscReal)*2*snx*sny,&scoord);
		PetscMalloc(sizeof(double)*snx*sny,&sheight);
		
		Lx = gmax[0] - gmin[0];
		Ly = gmax[2] - gmin[2];
		
		dx = Lx/((double)snx);
		dy = Ly/((double)sny);
		
		for (jj=0; jj<sny; jj++) {
			for (ii=0; ii<snx; ii++) {
				scoord[2*(ii+snx*jj)+0] = 0.5*dx + ii * dx;
				scoord[2*(ii+snx*jj)+1] = 0.5*dy + jj * dy;
				sheight[ii+snx*jj] = 0.0;
			}
		}
		ierr = InterpolateMSurf0ToSPMSurfIKGrid(dm_spmsurf0,(PetscInt)smx,(PetscInt)smy,scoord,sheight);CHKERRQ(ierr);
		
		/* interpolate topo: da_spmsurf0 -> 2d mesh */
		
		PetscGetTime(&t0);
#ifdef PTATIN_HAVE_FASTSCAPE_V3
		/* scale topo */
		for (k=0; k<snx*sny; k++) {
			sheight[k] = sheight[k] * Lstar;
		}
		/* scale geometry spacing */
		dx = dx * Lstar;
		dy = dy * Lstar;

		/* scale dt */
		
		
		dt = 0.25 * dt_mechanical * (Lstar/Vstar);
		m = 1.33;
		kf = 1.0;
		kd = 1.1;
		bc = 1;
		law = 1;
		fastscape_(sheight,&snx,&sny,&dx,&dy,&dt,&law,&m,&kf,&kd,&bc);
		
		/* unscale topo */
		for (k=0; k<snx*sny; k++) {
			sheight[k] = sheight[k] / Lstar;
		}
		
#endif
		PetscGetTime(&t1);

		/* interpolate topo: 2d mesh -> da_spmsurf0  */
		ierr = InterpolateSPMSurfIKGridToMSurf0((PetscInt)smx,(PetscInt)smy,scoord,sheight,dm_spmsurf0);CHKERRQ(ierr);
		if (debug) {
			ierr = DMDAViewPetscVTK(dm_spmsurf0,PETSC_NULL,"surf_extraction_interp.vtk");CHKERRQ(ierr);
		}
		
		PetscFree(sheight);
		PetscFree(scoord);
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

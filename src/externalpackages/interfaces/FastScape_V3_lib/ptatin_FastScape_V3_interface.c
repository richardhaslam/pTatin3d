
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
																																	 PetscReal dt_mechanical,PetscReal dt_spm,
																																	 PetscInt _law,PetscReal _m,PetscReal _kf,PetscReal _kd,PetscInt _bc);
#ifdef PTATIN_HAVE_FASTSCAPE_V3
extern void fastscape_(double *sheight,int* snx,int* sny,double* dx,double* dy,int* nsteps,int* nfreq,double* dt,int* law,double* m,double* kf,double* kd,int* bc);
#endif



#undef __FUNCT__
#define __FUNCT__ "ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3"
PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X,
																																	PetscInt refinement_factor,
																																	PetscReal Lstar,PetscReal Vstar,
																																	PetscReal dt_mechanical,PetscReal dt_spm,
																																	PetscInt _law,PetscReal _m,PetscReal _kf,PetscReal _kd,PetscInt _bc)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
#ifdef PTATIN_HAVE_FASTSCAPE_V3
	ierr = _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pctx,X,refinement_factor,Lstar,Vstar,dt_mechanical,dt_spm,_law,_m,_kf,_kd,_bc);CHKERRQ(ierr);
#else
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatind3D must be compiled with external package <FastScape_V3_lib: Landscape evolution model of Jean Braun>");
#endif
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3"
PetscErrorCode _ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X,
																																	 PetscInt refinement_factor,
																																	 PetscReal Lstar,PetscReal Vstar,
																																	 PetscReal dt_mechanical,PetscReal dt_spm,
																																	 PetscInt _law,PetscReal _m,PetscReal _kf,PetscReal _kd,PetscInt _bc)
{
	PetscBool       debug = PETSC_TRUE;
	PhysCompStokes  stokes;
	DM              stokes_pack,dav,dap;
	DM              dm_spmsurf0;
	PetscInt        k,mx,my,mz,JMAX;
	PetscReal       gmin[3],gmax[3];
	PetscLogDouble  t0,t1;
	PetscBool        flg;
	PetscErrorCode  ierr;
	
	
	PetscFunctionBegin;
	
	/* parse parameters */
	flg = PETSC_FALSE;
	PetscOptionsGetInt(NULL,"-fastscape_refinement_factor",&refinement_factor,&flg);
	if (flg) { PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] Using refinement factor %D\n",refinement_factor); }

	flg = PETSC_FALSE;
	PetscOptionsGetReal(NULL,"-fastscape_dt",&dt_spm,&flg);
	if (flg) { PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] Using dt %1.4e\n",dt_spm); }
	
	if (dt_spm > dt_mechanical) {
		dt_spm = dt_mechanical;
	}
	
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
			ierr = DMDAViewPetscVTK(dm_spmsurf0,NULL,"surf_extraction_ic.vtk");CHKERRQ(ierr);
		}
	}

	PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] Mechanical mx: %D \n",mx);
	PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] Mechanical my: %D \n",my);
	PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] Mechanical mz: %D \n",mz);
	PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] Mechanical dt: %1.4e \n",dt_mechanical);
	
	if (dm_spmsurf0) {
		double *sheight;
		int    snx,sny,law,bc,nsteps,nfreq;
		double dx,dy,dt,m,kf,kd;
		int    ii,jj,smx,smy;
		double *scoord;
		double Lx,Ly;
		
		/* generate regular 2d mesh */
		smx = mx * refinement_factor;
		smy = mz * refinement_factor;
		
		snx = (int)smx + 1;
		sny = (int)smy + 1;
		
		/* FastScape requires nx,ny be a multiple of 4.... Beware! */
		if ((snx/4)*4 != snx) {
			snx=(snx/4)*4;
		}
		if ((sny/4)*4 != sny) {
			sny=(sny/4)*4;
		}
		smx = snx - 1;
		smy = sny - 1;
		
		PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] SPM snx: %d \n",snx);
		PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] SPM sny: %d \n",sny);

		PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] SPM smx: %d \n",smx);
		PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] SPM smy: %d \n",smy);
		PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] SPM dt: %1.4e \n",dt_spm);
		
		PetscMalloc(sizeof(PetscReal)*2*snx*sny,&scoord);
		PetscMalloc(sizeof(double)*snx*sny,&sheight);
		
		Lx = gmax[0] - gmin[0];
		Ly = gmax[2] - gmin[2];
		
		dx = Lx/((double)snx);
		dy = Ly/((double)sny);
		
		/* generate coordinates for interpolation */
		for (jj=0; jj<sny; jj++) {
			for (ii=0; ii<snx; ii++) {
				scoord[2*(ii+snx*jj)+0] = 0.5*dx + ii * dx;
				scoord[2*(ii+snx*jj)+1] = 0.5*dy + jj * dy;
				sheight[ii+snx*jj] = 0.0;
			}
		}

		/* interpolate topo: da_spmsurf0 -> 2d mesh */
		ierr = InterpolateMSurf0ToSPMSurfIKGrid(dm_spmsurf0,(PetscInt)smx,(PetscInt)smy,scoord,sheight);CHKERRQ(ierr);
		
#ifdef PTATIN_HAVE_FASTSCAPE_V3
		PetscTime(&t0);
		/* scale topo */
		for (k=0; k<snx*sny; k++) {
			sheight[k] = sheight[k] * Lstar;
		}
		/* scale geometry spacing */
		dx = dx * Lstar;
		dy = dy * Lstar;

		nsteps = dt_mechanical / dt_spm;
		nsteps++;
		dt_spm = dt_mechanical/(PetscReal)nsteps;
		/* scale dt */
		dt = dt_spm * (Lstar/Vstar);
		
		/* assign */
		nsteps = dt_mechanical / dt_spm;
		nfreq = 1;
		m   = (double)_m;
		kf  = (double)_kf;
		kd  = (double)_kd; /* not used */
		bc  = (int)_bc;   /* four digit number 0101 */
		law = (int)_law;
		fastscape_(sheight,&snx,&sny,&dx,&dy,&nsteps,&nfreq,&dt,&law,&m,&kf,&kd,&bc);
		
		/* unscale topo */
		for (k=0; k<snx*sny; k++) {
			sheight[k] = sheight[k] / Lstar;
		}
		PetscTime(&t1);
		PetscPrintf(PETSC_COMM_WORLD,"  [libFastScape] Compute time %1.4e (sec)\n",t1-t0);
#endif

		/* interpolate topo: 2d mesh -> da_spmsurf0  */
		ierr = InterpolateSPMSurfIKGridToMSurf0((PetscInt)smx,(PetscInt)smy,scoord,sheight,dm_spmsurf0);CHKERRQ(ierr);
		if (debug) {
			ierr = DMDAViewPetscVTK(dm_spmsurf0,NULL,"surf_extraction_interp.vtk");CHKERRQ(ierr);
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

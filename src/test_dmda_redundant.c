
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_update_coords.h"
#include "dmda_redundant.h"


#undef __FUNCT__
#define __FUNCT__ "DMDAPerturbCoordinates"
PetscErrorCode DMDAPerturbCoordinates(DM da,PetscScalar perturbA)
{
	PetscScalar perturb;
	PetscBool flg;
	PetscReal avgdx,avgdy,avgdz;
	PetscInt i,j,k,si,sj,sk,nx,ny,nz,M,N,P;
	DM cda;
	Vec coord;
	DMDACoor3d ***_coord;
	PetscReal gmin[3],gmax[3],v0,v1,v2;
	PetscRandom rand;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	perturb = perturbA;
	ierr = PetscOptionsGetScalar(PETSC_NULL,"-perturb",&perturb,&flg);CHKERRQ(ierr);
	
	/* get average cell sizes */
	ierr = DMDAGetBoundingBox(da,gmin,gmax);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	avgdx = (gmax[0]-gmin[0])/( (PetscReal)(M-1) );
	avgdy = (gmax[1]-gmin[1])/( (PetscReal)(N-1) );
	avgdz = (gmax[2]-gmin[2])/( (PetscReal)(P-1) );
	
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	
	/* create random numbers */
	ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rand);CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rand);CHKERRQ(ierr);
	ierr = PetscRandomSetInterval(rand,-1.0,1.0);CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&_coord);CHKERRQ(ierr);
	for(k=sk;k<sk+nz;k++) {
		for(j=sj;j<sj+ny;j++) {
			for(i=si;i<si+nx;i++) {
				ierr = PetscRandomGetValueReal(rand,&v0);CHKERRQ(ierr);
				ierr = PetscRandomGetValueReal(rand,&v1);CHKERRQ(ierr);
				ierr = PetscRandomGetValueReal(rand,&v2);CHKERRQ(ierr);
				
				_coord[k][j][i].x = _coord[k][j][i].x + v0 * perturb * avgdx;
				_coord[k][j][i].y = _coord[k][j][i].y + v1 * perturb * avgdy;
				_coord[k][j][i].z = _coord[k][j][i].z + v2 * perturb * avgdz;
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda,coord,&_coord);CHKERRQ(ierr);
	
	
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test_DMDACreate3dRedundant"
PetscErrorCode test_DMDACreate3dRedundant(PetscInt nx,PetscInt ny,PetscInt nz)
{
	PetscErrorCode ierr;
	PetscReal x0,x1,y0,y1,z0,z1;
	DM da,da_red;
	Vec x;
	PetscViewer vv;
	PetscReal avgdx,avgdy,avgdz,dx;
	PetscReal gmin[3],gmax[3];
	PetscInt M,N,P;
	PetscInt nxs,nys,nzs,si,sj,sk;
	
	PetscFunctionBegin;
	
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 3,1, 0,0,0,&da);CHKERRQ(ierr);
	
	x0 = y0 = z0 = -1.0;
	x1 = y1 = z1 = 1.0;
	ierr = DMDASetUniformCoordinates(da, x0,x1, y0,y1, z0,z1);CHKERRQ(ierr);

	/* get average cell sizes */
	ierr = DMDAGetBoundingBox(da,gmin,gmax);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	avgdx = (gmax[0]-gmin[0])/( (PetscReal)(M-1) );
	avgdy = (gmax[1]-gmin[1])/( (PetscReal)(N-1) );
	avgdz = (gmax[2]-gmin[2])/( (PetscReal)(P-1) );
	dx = PetscMax(avgdx,avgdy);
	dx = PetscMax(dx,avgdz);
	ierr = DMDAPerturbCoordinates(da,0.1*dx);CHKERRQ(ierr);
	
	/*
	ierr = DMDACreate3dRedundant(da, 0,3, 5,10, 4,5, 3, &da_red );CHKERRQ(ierr);
	*/

	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&nxs,&nys,&nzs);CHKERRQ(ierr);
	ierr = DMDACreate3dRedundant(da, si,si+nxs, sj,sj+nys, P-1,P, 3, &da_red );CHKERRQ(ierr);
	
	/* output */
	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, "test_dmda_redundant_in.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr  = VecDestroy(&x);CHKERRQ(ierr);

	{
		char *name;
		PetscMPIInt rank;
		MPI_Comm_rank(((PetscObject)(da))->comm,&rank);
		asprintf(&name,"test_dmda_redundant_out_%d.vtk",rank);
		ierr = PetscViewerASCIIOpen(((PetscObject)(da_red))->comm, name, &vv);CHKERRQ(ierr);
		free(name);
	}
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_red,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da_red, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr  = VecDestroy(&x);CHKERRQ(ierr);
	
	
	ierr = DMDestroy(&da);CHKERRQ(ierr);
	ierr = DMDestroy(&da_red);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}


int main( int argc,char **argv )
{
	PetscErrorCode ierr;
	PetscInt mx,my,mz;
	
	PetscInitialize(&argc,&argv,(char *)0,0);
	
	mx = my = mz = 10;
	PetscOptionsGetInt( PETSC_NULL, "-mx", &mx, 0 );
	PetscOptionsGetInt( PETSC_NULL, "-my", &my, 0 );
	PetscOptionsGetInt( PETSC_NULL, "-mz", &mz, 0 );
	
	ierr = test_DMDACreate3dRedundant(mx,my,mz);CHKERRQ(ierr);

	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}

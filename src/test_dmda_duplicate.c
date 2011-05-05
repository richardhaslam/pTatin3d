
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define _GNU_SOURCE

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_update_coords.h"
#include "dmda_duplicate.h"
#include "dmda_view_petscvtk.h"


#undef __FUNCT__
#define __FUNCT__ "DMDAPerturbCoordinates"
PetscErrorCode DMDAPerturbCoordinates(DM da,PetscScalar perturb)
{
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
	ierr = PetscOptionsGetScalar(PETSC_NULL,"-perturb",&perturb,&flg);CHKERRQ(ierr);
	
	/* get average cell sizes */
	ierr = DMDAGetBoundingBox(da,gmin,gmax);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,0,0,0);CHKERRQ(ierr);
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
	ierr = PetscRandomDestroy(rand);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test_DMDADuplicateLayout"
PetscErrorCode test_DMDADuplicateLayout(PetscInt nx,PetscInt ny,PetscInt nz)
{
	PetscErrorCode ierr;
	PetscReal x0,x1,y0,y1,z0,z1;
	DM da,da2;
	Vec coords;
	PetscReal gmin[3],gmax[3];
	PetscReal avgdx,avgdy,avgdz,dx;
	PetscInt M,N,P;
	
	
	PetscFunctionBegin;
	
	/* create da1 */
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_NONPERIODIC,DMDA_STENCIL_BOX,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 3,1, 0,0,0,&da);CHKERRQ(ierr);

	/* set an aspect ratio box, then perturb the coords */	
	x0 = y0 = z0 = -1.0;
	x1 = y1 = z1 = 1.0;
	ierr = DMDASetUniformCoordinates(da, x0,x1, y0,y1, z0,z1);CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinates(da,&coords);CHKERRQ(ierr);

	ierr = VecStrideScale(coords,0,10.0);CHKERRQ(ierr);
	ierr = VecStrideScale(coords,1,20.0);CHKERRQ(ierr);
	ierr = VecStrideScale(coords,2,30.0);CHKERRQ(ierr);

	
	/* get average cell sizes */
	ierr = DMDAGetBoundingBox(da,gmin,gmax);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,0,0,0);CHKERRQ(ierr);
	avgdx = (gmax[0]-gmin[0])/( (PetscReal)(M-1) );
	avgdy = (gmax[1]-gmin[1])/( (PetscReal)(N-1) );
	avgdz = (gmax[2]-gmin[2])/( (PetscReal)(P-1) );
	dx = PetscMax(avgdx,avgdy);
	dx = PetscMax(dx,avgdz);
	ierr = DMDAPerturbCoordinates(da,0.1*dx);CHKERRQ(ierr);

	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	/* copy layout */
	ierr = DMDADuplicateLayout(da,5,PETSC_DECIDE,PETSC_DECIDE,&da2);CHKERRQ(ierr);
	
	/* output */
	ierr = DMDAViewPetscVTK(da, PETSC_NULL, "test_dmda_dup_1.vtk");CHKERRQ(ierr);
	ierr = DMDAViewPetscVTK(da2, PETSC_NULL, "test_dmda_dup_2.vtk");CHKERRQ(ierr);

	ierr = DMDestroy(da);CHKERRQ(ierr);
	ierr = DMDestroy(da2);CHKERRQ(ierr);
	
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
	
	ierr = test_DMDADuplicateLayout(mx,my,mz);CHKERRQ(ierr);
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}

/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   test_dmda_update_coords.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define _GNU_SOURCE

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "ptatin_init.h"
#include "dmda_update_coords.h"


#undef __FUNCT__
#define __FUNCT__ "test_DMDAUpdateGhostedCoordinates"
PetscErrorCode test_DMDAUpdateGhostedCoordinates(PetscInt nx,PetscInt ny,PetscInt nz)
{
	PetscErrorCode ierr;
	PetscReal x0,x1,y0,y1,z0,z1;
	DM da;
	Vec coords,x;
	PetscViewer vv;
	
	
	PetscFunctionBegin;
	
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 3,1, 0,0,0,&da);CHKERRQ(ierr);
	
	x0 = y0 = z0 = -1.0;
	x1 = y1 = z1 = 1.0;
	ierr = DMDASetUniformCoordinates(da, x0,x1, y0,y1, z0,z1);CHKERRQ(ierr);
	
	ierr = DMGetCoordinates(da,&coords);CHKERRQ(ierr);

//	ierr = VecScale(coords,10.0);CHKERRQ(ierr);
	ierr = VecStrideScale(coords,0,10.0);CHKERRQ(ierr);
	ierr = VecStrideScale(coords,1,20.0);CHKERRQ(ierr);
	ierr = VecStrideScale(coords,2,30.0);CHKERRQ(ierr);

	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	/* output */
	ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)da), "test1.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerPushFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	
	ierr = DMDestroy(&da);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}


int main( int argc,char **argv )
{
	PetscErrorCode ierr;
	PetscInt mx,my,mz;
	
	ierr = pTatinInitialize(&argc,&argv,(char *)0,NULL);CHKERRQ(ierr);
	
	mx = my = mz = 10;
	PetscOptionsGetInt(NULL, NULL, "-mx", &mx, 0 );
	PetscOptionsGetInt(NULL, NULL, "-my", &my, 0 );
	PetscOptionsGetInt(NULL, NULL, "-mz", &mz, 0 );
	
	ierr = test_DMDAUpdateGhostedCoordinates(mx,my,mz);CHKERRQ(ierr);
	/* todo */
	/*
	ierr = test_DMDASetCoordinatesFromLocalVector(mx,my,mz);CHKERRQ(ierr);
	ierr = test_DMDASetCoordinatesU(mx,my,mz);CHKERRQ(ierr);
	ierr = test_DMDACloneCoordinates(mx,my,mz);CHKERRQ(ierr);
	*/
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

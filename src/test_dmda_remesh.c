
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_update_coords.h"
#include "dmda_remesh.h"
#include "mesh_deformation.h"

#undef __FUNCT__
#define __FUNCT__ "test_DMDARemeshSetUniformCoordinatesInPlane_IJ"
PetscErrorCode test_DMDARemeshSetUniformCoordinatesInPlane_IJ(PetscInt nx,PetscInt ny,PetscInt nz)
{
	PetscErrorCode ierr;
	PetscReal x0,x1,y0,y1,z0,z1;
	DM da;
	Vec x;
	DMDACoor3d plane[4];
	PetscViewer vv;
	PetscInt M,N,P;
	PetscInt nxs,nys,nzs,si,sj,sk;
	
	PetscFunctionBegin;
	
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 3,1, 0,0,0,&da);CHKERRQ(ierr);
	
	x0 = y0 = z0 = -1.0;
	x1 = y1 = z1 = 1.0;
	ierr = DMDASetUniformCoordinates(da, x0,x1, y0,y1, z0,z1);CHKERRQ(ierr);

	/* output */
	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, "test_dmda_remesh_in.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);

	/* remesh */
	plane[0].x = -1.5;   plane[0].y = -1.1;   plane[0].z = -2.0;
	plane[1].x = -1.1;   plane[1].y = 1.1;  	plane[1].z = -1.8;
	plane[2].x = 1.3;    plane[2].y = 1.2;  	plane[2].z = -1.6;
	plane[3].x = 0.9;    plane[3].y = -1.2;  	plane[3].z = -1.4;

	ierr = DMDARemeshSetUniformCoordinatesInPlane_IJ(da, 0, plane );CHKERRQ(ierr);

	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, "test_dmda_remesh_out.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	
	ierr = DMDestroy(&da);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test_DMDARemeshSetUniformCoordinatesBetweenKLayers3d"
PetscErrorCode test_DMDARemeshSetUniformCoordinatesBetweenKLayers3d(PetscInt nx,PetscInt ny,PetscInt nz)
{
	PetscErrorCode ierr;
	PetscReal x0,x1,y0,y1,z0,z1;
	DM da;
	Vec x;
	DMDACoor3d plane[4];
	PetscViewer vv;
	PetscInt M,N,P;
	
	PetscFunctionBegin;
	
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 3,1, 0,0,0,&da);CHKERRQ(ierr);
	
	x0 = y0 = z0 = -1.0;
	x1 = y1 = z1 = 1.0;
	ierr = DMDASetUniformCoordinates(da, x0,x1, y0,y1, z0,z1);CHKERRQ(ierr);
	
	/* output */
	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, "test_dmda_remesh_in.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr  = VecDestroy(&x);CHKERRQ(ierr);
	
	/* remesh */
	plane[0].x = -1.5;   plane[0].y = -1.1;   plane[0].z = -2.0;
	plane[1].x = -1.1;   plane[1].y = 1.1;  	plane[1].z = -1.8;
	plane[2].x = 1.3;    plane[2].y = 1.2;  	plane[2].z = -1.6;
	plane[3].x = 0.9;    plane[3].y = -1.2;  	plane[3].z = -1.4;
	
	ierr = DMDARemeshSetUniformCoordinatesInPlane_IJ(da, 0, plane );CHKERRQ(ierr);

	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, "test_dmda_remesh_in1.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr  = VecDestroy(&x);CHKERRQ(ierr);
	
	
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenKLayers3d(da,0,P);CHKERRQ(ierr);
	
	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, "test_dmda_remesh_out.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr  = VecDestroy(&x);CHKERRQ(ierr);
	
	ierr = DMDestroy(&da);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test_DMDARemeshSetUniformCoordinatesBetweenKLayers3d_b"
PetscErrorCode test_DMDARemeshSetUniformCoordinatesBetweenKLayers3d_b(PetscInt nx,PetscInt ny,PetscInt nz)
{
	PetscErrorCode ierr;
	PetscReal x0,x1,y0,y1,z0,z1;
	DM da;
	Vec x;
	DMDACoor3d plane[4];
	PetscViewer vv;
	PetscInt M,N,P;
	
	PetscFunctionBegin;
	
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 3,1, 0,0,0,&da);CHKERRQ(ierr);
	
	x0 = y0 = z0 = 0.0;
	x1 = y1 = z1 = 1.0;
	ierr = DMDASetUniformCoordinates(da, x0,x1, y0,y1, z0,z1);CHKERRQ(ierr);
	
	/* output */
	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, "test_dmda_remesh_in.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr  = VecDestroy(&x);CHKERRQ(ierr);
	
	ierr = MeshDeformation_ShearXY(da);CHKERRQ(ierr);
	
	/* output sheared mesh */
	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, "test_dmda_remesh_out.vtk", &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
	ierr = DMView(da, vv);CHKERRQ(ierr);
	ierr = VecView(x, vv);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
	ierr  = VecDestroy(&x);CHKERRQ(ierr);
	
	ierr = DMDestroy(&da);CHKERRQ(ierr);
	
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
	
//	ierr = test_DMDARemeshSetUniformCoordinatesInPlane_IJ(mx,my,mz);CHKERRQ(ierr);

//	ierr = test_DMDARemeshSetUniformCoordinatesBetweenKLayers3d(mx,my,mz);CHKERRQ(ierr);

	ierr = test_DMDARemeshSetUniformCoordinatesBetweenKLayers3d_b(mx,my,mz);CHKERRQ(ierr);
	
	ierr = PetscFinalize();CHKERRQ(ierr);

	return 0;
}

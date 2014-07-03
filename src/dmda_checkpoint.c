/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      dmda_checkpoint.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_update_coords.h"
#include "dmda_checkpoint.h"

typedef enum { DMDA_DIM=0, DMDA_M,DMDA_N,DMDA_P, DMDA_DOF, DMDA_SW, DMDA_WRAPX,DMDA_WRAPY,DMDA_WRAPZ, DMDA_ST, DMDA_RX,DMDA_RY,DMDA_RZ, DMDA_COORD, DMDA_ENDFLAG } DMDAFields;

#undef __FUNCT__  
#define __FUNCT__ "DMDAPackDataToFile"
PetscErrorCode DMDAPackDataToFile(DM da,const char name[])
{
	PetscErrorCode ierr;
	PetscViewer    v;
	MPI_Comm       comm;
	PetscMPIInt    rank;
	PetscInt       i,L, dim, M,N,P, dof, sw;
	PetscInt       refine_x, refine_y, refine_z;
	DMBoundaryType wrap[3];
	DMDAStencilType  st;
	PetscScalar      val;
	Vec              dd,coords;
	PetscBool        has_coords;
	
	
	PetscFunctionBegin;
	if (da == NULL) SETERRQ( PETSC_COMM_WORLD,PETSC_ERR_USER, "da is NULL" );

	ierr = DMGetCoordinates(da,&coords);CHKERRQ(ierr);
	if (coords)  { has_coords = PETSC_TRUE;  }
	if (!coords) { has_coords = PETSC_FALSE; }
	
	/* write coordinates out to disk */
	if (has_coords == PETSC_TRUE) {
		char coord_file[PETSC_MAX_PATH_LEN];

		sprintf(coord_file,"%s.coords",name);
		
		ierr = PetscViewerCreate(PetscObjectComm((PetscObject)da),&v);CHKERRQ(ierr);
		ierr = PetscViewerSetType(v,PETSCVIEWERBINARY);CHKERRQ(ierr);
		ierr = PetscViewerFileSetMode(v,FILE_MODE_WRITE);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO
		ierr = PetscViewerBinarySetMPIIO(v);CHKERRQ(ierr);
#endif
		ierr = PetscViewerFileSetName(v,coord_file);CHKERRQ(ierr);
		ierr = VecView(coords,v);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&v);CHKERRQ(ierr);
	}
	
	PetscObjectGetComm( (PetscObject)da, &comm );
	ierr = MPI_Comm_rank( comm, &rank );CHKERRQ(ierr);
	if (rank != 0) { PetscFunctionReturn(0); }
	
	ierr = DMDAGetInfo( da, &dim, &M,&N,&P, 0,0,0, &dof, &sw, &wrap[0],&wrap[1],&wrap[2], &st );CHKERRQ(ierr);
	ierr = DMDAGetRefinementFactor( da, &refine_x, &refine_y, &refine_z );CHKERRQ(ierr);
	
	// dim , M,N,P , dof , sw , wrap , stencil , refx,refy,refz has_coords = 12
	//  0    1,2,3    4     5    6        7          8,9,10     11
	ierr = VecCreate( PETSC_COMM_SELF, &dd );CHKERRQ(ierr);
	L = (PetscInt)(DMDA_ENDFLAG);

	ierr = VecSetSizes( dd, PETSC_DECIDE, L );CHKERRQ(ierr);
	ierr = VecSetType( dd, VECSEQ );CHKERRQ(ierr);
	
	val = (PetscScalar)dim + 0.1;		VecSetValue( dd, (PetscInt)DMDA_DIM, val, INSERT_VALUES );
	
	val = (PetscScalar)M + 0.1;			VecSetValue( dd, (PetscInt)DMDA_M, val, INSERT_VALUES );
	val = (PetscScalar)N + 0.1;			VecSetValue( dd, (PetscInt)DMDA_N, val, INSERT_VALUES );
	val = (PetscScalar)P + 0.1;			VecSetValue( dd, (PetscInt)DMDA_P, val, INSERT_VALUES );
	
	val = (PetscScalar)dof + 0.1;		VecSetValue( dd, (PetscInt)DMDA_DOF, val, INSERT_VALUES );
	val = (PetscScalar)sw + 0.1;		VecSetValue( dd, (PetscInt)DMDA_SW, val, INSERT_VALUES );
	
	//
	for (i=0; i<3; i++) {
		switch (wrap[i]) {
			case DM_BOUNDARY_NONE:
				val = 0.1;
				break;
			case DM_BOUNDARY_GHOSTED:
				val = 1.1;
				break;
			case DM_BOUNDARY_MIRROR:
				val = 2.1;
				break;
			case DM_BOUNDARY_PERIODIC:
				val = 3.1;
				break;
		}
		if (i == 0) { ierr = VecSetValue( dd, (PetscInt)DMDA_WRAPX, val, INSERT_VALUES );CHKERRQ(ierr); }
		if (i == 1) { ierr = VecSetValue( dd, (PetscInt)DMDA_WRAPY, val, INSERT_VALUES );CHKERRQ(ierr); }
		if (i == 2) { ierr = VecSetValue( dd, (PetscInt)DMDA_WRAPZ, val, INSERT_VALUES );CHKERRQ(ierr); }
	}
	//
	
	switch (st) {
		case DMDA_STENCIL_STAR:
			val = 0.1;
			break;
		case DMDA_STENCIL_BOX:
			val = 1.1;
			break;
	}			
	ierr = VecSetValue( dd, (PetscInt)DMDA_ST, val, INSERT_VALUES );CHKERRQ(ierr);
	
	/* ref x,y,z */
	val = (PetscScalar)refine_x + 0.1;			ierr = VecSetValue( dd, (PetscInt)DMDA_RX, val, INSERT_VALUES );CHKERRQ(ierr);
	val = (PetscScalar)refine_y + 0.1;			ierr = VecSetValue( dd, (PetscInt)DMDA_RY, val, INSERT_VALUES );CHKERRQ(ierr);
	val = (PetscScalar)refine_z + 0.1;			ierr = VecSetValue( dd, (PetscInt)DMDA_RZ, val, INSERT_VALUES );CHKERRQ(ierr);
	
	if (has_coords == PETSC_TRUE) {
		ierr = VecSetValue( dd, (PetscInt)DMDA_COORD, 1, INSERT_VALUES );CHKERRQ(ierr);
	} else {
		ierr = VecSetValue( dd, (PetscInt)DMDA_COORD, 0, INSERT_VALUES );CHKERRQ(ierr);
	}
	
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_WRITE, &v );CHKERRQ(ierr);
	ierr = VecView(dd,v);CHKERRQ(ierr);
	
	ierr = VecDestroy(&dd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&v);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DMDACreateFromPackDataToFile"
PetscErrorCode DMDACreateFromPackDataToFile(MPI_Comm comm,const char name[],DM *da)
{
	PetscErrorCode ierr;
	PetscViewer    v;
	Vec            dd;
	PetscScalar    *data;
	PetscInt       dim, M,N,P, dof, sw;
	PetscInt       refine_x, refine_y, refine_z;
	DMBoundaryType wrap[3];
	DMDAStencilType  st;
	PetscInt         i,L,convert;
	PetscBool        has_coords;
	
	
	PetscFunctionBegin;
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_READ, &v );CHKERRQ(ierr);

	ierr = VecCreate( PETSC_COMM_SELF, &dd );CHKERRQ(ierr);
	L = (PetscInt)(DMDA_ENDFLAG);
	ierr = VecSetSizes( dd, PETSC_DECIDE, L );CHKERRQ(ierr);
	ierr = VecSetType( dd, VECSEQ );CHKERRQ(ierr);

	ierr = VecLoad(dd,v);CHKERRQ(ierr);
	/* 
	 putain - VecLoadIntoVector inserts the option below into the command line.
	 This will screw shit up if you load in vectors with different block sizes.
	 */
	ierr = PetscOptionsClearValue("-vecload_block_size");CHKERRQ(ierr);
	
	ierr = PetscViewerDestroy(&v);CHKERRQ(ierr);
	
	ierr = VecGetArray( dd, &data );CHKERRQ(ierr);
	
	//
	dim = (PetscInt)data[DMDA_DIM];
	
	M = (PetscInt)data[DMDA_M];
	N = (PetscInt)data[DMDA_N];
	P = (PetscInt)data[DMDA_P];
	
	dof = (PetscInt)data[DMDA_DOF];
	sw = (PetscInt)data[DMDA_SW];
	//
	for (i=0; i<3; i++) {
		if (i == 0) { convert = (PetscInt)data[DMDA_WRAPX]; }
		if (i == 1) { convert = (PetscInt)data[DMDA_WRAPY]; }
		if (i == 2) { convert = (PetscInt)data[DMDA_WRAPZ]; }

		switch (convert) {
			case 0:
				wrap[i] = DM_BOUNDARY_NONE;
				break;
			case 1:
				wrap[i] = DM_BOUNDARY_GHOSTED;
				break;
			case 2:
				wrap[i] = DM_BOUNDARY_MIRROR;
				break;
			case 3:
				wrap[i] = DM_BOUNDARY_PERIODIC;
				break;
		}
	}
	convert = (PetscInt)data[DMDA_ST];
	st = DMDA_STENCIL_STAR;
	switch (convert) {
		case 0:
			st = DMDA_STENCIL_STAR;
			break;
		case 1:
			st = DMDA_STENCIL_BOX;
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unknown stenctil type detected");
			break;
	}			
	
	/* ref x,y,z */
	refine_x = (PetscInt)data[DMDA_RX];
	refine_y = (PetscInt)data[DMDA_RY];
	refine_z = (PetscInt)data[DMDA_RZ];

	has_coords = PETSC_FALSE;
	if ((PetscInt)data[DMDA_COORD] == 1) {
		has_coords = PETSC_TRUE;
	}
	
	ierr = VecRestoreArray( dd, &data );CHKERRQ(ierr);
	ierr = VecDestroy(&dd);CHKERRQ(ierr);
	
    if (dim == 3) {
        ierr = DMDACreate3d( comm, wrap[0],wrap[1],wrap[2], st, M,N,P, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, dof,sw, 0,0,0, da );CHKERRQ(ierr);
    } else {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only dimension=3 supported");
    }
	ierr = DMDASetRefinementFactor( *da, refine_x, refine_y, refine_z );CHKERRQ(ierr);
	
	/* write coordinates out to disk */
	if (has_coords == PETSC_TRUE) {
		char coord_file[PETSC_MAX_PATH_LEN];
		Vec  da_coords;
		
		ierr = DMDASetUniformCoordinates(*da, 0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
		ierr = DMGetCoordinates(*da, &da_coords);CHKERRQ(ierr);
		
		sprintf(coord_file,"%s.coords",name);

		ierr = PetscViewerCreate(PetscObjectComm((PetscObject)*da),&v);CHKERRQ(ierr);
		ierr = PetscViewerSetType(v,PETSCVIEWERBINARY);CHKERRQ(ierr);
		ierr = PetscViewerFileSetMode(v,FILE_MODE_READ);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO
		ierr = PetscViewerBinarySetMPIIO(v);CHKERRQ(ierr);
#endif
		ierr = PetscViewerFileSetName(v,coord_file);CHKERRQ(ierr);
		ierr = VecLoad(da_coords,v);CHKERRQ(ierr);

		/* 
		 putain - VecLoadIntoVector inserts the option below into the command line.
		 This will screw shit up if you load in vectors with different block sizes.
		 */
		ierr = PetscOptionsClearValue("-vecload_block_size"); CHKERRQ(ierr);
		
		ierr = PetscViewerDestroy(&v);CHKERRQ(ierr);		
		ierr = DMDAUpdateGhostedCoordinates(*da);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DMDALoadGlobalVectorFromFile"
PetscErrorCode DMDALoadGlobalVectorFromFile(DM da,const char name[],Vec *da_x)
{
	PetscErrorCode ierr;
	PetscViewer    v;
	MPI_Comm       comm;
	Vec            xn;
	

	PetscFunctionBegin;
	PetscObjectGetComm( (PetscObject)da, &comm );
	if (da == NULL) SETERRQ(comm,PETSC_ERR_USER, "da is NULL");
	
	ierr = DMCreateGlobalVector(da, &xn);CHKERRQ(ierr);

	ierr = PetscViewerCreate(comm,&v);CHKERRQ(ierr);
	ierr = PetscViewerSetType(v,PETSCVIEWERBINARY);CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(v,FILE_MODE_READ);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO
	ierr = PetscViewerBinarySetMPIIO(v);CHKERRQ(ierr);
#endif
	ierr = PetscViewerFileSetName(v,name);CHKERRQ(ierr);
	ierr = VecLoad(xn,v); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
	
	/* 
	 putain - VecLoadIntoVector inserts the option below into the command line.
	 This will screw shit up if you load in vectors with different block sizes.
	 */
	ierr = PetscOptionsClearValue("-vecload_block_size");CHKERRQ(ierr);
	
	*da_x = xn;
	
	PetscFunctionReturn(0);	
}

#undef __FUNCT__  
#define __FUNCT__ "DMDALoadCoordinatesFromFile"
PetscErrorCode DMDALoadCoordinatesFromFile(DM da,const char name[])
{
	PetscErrorCode ierr;
	DM             cda;
	Vec            coords,da_coords;
	
	
	PetscFunctionBegin;
	if (da == NULL) SETERRQ( PetscObjectComm((PetscObject)da),PETSC_ERR_USER, "da is NULL" );
	
	/* make sure the vector is present */
	ierr = DMDASetUniformCoordinates(da, 0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMDALoadGlobalVectorFromFile(cda,name,&coords);CHKERRQ(ierr);
	
	/* set the global coordinates */
	ierr = DMGetCoordinates(da,&da_coords);CHKERRQ(ierr);
	ierr = VecCopy(coords,da_coords);CHKERRQ(ierr);
	
	/* make sure the local coordinates are upto date */
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	ierr = VecDestroy(&coords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);		
}

#undef __FUNCT__  
#define __FUNCT__ "DMDAWriteVectorToFile"
PetscErrorCode DMDAWriteVectorToFile(Vec x,const char name[],PetscBool zip_file)
{
	char fieldname[PETSC_MAX_PATH_LEN];
	PetscViewer viewer;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	if (zip_file) {
		sprintf(fieldname,"%s.gz",name);
	} else {
		sprintf(fieldname,"%s",name);
	}
	
  ierr = PetscViewerCreate(PetscObjectComm((PetscObject)x),&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO	
	ierr = PetscViewerBinarySetMPIIO(viewer);CHKERRQ(ierr);
#endif
	ierr = PetscViewerFileSetName(viewer,fieldname);CHKERRQ(ierr);
	
	ierr = VecView(x,viewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecLoadFromFile"
PetscErrorCode VecLoadFromFile(Vec x,const char name[])
{
	PetscErrorCode ierr;
	PetscViewer    v;
	
	
	PetscFunctionBegin;
	
	ierr = PetscViewerCreate(PetscObjectComm((PetscObject)x),&v);CHKERRQ(ierr);
	ierr = PetscViewerSetType(v,PETSCVIEWERBINARY);CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(v,FILE_MODE_READ);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO
	ierr = PetscViewerBinarySetMPIIO(v);CHKERRQ(ierr);
#endif
	ierr = PetscViewerFileSetName(v,name);CHKERRQ(ierr);
	ierr = VecLoad(x,v); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
	
	/* 
	 putain - VecLoadIntoVector inserts the option below into the command line.
	 This will screw shit up if you load in vectors with different block sizes.
	 */
	ierr = PetscOptionsClearValue("-vecload_block_size");CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}


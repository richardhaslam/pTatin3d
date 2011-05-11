

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
	PetscViewer v;
	MPI_Comm comm;
	PetscMPIInt rank;
	PetscInt i,L, dim, M,N,P, dof, sw;
	PetscInt refine_x, refine_y, refine_z;
	DMDABoundaryType wrap[3];
	DMDAStencilType st;
	PetscScalar val;
	Vec dd;
	Vec coords;
	PetscBool has_coords;
	
	PetscFunctionBegin;
	if( da == PETSC_NULL ) SETERRQ( PETSC_COMM_WORLD,PETSC_ERR_USER, "da is NULL" );

	ierr = DMDAGetCoordinates(da,&coords);CHKERRQ(ierr);
	if (coords)  { has_coords = PETSC_TRUE;  }
	if (!coords) { has_coords = PETSC_FALSE; }
	
	/* write coordinates out to disk */
	if (has_coords==PETSC_TRUE) {
		char coord_file[1000];

		sprintf(coord_file,"coords_%s",name);
		ierr = PetscViewerBinaryOpen( ((PetscObject)da)->comm, coord_file, FILE_MODE_WRITE, &v );CHKERRQ(ierr);
		ierr = VecView( coords, v );CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&v);CHKERRQ(ierr);		
	}
	
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_rank( comm, &rank );
	if( rank != 0 ) { PetscFunctionReturn(0); }
	
	
	ierr = DMDAGetInfo( da, &dim, &M,&N,&P, 0,0,0, &dof, &sw, &wrap[0],&wrap[1],&wrap[2], &st ); CHKERRQ(ierr);
	ierr = DMDAGetRefinementFactor( da, &refine_x, &refine_y, &refine_z ); CHKERRQ(ierr);
	
	// dim , M,N,P , dof , sw , wrap , stencil , refx,refy,refz has_coords = 12
	//  0    1,2,3    4     5    6        7          8,9,10     11
	ierr = VecCreate( PETSC_COMM_SELF, &dd ); CHKERRQ(ierr);
	L = (PetscInt)(DMDA_ENDFLAG);
	printf("L=%d \n", L );
	ierr = VecSetSizes( dd, PETSC_DECIDE, L ); CHKERRQ(ierr);
	ierr = VecSetType( dd, VECSEQ ); CHKERRQ(ierr);
	
	val = (PetscScalar)dim + 0.1;		VecSetValue( dd, DMDA_DIM, val, INSERT_VALUES );
	
	val = (PetscScalar)M + 0.1;			VecSetValue( dd, DMDA_M, val, INSERT_VALUES );
	val = (PetscScalar)N + 0.1;			VecSetValue( dd, DMDA_N, val, INSERT_VALUES );
	val = (PetscScalar)P + 0.1;			VecSetValue( dd, DMDA_P, val, INSERT_VALUES );
	
	val = (PetscScalar)dof + 0.1;		VecSetValue( dd, DMDA_DOF, val, INSERT_VALUES );
	val = (PetscScalar)sw + 0.1;		VecSetValue( dd, DMDA_SW, val, INSERT_VALUES );
	
	//
	for (i=0; i<3; i++) {
	switch(wrap[i]) {
		case DMDA_BOUNDARY_NONE:
			val = 0.1;
			break;
		case DMDA_BOUNDARY_GHOSTED:
			val = 1.1;
			break;
		case DMDA_BOUNDARY_MIRROR:
			val = 2.1;
			break;
		case DMDA_BOUNDARY_PERIODIC:
			val = 3.1;
			break;
/*
		case DMDA_XPERIODIC:
			val = 1.1;
			break;
		case DMDA_YPERIODIC:
			val = 2.1;
			break;
		case DMDA_XYPERIODIC:
			val = 3.1;
			break;
		case DMDA_XYZPERIODIC:
			val = 4.1;
			break;
		case DMDA_XZPERIODIC:
			val = 5.1;
			break;
		case DMDA_YZPERIODIC:
			val = 6.1;
			break;
		case DMDA_ZPERIODIC:
			val = 7.1;
			break;
		case DMDA_XYZGHOSTED:
			val = 8.1;
			break;
*/
	}
	if (i==0) { ierr = VecSetValue( dd, DMDA_WRAPX, val, INSERT_VALUES ); CHKERRQ(ierr); }
	if (i==1) { ierr = VecSetValue( dd, DMDA_WRAPY, val, INSERT_VALUES ); CHKERRQ(ierr); }
	if (i==2) { ierr = VecSetValue( dd, DMDA_WRAPZ, val, INSERT_VALUES ); CHKERRQ(ierr); }
	}
	//
	
	switch(st) {
		case DMDA_STENCIL_STAR:
			val = 0.1;
			break;
		case DMDA_STENCIL_BOX:
			val = 1.1;
			break;
	}			
	ierr = VecSetValue( dd, DMDA_ST, val, INSERT_VALUES ); CHKERRQ(ierr);
	
	/* ref x,y,z */
	val = (PetscScalar)refine_x + 0.1;			ierr = VecSetValue( dd, DMDA_RX, val, INSERT_VALUES );CHKERRQ(ierr);
	val = (PetscScalar)refine_y + 0.1;			ierr = VecSetValue( dd, DMDA_RY, val, INSERT_VALUES );CHKERRQ(ierr);
	val = (PetscScalar)refine_z + 0.1;			ierr = VecSetValue( dd, DMDA_RZ, val, INSERT_VALUES );CHKERRQ(ierr);
	
	if (has_coords==PETSC_TRUE) {
		ierr = VecSetValue( dd, DMDA_COORD, 1, INSERT_VALUES );CHKERRQ(ierr);
	} else {
		ierr = VecSetValue( dd, DMDA_COORD, 0, INSERT_VALUES );CHKERRQ(ierr);
	}
	
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_WRITE, &v ); CHKERRQ(ierr);
	ierr = VecView(dd,v); CHKERRQ(ierr);
	
	ierr = VecDestroy(&dd); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DMDACreateFromPackDataToFile"
PetscErrorCode DMDACreateFromPackDataToFile(MPI_Comm comm,const char name[],DM *da)
{
	PetscErrorCode ierr;
	PetscViewer v;
	Vec dd;
	PetscScalar *data;
	PetscInt dim, M,N,P, dof, sw;
	PetscInt refine_x, refine_y, refine_z;
	DMDABoundaryType wrap[3];
	DMDAStencilType st;
	PetscInt i,L,convert;
	PetscBool has_coords;
	
	PetscFunctionBegin;
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, name, FILE_MODE_READ, &v ); CHKERRQ(ierr);

	ierr = VecCreate( PETSC_COMM_SELF, &dd ); CHKERRQ(ierr);
	L = (PetscInt)(DMDA_ENDFLAG);
	ierr = VecSetSizes( dd, PETSC_DECIDE, L ); CHKERRQ(ierr);
	ierr = VecSetType( dd, VECSEQ ); CHKERRQ(ierr);

	ierr = VecLoad(dd,v); CHKERRQ(ierr);

	ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
	
	ierr = VecGetArray( dd, &data ); CHKERRQ(ierr);
	
	////
	dim = (PetscInt)data[DMDA_DIM];
	
	M = (PetscInt)data[DMDA_M];
	N = (PetscInt)data[DMDA_N];
	P = (PetscInt)data[DMDA_P];
	
	dof = (PetscInt)data[DMDA_DOF];
	sw = (PetscInt)data[DMDA_SW];
	
	
	//
	for (i=0; i<3; i++) {
	if (i==0) { convert = (PetscInt)data[DMDA_WRAPX]; }
	if (i==1) { convert = (PetscInt)data[DMDA_WRAPY]; }
	if (i==2) { convert = (PetscInt)data[DMDA_WRAPZ]; }

	switch(convert) {
		case 0:
			wrap[i] = DMDA_BOUNDARY_NONE;
			break;
		case 1:
			wrap[i] = DMDA_BOUNDARY_GHOSTED;
			break;
		case 2:
			wrap[i] = DMDA_BOUNDARY_MIRROR;
			break;
		case 3:
			wrap[i] = DMDA_BOUNDARY_PERIODIC;
			break;
	}
	}
	convert = (PetscInt)data[DMDA_ST];
	switch(convert) {
		case 0:
			st = DMDA_STENCIL_STAR;
			break;
		case 1:
			st = DMDA_STENCIL_BOX;
			break;
	}			
	
	/* ref x,y,z */
	refine_x = (PetscInt)data[DMDA_RX];
	refine_y = (PetscInt)data[DMDA_RY];
	refine_z = (PetscInt)data[DMDA_RZ];

	has_coords = PETSC_FALSE;
	if ( (PetscInt)data[DMDA_COORD]==1 ) {
		has_coords = PETSC_TRUE;
	}
	
	ierr = VecRestoreArray( dd, &data ); CHKERRQ(ierr);
	ierr = VecDestroy(&dd); CHKERRQ(ierr);
	
	ierr = DMDACreate3d( comm, wrap[0],wrap[1],wrap[2], st, M,N,P, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, dof,sw, 0,0,0, da ); CHKERRQ(ierr);
	ierr = DMDASetRefinementFactor( *da, refine_x, refine_y, refine_z ); CHKERRQ(ierr);
	
	/* write coordinates out to disk */
	if (has_coords==PETSC_TRUE) {
		char coord_file[1000];
		DM cda;
		Vec da_coords;
		
		ierr = DMDASetUniformCoordinates(*da, 0.0,1.0,0.0,1.0,0.0,1.0); CHKERRQ(ierr);
		ierr = DMDAGetCoordinates(*da, &da_coords);CHKERRQ(ierr);
		
		sprintf(coord_file,"coords_%s",name);
		ierr = PetscViewerBinaryOpen( ((PetscObject)(*da))->comm, coord_file, FILE_MODE_READ, &v );CHKERRQ(ierr);
		ierr = VecLoad(da_coords,v);CHKERRQ(ierr);
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
	PetscViewer v;
	MPI_Comm comm;
	Vec xn;
	

	PetscFunctionBegin;
	PetscObjectGetComm( (PetscObject)da, &comm );
	if (da==PETSC_NULL) SETERRQ(comm,PETSC_ERR_USER, "da is NULL");
	
	ierr = DMCreateGlobalVector(da, &xn); CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(comm, name, FILE_MODE_READ, &v); CHKERRQ(ierr);
	ierr = VecLoad(xn,v); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
	
	/* 
	 putain - VecLoadIntoVector inserts the option below into the command line.
	 This will screw shit up if you load in vectors with different block sizes.
	 */
	ierr = PetscOptionsClearValue("-vecload_block_size"); CHKERRQ(ierr);
	
	*da_x = xn;
	
	PetscFunctionReturn(0);	
}

#undef __FUNCT__  
#define __FUNCT__ "DMDALoadCoordinatesFromFile"
PetscErrorCode DMDALoadCoordinatesFromFile(DM da,const char name[])
{
	PetscErrorCode ierr;
	DM cda;
	Vec coords, da_coords;
	
	PetscFunctionBegin;
	if( da == PETSC_NULL ) SETERRQ( ((PetscObject)da)->comm,PETSC_ERR_USER, "da is NULL" );
	
	/* make sure the vector is present */
	ierr = DMDASetUniformCoordinates(da, 0.0,1.0,0.0,1.0,0.0,1.0); CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinateDA(da, &cda); CHKERRQ(ierr);
	ierr = DMDALoadGlobalVectorFromFile(cda, name, &coords); CHKERRQ(ierr);
	
	/* set the global coordinates */
	ierr = DMDAGetCoordinates(da, &da_coords);CHKERRQ(ierr);
	ierr = VecCopy(coords,da_coords);CHKERRQ(ierr);
	
	/* make sure the local coordinates are upto date */
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	ierr = VecDestroy(&coords); CHKERRQ(ierr);
	
	PetscFunctionReturn(0);		
}


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
 **    filename:   dmda_redundant.c
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
#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "petsc/private/dmdaimpl.h"    /*I   "petscdm.h"   I*/
#include "sub_comm.h"
#include "dmda_update_coords.h"
#include "dmda_redundant.h"



/* This is only needed as the DMDACreate3d() forces a call to DMSetFromOptions().
 This parses the options -da_processors_x,y,z even when the comm size is 1, which screws up my seq surface DM.
 This function is a copy of DMDACreate3d() without the call to DMSetFromOptions()
  Note: this can now be removed, as this call was taken out in PETSc 3.8
*/
PetscErrorCode  x_DMDACreate3d(
  MPI_Comm comm,DMBoundaryType wrap[],DMDAStencilType stencil_type,
  PetscInt M, PetscInt N,PetscInt P,PetscInt m,PetscInt n,PetscInt p,
  PetscInt dof,PetscInt s,const PetscInt lx[],const PetscInt ly[],const PetscInt lz[],DM *da)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMDACreate(comm, da);CHKERRQ(ierr);
  ierr = DMSetDimension(*da, 3);CHKERRQ(ierr);
  ierr = DMDASetSizes(*da, M, N, P);CHKERRQ(ierr);
  ierr = DMDASetNumProcs(*da, m, n, p);CHKERRQ(ierr);
  ierr = DMDASetBoundaryType(*da, wrap[0],wrap[1],wrap[2]);CHKERRQ(ierr);
  ierr = DMDASetDof(*da, dof);CHKERRQ(ierr);
  ierr = DMDASetStencilType(*da, stencil_type);CHKERRQ(ierr);
  ierr = DMDASetStencilWidth(*da, s);CHKERRQ(ierr);
  ierr = DMDASetOwnershipRanges(*da, lx, ly, lz);CHKERRQ(ierr);
  /* This violates the behavior for other classes, but right now users expect negative dimensions to be handled this way */
/*  ierr = DMSetFromOptions(*da);CHKERRQ(ierr); */
  ierr = DMSetUp(*da);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
 Given a DA (possibly distributed, we form a redundant DMDA which spans
 the nodes defined via (si,ei), (sj,ej) and (sk,ek).

 The nodes which are retrieved are:
 si , si + 1, si + 2, ... ei-1
 */
PetscErrorCode DMDACreate3dRedundant(DM da,PetscInt si, PetscInt ei, PetscInt sj, PetscInt ej, PetscInt sk, PetscInt ek,PetscInt n_dofs,DM *_seq_DA )
{
	PetscInt M,N,P,dof;
	PetscInt sw;
	DMBoundaryType wrap[3];
	DMBoundaryType wrapNP[] = {DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE};
	DMDAStencilType st;
	DM _sda;
	PetscInt i,j,k,c;
	PetscInt Ml,Nl,Pl;
	PetscInt *fine_indices;
	IS is_local, is_fine;
	VecScatter sctx;
	DM cda_fine;
	Vec global_coords_fine, global_coords_fine_NATURAL;
	Vec new_local_coords;
	PetscErrorCode ierr;


	ierr = DMDAGetInfo( da, 0, &M,&N,&P, 0,0,0, &dof,&sw,&wrap[0],&wrap[1],&wrap[2],&st );CHKERRQ(ierr);
	/* check ranges */
	if( si < 0 || (ei-1) > M-1 ) {
		SETERRQ( PetscObjectComm((PetscObject)da), PETSC_ERR_USER, "Specified range of nodes in i-direction to gather is too large" );
	}
	if( sj < 0 || (ej-1) > N-1 ) {
		SETERRQ( PetscObjectComm((PetscObject)da), PETSC_ERR_USER, "Specified range of nodes in j-direction to gather is too large" );
	}
	if( sk < 0 || (ek-1) > P-1 ) {
		SETERRQ( PetscObjectComm((PetscObject)da), PETSC_ERR_USER, "Specified range of nodes in k-direction to gather is too large" );
	}

	/* check for cases where user wanted a plane */
	if (si == ei) {
		SETERRQ( PetscObjectComm((PetscObject)da), PETSC_ERR_USER, "Specified range of nodes in i-direction too. If you wanted a plane, indicate si,si+1 as your range" );
	}
	if (sj == ej) {
		SETERRQ( PetscObjectComm((PetscObject)da), PETSC_ERR_USER, "Specified range of nodes in j-direction too. If you wanted a plane, indicate sj,sj+1 as your range" );
	}
	if (sk == ek) {
		SETERRQ( PetscObjectComm((PetscObject)da), PETSC_ERR_USER, "Specified range of nodes in k-direction too. If you wanted a plane, indicate sk,sk+1 as your range" );
	}


	/* modify the bounds if they go outside */
	/* why do we want to do this??? */
	/*
	if( si < 0 ) { si = 0; }
	if( sj < 0 ) { sj = 0; }
	if( sk < 0 ) { sk = 0; }

	if( ei > M ) { ei = M; }
	if( ej > N ) { ej = N; }
	if( ek > P ) { ek = P; }
	*/

	/* Get the coordinate vector from the distributed array */
	ierr = DMGetCoordinateDM(da,&cda_fine);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&global_coords_fine);CHKERRQ(ierr); /* a global vector */
	ierr = DMDACreateNaturalVector( cda_fine, &global_coords_fine_NATURAL );CHKERRQ(ierr);

	ierr = DMDAGlobalToNaturalBegin( cda_fine,global_coords_fine,INSERT_VALUES, global_coords_fine_NATURAL );CHKERRQ(ierr);
	ierr = DMDAGlobalToNaturalEnd( cda_fine,global_coords_fine,INSERT_VALUES, global_coords_fine_NATURAL );CHKERRQ(ierr);


	/* get indices of the guys I want to grab */
	Ml = ei-si;
	Nl = ej-sj;
	Pl = ek-sk;
	ierr = PetscMalloc( sizeof(PetscInt)*Ml*Nl*Pl*3, &fine_indices );CHKERRQ(ierr);

	c = 0;
	for( k=sk; k<ek; k++ ){
		for( j=sj; j<ej; j++ ){
			for( i=si; i<ei; i++ ){
				PetscInt nidx;

				nidx = (i) + (j)*M + (k)*M*N;
				fine_indices[c  ] = 3 * nidx     ;
				fine_indices[c+1] = 3 * nidx + 1 ;
				fine_indices[c+2] = 3 * nidx + 2 ;

				c = c + 3;
			}
		}
	}

	/* generate a local da to store the coords on */
  /* Note, we cannot use
	   ierr = DMDACreate3d( PETSC_COMM_SELF, DMDA_NONPERIODIC, st, (ei-si),(ej-sj),(ek-sk), PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, n_dofs, sw, 0,0,0, &_sda );CHKERRQ(ierr);
	   ierr = DMSetUp(_sda);CHKERRQ(ierr);
	   ierr = DMSetOptionsPrefix(_sda,"redundant_");CHKERRQ(ierr);

	   We have to call this, as flags like -da_processors_z 3 get filtered down into DMSetUp_DA_3D() even when
	   1) DMOptionsPrefix is set
	   2) The comm is serial
	   3) 1 is passed in for the processor size
	*/

	ierr = x_DMDACreate3d( PETSC_COMM_SELF, wrapNP, st, (ei-si),(ej-sj),(ek-sk), PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, n_dofs, sw, 0,0,0, &_sda );CHKERRQ(ierr);
	/* more hacky shit due to DMSetFromOptions() in constructor */
  /* Note: DMSetFromOptions() was removed in PETSc 3.8 so this hack could be cleaned up */
	{
		DM dd_da_coordinates;
		PetscInt _m,_n,_p,_M,_N,_P;
		ierr = DMDAGetInfo(_sda,0,&_m,&_n,&_p,&_M,&_N,&_P,0,0,0,0,0,0);CHKERRQ(ierr);
		ierr = x_DMDACreate3d(PetscObjectComm((PetscObject)_sda),wrapNP,st,_m,_n,_p,_M,_N,_P,3,sw,0,0,0,&dd_da_coordinates);CHKERRQ(ierr);
		ierr = DMSetCoordinateDM(_sda,dd_da_coordinates);CHKERRQ(ierr);
        ierr = DMDestroy(&dd_da_coordinates);CHKERRQ(ierr); /* hand back destruction to the _sda object */
	}
	ierr = DMDASetUniformCoordinates( _sda, 0.0,1.0, 0.0,1.0, 0.0,1.0 );CHKERRQ(ierr);

	ierr = DMGetCoordinates(_sda,&new_local_coords);CHKERRQ(ierr);

	/* generate scatter */
	ierr = ISCreateGeneral( PetscObjectComm((PetscObject)da), Ml*Nl*Pl*3, fine_indices, PETSC_USE_POINTER, &is_fine );CHKERRQ(ierr);
	ierr = ISCreateStride( PETSC_COMM_SELF,Ml*Nl*Pl*3,0,1,&is_local);CHKERRQ(ierr);

	ierr = VecScatterCreate( global_coords_fine_NATURAL,is_fine, new_local_coords,is_local, &sctx );CHKERRQ(ierr);

	ierr = VecScatterBegin( sctx, global_coords_fine_NATURAL,new_local_coords,INSERT_VALUES, SCATTER_FORWARD );CHKERRQ(ierr);
	ierr = VecScatterEnd(   sctx, global_coords_fine_NATURAL,new_local_coords,INSERT_VALUES, SCATTER_FORWARD );CHKERRQ(ierr);

	/* We scattered the parallel coordinates into the seq coords on the new da */
	/* Since the new da is sequential, we don't need to call DMDAUpdateGhostedCoordinates() */

	/* tidy */
	ierr = VecScatterDestroy(&sctx);CHKERRQ(ierr);
	ierr = ISDestroy(&is_fine);CHKERRQ(ierr);
	ierr = PetscFree(fine_indices);CHKERRQ(ierr);
	ierr = ISDestroy(&is_local);CHKERRQ(ierr);
	ierr = VecDestroy(&global_coords_fine_NATURAL);CHKERRQ(ierr);

	*_seq_DA = _sda;

	PetscFunctionReturn(0);
}


PetscErrorCode DMDACreate3dSemiRedundant(DM da,PetscInt nred,PetscMPISubComm *sub,DM *sda)
{
	MPI_Comm comm,subcomm;
	PetscMPISubComm _sub;
	DM _sda,seq_da;
	PetscBool active;
	PetscInt M,N,P,dof,sw;
	DMBoundaryType wrap[3];
	DMDAStencilType st;
	PetscInt si[3],gnx[3];
	PetscErrorCode ierr;

	comm = PetscObjectComm((PetscObject)da);

	ierr = PetscMPISubCommCreate(comm,nred,&_sub);CHKERRQ(ierr);
	ierr = PetscMPISubCommGetComm(_sub,&subcomm);CHKERRQ(ierr);
	ierr = PetscMPISubCommGetActive(_sub,&active);CHKERRQ(ierr);

	/* get properties from original da, create sub da letting petsc determine distribution */
	ierr = DMDAGetInfo( da, 0, &M,&N,&P, 0,0,0, &dof,&sw,&wrap[0],&wrap[1],&wrap[2],&st );CHKERRQ(ierr);
	if (active) {
		ierr = DMDACreate3d(subcomm,wrap[0],wrap[1],wrap[2], st, M,N,P, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, dof,sw, NULL,NULL,NULL, &_sda );CHKERRQ(ierr);
    ierr = DMSetUp(_sda);CHKERRQ(ierr);
	} else {
		_sda = NULL;
	}

	/* fetch distribution information about sub da */
	/* active ranks fetch a range of entries, non active ranks fetch a single local node */
	si[0]  = si[1]  = si[2]  = -1;
	gnx[0] = gnx[1] = gnx[2] = -1;
	if (active) {
		ierr = DMDAGetCorners(_sda,&si[0],&si[1],&si[2],&gnx[0],&gnx[1],&gnx[2]);CHKERRQ(ierr);
	} else {
		ierr = DMDAGetCorners(da,&si[0],&si[1],&si[2],NULL,NULL,NULL);CHKERRQ(ierr);
		gnx[0] = gnx[1] = gnx[2] = 1;
	}

	ierr = DMDACreate3dRedundant(da,si[0],si[0]+gnx[0], si[1],si[1]+gnx[1], si[2],si[2]+gnx[2], dof, &seq_da );CHKERRQ(ierr);

	/* now copy the coordiantes */
	if (active) {
		PetscInt s_si[3],seq_si[3],s_gnx[3],seq_gnx[3];
		DM seq_cda,s_cda;
		Vec seq_coor,s_coor;
		DMDACoor3d ***LA_seq_coor,***LA_s_coor;
		PetscInt i,j,k;

		ierr = DMDASetUniformCoordinates( _sda, 0.0,1.0, 0.0,1.0, 0.0,1.0 );CHKERRQ(ierr);

		ierr = DMGetCoordinateDM(seq_da,&seq_cda);CHKERRQ(ierr);
		ierr = DMGetCoordinates(seq_da,&seq_coor);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(seq_cda,seq_coor,&LA_seq_coor);CHKERRQ(ierr);

		ierr = DMGetCoordinateDM(_sda,&s_cda);CHKERRQ(ierr);
		ierr = DMGetCoordinates(_sda,&s_coor);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(s_cda,s_coor,&LA_s_coor);CHKERRQ(ierr);

		/* check sizes */
		ierr = DMDAGetCorners(_sda,  &s_si[0],  &s_si[1],  &s_si[2],  &s_gnx[0],  &s_gnx[1],  &s_gnx[2]);CHKERRQ(ierr);
		ierr = DMDAGetCorners(seq_da,&seq_si[0],&seq_si[1],&seq_si[2],&seq_gnx[0],&seq_gnx[1],&seq_gnx[2]);CHKERRQ(ierr);

		if (s_gnx[0] != seq_gnx[0]) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"i range doesn't match"); }
		if (s_gnx[1] != seq_gnx[1]) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"j range doesn't match"); }
		if (s_gnx[2] != seq_gnx[2]) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"k range doesn't match"); }

		for (k=s_si[2]; k<s_si[2]+s_gnx[2]; k++) {
			for (j=s_si[1]; j<s_si[1]+s_gnx[1]; j++) {
				for (i=s_si[0]; i<s_si[0]+s_gnx[0]; i++) {

					LA_s_coor[k][j][i].x = LA_seq_coor[ k - s_si[2] ][ j - s_si[1] ][ i - s_si[0] ].x;
					LA_s_coor[k][j][i].y = LA_seq_coor[ k - s_si[2] ][ j - s_si[1] ][ i - s_si[0] ].y;
					LA_s_coor[k][j][i].z = LA_seq_coor[ k - s_si[2] ][ j - s_si[1] ][ i - s_si[0] ].z;

				}
			}
		}

		ierr = DMDAVecRestoreArray(s_cda,s_coor,&LA_s_coor);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(seq_cda,seq_coor,&LA_seq_coor);CHKERRQ(ierr);


		ierr = DMDAUpdateGhostedCoordinates(_sda);CHKERRQ(ierr);
	}

	ierr = DMDestroy(&seq_da);CHKERRQ(ierr);


	*sub = _sub;
	*sda = _sda;

	PetscFunctionReturn(0);
}


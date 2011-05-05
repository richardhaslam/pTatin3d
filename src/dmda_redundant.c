#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "private/daimpl.h"    /*I   "petscdm.h"   I*/
#include "dmda_redundant.h"

extern PetscErrorCode DMView_DA_Private(DM da);

/* This is only needed as the DMDACreate3d() forces a call to DMSetFromOptions().
 This parses the options -da_processors_x,y,z even when the comm size is 1, which screws up my seq surface DM.
 This function is a copy of DMDACreate3d() without the call to DMSetFromOptions()
*/
#undef __FUNCT__
#define __FUNCT__ "x_DMDACreate3d"
PetscErrorCode  x_DMDACreate3d(MPI_Comm comm,DMDAPeriodicType wrap,DMDAStencilType stencil_type,PetscInt M,
														 PetscInt N,PetscInt P,PetscInt m,PetscInt n,PetscInt p,PetscInt dof,PetscInt s,const PetscInt lx[],const PetscInt ly[],const PetscInt lz[],DM *da)
{
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  ierr = DMDACreate(comm, da);CHKERRQ(ierr);
  ierr = DMDASetDim(*da, 3);CHKERRQ(ierr);
  ierr = DMDASetSizes(*da, M, N, P);CHKERRQ(ierr);
  ierr = DMDASetNumProcs(*da, m, n, p);CHKERRQ(ierr);
  ierr = DMDASetPeriodicity(*da, wrap);CHKERRQ(ierr);
  ierr = DMDASetDof(*da, dof);CHKERRQ(ierr);
  ierr = DMDASetStencilType(*da, stencil_type);CHKERRQ(ierr);
  ierr = DMDASetStencilWidth(*da, s);CHKERRQ(ierr);
  ierr = DMDASetOwnershipRanges(*da, lx, ly, lz);CHKERRQ(ierr);
  /* This violates the behavior for other classes, but right now users expect negative dimensions to be handled this way */
/*  ierr = DMSetFromOptions(*da);CHKERRQ(ierr); */
  ierr = DMSetUp(*da);CHKERRQ(ierr);
  ierr = DMView_DA_Private(*da);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
 Given a DA (possibly distributed, we form a redundant DMDA which spans 
 the nodes defined via (si,ei), (sj,ej) and (sk,ek).
 
 The nodes which are retrieved are:
 si , si + 1, si + 2, ... ei-1 
 */
#undef __FUNCT__
#define __FUNCT__ "DMDACreate3dRedundant"
PetscErrorCode DMDACreate3dRedundant(DM da,
																	 PetscInt si, PetscInt ei, PetscInt sj, PetscInt ej, PetscInt sk, PetscInt ek,
																	 PetscInt n_dofs,
																	 DM *_seq_DA )
{
	PetscInt M,N,P,dof;
	PetscInt sw;
	DMDAPeriodicType wrap;
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
	
	
	ierr = DMDAGetInfo( da, 0, &M,&N,&P, 0,0,0, &dof,&sw,&wrap,&st );CHKERRQ(ierr);
	/* check ranges */
	if( si < 0 || ei > M ) {
		SETERRQ( ((PetscObject)(da))->comm, PETSC_ERR_USER, "Specified range of nodes in i-direction to gather is too large" );
	}
	if( sj < 0 || ej > N ) {
		SETERRQ( ((PetscObject)(da))->comm, PETSC_ERR_USER, "Specified range of nodes in j-direction to gather is too large" );
	}
	if( sk < 0 || ek > P ) {
		SETERRQ( ((PetscObject)(da))->comm, PETSC_ERR_USER, "Specified range of nodes in k-direction to gather is too large" );
	}

	/* check for cases where user wanted a plane */
	if (si == ei) {
		SETERRQ( ((PetscObject)(da))->comm, PETSC_ERR_USER, "Specified range of nodes in i-direction too. If you wanted a plane, indicate si,si+1 as your range" );
	}
	if (sj == ej) {
		SETERRQ( ((PetscObject)(da))->comm, PETSC_ERR_USER, "Specified range of nodes in j-direction too. If you wanted a plane, indicate sj,sj+1 as your range" );
	}
	if (sk == ek) {
		SETERRQ( ((PetscObject)(da))->comm, PETSC_ERR_USER, "Specified range of nodes in k-direction too. If you wanted a plane, indicate sk,sk+1 as your range" );
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
	ierr = DMDAGetCoordinateDA(da,&cda_fine);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&global_coords_fine);CHKERRQ(ierr); /* a global vector */
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
#if 0
	ierr = DMDACreate3d( PETSC_COMM_SELF, DMDA_NONPERIODIC, st, (ei-si),(ej-sj),(ek-sk), PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, n_dofs, sw, 0,0,0, &_sda );CHKERRQ(ierr);
	/* this doesn't do anything :( */
	ierr = DMSetOptionsPrefix(_sda,"redundant_");CHKERRQ(ierr);
	/* we have to call this, as flags like -da_processors_z 3 get filtered down into DMSetUp_DA_3D() even when 
	 1) DMOptionsPrefix is set
	 2) The comm is serial
	 3) 1 is passed in for the processor size
	*/
#endif
	
	ierr = x_DMDACreate3d( PETSC_COMM_SELF, DMDA_NONPERIODIC, st, (ei-si),(ej-sj),(ek-sk), PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, n_dofs, sw, 0,0,0, &_sda );CHKERRQ(ierr);
	/* more hacky shit due to DMSetFromOptions() in constructor */
	{
		DM_DA          *dd = (DM_DA*)_sda->data;
		DM dd_da_coordinates;
		PetscInt _m,_n,_p,_M,_N,_P;
		ierr = DMDAGetInfo(_sda,0,&_m,&_n,&_p,&_M,&_N,&_P,0,0,0,0);CHKERRQ(ierr);
		ierr = x_DMDACreate3d(((PetscObject)_sda)->comm,st,DMDA_STENCIL_BOX,_m,_n,_p,_M,_N,_P,3,sw,0,0,0,&dd_da_coordinates);CHKERRQ(ierr);
		dd->da_coordinates = dd_da_coordinates;
	
	}
	ierr = DMDASetUniformCoordinates( _sda, 0.0,1.0, 0.0,1.0, 0.0,1.0 );CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinates(_sda,&new_local_coords);CHKERRQ(ierr);
	
	/* generate scatter */
	ierr = ISCreateGeneral( PETSC_COMM_WORLD, Ml*Nl*Pl*3, fine_indices, PETSC_USE_POINTER, &is_fine );CHKERRQ(ierr);
	ierr = ISCreateStride( PETSC_COMM_SELF,Ml*Nl*Pl*3,0,1,&is_local);CHKERRQ(ierr);
	
	ierr = VecScatterCreate( global_coords_fine_NATURAL,is_fine, new_local_coords,is_local, &sctx );CHKERRQ(ierr);
	
	ierr = VecScatterBegin( sctx, global_coords_fine_NATURAL,new_local_coords,INSERT_VALUES, SCATTER_FORWARD );CHKERRQ(ierr);
	ierr = VecScatterEnd(   sctx, global_coords_fine_NATURAL,new_local_coords,INSERT_VALUES, SCATTER_FORWARD );CHKERRQ(ierr);
	
	/* We scattered the parallel coordinates into the seq coords on the new da */
	/* Since the new da is sequential, we don't need to call DMDAUpdateGhostedCoordinates() */
	
	/* tidy */
	ierr = VecScatterDestroy(sctx);CHKERRQ(ierr);
	ierr = ISDestroy(is_fine);CHKERRQ(ierr);
	ierr = PetscFree(fine_indices);CHKERRQ(ierr);
	ierr = ISDestroy(is_local);CHKERRQ(ierr);
	ierr = VecDestroy(global_coords_fine_NATURAL);CHKERRQ(ierr);
	
	*_seq_DA = _sda;
	
	PetscFunctionReturn(0);
}

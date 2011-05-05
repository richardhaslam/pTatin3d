
#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_compare.h"


#define __DMDACompareStructures(var,v1,v2)PetscSynchronizedPrintf( comm, "[%d]: DMDA->%s are different { %d , %d } \n", rank, (var), (int)((v1)),(int)((v2)) )

#undef __FUNCT__
#define __FUNCT__ "DMDACompareStructures"
PetscErrorCode DMDACompareStructures(DM da1,DM da2,PetscBool *flg)
{
	PetscInt si1,sj1,sk1 , si2,sj2,sk2;
	PetscInt mx1,my1,mz1 , mx2,my2,mz2;
	PetscInt M1,N1,P1 , M2,N2,P2;
	PetscInt cx1,cy1,cz1 , cx2,cy2,cz2;
	PetscInt dim1 , dim2;
	PetscInt sw1 , sw2;
	DMDAPeriodicType wrap1 , wrap2;
	DMDAStencilType st1 , st2;
	MPI_Comm comm;
	PetscMPIInt rank;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	*flg = PETSC_TRUE;

	ierr = DMDAGetInfo( da1, &dim1, &M1,&N1,&P1, &cx1,&cy1,&cz1, 0, &sw1, &wrap1, &st1 );CHKERRQ(ierr);
	ierr = DMDAGetInfo( da2, &dim2, &M2,&N2,&P2, &cx2,&cy2,&cz2, 0, &sw2, &wrap2, &st2 );CHKERRQ(ierr);
	
	ierr = DMDAGetCorners( da1, &si1,&sj1,&sk1 , &mx1,&my1,&mz1 );CHKERRQ(ierr);
	ierr = DMDAGetCorners( da2, &si2,&sj2,&sk2 , &mx2,&my2,&mz2 );CHKERRQ(ierr);
	
	
	PetscObjectGetComm( (PetscObject)da1, &comm );
	MPI_Comm_rank( comm, &rank );
	
	if( dim1 != dim2 ) {	PetscSynchronizedPrintf( comm, "[%d]: DA->dim are different { %d , %d } \n", dim1,dim2 );		*flg = PETSC_FALSE;	}
	
	if( M1 != M2 ) {	__DMDACompareStructures( "M",M1,M2 );		*flg = PETSC_FALSE;	}
	if( N1 != N2 ) {	__DMDACompareStructures( "N",N1,N2 );		*flg = PETSC_FALSE;	}
	if( P1 != P2 ) {	__DMDACompareStructures( "P",P1,P2 );		*flg = PETSC_FALSE;	}
	
	if( cx1 != cx2 ) {	__DMDACompareStructures( "px",cx1,cx2 );		*flg = PETSC_FALSE;	}
	if( cy1 != cy2 ) {	__DMDACompareStructures( "py",cy1,cy2 );		*flg = PETSC_FALSE;	}
	if( cz1 != cz2 ) {	__DMDACompareStructures( "pz",cz1,cz2 );		*flg = PETSC_FALSE;	}
	
	if( sw1 != sw2 ) {		__DMDACompareStructures( "stencil_width",sw1,sw2 );		*flg = PETSC_FALSE;	}
	if( wrap1 != wrap2 ) {	__DMDACompareStructures( "wrap",wrap1,wrap2 );			*flg = PETSC_FALSE;	}
	
	if( si1 != si2 ) {	__DMDACompareStructures( "si",si1,si2 );		*flg = PETSC_FALSE;	}
	if( sj1 != sj2 ) {	__DMDACompareStructures( "sj",sj1,sj2 );		*flg = PETSC_FALSE;	}
	if( sk1 != sk2 ) {	__DMDACompareStructures( "sk",sk1,sk2 );		*flg = PETSC_FALSE;	}
	
	if( mx1 != mx2 ) {	__DMDACompareStructures( "mx",mx1,mx2 );		*flg = PETSC_FALSE;	}
	if( my1 != my2 ) {	__DMDACompareStructures( "my",my1,my2 );		*flg = PETSC_FALSE;	}
	if( mz1 != mz2 ) {	__DMDACompareStructures( "mz",mz1,mz2 );		*flg = PETSC_FALSE;	}
	
	
	if( *flg == PETSC_FALSE ) {
		PetscSynchronizedFlush( comm );
	}
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMDA_CheckNodeIndex2d"
PetscErrorCode DMDA_CheckNodeIndex2d(DM da,PetscBool ghosted,PetscInt i,PetscInt j,PetscInt k )
{
	PetscInt si,sj,nx,ny;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDAGetCorners(da,&si,&sj,0,&nx,&ny,0);CHKERRQ(ierr);
	if(ghosted==PETSC_TRUE) {
		ierr = DMDAGetGhostCorners(da,&si,&sj,0,&nx,&ny,0);CHKERRQ(ierr);
	}
	if( i<si ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "i=%D < start_index_i(%D)", i,si );
	if( j<sj ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "j=%D < start_index_j(%D)", j,sj );
	
	if( i>=(si+nx) ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "i=%D >= end_index_i(%D)", i,si+nx );
	if( j>=(sj+ny) ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "j=%D >= end_index_j(%D)", j,sj+ny );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDA_CheckNodeIndex3d"
PetscErrorCode DMDA_CheckNodeIndex3d(DM da,PetscBool ghosted,PetscInt i,PetscInt j,PetscInt k )
{
	PetscErrorCode ierr;
	PetscInt si,sj,sk,nx,ny,nz;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	if(ghosted==PETSC_TRUE) {
		ierr = DMDAGetGhostCorners(da,&si,&sj,0,&nx,&ny,0);CHKERRQ(ierr);
	}
	if( i<si ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "i=%D < start_index_i(%D)", i,si );
	if( j<sj ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "j=%D < start_index_j(%D)", j,sj );
	if( k<sk ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "k=%D < start_index_k(%D)", k,sk );
	
	if( i>=(si+nx) ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "i=%D >= end_index_i(%D)", i,si+nx );
	if( j>=(sj+ny) ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "j=%D >= end_index_j(%D)", j,sj+ny );
	if( k>=(sk+nz) ) SETERRQ2( ((PetscObject)da)->comm, PETSC_ERR_USER, "k=%D >= end_index_k(%D)", k,sk+nz );
	
	PetscFunctionReturn(0);
}


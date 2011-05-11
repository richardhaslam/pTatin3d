
#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_compare.h"
#include "dmda_update_coords.h"
#include "dmda_redundant.h"
#include "dmda_remesh.h"

#undef __FUNCT__
#define __FUNCT__ "DMDAGetCornerCoordinatesInPlane_IJ"
PetscErrorCode DMDAGetCornerCoordinatesInPlane_IJ(DM da,PetscInt K,DMDACoor3d coords[])
{
	PetscInt ii,jj,si,sj,sk;
	PetscInt MX,MY,MZ,nx,ny,nz;
	PetscReal corners[4][3];
	Vec coordinates;
	MPI_Comm comm;
	DM cda;
	DMDACoor3d ***LA_coords;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&MX,&MY,&MZ, 0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da, &si,&sj,&sk, &nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coordinates);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
	
	comm = ((PetscObject)da)->comm;
	if ( (K>=sk) && (K<sk+nz) ) {
		if ( (si==0) && (sj==0) ) {
			ii = 0;
			jj = 0;
			corners[0][0] = LA_coords[K][jj][ii].x;
			corners[0][1] = LA_coords[K][jj][ii].y;
			corners[0][2] = LA_coords[K][jj][ii].z;
			ierr = MPI_Bcast(corners[0],3,MPIU_REAL,0,comm);CHKERRQ(ierr);
		}}
	
	if ( (K>=sk) && (K<sk+nz) ) {
		if ( (si==0) && (sj==MY-1) ) {
			ii = 0;
			jj = MY-1;
			corners[1][0] = LA_coords[K][jj][ii].x;
			corners[1][1] = LA_coords[K][jj][ii].y;
			corners[1][2] = LA_coords[K][jj][ii].z;
			ierr = MPI_Bcast(corners[1],3,MPIU_REAL,0,comm);CHKERRQ(ierr);
		}}
	
	if ( (K>=sk) && (K<sk+nz) ) {
		if ( (si==MX-1) && (sj==MY-1) ) {
			ii = MX-1;
			jj = MY-1;
			corners[2][0] = LA_coords[K][jj][ii].x;
			corners[2][1] = LA_coords[K][jj][ii].y;
			corners[2][2] = LA_coords[K][jj][ii].z;
			ierr = MPI_Bcast(corners[2],3,MPIU_REAL,0,comm);CHKERRQ(ierr);
		}}
	
	if ( (K>=sk) && (K<sk+nz) ) {
		if ( (si==MX-1) && (sj==0) ) {
			ii = MX-1;
			jj = 0;
			corners[3][0] = LA_coords[K][jj][ii].x;
			corners[3][1] = LA_coords[K][jj][ii].y;
			corners[3][2] = LA_coords[K][jj][ii].z;
			ierr = MPI_Bcast(corners[3],3,MPIU_REAL,0,comm);CHKERRQ(ierr);
		}}
	
	ierr = DMDAVecRestoreArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
	
	coords[0].x = corners[0][0];    coords[0].y = corners[0][1];    	coords[0].z = corners[0][2];
	coords[1].x = corners[1][0];    coords[1].y = corners[1][1];    	coords[1].z = corners[1][2];
	coords[2].x = corners[2][0];    coords[2].y = corners[2][1];    	coords[2].z = corners[2][2];
	coords[3].x = corners[3][0];    coords[3].y = corners[3][1];    	coords[3].z = corners[3][2];
	
	PetscFunctionReturn(0);
}

/*
 Remesh all nodes in plane K, in a regular fashion.
 coords[] defines the vertices of the bilinear plane which defines the surface we remesh over.
 The order should be
 (i,j+1)[1]  (i+1,j+1)[2]
 (i,j)  [0]  (i+1,j)  [3]
*/
#undef __FUNCT__
#define __FUNCT__ "DMDARemeshSetUniformCoordinatesInPlane_IJ"
PetscErrorCode DMDARemeshSetUniformCoordinatesInPlane_IJ(DM da,PetscInt K,DMDACoor3d coords[])
{
	PetscErrorCode ierr;
	PetscInt si,sj,sk,nx,ny,nz,i,j;
	PetscInt MX,MY;
	PetscInt n;
	PetscReal dxi,deta,Ni[4];
	DM cda;
	Vec coord;
	DMDACoor3d ***_coord;

	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&MX,&MY,0, 0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&_coord);CHKERRQ(ierr);
	
	if ( (K>=sk) && (K<sk+nz) ) {
		dxi = 2.0/(PetscReal)(MX-1);
		deta = 2.0/(PetscReal)(MY-1);
		for( j=sj; j<sj+ny; j++ ) {
			for( i=si; i<si+nx; i++ ) {
				PetscReal xi,eta, xn,yn,zn;
				
				xi  = -1.0 + i*dxi;
				eta = -1.0 + j*deta;
				
				Ni[0] = 0.25 * (1.0-xi) * (1.0-eta);
				Ni[1] = 0.25 * (1.0-xi) * (1.0+eta);
				Ni[2] = 0.25 * (1.0+xi) * (1.0+eta);
				Ni[3] = 0.25 * (1.0+xi) * (1.0-eta);
				
				xn = yn = zn = 0.0;
				for (n=0; n<4; n++) {
					xn = xn + Ni[n] * coords[n].x;
					yn = yn + Ni[n] * coords[n].y;
					zn = zn + Ni[n] * coords[n].z;
				}
				_coord[K][j][i].x = xn;
				_coord[K][j][i].y = yn;
				_coord[K][j][i].z = zn;
			}
		}
	}
	/* tidy up */
	ierr = DMDAVecRestoreArray(cda,coord,&_coord);CHKERRQ(ierr);

	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMDARemeshSetUniformCoordinatesBetweenKLayers3d_MPI"
PetscErrorCode DMDARemeshSetUniformCoordinatesBetweenKLayers3d_MPI( DM da, PetscInt startK, PetscInt endK )
{
	PetscInt si,sj,sk,nx,ny,nz,i,j,k;
	PetscInt s_si,s_sj,s_sk,s_nx,s_ny,s_nz;
	DM surface1_da,surface2_da;
	DM surface1_cda,surface2_cda;
	Vec surface1_coords,surface2_coords;
	PetscScalar *surface1_nodes,*surface2_nodes;
	DM cda;
	Vec coords;
	DMDACoor3d ***nodes;
	DMDACoor3d DX;
	PetscInt start,end, DL, RANGE;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* build top and bottom surface da on my processor */
	ierr = DMDAGetCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
	ierr = DMDACreate3dRedundant( da, si,si+nx, sj,sj+ny, startK,  startK+1, 1, &surface1_da );CHKERRQ(ierr);
	ierr = DMDACreate3dRedundant( da, si,si+nx, sj,sj+ny, endK-1,  endK,     1, &surface2_da );CHKERRQ(ierr);
	
	/*
	{
		PetscViewer vv;
		Vec x;
		ierr = PetscViewerASCIIOpen(((PetscObject)(surface1_da))->comm, "test_dmda_remesh_s1.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(surface1_da,&x);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
		ierr = DMView(surface1_da, vv);CHKERRQ(ierr);
		ierr = VecView(x, vv);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(vv);CHKERRQ(ierr);
		ierr  = VecDestroy(x);CHKERRQ(ierr);
	}
	{
		PetscViewer vv;
		Vec x;
		ierr = PetscViewerASCIIOpen(((PetscObject)(surface2_da))->comm, "test_dmda_remesh_s2.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(surface2_da,&x);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
		ierr = DMView(surface2_da, vv);CHKERRQ(ierr);
		ierr = VecView(x, vv);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(vv);CHKERRQ(ierr);
		ierr  = VecDestroy(x);CHKERRQ(ierr);
	}
	*/
	

//	ierr = DMDAGetCoordinateDA( surface1_da, &surface1_cda);CHKERRQ(ierr); /* don't access coordinate da's on these guys */
	ierr = DMDAGetCoordinates( surface1_da,&surface1_coords );CHKERRQ(ierr);
//	ierr = DMDAVecGetArray(surface1_cda,surface1_coords,&surface1_nodes);CHKERRQ(ierr);
	ierr = VecGetArray(surface1_coords,&surface1_nodes);CHKERRQ(ierr);
	
//	ierr = DMDAGetCoordinateDA( surface2_da, &surface2_cda);CHKERRQ(ierr); /* don't access coordinate da's on these guys! */
	ierr = DMDAGetCoordinates( surface2_da,&surface2_coords );CHKERRQ(ierr);
//	ierr = DMDAVecGetArray(surface2_cda,surface2_coords,&surface2_nodes);CHKERRQ(ierr);
	ierr = VecGetArray(surface2_coords,&surface2_nodes);CHKERRQ(ierr);

	
	ierr = DMDAGetCoordinateDA( da, &cda);CHKERRQ(ierr);

	ierr = DMDAGetCoordinates( da,&coords );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coords,&nodes);CHKERRQ(ierr);

	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
//	ierr = DMDAGetCorners(surface2_cda,&s_si,&s_sj,&s_sk,&s_nx,&s_ny,&s_nz);CHKERRQ(ierr);
	ierr = DMDAGetCorners(surface2_da,&s_si,&s_sj,&s_sk,&s_nx,&s_ny,&s_nz);CHKERRQ(ierr);

	/* starts won't match as the the surface meshes are sequential */
/*
	if( si != s_si ) {  SETERRQ( PETSC_COMM_SELF, PETSC_ERR_USER,"si on da must match surface da (s1)" );  }
        if( sj != s_sj ) {  SETERRQ( PETSC_COMM_SELF, PETSC_ERR_USER,"sj on da must match surface da (s1)" );  }
*/
        if( nx != s_nx ) {  SETERRQ( PETSC_COMM_SELF, PETSC_ERR_USER,"nx on da must match surface da (s1)" );  }
        if( ny != s_ny ) {  SETERRQ( PETSC_COMM_SELF, PETSC_ERR_USER,"ny on da must match surface da (s1)" );  }

        if( s_sk != 0  ) {  SETERRQ( PETSC_COMM_SELF, PETSC_ERR_USER,"s_sk on surface da should be 0 (s1)" );  }
        if( s_nz != 1  ) {  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"s_nz on surface da should be 1 (s1)" );  }



#if 0	
	/* Some portion of grid may overlap sub-domains */
	/* Figure out the range of k indices I need to traverse */
	start = 0;
	end   = 0;
	if( (endK+1) >= (sk+nz) ) {
		printf("sk=%d,nz=%d\n",sk,nz);
		end = PetscMin( sk+nz, endK+1 );
	}
	if( (startK) > (sk) ) {
		start = PetscMax( sk, startK );
	}
	DL = end - start;
	if( DL < 0 ) {
		SETERRQ( ((PetscObject)da)->comm, PETSC_ERR_USER, "DL cannot be negative" );
	}
	
	/* total range of k indices to span */
	RANGE = endK - startK;
	
	printf("DL = %d \n", DL );
	
	if( DL!=0 ) {
		
		for( j=sj; j<sj+ny; j++ ) {
			for( i=si; i<si+nx; i++ ) {
				
				DX.x = (  surface2_nodes[0][j][i].x  -  surface1_nodes[0][j][i].x  )/( (PetscScalar)(RANGE) );
				DX.y = (  surface2_nodes[0][j][i].y  -  surface1_nodes[0][j][i].y  )/( (PetscScalar)(RANGE) );
				DX.z = (  surface2_nodes[0][j][i].z  -  surface1_nodes[0][j][i].z  )/( (PetscScalar)(RANGE) );
				
				for( k=start; k<end; k++ ) {
					ierr = DMDA_CheckNodeIndex3d(da,PETSC_FALSE,i,j,k); CHKERRQ(ierr);
					
					nodes[k][j][i].x = surface1_nodes[0+nx*j+i].x + ((PetscScalar)k) * DX.x;
					nodes[k][j][i].y = surface1_nodes[0+nx*j+i].y + ((PetscScalar)k) * DX.y;
					nodes[k][j][i].z = surface1_nodes[0+nx*j+i].z + ((PetscScalar)k) * DX.z;
				}
			}
		}
		
	}
#endif

	RANGE = endK - startK;
	for( j=0; j<ny; j++ ) {
		for( i=0; i<nx; i++ ) {
			
			DX.x = (  surface2_nodes[3*(0+nx*j+i)  ]  -  surface1_nodes[3*(0+nx*j+i)  ]  )/( (PetscScalar)(RANGE-1) );
			DX.y = (  surface2_nodes[3*(0+nx*j+i)+1]  -  surface1_nodes[3*(0+nx*j+i)+1]  )/( (PetscScalar)(RANGE-1) );
			DX.z = (  surface2_nodes[3*(0+nx*j+i)+2]  -  surface1_nodes[3*(0+nx*j+i)+2]  )/( (PetscScalar)(RANGE-1) );

			for( k=sk; k<sk+nz; k++ ) {
				PetscInt ii = i + si;
				PetscInt jj = j + sj;
				ierr = DMDA_CheckNodeIndex3d(da,PETSC_FALSE,ii,jj,k); CHKERRQ(ierr);
				if ( (k>=startK) && (k<endK) ){
					nodes[k][jj][ii].x = surface1_nodes[3*(0+nx*j+i)  ] + ((PetscScalar)k) * DX.x;
					nodes[k][jj][ii].y = surface1_nodes[3*(0+nx*j+i)+1] + ((PetscScalar)k) * DX.y;
					nodes[k][jj][ii].z = surface1_nodes[3*(0+nx*j+i)+2] + ((PetscScalar)k) * DX.z;
				}
			}
		}
	}
	
	/* tidy up */
	ierr = DMDAVecRestoreArray(cda,coords,&nodes);CHKERRQ(ierr);

	ierr = VecRestoreArray(surface1_coords,&surface1_nodes);CHKERRQ(ierr);
	ierr = VecRestoreArray(surface2_coords,&surface2_nodes);CHKERRQ(ierr);

	ierr = DMDestroy(&surface1_da);CHKERRQ(ierr);
	ierr = DMDestroy(&surface2_da);CHKERRQ(ierr);
	
	/* update */
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDARemeshSetUniformCoordinatesBetweenKLayers3d"
PetscErrorCode DMDARemeshSetUniformCoordinatesBetweenKLayers3d( DM da, PetscInt startK, PetscInt endK )
{
	MPI_Comm comm;
	PetscMPIInt size;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_size( comm, &size );
	if(size==1) {
		ierr = DMDARemeshSetUniformCoordinatesBetweenKLayers3d_MPI( da, startK, endK );CHKERRQ(ierr);
	}
	else {
		ierr = DMDARemeshSetUniformCoordinatesBetweenKLayers3d_MPI( da, startK, endK );CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}




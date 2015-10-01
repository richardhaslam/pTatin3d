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
 **    filename:   dmda_remesh.c
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

#include "dmda_compare.h"
#include "dmda_update_coords.h"
#include "dmda_redundant.h"
#include "dmda_remesh.h"
#include "mesh_deformation.h"

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
	ierr = DMGetCoordinates(da,&coordinates);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coordinates,&LA_coords);CHKERRQ(ierr);
	
	comm = PetscObjectComm((PetscObject)da);
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
	
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&coord);CHKERRQ(ierr);
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
#define __FUNCT__ "DMDARemeshSetUniformCoordinatesInPlane_IK"
PetscErrorCode DMDARemeshSetUniformCoordinatesInPlane_IK(DM da,PetscInt J,DMDACoor3d coords[])
{
	PetscErrorCode ierr;
	PetscInt si,sj,sk,nx,ny,nz,i,k;
	PetscInt MX,MZ;
	PetscInt n;
	PetscReal dxi,dzeta,Ni[4];
	DM cda;
	Vec coord;
	DMDACoor3d ***_coord;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&MX,0,&MZ, 0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&_coord);CHKERRQ(ierr);
	
	if ( (J>=sj) && (J<sj+ny) ) {
		dxi   = 2.0/(PetscReal)(MX-1);
		dzeta = 2.0/(PetscReal)(MZ-1);
		for( k=sk; k<sk+nz; k++ ) {
			for( i=si; i<si+nx; i++ ) {
				PetscReal xi,zeta, xn,yn,zn;
				
				xi   = -1.0 + i*dxi;
				zeta = -1.0 + k*dzeta;
				
				Ni[0] = 0.25 * (1.0-xi) * (1.0-zeta);
				Ni[1] = 0.25 * (1.0-xi) * (1.0+zeta);
				Ni[2] = 0.25 * (1.0+xi) * (1.0+zeta);
				Ni[3] = 0.25 * (1.0+xi) * (1.0-zeta);
				
				xn = yn = zn = 0.0;
				for (n=0; n<4; n++) {
					xn = xn + Ni[n] * coords[n].x;
					yn = yn + Ni[n] * coords[n].y;
					zn = zn + Ni[n] * coords[n].z;
				}
				_coord[k][J][i].x = xn;
				_coord[k][J][i].y = yn;
				_coord[k][J][i].z = zn;
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
	Vec surface1_coords,surface2_coords;
	PetscScalar *surface1_nodes,*surface2_nodes;
	DM cda;
	Vec coords;
	DMDACoor3d ***nodes;
	DMDACoor3d DX;
	PetscInt RANGE;
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
		ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)surface1_da), "test_dmda_remesh_s1.vtk", &vv);CHKERRQ(ierr);
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
		ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)surface2_da), "test_dmda_remesh_s2.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(surface2_da,&x);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)x, "phi" );CHKERRQ(ierr);
		ierr = DMView(surface2_da, vv);CHKERRQ(ierr);
		ierr = VecView(x, vv);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(vv);CHKERRQ(ierr);
		ierr  = VecDestroy(x);CHKERRQ(ierr);
	}
	*/
	

//	ierr = DMGetCoordinateDM( surface1_da, &surface1_cda);CHKERRQ(ierr); /* don't access coordinate da's on these guys */
	ierr = DMGetCoordinates( surface1_da,&surface1_coords );CHKERRQ(ierr);
//	ierr = DMDAVecGetArray(surface1_cda,surface1_coords,&surface1_nodes);CHKERRQ(ierr);
	ierr = VecGetArray(surface1_coords,&surface1_nodes);CHKERRQ(ierr);
	
//	ierr = DMGetCoordinateDM( surface2_da, &surface2_cda);CHKERRQ(ierr); /* don't access coordinate da's on these guys! */
	ierr = DMGetCoordinates( surface2_da,&surface2_coords );CHKERRQ(ierr);
//	ierr = DMDAVecGetArray(surface2_cda,surface2_coords,&surface2_nodes);CHKERRQ(ierr);
	ierr = VecGetArray(surface2_coords,&surface2_nodes);CHKERRQ(ierr);

	
	ierr = DMGetCoordinateDM( da, &cda);CHKERRQ(ierr);

	ierr = DMGetCoordinates( da,&coords );CHKERRQ(ierr);
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
		SETERRQ( PetscObjectComm((PetscObject)da), PETSC_ERR_USER, "DL cannot be negative" );
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
					nodes[k][jj][ii].x = surface1_nodes[3*(0+nx*j+i)  ] + ((PetscScalar)(k-startK)) * DX.x;
					nodes[k][jj][ii].y = surface1_nodes[3*(0+nx*j+i)+1] + ((PetscScalar)(k-startK)) * DX.y;
					nodes[k][jj][ii].z = surface1_nodes[3*(0+nx*j+i)+2] + ((PetscScalar)(k-startK)) * DX.z;
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
#define __FUNCT__ "DMDARemeshSetUniformCoordinatesBetweenJLayers3d_MPI"
PetscErrorCode DMDARemeshSetUniformCoordinatesBetweenJLayers3d_MPI( DM da, PetscInt startJ, PetscInt endJ )
{
	PetscInt si,sj,sk,nx,ny,nz,i,j,k;
	PetscInt s_si,s_sj,s_sk,s_nx,s_ny,s_nz;
	DM surface1_da,surface2_da;
	Vec surface1_coords,surface2_coords;
	PetscScalar *surface1_nodes,*surface2_nodes;
	DM cda;
	Vec coords;
	DMDACoor3d ***nodes;
	DMDACoor3d DX;
	PetscInt RANGE;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* build top and bottom surface da on my processor */
	ierr = DMDAGetCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
	ierr = DMDACreate3dRedundant( da, si,si+nx, startJ,  startJ+1, sk,sk+nz, 1, &surface1_da );CHKERRQ(ierr);
	ierr = DMDACreate3dRedundant( da, si,si+nx, endJ-1,  endJ,     sk,sk+nz, 1, &surface2_da );CHKERRQ(ierr);
	
	ierr = DMGetCoordinates( surface1_da,&surface1_coords );CHKERRQ(ierr);
	ierr = VecGetArray(surface1_coords,&surface1_nodes);CHKERRQ(ierr);
	
	ierr = DMGetCoordinates( surface2_da,&surface2_coords );CHKERRQ(ierr);
	ierr = VecGetArray(surface2_coords,&surface2_nodes);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM( da, &cda);CHKERRQ(ierr);
	
	ierr = DMGetCoordinates( da,&coords );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coords,&nodes);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMDAGetCorners(surface2_da,&s_si,&s_sj,&s_sk,&s_nx,&s_ny,&s_nz);CHKERRQ(ierr);
	
	/* starts won't match as the the surface meshes are sequential */
	if( nx != s_nx ) {  SETERRQ( PETSC_COMM_SELF, PETSC_ERR_USER,"nx on da must match surface da (s1)" );  }
	if( nz != s_nz ) {  SETERRQ( PETSC_COMM_SELF, PETSC_ERR_USER,"nz on da must match surface da (s1)" );  }
	
	if( s_sj != 0  ) {  SETERRQ( PETSC_COMM_SELF, PETSC_ERR_USER,"s_sj on surface da should be 0 (s1)" );  }
	if( s_ny != 1  ) {  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"s_ny on surface da should be 1 (s1)" );  }
	
	RANGE = endJ - startJ;
	for( k=0; k<nz; k++ ) {
		for( i=0; i<nx; i++ ) {
			
			DX.x = (  surface2_nodes[3*(0+nx*k+i)  ]  -  surface1_nodes[3*(0+nx*k+i)  ]  )/( (PetscScalar)(RANGE-1) );
			DX.y = (  surface2_nodes[3*(0+nx*k+i)+1]  -  surface1_nodes[3*(0+nx*k+i)+1]  )/( (PetscScalar)(RANGE-1) );
			DX.z = (  surface2_nodes[3*(0+nx*k+i)+2]  -  surface1_nodes[3*(0+nx*k+i)+2]  )/( (PetscScalar)(RANGE-1) );
			
			for( j=sj; j<sj+ny; j++ ) {
				PetscInt ii = i + si;
				PetscInt kk = k + sk;
				ierr = DMDA_CheckNodeIndex3d(da,PETSC_FALSE,ii,j,kk); CHKERRQ(ierr);
				if ( (j>=startJ) && (j<endJ) ){
					nodes[kk][j][ii].x = surface1_nodes[3*(0+nx*k+i)  ] + ((PetscScalar)(j-startJ)) * DX.x;
					nodes[kk][j][ii].y = surface1_nodes[3*(0+nx*k+i)+1] + ((PetscScalar)(j-startJ)) * DX.y;
					nodes[kk][j][ii].z = surface1_nodes[3*(0+nx*k+i)+2] + ((PetscScalar)(j-startJ)) * DX.z;
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
	ierr = MPI_Comm_size( comm, &size );CHKERRQ(ierr);
	if(size==1) {
		ierr = DMDARemeshSetUniformCoordinatesBetweenKLayers3d_MPI( da, startK, endK );CHKERRQ(ierr);
	}
	else {
		ierr = DMDARemeshSetUniformCoordinatesBetweenKLayers3d_MPI( da, startK, endK );CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDARemeshSetUniformCoordinatesBetweenJLayers3d"
PetscErrorCode DMDARemeshSetUniformCoordinatesBetweenJLayers3d( DM da, PetscInt startJ, PetscInt endJ )
{
	MPI_Comm comm;
	PetscMPIInt size;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscObjectGetComm( (PetscObject)da, &comm );
	ierr = MPI_Comm_size( comm, &size );CHKERRQ(ierr);
	if(size==1) {
		ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d_MPI( da, startJ, endJ );CHKERRQ(ierr);
	}
	else {
		ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d_MPI( da, startJ, endJ );CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDARemeshJMAX_UpdateHeightsFromInterior"
PetscErrorCode DMDARemeshJMAX_UpdateHeightsFromInterior(DM da)
{
	PetscErrorCode ierr;
	PetscInt M,N,P,si,sj,sk,nx,ny,nz,i,k;
	DM cda;
	Vec coords,gcoords;
	DMDACoor3d ***LA_coords,***LA_gcoords;
	
	PetscFunctionBegin;

	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&coords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coords,&LA_coords);CHKERRQ(ierr);

	ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	if (sj+ny == N) { /* we contain the surface mesh */
		
		/* check IMIN */
		if (si == 0) {
			i = 0;
			for (k=sk; k<sk+nz; k++) {
				LA_coords[k][N-1][i].y = LA_gcoords[k][N-1][i+1].y;
			}
		}
		
		/* check IMAX */
		if (si + nx == M) {
			i = M-1;
			for (k=sk; k<sk+nz; k++) {
				LA_coords[k][N-1][i].y = LA_gcoords[k][N-1][i-1].y;
			}
		}
		
		/* check KMIN */
		if (sk == 0) {
			k = 0;
			for (i=si; i<si+nx; i++) {
				LA_coords[k][N-1][i].y = LA_gcoords[k+1][N-1][i].y;
			}
		}
		
		/* check KMAX */
		if (sk + nz == P) {
			k = P-1;
			for (i=si; i<si+nx; i++) {
				LA_coords[k][N-1][i].y = LA_gcoords[k-1][N-1][i].y;
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,coords,&LA_coords);CHKERRQ(ierr);
	
	/* update */
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);

	
	PetscFunctionReturn(0);
}

/*
 PetscInt dir     : Value of 0, 1, or 2 indicating in which diection you wish refiniment to occur along 
 PetscReal factor : Controls aggressiveness of refinement in central section of domain. Values larger than one incur refinement.
 PetscReal x1_frac,x2_frac : Define the start end fraction of the sector in the domain here refinement will occur
 Domain is mapped like this
 
 xprime = slope * (x - x_ref) + xprime_ref
 
 */
#undef __FUNCT__
#define __FUNCT__ "DMDASetCoordinatesColumnRefinement"
PetscErrorCode DMDASetCoordinatesColumnRefinement(DM da,PetscInt dir,PetscReal factor,PetscReal x1_frac,PetscReal x2_frac)
{
	PetscErrorCode ierr;
	PetscInt si,sj,sk,nx,ny,nz,i,j,k,M,N,P;
	DM cda,da_min,da_max,cda_min,cda_max;
	Vec coord,coord_min,coord_max;
	DMDACoor3d ***LA_coords,***LA_coords_da_min,***LA_coords_da_max;
	PetscReal x0prime,x1prime,x2prime,x3prime;
	PetscReal x0,x1,x2,x3;
	
	
	PetscFunctionBegin;
	
	if ((dir < 0) || (dir > 3)) {
		SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Value \"dir\" must be one of {0,1,2}");
	}
	if (factor < 1.0) {
		SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Value \"factor\" must be >= 1.0");
	}
	
	ierr = DMDAGetInfo(da,0,&M,&N,&P, 0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	
	switch (dir) {
		case 0:
			ierr = DMDACreate3dRedundant(da,M-1,M,sj,sj+ny,sk,sk+nz, 1, &da_max);CHKERRQ(ierr);
			ierr = DMDACreate3dRedundant(da,0,1,  sj,sj+ny,sk,sk+nz, 1, &da_min);CHKERRQ(ierr);
			break;
		case 1:
			ierr = DMDACreate3dRedundant(da,si,si+nx,N-1,N,sk,sk+nz, 1, &da_max);CHKERRQ(ierr);
			ierr = DMDACreate3dRedundant(da,si,si+nx,0,1,  sk,sk+nz, 1, &da_min);CHKERRQ(ierr);
			break;
		case 2:
			ierr = DMDACreate3dRedundant(da,si,si+nx,sj,sj+ny,P-1,P, 1, &da_max);CHKERRQ(ierr);
			ierr = DMDACreate3dRedundant(da,si,si+nx,sj,sj+ny,0,1,   1, &da_min);CHKERRQ(ierr);
			break;
	}
	
	
	
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coords);CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(da_min,&cda_min);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_min,&coord_min);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_min,coord_min,&LA_coords_da_min);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(da_max,&cda_max);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_max,&coord_max);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_max,coord_max,&LA_coords_da_max);CHKERRQ(ierr);

	/* uniformily set coordinates between [x0prime,x3prime] */
	/*
	switch (dir) {
		case 0:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"dir = I unsupported: Requires DMDARemeshSetUniformCoordinatesBetweenILayers3d()");
			break;
		case 1:
			ierr = DMDARemeshSetUniformCoordinatesBetweenJLayers3d(da,0,N-1);CHKERRQ(ierr);
			break;
		case 2:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"dir = K unsupported: Requires DMDARemeshSetUniformCoordinatesBetweenKLayers3d()");
			break;
	}	
	*/
	
	// ierr = DMDASetUniformCoordinates1D(da,dir,x0prime,x3prime);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (j=sj; j<sj+ny; j++) {
			for (i=si; i<si+nx; i++) {
				PetscScalar pos[3],coord_prime,dl;
				PetscInt ii;
				
				if (dir == 0) {
					/* fetch x0,x3 from surface mesh */
					x0 = LA_coords_da_min[k-sk][j-sj][0].x;
					x3 = LA_coords_da_max[k-sk][j-sj][0].x;
					ii = i;
					dl = (x3prime-x0prime)/((PetscReal)(M-1));
				} else if (dir == 1) {
					/* fetch x0,x3 from surface mesh */
					x0 = LA_coords_da_min[k-sk][0][i-si].y;
					x3 = LA_coords_da_max[k-sk][0][i-si].y;
					ii = j;
					dl = (x3prime-x0prime)/((PetscReal)(N-1));
				} else {
					/* fetch x0,x3 from surface mesh */
					x0 = LA_coords_da_min[0][j-sj][i-si].z;
					x3 = LA_coords_da_max[0][j-sj][i-si].z;
					ii = k;
					dl = (x3prime-x0prime)/((PetscReal)(P-1));
				}
				
				/* compute x1,x2 from x0,x3 */
				x1 = x0 + (x3-x0) * x1_frac;
				x2 = x0 + (x3-x0) * x2_frac;
				
				x0prime = x0;
				x1prime = x1;
				// x = [(x2 - x1)/fac] * (xprime - x1prime) + x1
				// x2prime
				//x2prime = factor + x1prime;
				x2prime = (x2-x1)*factor + x1;
				// x3prime
				x3prime = x3 - x2 + x2prime;
				
				if (dir == 0) {
					/* fetch x0,x3 from surface mesh */
					ii = i;
					dl = (x3prime-x0prime)/((PetscReal)(M-1));
				} else if (dir == 1) {
					/* fetch x0,x3 from surface mesh */
					ii = j;
					dl = (x3prime-x0prime)/((PetscReal)(N-1));
				} else {
					/* fetch x0,x3 from surface mesh */
					ii = k;
					dl = (x3prime-x0prime)/((PetscReal)(P-1));
				}
				
				pos[0] = LA_coords[k][j][i].x;
				pos[1] = LA_coords[k][j][i].y;
				pos[2] = LA_coords[k][j][i].z;
				
				coord_prime = x0prime + ii*dl;
				
				pos[dir] = coord_prime;
				
				LA_coords[k][j][i].x = pos[0];
				LA_coords[k][j][i].y = pos[1];
				LA_coords[k][j][i].z = pos[2];
			}
		}
	}
	
	
	for (k=sk; k<sk+nz; k++) {
		for (j=sj; j<sj+ny; j++) {
			for (i=si; i<si+nx; i++) {
				PetscScalar pos[3],coord_prime,new_coord;

				if (dir == 0) {
					/* fetch x0,x3 from surface mesh */
					x0 = LA_coords_da_min[k-sk][j-sj][0].x;
					x3 = LA_coords_da_max[k-sk][j-sj][0].x;
				} else if (dir == 1) {
					/* fetch x0,x3 from surface mesh */
					x0 = LA_coords_da_min[k-sk][0][i-si].y;
					x3 = LA_coords_da_max[k-sk][0][i-si].y;
				} else {
					/* fetch x0,x3 from surface mesh */
					x0 = LA_coords_da_min[0][j-sj][i-si].z;
					x3 = LA_coords_da_max[0][j-sj][i-si].z;
				}

				/* compute x1,x2 from x0,x3 */
				x1 = x0 + (x3-x0) * x1_frac;
				x2 = x0 + (x3-x0) * x2_frac;
				
				x0prime = x0;
				x1prime = x1;
				// x = [(x2 - x1)/fac] * (xprime - x1prime) + x1
				// x2prime
				//x2prime = factor + x1prime;
				x2prime = (x2-x1)*factor + x1;
				// x3prime
				x3prime = x3 - x2 + x2prime;
				
				
				pos[0] = LA_coords[k][j][i].x;
				pos[1] = LA_coords[k][j][i].y;
				pos[2] = LA_coords[k][j][i].z;
				
				coord_prime = pos[dir];
				
				if (pos[dir] <= x1prime) {
					new_coord = coord_prime;
				} else if (pos[dir] > x1prime && pos[dir] < x2prime) {
					//new_coord = (x2-x1) * (coord_prime - x1prime)/(factor) + x1;
					new_coord = (1.0) * (coord_prime - x1prime)/(factor) + x1;
				} else {
					new_coord = 1.0 * (coord_prime - x2prime) + x2;
				}
				
				pos[dir] = new_coord;
				
				LA_coords[k][j][i].x = pos[0];
				LA_coords[k][j][i].y = pos[1];
				LA_coords[k][j][i].z = pos[2];
			}
		}
	}
	ierr = DMDAVecRestoreArray(cda_max,coord_max,&LA_coords_da_max);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda_min,coord_min,&LA_coords_da_min);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coords);CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDACoordinateRefinementTransferFunction_OrthogonalFaces"
PetscErrorCode DMDACoordinateRefinementTransferFunction_OrthogonalFaces(DM da,PetscInt dir,PetscInt npoints,PetscReal xref[],PetscReal xnatural[])
{
    PetscErrorCode ierr;
    DM cda;
    Vec coord;
    PetscReal gmin[3],gmax[3];
    PetscInt nc,M,N,P,si,sj,sk,nx,ny,nz,i,j,k;
    PetscScalar *LA_coords;
    MPI_Comm comm;
    
    ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
    ierr = DMDAGetBoundingBox(da,gmin,gmax);CHKERRQ(ierr);
    
    if (xnatural[0] < 0.0) SETERRQ1(comm,PETSC_ERR_USER,"xnatural[0] must be >= 0.0 (%1.4e)",xnatural[0]);
    if (xnatural[npoints-1] > 1.0) SETERRQ1(comm,PETSC_ERR_USER,"mesh_xmax > xnatural[last] must be < 1.0 (%1.4e)",xnatural[npoints-1]);
    if ((dir < 0) || (dir > 3)) SETERRQ(comm,PETSC_ERR_SUP,"Direction must be {0,1,2}");

    
    /* reset mesh to be mapped to the reference coordinate system [0,1] */
    ierr = DMDASetUniformCoordinates1D(da,dir,0.0,1.0);CHKERRQ(ierr);
    
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = VecGetArray(coord,&LA_coords);CHKERRQ(ierr);
    
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                PetscReal xc_nat,xc_ref,xnatural0,xnatural1,xref0,xref1;
                PetscInt nid,region;
                
                nid = (i) + (j) * nx + (k) * nx * ny;
                xc_ref = LA_coords[3*nid + dir]; /* dir = 0,1,2 */
                
                region = -1;
                for (nc=0; nc<npoints-1; nc++) {
                    if ((xc_ref >= xref[nc]) && (xc_ref <= xref[nc+1])) {
                        region = nc;
                        break;
                    }
                }
                if (region == -1) {
                    switch (dir) {
                        case 0:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index i=%D",i);
                            break;
                        case 1:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index j=%D",j);
                            break;
                        case 2:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index k=%D",k);
                            break;
                        default:
                            break;
                    }
                }
                
                /* perform interpolation */
                xc_nat = 0.0;
                
                xnatural0 = xnatural[region];
                xnatural1 = xnatural[region+1];
                
                xref0 = xref[region];
                xref1 = xref[region+1];
                
                xc_nat = (xc_ref - xref0) * (xnatural1 - xnatural0)/(xref1 - xref0) + xnatural0;
                
                /* set new coordinate */
                LA_coords[3*nid + dir] = gmin[dir] + xc_nat * (gmax[dir] - gmin[dir]);  /* dir = 0,1,2 */
            }
        }
    }
    
	ierr = VecRestoreArray(coord,&LA_coords);CHKERRQ(ierr);
    
    ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDACoordinateRefinementTransferFunction_PreserveFaceGeometry"
PetscErrorCode DMDACoordinateRefinementTransferFunction_PreserveFaceGeometry(DM da,PetscInt dir,PetscInt npoints,PetscReal xref[],PetscReal xnatural[])
{
    PetscErrorCode ierr;
    DM da_min,da_max,cda,cda_min,cda_max;
	Vec coord,coord_min,coord_max;
	DMDACoor3d ***LA_coords_da_min,***LA_coords_da_max;
    PetscReal gmin[3],gmax[3];
    PetscInt nc,M,N,P,si,sj,sk,nx,ny,nz,i,j,k;
    PetscScalar *LA_coords;
    MPI_Comm comm;
    
    ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
    ierr = DMDAGetBoundingBox(da,gmin,gmax);CHKERRQ(ierr);
    
    if (xnatural[0] < 0.0) SETERRQ1(comm,PETSC_ERR_USER,"xnatural[0] must be >= 0.0 (%1.4e)",xnatural[0]);
    if (xnatural[npoints-1] > 1.0) SETERRQ1(comm,PETSC_ERR_USER,"mesh_xmax > xnatural[last] must be < 1.0 (%1.4e)",xnatural[npoints-1]);
    
	switch (dir) {
		case 0:
			ierr = DMDACreate3dRedundant(da,M-1,M,sj,sj+ny,sk,sk+nz, 1, &da_max);CHKERRQ(ierr);
			ierr = DMDACreate3dRedundant(da,0,1,  sj,sj+ny,sk,sk+nz, 1, &da_min);CHKERRQ(ierr);
			break;
		case 1:
			ierr = DMDACreate3dRedundant(da,si,si+nx,N-1,N,sk,sk+nz, 1, &da_max);CHKERRQ(ierr);
			ierr = DMDACreate3dRedundant(da,si,si+nx,0,1,  sk,sk+nz, 1, &da_min);CHKERRQ(ierr);
			break;
		case 2:
			ierr = DMDACreate3dRedundant(da,si,si+nx,sj,sj+ny,P-1,P, 1, &da_max);CHKERRQ(ierr);
			ierr = DMDACreate3dRedundant(da,si,si+nx,sj,sj+ny,0,1,   1, &da_min);CHKERRQ(ierr);
			break;
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Direction must be {0,1,2}");
            break;
	}
    
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = VecGetArray(coord,&LA_coords);CHKERRQ(ierr);
    
	ierr = DMGetCoordinateDM(da_min,&cda_min);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_min,&coord_min);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_min,coord_min,&LA_coords_da_min);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(da_max,&cda_max);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_max,&coord_max);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_max,coord_max,&LA_coords_da_max);CHKERRQ(ierr);
    
    /* define uniform coordinate spacing in direction */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                PetscReal gmin_col,ds;
                PetscInt  nid;
                
                nid = (i) + (j) * nx + (k) * nx * ny;

                switch (dir) {
                    case 0:
                        gmin_col = LA_coords_da_min[k][j][0].x;
                        ds = (LA_coords_da_max[k][j][0].x - LA_coords_da_min[k][j][0].x)/((PetscReal)M);
                        LA_coords[3*nid + dir] = gmin_col + (i+si)*ds;
                        break;
                    case 1:
                        gmin_col = LA_coords_da_min[k][0][i].y;
                        ds = (LA_coords_da_max[k][0][i].y - LA_coords_da_min[k][0][i].y)/((PetscReal)N);
                        LA_coords[3*nid + dir] = gmin_col + (j+sj)*ds;
                        break;
                    case 2:
                        gmin_col = LA_coords_da_min[0][j][i].z;
                        ds = (LA_coords_da_max[0][j][i].z - LA_coords_da_min[0][j][i].z)/((PetscReal)P);
                        LA_coords[3*nid + dir] = gmin_col + (k+sk)*ds;
                        break;
                    default:
                        gmin_col = 0.0;
                        ds = 0.0;
                        break;
                }
            }
        }
    }
    
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                PetscReal xc_nat,xc_ref,xnatural0,xnatural1,xref0,xref1,gmin_col,gmax_gmin_col;
                PetscInt  nid,region;
                
                nid = (i) + (j) * nx + (k) * nx * ny;

                switch (dir) {
                    case 0:
                        gmin_col      = LA_coords_da_min[k][j][0].x;
                        gmax_gmin_col = LA_coords_da_max[k][j][0].x - LA_coords_da_min[k][j][0].x;
                        break;
                    case 1:
                        gmin_col      = LA_coords_da_min[k][0][i].y;
                        gmax_gmin_col = LA_coords_da_max[k][0][i].y - LA_coords_da_min[k][0][i].y;
                        break;
                    case 2:
                        gmin_col      = LA_coords_da_min[0][j][i].z;
                        gmax_gmin_col = LA_coords_da_max[0][j][i].z - LA_coords_da_min[0][j][i].z;
                        break;
                    default:
                        gmin_col = 0.0;
                        gmax_gmin_col = 0.0;
                        break;
                }
                
                /* normalize mesh coord(dir) to be mapped to the reference coordinate system [0,1] */
                xc_ref = (LA_coords[3*nid + dir] - gmin_col)/gmax_gmin_col; /* dir = 0,1,2 */
                
                region = -1;
                for (nc=0; nc<npoints-1; nc++) {
                    if ((xc_ref >= xref[nc]) && (xc_ref <= xref[nc+1])) {
                        region = nc;
                        break;
                    }
                }
                if (region == -1) {
                    switch (dir) {
                        case 0:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index i=%D",i);
                            break;
                        case 1:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index j=%D",j);
                            break;
                        case 2:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index k=%D",k);
                            break;
                        default:
                            break;
                    }
                }
                
                /* perform interpolation */
                xc_nat = 0.0;
                
                xnatural0 = xnatural[region];
                xnatural1 = xnatural[region+1];
                
                xref0 = xref[region];
                xref1 = xref[region+1];
                
                /* set new coordinate */
                xc_nat = (xc_ref - xref0) * (xnatural1 - xnatural0)/(xref1 - xref0) + xnatural0;
                LA_coords[3*nid + dir] = gmin_col + xc_nat * (gmax_gmin_col);  /* dir = 0,1,2 */
            }
        }
    }
    
	ierr = VecRestoreArray(coord,&LA_coords);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(cda_max,coord_max,&LA_coords_da_max);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda_min,coord_min,&LA_coords_da_min);CHKERRQ(ierr);
    
    ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);

    ierr = DMDestroy(&da_min);CHKERRQ(ierr);
    ierr = DMDestroy(&da_max);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

/*
 Maps xref[] --> xnatural[] in any i,j,k direction
 - xref[] represents a point distribution
 - xnatural[] represents the graduated point distribution
 
 A reference coordinate system is used to allow using the mesh graduating functions.
 Thus, each xref[i] and xnatural[i] must be in the range [0,1].
 This formulation allows one to defined arbitrary mesh space under the following conditions:
    (i) for any domain size
   (ii) for case when you wish to preserve an existing geometry on two faces
 To utilize case (ii), set preserve_face_geometry = PETSC_TRUE
*/
#undef __FUNCT__
#define __FUNCT__ "DMDACoordinateRefinementTransferFunction"
PetscErrorCode DMDACoordinateRefinementTransferFunction(DM da,PetscInt dir,PetscBool preserve_face_geometry,PetscInt npoints,PetscReal xref[],PetscReal xnatural[])
{
    PetscErrorCode ierr;

    if (!preserve_face_geometry) {
        ierr = DMDACoordinateRefinementTransferFunction_OrthogonalFaces(da,dir,npoints,xref,xnatural);CHKERRQ(ierr);
    } else {
        ierr = DMDACoordinateRefinementTransferFunction_PreserveFaceGeometry(da,dir,npoints,xref,xnatural);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

/* Depreciated */
#undef __FUNCT__
#define __FUNCT__ "_DMDACoordinateRefinementTransferFunction"
PetscErrorCode _DMDACoordinateRefinementTransferFunction(DM da,PetscInt dir,PetscInt npoints,PetscReal xref[],PetscReal xnatural[])
{
    PetscErrorCode ierr;
    DM cda;
    Vec coord;
    PetscReal gmin[3],gmax[3];
    PetscInt nc,M,N,P,si,sj,sk,nx,ny,nz,i,j,k;
    PetscScalar *LA_coords;
    MPI_Comm comm;
    
    ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
    ierr = DMDAGetBoundingBox(da,gmin,gmax);CHKERRQ(ierr);
    
    switch (dir) {
        case 0:
            if (gmin[0] < xnatural[0]) SETERRQ2(comm,PETSC_ERR_USER,"mesh_xmin[0] < xnatural[0] (%1.4e < %1.4e)",gmin[0],xnatural[0]);
            if (gmax[0] > xnatural[npoints-1]) SETERRQ2(comm,PETSC_ERR_USER,"mesh_xmax > xnatural[last] (%1.4e > %1.4e)",gmax[0],xnatural[npoints-1]);
            break;
        case 1:
            if (gmin[1] < xnatural[0]) SETERRQ2(comm,PETSC_ERR_USER,"mesh_xmin[0] < xnatural[0] (%1.4e < %1.4e)",gmin[1],xnatural[0]);
            if (gmax[1] > xnatural[npoints-1]) SETERRQ2(comm,PETSC_ERR_USER,"mesh_xmax > xnatural[last] (%1.4e > %1.4e)",gmax[1],xnatural[npoints-1]);
            break;
        case 2:
            if (gmin[2] < xnatural[0]) SETERRQ2(comm,PETSC_ERR_USER,"mesh_xmin[0] < xnatural[0] (%1.4e < %1.4e)",gmin[2],xnatural[0]);
            if (gmax[2] > xnatural[npoints-1]) SETERRQ2(comm,PETSC_ERR_USER,"mesh_xmax > xnatural[last] (%1.4e > %1.4e)",gmax[2],xnatural[npoints-1]);
            break;
        default:
            SETERRQ(comm,PETSC_ERR_SUP,"Default direction no supported. dir must be [0,1,2]");
            break;
    }
    
    /* reset mesh to be mapped to the reference coordinate system [0,1] */
    ierr = DMDASetUniformCoordinates1D(da,dir,0.0,1.0);CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&coord);CHKERRQ(ierr);
	ierr = VecGetArray(coord,&LA_coords);CHKERRQ(ierr);
    
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                PetscReal xc_nat,xc_ref,xnatural0,xnatural1,xref0,xref1;
                PetscInt nid,region;
                
                nid = (i) + (j) * nx + (k) * nx * ny;
                xc_ref = LA_coords[3*nid + dir]; /* dir = 0,1,2 */
                
                region = -1;
                for (nc=0; nc<npoints-1; nc++) {
                    if ((xc_ref >= xref[nc]) && (xc_ref <= xref[nc+1])) {
                        region = nc;
                        break;
                    }
                }
                if (region == -1) {
                    switch (dir) {
                        case 0:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index i=%D",i);
                            break;
                        case 1:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index j=%D",j);
                            break;
                        case 2:
                            SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to determine which sector contains node index k=%D",k);
                            break;
                    }
                }
                
                /* perform interpolation */
                xc_nat = 0.0;
                
                xnatural0 = xnatural[region];
                xnatural1 = xnatural[region+1];

                xref0 = xref[region];
                xref1 = xref[region+1];
                
                xc_nat = (xc_ref - xref0) * (xnatural1 - xnatural0)/(xref1 - xref0) + xnatural0;
                
                /* set new coordinate */
                LA_coords[3*nid + dir] = xc_nat;  /* dir = 0,1,2 */
            }
        }
    }

	ierr = VecRestoreArray(coord,&LA_coords);CHKERRQ(ierr);
    
    ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


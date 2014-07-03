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
 **    Filename:      dmda_element_q2p1.c
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


#include "petsc.h"
#include "petscdm.h"
#include "ptatin3d_defs.h"
#include "dmda_element_q2p1.h"


/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetSizeElementQ2"
PetscErrorCode DMDAGetSizeElementQ2(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ)
{
	const PetscInt order = 2;
	PetscInt M,N,P,mx,my,mz,width;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width!=2) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 2 for Q2");
	}
	
	/* M = order*mx+1 */
	if ( (M-1)%order != 0 ) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"DMDA is not compatible with Q2 elements in x direction");
	}
	if (MX) { *MX = (M-1)/order; }
	
	//
	if ( (N-1)%order != 0 ) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"DMDA is not compatible with Q2 elements in y direction");
	}
	if (MY) { *MY = (N-1)/order; }
	
	//
	if ( (P-1)%order != 0 ) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"DMDA is not compatible with Q2 elements in z direction");
	}
	if (MZ) { *MZ = (P-1)/order; }
	
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetLocalSizeElementQ2"
PetscErrorCode DMDAGetLocalSizeElementQ2(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	const PetscInt order = 2;
	PetscInt i,j,k,start;
	PetscInt cntx,cnty,cntz;
	PetscInt si,sj,sk,m,n,p,M,N,P,width;
	PetscInt sig,sjg,skg,mg,ng,pg;
	PetscErrorCode ierr;
	PetscMPIInt rank;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sig,&sjg,&skg,&mg,&ng,&pg);CHKERRQ(ierr);
	if (width!=2) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 2 for Q2");
	}
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	//printf("[%d]: i(%d->%d) : j(%d->%d) \n", rank,si,si+m,sj,sj+n);
	
	cntx = cnty = cntz = 0;

	
	/* ======================================================================================== */
	// x
	start = -1;
	for (i=si; i<si+m; i++) {
		if (i%2==0 && i==si && i!=0) { continue; } /* reject first ghost if its's even */
		if (i%2==0) {
			start = i;
			break;
		}
	}
	if (start == -1) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot determine start index in I"); }
	while (start + 2 * cntx < si+m) {
		PetscInt n0,n2;
		
		n0 = start + 2 * cntx;
		n2 = n0 + 2;
		
		if (n2<sig+mg) {
//			PetscPrintf(PETSC_COMM_SELF,"ELEMENT(i) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,si,si+m-1 );
			cntx++;
			continue;
		}		
		
		if (si+m-n2>1) {
//			PetscPrintf(PETSC_COMM_SELF,"ELEMENT(i) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,si,si+m-1 );
			cntx++;
			continue;
		}		
		
		if (si+m-n2<=1) {
			break;
		}		
	}
	
	/* ======================================================================================== */
	// y
	start = -1;
	for (j=sj; j<sj+n; j++) {
		if (j%2==0 && j==sj && j!=0) { continue; } /* reject first ghost if its's even */
		if (j%2==0) {
			start = j;
			break;
		}
	}
	if (start == -1) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot determine start index in J"); }	
	while (start + 2 * cnty < sj+n) {
		PetscInt n0,n2;
		
		n0 = start + 2 * cnty;
		n2 = n0 + 2;
		
		/* if start and end of element are inside global range - keep it */
		if (n2<sjg+ng) {
//			PetscPrintf(PETSC_COMM_SELF,"ELEMENT(j) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,sj,sj+n-1 );
			cnty++;
			continue;
		}		
		
		if (sj+n-n2>1) {
//			PetscPrintf(PETSC_COMM_SELF,"ELEMENT(j) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,sj,sj+n-1 );
			cnty++;
			continue;
		}		
		
		if (sj+n-n2<=1) {
			break;
		}		
	}
	
	/* ======================================================================================== */
	// z
	start = -1;
	for (k=sk; k<sk+p; k++) {
		if (k%2==0 && k==sk && k!=0) { continue; } /* reject first ghost if its's even */
		if (k%2==0) {
			start = k;
			break;
		}
	}
	if (start == -1) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot determine start index in K"); }
	while (start + 2 * cntz < p+sk) {
		PetscInt n0,n2;
		
		n0 = start + 2 * cntz;
		n2 = n0 + 2;
		
		//		printf("last-n2 = %d \n", sk+p-n2 );
		/* if start and end of element are inside global range - keep it */
		if (n2<skg+pg) {
//			PetscPrintf(PETSC_COMM_SELF,"ELEMENT(k) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,sk,sk+p-1 );
			//			printf("[GLOBAL] usng start id [k] %d [%d-%d]l [%d-%d]g\n", n0,sk,sk+p-1,skg,skg+pg-1);
			cntz++;
			continue;
		}		
		
		if (sk+p-n2>1) {
//			PetscPrintf(PETSC_COMM_SELF,"ELEMENT(k) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,sk,sk+p-1 );
			//			printf("[LOCAL] usng start id [k] %d [%d-%d]l [%d-%d]g\n", n0,sk,sk+p-1,skg,skg+pg-1);
			cntz++;
			continue;
		}		
		
		if (sk+p-n2<=1) {
			/* this means the element is taking two entries from the ghost cells */
			//			printf("[OUTSIDE] element (%d .. %d) [%d-%d]l [%d-%d]g\n", n0,n2,sk,sk+p-1,skg,skg+pg-1);
			break;
		}		
	}
	
	/* ======================================================================================== */
	
#if 0	
	// z
	for (k=sk; k<sk+p; k++) {
		if (k%2==0 && k==sk && k!=0) { continue; } /* reject first ghost if its's even */

		if (k%2==0) {
			start = k;
			break;
		}
	}
	while (start + 2 * cntz < p+sk) {
		PetscInt n0,n2;
		
		n0 = start + 2 * cntz;
		n2 = n0 + 2;
		
//		printf("last-n2 = %d \n", sk+p-n2 );
		/* if start and end of element are inside global range - keep it */
		if (n2<skg+pg) {
			PetscPrintf(PETSC_COMM_SELF,"ELEMENT (%4d - %4d) inside range l[%4d-%4d]\n", n0,n2,sk,sk+p-1 );
//			printf("[GLOBAL] usng start id [k] %d [%d-%d]l [%d-%d]g\n", n0,sk,sk+p-1,skg,skg+pg-1);
			cntz++;
			continue;
		}		
		
		if (sk+p-n2>1) {
			PetscPrintf(PETSC_COMM_SELF,"ELEMENT (%4d - %4d) inside range l[%4d-%4d]\n", n0,n2,sk,sk+p-1 );
//			printf("[LOCAL] usng start id [k] %d [%d-%d]l [%d-%d]g\n", n0,sk,sk+p-1,skg,skg+pg-1);
			cntz++;
			continue;
		}		

		if (sk+p-n2<=1) {
			/* this means the element is taking two entries from the ghost cells */
//			printf("[OUTSIDE] element (%d .. %d) [%d-%d]l [%d-%d]g\n", n0,n2,sk,sk+p-1,skg,skg+pg-1);
			break;
		}		
		
	}
	/*printf("mx,my,mz = %d %d %d \n", cntx,cnty,cntz);*/
#endif
	
	if (mx) { *mx = cntx; }
	if (my) { *my = cnty; }
	if (mz) { *mz = cntz; }
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetCornersElementQ2"
PetscErrorCode DMDAGetCornersElementQ2(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	const PetscInt order = 2;
	PetscInt i,j,k;
	PetscInt cntx,cnty,cntz;
	PetscInt si,sj,sk,m,n,p,M,N,P,width;
	PetscInt sig,sjg,skg,mg,ng,pg;
	PetscErrorCode ierr;
	PetscMPIInt rank;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sig,&sjg,&skg,&mg,&ng,&pg);CHKERRQ(ierr);
	if (width!=2) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 2 for Q2");
	}
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	/*PetscPrintf(PETSC_COMM_SELF,"[%d]: %d->%d : %d->%d \n", rank,si,si+m,sj,sj+n);*/
	
	cntx = cnty = cntz = 0;
	// x
	for (i=si; i<si+m; i++) {
		if (i%2==0 && i==si && i!=0) { continue; } /* reject first ghost if its's even */
		if (i%2==0) {
			*sei = i;
			break;
		}
	}

	// y
	for (j=sj; j<sj+n; j++) {
		if (j%2==0 && j==sj && j!=0) { continue; } /* reject first ghost if its's even */
		if (j%2==0) {
			*sej = j;
			break;
		}
	}
	
	// z
	for (k=sk; k<sk+p; k++) {
		if (k%2==0 && k==sk && k!=0) { continue; } /* reject first ghost if its's even */
		if (k%2==0) {
			*sek = k;
			break;
		}
	}
	/*PetscPrintf(PETSC_COMM_SELF,"si,sj,sk = %d %d %d \n", *sei,*sej,*sek);*/
	
	ierr = DMDAGetLocalSizeElementQ2(da,mx,my,mz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetOwnershipRangesElementQ2"
PetscErrorCode DMDAGetOwnershipRangesElementQ2(DM da,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt **si,PetscInt **sj,PetscInt **sk,PetscInt **_mx,PetscInt **_my,PetscInt **_mz)
{
	PetscMPIInt nproc,rank;
	MPI_Comm comm;
	PetscInt M,N,P,pM,pN,pP;
	PetscInt i,j,k,II,dim,esi,esj,esk,mx,my,mz;
	PetscInt *olx,*oly,*olz;
	PetscInt *lmx,*lmy,*lmz,*tmp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* create file name */
	PetscObjectGetComm( (PetscObject)da, &comm );
	ierr = MPI_Comm_size( comm, &nproc );
	ierr = MPI_Comm_rank( comm, &rank );
	
	ierr = DMDAGetInfo( da, &dim, &M,&N,&P, &pM,&pN,&pP, 0, 0, 0,0,0, 0 );CHKERRQ(ierr);
	ierr = DMDAGetCornersElementQ2(da,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
	if (dim == 1) {
		pN = 1;
		pP = 1;
	}
	if (dim == 2) {
		pP = 1;
	}
	
	ierr = PetscMalloc( sizeof(PetscInt)*(nproc), &tmp );CHKERRQ(ierr);
	
	ierr = PetscMalloc( sizeof(PetscInt)*(pM+1), &olx );CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*(pN+1), &oly );CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*(pP+1), &olz );CHKERRQ(ierr);
	
	ierr = PetscMalloc( sizeof(PetscInt)*(pM+1), &lmx );CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*(pN+1), &lmy );CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*(pP+1), &lmz );CHKERRQ(ierr);
	
	if (dim >= 1) {
		ierr = MPI_Allgather ( &esi, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		j = k = 0;
		for( i=0; i<pM; i++ ) {
			PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
			olx[i] = tmp[procid];
		}
		
		ierr = MPI_Allgather ( &mx, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		j = k = 0;
		for( i=0; i<pM; i++ ) {
			PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
			lmx[i] = tmp[procid];
		}
	}
	
	if (dim >= 2 ) {
		ierr = MPI_Allgather ( &esj, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		i = k = 0;
		for( j=0; j<pN; j++ ) {
			PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
			oly[j] = tmp[procid];
		}
		
		ierr = MPI_Allgather ( &my, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		i = k = 0;
		for( j=0; j<pN; j++ ) {
			PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
			lmy[j] = tmp[procid];
		}
	}
	
	if (dim == 3 ) {
		ierr = MPI_Allgather ( &esk, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		i = j = 0;
		for( k=0; k<pP; k++ ) {
			PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
			olz[k] = tmp[procid];
		}
		
		ierr = MPI_Allgather ( &mz, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		i = j = 0;
		for( k=0; k<pP; k++ ) {
			PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
			lmz[k] = tmp[procid];
		}
	}
	
	if(m) { *m = pM; }
	if(n) { *n = pN; }
	if(p) { *p = pP; }
	
	if(si) { *si = olx; } else { ierr = PetscFree(olx);CHKERRQ(ierr); }
	if(sj) { *sj = oly; } else { ierr = PetscFree(oly);CHKERRQ(ierr); }
	if(sk) { *sk = olz; } else { ierr = PetscFree(olz);CHKERRQ(ierr); }
	
	if(_mx) { *_mx = lmx; } else { ierr = PetscFree(lmx);CHKERRQ(ierr); }
	if(_my) { *_my = lmy; } else { ierr = PetscFree(lmy);CHKERRQ(ierr); }
	if(_mz) { *_mz = lmz; } else { ierr = PetscFree(lmz);CHKERRQ(ierr); }
	
	ierr = PetscFree(tmp);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAGetElements_DA_Q2_3D"
PetscErrorCode DMDAGetElements_DA_Q2_3D(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx)
{
  DM_DA          *da = (DM_DA*)dm->data;
	const PetscInt order = 2;
	PetscErrorCode ierr;
	PetscInt *idx,mx,my,mz,_npe, M,N,P;
	PetscInt ei,ej,ek,i,j,k,elcnt,esi,esj,esk,gsi,gsj,gsk,nid[27],n,d,X,Y,Z,width;
	PetscInt *el;
	PetscInt dof;
	PetscMPIInt rank;
	PetscFunctionBegin;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width!=2) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 2 for Q2");
	}
	
	_npe = (order + 1)*(order + 1)*(order + 1);
  if (!da->e) {
		ierr = DMDAGetCornersElementQ2(dm,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
		ierr = PetscMalloc(sizeof(PetscInt)*(mx*my*mz*_npe+1),&idx);CHKERRQ(ierr);
		ierr = DMDAGetGhostCorners(dm,&gsi,&gsj,&gsk, &X,&Y,&Z);CHKERRQ(ierr);
		ierr = DMDAGetInfo(dm,0, 0,0,0, 0,0,0, &dof,0, 0,0,0, 0);CHKERRQ(ierr);
		
		elcnt = 0;

		for (ek=0; ek<mz; ek++) {
			k = esk-gsk + 2*ek;
			for (ej=0; ej<my; ej++) {
				j = esj-gsj + 2*ej;
				for (ei=0; ei<mx; ei++) {
					i = esi-gsi + 2*ei;

					el = &idx[_npe*elcnt];
					
					//				printf("[%d]: i=%d \n", rank,i );
					nid[ 0] = (i  ) + (j  ) *X  + (k  ) *X*Y; 
					nid[ 1] = (i+1) + (j  ) *X  + (k  ) *X*Y; 
					nid[ 2] = (i+2) + (j  ) *X  + (k  ) *X*Y; 
					
					nid[ 3] = (i  ) + (j+1) *X  + (k  ) *X*Y; 
					nid[ 4] = (i+1) + (j+1) *X  + (k  ) *X*Y; 
					nid[ 5] = (i+2) + (j+1) *X  + (k  ) *X*Y; 
					
					nid[ 6] = (i  ) + (j+2) *X  + (k  ) *X*Y; 
					nid[ 7] = (i+1) + (j+2) *X  + (k  ) *X*Y; 
					nid[ 8] = (i+2) + (j+2) *X  + (k  ) *X*Y; 
					//
					nid[ 9] = (i  ) + (j  ) *X  + (k+1) *X*Y; 
					nid[10] = (i+1) + (j  ) *X  + (k+1) *X*Y; 
					nid[11] = (i+2) + (j  ) *X  + (k+1) *X*Y; 
					
					nid[12] = (i  ) + (j+1) *X  + (k+1) *X*Y; 
					nid[13] = (i+1) + (j+1) *X  + (k+1) *X*Y; 
					nid[14] = (i+2) + (j+1) *X  + (k+1) *X*Y; 
					
					nid[15] = (i  ) + (j+2) *X  + (k+1) *X*Y; 
					nid[16] = (i+1) + (j+2) *X  + (k+1) *X*Y; 
					nid[17] = (i+2) + (j+2) *X  + (k+1) *X*Y; 
					//
					nid[18] = (i  ) + (j  ) *X  + (k+2) *X*Y; 
					nid[19] = (i+1) + (j  ) *X  + (k+2) *X*Y; 
					nid[20] = (i+2) + (j  ) *X  + (k+2) *X*Y; 
					
					nid[21] = (i  ) + (j+1) *X  + (k+2) *X*Y; 
					nid[22] = (i+1) + (j+1) *X  + (k+2) *X*Y; 
					nid[23] = (i+2) + (j+1) *X  + (k+2) *X*Y; 
					
					nid[24] = (i  ) + (j+2) *X  + (k+2) *X*Y; 
					nid[25] = (i+1) + (j+2) *X  + (k+2) *X*Y; 
					nid[26] = (i+2) + (j+2) *X  + (k+2) *X*Y; 

/*					
					if(rank==1) printf("%d,%d,%d %d %d %d, %d %d %d , %d %d %d \n", i,j,k,
											 nid[0],nid[1],nid[2],nid[3],nid[4],nid[5],nid[6],nid[7],nid[8] );
*/					
					for (n=0; n<_npe; n++) {
						if (nid[n]>M*N*P) { 
							SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of global nodes");
						}
						if (nid[n]>X*Y*Z) { 
							SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of local nodes");
						}
						el[n] = nid[n]; //gidx[dof*nid[n]+0]/dof;
					}
					
					elcnt++;
				}
			}
		}
		
		da->e  = idx;
		da->ne = elcnt;
	}
	
	*eidx = da->e;
	*npe = _npe;
	*nel = da->ne;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAGetElements_DA_Q2"
PetscErrorCode DMDAGetElements_DA_Q2(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[])
{
  DM_DA          *da = (DM_DA*)dm->data;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (da->dim==-1) {
    *nel = 0; *nen = 0; *e = NULL;
  } else if (da->dim==1) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DMDA doesn't support Q2 in 1D");
  } else if (da->dim==2) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DMDA doesn't support Q2 in 2D");
  } else if (da->dim==3) {
    ierr = DMDAGetElements_DA_Q2_3D(dm,nel,nen,e);CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"DMDA dimension not 1, 2, or 3, it is %D\n",da->dim);
  }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAGetElements_DA_P0MD_3D"
PetscErrorCode DMDAGetElements_DA_P0MD_3D(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx)
{
  DM_DA          *da = (DM_DA*)dm->data;
	PetscErrorCode ierr;
	PetscInt *idx,mx,my,mz,_npe;
	PetscInt ei,ej,ek,i,j,k,elcnt,cnt,esi,esj,esk,gsi,gsj,gsk,nid[100],d,X,width;
	PetscInt *el,M,N,P,dof;
	PetscMPIInt rank;
	PetscFunctionBegin;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width!=0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 0 for P0 with multi-dofs");
	}
	
	_npe = 1;
	ierr = DMDAGetInfo(dm,0, 0,0,0, 0,0,0, &dof,0, 0,0,0, 0);CHKERRQ(ierr);
  if (dof>100) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"dof > 100"); }
	
	if (!da->e) {
		ierr = DMDAGetCorners(dm,&esi,&esj,0,&mx,&my,&mz);CHKERRQ(ierr);
		ierr = PetscMalloc(sizeof(PetscInt)*(mx*my*mz*_npe*dof+1),&idx);CHKERRQ(ierr); /* we add one extra space to allow for cases when a proc has zero elements, and we don't want to malloc 0 bytes */
		
		elcnt = 0;
		cnt = 0;
		for (ek=0; ek<mz; ek++) {
			for (ej=0; ej<my; ej++) {
				for (ei=0; ei<mx; ei++) {
					el = &idx[_npe*dof*elcnt];
					
					for (d=0; d<_npe*dof; d++) {
						nid[d] = (_npe*dof) * elcnt + d;
					}

					for (d=0; d<_npe*dof; d++) {
						if (nid[d]>M*N*P*dof) {
							SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of global dofs");
						}
						el[d] = nid[d];
					}
					
					elcnt++;
				}
			}
		}
		
		da->e  = idx;
		da->ne = elcnt;
	}
	
	*eidx = da->e;
	*npe  = _npe * dof;
	*nel  = da->ne;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAGetElements_DA_P1"
PetscErrorCode DMDAGetElements_DA_P1(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[])
{
  DM_DA          *da = (DM_DA*)dm->data;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (da->dim==-1) {
    *nel = 0; *nen = 0; *e = NULL;
  } else if (da->dim==1) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DMDA doesn't support P1 in 1D");
  } else if (da->dim==2) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DMDA doesn't support P1 in 2D");
  } else if (da->dim==3) {
    ierr = DMDAGetElements_DA_P0MD_3D(dm,nel,nen,e);CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"DMDA dimension not 1, 2, or 3, it is %D\n",da->dim);
  }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAGetElements_pTatinQ2P1"
PetscErrorCode DMDAGetElements_pTatinQ2P1(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[])
{
  PetscErrorCode ierr;
  PetscInt dof,sw;
  static PetscInt beenHereQ2 = 0;
  static PetscInt beenHereP0 = 0;
  
	ierr = DMDAGetInfo(dm, 0, 0,0,0, 0,0,0, &dof,&sw, 0,0,0, 0);CHKERRQ(ierr);
  
  if (sw == 2) {
    //if (!beenHereQ2) { PetscPrintf(PETSC_COMM_WORLD,"ASSUMING stencil width=2 implies Q2 basis...\n"); }
    ierr = DMDAGetElements_DA_Q2(dm,nel,nen,e);CHKERRQ(ierr);
    beenHereQ2 = 1;
  } else if (sw == 0) {
    //if (!beenHereP0) { PetscPrintf(PETSC_COMM_WORLD,"ASSUMING stencil width=0 implies P0 basis...\n"); }
    ierr = DMDAGetElements_DA_P1(dm,nel,nen,e);CHKERRQ(ierr);
    beenHereP0 = 1;
  } else {
    ierr = DMDAGetElements(dm,nel,nen,e);CHKERRQ(ierr);
  }
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetElementType_Q2"
PetscErrorCode  DMDASetElementType_Q2(DM da)
{
  DM_DA          *dd = (DM_DA*)da->data;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscValidHeaderSpecific(da,DM_CLASSID,1);
	//  PetscValidLogicalCollectiveEnum(da,etype,2);
  if (dd->elementtype) {
    ierr = PetscFree(dd->e);CHKERRQ(ierr);
		//    dd->elementtype = etype;
    dd->ne          = 0; 
    dd->e           = NULL;
  }
  //da->ops->getelements = DMGetElements_DA_Q2;
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetElementType_P1"
PetscErrorCode  DMDASetElementType_P1(DM da)
{
  DM_DA          *dd = (DM_DA*)da->data;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscValidHeaderSpecific(da,DM_CLASSID,1);
	//  PetscValidLogicalCollectiveEnum(da,etype,2);
  if (dd->elementtype) {
    ierr = PetscFree(dd->e);CHKERRQ(ierr);
		//    dd->elementtype = etype;
    dd->ne          = 0; 
    dd->e           = NULL;
  }
	//da->ops->getelements = DMGetElements_DA_P1;
	
  PetscFunctionReturn(0);
}


/* element helpers */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetElementCoordinatesQ2_3D"
PetscErrorCode DMDAGetElementCoordinatesQ2_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		elcoords[3*n  ] = LA_gcoords[3*elnid[n]  ];
		elcoords[3*n+1] = LA_gcoords[3*elnid[n]+1];
		elcoords[3*n+2] = LA_gcoords[3*elnid[n]+2];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAGetScalarElementFieldQ2_3D"
PetscErrorCode DMDAGetScalarElementFieldQ2_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		elfield[n] = LA_gfield[elnid[n]];
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DMDAGetVectorElementFieldQ2_3D"
PetscErrorCode DMDAGetVectorElementFieldQ2_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		elfield[3*n  ] = LA_gfield[3*elnid[n]  ];
		elfield[3*n+1] = LA_gfield[3*elnid[n]+1];
		elfield[3*n+2] = LA_gfield[3*elnid[n]+2];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAGetScalarElementField"
PetscErrorCode DMDAGetScalarElementField(PetscScalar elfield[],PetscInt npe,PetscInt elnid[],PetscScalar LA_gfield[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<npe; n++) {
		elfield[n] = LA_gfield[elnid[n]];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetValuesLocalStencil_AddValues_DOF"
PetscErrorCode DMDASetValuesLocalStencil_AddValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[])
{
  PetscInt n,d,el_idx,idx;
	
  PetscFunctionBegin;
	for (d=0; d<ndof; d++) {
		for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
			el_idx = ndof*n + d;
			idx    = eqn[el_idx];
			fields_F[idx] += Fe[el_idx];
		}
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetValuesLocalStencil_SetValues_DOF"
PetscErrorCode DMDASetValuesLocalStencil_SetValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[])
{
  PetscInt n,d,el_idx,idx;
	
  PetscFunctionBegin;
	for (d=0; d<ndof; d++) {
		for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
			el_idx = ndof*n + d;
			idx    = eqn[el_idx];
			fields_F[idx] = Fe[el_idx];
		}
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Q2GetElementLocalIndicesDOF"
PetscErrorCode Q2GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[])
{
	PetscInt n,d;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (d=0; d<ndof; d++) {
		for (n=0; n<U_BASIS_FUNCTIONS; n++) {
			el_localIndices[ndof*n+d] = ndof*elnid[n]+d;
		}
	}		
	PetscFunctionReturn(0);
}


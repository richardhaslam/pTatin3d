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
 **    filename:   dmda_element_q1macrop1.c
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
#include <petscdm.h>
#include <petsc/private/dmdaimpl.h>
#include "ptatin3d_defs.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1macrop1.h"

typedef struct _p_DMDAEQ1MacroCtx *DMDAEQ1MacroCtx;

struct _p_DMDAEQ1MacroCtx {
	PetscInt mx_mixed_space,  my_mixed_space,  mz_mixed_space;
	PetscInt mx_natural_space,my_natural_space,mz_natural_space;
	PetscInt Mx_mixed_space,  My_mixed_space,  Mz_mixed_space;
	PetscInt Mx_natural_space,My_natural_space,Mz_natural_space;
	PetscInt *element_macro_index;            /* (mx.my.mz)_natural_space : maps to macro index */
	PetscInt *element_node_map_natural_space; /* (mx.my.mz)_natural_space * 8 */
	PetscInt *element_node_map_mixed_space;   /* (mx.my.mz)_mixed_space * 27 */
	PetscInt *element_macro_local_index;      /* (mx.my.mz)_natural_space : maps to local index in macro cell (for quadrature) */
	
//	fp_GetElementsNatural();
//	fp_GetElementsMixed();
//	fp_GetBasisFunction();
//	fp_GetBasisFunctionDerivatives();
};

#define DMDAEQ1MacroCtxName "DMDAEQ1Macro"

PetscErrorCode _DMDAEQ1Macro_MixedSpace_GetSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ);
PetscErrorCode _DMDAEQ1Macro_NaturalSpace_GetSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ);
PetscErrorCode _DMDAEQ1Macro_MixedSpace_GetLocalSizeElement(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz);
PetscErrorCode _DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz);
PetscErrorCode _DMDAEQ1Macro_MixedSpace_GetOwnershipRangesElement(DM da,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt **si,PetscInt **sj,PetscInt **sk,PetscInt **_mx,PetscInt **_my,PetscInt **_mz);
PetscErrorCode _DMDAEQ1Macro_MixedSpace_GetElements3D(DM dm,PetscInt *nel,PetscInt *npe,PetscInt **eidx);
PetscErrorCode _DMDAEQ1Macro_NaturalSpace_GetElements3D(DM dm,PetscInt *nel,PetscInt *npe,PetscInt **eidx);
PetscErrorCode _DMDAEQ1Macro_NaturalSpaceToMixedSpace3D(DM dm,PetscInt **natural_mixed_map);
PetscErrorCode _DMDAEQ1Macro_NaturalSpaceToMixedLocalSpace3D(DM dm,PetscInt **natural_mixed_map);

PetscErrorCode DMDAEQ1Macro_MixedSpace_GetCornersElement(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz);


#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_FetchContext"
PetscErrorCode DMDAEQ1Macro_FetchContext(DM da,DMDAEQ1MacroCtx *ctx)
{
	PetscContainer container;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	container = NULL;
	ierr = PetscObjectQuery((PetscObject)da,DMDAEQ1MacroCtxName,(PetscObject*)&container);CHKERRQ(ierr);
	if (!container) SETERRQ1(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONG,"No data with name \"%s\" was composed with this DAE",DMDAEQ1MacroCtxName);
	ierr = PetscContainerGetPointer(container,(void**)ctx);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}


/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "_DMDAEQ1Macro_MixedSpace_GetSizeElement"
PetscErrorCode _DMDAEQ1Macro_MixedSpace_GetSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ)
{
	const PetscInt order = 2;
	PetscInt M,N,P,width;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (MX) { *MX = (M-1)/order; }
	if (MY) { *MY = (N-1)/order; }
	if (MZ) { *MZ = (P-1)/order; }

	if (width!=1) {
		SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Stencil width must be 1 for Q1Macro");
	}
	
	/* M = order*mx+1 */
	if ( (M-1)%order != 0 ) {
		SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"DMDA is not compatible with Q1Macro elements in x direction");
	}
	
	//
	if ( (N-1)%order != 0 ) {
		SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"DMDA is not compatible with Q1Macro elements in y direction");
	}
	
	//
	if ( (P-1)%order != 0 ) {
		SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"DMDA is not compatible with Q1Macro elements in z direction");
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DMDAEQ1Macro_NaturalSpace_GetSizeElement"
PetscErrorCode _DMDAEQ1Macro_NaturalSpace_GetSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ)
{
	PetscInt m,n,p;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
  
  m = 0; n = 0; p = 0;
	ierr = _DMDAEQ1Macro_MixedSpace_GetSizeElement(da,&m,&n,&p);CHKERRQ(ierr);
	if (MX) { *MX = 2 * m; }
	if (MY) { *MY = 2 * n; }
	if (MZ) { *MZ = 2 * p; }
	
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "_DMDAEQ1Macro_MixedSpace_GetLocalSizeElement"
PetscErrorCode _DMDAEQ1Macro_MixedSpace_GetLocalSizeElement(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
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
	ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)da),&rank);CHKERRQ(ierr);
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
	
	
	if (mx) { *mx = cntx; }
	if (my) { *my = cnty; }
	if (mz) { *mz = cntz; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement"
PetscErrorCode _DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	PetscInt m,n,p;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = _DMDAEQ1Macro_MixedSpace_GetLocalSizeElement(da,&m,&n,&p);CHKERRQ(ierr);
	if (mx) { *mx = 2 * m; }
	if (my) { *my = 2 * n; }
	if (mz) { *mz = 2 * p; }
	
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_GetCornersElement"
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetCornersElement(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	PetscInt i,j,k;
	PetscInt si,sj,sk,m,n,p,M,N,P,width;
	PetscInt sig,sjg,skg,mg,ng,pg;
	PetscErrorCode ierr;
	int rank;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sig,&sjg,&skg,&mg,&ng,&pg);CHKERRQ(ierr);
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	/*PetscPrintf(PETSC_COMM_SELF,"[%d]: %d->%d : %d->%d \n", rank,si,si+m,sj,sj+n);*/

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
	
	ierr = _DMDAEQ1Macro_MixedSpace_GetLocalSizeElement(da,mx,my,mz);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "_DMDAEQ1Macro_MixedSpace_GetOwnershipRangesElement"
PetscErrorCode _DMDAEQ1Macro_MixedSpace_GetOwnershipRangesElement(DM da,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt **si,PetscInt **sj,PetscInt **sk,PetscInt **_mx,PetscInt **_my,PetscInt **_mz)
{
	PetscMPIInt nproc,rank;
	MPI_Comm comm;
	PetscInt M,N,P,pM,pN,pP;
	PetscInt i,j,k,dim,esi,esj,esk,mx,my,mz;
	PetscInt *olx,*oly,*olz;
	PetscInt *lmx,*lmy,*lmz,*tmp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* create file name */
	PetscObjectGetComm( (PetscObject)da, &comm );
	ierr = MPI_Comm_size( comm, &nproc );CHKERRQ(ierr);
	ierr = MPI_Comm_rank( comm, &rank );CHKERRQ(ierr);
	
	ierr = DMDAGetInfo( da, &dim, &M,&N,&P, &pM,&pN,&pP, 0, 0, 0,0,0, 0 );CHKERRQ(ierr);
	ierr = DMDAEQ1Macro_MixedSpace_GetCornersElement(da,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
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
#define __FUNCT__ "_DMDAEQ1Macro_MixedSpace_GetElements3D"
PetscErrorCode _DMDAEQ1Macro_MixedSpace_GetElements3D(DM dm,PetscInt *nel,PetscInt *npe,PetscInt **eidx)
{
	PetscErrorCode ierr;
	PetscInt *idx,mx,my,mz,_npe, M,N,P;
	PetscInt ei,ej,ek,i,j,k,elcnt,esi,esj,esk,gsi,gsj,gsk,nid[27],n,X,Y,Z,width;
	PetscInt *el;
	PetscInt dof;
	int rank;
	PetscFunctionBegin;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width!=1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 1 for Q1Macro");
	}
	
	_npe = 3 * 3 * 3;
	ierr = DMDAEQ1Macro_MixedSpace_GetCornersElement(dm,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
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

				for (n=0; n<_npe; n++) {
					if (nid[n]>M*N*P) { 
						SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of global nodes");
					}
					if (nid[n]>X*Y*Z) { 
						SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of local nodes");
					}
					el[n] = nid[n];
				}
				
				elcnt++;
			}
		}
	}
	
	*eidx = idx;
	*npe  = _npe;
	*nel  = elcnt;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DMDAEQ1Macro_NaturalSpace_GetElements3D"
PetscErrorCode _DMDAEQ1Macro_NaturalSpace_GetElements3D(DM dm,PetscInt *nel,PetscInt *npe,PetscInt **eidx)
{
	PetscErrorCode ierr;
	PetscInt *idx,mx,my,mz,_npe, M,N,P;
	PetscInt ei,ej,ek,i,j,k,ii,jj,kk,elcnt,esi,esj,esk,gsi,gsj,gsk,nid[8],n,X,Y,Z,width;
	PetscInt mx_natural,my_natural,mz_natural;
	PetscInt *el;
	PetscInt dof;
	int rank;
	PetscFunctionBegin;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width!=1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 1 for Q1Macro");
	}
	
	_npe = 2 * 2 * 2;

	ierr = DMDAEQ1Macro_MixedSpace_GetCornersElement(dm,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
  mx_natural = 0; my_natural = 0; mz_natural = 0;
	ierr = _DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(dm,&mx_natural,&my_natural,&mz_natural);CHKERRQ(ierr);
	//ierr = DMDAEQ1Macro_NaturalSpace_GetCornersElement(dm,&esi_natural,&esj_natural,&esk_natural,);CHKERRQ(ierr);

	ierr = PetscMalloc(sizeof(PetscInt)*(mx_natural*my_natural*mz_natural*_npe+1),&idx);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(dm,&gsi,&gsj,&gsk, &X,&Y,&Z);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dm,0, 0,0,0, 0,0,0, &dof,0, 0,0,0, 0);CHKERRQ(ierr);
	
	elcnt = 0;
	
	/* loop over macro elements */
	for (ek=0; ek<mz; ek++) {
		k = esk-gsk + 2*ek;
		for (ej=0; ej<my; ej++) {
			j = esj-gsj + 2*ej;
			for (ei=0; ei<mx; ei++) {
				i = esi-gsi + 2*ei;
				
				/* loop over sub elements */
				for (kk=0; kk<2; kk++) {
					for (jj=0; jj<2; jj++) {
						for (ii=0; ii<2; ii++) {
							PetscInt natural_idx;
							PetscInt ni,nj,nk;
							
							//macro_idx   = (ei)      + (ej)*mx              + (ek)*mx*my;
							natural_idx = (2*ei+ii) + (2*ej+jj)*mx_natural + (2*ek+kk)*mx_natural*my_natural;
							
							el = &idx[_npe*natural_idx];
							
							ni = i + ii;
							nj = j + jj;
							nk = k + kk;
							
							//				printf("[%d]: i=%d \n", rank,i );
							nid[0] = (ni  ) + (nj  ) *X  + (nk  ) *X*Y; 
							nid[1] = (ni+1) + (nj  ) *X  + (nk  ) *X*Y; 
							
							nid[2] = (ni  ) + (nj+1) *X  + (nk  ) *X*Y; 
							nid[3] = (ni+1) + (nj+1) *X  + (nk  ) *X*Y; 
							//
							nid[4] = (ni  ) + (nj  ) *X  + (nk+1) *X*Y; 
							nid[5] = (ni+1) + (nj  ) *X  + (nk+1) *X*Y; 
							
							nid[6] = (ni  ) + (nj+1) *X  + (nk+1) *X*Y; 
							nid[7] = (ni+1) + (nj+1) *X  + (nk+1) *X*Y; 
							
							//printf("[%d]: macro[%d,%d,%d]:sub[%d,%d,%d]:ijk[%d,%d,%d] \n", rank,ei,ej,ek,ii,jj,kk,i,j,k );
							//printf("  nid: {%d , %d , %d , %d ; %d , %d , %d , %d} \n", nid[0],nid[1],nid[2],nid[3],nid[4],nid[5],nid[6],nid[7] );
							for (n=0; n<_npe; n++) {
								if (nid[n]>M*N*P) { 
									SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of global nodes");
								}
								if (nid[n]>X*Y*Z) { 
									SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of local nodes");
								}
								el[n] = nid[n];
							}
							
							elcnt++;
							
						}
					}
				}
				
				
			}
		}
	}

	*eidx = idx;
	*npe  = _npe;
	*nel  = elcnt;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DMDAEQ1Macro_NaturalSpaceToMixedSpace3D"
PetscErrorCode _DMDAEQ1Macro_NaturalSpaceToMixedSpace3D(DM dm,PetscInt **natural_mixed_map)
{
	PetscInt *idx,mx,my,mz, M,N,P;
	PetscInt esi,esj,esk,ei,ej,ek,ii,jj,kk,width;
	PetscInt mx_natural,my_natural,mz_natural;
	PetscErrorCode ierr;

	PetscFunctionBegin;

	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width!=1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 1 for Q1Macro");
	}
	
	ierr = DMDAEQ1Macro_MixedSpace_GetCornersElement(dm,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);

  mx_natural = 0; my_natural = 0; mz_natural = 0;
	ierr = _DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(dm,&mx_natural,&my_natural,&mz_natural);CHKERRQ(ierr);
	//ierr = DMDAEQ1Macro_NaturalSpace_GetCornersElement(dm,&esi_natural,&esj_natural,&esk_natural,);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(mx_natural*my_natural*mz_natural),&idx);CHKERRQ(ierr);
	
	/* loop over macro elements */
	for (ek=0; ek<mz; ek++) {
		/* loop over sub elements */
		for (kk=0; kk<2; kk++) {
			
			for (ej=0; ej<my; ej++) {
				for (jj=0; jj<2; jj++) {
					
					for (ei=0; ei<mx; ei++) {
						for (ii=0; ii<2; ii++) {
							PetscInt macro_idx,natural_idx;
							
							macro_idx   = (ei)      + (ej)*mx              + (ek)*mx*my;
							natural_idx = (2*ei+ii) + (2*ej+jj)*mx_natural + (2*ek+kk)*mx_natural*my_natural;
							
							idx[natural_idx] = macro_idx;
							
						}
					}
					
				}		
			}
			
		}
	}
	
	*natural_mixed_map = idx;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DMDAEQ1Macro_NaturalSpaceToMixedLocalSpace3D"
PetscErrorCode _DMDAEQ1Macro_NaturalSpaceToMixedLocalSpace3D(DM dm,PetscInt **natural_mixed_map)
{
	PetscInt *idx,mx,my,mz, M,N,P;
	PetscInt esi,esj,esk,ei,ej,ek,ii,jj,kk,width;
	PetscInt mx_natural,my_natural,mz_natural;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width!=1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stencil width must be 1 for Q1Macro");
	}
	
	ierr = DMDAEQ1Macro_MixedSpace_GetCornersElement(dm,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
  mx_natural = 0; my_natural = 0; mz_natural = 0;
	ierr = _DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(dm,&mx_natural,&my_natural,&mz_natural);CHKERRQ(ierr);
	//ierr = DMDAEQ1Macro_NaturalSpace_GetCornersElement(dm,&esi_natural,&esj_natural,&esk_natural,);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(mx_natural*my_natural*mz_natural),&idx);CHKERRQ(ierr);
	
	/* loop over macro elements */
	for (ek=0; ek<mz; ek++) {
		/* loop over sub elements */
		for (kk=0; kk<2; kk++) {
			
			for (ej=0; ej<my; ej++) {
				for (jj=0; jj<2; jj++) {
					
					for (ei=0; ei<mx; ei++) {
						for (ii=0; ii<2; ii++) {
							PetscInt natural_idx;
							
							//macro_idx   = (ei)      + (ej)*mx              + (ek)*mx*my;
							natural_idx = (2*ei+ii) + (2*ej+jj)*mx_natural + (2*ek+kk)*mx_natural*my_natural;
							
							idx[natural_idx] = ii + jj*2 + kk*2*2;
							
						}
					}
					
				}		
			}
			
		}
	}
	
	*natural_mixed_map = idx;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetElementType_Q1Macro"
PetscErrorCode DMDASetElementType_Q1Macro(DM da)
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
  //da->ops->getelements = DMGetElements_DA_Q1;
	
  PetscFunctionReturn(0);
}

/* getters */
#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_GetSizeElement"
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ)
{
	DMDAEQ1MacroCtx ctx;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAEQ1Macro_FetchContext(da,&ctx);CHKERRQ(ierr);
	if(MX) { *MX = ctx->Mx_mixed_space; }
	if(MY) { *MY = ctx->My_mixed_space; }
	if(MZ) { *MZ = ctx->Mz_mixed_space; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_GetLocalSizeElement"
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetLocalSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ)
{
	DMDAEQ1MacroCtx ctx;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAEQ1Macro_FetchContext(da,&ctx);CHKERRQ(ierr);
	if(MX) { *MX = ctx->mx_mixed_space; }
	if(MY) { *MY = ctx->my_mixed_space; }
	if(MZ) { *MZ = ctx->mz_mixed_space; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_NaturalSpace_GetSizeElement"
PetscErrorCode DMDAEQ1Macro_NaturalSpace_GetSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ)
{
	DMDAEQ1MacroCtx ctx;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAEQ1Macro_FetchContext(da,&ctx);CHKERRQ(ierr);
	if(MX) { *MX = ctx->Mx_natural_space; }
	if(MY) { *MY = ctx->My_natural_space; }
	if(MZ) { *MZ = ctx->Mz_natural_space; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement"
PetscErrorCode DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ)
{
	DMDAEQ1MacroCtx ctx;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAEQ1Macro_FetchContext(da,&ctx);CHKERRQ(ierr);
	if(MX) { *MX = ctx->mx_natural_space; }
	if(MY) { *MY = ctx->my_natural_space; }
	if(MZ) { *MZ = ctx->mz_natural_space; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEGetElements_Q1MacroMixedSpace"
PetscErrorCode DMDAEGetElements_Q1MacroMixedSpace(DM da,PetscInt *nel,PetscInt *nen,PetscInt *e[])
{
	DMDAEQ1MacroCtx ctx;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAEQ1Macro_FetchContext(da,&ctx);CHKERRQ(ierr);
	
	if (e)   { *e = ctx->element_node_map_mixed_space; }
	if (nen) { *nen = 27; }
	if (nel) { *nel = ctx->mx_mixed_space * ctx->my_mixed_space * ctx->mz_mixed_space; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEGetElements_Q1MacroNaturalSpace"
PetscErrorCode DMDAEGetElements_Q1MacroNaturalSpace(DM da,PetscInt *nel,PetscInt *nen,PetscInt *e[])
{
	DMDAEQ1MacroCtx ctx;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAEQ1Macro_FetchContext(da,&ctx);CHKERRQ(ierr);
	
	if (e)   { *e = ctx->element_node_map_natural_space; }
	if (nen) { *nen = 8; }
	if (nel) { *nel = ctx->mx_natural_space * ctx->my_natural_space * ctx->mz_natural_space; }

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEGetElementMap_Q1MacroNaturalToMixedSpace"
PetscErrorCode DMDAEGetElementMap_Q1MacroNaturalToMixedSpace(DM da,PetscInt *nel,PetscInt *e[])
{
	DMDAEQ1MacroCtx ctx;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAEQ1Macro_FetchContext(da,&ctx);CHKERRQ(ierr);
	
	if (e)   { *e = ctx->element_macro_index; }
	if (nel) { *nel = ctx->mx_natural_space * ctx->my_natural_space * ctx->mz_natural_space; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEGetElementMap_Q1MacroNaturalToMixedLocalSpace"
PetscErrorCode DMDAEGetElementMap_Q1MacroNaturalToMixedLocalSpace(DM da,PetscInt *nel,PetscInt *e[])
{
	DMDAEQ1MacroCtx ctx;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAEQ1Macro_FetchContext(da,&ctx);CHKERRQ(ierr);
	
	if (e)   { *e = ctx->element_macro_local_index; }
	if (nel) { *nel = ctx->mx_natural_space * ctx->my_natural_space * ctx->mz_natural_space; }
	
	PetscFunctionReturn(0);
}

/* functions promoted to be public */
#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_GetOwnershipRangesElement"
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetOwnershipRangesElement(DM da,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt **si,PetscInt **sj,PetscInt **sk,PetscInt **_mx,PetscInt **_my,PetscInt **_mz)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = _DMDAEQ1Macro_MixedSpace_GetOwnershipRangesElement(da,m,n,p,si,sj,sk,_mx,_my,_mz);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


/* constructor */
#undef __FUNCT__
#define __FUNCT__ "DMDAESetType_Q1Macro"
PetscErrorCode  DMDAESetType_Q1Macro(DM da)
{
	DMDAEQ1MacroCtx ctx;
	PetscInt        nel,npe;
	PetscContainer  container;
	PetscErrorCode  ierr;
	
	/* reset all existing element information */
	ierr = DMDASetElementType_Q1Macro(da);CHKERRQ(ierr);
	

	/* create data bucket */
	ierr = PetscMalloc(sizeof(struct _p_DMDAEQ1MacroCtx),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_DMDAEQ1MacroCtx));CHKERRQ(ierr);
	
	/* fill it out */
	ierr = _DMDAEQ1Macro_MixedSpace_GetSizeElement(da,&ctx->Mx_mixed_space,&ctx->My_mixed_space,&ctx->Mz_mixed_space);CHKERRQ(ierr);
	ierr = _DMDAEQ1Macro_NaturalSpace_GetSizeElement(da,&ctx->Mx_natural_space,&ctx->My_natural_space,&ctx->Mz_natural_space);CHKERRQ(ierr);
	ierr = _DMDAEQ1Macro_MixedSpace_GetLocalSizeElement(da,&ctx->mx_mixed_space,&ctx->my_mixed_space,&ctx->mz_mixed_space);CHKERRQ(ierr);
	ierr = _DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(da,&ctx->mx_natural_space,&ctx->my_natural_space,&ctx->mz_natural_space);CHKERRQ(ierr);

	//PetscInt *element_macro_index;            /* (mx.my.mz)_natural_space : maps to macro index */
	//PetscInt *element_macro_local_index;      /* (mx.my.mz)_natural_space : maps to local index in macro cell (for quadrature) */
	//PetscInt *element_node_map_natural_space; /* (mx.my.mz)_natural_space * 8 */
	//PetscInt *element_node_map_mixed_space;   /* (mx.my.mz)_mixed_space * 27 */
	ierr = _DMDAEQ1Macro_NaturalSpaceToMixedSpace3D(da,&ctx->element_macro_index);CHKERRQ(ierr);
	ierr = _DMDAEQ1Macro_NaturalSpaceToMixedLocalSpace3D(da,&ctx->element_macro_local_index);CHKERRQ(ierr);

	ierr = _DMDAEQ1Macro_NaturalSpace_GetElements3D(da,&nel,&npe,&ctx->element_node_map_natural_space);CHKERRQ(ierr);
	ierr = _DMDAEQ1Macro_MixedSpace_GetElements3D(da,&nel,&npe,&ctx->element_node_map_mixed_space);CHKERRQ(ierr);
	
	/* attach bucket to dm */
  ierr = PetscContainerCreate(PetscObjectComm((PetscObject)da),&container);CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(container,(void*)ctx);CHKERRQ(ierr);
	
	ierr = PetscObjectCompose((PetscObject)da,DMDAEQ1MacroCtxName,(PetscObject)container);CHKERRQ(ierr);
	
	
  PetscFunctionReturn(0);
}

/* element helpers */
#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_GetElementCoordinates_3D"
PetscErrorCode DMDAEQ1Macro_GetElementCoordinates_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[])
{
	PetscInt n;

	PetscFunctionBegin;
	for (n=0; n<8; n++) {
		elcoords[3*n  ] = LA_gcoords[3*elnid[n]  ];
		elcoords[3*n+1] = LA_gcoords[3*elnid[n]+1];
		elcoords[3*n+2] = LA_gcoords[3*elnid[n]+2];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_GetElementCoordinates_3D"
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetElementCoordinates_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[])
{
	PetscInt n;

	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		elcoords[3*n  ] = LA_gcoords[3*elnid[n]  ];
		elcoords[3*n+1] = LA_gcoords[3*elnid[n]+1];
		elcoords[3*n+2] = LA_gcoords[3*elnid[n]+2];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_GetScalarElementField_3D"
PetscErrorCode DMDAEQ1Macro_GetScalarElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
{
	PetscInt n;

	PetscFunctionBegin;
	for (n=0; n<8; n++) {
		elfield[n] = LA_gfield[elnid[n]];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_GetScalarElementField_3D"
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetScalarElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
{
	PetscInt n;

	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		elfield[n] = LA_gfield[elnid[n]];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_GetVectorElementField_3D"
PetscErrorCode DMDAEQ1Macro_GetVectorElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
{
	PetscInt n;

	PetscFunctionBegin;
	for (n=0; n<8; n++) {
		elfield[3*n  ] = LA_gfield[3*elnid[n]  ];
		elfield[3*n+1] = LA_gfield[3*elnid[n]+1];
		elfield[3*n+2] = LA_gfield[3*elnid[n]+2];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_GetVectorElementField_3D"
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetVectorElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
{
	PetscInt n;

	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		elfield[3*n  ] = LA_gfield[3*elnid[n]  ];
		elfield[3*n+1] = LA_gfield[3*elnid[n]+1];
		elfield[3*n+2] = LA_gfield[3*elnid[n]+2];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_SetValuesLocalStencil_AddValues_DOF"
PetscErrorCode DMDAEQ1Macro_SetValuesLocalStencil_AddValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[])
{
  PetscInt n,d,el_idx,idx;
	
  PetscFunctionBegin;
	for (d=0; d<ndof; d++) {
		for (n = 0; n<8; n++) {
			el_idx = ndof*n + d;
			idx    = eqn[el_idx];
			fields_F[idx] += Fe[el_idx];
		}
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_SetValuesLocalStencil_AddValues_DOF"
PetscErrorCode DMDAEQ1Macro_MixedSpace_SetValuesLocalStencil_AddValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[])
{
  PetscInt n,d,el_idx,idx;
	
  PetscFunctionBegin;
	for (d=0; d<ndof; d++) {
		for (n = 0; n<27; n++) {
			el_idx = ndof*n + d;
			idx    = eqn[el_idx];
			fields_F[idx] += Fe[el_idx];
		}
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_GetElementLocalIndicesDOF"
PetscErrorCode DMDAEQ1Macro_GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[])
{
	PetscInt n,d;

	PetscFunctionBegin;
	for (d=0; d<ndof; d++) {
		for (n=0; n<8; n++) {
			el_localIndices[ndof*n+d] = ndof*elnid[n]+d;
		}
	}		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1Macro_MixedSpace_GetElementLocalIndicesDOF"
PetscErrorCode DMDAEQ1Macro_MixedSpace_GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[])
{
	PetscInt n,d;

	PetscFunctionBegin;
	for (d=0; d<ndof; d++) {
		for (n=0; n<27; n++) {
			el_localIndices[ndof*n+d] = ndof*elnid[n]+d;
		}
	}		
	PetscFunctionReturn(0);
}


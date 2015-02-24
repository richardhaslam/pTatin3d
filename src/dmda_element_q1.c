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
 **    filename:   dmda_element_q1.c
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


#include "petsc.h"
#include "petscdm.h"
#include "ptatin3d_defs.h"
#include "dmda_update_coords.h"
#include "dmda_compare.h"
#include "dmdae.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1.h"


/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetLocalSizeElementQ1"
PetscErrorCode DMDAGetLocalSizeElementQ1(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	PetscInt cntx,cnty,cntz;
	PetscInt si,sj,sk,M,N,P,m,n,p,width;
	PetscInt sig,sjg,skg,mg,ng,pg;
	PetscErrorCode ierr;
	PetscMPIInt rank;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sig,&sjg,&skg,&mg,&ng,&pg);CHKERRQ(ierr);
	if (width != 1) {
		SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Stencil width must be 1 for Q1");
	}
	ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)da),&rank);CHKERRQ(ierr);
	//printf("[%d]: i(%d->%d) : j(%d->%d) \n", rank,si,si+m,sj,sj+n);
	
	cntx = cnty = cntz = 0;
	
	/* ======================================================================================== */
	// x
	cntx = mg + 1;
	/* if on right wall */
	if ( (sig+mg) == (M) ) {
		cntx--;
	}

	// y
	cnty = ng + 1;
	if ( (sjg+ng) == (N) ) {
		cnty--;
	}

	// z
	cntz = pg + 1;
	if ( (skg+pg) == (P) ) {
		cntz--;
	}
	
	
	if (mx) { *mx = cntx; }
	if (my) { *my = cnty; }
	if (mz) { *mz = cntz; }
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetCornersElementQ1"
PetscErrorCode DMDAGetCornersElementQ1(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	PetscInt sig,sjg,skg,M,N,P,width;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sig,&sjg,&skg,0,0,0);CHKERRQ(ierr);
	if (width != 1) {
		SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Stencil width must be 1 for Q1");
	}
	
	// x
	*sei = sig;
	// y
	*sej = sjg;
	// z
	*sek = skg;
	/*PetscPrintf(PETSC_COMM_SELF,"si,sj,sk = %d %d %d \n", *sei,*sej,*sek);*/
	
	ierr = DMDAGetLocalSizeElementQ1(da,mx,my,mz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAGetElements_DA_Q1_3D"
PetscErrorCode DMDAGetElements_DA_Q1_3D(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx)
{
  DM_DA          *da = (DM_DA*)dm->data;
	const PetscInt order = 1;
	PetscErrorCode ierr;
	PetscInt *idx,mx,my,mz,_npe, M,N,P;
	PetscInt ei,ej,ek,i,j,k,elcnt,esi,esj,esk,gsi,gsj,gsk,nid[8],n,X,Y,Z,width;
	PetscInt *el;
	PetscMPIInt rank;
	PetscFunctionBegin;
	
	ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width!=1) {
		SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Stencil width must be 1 for Q1");
	}
	
	_npe = (order + 1)*(order + 1)*(order + 1);
  if (!da->e) {
		DMDAE dmdae;
		
		ierr = DMGetDMDAE(dm,&dmdae);CHKERRQ(ierr);
		mx = dmdae->lmx;
		my = dmdae->lmy;
		mz = dmdae->lmz;
		
		esi = dmdae->si;
		esj = dmdae->sj;
		esk = dmdae->sk;
		
		ierr = PetscMalloc(sizeof(PetscInt)*(mx*my*mz*_npe+1),&idx);CHKERRQ(ierr);
		ierr = DMDAGetGhostCorners(dm,&gsi,&gsj,&gsk, &X,&Y,&Z);CHKERRQ(ierr);
		
		elcnt = 0;
		for (ek=0; ek<mz; ek++) {
			k = (esk-gsk) + ek;

			for (ej=0; ej<my; ej++) {
				j = (esj-gsj) + ej;
			
				for (ei=0; ei<mx; ei++) {
					i = (esi-gsi) + ei;
				
					el = &idx[_npe*elcnt];
					
					nid[0] = (i  ) + (j  ) *X  + (k  ) *X*Y;
					nid[1] = (i+1) + (j  ) *X  + (k  ) *X*Y;
					
					nid[2] = (i  ) + (j+1) *X  + (k  ) *X*Y;
					nid[3] = (i+1) + (j+1) *X  + (k  ) *X*Y;
					
					nid[4] = (i  ) + (j  ) *X  + (k+1) *X*Y;
					nid[5] = (i+1) + (j  ) *X  + (k+1) *X*Y;
					
					nid[6] = (i  ) + (j+1) *X  + (k+1) *X*Y;
					nid[7] = (i+1) + (j+1) *X  + (k+1) *X*Y;

					//printf("e=%d: [%d %d %d %d] [%d %d %d %d] \n",elcnt,nid[0],nid[1],nid[2],nid[3],nid[0+4],nid[1+4],nid[2+4],nid[3+4] );
					
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
#define __FUNCT__ "DMDAGetElementsQ1"
PetscErrorCode DMDAGetElementsQ1(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx)
{
	PetscInt dim,sw;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(dm, &dim, 0,0,0, 0,0,0, 0,&sw, 0,0,0, 0);CHKERRQ(ierr);
	if (sw != 1) {
		SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Stencil width must equal 1");
	}
	switch (dim) {
		case 1:
			SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Dimension width must equal 3. Your DM has dim=1");
			break;
		case 2:
			SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Dimension width must equal 3. Your DM has dim=2");
			break;
		case 3:
			ierr = DMDAGetElements_DA_Q1_3D(dm,nel,npe,eidx);CHKERRQ(ierr);
			break;
	}
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMDASetElementType_Q1"
PetscErrorCode  DMDASetElementType_Q1(DM da)
{
  DM_DA          *dd = (DM_DA*)da->data;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscValidHeaderSpecific(da,DM_CLASSID,1);
  if (dd->elementtype) {
    ierr = PetscFree(dd->e);CHKERRQ(ierr);
    dd->ne          = 0; 
    dd->e           = NULL;
  }
	
  PetscFunctionReturn(0);
}

/* constructors */
#undef __FUNCT__
#define __FUNCT__ "DMDAProjectCoordinatesQ2toOverlappingQ1_3d"
PetscErrorCode DMDAProjectCoordinatesQ2toOverlappingQ1_3d(DM daq2,DM daq1)
{
	PetscErrorCode ierr;
	Vec coordsQ2, coordsQ1;
	DMDACoor3d ***LA_coordsQ2;
	DMDACoor3d ***LA_coordsQ1;
	PetscInt si1,sj1,sk1,nx1,ny1,nz1,i,j,k;
	PetscInt si2,sj2,sk2,nx2,ny2,nz2;
	DM cdaQ2,cdaQ1;
	
	PetscFunctionBegin;
	
	ierr = DMGetCoordinatesLocal(daq2,&coordsQ2);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(daq2,&cdaQ2);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cdaQ2,coordsQ2,&LA_coordsQ2);CHKERRQ(ierr);
	
	ierr = DMDASetUniformCoordinates(daq1,-2.0,2.0,-2.0,2.0,-2.0,2.0);CHKERRQ(ierr);
	ierr = DMGetCoordinates(daq1,&coordsQ1);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(daq1,&cdaQ1);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cdaQ1,coordsQ1,&LA_coordsQ1);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(     daq1,&si1,&sj1,&sk1 , &nx1,&ny1,&nz1);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(daq2,&si2,&sj2,&sk2 , &nx2,&ny2,&nz2);CHKERRQ(ierr);

	/*
	if ( (2*si1<si2) || (si1+nx1>si2+nx2) ) {
		printf("si1=%d, si2=%d, : si1+nx1=%d, si2+nx=%d \n", si1,si2, si1+nx1,si2+nx2 );
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-i) must overlap DA(Q1) in global space");
	}
	if ( (2*sj1<sj2) || (sj1+ny1>sj2+ny2) ) {
		printf("sj1=%d, sj2=%d, : sj1+ny1=%d, sj2+ny=%d \n", sj1,sj2, sj1+ny1,sj2+ny2 );
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-j) must overlap DA(Q1) in global space");
	}
	if ( (2*sk1<sk2) || (sk1+nz1>sk2+nz2) ) {
		printf("sk1=%d, sk2=%d, : sk1+nz1=%d, sk2+nz=%d \n", sk1,sk2, sk1+nz1,sk2+nz2 );
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-k) must overlap DA(Q1) in global space");
	}
*/	
	for( k=sk1; k<sk1+nz1; k++ ) {
		if ( (2*k<sk2) || (2*k>sk2+nz2) ) {
			PetscPrintf(PETSC_COMM_SELF,"sk1=%D, sk2=%D, : sk1+nz1=%D, sk2+nz=%D \n", sk1,sk2, sk1+nz1,sk2+nz2 );
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-k) must overlap DA(Q1) in global space");
		}
		
		for( j=sj1; j<sj1+ny1; j++ ) {
			if ( (2*j<sj2) || (2*j>sj2+ny2) ) {
				PetscPrintf(PETSC_COMM_SELF,"sj1=%D, sj2=%D, : sj1+ny1=%D, sj2+ny=%D \n", sj1,sj2, sj1+ny1,sj2+ny2 );
				SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-j) must overlap DA(Q1) in global space");
			}
			
			for( i=si1; i<si1+nx1; i++ ) {
				if ( (2*i<si2) || (2*i>si2+nx2) ) {
					PetscPrintf(PETSC_COMM_SELF,"si1=%D, si2=%D, : si1+nx1=%D, si2+nx=%D \n", si1,si2, si1+nx1,si2+nx2 );
					SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-i) must overlap DA(Q1) in global space");
				}
				
				//printf("(i,j,k-%d,%d,%d) => (I,J,K-%d,%d,%d)\n",i,j,k,2*i,2*j,2*k);
				LA_coordsQ1[k][j][i].x = LA_coordsQ2[2*k][2*j][2*i].x;
				LA_coordsQ1[k][j][i].y = LA_coordsQ2[2*k][2*j][2*i].y;
				LA_coordsQ1[k][j][i].z = LA_coordsQ2[2*k][2*j][2*i].z;
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(cdaQ1,coordsQ1,&LA_coordsQ1);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cdaQ2,coordsQ2,&LA_coordsQ2);CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(daq1);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DMDACreateOverlappingQ1FromQ2"
PetscErrorCode DMDACreateOverlappingQ1FromQ2(DM dmq2,PetscInt ndofs,DM *dmq1)
{
	DM dm;
	DMDAE dae;
	PetscInt stencilq2,sei,sej,sek,lmx,lmy,lmz,MX,MY,MZ,Mp,Np,Pp,i,j,k;
	PetscInt *siq2,*sjq2,*skq2,*lmxq2,*lmyq2,*lmzq2,*lxq1,*lyq1,*lzq1;
	PetscInt *lsip,*lsjp,*lskp;
	PetscErrorCode ierr;
	PetscMPIInt rank;
	
	PetscFunctionBegin;
	
	ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dmq2),&rank);CHKERRQ(ierr);
	
	ierr = DMDAGetCornersElementQ2(dmq2,&sei,&sej,&sek,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	ierr = DMDAGetSizeElementQ2(dmq2,&MX,&MY,&MZ);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dmq2,0, 0,0,0, &Mp,&Np,&Pp,0,&stencilq2, 0,0,0, 0);CHKERRQ(ierr);
	if (stencilq2 != 2) {
		SETERRQ(PetscObjectComm((PetscObject)dmq2),PETSC_ERR_SUP,"Attempting to create an overlapping DMDA from a DMDA which doesn't have a stencil width of 2... probably the generator isn't a Q2");
	}
	ierr = DMDAGetOwnershipRangesElementQ2(dmq2,0,0,0,&siq2,&sjq2,&skq2,&lmxq2,&lmyq2,&lmzq2);CHKERRQ(ierr);
	
	
	ierr = PetscMalloc(sizeof(PetscInt)*Mp,&lsip);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Np,&lsjp);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Pp,&lskp);CHKERRQ(ierr);
	
	for (i=0; i<Mp; i++) {
		lsip[i] = siq2[i]/2;
	}
	
	for (j=0; j<Np; j++) {
		lsjp[j] = sjq2[j]/2;
	}

	for (k=0; k<Pp; k++) {
		lskp[k] = skq2[k]/2;
	}
	
	ierr = PetscMalloc(sizeof(PetscInt)*Mp,&lxq1);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Np,&lyq1);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Pp,&lzq1);CHKERRQ(ierr);
	
	for (i=0; i<Mp-1; i++) {
		lxq1[i] = siq2[i+1]/2 - siq2[i]/2;
	}
	lxq1[Mp-1] = (2*MX - siq2[Mp-1])/2;
	lxq1[Mp-1]++;
	
	for (j=0; j<Np-1; j++) {
		lyq1[j] = sjq2[j+1]/2 - sjq2[j]/2;
	}
	lyq1[Np-1] = (2*MY - sjq2[Np-1])/2;
	lyq1[Np-1]++;

	for (k=0; k<Pp-1; k++) {
		lzq1[k] = skq2[k+1]/2 - skq2[k]/2;
	}
	lzq1[Pp-1] = (2*MZ - skq2[Pp-1])/2;
	lzq1[Pp-1]++;
	
	/*
	if (rank==0) {
		for (i=0; i<Mp; i++) {
			PetscPrintf(PETSC_COMM_SELF,"rank[%D] startI = %D: lmxq1 = %D \n",i,lsip[i],lxq1[i]);
		}
		for (j=0; j<Np; j++) {
			PetscPrintf(PETSC_COMM_SELF,"rank[%D] startJ = %D: lmyq1 = %D \n",j,lsjp[j],lyq1[j]);
		}
		for (k=0; k<Pp; k++) {
			PetscPrintf(PETSC_COMM_SELF,"rank[%D] startK = %D: lmzq1 = %D \n",k,lskp[k],lzq1[k]);
		}
	}
	*/
	
	ierr = DMDACreate3d(PetscObjectComm((PetscObject)dmq2),DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX+1,MY+1,MZ+1, Mp,Np,Pp, ndofs,1, lxq1,lyq1,lzq1, &dm );CHKERRQ(ierr);
	
	//ierr = DMView(dmq2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	//ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	/* add the space for the data structure */
	ierr = DMAttachDMDAE(dm);CHKERRQ(ierr);
	/* fetch the data structure */
	ierr = DMGetDMDAE(dm,&dae);CHKERRQ(ierr);
	dae->ne = MX * MY * MZ;
	dae->lne = lmx * lmy * lmz;
	
	dae->mx = MX;
	dae->my = MY;
	dae->mz = MZ;
	dae->lmx = lmx;
	dae->lmy = lmy;
	dae->lmz = lmz;
	
	dae->si = sei/2;
	dae->sj = sej/2;
	dae->sk = sek/2;
	
	dae->npe = 8;
	dae->nps = 2;
	dae->overlap = 0;
	
	for (i=0; i<Mp; i++) {
		lxq1[i] = lmxq2[i];
	}
	for (j=0; j<Np; j++) {
		lyq1[j] = lmyq2[j];
	}
	for (k=0; k<Pp; k++) {
		lzq1[k] = lmzq2[k];
	}
	dae->lmxp = lxq1;
	dae->lmyp = lyq1;
	dae->lmzp = lzq1;
	
	dae->lsip = lsip;
	dae->lsjp = lsjp;
	dae->lskp = lskp;
	
	PetscFree(siq2);
	PetscFree(sjq2);
	PetscFree(skq2);
	PetscFree(lmxq2);
	PetscFree(lmyq2);
	PetscFree(lmzq2);
	
	/* force element creation using MY numbering */
	{
		PetscInt nel,nen;
		const PetscInt *els;

		ierr = DMDASetElementType_Q1(dm);CHKERRQ(ierr);
		ierr = DMDAGetElements_DA_Q1_3D(dm,&nel,&nen,&els);CHKERRQ(ierr);
	}	
	
	/* force coordinate copy */
	ierr = DMDAProjectCoordinatesQ2toOverlappingQ1_3d(dmq2,dm);CHKERRQ(ierr);

	*dmq1 = dm;
	
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DMDACreateNestedQ1FromQ2"
PetscErrorCode DMDACreateNestedQ1FromQ2(DM dmq2,PetscInt ndofs,DM *dmq1)
{
	DM dm;
	DMDAE dae;
	PetscInt stencilq2,sei,sej,sek,lmx,lmy,lmz,MX,MY,MZ,Mp,Np,Pp,i,j,k;
	PetscInt *siq2,*sjq2,*skq2,*lmxq2,*lmyq2,*lmzq2,*lxq1,*lyq1,*lzq1;
	const PetscInt *lxq2,*lyq2,*lzq2;
	PetscInt *lsip,*lsjp,*lskp;
	PetscErrorCode ierr;
	PetscMPIInt rank;
	
	PetscFunctionBegin;
	
	ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dmq2),&rank);CHKERRQ(ierr);
	
	ierr = DMDAGetCornersElementQ2(dmq2,&sei,&sej,&sek,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	ierr = DMDAGetSizeElementQ2(dmq2,&MX,&MY,&MZ);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dmq2,0, 0,0,0, &Mp,&Np,&Pp,0,&stencilq2, 0,0,0, 0);CHKERRQ(ierr);
	if (stencilq2 != 2) {
		SETERRQ(PetscObjectComm((PetscObject)dmq2),PETSC_ERR_SUP,"Attempting to create an overlapping DMDA from a DMDA which doesn't have a stencil width of 2... probably the generator isn't a Q2");
	}
	
	ierr = DMDAGetOwnershipRangesElementQ2(dmq2,0,0,0,&siq2,&sjq2,&skq2,&lmxq2,&lmyq2,&lmzq2);CHKERRQ(ierr);
	
	
	ierr = PetscMalloc(sizeof(PetscInt)*Mp,&lsip);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Np,&lsjp);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Pp,&lskp);CHKERRQ(ierr);
	
	for (i=0; i<Mp; i++) {
		lsip[i] = siq2[i];
		PetscPrintf(PETSC_COMM_SELF,"[%D]: si[%D] %D \n", rank,i,lsip[i]);
	}
	
	for (j=0; j<Np; j++) {
		lsjp[j] = sjq2[j];
		PetscPrintf(PETSC_COMM_SELF,"[%D]: sj[%D] %D \n", rank,j,lsjp[j]);
	}
	
	for (k=0; k<Pp; k++) {
		lskp[k] = skq2[k];
		PetscPrintf(PETSC_COMM_SELF,"[%D]: sk[%D] %D \n", rank,k,lskp[k]);
	}

/*
	ierr = PetscMalloc(sizeof(PetscInt)*Mp,&lxq1);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Np,&lyq1);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Pp,&lzq1);CHKERRQ(ierr);
	
	for (i=0; i<Mp-1; i++) {
		lxq1[i] = siq2[i+1] - siq2[i];
	}
	lxq1[Mp-1] = 2*MX+1 - siq2[Mp-1];
	
	for (j=0; j<Np-1; j++) {
		lyq1[j] = sjq2[j+1] - sjq2[j];
	}
	lyq1[Np-1] = 2*MY+1 - sjq2[Np-1];
	
	for (k=0; k<Pp-1; k++) {
		lzq1[k] = skq2[k+1] - skq2[k];
		printf("[%d]: lz_k[%d] %d \n", rank,k,lzq1[k]);
	}
	lzq1[Pp-1] = 2*MZ+1 - skq2[Pp-1];
	printf("[%d]: lz_k[%d] %d \n", rank,Pp-1,lzq1[Pp-1]);
*/	
	ierr = DMDAGetOwnershipRanges(dmq2,&lxq2,&lyq2,&lzq2);CHKERRQ(ierr);

	ierr = PetscMalloc(sizeof(PetscInt)*Mp,&lxq1);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Np,&lyq1);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Pp,&lzq1);CHKERRQ(ierr);
	for (i=0; i<Mp; i++) {
		lxq1[i] = lxq2[i];
	}
	
	for (j=0; j<Np; j++) {
		lyq1[j] = lyq2[j];
	}
	
	for (k=0; k<Pp; k++) {
		lzq1[k] = lzq2[k];
	}
	
	if (rank==0) {
		for (i=0; i<Mp; i++) {
			PetscPrintf(PETSC_COMM_SELF,"rank[%D] startI = %D: lmxq1 = %D \n",i,lsip[i],lxq1[i]);
		}
		for (j=0; j<Np; j++) {
			PetscPrintf(PETSC_COMM_SELF,"rank[%D] startJ = %D: lmyq1 = %D \n",j,lsjp[j],lyq1[j]);
		}
		for (k=0; k<Pp; k++) {
			PetscPrintf(PETSC_COMM_SELF,"rank[%D] startK = %D: lmzq1 = %D \n",k,lskp[k],lzq1[k]);
		}
	}
	
	ierr = DMDACreate3d(PetscObjectComm((PetscObject)dmq2),DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, 2*MX+1,2*MY+1,2*MZ+1, Mp,Np,Pp, ndofs,1, lxq2,lyq2,lzq2, &dm );CHKERRQ(ierr);
	
	//ierr = DMView(dmq2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	//ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	/* add the space for the data structure */
	ierr = DMAttachDMDAE(dm);CHKERRQ(ierr);
	/* fetch the data structure */
	ierr = DMGetDMDAE(dm,&dae);CHKERRQ(ierr);
	dae->ne = 8 * MX * MY * MZ;
	dae->lne = 8 * lmx * lmy * lmz;
	
	dae->mx = 2*MX;
	dae->my = 2*MY;
	dae->mz = 2*MZ;
	dae->lmx = 2*lmx;
	dae->lmy = 2*lmy;
	dae->lmz = 2*lmz;
	
	dae->si = sei;
	dae->sj = sej;
	dae->sk = sek;
	
	dae->npe = 8;
	dae->nps = 2;
	dae->overlap = 0;
	
	/* for convience we scale these and re-use the data */
	for (i=0; i<Mp; i++) {
		lxq1[i] = 2*lmxq2[i];
	}
	for (j=0; j<Np; j++) {
		lyq1[j] = 2*lmyq2[j];
	}
	for (k=0; k<Pp; k++) {
		lzq1[k] = 2*lmzq2[k];
	}
	dae->lmxp = lxq1;
	dae->lmyp = lyq1;
	dae->lmzp = lzq1;
	
	dae->lsip = lsip;
	dae->lsjp = lsjp;
	dae->lskp = lskp;
	
	PetscFree(siq2);
	PetscFree(sjq2);
	PetscFree(skq2);
	PetscFree(lmxq2);
	PetscFree(lmyq2);
	PetscFree(lmzq2);
	
	/* force element creation using MY numbering */
	{
		PetscInt nel,nen;
		const PetscInt *els;

		ierr = DMDASetElementType_Q1(dm);CHKERRQ(ierr);
		ierr = DMDAGetElements_DA_Q1_3D(dm,&nel,&nen,&els);CHKERRQ(ierr);
	}	
	
	/* force coordinate copy */
	{
		Vec coordq1,coordq2;
		
		ierr = DMDASetUniformCoordinates(dm,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);

		ierr = DMGetCoordinates(dmq2,&coordq2);CHKERRQ(ierr);
		ierr = DMGetCoordinates(dm,&coordq1);CHKERRQ(ierr);
		ierr = VecCopy(coordq2,coordq1);CHKERRQ(ierr);
		ierr = DMDAUpdateGhostedCoordinates(dm);CHKERRQ(ierr);
	}

	*dmq1 = dm;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DMDACreateQ1"
PetscErrorCode DMDACreateQ1(MPI_Comm comm,PetscInt MX,PetscInt MY,PetscInt MZ,PetscInt ndofs,DM *dmq1)
{
	DM dm;
	DMDAE dae;
	PetscInt sei,sej,sek,lmx,lmy,lmz;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDACreate3d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX+1,MY+1,MZ+1, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, ndofs,1, 0,0,0, &dm );CHKERRQ(ierr);
	
	/* add the space for the data structure */
	ierr = DMAttachDMDAE(dm);CHKERRQ(ierr);
	/* fetch the data structure */
	ierr = DMGetDMDAE(dm,&dae);CHKERRQ(ierr);

	dae->ne = MX * MY * MZ;
	dae->mx = MX;
	dae->my = MY;
	dae->mz = MZ;

	ierr = DMDAGetLocalSizeElementQ1(dm,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	dae->lne = lmx * lmy * lmz;	
	dae->lmx = lmx;
	dae->lmy = lmy;
	dae->lmz = lmz;
	
	ierr = DMDAGetCornersElementQ1(dm,&sei,&sej,&sek,0,0,0);CHKERRQ(ierr);
	dae->si = sei;
	dae->sj = sej;
	dae->sk = sek;
	
	dae->npe = 8;
	dae->nps = 2;
	dae->overlap = 0;
	
	/* for convience we scale these and re-use the data */
	dae->lmxp = NULL;
	dae->lmyp = NULL;
	dae->lmzp = NULL;
	
	dae->lsip = NULL;
	dae->lsjp = NULL;
	dae->lskp = NULL;
	
	/* force element creation using MY numbering */
	{
		PetscInt nel,nen;
		const PetscInt *els;
		
		ierr = DMDASetElementType_Q1(dm);CHKERRQ(ierr);
		ierr = DMDAGetElements_DA_Q1_3D(dm,&nel,&nen,&els);CHKERRQ(ierr);
	}	
	
	/* force coordinate copy */
	ierr = DMDASetUniformCoordinates(dm,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dm);CHKERRQ(ierr);
	
	*dmq1 = dm;
	
	PetscFunctionReturn(0);
}

/* element helpers */
#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1_GetElementCoordinates_3D"
PetscErrorCode DMDAEQ1_GetElementCoordinates_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[])
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
#define __FUNCT__ "DMDAEQ1_GetScalarElementField_3D"
PetscErrorCode DMDAEQ1_GetScalarElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
{
	PetscInt n;

	PetscFunctionBegin;
	for (n=0; n<8; n++) {
		elfield[n] = LA_gfield[elnid[n]];
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1_GetVectorElementField_3D"
PetscErrorCode DMDAEQ1_GetVectorElementField_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
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
#define __FUNCT__ "DMDAEQ1_SetValuesLocalStencil_AddValues_DOF"
PetscErrorCode DMDAEQ1_SetValuesLocalStencil_AddValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[])
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
#define __FUNCT__ "DMDAEQ1_GetElementLocalIndicesDOF"
PetscErrorCode DMDAEQ1_GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[])
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
#define __FUNCT__ "DMDAProjectVectorQ2toOverlappingQ1_3d"
PetscErrorCode DMDAProjectVectorQ2toOverlappingQ1_3d(DM daq2,Vec x2,DM daq1,Vec x1)
{
	PetscErrorCode ierr;
	PetscInt       d,ndofq2,ndofq1;
	PetscScalar    ****LA_fieldQ2;
	PetscScalar    ****LA_fieldQ1;
	PetscInt       si1,sj1,sk1,nx1,ny1,nz1,i,j,k;
	PetscInt       si2,sj2,sk2,nx2,ny2,nz2;
	Vec            x2_local;
	
	PetscFunctionBegin;

	ierr = DMDAGetInfo(daq2,0, 0,0,0, 0,0,0, &ndofq2,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetInfo(daq1,0, 0,0,0, 0,0,0, &ndofq1,0, 0,0,0, 0);CHKERRQ(ierr);
	
	/* scatter Q2 (ghosted field) */
	ierr = DMGetLocalVector(daq2,&x2_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(daq2,x2,INSERT_VALUES,x2_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (daq2,x2,INSERT_VALUES,x2_local);CHKERRQ(ierr);
	
	ierr = DMDAVecGetArrayDOF(daq2,x2_local,&LA_fieldQ2);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(daq1,x1,&LA_fieldQ1);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(     daq1,&si1,&sj1,&sk1 , &nx1,&ny1,&nz1);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(daq2,&si2,&sj2,&sk2 , &nx2,&ny2,&nz2);CHKERRQ(ierr);
	
	for( k=sk1; k<sk1+nz1; k++ ) {
		if ( (2*k<sk2) || (2*k>sk2+nz2) ) {
			PetscPrintf(PETSC_COMM_SELF,"sk1=%D, sk2=%D, : sk1+nz1=%D, sk2+nz=%D \n", sk1,sk2, sk1+nz1,sk2+nz2 );
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-k) must overlap DA(Q1) in global space");
		}
		
		for( j=sj1; j<sj1+ny1; j++ ) {
			if ( (2*j<sj2) || (2*j>sj2+ny2) ) {
				PetscPrintf(PETSC_COMM_SELF,"sj1=%D, sj2=%D, : sj1+ny1=%D, sj2+ny=%D \n", sj1,sj2, sj1+ny1,sj2+ny2 );
				SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-j) must overlap DA(Q1) in global space");
			}
			
			for( i=si1; i<si1+nx1; i++ ) {
				if ( (2*i<si2) || (2*i>si2+nx2) ) {
					PetscPrintf(PETSC_COMM_SELF,"si1=%D, si2=%D, : si1+nx1=%D, si2+nx=%D \n", si1,si2, si1+nx1,si2+nx2 );
					SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DA(Q2) (ghosted-i) must overlap DA(Q1) in global space");
				}
				
				//printf("(i,j,k-%d,%d,%d) => (I,J,K-%d,%d,%d)\n",i,j,k,2*i,2*j,2*k);
				for (d=0; d<ndofq2; d++) {
					LA_fieldQ1[k][j][i][d] = LA_fieldQ2[2*k][2*j][2*i][d];
				}
			}
		}
	}
	
	ierr = DMDAVecRestoreArrayDOF(daq1,x1,&LA_fieldQ1);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(daq2,x2_local,&LA_fieldQ2);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(daq2,&x2_local);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAProjectVectorQ2toNestedQ1_3d"
PetscErrorCode DMDAProjectVectorQ2toNestedQ1_3d(DM daq2,Vec x2,DM daq1,Vec x1)
{
	PetscErrorCode ierr;
	PetscInt si1,sj1,sk1,nx1,ny1,nz1;
	PetscInt si2,sj2,sk2,nx2,ny2,nz2;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetCorners(daq1,&si1,&sj1,&sk1 , &nx1,&ny1,&nz1);CHKERRQ(ierr);
	ierr = DMDAGetCorners(daq2,&si2,&sj2,&sk2 , &nx2,&ny2,&nz2);CHKERRQ(ierr);

	if (si1 != si2){ SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"index does not match (si)"); }
	if (sj1 != sj2){ SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"index does not match (sj)"); }
	if (sk1 != sk2){ SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"index does not match (sk)"); }

	if (nx1 != nz2){ SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"index does not match (nx)"); }
	if (ny1 != ny2){ SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"index does not match (ny)"); }
	if (nz1 != nz2){ SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"index does not match (nz)"); }
		
	ierr = VecCopy(x2,x1);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAProjectVectorQ2toQ1"
PetscErrorCode DMDAProjectVectorQ2toQ1(DM daq2,Vec x2,DM daq1,Vec x1,PetscInt mesh_type)
{
	PetscErrorCode ierr;
	PetscInt dim,ndofq2,ndofq1;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(daq2,&dim, 0,0,0, 0,0,0, &ndofq2,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetInfo(daq1,0, 0,0,0, 0,0,0, &ndofq1,0, 0,0,0, 0);CHKERRQ(ierr);
	if (ndofq1 != ndofq2) {
		SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"ndof on DM's do not match");
	}
	
	switch (dim) {
		case 1:
			SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"Only 3d supported");
			break;
		case 2:
			SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"Only 3d supported");
			break;
		case 3:
			if (mesh_type == 0) {
				SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"mesh type must be 1 or 2");
			} else if (mesh_type == 1) {
				ierr = DMDAProjectVectorQ2toOverlappingQ1_3d(daq2,x2,daq1,x1);CHKERRQ(ierr);
			} else if (mesh_type == 2) {
				ierr = DMDAProjectVectorQ2toNestedQ1_3d(daq2,x2,daq1,x1);CHKERRQ(ierr);
			} else {
				SETERRQ(PetscObjectComm((PetscObject)daq2),PETSC_ERR_SUP,"Unknown mesh type");
			}
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAProjectCoordinatesQ2toQ1"
PetscErrorCode DMDAProjectCoordinatesQ2toQ1(DM daq2,DM daq1,PetscInt mesh_type)
{
	PetscErrorCode ierr;
	DM cdaq2,cdaq1;
	Vec coordq2,coordq1;
	
	PetscFunctionBegin;
	
	ierr = DMGetCoordinateDM(daq2,&cdaq2);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(daq1,&cdaq1);CHKERRQ(ierr);

	ierr = DMGetCoordinates(daq2,&coordq2);CHKERRQ(ierr);
	ierr = DMGetCoordinates(daq1,&coordq1);CHKERRQ(ierr);
	
	ierr = DMDAProjectVectorQ2toQ1(cdaq2,coordq2,cdaq1,coordq1,mesh_type);CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(daq1);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

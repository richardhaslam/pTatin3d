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
 **    filename:   dmda_element_q2p1.c
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
#include <petsc/private/dmimpl.h>
#include <petsc/private/dmdaimpl.h>
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin3d_defs.h"
#include "ptatin3d_stokes.h"
#include "dmda_element_q2p1.h"
#include "ptable.h"


/* DA Q2 1D,2D,3D */
PetscErrorCode DMDAGetSizeElementQ2(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ)
{
  const PetscInt order = 2;
  PetscInt M,N,P,width;
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
PetscErrorCode DMDAGetLocalSizeElementQ2(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz)
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
//      PetscPrintf(PETSC_COMM_SELF,"ELEMENT(i) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,si,si+m-1 );
      cntx++;
      continue;
    }

    if (si+m-n2>1) {
//      PetscPrintf(PETSC_COMM_SELF,"ELEMENT(i) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,si,si+m-1 );
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
//      PetscPrintf(PETSC_COMM_SELF,"ELEMENT(j) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,sj,sj+n-1 );
      cnty++;
      continue;
    }

    if (sj+n-n2>1) {
//      PetscPrintf(PETSC_COMM_SELF,"ELEMENT(j) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,sj,sj+n-1 );
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

    //    printf("last-n2 = %d \n", sk+p-n2 );
    /* if start and end of element are inside global range - keep it */
    if (n2<skg+pg) {
//      PetscPrintf(PETSC_COMM_SELF,"ELEMENT(k) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,sk,sk+p-1 );
      //      printf("[GLOBAL] usng start id [k] %d [%d-%d]l [%d-%d]g\n", n0,sk,sk+p-1,skg,skg+pg-1);
      cntz++;
      continue;
    }

    if (sk+p-n2>1) {
//      PetscPrintf(PETSC_COMM_SELF,"ELEMENT(k) (%6d - %6d) inside range l[%6d - %6d]\n", n0,n2,sk,sk+p-1 );
      //      printf("[LOCAL] usng start id [k] %d [%d-%d]l [%d-%d]g\n", n0,sk,sk+p-1,skg,skg+pg-1);
      cntz++;
      continue;
    }

    if (sk+p-n2<=1) {
      /* this means the element is taking two entries from the ghost cells */
      //      printf("[OUTSIDE] element (%d .. %d) [%d-%d]l [%d-%d]g\n", n0,n2,sk,sk+p-1,skg,skg+pg-1);
      break;
    }
  }

  /* ======================================================================================== */

  if (mx) { *mx = cntx; }
  if (my) { *my = cnty; }
  if (mz) { *mz = cntz; }
  PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
PetscErrorCode DMDAGetCornersElementQ2(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
  PetscInt i,j,k;
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
PetscErrorCode DMDAGetOwnershipRangesElementQ2(DM da,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt **si,PetscInt **sj,PetscInt **sk,PetscInt **_mx,PetscInt **_my,PetscInt **_mz)
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
  ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

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

PetscErrorCode DMDAGetElements_DA_Q2_3D(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx)
{
  DM_DA          *da = (DM_DA*)dm->data;
  const PetscInt order = 2;
  PetscErrorCode ierr;
  PetscInt *idx,mx,my,mz,_npe,M,N,P;
  PetscInt ei,ej,ek,i,j,k,elcnt,esi,esj,esk,gsi,gsj,gsk,nid[27],n,X,Y,Z,width;
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

          //        printf("[%d]: i=%d \n", rank,i );
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

  if (eidx) { *eidx = da->e; }
  if (npe)  {  *npe = _npe; }
  if (nel)  {  *nel = da->ne; }

  PetscFunctionReturn(0);
}

PetscErrorCode DMDAGetElements_DA_Q2(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[])
{
  PetscErrorCode ierr;
  PetscInt       dim;
  PetscFunctionBegin;
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  if (dim==-1) {
    *nel = 0; *nen = 0; *e = NULL;
  } else if (dim==1) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DMDA doesn't support Q2 in 1D");
  } else if (dim==2) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DMDA doesn't support Q2 in 2D");
  } else if (dim==3) {
    ierr = DMDAGetElements_DA_Q2_3D(dm,nel,nen,e);CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"DMDA dimension not 1, 2, or 3, it is %D\n",dim);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMDAGetElements_DA_P0MD_3D(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx)
{
  DM_DA          *da = (DM_DA*)dm->data;
  PetscErrorCode ierr;
  PetscInt *idx,mx,my,mz,_npe;
  PetscInt d,ei,ej,ek,elcnt,esi,esj,nid[100],width;
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

PetscErrorCode DMDAGetElements_DA_P1(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[])
{
  PetscErrorCode ierr;
  PetscInt       dim;
  PetscFunctionBegin;
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  if (dim==-1) {
    *nel = 0; *nen = 0; *e = NULL;
  } else if (dim==1) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DMDA doesn't support P1 in 1D");
  } else if (dim==2) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"DMDA doesn't support P1 in 2D");
  } else if (dim==3) {
    ierr = DMDAGetElements_DA_P0MD_3D(dm,nel,nen,e);CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"DMDA dimension not 1, 2, or 3, it is %D\n",dim);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMDAGetElements_pTatinQ2P1(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[])
{
  PetscErrorCode ierr;
  PetscInt dof,sw;

  ierr = DMDAGetInfo(dm, 0, 0,0,0, 0,0,0, &dof,&sw, 0,0,0, 0);CHKERRQ(ierr);

  if (sw == 2) {
    ierr = DMDAGetElements_DA_Q2(dm,nel,nen,e);CHKERRQ(ierr);
  } else if (sw == 0) {
    ierr = DMDAGetElements_DA_P1(dm,nel,nen,e);CHKERRQ(ierr);
  } else {
    ierr = DMDAGetElements(dm,nel,nen,e);CHKERRQ(ierr);
  }

    PetscFunctionReturn(0);
}

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
PetscErrorCode DMDAGetElementCoordinatesQ2_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[])
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

PetscErrorCode DMDAGetScalarElementFieldQ2_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
{
  PetscInt n;

  PetscFunctionBegin;
  for (n=0; n<27; n++) {
    elfield[n] = LA_gfield[elnid[n]];
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DMDAGetVectorElementFieldQ2_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[])
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

PetscErrorCode DMDAGetScalarElementField(PetscScalar elfield[],PetscInt npe,PetscInt elnid[],PetscScalar LA_gfield[])
{
  PetscInt n;

  PetscFunctionBegin;
  for (n=0; n<npe; n++) {
    elfield[n] = LA_gfield[elnid[n]];
  }
  PetscFunctionReturn(0);
}

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

PetscErrorCode Q2GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[])
{
  PetscInt n,d;
  PetscFunctionBegin;
  for (d=0; d<ndof; d++) {
    for (n=0; n<U_BASIS_FUNCTIONS; n++) {
      el_localIndices[ndof*n+d] = ndof*elnid[n]+d;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateMatrix_DAQ2(DM dm,Mat *B)
{
  PetscErrorCode         ierr;
  PTable                 table;
  PetscInt               start,end,M,m,e,i,j,nelements,nen_u;
  const PetscInt         *elnidx_u;
  PetscInt               num_gindices,ge_eqnums[3*Q2_NODES_PER_EL_3D];
  PetscInt               vel_el_lidx[3*U_BASIS_FUNCTIONS];
  Vec                    xg;
  PetscInt               *dnnz = NULL,*onnz = NULL;
  ISLocalToGlobalMapping ltog;
  const PetscInt         *gindices;
  Mat                    A;
  
  ierr = DMCreateGlobalVector(dm,&xg);CHKERRQ(ierr);
  ierr = VecGetSize(xg,&M);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(xg,&start,&end);CHKERRQ(ierr);
  m = end - start;
  ierr = VecDestroy(&xg);CHKERRQ(ierr);
  
  ierr = PTableCreate(PETSC_COMM_WORLD,&table);CHKERRQ(ierr);
  ierr = PTableSetRange(table,start,end);CHKERRQ(ierr);
  ierr = PTableSetType(table,PTABLE_DENSE);CHKERRQ(ierr);
  ierr = PTableSetup(table);CHKERRQ(ierr);
  
  ierr = DMDAGetElements_pTatinQ2P1(dm,&nelements,&nen_u,&elnidx_u);CHKERRQ(ierr);
  ierr = DMGetLocalToGlobalMapping(dm, &ltog);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingGetSize(ltog, &num_gindices);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingGetIndices(ltog, &gindices);CHKERRQ(ierr);
  
  for (e=0; e<nelements; e++) {
    
    /* get local indices */
    ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
    
    /* get global indices */
    for (i=0; i<nen_u; i++) {
      const int nid = elnidx_u[nen_u*e + i];
      
      ge_eqnums[3*i  ] = gindices[ 3*nid   ];
      ge_eqnums[3*i+1] = gindices[ 3*nid+1 ];
      ge_eqnums[3*i+2] = gindices[ 3*nid+2 ];
    }
    
    for (i=0; i<3*nen_u; i++) {
      for (j=0; j<3*nen_u; j++) {
        PetscInt row,col;
        
        row = ge_eqnums[i];
        col = ge_eqnums[j];
        ierr = PTableSetValue(table,(long int)row,(long int)col);CHKERRQ(ierr);
      }
    }
    
  }
  ierr = PTableSynchronize(table);CHKERRQ(ierr);
  
  /* determine non-zero structure using PTable */
  
  ierr = PetscMalloc1((end-start),&dnnz);CHKERRQ(ierr);
  ierr = PetscMalloc1((end-start),&onnz);CHKERRQ(ierr);
  
  for (i=start; i<end; i++) {
    long int rl;
    const long int *vals;
    
    dnnz[i-start] = 0;
    onnz[i-start] = 0;
    
    ierr = PTableGetValues(table,i,&rl,&vals);CHKERRQ(ierr);
    
    for (j=0; j<rl; j++) {
      PetscInt col = (PetscInt)vals[j];
      
      if (col >= start && col < end) {
        dnnz[i-start]++;
      } else {
        onnz[i-start]++;
      }
    }
  }
  
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,m,m,M,M);CHKERRQ(ierr);
  ierr = MatSetBlockSize(A,3);CHKERRQ(ierr);
  ierr = MatSetType(A,MATAIJ);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  
  ierr = MatSeqAIJSetPreallocation(A,0,dnnz);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,0,dnnz,0,onnz);CHKERRQ(ierr);
  
  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  
  /* push non-zero structure using PTable */
  for (i=start; i<end; i++) {
    long int rl;
    const long int *vals;
    
    dnnz[i-start] = 0;
    onnz[i-start] = 0;
    
    ierr = PTableGetValues(table,i,&rl,&vals);CHKERRQ(ierr);
    
    for (j=0; j<rl; j++) {
      PetscInt col = (PetscInt)vals[j];
      
      ierr = MatSetValue(A,i,col,0.0,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  
  ierr = ISLocalToGlobalMappingRestoreIndices(ltog, &gindices);CHKERRQ(ierr);
  ierr = PetscFree(dnnz);CHKERRQ(ierr);
  ierr = PetscFree(onnz);CHKERRQ(ierr);
  ierr = PTableDestroy(&table);CHKERRQ(ierr);
  
  *B = A;
  
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateMatrixAuuHelper(DM dm,PetscBool use_custom_mapping,Mat *B)
{
  PetscErrorCode ierr;
  if (use_custom_mapping) {
    ierr = DMCreateMatrix_DAQ2(dm,B);CHKERRQ(ierr);
  } else {
    ierr = DMCreateMatrix(dm,B);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

typedef struct {
  PetscInt   refcnt;
  VecScatter global2local;
  IS         from,to;
  PetscInt   *elements_default;
  PetscInt   *elements_custom;
  PetscInt   nlocal,nelements,mx,my,mz;
  /* reference function pointers to the standard / default DMDA implementations */
  PetscErrorCode (*createlocalvector_default)(DM,Vec*);
  PetscErrorCode (*creatematrix_default)(DM, Mat*);
  PetscErrorCode (*globaltolocalbegin_default)(DM,Vec,InsertMode,Vec);
  PetscErrorCode (*globaltolocalend_default)(DM,Vec,InsertMode,Vec);
  PetscErrorCode (*localtoglobalbegin_default)(DM,Vec,InsertMode,Vec);
  PetscErrorCode (*localtoglobalend_default)(DM,Vec,InsertMode,Vec);
  /* reference function pointers to the standard / default Coordinate DMDA implementations */
  PetscErrorCode (*createlocalvector_default_c)(DM,Vec*);
  PetscErrorCode (*creatematrix_default_c)(DM, Mat*);
  PetscErrorCode (*globaltolocalbegin_default_c)(DM,Vec,InsertMode,Vec);
  PetscErrorCode (*globaltolocalend_default_c)(DM,Vec,InsertMode,Vec);
  PetscErrorCode (*localtoglobalbegin_default_c)(DM,Vec,InsertMode,Vec);
  PetscErrorCode (*localtoglobalend_default_c)(DM,Vec,InsertMode,Vec);
} DMDAQ2Mapping;

PetscErrorCode DMDAQ2MappingDestroy(DMDAQ2Mapping **m)
{
  DMDAQ2Mapping  *map;
  PetscErrorCode ierr;
  
  if (!m) PetscFunctionReturn(0);
  map = *m;
  if (--(map->refcnt) > 0) { *m = NULL; PetscFunctionReturn(0); }
  
  ierr = VecScatterDestroy(&map->global2local);CHKERRQ(ierr);
  ierr = ISDestroy(&map->from);CHKERRQ(ierr);
  ierr = ISDestroy(&map->to);CHKERRQ(ierr);
  ierr = PetscFree(map->elements_default);CHKERRQ(ierr);
  ierr = PetscFree(map->elements_custom);CHKERRQ(ierr);
  ierr = PetscFree(map);CHKERRQ(ierr);
  *m = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode PetscContainer_DMDAQ2MappingDestroy(void *ctx)
{
  PetscErrorCode ierr;
  DMDAQ2Mapping  *map = (DMDAQ2Mapping*)ctx;
  ierr = DMDAQ2MappingDestroy(&map);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* dmda over-rides */
PetscErrorCode DMGlobalToLocalBegin_DAQ2(DM da,Vec g,InsertMode mode,Vec l)
{
  PetscErrorCode ierr;
  PetscContainer c = NULL;
  DMDAQ2Mapping  *mapping = NULL;
  
  PetscFunctionBegin;
  PetscObjectQuery((PetscObject)da,"DMDACustomMapping",(PetscObject*)&c);
  PetscContainerGetPointer(c,(void**)&mapping);
  //printf("custom g2l-b\n");
  ierr = VecScatterBegin(mapping->global2local,g,l,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalEnd_DAQ2(DM da,Vec g,InsertMode mode,Vec l)
{
  PetscErrorCode ierr;
  PetscContainer c = NULL;
  DMDAQ2Mapping  *mapping = NULL;
  
  PetscFunctionBegin;
  PetscObjectQuery((PetscObject)da,"DMDACustomMapping",(PetscObject*)&c);
  PetscContainerGetPointer(c,(void**)&mapping);
  //printf("custom g2l-e\n");
  ierr = VecScatterEnd(mapping->global2local,g,l,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalBegin_DAQ2(DM da,Vec l,InsertMode mode,Vec g)
{
  PetscErrorCode ierr;
  PetscContainer c = NULL;
  DMDAQ2Mapping  *mapping = NULL;
  
  PetscFunctionBegin;
  PetscObjectQuery((PetscObject)da,"DMDACustomMapping",(PetscObject*)&c);
  PetscContainerGetPointer(c,(void**)&mapping);
  //printf("custom l2g-b\n");
  if (mode == ADD_VALUES) {
    ierr = VecScatterBegin(mapping->global2local,l,g,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  } else if (mode == INSERT_VALUES) {
    ierr = VecScatterBegin(mapping->global2local,l,g,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Not yet implemented");
  PetscFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalEnd_DAQ2(DM da,Vec l,InsertMode mode,Vec g)
{
  PetscErrorCode ierr;
  PetscContainer c = NULL;
  DMDAQ2Mapping  *mapping = NULL;
  
  PetscFunctionBegin;
  PetscObjectQuery((PetscObject)da,"DMDACustomMapping",(PetscObject*)&c);
  PetscContainerGetPointer(c,(void**)&mapping);
  //printf("custom l2g-2\n");
  if (mode == ADD_VALUES) {
    ierr = VecScatterEnd(mapping->global2local,l,g,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  } else if (mode == INSERT_VALUES) {
    ierr = VecScatterEnd(mapping->global2local,l,g,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Not yet implemented");
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateLocalVector_DAQ2(DM da,Vec *g)
{
  PetscErrorCode ierr;
  PetscContainer c = NULL;
  DMDAQ2Mapping  *mapping = NULL;
  PetscInt       length;
  
  PetscFunctionBegin;
  PetscObjectQuery((PetscObject)da,"DMDACustomMapping",(PetscObject*)&c);
  PetscContainerGetPointer(c,(void**)&mapping);
  
  //printf("custom CreateLocal\n");
  length = mapping->nlocal * 3;
  ierr = VecCreate(PETSC_COMM_SELF,g);CHKERRQ(ierr);
  ierr = VecSetSizes(*g,length,length);CHKERRQ(ierr);
  ierr = VecSetBlockSize(*g,3);CHKERRQ(ierr);
  ierr = VecSetType(*g,VECSEQ);CHKERRQ(ierr);
  ierr = VecSetDM(*g, da);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMDAQ2ElementLayout(PetscInt eidx[],PetscInt mx,PetscInt my,PetscInt mz)
{
  PetscInt       elcnt,n,i,j,k,ei,ej,ek,X,Y,nid[27];
  PetscInt       *el;
  PetscErrorCode ierr;
  
  ierr = PetscMemzero(eidx,sizeof(PetscInt)*(mx*my*mz)*27);CHKERRQ(ierr);
  
  elcnt = 0;
  X = 2*mx + 1;
  Y = 2*my + 1;
  
  for (ek=0; ek<mz; ek++) {
    k = 2*ek;
    for (ej=0; ej<my; ej++) {
      j = 2*ej;
      for (ei=0; ei<mx; ei++) {
        i = 2*ei;
        
        el = &eidx[27*elcnt];
        
        nid[ 0] = (i  ) + (j  ) *X  + (k  ) *X*Y;
        nid[ 1] = (i+1) + (j  ) *X  + (k  ) *X*Y;
        nid[ 2] = (i+2) + (j  ) *X  + (k  ) *X*Y;
        
        nid[ 3] = (i  ) + (j+1) *X  + (k  ) *X*Y;
        nid[ 4] = (i+1) + (j+1) *X  + (k  ) *X*Y;
        nid[ 5] = (i+2) + (j+1) *X  + (k  ) *X*Y;
        
        nid[ 6] = (i  ) + (j+2) *X  + (k  ) *X*Y;
        nid[ 7] = (i+1) + (j+2) *X  + (k  ) *X*Y;
        nid[ 8] = (i+2) + (j+2) *X  + (k  ) *X*Y;
        
        nid[ 9] = (i  ) + (j  ) *X  + (k+1) *X*Y;
        nid[10] = (i+1) + (j  ) *X  + (k+1) *X*Y;
        nid[11] = (i+2) + (j  ) *X  + (k+1) *X*Y;
        
        nid[12] = (i  ) + (j+1) *X  + (k+1) *X*Y;
        nid[13] = (i+1) + (j+1) *X  + (k+1) *X*Y;
        nid[14] = (i+2) + (j+1) *X  + (k+1) *X*Y;
        
        nid[15] = (i  ) + (j+2) *X  + (k+1) *X*Y;
        nid[16] = (i+1) + (j+2) *X  + (k+1) *X*Y;
        nid[17] = (i+2) + (j+2) *X  + (k+1) *X*Y;
        
        nid[18] = (i  ) + (j  ) *X  + (k+2) *X*Y;
        nid[19] = (i+1) + (j  ) *X  + (k+2) *X*Y;
        nid[20] = (i+2) + (j  ) *X  + (k+2) *X*Y;
        
        nid[21] = (i  ) + (j+1) *X  + (k+2) *X*Y;
        nid[22] = (i+1) + (j+1) *X  + (k+2) *X*Y;
        nid[23] = (i+2) + (j+1) *X  + (k+2) *X*Y;
        
        nid[24] = (i  ) + (j+2) *X  + (k+2) *X*Y;
        nid[25] = (i+1) + (j+2) *X  + (k+2) *X*Y;
        nid[26] = (i+2) + (j+2) *X  + (k+2) *X*Y;
        
        for (n=0; n<27; n++) {
          el[n] = nid[n];
        }
        
        elcnt++;
      }
    }
  }
  PetscFunctionReturn(0);
}

/*
 This might be a good idea for clean switching between standard DMDA
 and the custom methods associated with Q2 finite elements.
*/
PetscErrorCode DMDAUseQ2Mappings(DM dm,PetscBool useQ2)
{
  DMDAQ2Mapping  *mapping = NULL;
  PetscContainer myg2lctx = NULL;
  DM             cdm;
  PetscErrorCode ierr;
  
  ierr = DMGetCoordinateDM(dm,&cdm);CHKERRQ(ierr);
  
  // query for object
  ierr = PetscObjectQuery((PetscObject)dm,"DMDACustomMapping",(PetscObject*)&myg2lctx);CHKERRQ(ierr);
  if (!myg2lctx) {
    PetscContainer c = NULL;
    PetscInt e,nelements,nen_u,i;
    const PetscInt *elnidx_u;
    PTable table;
    PetscInt               num_gindices,ge_eqnums[3*Q2_NODES_PER_EL_3D];
    PetscInt               vel_el_lidx[3*U_BASIS_FUNCTIONS];
    ISLocalToGlobalMapping ltog;
    const PetscInt         *gindices;
    
    
    // If the container is not found,
    // (i) initialize DMDAQ2Mapping mapping
    // (ii) copy original function pointers, store them
    // (iii) create containers
    // (iv) stuff mapping in to the container and attach it to the dm
    
    // (i)
    ierr = PetscMalloc1(1,&mapping);CHKERRQ(ierr);
    ierr = PetscMemzero(mapping,sizeof(DMDAQ2Mapping));CHKERRQ(ierr);
    mapping->refcnt = 1;
    
    // (iii)
    ierr = PetscContainerCreate(PetscObjectComm((PetscObject)dm),&c);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(c,(void*)mapping);CHKERRQ(ierr);
    ierr = PetscContainerSetUserDestroy(c,PetscContainer_DMDAQ2MappingDestroy);CHKERRQ(ierr);
    
    // (iv) - compose early so DMCreateLocalVector_DAQ2() can be called
    ierr = PetscObjectCompose((PetscObject)dm,"DMDACustomMapping",(PetscObject)c);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)cdm,"DMDACustomMapping",(PetscObject)c);CHKERRQ(ierr);
    
    // setup sizes
    ierr = DMDAGetLocalSizeElementQ2(dm,&mapping->mx,&mapping->my,&mapping->mz);CHKERRQ(ierr);
    mapping->nelements = mapping->mx * mapping->my * mapping->mz;
    mapping->nlocal = (2*mapping->mx+1) * (2*mapping->my+1) * (2*mapping->mz+1);
    
    // setup element arrays
    ierr = PetscMalloc1(mapping->nelements*27,&mapping->elements_default);CHKERRQ(ierr);
    ierr = PetscMalloc1(mapping->nelements*27,&mapping->elements_custom);CHKERRQ(ierr);
    
    ierr = DMDAGetElements_pTatinQ2P1(dm,&nelements,&nen_u,&elnidx_u);CHKERRQ(ierr);
    
    ierr = PetscMemcpy(mapping->elements_default,(PetscInt*)elnidx_u,sizeof(PetscInt)*mapping->nelements*27);CHKERRQ(ierr);
    ierr = DMDAQ2ElementLayout(mapping->elements_custom,mapping->mx,mapping->my,mapping->mz);CHKERRQ(ierr);
    
    // setup table
    ierr = PTableCreate(PETSC_COMM_SELF,&table);CHKERRQ(ierr);
    ierr = PTableSetRange(table,0,mapping->nlocal * 3);CHKERRQ(ierr);
    ierr = PTableSetType(table,PTABLE_DENSE);CHKERRQ(ierr);
    ierr = PTableSetup(table);CHKERRQ(ierr);
    
    ierr = DMGetLocalToGlobalMapping(dm,&ltog);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetSize(ltog,&num_gindices);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetIndices(ltog,&gindices);CHKERRQ(ierr);
    
    for (e=0; e<nelements; e++) {
      
      /* get local indices */
      ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
      
      /* get global indices */
      for (i=0; i<nen_u; i++) {
        const int nid = elnidx_u[nen_u*e + i]; /* use default map to fetch global dof indices */
        
        ge_eqnums[3*i  ] = gindices[ 3*nid   ];
        ge_eqnums[3*i+1] = gindices[ 3*nid+1 ];
        ge_eqnums[3*i+2] = gindices[ 3*nid+2 ];
      }
      
      for (i=0; i<nen_u; i++) {
        PetscInt basis_id,local_dof_id;
        
        basis_id     = mapping->elements_custom[nen_u*e + i]; /* use custom map to fetch global dof indices for row */
        
        local_dof_id = 3 * basis_id + 0;
        ierr = PTableSetValue(table,local_dof_id,ge_eqnums[3*i  ]);CHKERRQ(ierr);
        
        local_dof_id = 3 * basis_id + 1;
        ierr = PTableSetValue(table,local_dof_id,ge_eqnums[3*i+1]);CHKERRQ(ierr);
        
        local_dof_id = 3 * basis_id + 2;
        ierr = PTableSetValue(table,local_dof_id,ge_eqnums[3*i+2]);CHKERRQ(ierr);
      }
    }
    ierr = ISLocalToGlobalMappingRestoreIndices(ltog,&gindices);CHKERRQ(ierr);
    ierr = PTableSynchronize(table);CHKERRQ(ierr);
    
    // setup IS
    ierr = PTableFlattenIntoIS(table,PETSC_TRUE,&mapping->to,&mapping->from);CHKERRQ(ierr);
    
    // setup VecScatter
    {
      Vec xg,xl;
      
      ierr = DMCreateGlobalVector(dm,&xg);CHKERRQ(ierr);
      ierr = DMCreateLocalVector_DAQ2(dm,&xl);CHKERRQ(ierr);
      
      ierr = VecScatterCreate(xg,mapping->from,xl,mapping->to,&mapping->global2local);CHKERRQ(ierr);
      
      ierr = VecDestroy(&xl);CHKERRQ(ierr);
      ierr = VecDestroy(&xg);CHKERRQ(ierr);
    }
    
    // (ii)
    mapping->createlocalvector_default  = dm->ops->createlocalvector;
    mapping->creatematrix_default       = dm->ops->creatematrix;
    mapping->globaltolocalbegin_default = dm->ops->globaltolocalbegin;
    mapping->globaltolocalend_default   = dm->ops->globaltolocalend;
    mapping->localtoglobalbegin_default = dm->ops->localtoglobalbegin;
    mapping->localtoglobalend_default   = dm->ops->localtoglobalend;
    
    mapping->createlocalvector_default_c  = cdm->ops->createlocalvector;
    mapping->creatematrix_default_c       = cdm->ops->creatematrix;
    mapping->globaltolocalbegin_default_c = cdm->ops->globaltolocalbegin;
    mapping->globaltolocalend_default_c   = cdm->ops->globaltolocalend;
    mapping->localtoglobalbegin_default_c = cdm->ops->localtoglobalbegin;
    mapping->localtoglobalend_default_c   = cdm->ops->localtoglobalend;
    
    
    ierr = PTableDestroy(&table);CHKERRQ(ierr);
    
    myg2lctx = c;
    mapping = NULL;
  }
  
  ierr = PetscContainerGetPointer(myg2lctx,(void**)&mapping);CHKERRQ(ierr);
  if (!useQ2) {
    DM_DA *da_context = (DM_DA*)dm->data;
    
    // flush exisiting cached vectors //
    ierr = DMClearGlobalVectors(dm);CHKERRQ(ierr);
    ierr = DMClearLocalVectors(dm);CHKERRQ(ierr);
    
    // restore original function pointers
    dm->ops->createlocalvector  = mapping->createlocalvector_default;
    dm->ops->creatematrix       = mapping->creatematrix_default;
    dm->ops->globaltolocalbegin = mapping->globaltolocalbegin_default;
    dm->ops->globaltolocalend   = mapping->globaltolocalend_default;
    dm->ops->localtoglobalbegin = mapping->localtoglobalbegin_default;
    dm->ops->localtoglobalend   = mapping->localtoglobalend_default;
    
    // over-ride elements
    ierr = PetscMemcpy(da_context->e,mapping->elements_default,sizeof(PetscInt)*mapping->nelements*27);CHKERRQ(ierr);
    
    ierr = DMGetCoordinateDM(dm,&cdm);CHKERRQ(ierr);
    // restore original function pointers
    cdm->ops->createlocalvector  = mapping->createlocalvector_default_c;
    cdm->ops->creatematrix       = mapping->creatematrix_default_c;
    cdm->ops->globaltolocalbegin = mapping->globaltolocalbegin_default_c;
    cdm->ops->globaltolocalend   = mapping->globaltolocalend_default_c;
    cdm->ops->localtoglobalbegin = mapping->localtoglobalbegin_default_c;
    cdm->ops->localtoglobalend   = mapping->localtoglobalend_default_c;
    
    ierr = VecDestroy(&dm->coordinatesLocal);CHKERRQ(ierr);
    dm->coordinatesLocal = NULL;
  } else {
    DM_DA *da_context = (DM_DA*)dm->data;
    
    PetscPrintf(PETSC_COMM_WORLD,"**** using custom Q2 DMDA mappings ****\n");
    
    // flush exisiting cached vectors //
    ierr = DMClearGlobalVectors(dm);CHKERRQ(ierr);
    ierr = DMClearLocalVectors(dm);CHKERRQ(ierr);
    
    // swap in new methods //
    dm->ops->createlocalvector  = DMCreateLocalVector_DAQ2;
    //da->ops->creatematrix       = DMCreateMatrix_DAQ2;
    dm->ops->globaltolocalbegin = DMGlobalToLocalBegin_DAQ2;
    dm->ops->globaltolocalend   = DMGlobalToLocalEnd_DAQ2;
    dm->ops->localtoglobalbegin = DMLocalToGlobalBegin_DAQ2;
    dm->ops->localtoglobalend   = DMLocalToGlobalEnd_DAQ2;
    
    // over-ride elements
    ierr = PetscMemcpy(da_context->e,mapping->elements_custom,sizeof(PetscInt)*mapping->nelements*27);CHKERRQ(ierr);
    
    ierr = DMGetCoordinateDM(dm,&cdm);CHKERRQ(ierr);
    // swap in new methods - under the assumption that dm.dof = cdm.dof //
    cdm->ops->createlocalvector  = DMCreateLocalVector_DAQ2;
    cdm->ops->globaltolocalbegin = DMGlobalToLocalBegin_DAQ2;
    cdm->ops->globaltolocalend   = DMGlobalToLocalEnd_DAQ2;
    cdm->ops->localtoglobalbegin = DMLocalToGlobalBegin_DAQ2;
    cdm->ops->localtoglobalend   = DMLocalToGlobalEnd_DAQ2;
    
    ierr = VecDestroy(&dm->coordinatesLocal);CHKERRQ(ierr);
    dm->coordinatesLocal = NULL;
  }
  PetscFunctionReturn(0);
}

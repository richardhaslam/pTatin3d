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
 **    filename:   dmdae.c
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

#include "dmdae.h"


PetscErrorCode DMDAECreate(DMDAE *dae)
{
  DMDAE d;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *dae = NULL;
  ierr = PetscMalloc(sizeof(struct _p_DMDAE),&d);CHKERRQ(ierr);
  ierr = PetscMemzero(d,sizeof(struct _p_DMDAE));CHKERRQ(ierr);
  *dae = d;

  PetscFunctionReturn(0);
}

PetscErrorCode DMDAEDestroy(DMDAE *dae)
{
  DMDAE d;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!dae){ PetscFunctionReturn(0); }
  d = (*dae);
  ierr = PetscFree(d->lsip);CHKERRQ(ierr);
  ierr = PetscFree(d->lsjp);CHKERRQ(ierr);
  ierr = PetscFree(d->lskp);CHKERRQ(ierr);
  ierr = PetscFree(d->lmxp);CHKERRQ(ierr);
  ierr = PetscFree(d->lmyp);CHKERRQ(ierr);
  ierr = PetscFree(d->lmzp);CHKERRQ(ierr);
  ierr = PetscFree(d);CHKERRQ(ierr);
  *dae = NULL;

  PetscFunctionReturn(0);
}

PetscErrorCode DMDAEDeepCopy(DMDAE dae1,PetscInt NP[],DMDAE dae2)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dae2->ne  = dae1->ne;
  dae2->lne = dae1->lne;

  dae2->mx = dae1->mx;
  dae2->my = dae1->my;
  dae2->mz = dae1->mz;

  dae2->lmx = dae1->lmx;
  dae2->lmy = dae1->lmy;
  dae2->lmz = dae1->lmz;

  dae2->si = dae1->si;
  dae2->sj = dae1->sj;
  dae2->sk = dae1->sk;

  dae2->npe = dae1->npe;
  dae2->nps = dae1->nps;
  dae2->overlap = dae1->overlap;

  dae2->sgi = dae1->sgi;
  dae2->sgj = dae1->sgj;
  dae2->sgk = dae1->sgk;

  ierr = PetscMalloc(sizeof(PetscInt)*NP[0],&dae2->lsip);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscInt)*NP[1],&dae2->lsjp);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscInt)*NP[2],&dae2->lskp);CHKERRQ(ierr);
  ierr = PetscMemcpy(dae2->lsip,dae1->lsip,sizeof(PetscInt)*NP[0]);CHKERRQ(ierr);
  ierr = PetscMemcpy(dae2->lsjp,dae1->lsjp,sizeof(PetscInt)*NP[1]);CHKERRQ(ierr);
  ierr = PetscMemcpy(dae2->lskp,dae1->lskp,sizeof(PetscInt)*NP[2]);CHKERRQ(ierr);

  ierr = PetscMalloc(sizeof(PetscInt)*NP[0],&dae2->lmxp);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscInt)*NP[1],&dae2->lmyp);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscInt)*NP[2],&dae2->lmzp);CHKERRQ(ierr);
  ierr = PetscMemcpy(dae2->lsip,dae1->lmxp,sizeof(PetscInt)*NP[0]);CHKERRQ(ierr);
  ierr = PetscMemcpy(dae2->lsjp,dae1->lmyp,sizeof(PetscInt)*NP[1]);CHKERRQ(ierr);
  ierr = PetscMemcpy(dae2->lskp,dae1->lmzp,sizeof(PetscInt)*NP[2]);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMDAECopy(DMDAE dae1,DMDAE dae2)
{

  PetscFunctionBegin;
  dae2->ne  = dae1->ne;
  dae2->lne = dae1->lne;

  dae2->mx = dae1->mx;
  dae2->my = dae1->my;
  dae2->mz = dae1->mz;

  dae2->lmx = dae1->lmx;
  dae2->lmy = dae1->lmy;
  dae2->lmz = dae1->lmz;

  dae2->si = dae1->si;
  dae2->sj = dae1->sj;
  dae2->sk = dae1->sk;

  dae2->npe = dae1->npe;
  dae2->nps = dae1->nps;
  dae2->overlap = dae1->overlap;

  dae2->sgi = dae1->sgi;
  dae2->sgj = dae1->sgj;
  dae2->sgk = dae1->sgk;

  dae2->lsip = dae1->lsip;
  dae2->lsjp = dae1->lsjp;
  dae2->lskp = dae1->lskp;

  dae2->lmxp = dae1->lmxp;
  dae2->lmyp = dae1->lmyp;
  dae2->lmzp = dae1->lmzp;

  PetscFunctionReturn(0);
}

PetscErrorCode DMAttachDMDAE(DM dm)
{
  PetscContainer container;
  DMDAE dae;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMDAECreate(&dae);CHKERRQ(ierr);

  ierr = PetscContainerCreate(PETSC_COMM_WORLD,&container);CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(container,(void*)dae);CHKERRQ(ierr);

  ierr = PetscObjectCompose((PetscObject)dm,"DMDAEobject",(PetscObject)container);CHKERRQ(ierr);
  ierr = PetscContainerDestroy(&container);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMGetDMDAE(DM dm,DMDAE *dae)
{
  DMDAE          d;
  PetscContainer container;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *dae = NULL;
  ierr = PetscObjectQuery((PetscObject)dm,"DMDAEobject",(PetscObject*)&container);CHKERRQ(ierr);
  if (!container) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"No data with name \"DMDAEobject\" was composed with DM");
  ierr = PetscContainerGetPointer(container,(void**)&d);CHKERRQ(ierr);
  *dae = d;

  PetscFunctionReturn(0);
}

PetscErrorCode DMDestroyDMDAE(DM dm)
{
  DMDAE d;
  PetscContainer container;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)dm,"DMDAEobject",(PetscObject*)&container);CHKERRQ(ierr);
  if (!container) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"No data with name \"DMDAEobject\" was composed with DM");
  ierr = PetscContainerGetPointer(container,(void**)&d);CHKERRQ(ierr);

  ierr = DMDAEDestroy(&d);CHKERRQ(ierr);
  //ierr = PetscContainerDestroy(&container);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMDAEGetOwnershipRanges(DM dm,
    PetscInt *m,PetscInt *n,PetscInt *p,
    PetscInt **si,PetscInt **sj,PetscInt **sk,
    PetscInt **mx,PetscInt **my,PetscInt **mz)
{
  DMDAE dae;
  PetscInt pM,pN,pP;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetDMDAE(dm,&dae);CHKERRQ(ierr);

  ierr = DMDAGetInfo(dm,0,0,0,0,&pM,&pN,&pP, 0, 0, 0,0,0, 0 );CHKERRQ(ierr);
  if (m) { *m = pM; }
  if (n) { *n = pN; }
  if (p) { *p = pP; }

  if (si) { *si = dae->lsip; }
  if (sj) { *sj = dae->lsjp; }
  if (sk) { *sk = dae->lskp; }

  if (mx) { *mx = dae->lmxp; }
  if (my) { *my = dae->lmyp; }
  if (mz) { *mz = dae->lmzp; }

  PetscFunctionReturn(0);
}

PetscErrorCode DMDAEGetCornersElement(DM dm,PetscInt *esi,PetscInt *esj,PetscInt *esk,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
  DMDAE dae;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetDMDAE(dm,&dae);CHKERRQ(ierr);
  if (esi) { *esi = dae->si; }
  if (esj) { *esj = dae->sj; }
  if (esk) { *esk = dae->sk; }

  if (mx) { *mx = dae->lmx; }
  if (my) { *my = dae->lmy; }
  if (mz) { *mz = dae->lmz; }

  PetscFunctionReturn(0);
}


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
 **    Filename:      dmdae.h
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

#ifndef __ptatin_DMDAE_h__
#define __ptatin_DMDAE_h__

#include "petsc.h"
#include "petscdm.h"

/* little data bucket for all things needed to define elements */
typedef struct _p_DMDAE *DMDAE;

struct _p_DMDAE {
	PetscInt ne,lne;
	PetscInt mx,my,mz;
	PetscInt lmx,lmy,lmz;
	PetscInt *lsip,*lsjp,*lskp;
	PetscInt *lmxp,*lmyp,*lmzp;
	PetscInt si,sj,sk;
	PetscInt npe; /* nodes per element */
	PetscInt nps; /* nodes per side */
	PetscInt overlap;
	PetscInt sgi,sgj,sgk;
};

PetscErrorCode DMDAEDeepCopy(DMDAE dae1,PetscInt NP[],DMDAE dae2);
PetscErrorCode DMDAECopy(DMDAE dae1,DMDAE dae2);
PetscErrorCode DMDAECreate(DMDAE *dae);
PetscErrorCode DMDAEDestroy(DMDAE *dae);
PetscErrorCode DMAttachDMDAE(DM dm);
PetscErrorCode DMGetDMDAE(DM dm,DMDAE *dae);
PetscErrorCode DMDestroyDMDAE(DM dm);

/* helpers */
PetscErrorCode DMDAEGetOwnershipRanges(DM dm,
																			 PetscInt *m,PetscInt *n,PetscInt *p,
																			 PetscInt **si,PetscInt **sj,PetscInt **sk,
																			 PetscInt **mx,PetscInt **my,PetscInt **mz);
PetscErrorCode DMDAEGetCornersElement(DM dm,PetscInt *esi,PetscInt *esj,PetscInt *esk,PetscInt *mx,PetscInt *my,PetscInt *mz);

#endif


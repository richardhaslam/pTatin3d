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
 **    filename:   dmda_project_coords.c
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
#include <petsc/private/dmdaimpl.h>

#include "dmda_update_coords.h"
#include "dmda_project_coords.h"

PetscErrorCode DMDARestrictCoordinatesHierarchy(DM da[],PetscInt nlevels)
{
  DM daf,dac;
  DM cdaf,cdac;
  Vec coordsc,coordsf;
  Mat inject;
  PetscInt L;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  for (L=nlevels-1; L>0; L--) {
    daf = da[L];
    dac = da[L-1];

    ierr = DMGetCoordinateDM(dac,&cdac);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(daf,&cdaf);CHKERRQ(ierr);

    ierr = DMGetCoordinates(dac,&coordsc);CHKERRQ(ierr);
    ierr = DMGetCoordinates(daf,&coordsf);CHKERRQ(ierr);

    ierr = DMCreateInjection(cdac,cdaf,&inject);CHKERRQ(ierr);
    ierr = MatMult(inject,coordsf,coordsc);CHKERRQ(ierr);
    ierr = MatDestroy(&inject);CHKERRQ(ierr);
  }

  /* update ghost coordinates for all levels except the finest */
  for (L=0; L<nlevels-1; L++) {
    ierr = DMDAUpdateGhostedCoordinates(da[L]);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DMDARestrictCoordinates(DM daf,DM dac)
{
  DM cdaf,cdac;
  Vec coordsc,coordsf;
  Mat inject;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = DMGetCoordinateDM(dac,&cdac);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(daf,&cdaf);CHKERRQ(ierr);

  ierr = DMGetCoordinates(dac,&coordsc);CHKERRQ(ierr);
  ierr = DMGetCoordinates(daf,&coordsf);CHKERRQ(ierr);

  ierr = DMCreateInjection(cdac,cdaf,&inject);CHKERRQ(ierr);
  ierr = MatMult(inject,coordsf,coordsc);CHKERRQ(ierr);
  ierr = MatDestroy(&inject);CHKERRQ(ierr);

  /* update ghost coordinates */
  ierr = DMDAUpdateGhostedCoordinates(dac);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode  DMDASetCoarseningFactor(DM da,PetscInt cx,PetscInt cy,PetscInt cz)
{
  DM_DA *dd = (DM_DA*)da->data;
  if (cx > 0) dd->coarsen_x = cx;
  if (cy > 0) dd->coarsen_y = cy;
  if (cz > 0) dd->coarsen_z = cz;

  PetscFunctionReturn(0);
}


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
 **    filename:   dmda_view_petscvtk.c
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

#include "dmda_view_petscvtk.h"


#undef __FUNCT__
#define __FUNCT__ "DMDAViewPetscVTK"
PetscErrorCode DMDAViewPetscVTK(DM da,Vec field,const char name[])
{
	Vec x;
	PetscViewer vv;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)da), name, &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);

	/* view mesh */
	ierr = DMView(da, vv);CHKERRQ(ierr);

	/* view field */
	/* if the vector is valid, write it out - else create an empty field */
	if (field) {
		const char *name;
		name = NULL;
		/* temp work around - calling GetName forces a name to be inserted if you isn't there 
		- in parallel an error will occur if [1]PETSC ERROR: VecView_MPI_DA() line 464 in src/dm/impls/da/gr2.c
		if the name is null
		*/
		ierr = PetscObjectGetName( (PetscObject)field,&name);CHKERRQ(ierr);
		ierr = VecView(field,vv);CHKERRQ(ierr);
	} else {
		ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)x, "empty_field" );CHKERRQ(ierr);
		ierr = VecView(x,vv);CHKERRQ(ierr);
		ierr  = VecDestroy(&x);CHKERRQ(ierr);
	}
	
	ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

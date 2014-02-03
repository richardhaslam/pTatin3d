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
 **    Filename:      dmda_project_coords.c
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


#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_update_coords.h"
#include "dmda_project_coords.h"

#undef __FUNCT__
#define __FUNCT__ "DMDARestrictCoordinatesHierarchy"
PetscErrorCode DMDARestrictCoordinatesHierarchy(DM da[],PetscInt nlevels)
{
	DM daf,dac;
	DM cdaf,cdac;
	Vec coordsc,coordsf;
	VecScatter inject;
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
		ierr = VecScatterBegin(inject,coordsf,coordsc,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
		ierr = VecScatterEnd(inject  ,coordsf,coordsc,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
		ierr = VecScatterDestroy(&inject);CHKERRQ(ierr);
	}
	
	/* update ghost coordinates for all levels except the finest */
	for (L=0; L<nlevels-1; L++) {
		ierr = DMDAUpdateGhostedCoordinates(da[L]);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMDARestrictCoordinates"
PetscErrorCode DMDARestrictCoordinates(DM daf,DM dac)
{
	DM cdaf,cdac;
	Vec coordsc,coordsf;
	VecScatter inject;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMGetCoordinateDM(dac,&cdac);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(daf,&cdaf);CHKERRQ(ierr);
	
	ierr = DMGetCoordinates(dac,&coordsc);CHKERRQ(ierr);
	ierr = DMGetCoordinates(daf,&coordsf);CHKERRQ(ierr);
	
	ierr = DMCreateInjection(cdac,cdaf,&inject);CHKERRQ(ierr);
	ierr = VecScatterBegin(inject,coordsf,coordsc,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterEnd(inject  ,coordsf,coordsc,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterDestroy(&inject);CHKERRQ(ierr);
	
	/* update ghost coordinates */
	ierr = DMDAUpdateGhostedCoordinates(dac);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

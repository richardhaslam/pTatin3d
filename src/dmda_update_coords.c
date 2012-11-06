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
 **    Filename:      dmda_update_coords.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_update_coords.h"

#undef __FUNCT__
#define __FUNCT__ "DMDAUpdateGhostedCoordinates"
PetscErrorCode DMDAUpdateGhostedCoordinates(DM da)
{
	PetscErrorCode ierr;
	Vec da_coordinates, gcoords;
	DM _dac;
	
	PetscFunctionBegin;
	ierr = DMDAGetCoordinateDA(da,&_dac);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(da,&da_coordinates);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(_dac,da_coordinates,INSERT_VALUES,gcoords);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(_dac,da_coordinates,INSERT_VALUES,gcoords);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetCoordinatesFromLocalVector"
PetscErrorCode DMDASetCoordinatesFromLocalVector(DM da,Vec local_coords)
{
	PetscErrorCode ierr;
	Vec da_coordinates;
	DM dac;
	
	PetscFunctionBegin;
	ierr = DMDAGetCoordinateDA(da,&dac);CHKERRQ(ierr);
	
	/* scatter new existing coords into global_coords */
	ierr = DMDAGetCoordinates( da, &da_coordinates );CHKERRQ(ierr);
	ierr = VecZeroEntries(da_coordinates);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin( dac, local_coords, INSERT_VALUES, da_coordinates );CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  ( dac, local_coords, INSERT_VALUES, da_coordinates );CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetCoordinatesU"
PetscErrorCode DMDASetCoordinatesU(DM da,Vec coords)
{
	PetscErrorCode ierr;
	Vec da_coords;
	
	PetscFunctionBegin;
	ierr = DMDAGetCoordinates(da,&da_coords);CHKERRQ(ierr);
	ierr = VecCopy(coords,da_coords);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDACloneCoordinates"
PetscErrorCode DMDACloneCoordinates(DM da,DM da_clone)
{
	PetscErrorCode ierr;
	Vec coords, coords_clone;
	
	PetscFunctionBegin;
	ierr = DMDAGetCoordinates( da, &coords ); CHKERRQ(ierr);
	ierr = DMDAGetCoordinates( da_clone, &coords_clone ); CHKERRQ(ierr);
	ierr = VecCopy( coords, coords_clone ); CHKERRQ(ierr);

	ierr = DMDAGetGhostedCoordinates( da, &coords ); CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( da_clone, &coords_clone ); CHKERRQ(ierr);
	ierr = VecCopy( coords, coords_clone ); CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

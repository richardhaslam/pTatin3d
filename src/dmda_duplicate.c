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
 **    filename:   dmda_duplicate.c
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

#include "dmda_update_coords.h"
#include "dmda_duplicate.h"

#undef __FUNCT__
#define __FUNCT__ "DMDADuplicateLayout"
PetscErrorCode DMDADuplicateLayout(DM da1,PetscInt dof2,PetscInt sw2,DMDAStencilType st2,DM *da2)
{
    PetscInt        si1,sj1,sk1;
    PetscInt        mx1,my1,mz1;
    PetscInt        M1,N1,P1;
    PetscInt        cx1,cy1,cz1;
    PetscInt        dim1;
    PetscInt        dof1,sw1;
    DMDAStencilType st1;
    DMBoundaryType  wrap1[3];
    const PetscInt  *lx,*ly,*lz;
    Vec             coords;
    PetscErrorCode  ierr;

    PetscFunctionBegin;
    ierr = DMDAGetInfo( da1, &dim1, &M1,&N1,&P1, &cx1,&cy1,&cz1, &dof1, &sw1, &wrap1[0],&wrap1[1],&wrap1[2], &st1 );CHKERRQ(ierr);
    ierr = DMDAGetCorners( da1, &si1,&sj1,&sk1 , &mx1,&my1,&mz1 );CHKERRQ(ierr);

    if (dof2==PETSC_DECIDE) { dof2 = dof1; }
    if (sw2==PETSC_DECIDE)  {  sw2 = sw1;  }

    if (dim1==1) {
        ierr = DMDAGetOwnershipRanges(da1,&lx,NULL,NULL);CHKERRQ(ierr);
        ierr = DMDACreate1d(PetscObjectComm((PetscObject)da1), wrap1[0],M1, dof2,sw2, lx, da2);CHKERRQ(ierr);
    } else if (dim1==2) {
        ierr = DMDAGetOwnershipRanges(da1,&lx,&ly,NULL);CHKERRQ(ierr);
        ierr = DMDACreate2d(PetscObjectComm((PetscObject)da1), wrap1[0],wrap1[1],st2, M1,N1,cx1,cy1, dof2,sw2, lx,ly, da2);CHKERRQ(ierr);
    } else if (dim1==3) {
        ierr = DMDAGetOwnershipRanges(da1,&lx,&ly,&lz);CHKERRQ(ierr);
        ierr = DMDACreate3d(PetscObjectComm((PetscObject)da1), wrap1[0],wrap1[1],wrap1[2],st2, M1,N1,P1,cx1,cy1,cz1, dof2,sw2, lx,ly,lz, da2);CHKERRQ(ierr);
    } else {
        SETERRQ(PetscObjectComm((PetscObject)da1),PETSC_ERR_USER,"Unknown dimension for DMDA");
    }
    ierr = DMSetUp(*da2);CHKERRQ(ierr);

    ierr = DMGetCoordinates(da1,&coords);CHKERRQ(ierr);
    if (coords) {
        if (dim1==1) {
            ierr = DMDASetUniformCoordinates(*da2,-1.0,1.0, 0.,0., 0.,0.);CHKERRQ(ierr);
        } else if (dim1==2) {
            ierr = DMDASetUniformCoordinates(*da2,-1.0,1.0, -1.0,1.0, 0.,0.);CHKERRQ(ierr);
        } else {
            ierr = DMDASetUniformCoordinates(*da2,-1.0,1.0, -1.0,1.0, -1.0,1.0);CHKERRQ(ierr);
        }
        ierr = DMDASetCoordinatesU(*da2,coords);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


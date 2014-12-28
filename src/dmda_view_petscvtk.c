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
 **    Filename:      dmda_view_petscvtk.c
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

#include "dmda_view_petscvtk.h"


#undef __FUNCT__
#define __FUNCT__ "DMDAViewPetscLegacyVTK"
PetscErrorCode DMDAViewPetscLegacyVTK(DM da,Vec field,const char filename[])
{
    PetscViewer    vv;
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)da),filename,&vv);CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(vv,PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);

    /* view mesh */
    ierr = DMView(da,vv);CHKERRQ(ierr);

    /* view field */
    /* if the vector is valid, write it out - else create an empty field */
    if (field) {
        const char *name;
		
        name = NULL;
        /* temp work around - calling GetName forces a name to be inserted if you isn't there
         - in parallel an error will occur if [1]PETSC ERROR: VecView_MPI_DA() line 464 in src/dm/impls/da/gr2.c
         if the name is null
        */
        ierr = PetscObjectGetName((PetscObject)field,&name);CHKERRQ(ierr);
        ierr = VecView(field,vv);CHKERRQ(ierr);
    } else {
        Vec x;
        
        ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)x,"empty_field");CHKERRQ(ierr);
        ierr = VecView(x,vv);CHKERRQ(ierr);
        ierr = VecDestroy(&x);CHKERRQ(ierr);
    }
    ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAViewPetscVTS"
PetscErrorCode DMDAViewPetscVTS(DM dm,Vec field,const char filename[])
{
    MPI_Comm       comm;
    PetscViewer    viewer;
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    PetscObjectGetComm((PetscObject)dm,&comm);
    ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERVTK);CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);
    
    if (field) {
        ierr = VecView(field,viewer);CHKERRQ(ierr);
    } else {
        Vec x;
        
        ierr = DMCreateGlobalVector(dm,&x);CHKERRQ(ierr);
        ierr = VecView(x,viewer);CHKERRQ(ierr);
        ierr = VecDestroy(&x);CHKERRQ(ierr);
    }

    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAViewFieldsPetscVTS"
PetscErrorCode DMDAViewFieldsPetscVTS(DM dm,PetscInt nf,Vec fields[],const char filename[])
{
    MPI_Comm       comm;
    PetscInt       n;
    PetscViewer    viewer;
    Vec            dummy;
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    
    if (nf == 1) {
        ierr = DMDAViewPetscVTS(dm,fields[0],filename);CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
    
    PetscObjectGetComm((PetscObject)dm,&comm);
    ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERVTK);CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr);

    /* create a dummy vector if required */
    dummy = NULL;
    for (n=0; n<nf; n++) {
        if (!fields[n]) {
            ierr = DMCreateGlobalVector(dm,&dummy);CHKERRQ(ierr);
            break;
        }
    }
    
    /* sweep through all fields[], checking for null entries */
    for (n=0; n<nf; n++) {
        Vec field = fields[n];
        
        if (field) {
            ierr = VecView(field,viewer);CHKERRQ(ierr);
        } else {
            ierr = VecView(dummy,viewer);CHKERRQ(ierr);
        }
    }
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    if (dummy) {
        ierr = VecDestroy(&dummy);CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

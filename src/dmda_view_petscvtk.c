
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
	ierr = PetscViewerASCIIOpen(((PetscObject)(da))->comm, name, &vv);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);

	/* view mesh */
	ierr = DMView(da, vv);CHKERRQ(ierr);

	/* view field */
	/* if the vector is valid, write it out - else create an empty field */
	if (field) {
		const char *name;
		name = PETSC_NULL;
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

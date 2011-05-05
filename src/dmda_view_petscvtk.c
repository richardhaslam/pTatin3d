
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
		ierr = VecView(field,vv);CHKERRQ(ierr);
	} else {
		ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)x, "empty_field" );CHKERRQ(ierr);
		ierr = VecView(x,vv);CHKERRQ(ierr);
		ierr  = VecDestroy(x);CHKERRQ(ierr);
	}
	
	ierr = PetscViewerDestroy(vv);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
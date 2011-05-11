

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
		
		ierr = DMDAGetCoordinateDA(dac,&cdac);CHKERRQ(ierr);
		ierr = DMDAGetCoordinateDA(daf,&cdaf);CHKERRQ(ierr);
		
		ierr = DMDAGetCoordinates(dac,&coordsc);CHKERRQ(ierr);
		ierr = DMDAGetCoordinates(daf,&coordsf);CHKERRQ(ierr);
		
		ierr = DMGetInjection(cdac,cdaf,&inject);CHKERRQ(ierr);
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
	
	ierr = DMDAGetCoordinateDA(dac,&cdac);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(daf,&cdaf);CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinates(dac,&coordsc);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(daf,&coordsf);CHKERRQ(ierr);
	
	ierr = DMGetInjection(cdac,cdaf,&inject);CHKERRQ(ierr);
	ierr = VecScatterBegin(inject,coordsf,coordsc,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterEnd(inject  ,coordsf,coordsc,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterDestroy(&inject);CHKERRQ(ierr);
	
	/* update ghost coordinates */
	ierr = DMDAUpdateGhostedCoordinates(dac);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

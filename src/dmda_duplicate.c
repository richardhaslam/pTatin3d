
#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "dmda_update_coords.h"
#include "dmda_duplicate.h"

#undef __FUNCT__
#define __FUNCT__ "DMDADuplicateLayout"
PetscErrorCode DMDADuplicateLayout(DM da1,PetscInt dof2,PetscInt sw2,DMDAStencilType st2,DM *da2)
{
	PetscInt si1,sj1,sk1;
	PetscInt mx1,my1,mz1;
	PetscInt M1,N1,P1;
	PetscInt cx1,cy1,cz1;
	PetscInt dim1;
	PetscInt dof1,sw1;
	DMDAStencilType st1;
	DMDAPeriodicType wrap1;
	const PetscInt *lx,*ly,*lz;
	Vec coords;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo( da1, &dim1, &M1,&N1,&P1, &cx1,&cy1,&cz1, &dof1, &sw1, &wrap1, &st1 );CHKERRQ(ierr);
	ierr = DMDAGetCorners( da1, &si1,&sj1,&sk1 , &mx1,&my1,&mz1 );CHKERRQ(ierr);

	if (dof2==PETSC_DECIDE) { dof2 = dof1; }
	if (sw2==PETSC_DECIDE)  {  sw2 = sw1;  }
	if (st2==PETSC_DECIDE)  {  st2 = st1;  }

	if (dim1==1) {
		ierr = DMDAGetOwnershipRanges(da1,&lx,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		ierr = DMDACreate1d(((PetscObject)da1)->comm, wrap1,M1, dof2,sw2, lx, da2);CHKERRQ(ierr);
	} else if (dim1==2) {
		ierr = DMDAGetOwnershipRanges(da1,&lx,&ly,PETSC_NULL);CHKERRQ(ierr);
		ierr = DMDACreate2d(((PetscObject)da1)->comm, wrap1,st2, M1,N1,cx1,cy1, dof2,sw2, lx,ly, da2);CHKERRQ(ierr);
	} else if (dim1==3) {
		ierr = DMDAGetOwnershipRanges(da1,&lx,&ly,&lz);CHKERRQ(ierr);
		ierr = DMDACreate3d(((PetscObject)da1)->comm, wrap1,st2, M1,N1,P1,cx1,cy1,cz1, dof2,sw2, lx,ly,lz, da2);CHKERRQ(ierr);
	} else {
		SETERRQ(((PetscObject)da1)->comm,PETSC_ERR_USER,"Unknown dimension for DMDA");
	}

	ierr = DMDAGetCoordinates(da1,&coords);CHKERRQ(ierr);
	if (coords) {
		if (dim1==1) {
			ierr = DMDASetUniformCoordinates(*da2,-1.0,1.0, PETSC_NULL,PETSC_NULL, PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		} else if (dim1==2) {
			ierr = DMDASetUniformCoordinates(*da2,-1.0,1.0, -1.0,1.0, PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		} else {
			ierr = DMDASetUniformCoordinates(*da2,-1.0,1.0, -1.0,1.0, -1.0,1.0);CHKERRQ(ierr);
		}
		ierr = DMDASetCoordinatesU(*da2,coords);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}



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
/*@
   DMDAUpdateGhostedCoordinates - Performs the scatter from the DA's global coodinate vector
      to the local coordinate vector. The local coordinate vector includes the ghost value.

   Collective on DA

   Input Parameter:
.  da - the distributed array

  Level: intermediate

.keywords: distributed array, get, corners, nodes, local indices, coordinates

.seealso: DAGetGhostCorners(), DASetCoordinates(), DASetUniformCoordinates(), DAGetCoordinates(), DAGetCoordinateDA()
@*/
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
/*@
   DMDASetCoordinatesFromLocalVector - Sets into the DA's global coordinate vector
      the coordinate value stored in a local vector (ghost nodes are ignored).
      The local vector has the same parallel layout as the one obtained from
      calling DAGetGhostedCoordinates().

   Not Collective

   Input Parameter:
+  da - the distributed array
-  local_coords - local coordinate vector

   Note:
    The coordinates may include those for all ghost points, however ghost values 
    will be ignored when the local coords are inserted into the global coordinate vector.

     The user is responsible for destroying the vector local_coords.

  Level: intermediate

.keywords: distributed array, get, corners, nodes, local indices, coordinates

.seealso: DAGetGhostCorners(), DAGetCoordinates(), DAGetGhostCoordinates(), DASetUniformCoordinates(), DAGetCoordinateDA()
@*/
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
/*@
   DMDASetCoordinatesU - Sets into the DA a vector that indicates the 
      coordinates of the local nodes (NOT including ghost nodes).
      Following the setting of the new coordinates, the ghost values 
      stored on the local vector are also updated.

   Not Collective

   Input Parameter:
+  da - the distributed array
-  c - coordinate vector

   Note:
    The coordinates should NOT include those for all ghost points

     Does NOT increase the reference count of this vector, so caller should NOT
  destroy the vector.

  Level: intermediate

.keywords: distributed array, get, corners, nodes, local indices, coordinates

.seealso: DAGetGhostCorners(), DAGetCoordinates(), DASetUniformCoordinates(). DAGetGhostCoordinates(), DAGetCoordinateDA()
@*/
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

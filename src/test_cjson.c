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
 **    Filename:      test_cjson.c
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

#include "petsc.h"
#include "ptatin3d.h"
#include "ptatin_init.h"
#include "geometry_object.h"
#include "geometry_object_parse.h"


#undef __FUNCT__
#define __FUNCT__ "test_GeometryObjectParse_cJSON"
PetscErrorCode test_GeometryObjectParse_cJSON(void)
{
	PetscErrorCode ierr;
    PetscInt k,ngo;
    GeometryObject *golist;
    
    PetscPrintf(PETSC_COMM_WORLD,"** test_GeometryObjectParse_cJSON **\n");
    ierr = GeometryObjectLoadJSON("src/test_option_files/test_geom.json",&ngo,&golist);CHKERRQ(ierr);

    for (k=0; k<ngo; k++) {
        ierr = GeometryObjectView(golist[k]);CHKERRQ(ierr);
    }
    
    for (k=0; k<ngo; k++) {
        ierr = GeometryObjectDestroy(&golist[k]);CHKERRQ(ierr);
    }
    free(golist);
    
	PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	
	ierr = pTatinInitialize(&argc,&argv,(char *)0,NULL);CHKERRQ(ierr);

    //for (int i=0; i<100000000; i++) {
        ierr = test_GeometryObjectParse_cJSON();CHKERRQ(ierr);
	//}
    
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

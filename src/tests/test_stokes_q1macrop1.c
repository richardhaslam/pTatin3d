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
 **    filename:   test_stokes_q1macrop1.c
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

static const char help[] = "Test mesh generation of Q1MacroP1\n\n";


#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"

#include "ptatin_utils.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1macrop1.h"


extern PetscErrorCode PhysCompCreateMesh_Stokes3d_Q1MacroP1(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx);



PetscErrorCode test_q1macrop1_a(void)
{
	PetscErrorCode ierr;
	PhysCompStokes stokes;
	PetscInt mx,my,mz,m,n,p;
	PetscInt nels,*map,e;
	DM dav;
	PetscMPIInt rank;
	PetscFunctionBegin;

	mx = my = mz = 4;
	PetscOptionsGetInt(NULL,NULL,"-mx",&mx,0);
	PetscOptionsGetInt(NULL,NULL,"-my",&my,0);
	PetscOptionsGetInt(NULL,NULL,"-mz",&mz,0);

	ierr = PhysCompCreate_Stokes(&stokes);CHKERRQ(ierr);
	ierr = PhysCompCreateMesh_Stokes3d_Q1MacroP1(mx,my,mz,stokes);CHKERRQ(ierr);

	dav = stokes->dav;
	MPI_Comm_rank(PetscObjectComm((PetscObject)dav),&rank);

	ierr = DMDAEQ1Macro_MixedSpace_GetSizeElement(dav,&m,&n,&p);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"MixedSpace size: %D x %D x %D \n", m,n,p);

	ierr = DMDAEQ1Macro_MixedSpace_GetLocalSizeElement(dav,&m,&n,&p);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"[%d] MixedSpace local size: %D x %D x %D \n", (int)rank,m,n,p);

	ierr = DMDAEQ1Macro_NaturalSpace_GetSizeElement(dav,&m,&n,&p);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"NaturalSpace size: %D x %D x %D \n", m,n,p);

	ierr = DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(dav,&m,&n,&p);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"[%d] NaturalSpace local size: %D x %D x %D \n", (int)rank,m,n,p);


	ierr = DMDAEGetElementMap_Q1MacroNaturalToMixedSpace(dav,&nels,&map);CHKERRQ(ierr);
	printf("Map(Natural->Mixed): \n");
	for (e=0; e<nels; e++) {
//		printf(" [e=%d]: %d \n", e,map[e] );
	}

	ierr = DMDAEGetElementMap_Q1MacroNaturalToMixedLocalSpace(stokes->dav,&nels,&map);CHKERRQ(ierr);
	printf("Map(Natural->MixedLocal): \n");
	for (e=0; e<nels; e++) {
//		printf(" [e=%d]: %d \n",e, map[e] );
	}



	PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
	PetscErrorCode ierr;

	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);

	ierr = test_q1macrop1_a();CHKERRQ(ierr);

	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

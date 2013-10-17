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
 **    Filename:      ptatin_init.c
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
#include "stdio.h"
#include "stdlib.h"

#define PTATIN_SVN_REVISION "r4831 2013-10-13 19:54:50 +0200 dmay"

#define STR_VALUE(arg)      #arg
#define STRINGIFY_ARG(name) STR_VALUE(name)

#undef __FUNCT__
#define __FUNCT__ "pTatinCheckCompilationFlags"
PetscErrorCode pTatinCheckCompilationFlags(const char flags[])
{
	char *loc = NULL;
	int throw_warning = 0;
	
	PetscFunctionBegin;
	
	loc = strstr(flags,"-O0");
	if (loc != NULL) {
		throw_warning = 1;
	}
	if (throw_warning == 1) {
		PetscPrintf(PETSC_COMM_WORLD,"**                                                                       \n");
		PetscPrintf(PETSC_COMM_WORLD,"** pTatin3d is configured with debug options \n");
		PetscPrintf(PETSC_COMM_WORLD,"**   TATIN_CFLAGS = %s\n",flags);
		PetscPrintf(PETSC_COMM_WORLD,"** For signifcant performance improvements, please consult the makefile  \n");
		PetscPrintf(PETSC_COMM_WORLD,"** and set TATIN_CFLAGS to an optimized value suitable for your machine. \n");
		PetscPrintf(PETSC_COMM_WORLD,"**                                                                       \n");
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinWritePreamble"
PetscErrorCode pTatinWritePreamble(void)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"** ====================================================================================== \n");
	PetscPrintf(PETSC_COMM_WORLD,"**\n");	
	PetscPrintf(PETSC_COMM_WORLD,"**             ___________                          _______\n");
	PetscPrintf(PETSC_COMM_WORLD,"**     _______/          /_____ ________ __ _   ___/       \\ ____\n");
	PetscPrintf(PETSC_COMM_WORLD,"**    /      /___    ___/      /__   __/  /  \\ /  /\\__ _   /     \\\n");
	PetscPrintf(PETSC_COMM_WORLD,"**   /  //  /   /   /  /  //  /  /  / /  /    /  /___/_   /  //  /\n");
	PetscPrintf(PETSC_COMM_WORLD,"**  /  ___ /   /   /  /  _   /  /  / /  /  /    //       /  //  /\n");
	PetscPrintf(PETSC_COMM_WORLD,"** /__/       /___/  /__//__/  /__/ /__/__/ \\__//_______/______/\n");
	PetscPrintf(PETSC_COMM_WORLD,"**\n");	
	PetscPrintf(PETSC_COMM_WORLD,"** Authors:  Dave A. May          (dave.may@erdw.ethz.ch)           \n");
	PetscPrintf(PETSC_COMM_WORLD,"**           Laetitia Le Pourhiet (laetitia.le_pourhiet@upmc.fr)    \n");
	PetscPrintf(PETSC_COMM_WORLD,"** Revision: %s #\n", PTATIN_SVN_REVISION);
	PetscPrintf(PETSC_COMM_WORLD,"**\n");	
	
#ifdef COMPFLAGS
	#define STR_ARG_NAME STRINGIFY_ARG(COMPFLAGS)
	ierr = pTatinCheckCompilationFlags(STR_ARG_NAME);CHKERRQ(ierr);
#endif

	PetscPrintf(PETSC_COMM_WORLD,"** ====================================================================================== \n");
				 
	PetscFunctionReturn(0);
}


extern PetscErrorCode KSPCreate_ChebychevRN(KSP ksp);
extern PetscErrorCode PCCreate_SemiRedundant(PC pc);

#undef __FUNCT__
#define __FUNCT__ "pTatinInitialize"
PetscErrorCode pTatinInitialize(int *argc,char ***args,const char file[],const char help[])
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscInitialize(argc,args,file,help);CHKERRQ(ierr);

	ierr = KSPRegisterDynamic("chebychevrn","./","KSPCreate_ChebychevRN",KSPCreate_ChebychevRN);CHKERRQ(ierr);
	ierr = PCRegisterDynamic("semiredundant","./","PCCreate_SemiRedundant",PCCreate_SemiRedundant);CHKERRQ(ierr);

	ierr = pTatinWritePreamble();CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinFinalize"
PetscErrorCode pTatinFinalize(void)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

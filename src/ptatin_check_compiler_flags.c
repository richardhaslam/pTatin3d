
#include "petsc.h"
#include "stdio.h"
#include "stdlib.h"

#define PTATIN_SVN_REVISION "3636 2012-11-06 19:14:24Z dmay"

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
	pTatinCheckCompilationFlags(STR_ARG_NAME);
#endif

	PetscPrintf(PETSC_COMM_WORLD,"** ====================================================================================== \n");
				 
	PetscFunctionReturn(0);
}



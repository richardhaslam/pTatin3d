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
 **    filename:   ptatin_init.c
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


#include "petsc.h"
#include "stdio.h"
#include "stdlib.h"

#include "ptatin_version_info.h"
#include "ptatin3d.h"
#include "ptatin_models.h"

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
		PetscPrintf(PETSC_COMM_WORLD,"** WARNING pTatin3d appears to have been compiled with debug options \n");
		//PetscPrintf(PETSC_COMM_WORLD,"**   TATIN_CFLAGS = %s\n",flags);
		PetscPrintf(PETSC_COMM_WORLD,"** For significant performance improvements, please consult the file makefile.arch  \n");
		PetscPrintf(PETSC_COMM_WORLD,"** Adjust TATIN_CFLAGS to include aggressive compiler optimizations \n");
		PetscPrintf(PETSC_COMM_WORLD,"**                                                                       \n");
	}

#if defined(__AVX__) && !defined(TATIN_HAVE_AVX)
		PetscPrintf(PETSC_COMM_WORLD,"** WARNING pTatin3d appears to have been compiled without AVX support, on a system with AVX support. \n");
		PetscPrintf(PETSC_COMM_WORLD,"** For significant performance gains, set CONFIG_AVX=y in makefile.arch\n");
		PetscPrintf(PETSC_COMM_WORLD,"**                                                                     \n");
#endif
	
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
	PetscPrintf(PETSC_COMM_WORLD,"**           Jed Brown            (jedbrown@mcs.anl.gov)            \n");
	PetscPrintf(PETSC_COMM_WORLD,"**\n");
	
#if(PTATIN_VERSION_INFO == 1)
	PetscPrintf(PETSC_COMM_WORLD,"** %s \n", PTATIN_VERSION_CNTR_REPO);
	PetscPrintf(PETSC_COMM_WORLD,"** %s \n", PTATIN_VERSION_CNTR_REVISION);
	PetscPrintf(PETSC_COMM_WORLD,"** %s \n", PTATIN_VERSION_CNTR_LOG);
#endif
    
#ifdef COMPFLAGS
	#define STR_ARG_NAME STRINGIFY_ARG(COMPFLAGS)
	PetscPrintf(PETSC_COMM_WORLD,"**                                                                       \n");
	PetscPrintf(PETSC_COMM_WORLD,"** TATIN_CFLAGS = %s\n",STR_ARG_NAME);
	PetscPrintf(PETSC_COMM_WORLD,"**                                                                       \n");
	ierr = pTatinCheckCompilationFlags(STR_ARG_NAME);CHKERRQ(ierr);
	#undef STR_ARG_NAME
#endif
	PetscPrintf(PETSC_COMM_WORLD,"** ====================================================================================== \n");
				 
	PetscFunctionReturn(0);
}


extern PetscErrorCode KSPCreate_ChebychevRN(KSP ksp);
extern PetscErrorCode PCCreate_SemiRedundant(PC pc);
extern PetscErrorCode PCCreate_WSMP(PC pc);
extern PetscErrorCode PCCreate_DMDARepart(PC pc);
extern PetscLogEvent MAT_MultMFA11;
extern PetscLogEvent MAT_MultMFA11_setup;
extern PetscLogEvent MAT_MultMFA11_copyto;
extern PetscLogEvent MAT_MultMFA11_kernel;
extern PetscLogEvent MAT_MultMFA11_copyfrom;
extern PetscLogEvent MAT_MultMFA11_merge;

extern PetscLogEvent MAT_MultMFA; /* stokes operator */
extern PetscLogEvent MAT_MultMFA12; /* stokes operator - gradient operator */
extern PetscLogEvent MAT_MultMFA21; /* stokes operator - divergence operator */

extern PetscLogEvent MAT_MultMFA_QuasiNewtonX; 
extern PetscLogEvent MAT_MultMFA11_QuasiNewtonX;
extern PetscLogEvent MAT_MultMFA12_QuasiNewtonX;
extern PetscLogEvent MAT_MultMFA21_QuasiNewtonX;

PetscClassId PTATIN_CLASSID;
extern PetscLogEvent PTATIN_DataExchangerTopologySetup;
extern PetscLogEvent PTATIN_DataExchangerBegin;
extern PetscLogEvent PTATIN_DataExchangerEnd;

extern PetscLogEvent PTATIN_MaterialPointAdvGlobalCoordUpdate;
extern PetscLogEvent PTATIN_MaterialPointAdvLocalCoordUpdate;
extern PetscLogEvent PTATIN_MaterialPointAdvCommunication;
extern PetscLogEvent PTATIN_MaterialPointAdvRemoval;

//extern PetscLogEvent PTATIN_MaterialPointProjection;

extern PetscLogEvent PTATIN_MaterialPointPopulationControlInsert;
extern PetscLogEvent PTATIN_MaterialPointPopulationControlRemove;

extern PetscLogEvent PTATIN_ModelInitialize;
extern PetscLogEvent PTATIN_ModelApplyInitialSolution;
extern PetscLogEvent PTATIN_ModelApplyInitialMeshGeometry;
extern PetscLogEvent PTATIN_ModelApplyInitialMaterialGeometry;
extern PetscLogEvent PTATIN_ModelApplyInitialStokesVariableMarkers;
extern PetscLogEvent PTATIN_ModelApplyBoundaryCondition;
extern PetscLogEvent PTATIN_ModelApplyBoundaryConditionMG;
extern PetscLogEvent PTATIN_ModelApplyMaterialBoundaryCondition;
extern PetscLogEvent PTATIN_ModelUpdateMeshGeometry;
extern PetscLogEvent PTATIN_ModelOutput;


#undef __FUNCT__
#define __FUNCT__ "pTatinInitialize"
PetscErrorCode pTatinInitialize(int *argc,char ***args,const char file[],const char help[])
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = PetscInitialize(argc,args,file,help);CHKERRQ(ierr);

  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);

	ierr = KSPRegister("chebychevrn",KSPCreate_ChebychevRN);CHKERRQ(ierr);
	ierr = PCRegister("semiredundant",PCCreate_SemiRedundant);CHKERRQ(ierr);
	ierr = PCRegister("wsmp",PCCreate_WSMP);CHKERRQ(ierr);
	ierr = PCRegister("dmdarepart",PCCreate_DMDARepart);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA11",MAT_CLASSID,&MAT_MultMFA11);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA11_setup",MAT_CLASSID,&MAT_MultMFA11_setup);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA11_copyto",MAT_CLASSID,&MAT_MultMFA11_copyto);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA11_kernel",MAT_CLASSID,&MAT_MultMFA11_kernel);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA11_copy_from",MAT_CLASSID,&MAT_MultMFA11_copyfrom);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA11_merge",MAT_CLASSID,&MAT_MultMFA11_merge);CHKERRQ(ierr);

	ierr = PetscLogEventRegister("MatMultMFA",  MAT_CLASSID,&MAT_MultMFA  );CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA12",MAT_CLASSID,&MAT_MultMFA12);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA21",MAT_CLASSID,&MAT_MultMFA21);CHKERRQ(ierr);

	ierr = PetscLogEventRegister("MatMultMFA_QuasiNewtonX",  MAT_CLASSID,&MAT_MultMFA_QuasiNewtonX  );CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA11_QuasiNewtonX",MAT_CLASSID,&MAT_MultMFA11_QuasiNewtonX);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA12_QuasiNewtonX",MAT_CLASSID,&MAT_MultMFA12_QuasiNewtonX);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MatMultMFA21_QuasiNewtonX",MAT_CLASSID,&MAT_MultMFA21_QuasiNewtonX);CHKERRQ(ierr);

	ierr = PetscClassIdRegister("ptatin",&PTATIN_CLASSID);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("DataExTopoSetup",PTATIN_CLASSID,&PTATIN_DataExchangerTopologySetup);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("DataExBegin",    PTATIN_CLASSID,&PTATIN_DataExchangerBegin);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("DataExEnd",      PTATIN_CLASSID,&PTATIN_DataExchangerEnd);CHKERRQ(ierr);

	ierr = PetscLogEventRegister("ModelInit",      PTATIN_CLASSID,&PTATIN_ModelInitialize);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelInitSoln",  PTATIN_CLASSID,&PTATIN_ModelApplyInitialSolution);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelInitMesh",  PTATIN_CLASSID,&PTATIN_ModelApplyInitialMeshGeometry);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelInitMat",   PTATIN_CLASSID,&PTATIN_ModelApplyInitialMaterialGeometry);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelInitStkVar",PTATIN_CLASSID,&PTATIN_ModelApplyInitialStokesVariableMarkers);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelBC",        PTATIN_CLASSID,&PTATIN_ModelApplyBoundaryCondition);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelBCMG",      PTATIN_CLASSID,&PTATIN_ModelApplyBoundaryConditionMG);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelMatBC",     PTATIN_CLASSID,&PTATIN_ModelApplyMaterialBoundaryCondition);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelUpdateMesh",PTATIN_CLASSID,&PTATIN_ModelUpdateMeshGeometry);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("ModelOutput",    PTATIN_CLASSID,&PTATIN_ModelOutput);CHKERRQ(ierr);
	
	ierr = PetscLogEventRegister("MPAdvGCoord", PTATIN_CLASSID,&PTATIN_MaterialPointAdvGlobalCoordUpdate);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MPAdvLCoord", PTATIN_CLASSID,&PTATIN_MaterialPointAdvLocalCoordUpdate);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MPAdvComm",   PTATIN_CLASSID,&PTATIN_MaterialPointAdvCommunication);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MPAdvRemove", PTATIN_CLASSID,&PTATIN_MaterialPointAdvRemoval);CHKERRQ(ierr);
	
	ierr = PetscLogEventRegister("MPPCInsert", PTATIN_CLASSID,&PTATIN_MaterialPointPopulationControlInsert);CHKERRQ(ierr);
	ierr = PetscLogEventRegister("MPPCRemove", PTATIN_CLASSID,&PTATIN_MaterialPointPopulationControlRemove);CHKERRQ(ierr);
	
	ierr = pTatinWritePreamble();CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinFinalize"
PetscErrorCode pTatinFinalize(void)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
  ierr = pTatinModelDeRegisterAll();CHKERRQ(ierr);
	ierr = PetscFinalize();CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#include "../ptatin_models.h"

/* Users add prototypes here */
extern PetscErrorCode pTatinModelRegister_Template(void);
extern PetscErrorCode pTatinModelRegister_ViscousSinker(void);
extern PetscErrorCode pTatinModelRegister_Gene3D(void);
extern PetscErrorCode pTatinModelRegister_Gene3DNueve(void);
extern PetscErrorCode pTatinModelRegister_Indentor(void);
extern PetscErrorCode pTatinModelRegister_Rift3D(void);
extern PetscErrorCode pTatinModelRegister_Rift3D_T(void);
extern PetscErrorCode pTatinModelRegister_Sierra(void);
extern PetscErrorCode pTatinModelRegister_Folding(void);
extern PetscErrorCode pTatinModelRegister_Folding2d(void);
extern PetscErrorCode pTatinModelRegister_AdvDiffExample(void);
extern PetscErrorCode pTatinModelRegister_BasinComp(void);
extern PetscErrorCode pTatinModelRegister_FaultFold(void);
extern PetscErrorCode pTatinModelRegister_WrenchFold(void);
extern PetscErrorCode pTatinModelRegister_Delamination(void);
extern PetscErrorCode pTatinModelRegister_Riftrh(void);
extern PetscErrorCode pTatinModelRegister_GeoMod2008(void);
extern PetscErrorCode pTatinModelRegister_FaultFoldPlastic(void);
extern PetscErrorCode pTatinModelRegister_MultilayerFolding(void);
extern PetscErrorCode pTatinModelRegister_SubmarineLavaFlow(void);
extern PetscErrorCode pTatinModelRegister_ExSubduction(void);
extern PetscErrorCode pTatinModelRegister_iPLUS(void);
extern PetscErrorCode pTatinModelRegister_Subduction_Initiation2d(void);
extern PetscErrorCode pTatinModelRegister_Thermal_Convection2d(void);

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegisterAll"
PetscErrorCode pTatinModelRegisterAll(void)
{
	PetscErrorCode ierr;

	PetscFunctionBegin;
	/* call registration functions for each model here */
	ierr = pTatinModelRegister_Template();CHKERRQ(ierr);
	ierr = pTatinModelRegister_ViscousSinker();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Gene3D();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Gene3DNueve();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Indentor();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Rift3D();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Rift3D_T();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Sierra();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Folding();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Folding2d();CHKERRQ(ierr);
	ierr = pTatinModelRegister_AdvDiffExample();CHKERRQ(ierr);
	ierr = pTatinModelRegister_BasinComp();CHKERRQ(ierr);
	ierr = pTatinModelRegister_FaultFold();CHKERRQ(ierr);
	ierr = pTatinModelRegister_WrenchFold();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Delamination();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Riftrh();CHKERRQ(ierr);
	ierr = pTatinModelRegister_GeoMod2008();CHKERRQ(ierr);
	ierr = pTatinModelRegister_FaultFoldPlastic();CHKERRQ(ierr);
	ierr = pTatinModelRegister_MultilayerFolding();CHKERRQ(ierr);
	ierr = pTatinModelRegister_SubmarineLavaFlow();CHKERRQ(ierr);
	ierr = pTatinModelRegister_ExSubduction();CHKERRQ(ierr);
	ierr = pTatinModelRegister_iPLUS();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Subduction_Initiation2d();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Thermal_Convection2d();CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

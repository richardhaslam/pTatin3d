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
 **    filename:   ptatin_models_reg.c
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
#include "../ptatin_models.h"

/* Users add prototypes here */
extern PetscErrorCode pTatinModelRegister_Template(void);
extern PetscErrorCode pTatinModelRegister_ViscousSinker(void);
extern PetscErrorCode pTatinModelRegister_Gene3D(void);
extern PetscErrorCode pTatinModelRegister_Indentor(void);
extern PetscErrorCode pTatinModelRegister_Rift3D_T(void);
extern PetscErrorCode pTatinModelRegister_AdvDiffExample(void);
extern PetscErrorCode pTatinModelRegister_Delamination(void);
extern PetscErrorCode pTatinModelRegister_Riftrh(void);
extern PetscErrorCode pTatinModelRegister_Rift_oblique3d(void);
extern PetscErrorCode pTatinModelRegister_GeoMod2008(void);
extern PetscErrorCode pTatinModelRegister_MultilayerFolding(void);
extern PetscErrorCode pTatinModelRegister_SubmarineLavaFlow(void);
extern PetscErrorCode pTatinModelRegister_ExSubduction(void);
extern PetscErrorCode pTatinModelRegister_iPLUS(void);
extern PetscErrorCode pTatinModelRegister_Subduction_Initiation2d(void);
extern PetscErrorCode pTatinModelRegister_Thermal_Convection2d(void);
extern PetscErrorCode pTatinModelRegister_ThermalSB(void);
extern PetscErrorCode pTatinModelRegister_SD3D(void);
extern PetscErrorCode pTatinModelRegister_PAS(void);
extern PetscErrorCode pTatinModelRegister_PD(void);
extern PetscErrorCode pTatinModelCreate_StaticBox(pTatinModel);
extern PetscErrorCode pTatinModelCreate_StaticBoxTM(pTatinModel);
extern PetscErrorCode pTatinModelCreate_AnlVV(pTatinModel m);

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
	ierr = pTatinModelRegister_Indentor();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Rift3D_T();CHKERRQ(ierr);
	ierr = pTatinModelRegister_AdvDiffExample();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Delamination();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Riftrh();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Rift_oblique3d();CHKERRQ(ierr);
	ierr = pTatinModelRegister_GeoMod2008();CHKERRQ(ierr);
	ierr = pTatinModelRegister_MultilayerFolding();CHKERRQ(ierr);
	ierr = pTatinModelRegister_SubmarineLavaFlow();CHKERRQ(ierr);
	ierr = pTatinModelRegister_ExSubduction();CHKERRQ(ierr);
	ierr = pTatinModelRegister_iPLUS();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Subduction_Initiation2d();CHKERRQ(ierr);
	ierr = pTatinModelRegister_Thermal_Convection2d();CHKERRQ(ierr);
	ierr = pTatinModelRegister_ThermalSB();CHKERRQ(ierr);
	ierr = pTatinModelRegister_SD3D();CHKERRQ(ierr);
	ierr = pTatinModelRegister_PAS();CHKERRQ(ierr);
	ierr = pTatinModelRegister_PD();CHKERRQ(ierr);

  ierr = pTatinModelDynamicRegister("static_box",pTatinModelCreate_StaticBox);CHKERRQ(ierr);
  ierr = pTatinModelDynamicRegister("static_box_thermomech",pTatinModelCreate_StaticBoxTM);CHKERRQ(ierr);
  ierr = pTatinModelDynamicRegister("analytics_vv",pTatinModelCreate_AnlVV);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

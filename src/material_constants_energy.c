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
 **    filename:   material_constants_energy.c
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
#include "data_bucket.h"
#include "material_constants_energy.h"

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsEnergyInitialize"
PetscErrorCode MaterialConstantsEnergyInitialize(DataBucket db)
{
	DataBucketRegisterField(db,EnergyMaterialConstants_classname,       sizeof(EnergyMaterialConstants),NULL);
	DataBucketRegisterField(db,EnergySourceConst_classname,             sizeof(EnergySourceConst),NULL);
	DataBucketRegisterField(db,EnergySourceDecay_classname,             sizeof(EnergySourceDecay),NULL);
	DataBucketRegisterField(db,EnergySourceAdiabaticAdvection_classname,sizeof(EnergySourceAdiabaticAdvection),NULL);
	DataBucketFinalize(db);
	
  PetscFunctionReturn(0);
}



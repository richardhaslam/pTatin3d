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
 **    filename:   material_constants_energy.h
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

#ifndef __material_constants_energy_h__
#define __material_constants_energy_h__

/* add auto geneterated classes here */
#include "material_constants/EnergyMaterialConstants_def.h"
#include "material_constants/EnergySourceConst_def.h"
#include "material_constants/EnergySourceDecay_def.h"
#include "material_constants/EnergySourceAdiabaticAdvection_def.h"


typedef enum { 
	ENERGYSOURCE_NONE = 0,
	ENERGYSOURCE_CONSTANT ,
    ENERGYSOURCE_SHEAR_HEATING,
    ENERGYSOURCE_DECAY,
    ENERGYSOURCE_ADIABATIC,
    ENERGYSOURCE_ADIABATIC_ADVECTION
} EnergySourceType;

typedef enum {
    ENERGYCONDUCTIVITY_CONSTANT = 0,
    ENERGYCONDUCTIVITY_TEMP_DEP_THRESHOLD,
} EnergyConductivityType;


typedef enum {
	ENERGYDENSITY_CONSTANT = 0,
	ENERGYDENSITY_BOUSSINESQ
} EnergyDensityType;

PetscErrorCode MaterialConstantsEnergyInitialize(DataBucket db);

/*
PetscErrorCode MaterialConstantsSetDefaults(DataBucket db);

PetscErrorCode MaterialConstantsSetFromOptions_MaterialType(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_MaterialType(DataBucket db,const int region_id,PetscInt visc_t,PetscInt plast_t,PetscInt soft_t,PetscInt dens_t);
PetscErrorCode MaterialConstantsPrintValues_MaterialType(DataBucket db,const int region_id);

PetscErrorCode MaterialConstantsSetFromOptions_ViscosityConst(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_ViscosityConst(DataBucket db,const int region_id,PetscReal viscosity);
PetscErrorCode MaterialConstantsScaleValues_ViscosityConst(DataBucket db,const int region_id,PetscReal eta_star);
PetscErrorCode MaterialConstantsPrintValues_ViscosityConst(DataBucket db,const int region_id);
*/

#endif

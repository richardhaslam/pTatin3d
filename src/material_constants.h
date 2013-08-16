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
 **    Filename:      material_constants.h
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

#ifndef __material_constants_h__
#define __material_constants_h__

/* add auto geneterated classes here */
#include "material_constants/MaterialConst_MaterialType_def.h"
#include "material_constants/MaterialConst_ViscosityConst_def.h"
#include "material_constants/MaterialConst_DensityConst_def.h"
#include "material_constants/MaterialConst_DensityBoussinesq_def.h"
#include "material_constants/MaterialConst_ViscosityZ_def.h"
#include "material_constants/MaterialConst_ViscosityFK_def.h"
#include "material_constants/MaterialConst_ViscosityArrh_def.h"
#include "material_constants/MaterialConst_PlasticMises_def.h"
#include "material_constants/MaterialConst_PlasticDP_def.h"
#include "material_constants/MaterialConst_SoftLin_def.h"
#include "material_constants/MaterialConst_SoftExpo_def.h"

typedef enum { 
	VISCOUS_CONSTANT=0,
	VISCOUS_FRANKK,
    	VISCOUS_Z,
	VISCOUS_ARRHENIUS,
	VISCOUS_ARRHENIUS_2
} ViscousType;

typedef enum { 
	PLASTIC_NONE=0,
	PLASTIC_MISES,
	PLASTIC_DP
} PlasticType;

typedef enum { 
	SOFTENING_NONE=0,
	SOFTENING_LINEAR,
	SOFTENING_EXPONENTIAL
} SofteningType;

typedef enum { 
	DENSITY_CONSTANT=0,
	DENSITY_BOUSSINESQ
} DensityType;

PetscErrorCode MaterialConstantsInitialize(DataBucket *_db);
PetscErrorCode MaterialConstantsSetDefaults(DataBucket db);

PetscErrorCode MaterialConstantsSetFromOptions_MaterialType(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_MaterialType(DataBucket db,const int region_id,PetscInt visc_t,PetscInt plast_t,PetscInt soft_t,PetscInt dens_t);
PetscErrorCode MaterialConstantsPrintValues_MaterialType(DataBucket db,const int region_id);

PetscErrorCode MaterialConstantsSetFromOptions_ViscosityConst(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_ViscosityConst(DataBucket db,const int region_id,PetscReal viscosity);
PetscErrorCode MaterialConstantsScaleValues_ViscosityConst(DataBucket db,const int region_id,PetscReal eta_star);
PetscErrorCode MaterialConstantsPrintValues_ViscosityConst(DataBucket db,const int region_id);


PetscErrorCode MaterialConstantsSetFromOptions_ViscosityFK(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_ViscosityFK(DataBucket db,const int region_id,PetscReal eta0, PetscReal theta);
PetscErrorCode MaterialConstantsScaleValues_ViscosityFK(DataBucket db,const int region_id,PetscReal eta_star);
PetscErrorCode MaterialConstantsPrintValues_ViscosityFK(DataBucket db,const int region_id);

/* Viscosity Arrh has no scaling function */ 
PetscErrorCode MaterialConstantsSetFromOptions_ViscosityArrh(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_ViscosityArrh(DataBucket db,const int region_id,PetscReal preexpA,PetscReal Ascale,PetscReal entalpy,PetscReal Vmol,PetscReal nexp,PetscReal Tref);
PetscErrorCode MaterialConstantsPrintValues_ViscosityArrh(DataBucket db,const int region_id);


PetscErrorCode MaterialConstantsSetFromOptions_DensityConst(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_DensityConst(DataBucket db,const int region_id,PetscReal density);
PetscErrorCode MaterialConstantsScaleValues_DensityConst(DataBucket db,const int region_id,PetscReal rho_star);
PetscErrorCode MaterialConstantsPrintValues_DensityConst(DataBucket db,const int region_id);


PetscErrorCode MaterialConstantsSetFromOptions_DensityBoussinesq(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_DensityBoussinesq(DataBucket db,const int region_id,PetscReal density,PetscReal alpha,PetscReal beta);
PetscErrorCode MaterialConstantsPrintValues_DensityBoussinesq(DataBucket db,const int region_id);
PetscErrorCode MaterialConstantsScaleValues_DensityBoussinesq(DataBucket db,const int region_id,PetscReal rho_star,PetscReal sigma_star);

PetscErrorCode MaterialConstantsSetFromOptions_PlasticMises(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_PlasticMises(DataBucket db,const int region_id,PetscReal yield_stress,PetscReal yield_stress_inf);
PetscErrorCode MaterialConstantsScaleValues_PlasticMises(DataBucket db,const int region_id,PetscReal stress_star);
PetscErrorCode MaterialConstantsPrintValues_PlasticMises(DataBucket db,const int region_id);

PetscErrorCode MaterialConstantsSetFromOptions_PlasticDP(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_PlasticDP(DataBucket db,const int region_id,PetscReal friction,PetscReal friction_inf,PetscReal cohesion,PetscReal cohesion_inf,PetscReal tens_cutoff,PetscReal hst_cutoff);
PetscErrorCode MaterialConstantsScaleValues_PlasticDP(DataBucket db,const int region_id, PetscReal stress_star);
PetscErrorCode MaterialConstantsPrintValues_PlasticDP(DataBucket db,const int region_id);

PetscErrorCode MaterialConstantsSetFromOptions_ViscosityZ(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetValues_ViscosityZ(DataBucket db,const int region_id,PetscReal eta0,PetscReal zeta,PetscReal zref);
PetscErrorCode MaterialConstantsScaleValues_ViscosityZ(DataBucket db,const int region_id,PetscReal eta_star,PetscReal L_star);
PetscErrorCode MaterialConstantsPrintValues_ViscosityZ(DataBucket db,const int region_id);

PetscErrorCode MaterialConstantsSetDefault_SoftLin(DataBucket db);
PetscErrorCode MaterialConstantsSetFromOptions_SoftLin(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsPrintValues_SoftLin(DataBucket db,const int region_id);
PetscErrorCode MaterialConstantsSetValues_SoftLin(DataBucket db,const int region_id,PetscReal emin,PetscReal emax);



PetscErrorCode MaterialConstantsSetDefault_SoftExpo(DataBucket db);
PetscErrorCode MaterialConstantsSetFromOptions_SoftExpo(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsPrintValues_SoftExpo(DataBucket db,const int region_id);
PetscErrorCode MaterialConstantsSetValues_SoftExpo(DataBucket db,const int region_id,PetscReal emin,PetscReal efold);


PetscErrorCode MaterialConstantsReportParseError(const char model_name[],const char field_name[],const int region);

PetscErrorCode MaterialConstantsScaleAll(DataBucket db,const int region_id,PetscReal L_star, PetscReal U_star,PetscReal t_star,PetscReal eta_star,PetscReal rho_star,PetscReal P_star);
PetscErrorCode MaterialConstantsSetFromOptions(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsPrintAll(DataBucket db,const int region_id);
#endif


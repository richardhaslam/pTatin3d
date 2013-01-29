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
 **    Filename:      material_constants.c
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


#define _GNU_SOURCE
#include "petsc.h"
#include "swarm_fields.h"
#include "rheology.h"
#include "material_constants.h"


#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsInitialize"
PetscErrorCode MaterialConstantsInitialize(DataBucket *_db)
{
	DataBucket db;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	DataBucketCreate(&db);
	
	DataBucketRegisterField(db,MaterialConst_MaterialType_classname,      sizeof(MaterialConst_MaterialType),PETSC_NULL);
	DataBucketRegisterField(db,MaterialConst_ViscosityConst_classname,    sizeof(MaterialConst_ViscosityConst),PETSC_NULL);
    DataBucketRegisterField(db,MaterialConst_ViscosityZ_classname,    sizeof(MaterialConst_ViscosityZ),PETSC_NULL);
	DataBucketRegisterField(db,MaterialConst_PlasticMises_classname,      sizeof(MaterialConst_PlasticMises),PETSC_NULL);
    DataBucketRegisterField(db,MaterialConst_PlasticDP_classname,      sizeof(MaterialConst_PlasticDP),PETSC_NULL);

	DataBucketFinalize(db);
	

	DataBucketSetInitialSizes(db,200,0);
	ierr = MaterialConstantsSetDefaults(db);CHKERRQ(ierr);
	
	
	*_db = db;
	
	PetscFunctionReturn(0);
}

/* MaterialType */
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetDefault_MaterialType"
PetscErrorCode MaterialConstantsSetDefault_MaterialType(DataBucket db)
{
	int                        r,nregions;
	DataField                  PField;
	MaterialConst_MaterialType *data;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(db,&nregions,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(db,MaterialConst_MaterialType_classname,&PField);

	data = (MaterialConst_MaterialType*)PField->data; /* should write a function to do this */
	for (r=0; r<nregions; r++) {
		data[r].visc_type      = VISCOUS_CONSTANT;
		data[r].plastic_type   = PLASTIC_NONE;
		data[r].softening_type = SOFTENING_NONE;
		data[r].density_type   = DENSITY_CONSTANT;
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetFromOptions_MaterialType"
PetscErrorCode MaterialConstantsSetFromOptions_MaterialType(DataBucket db,const char model_name[],const int region_id,PetscBool essential)
{
	char                         opt_name[256];
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	PetscBool                    found;
	PetscInt                     value;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_MaterialType_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* visc_type */
	value = 1.0; /* default - not required is nothing happens if option not found */
	sprintf(opt_name,"-%s_%d",MaterialConst_MaterialType_member_names[0],region_id);
	ierr = PetscOptionsGetInt(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_MaterialTypeSetField_visc_type(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_MaterialType_member_names[0],region_id);CHKERRQ(ierr);
	}

    /* plastic_type */
	value = 1.0; /* default - not required is nothing happens if option not found */
	sprintf(opt_name,"-%s_%d",MaterialConst_MaterialType_member_names[1],region_id);
	ierr = PetscOptionsGetInt(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_MaterialTypeSetField_plastic_type(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_MaterialType_member_names[1],region_id);CHKERRQ(ierr);
	}

    /* softening_type */
	value = 1.0; /* default - not required is nothing happens if option not found */
	sprintf(opt_name,"-%s_%d",MaterialConst_MaterialType_member_names[2],region_id);
	ierr = PetscOptionsGetInt(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_MaterialTypeSetField_softening_type(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_MaterialType_member_names[2],region_id);CHKERRQ(ierr);
	}

    /* density_type */
	value = 1.0; /* default - not required is nothing happens if option not found */
	sprintf(opt_name,"-%s_%d",MaterialConst_MaterialType_member_names[3],region_id);
	ierr = PetscOptionsGetInt(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_MaterialTypeSetField_density_type(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_MaterialType_member_names[3],region_id);CHKERRQ(ierr);
	}

    
    
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}
/* DensityConst */
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetDefault_DensityConst"
PetscErrorCode MaterialConstantsSetDefault_ViscosityConst(DataBucket db)
{
	int                          r,nregions;
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(db,&nregions,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(db,MaterialConst_DensityConst_classname,&PField);
	
	data = (MaterialConst_DensityConst*)PField->data; /* should write a function to do this */
	for (r=0; r<nregions; r++) {
		data[r].density = 1.0;
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetFromOptions_DensityConst"
PetscErrorCode MaterialConstantsSetFromOptions_DensityConst(DataBucket db,const char model_name[],const int region_id,PetscBool essential)
{
	char                         opt_name[256];
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	PetscBool                    found;
	PetscReal                    value;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_DensityConst_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* eta0 */
	value = 1.0; /* default - not required is nothing happens if option not found */
	sprintf(opt_name,"-%s_%d",MaterialConst_DensityConst_member_names[0],region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_DensityConstSetField_density(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_DensityConst_member_names[0],region_id);CHKERRQ(ierr);
	}
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}



/* ViscosityConst */
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetDefault_ViscosityConst"
PetscErrorCode MaterialConstantsSetDefault_ViscosityConst(DataBucket db)
{
	int                          r,nregions;
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(db,&nregions,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityConst_classname,&PField);
	
	data = (MaterialConst_ViscosityConst*)PField->data; /* should write a function to do this */
	for (r=0; r<nregions; r++) {
		data[r].eta0 = 1.0;
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetFromOptions_ViscosityConst"
PetscErrorCode MaterialConstantsSetFromOptions_ViscosityConst(DataBucket db,const char model_name[],const int region_id,PetscBool essential)
{
	char                         opt_name[256];
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	PetscBool                    found;
	PetscReal                    value;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityConst_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* eta0 */
	value = 1.0; /* default - not required is nothing happens if option not found */
	sprintf(opt_name,"-%s_%d",MaterialConst_ViscosityConst_member_names[0],region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_ViscosityConstSetField_eta0(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_ViscosityConst_member_names[0],region_id);CHKERRQ(ierr);
	}
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

/* ViscosityZ */
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetDefault_ViscosityZ"
PetscErrorCode MaterialConstantsSetDefault_ViscosityZ(DataBucket db)
{
	int                          r,nregions;
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(db,&nregions,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityConst_classname,&PField);
	
	data = (MaterialConst_ViscosityConst*)PField->data; /* should write a function to do this */
	for (r=0; r<nregions; r++) {
		data[r].eta0 = 1.0;
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetFromOptions_ViscosityZ"
PetscErrorCode MaterialConstantsSetFromOptions_ViscosityZ(DataBucket db,const char model_name[],const int region_id,PetscBool essential)
{
	char                         opt_name[256];
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	PetscBool                    found;
	PetscReal                    value;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityZ_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* eta0 */
	value = 1.0; /* default - not required is nothing happens if option not found */

	sprintf(opt_name,"-%s_%d",MaterialConst_ViscosityZ_member_names[0],region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_ViscosityZSetField_eta0(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_ViscosityZ_member_names[0],region_id);CHKERRQ(ierr);
	}
    
    sprintf(opt_name,"-%s_%d",MaterialConst_ViscosityZ_member_names[1],region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_ViscosityZSetField_zeta(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_ViscosityZ_member_names[1],region_id);CHKERRQ(ierr);
	}
    
    sprintf(opt_name,"-%s_%d",MaterialConst_ViscosityZ_member_names[2],region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_ViscosityZSetField_zref(data,value);
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,MaterialConst_ViscosityZ_member_names[2],region_id);CHKERRQ(ierr);
	}
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

/* ViscosityConst */
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetDefault_PlasticMises"
PetscErrorCode MaterialConstantsSetDefault_PlasticMises(DataBucket db)
{
	int                          r,nregions;
	DataField                    PField;
	MaterialConst_PlasticMises   *data;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(db,&nregions,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticMises_classname,&PField);
	
	data = (MaterialConst_PlasticMises*)PField->data; /* should write a function to do this */
	for (r=0; r<nregions; r++) {
		data[r].tau_yield     = 1.0e32;
		data[r].tau_yield_inf = 1.0e32;
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetFromOptions_PlasticMises"
PetscErrorCode MaterialConstantsSetFromOptions_PlasticMises(DataBucket db,const char model_name[],const int region_id,PetscBool essential)
{
	char                         opt_name[256],*field_name;
	DataField                    PField;
	MaterialConst_PlasticMises   *data;
	PetscBool                    found;
	PetscReal                    value;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticMises_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* tau_yield */
	value      = 1.0; /* default - not required is nothing happens if option not found */
	field_name = MaterialConst_PlasticMises_member_names[0];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticMisesSetField_yield_stress(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}

	/* tau_yield_inf */
	value      = 1.0; /* default - not required is nothing happens if option not found */
	field_name = MaterialConst_PlasticMises_member_names[1];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticMisesSetField_yield_stress_inf(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}
	
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}




/*DP*/ 
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetDefault_PlasticDP"
PetscErrorCode MaterialConstantsSetDefault_PlasticDP(DataBucket db)
{
	int                          r,nregions;
	DataField                    PField;
	MaterialConst_PlasticDP      *data;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(db,&nregions,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticDP_classname,&PField);
	
	data = (MaterialConst_PlasticDP*)PField->data; /* should write a function to do this */
	for (r=0; r<nregions; r++) {
		data[r].Co      = 1.0e32;
		data[r].Co_inf  = 1.0e32;
        data[r].phi     = 0.6;
		data[r].phi_inf = 0.6;
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetFromOptions_PlasticDP"
PetscErrorCode MaterialConstantsSetFromOptions_PlasticDP(DataBucket db,const char model_name[],const int region_id,PetscBool essential)
{
	char                         opt_name[256],*field_name;
	DataField                    PField;
	MaterialConst_PlasticMises   *data;
	PetscBool                    found;
	PetscReal                    value;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticDP_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* friction */
	value      = 1.0; /* default - not required as nothing happens if option not found */
	field_name = MaterialConst_PlasticDP_member_names[0];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticDPSetField_friction(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}
    
	/* friction_inf */
	value      = 1.0; /* default - not required as nothing happens if option not found */
	field_name = MaterialConst_PlasticDP_member_names[1];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticDPSetField_friction_inf(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}
	/* Cohesion */
	value      = 1.0; /* default - not required as nothing happens if option not found */
	field_name = MaterialConst_PlasticDP_member_names[2];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticDPSetField_cohesion(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}
    
	/* cohesion_inf */
	value      = 1.0; /* default - not required as nothing happens if option not found */
	field_name = MaterialConst_PlasticDP_member_names[3];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticDPSetField_cohesion_inf(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}
	
    
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsReportParseError"
PetscErrorCode MaterialConstantsReportParseError(const char model_name[],const char field_name[],const int region)
{
	PetscFunctionBegin;
	if (model_name) {
		PetscPrintf(PETSC_COMM_WORLD,"Expected user to provide option for field (%s) in region (%D) via -%s%s_%D \n",field_name,region,  model_name,field_name,region);			
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Expected user to provide option for field (%s) in region (%D) via -%s_%D \n",field_name,region,    field_name,region);			
	}
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Option not found");
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetDefaults"
PetscErrorCode MaterialConstantsSetDefaults(DataBucket db)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MaterialConstantsSetDefault_MaterialType(db);CHKERRQ(ierr);
	ierr = MaterialConstantsSetDefault_ViscosityConst(db);CHKERRQ(ierr);
    ierr = MaterialConstantsSetDefault_ViscosityZ(db);CHKERRQ(ierr);
	ierr = MaterialConstantsSetDefault_PlasticMises(db);CHKERRQ(ierr);
    ierr = MaterialConstantsSetDefault_PlasticDP(db);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}




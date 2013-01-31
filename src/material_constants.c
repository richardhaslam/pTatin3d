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
    DataBucketRegisterField(db,MaterialConst_DensityConst_classname,      sizeof(MaterialConst_DensityConst),PETSC_NULL);
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
	MaterialConst_MaterialType   *data;
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

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsPrintValues_MaterialType"
PetscErrorCode MaterialConstantsPrintValues_MaterialType(DataBucket db,const int region_id)
{
	char                         opt_name[256];
	DataField                    PField;
	MaterialConst_MaterialType   *data;
	PetscInt                     value;
	
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_MaterialType_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	
	sprintf(opt_name,"-%s_%d",MaterialConst_MaterialType_member_names[0],region_id);
    MaterialConst_MaterialTypeGetField_visc_type(data,&value);
	PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %d  \n", opt_name ,value);

    sprintf(opt_name,"-%s_%d",MaterialConst_MaterialType_member_names[1],region_id);
	MaterialConst_MaterialTypeGetField_plastic_type(data,&value);
	PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %d  \n", opt_name ,value);	
    
	sprintf(opt_name,"-%s_%d",MaterialConst_MaterialType_member_names[2],region_id);
	MaterialConst_MaterialTypeGetField_softening_type(data,&value);
    PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %d  \n", opt_name ,value);
    
    sprintf(opt_name,"-%s_%d",MaterialConst_MaterialType_member_names[3],region_id);
	MaterialConst_MaterialTypeGetField_density_type(data,&value);
    PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %d  \n", opt_name ,value);
    
    
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetValues_MaterialType"
PetscErrorCode MaterialConstantsSetValues_MaterialType(DataBucket db,const int region_id,PetscInt visc_t,PetscInt plast_t,PetscInt soft_t,PetscInt dens_t)
{
	DataField                    PField;
	MaterialConst_MaterialType   *data;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_MaterialType_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* visc_type */
	if (visc_t != PETSC_DEFAULT) {
		MaterialConst_MaterialTypeSetField_visc_type(data,visc_t);
	}
	/* plastic_type */
	if (plast_t != PETSC_DEFAULT) {
		MaterialConst_MaterialTypeSetField_plastic_type(data,plast_t);
	}
	/* softening_type */
	if (soft_t != PETSC_DEFAULT) {
		MaterialConst_MaterialTypeSetField_softening_type(data,soft_t);
	}
	/* density_type */
	if (dens_t != PETSC_DEFAULT) {
		MaterialConst_MaterialTypeSetField_density_type(data,dens_t);
	}
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

/* DensityConst */
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetDefault_DensityConst"
PetscErrorCode MaterialConstantsSetDefault_DensityConst(DataBucket db)
{
	int                          r,nregions;
	DataField                    PField;
	MaterialConst_DensityConst   *data;
	
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
	MaterialConst_DensityConst   *data;
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

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsPrintValues_DensityConst"
PetscErrorCode MaterialConstantsPrintValues_DensityConst(DataBucket db,const int region_id)
{
	char                         opt_name[256];
	DataField                    PField;
	MaterialConst_DensityConst   *data;
	PetscReal                    value;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_DensityConst_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	
	sprintf(opt_name,"-%s_%d",MaterialConst_DensityConst_member_names[0],region_id);
    MaterialConst_DensityConstGetField_density(data,&value);
    PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
    	
	DataFieldRestoreAccess(PField);
  
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetValues_DensityConst"
PetscErrorCode MaterialConstantsSetValues_DensityConst(DataBucket db,const int region_id,PetscReal density)
{
	DataField                    PField;
	MaterialConst_DensityConst   *data;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_DensityConst_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	if (density != PETSC_DEFAULT) {
		MaterialConst_DensityConstSetField_density(data,density);
	}	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsScaleValues_DensityConst"
PetscErrorCode MaterialConstantsScaleValues_DensityConst(DataBucket db,const int region_id,PetscReal rho_star)
{
	DataField                    PField;
	MaterialConst_DensityConst   *data;
	PetscErrorCode               ierr;
	PetscReal density;
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_DensityConst_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
    MaterialConst_DensityConstGetField_density(data,&density);
    density=density/rho_star;
    MaterialConst_DensityConstSetField_density(data,density);
		
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

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsPrintValues_ViscosityConst"
PetscErrorCode MaterialConstantsPrintValues_ViscosityConst(DataBucket db,const int region_id)
{
	char                         opt_name[256];
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	PetscReal                    value;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityConst_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* eta0 */
	sprintf(opt_name,"-%s_%d",MaterialConst_ViscosityConst_member_names[0],region_id);
    MaterialConst_ViscosityConstGetField_eta0(data,&value);
    PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);

	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetValues_ViscosityConst"
PetscErrorCode MaterialConstantsSetValues_ViscosityConst(DataBucket db,const int region_id,PetscReal viscosity)
{
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityConst_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	if (viscosity != PETSC_DEFAULT) {
		MaterialConst_ViscosityConstSetField_eta0(data,viscosity);
	}
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsScaleValues_ViscosityConst"
PetscErrorCode MaterialConstantsScaleValues_ViscosityConst(DataBucket db,const int region_id,PetscReal eta_star)
{
	DataField                    PField;
	MaterialConst_ViscosityConst *data;
	PetscErrorCode               ierr;
	PetscReal viscosity;
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityConst_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
    MaterialConst_ViscosityConstGetField_eta0(data,&viscosity);
    viscosity=viscosity/eta_star;
    MaterialConst_ViscosityConstSetField_eta0(data,viscosity);
    
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
	MaterialConst_ViscosityZ     *data;
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
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsPrintValues_ViscosityZ"
PetscErrorCode MaterialConstantsPrintValues_ViscosityZ(DataBucket db,const int region_id)
{
	char                         opt_name[256];
	DataField                    PField;
	MaterialConst_ViscosityZ     *data;
	PetscReal                    value;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityZ_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* eta0 */
	
	sprintf(opt_name,"-%s_%d",MaterialConst_ViscosityZ_member_names[0],region_id);
		MaterialConst_ViscosityZGetField_eta0(data,&value);
	PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
	
    sprintf(opt_name,"-%s_%d",MaterialConst_ViscosityZ_member_names[1],region_id);
		MaterialConst_ViscosityZGetField_zeta(data,&value);
	PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
    
	sprintf(opt_name,"-%s_%d",MaterialConst_ViscosityZ_member_names[2],region_id);
		MaterialConst_ViscosityZGetField_zref(data,&value);
    PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetValues_ViscosityZ"
PetscErrorCode MaterialConstantsSetValues_ViscosityZ(DataBucket db,const int region_id,PetscReal eta0,PetscReal zeta,PetscReal zref)
{
	DataField                    PField;
	MaterialConst_ViscosityZ     *data;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityZ_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	if (eta0 != PETSC_DEFAULT) {
		MaterialConst_ViscosityZSetField_eta0(data,eta0);
	}
	if (zeta != PETSC_DEFAULT) {
		MaterialConst_ViscosityZSetField_zeta(data,zeta);
	}
	if (zref != PETSC_DEFAULT) {
		MaterialConst_ViscosityZSetField_zref(data,zref);
	}
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsScaleValues_ViscosityZ"
PetscErrorCode MaterialConstantsScaleValues_ViscosityZ(DataBucket db,const int region_id,PetscReal eta_star,PetscReal L_star)
{
	DataField                    PField;
	MaterialConst_ViscosityZ     *data;
	PetscErrorCode               ierr;
	PetscReal eta0,zeta,zref;
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_ViscosityZ_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
    MaterialConst_ViscosityZGetField_eta0(data,&eta0);
    eta0=eta0/eta_star;
	MaterialConst_ViscosityZSetField_eta0(data,eta0);
    
	MaterialConst_ViscosityZGetField_zeta(data,&zeta);
    zeta=zeta/L_star;
	MaterialConst_ViscosityZSetField_zeta(data,zeta);
	
	MaterialConst_ViscosityZGetField_zref(data,&zref);
	zref=zref/L_star;
	MaterialConst_ViscosityZSetField_zref(data,zref);	
	
    DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}



/* PlasticMises */
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
	field_name = (char*)MaterialConst_PlasticMises_member_names[0];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticMisesSetField_yield_stress(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}

	/* tau_yield_inf */
	value      = 1.0; /* default - not required is nothing happens if option not found */
	field_name = (char*)MaterialConst_PlasticMises_member_names[1];
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

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsPrintValues_PlasticMises"
PetscErrorCode MaterialConstantsPrintValues_PlasticMises(DataBucket db,const int region_id)
{
	char                         opt_name[256],*field_name;
	DataField                    PField;
	MaterialConst_PlasticMises   *data;
	PetscReal                    value;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticMises_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* insert options */
	/* tau_yield */
	field_name = (char*)MaterialConst_PlasticMises_member_names[0];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
    MaterialConst_PlasticMisesGetField_yield_stress(data,&value); /* use setter */
	PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
    
	field_name = (char*)MaterialConst_PlasticMises_member_names[1];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	MaterialConst_PlasticMisesGetField_yield_stress_inf(data,&value); /* use setter */
	PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
    
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetValues_PlasticMises"
PetscErrorCode MaterialConstantsSetValues_PlasticMises(DataBucket db,const int region_id,PetscReal yield_stress,PetscReal yield_stress_inf)
{
	DataField                    PField;
	MaterialConst_PlasticMises   *data;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticMises_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	if (yield_stress != PETSC_DEFAULT) {
		MaterialConst_PlasticMisesSetField_yield_stress(data,yield_stress); /* use setter */
	}
	if (yield_stress_inf != PETSC_DEFAULT) {
		MaterialConst_PlasticMisesSetField_yield_stress_inf(data,yield_stress_inf); /* use setter */
	} 
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsScaleValues_PlasticMises"
PetscErrorCode MaterialConstantsScaleValues_PlasticMises(DataBucket db,const int region_id,PetscReal stress_star)
{
	DataField                    PField;
	MaterialConst_PlasticMises   *data;
	PetscErrorCode               ierr;
	PetscReal yield_stress, yield_stress_inf;
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticMises_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	
    MaterialConst_PlasticMisesGetField_yield_stress(data,&yield_stress); 
    yield_stress = yield_stress / stress_star;
    MaterialConst_PlasticMisesSetField_yield_stress(data,yield_stress); 
	
    MaterialConst_PlasticMisesGetField_yield_stress_inf(data,&yield_stress_inf);
    yield_stress_inf = yield_stress_inf / stress_star;
    MaterialConst_PlasticMisesSetField_yield_stress_inf(data,yield_stress_inf);
    
	
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
	MaterialConst_PlasticDP      *data;
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
	field_name = (char*)MaterialConst_PlasticDP_member_names[0];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticDPSetField_friction(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}
    
	/* friction_inf */
	value      = 1.0; /* default - not required as nothing happens if option not found */
	field_name = (char*)MaterialConst_PlasticDP_member_names[1];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticDPSetField_friction_inf(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}
	/* Cohesion */
	value      = 1.0; /* default - not required as nothing happens if option not found */
	field_name = (char*)MaterialConst_PlasticDP_member_names[2];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	ierr = PetscOptionsGetReal(model_name,opt_name,&value,&found);CHKERRQ(ierr);
	if (found) {
		MaterialConst_PlasticDPSetField_cohesion(data,value); /* use setter */
	} else if ( (!found)  && (essential) ) {
		ierr = MaterialConstantsReportParseError(model_name,field_name,region_id);CHKERRQ(ierr);
	}
    
	/* cohesion_inf */
	value      = 1.0; /* default - not required as nothing happens if option not found */
	field_name = (char*)MaterialConst_PlasticDP_member_names[3];
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
#define __FUNCT__ "MaterialConstantsPrintValues_PlasticDP"
PetscErrorCode MaterialConstantsPrintValues_PlasticDP(DataBucket db,const int region_id)
{
	char                         opt_name[256],*field_name;
	DataField                    PField;
	MaterialConst_PlasticDP      *data;
	PetscReal                    value;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticDP_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	
	/* friction */
	field_name = (char*)MaterialConst_PlasticDP_member_names[0];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	MaterialConst_PlasticDPGetField_friction(data,&value); /* use setter */
PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
    
	/* friction_inf */
	field_name = (char*)MaterialConst_PlasticDP_member_names[1];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
    MaterialConst_PlasticDPGetField_friction_inf(data,&value); /* use setter */
PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
	/* Cohesion */
	field_name = (char*)MaterialConst_PlasticDP_member_names[2];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
		MaterialConst_PlasticDPGetField_cohesion(data,&value); /* use setter */
    PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
	/* cohesion_inf */
	field_name = (char*)MaterialConst_PlasticDP_member_names[3];
	sprintf(opt_name,"-%s_%d",field_name,region_id);
	MaterialConst_PlasticDPGetField_cohesion_inf(data,&value); /* use setter */
PetscPrintf(PETSC_COMM_WORLD,"Current Value %s   :  %1.4e  \n", opt_name ,value);
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetValues_PlasticDP"
PetscErrorCode MaterialConstantsSetValues_PlasticDP(DataBucket db,const int region_id,PetscReal friction,PetscReal friction_inf,PetscReal cohesion,PetscReal cohesion_inf)
{
	DataField                    PField;
	MaterialConst_PlasticDP      *data;
	PetscErrorCode               ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticDP_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	if (friction != PETSC_DEFAULT) {
		MaterialConst_PlasticDPSetField_friction(data,friction); /* use setter */
	}
	if (friction_inf != PETSC_DEFAULT) {
		MaterialConst_PlasticDPSetField_friction_inf(data,friction_inf); /* use setter */
	}
	if (cohesion != PETSC_DEFAULT) {
		MaterialConst_PlasticDPSetField_cohesion(data,cohesion); /* use setter */
	}
	if (cohesion_inf != PETSC_DEFAULT) {
		MaterialConst_PlasticDPSetField_cohesion_inf(data,cohesion_inf); /* use setter */
	}
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsScaleValues_PlasticDP"
PetscErrorCode MaterialConstantsScaleValues_PlasticDP(DataBucket db,const int region_id, PetscReal stress_star)
{
	DataField                    PField;
	MaterialConst_PlasticDP      *data;
	PetscErrorCode               ierr;
	PetscReal cohesion,cohesion_inf;
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(db,MaterialConst_PlasticDP_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	MaterialConst_PlasticDPGetField_cohesion(data,&cohesion); /* use setter */
    cohesion = cohesion/stress_star;
	MaterialConst_PlasticDPSetField_cohesion(data,cohesion); /* use setter */
    
	MaterialConst_PlasticDPGetField_cohesion_inf(data,&cohesion_inf); /* use setter */
	cohesion_inf = cohesion_inf/stress_star;
    MaterialConst_PlasticDPSetField_cohesion_inf(data,cohesion_inf); /* use setter */

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
    ierr = MaterialConstantsSetDefault_DensityConst(db);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsSetFromOptions"
PetscErrorCode MaterialConstantsSetFromOptions(DataBucket db,const char model_name[],const int region_id,PetscBool essential)
{
    DataField                    PField;
	MaterialConst_MaterialType   *data;
	PetscInt                     value;
	PetscErrorCode ierr;
    
	PetscFunctionBegin;
    
    ierr= MaterialConstantsSetFromOptions_MaterialType(db,model_name,region_id,essential);CHKERRQ(ierr);
    
    DataBucketGetDataFieldByName(db,MaterialConst_MaterialType_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
    
    MaterialConst_MaterialTypeGetField_visc_type(data,&value);
    
    switch (value) {
        case VISCOUS_CONSTANT:
            ierr= MaterialConstantsSetFromOptions_ViscosityConst(db,model_name,region_id,essential);CHKERRQ(ierr);
            break;
        case VISCOUS_FRANKK:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Rheology FRANKK NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;        
        case VISCOUS_Z:
            ierr= MaterialConstantsSetFromOptions_ViscosityZ(db,model_name,region_id,essential);CHKERRQ(ierr);
            break; 
        case VISCOUS_ARRHENIUS:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Rheology ARRHENIUS NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break; 
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"VISCOUS TYPE UNDEFINED");
            break;
    }
	MaterialConst_MaterialTypeGetField_plastic_type(data,&value);
    switch (value) {
        case PLASTIC_NONE:
            
            break;
        case PLASTIC_MISES:
            ierr= MaterialConstantsSetFromOptions_PlasticMises(db,model_name,region_id,essential);CHKERRQ(ierr);
            break;    
        case PLASTIC_DP:
            ierr= MaterialConstantsSetFromOptions_PlasticDP(db,model_name,region_id,essential);CHKERRQ(ierr);
            break;    
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"PLASTIC TYPE UNDEFINED");
            break;
    }
    
    MaterialConst_MaterialTypeGetField_softening_type(data,&value);
    switch (value) {
        case SOFTENING_NONE:
            
            break;
        case SOFTENING_LINEAR:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"SOFTENING_LINEAR NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;
        case SOFTENING_EXPONENTIAL:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"SOFTENING_EXPONENTIAL NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;        
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"DENSITY TYPE UNDEFINED");
            break;
    }
    
    MaterialConst_MaterialTypeGetField_density_type(data,&value);
    switch (value) {
        case DENSITY_CONSTANT:
            ierr= MaterialConstantsSetFromOptions_DensityConst(db,model_name,region_id,essential);CHKERRQ(ierr);            
            break;
        case DENSITY_BOUSINESQ:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Density BOUSSINESQ NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;    
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"DENSITY TYPE UNDEFINED");
            break;
    }
    
	DataFieldRestoreAccess(PField);
    
    
	PetscFunctionReturn(0);
	
}



#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsScaleAll"
PetscErrorCode MaterialConstantsScaleAll(DataBucket db,const int region_id,PetscReal L_star, PetscReal U_star,PetscReal t_star,PetscReal eta_star,PetscReal rho_star,PetscReal P_star)
{
    DataField                    PField;
	MaterialConst_MaterialType   *data;
	PetscInt                     value;
	PetscErrorCode ierr;
    
	PetscFunctionBegin;
    
    DataBucketGetDataFieldByName(db,MaterialConst_MaterialType_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
    
    MaterialConst_MaterialTypeGetField_visc_type(data,&value);
    switch (value) {
        case VISCOUS_CONSTANT:
            ierr = MaterialConstantsScaleValues_ViscosityConst(db,region_id,eta_star);CHKERRQ(ierr);
            break;
        case VISCOUS_FRANKK:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Rheology FRANKK NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;        
        case VISCOUS_Z:
            ierr = MaterialConstantsScaleValues_ViscosityZ(db,region_id,eta_star,L_star);CHKERRQ(ierr);
            break; 
        case VISCOUS_ARRHENIUS:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Rheology ARRHENIUS NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break; 
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"VISCOUS TYPE UNDEFINED");
            break;
    }
	MaterialConst_MaterialTypeGetField_plastic_type(data,&value);
    switch (value) {
        case PLASTIC_NONE:
            
            break;
        case PLASTIC_MISES:
            ierr = MaterialConstantsScaleValues_PlasticMises(db,region_id,P_star);CHKERRQ(ierr);
            break;    
        case PLASTIC_DP:
            ierr = MaterialConstantsScaleValues_PlasticDP(db,region_id,P_star);CHKERRQ(ierr);
            break;    
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"PLASTIC TYPE UNDEFINED");
            break;
    }
    
    MaterialConst_MaterialTypeGetField_softening_type(data,&value);
    switch (value) {
        case SOFTENING_NONE:
            
            break;
        case SOFTENING_LINEAR:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"SOFTENING_LINEAR NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;
        case SOFTENING_EXPONENTIAL:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"SOFTENING_EXPONENTIAL NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;        
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"DENSITY TYPE UNDEFINED");
            break;
    }
    
    MaterialConst_MaterialTypeGetField_density_type(data,&value);
    switch (value) {
        case DENSITY_CONSTANT:
            ierr = MaterialConstantsScaleValues_DensityConst(db,region_id,rho_star);CHKERRQ(ierr);            
            break;
        case DENSITY_BOUSINESQ:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Density BOUSSINESQ NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;    
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"DENSITY TYPE UNDEFINED");
            break;
    }
    
	DataFieldRestoreAccess(PField);
    
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialConstantsPrintAll"
PetscErrorCode MaterialConstantsPrintAll(DataBucket db,const int region_id)
{
    DataField                    PField;
	MaterialConst_MaterialType   *data;
	PetscInt                     value;
	PetscErrorCode ierr;
    
	PetscFunctionBegin;
    
    DataBucketGetDataFieldByName(db,MaterialConst_MaterialType_classname,&PField);
	
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
    
    MaterialConst_MaterialTypeGetField_visc_type(data,&value);
    switch (value) {
        case VISCOUS_CONSTANT:
            ierr = MaterialConstantsPrintValues_ViscosityConst(db,region_id);CHKERRQ(ierr);
            break;
        case VISCOUS_FRANKK:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Rheology FRANKK NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;        
        case VISCOUS_Z:
            ierr = MaterialConstantsPrintValues_ViscosityZ(db,region_id);CHKERRQ(ierr);
            break; 
        case VISCOUS_ARRHENIUS:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Rheology ARRHENIUS NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break; 
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"VISCOUS TYPE UNDEFINED");
            break;
    }
	MaterialConst_MaterialTypeGetField_plastic_type(data,&value);
    switch (value) {
        case PLASTIC_NONE:
            
            break;
        case PLASTIC_MISES:
            ierr = MaterialConstantsPrintValues_PlasticMises(db,region_id);CHKERRQ(ierr);
            break;    
        case PLASTIC_DP:
            ierr = MaterialConstantsPrintValues_PlasticDP(db,region_id);CHKERRQ(ierr);
            break;    
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"PLASTIC TYPE UNDEFINED");
            break;
    }
    
    MaterialConst_MaterialTypeGetField_softening_type(data,&value);
    switch (value) {
        case SOFTENING_NONE:
            
            break;
        case SOFTENING_LINEAR:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"SOFTENING_LINEAR NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;
        case SOFTENING_EXPONENTIAL:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"SOFTENING_EXPONENTIAL NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;        
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"DENSITY TYPE UNDEFINED");
            break;
    }
    
    MaterialConst_MaterialTypeGetField_density_type(data,&value);
    switch (value) {
        case DENSITY_CONSTANT:
    ierr = MaterialConstantsPrintValues_DensityConst(db,region_id);CHKERRQ(ierr);            
            break;
        case DENSITY_BOUSINESQ:
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Density BOUSSINESQ NOT IMPLEMENTED CHECK YOUR OPTIONS");
            break;    
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"DENSITY TYPE UNDEFINED");
            break;
    }
    
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}



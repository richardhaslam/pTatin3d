

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
	
	DataBucketRegisterField(db,MaterialConst_MaterialType_classname,    sizeof(MaterialConst_MaterialType),PETSC_NULL);
	DataBucketRegisterField(db,MaterialConst_ViscosityConst_classname,  sizeof(MaterialConst_ViscosityConst),PETSC_NULL);

	DataBucketFinalize(db);
	

	DataBucketSetInitialSizes(db,200,0);
	ierr = MaterialConstantsSetDefaults(db);CHKERRQ(ierr);
	
	
	*_db = db;
	
	PetscFunctionReturn(0);
}

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
#define __FUNCT__ "MaterialConstantsSetDefaults"
PetscErrorCode MaterialConstantsSetDefaults(DataBucket db)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MaterialConstantsSetDefault_MaterialType(db);CHKERRQ(ierr);
	ierr = MaterialConstantsSetDefault_ViscosityConst(db);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}




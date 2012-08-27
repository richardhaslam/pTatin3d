
#ifndef __material_constants_h__
#define __material_constants_h__

/* add auto geneterated classes here */
#include "material_constants/MaterialConst_MaterialType_def.h"
#include "material_constants/MaterialConst_ViscosityConst_def.h"
#include "material_constants/MaterialConst_PlasticMises_def.h"

PetscErrorCode MaterialConstantsInitialize(DataBucket *_db);
PetscErrorCode MaterialConstantsSetDefaults(DataBucket db);

PetscErrorCode MaterialConstantsSetFromOptions_ViscosityConst(DataBucket db,const char model_name[],const int region_id,PetscBool essential);
PetscErrorCode MaterialConstantsSetFromOptions_PlasticMises(DataBucket db,const char model_name[],const int region_id,PetscBool essential);

PetscErrorCode MaterialConstantsReportParseError(const char model_name[],const char field_name[],const int region);

#endif

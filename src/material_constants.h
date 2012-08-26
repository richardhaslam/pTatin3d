
#ifndef __material_constants_h__
#define __material_constants_h__

#include "material_constants/MaterialConst_MaterialType_def.h"
#include "material_constants/MaterialConst_ViscosityConst_def.h"

PetscErrorCode MaterialConstantsInitialize(DataBucket *_db);
PetscErrorCode MaterialConstantsSetDefaults(DataBucket db);

PetscErrorCode MaterialConstantsSetFromOptions_ViscosityConst(DataBucket db,const char model_name[],const int region_id,PetscBool essential);

#endif

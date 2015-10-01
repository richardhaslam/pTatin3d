
#ifndef __foundation_h__
#define __foundation_h__

#include "petsc.h"
#include "petscvec.h"
#include "cJSON.h"
#include "cJSON_helpers.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"


typedef struct _p_Foundation *Foundation;
typedef struct _p_FoundationUserVars *FoundationUserVars;
typedef struct _p_FoundationOutput *FoundationOutput;

typedef enum {
  FND_MeshGeomT_Cartesian=0,
  FND_MeshGeomT_Hex,
  FND_MeshGeomT_DMDAFile,
  FND_MeshGeomT_Lua,
  FND_MeshGeomT_NULL
} FND_MeshGeomType;
extern const char *FND_MeshGeomTypeNames[];

typedef enum {
  FND_MeshRefT_A=0,
  FND_MeshRefT_B,
  FND_MeshRefT_Lua,
  FND_MeshRefT_NULL
} FND_MeshRefinementType;
extern const char *FND_MeshRefinementTypeNames[];

typedef enum {
  FND_MatRegT_JSON=0,
  FND_MatRegT_DBFile,
  FND_MatRegT_Lua,
  FND_MatRegT_NULL
} FND_MaterialRegionType;
extern const char *FND_MaterialRegionTypeNames[];


typedef enum {
  FND_BCType_Dirichlet = 0,
  FND_BCType_Neumann
} FND_BoundaryConditionType;
extern const char *FND_BoundaryConditionTypeNames[];

typedef enum {
  FND_BCMode_Scalar = 0,
  FND_BCMode_Component,
  FND_BCMode_Vector,
  FND_BCMode_Tensor,
  FND_BCMode_User,
  FND_BCMode_Lua
} FND_BoundaryConditionMode;
extern const char *FND_BoundaryConditionModeNames[];

typedef enum {
  FND_OutputType_MeshDefault = 0,
  FND_OutputType_MeshLiteDefault,
  FND_OutputType_MaterialPointDefault,
  FND_OutputType_User,
} FND_OutputType;
extern const char *FND_OutputTypeNames[];


/* foundation parser helpers */
PetscErrorCode FoundationGetPtatinCtx(Foundation f,pTatinCtx *ctx);
PetscErrorCode FoundationGetUserVariables(Foundation f,FoundationUserVars *fv);
PetscErrorCode FoundationParseJSONGetItemEssential(cJSON *jroot,const char itemname[],cJSON **jitem);
PetscErrorCode FoundationParseJSONGetItemOptional(cJSON *jroot,const char itemname[],cJSON **jitem);
PetscErrorCode FoundationParseMaterialMetaData(pTatinCtx c,Foundation f);
PetscErrorCode ModelInitialMaterialState_Foundation(pTatinCtx c,void *ctx);

PetscErrorCode ModelInitialize_Foundation(pTatinCtx c,void *ctx);

PetscErrorCode ModelInitialMeshGeometry_Foundation(pTatinCtx c,void *ctx);
PetscErrorCode ModelInitialMaterialGeometry_Foundation(pTatinCtx c,void *ctx);

PetscErrorCode ModelInitialConditionMeshFields_Foundation(pTatinCtx c,Vec X,void *ctx);
PetscErrorCode ModelInitialConditonMarkerFields_Foundation(pTatinCtx c,Vec X,void *ctx);

PetscErrorCode ModelApplyBC_Foundation(pTatinCtx c,void *ctx);
PetscErrorCode ModelApplyBCMG_Foundation(PetscInt nl,BCList bclist[],DM dmv[],pTatinCtx c,void *ctx);

//PetscErrorCode ModelUpdateMeshGeometry_Foundation(pTatinCtx c,Vec X,void *ctx);

PetscErrorCode ModelOutput_Foundation(pTatinCtx c,Vec X,const char prefix[],void *ctx);
PetscErrorCode ModelDestroy_Foundation(pTatinCtx c,void *ctx);



PetscErrorCode FoundationUserRegisterFunctions(Foundation f);

#endif


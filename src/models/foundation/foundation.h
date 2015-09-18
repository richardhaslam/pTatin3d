
#ifndef __foundation_h__
#define __foundation_h__

#include "petsc.h"
#include "petscvec.h"
#include "cJSON.h"
#include "cJSON_helpers.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"


typedef struct _p_Foundation *Foundation;
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
  FND_MatRefT_JSON=0,
  FND_MatRegT_DBFile,
  FND_MatRegT_Lua,
  FND_MatRegT_NULL
} FND_MaterialRegionType;
extern const char *FND_MaterialRegionTypeNames[];


/* foundation parser helpers */
PetscErrorCode FoundationParseJSONGetItemEssential(cJSON *jroot,const char itemname[],cJSON **jitem);
PetscErrorCode FoundationParseJSONGetItemOptional(cJSON *jroot,const char itemname[],cJSON **jitem);
PetscErrorCode FoundationParseMaterialMetaData(pTatinCtx c,Foundation f);

PetscErrorCode ModelInitialize_Foundation(pTatinCtx c,void *ctx);
PetscErrorCode ModelInitialMeshGeometry_Foundation(pTatinCtx c,void *ctx);
PetscErrorCode ModelInitialMaterialGeometry_Foundation(pTatinCtx c,void *ctx);
PetscErrorCode ModelApplyBoundaryCondition_Foundation(pTatinCtx c,void *ctx);
PetscErrorCode ModelApplyBoundaryConditionMG_Foundation(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx);
PetscErrorCode ModelApplyUpdateMeshGeometry_Foundation(pTatinCtx c,Vec X,void *ctx);
PetscErrorCode ModelOutput_Foundation(pTatinCtx c,Vec X,const char prefix[],void *ctx);
PetscErrorCode ModelDestroy_Foundation(pTatinCtx c,void *ctx);

#endif


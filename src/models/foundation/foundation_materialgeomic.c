
#include "foundation.h"
#include "foundation_impl.h"
#include "geometry_object.h"
#include "geometry_object_parse.h"
#include "data_bucket.h"
#include "MPntStd_def.h"

#undef __FUNCT__
#define __FUNCT__ "MaterialICApply_JSON"
PetscErrorCode MaterialICApply_JSON(pTatinCtx c,Foundation f,cJSON *jitem,DataBucket matpoint)
{
  char *filename;
  int found,n_mp_points;
  PetscInt g,ng,p;
  GeometryObject *golist;
	DataField       PField_std;
  PetscErrorCode  ierr;
  
  cJSON_GetObjectValue_char(jitem,"RegionGeomFile",&found,&filename);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.RegionGeomFile missing");

  ierr = GeometryObjectLoadJSON(filename,&ng,&golist);CHKERRQ(ierr);

  DataBucketGetSizes(matpoint,&n_mp_points,0,0);

	DataBucketGetDataFieldByName(matpoint,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	for (p=0; p<n_mp_points; p++) {
		MPntStd *mp;
    double  *pos;
    int     inside,region_id;
    
		DataFieldAccessPoint(PField_std,p,   (void**)&mp);
		MPntStdGetField_global_coord(mp,&pos);
    
    region_id = -1;
    for (g=0; g<ng; g++) {
      ierr = GeometryObjectPointInside(golist[g],pos,&inside);CHKERRQ(ierr);
      if (inside) {
        region_id = g;
      }
    }
    if (region_id == -1) {
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"MatGeomIC.JSON: Failed to find geometry object containing material point");
    }
    
		MPntStdSetField_phase_index(mp,region_id);
    
  }
	DataFieldRestoreAccess(PField_std);
  
  
  /* tidy up geometry objects */
  for (g=0; g<ng; g++) {
    ierr = GeometryObjectDestroy(&golist[g]);CHKERRQ(ierr);
  }
  free(golist);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialMaterialGeometry_Foundation"
PetscErrorCode ModelInitialMaterialGeometry_Foundation(pTatinCtx c,void *ctx)
{
  Foundation      data = (Foundation)ctx;
  cJSON           *jroot,*jitem;
  int             i;
  DataBucket      matpoint;
  PetscErrorCode  ierr;
  
  /* look for MeshGeomICType */
  ierr = FoundationParseJSONGetItemEssential(data->root,"MatGeomICType",&jroot);CHKERRQ(ierr);

  /* parse type */
  jitem = NULL;
  for (i=0; i<(PetscInt)FND_MatRegT_NULL-1; i++) {
    ierr = FoundationParseJSONGetItemOptional(jroot,FND_MaterialRegionTypeNames[i],&jitem);CHKERRQ(ierr);
    if (jitem) {
      data->material_region_type = (FND_MaterialRegionType)i;
      break;
    }
  }
  
  /* apply method */
  ierr = pTatinGetMaterialPoints(c,&matpoint,NULL);CHKERRQ(ierr);
  switch (data->material_region_type) {
      
    case FND_MatRefT_JSON:
      ierr = MaterialICApply_JSON(c,data,jitem,matpoint);CHKERRQ(ierr);
      break;
      
    case FND_MatRegT_DBFile:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MatGeomICType.DBFile not implemented");
      break;
      
    case FND_MatRegT_Lua:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MatGeomICType.Lua not implemented");
      break;
      
    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"MatGeomICType.method missing: A valid constructor must be specified");
      break;
  }
  
  PetscFunctionReturn(0);
}






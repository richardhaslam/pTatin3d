
#include "foundation.h"
#include "foundation_impl.h"


typedef struct {
  PetscBool quasi_two_dimensional;
  double xx[2],yy[2],zz[2];
} MIC_Cartesian;

typedef struct {
  double xx[2],yy[2],zz[2];
} MIC_Hex;

typedef struct {
  char filename[PETSC_MAX_PATH_LEN];
} MIC_Lua;

typedef struct {
  MIC_Cartesian cart;
  MIC_Hex hex;
  MIC_Lua lua;
} MeshICCtx;


#undef __FUNCT__
#define __FUNCT__ "MeshICParse_Cartesian"
PetscErrorCode MeshICParse_Cartesian(cJSON *jobj,MIC_Cartesian *ctx)
{
  cJSON *jobjfield;
  int found,nv;
  PetscErrorCode ierr;
  
  /* essential */
  ierr = FoundationParseJSONGetItemEssential(jobj,"xrange",&jobjfield);CHKERRQ(ierr);
  if (cJSON_GetArraySize(jobjfield) != 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cartesian: Requires array \"xrange\" be of length 2");
  cJSON_GetObjectValue_doublearray(jobj,"xrange",&found,&nv,ctx->xx);
  
  ierr = FoundationParseJSONGetItemEssential(jobj,"yrange",&jobjfield);CHKERRQ(ierr);
  if (cJSON_GetArraySize(jobjfield) != 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cartesian: Requires array \"yrange\" be of length 2");
  cJSON_GetObjectValue_doublearray(jobj,"yrange",&found,&nv,ctx->yy);
  
  ierr = FoundationParseJSONGetItemEssential(jobj,"zrange",&jobjfield);CHKERRQ(ierr);
  if (cJSON_GetArraySize(jobjfield) != 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cartesian: Requires array \"zrange\" be of length 2");
  cJSON_GetObjectValue_doublearray(jobj,"zrange",&found,&nv,ctx->zz);
  
  /* optional */
  ctx->quasi_two_dimensional = PETSC_FALSE;
  cJSON_GetObjectValue_logical(jobj,"2dmode",&found,&nv);
  if (found) { ctx->quasi_two_dimensional = (PetscBool)nv; }
   
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MeshICApply_Cartesian"
PetscErrorCode MeshICApply_Cartesian(pTatinCtx c,Foundation f,MIC_Cartesian *cart,DM dav)
{
  PetscErrorCode ierr;

  ierr = DMDASetUniformCoordinates(dav,cart->xx[0],cart->xx[1],cart->yy[0],cart->yy[1],cart->zz[0],cart->zz[1]);CHKERRQ(ierr);
  if (cart->quasi_two_dimensional) {
    ierr = pTatin3d_DefineVelocityMeshGeometryQuasi2D(c);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialMeshGeometry_Foundation"
PetscErrorCode ModelInitialMeshGeometry_Foundation(pTatinCtx c,void *ctx)
{
  Foundation      data = (Foundation)ctx;
  cJSON           *jroot,*jitem;
  int             i;
  PhysCompStokes  stokes;
  DM              stokes_pack,dav,dap;
  MeshICCtx       mtype_data;
  PetscErrorCode  ierr;
  
  /* look for MeshGeomICType */
  ierr = FoundationParseJSONGetItemEssential(data->root,"MeshGeomICType",&jroot);CHKERRQ(ierr);

  /* parse type */
  jitem = NULL;
  for (i=0; i<(PetscInt)FND_MeshGeomT_NULL-1; i++) {
    ierr = FoundationParseJSONGetItemOptional(jroot,FND_MeshGeomTypeNames[i],&jitem);CHKERRQ(ierr);
    if (jitem) {
      data->mesh_geom_type = (FND_MeshGeomType)i;
      break;
    }
  }
  
  /* apply method */
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&stokes_pack);CHKERRQ(ierr);
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
  
  switch (data->mesh_geom_type) {
      
    case FND_MeshGeomT_Cartesian:
      ierr = MeshICParse_Cartesian(jitem,&mtype_data.cart);CHKERRQ(ierr);
      ierr = MeshICApply_Cartesian(c,data,&mtype_data.cart,dav);CHKERRQ(ierr);
      break;
      
    case FND_MeshGeomT_Hex:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MeshGeomICType.Hex not implemented");
      break;
      
    case FND_MeshGeomT_DMDAFile:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MeshGeomICType.DMDAFile not implemented");
      break;
      
    case FND_MeshGeomT_Lua:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MeshGeomICType.Lua not implemented");
      break;
      
    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"MeshGeomICType.method missing: A valid constructor must be specified");
      break;
  }
  
  PetscFunctionReturn(0);
}






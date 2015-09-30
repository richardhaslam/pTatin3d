
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"
#include "foundation_output.h"

#include "ptatin3d.h"
#include "material_point_utils.h"
#include "output_material_points.h"


#undef __FUNCT__
#define __FUNCT__ "FoundationOutput_DefaultMeshFields"
PetscErrorCode FoundationOutput_DefaultMeshFields(pTatinCtx ptatinctx,Foundation f,Vec X,const char prefix[])
{
  PetscErrorCode   ierr;
  PetscFunctionBegin;
  /* v and p are written */
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationOutput_DefaultMeshFieldsLite"
PetscErrorCode FoundationOutput_DefaultMeshFieldsLite(pTatinCtx ptatinctx,Foundation f,Vec X,const char prefix[])
{
  PetscErrorCode   ierr;
  PetscFunctionBegin;

  /* Light weight viewer: Only v is written. v and coords are expressed as floats */
  ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationOutput_DefaultMaterialPoints"
PetscErrorCode FoundationOutput_DefaultMaterialPoints(pTatinCtx ptatinctx,Foundation f,Vec X,const char prefix[])
{
  DataBucket       materialpoint;
  PetscErrorCode   ierr;
  
  PetscFunctionBegin;
  ierr = pTatinGetMaterialPoints(ptatinctx,&materialpoint,NULL);CHKERRQ(ierr);
  /* Customized viewer: User defines specific fields they want to view - NOTE not .pvd file will be created */
  {
    const int                 nf = 2;
    const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes };
    char                      mp_file_prefix[PETSC_MAX_PATH_LEN];
    
    PetscSNPrintf(mp_file_prefix,PETSC_MAX_PATH_LEN-1,"%s_mpoints",prefix);
    ierr = SwarmViewGeneric_ParaView(materialpoint,nf,mp_prop_list,ptatinctx->outputpath,mp_file_prefix);CHKERRQ(ierr);
  }
  
  /* Customized marker->cell viewer: Marker data is projected onto the velocity mesh. User defines specific fields */
  {
    const int                    nf = 3;
    const MaterialPointVariable  mp_prop_list[] = { MPV_region, MPV_viscosity, MPV_density };
    
    ierr = pTatin3d_ModelOutput_MarkerCellFields(ptatinctx,nf,mp_prop_list,prefix);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Foundation"
PetscErrorCode ModelOutput_Foundation(pTatinCtx ptatinctx,Vec X,const char prefix[],void *modelctx)
{
  Foundation       f = (Foundation)modelctx;
  PetscInt         i;
  PetscErrorCode   ierr;
  
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

  /*
  ierr = FoundationOutput_DefaultMeshFields(ptatinctx,f,X,prefix);CHKERRQ(ierr);
  ierr = FoundationOutput_DefaultMeshFieldsLite(ptatinctx,f,X,prefix);CHKERRQ(ierr);
  ierr = FoundationOutput_DefaultMaterialPoints(ptatinctx,f,X,prefix);CHKERRQ(ierr);
  ierr = FoundationOutput_DefaultMaterialPoints(ptatinctx,f,X,prefix);CHKERRQ(ierr);
  
  ierr = FoundationUserApply_Output(f,"MaterialPointStdViewer",X,prefix);CHKERRQ(ierr);
  */

  for (i=0; i<f->output->n_modes; i++) {
    switch (f->output->typelist[i]) {
      case FND_OutputType_MeshDefault:
        ierr = FoundationOutput_DefaultMeshFields(ptatinctx,f,X,prefix);CHKERRQ(ierr);
        break;

      case FND_OutputType_MeshLiteDefault:
        ierr = FoundationOutput_DefaultMeshFieldsLite(ptatinctx,f,X,prefix);CHKERRQ(ierr);
        break;

      case FND_OutputType_MaterialPointDefault:
        ierr = FoundationOutput_DefaultMaterialPoints(ptatinctx,f,X,prefix);CHKERRQ(ierr);
        break;

      case FND_OutputType_User:
        break;
    }
  }
  for (i=0; i<f->output->n_userfunc; i++) {
    ierr = f->output->userfuncs[i](f,X,prefix);CHKERRQ(ierr);
  }
  
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationOutputCreate"
PetscErrorCode FoundationOutputCreate(FoundationOutput *fo)
{
  FoundationOutput f;
  
  PetscMalloc(sizeof(struct _p_FoundationOutput),&f);
  PetscMemzero(f,sizeof(struct _p_FoundationOutput));
  f->n_userfunc = 0;
  *fo = f;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationOutputDestroy"
PetscErrorCode FoundationOutputDestroy(FoundationOutput *fo)
{
  FoundationOutput f;

  if (!fo) PetscFunctionReturn(0);
  f = *fo;
  PetscFree(f);
  *fo = NULL;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationOutputParse"
PetscErrorCode FoundationOutputParse(Foundation f,cJSON *root)
{
  cJSON *jroot,*typelist,*userlist,*jobj_k;
  int entry,nmodes;
  PetscBool match,user_found = PETSC_FALSE;
  PetscErrorCode ierr;
  
  /* look for MeshGeomICType */
  PetscPrintf(PETSC_COMM_WORLD,"[Foundation] Parsing \"ModelOutput\"\n");
  ierr = FoundationParseJSONGetItemEssential(root,"ModelOutput",&jroot);CHKERRQ(ierr);

  ierr = FoundationParseJSONGetItemEssential(jroot,"TypeList",&typelist);CHKERRQ(ierr);

  nmodes = cJSON_GetArraySize(typelist);
  printf("--------- Found %d modes \n",nmodes);
  f->output->n_modes = (PetscInt)nmodes;
  
  jobj_k = cJSON_GetArrayItemRoot(typelist);
  for (entry=0; entry<nmodes; entry++) {

    PetscStrcmp(jobj_k->valuestring,"MeshFields_Default",&match);
    if (match) { f->output->typelist[entry] = FND_OutputType_MeshDefault; }

    PetscStrcmp(jobj_k->valuestring,"MeshFieldsLite_Default",&match);
    if (match) { f->output->typelist[entry] = FND_OutputType_MeshLiteDefault; }

    PetscStrcmp(jobj_k->valuestring,"MaterialPoints_Default",&match);
    if (match) { f->output->typelist[entry] = FND_OutputType_MaterialPointDefault; }

    PetscStrcmp(jobj_k->valuestring,"User",&match);
    if (match) {
      f->output->typelist[entry] = FND_OutputType_User;
      user_found = PETSC_TRUE;
    }

    
    /*
    ierr = FoundationParseJSONGetItemOptional(jobj_k,"MeshFields_Default",&jdummy);CHKERRQ(ierr);
    if (match) { f->output->typelist[entry] = FND_OutputType_MeshDefault; continue; }

    ierr = FoundationParseJSONGetItemOptional(jobj_k,"MeshFieldsLite_Default",&jdummy);CHKERRQ(ierr);
    if (match) { f->output->typelist[entry] = FND_OutputType_MeshLiteDefault; continue; }
    
    ierr = FoundationParseJSONGetItemOptional(jobj_k,"MaterialPoints_Default",&jdummy);CHKERRQ(ierr);
    if (match) { f->output->typelist[entry] = FND_OutputType_MaterialPointDefault; continue; }

    ierr = FoundationParseJSONGetItemOptional(jobj_k,"User",&jdummy);CHKERRQ(ierr);
    if (match) {
      f->output->typelist[entry] = FND_OutputType_User;
      user_found = PETSC_TRUE;
      continue;
    }
    */
    jobj_k = cJSON_GetArrayItemNext(jobj_k);
  }

  if (user_found) {
    ierr = FoundationParseJSONGetItemEssential(jroot,"UserFunctionList",&userlist);CHKERRQ(ierr);
  
    nmodes = cJSON_GetArraySize(userlist);
    printf("--------- Found %d user functions \n",nmodes);
    f->output->n_userfunc = (PetscInt)nmodes;
    
    jobj_k = cJSON_GetArrayItemRoot(userlist);
    for (entry=0; entry<nmodes; entry++) {
      char                 *functionname;
      FoundationUserOutput fp;
      
      functionname = jobj_k->valuestring;
      
      ierr = PetscFunctionListFind(f->user->flist_output,(const char*)functionname,(void(**)(void))&fp);CHKERRQ(ierr);
      if (!fp) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserOutput: Method with name \"%s\" has not been registered",functionname);
      f->output->userfuncs[entry] = fp;
      
      jobj_k = cJSON_GetArrayItemNext(jobj_k);
    }
  }
  
  
  PetscFunctionReturn(0);
}

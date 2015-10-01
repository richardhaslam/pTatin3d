
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"
#include "foundation_output.h"

#include "ptatin3d.h"
#include "material_point_utils.h"
#include "output_material_points.h"
#include "energy_output.h"
#include "ptatin3d_energy.h"


#undef __FUNCT__
#define __FUNCT__ "FoundationOutput_DefaultMeshFields"
PetscErrorCode FoundationOutput_DefaultMeshFields(pTatinCtx ptatinctx,Foundation f,Vec X,const char prefix[])
{
  PetscErrorCode   ierr;
  PetscBool        active_energy;
  
  PetscFunctionBegin;
  /* v and p are written */
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);

	/* T written out */
	ierr = pTatinContextValid_Energy(ptatinctx,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
    PhysCompEnergy energy;
    Vec            temperature;
    
    ierr = pTatinGetContext_Energy(ptatinctx,&energy);CHKERRQ(ierr);
    ierr = pTatinPhysCompGetData_Energy(ptatinctx,&temperature,NULL);CHKERRQ(ierr);
    ierr = pTatin3d_ModelOutput_Temperature_Energy(ptatinctx,temperature,prefix);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationOutput_DefaultMeshFieldsLite"
PetscErrorCode FoundationOutput_DefaultMeshFieldsLite(pTatinCtx ptatinctx,Foundation f,Vec X,const char prefix[])
{
  PetscErrorCode   ierr;
  PetscBool        active_energy;

  PetscFunctionBegin;

  /* Light weight viewer: Only v is written. v and coords are expressed as floats */
  ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);

	/* T written out */
	ierr = pTatinContextValid_Energy(ptatinctx,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
    PhysCompEnergy energy;
    Vec            temperature;
    
    ierr = pTatinGetContext_Energy(ptatinctx,&energy);CHKERRQ(ierr);
    ierr = pTatinPhysCompGetData_Energy(ptatinctx,&temperature,NULL);CHKERRQ(ierr);
    ierr = pTatin3d_ModelOutput_Temperature_Energy(ptatinctx,temperature,prefix);CHKERRQ(ierr);
  }
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

  for (i=0; i<f->output->n_types; i++) {
    switch (f->output->type[i]) {
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
  for (i=0; i<f->output->n_userfuncs; i++) {
    ierr = f->output->userfunc[i](f,X,prefix);CHKERRQ(ierr);
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
  f->n_types = 0;
  f->n_userfuncs = 0;
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
  cJSON          *jroot,*typelist,*userlist,*jobj_k;
  int            entry,nmodes;
  PetscBool      match,user_found = PETSC_FALSE;
  PetscErrorCode ierr;
  
  /* look for MeshGeomICType */
  ierr = FoundationParseJSONGetItemEssential(root,"ModelOutput",&jroot);CHKERRQ(ierr);

  ierr = FoundationParseJSONGetItemEssential(jroot,"TypeList",&typelist);CHKERRQ(ierr);

  nmodes = cJSON_GetArraySize(typelist);
  PetscPrintf(PETSC_COMM_WORLD,"[foundation]--------- Found %d output types\n",nmodes);
  f->output->n_types = (PetscInt)nmodes;
  
  jobj_k = cJSON_GetArrayItemRoot(typelist);
  for (entry=0; entry<nmodes; entry++) {

    PetscStrcmp(jobj_k->valuestring,"MeshFields_Default",&match);
    if (match) { f->output->type[entry] = FND_OutputType_MeshDefault; }

    PetscStrcmp(jobj_k->valuestring,"MeshFieldsLite_Default",&match);
    if (match) { f->output->type[entry] = FND_OutputType_MeshLiteDefault; }

    PetscStrcmp(jobj_k->valuestring,"MaterialPoints_Default",&match);
    if (match) { f->output->type[entry] = FND_OutputType_MaterialPointDefault; }

    PetscStrcmp(jobj_k->valuestring,"User",&match);
    if (match) {
      f->output->type[entry] = FND_OutputType_User;
      user_found = PETSC_TRUE;
    }

    jobj_k = cJSON_GetArrayItemNext(jobj_k);
  }

  if (user_found) {
    ierr = FoundationParseJSONGetItemEssential(jroot,"UserFunctionList",&userlist);CHKERRQ(ierr);
  
    nmodes = cJSON_GetArraySize(userlist);
    PetscPrintf(PETSC_COMM_WORLD,"[foundation]--------- Found %d user functions\n",nmodes);
    f->output->n_userfuncs = (PetscInt)nmodes;
    
    jobj_k = cJSON_GetArrayItemRoot(userlist);
    for (entry=0; entry<nmodes; entry++) {
      char                 *functionname;
      FoundationUserOutput fp;
      
      functionname = jobj_k->valuestring;
      
      ierr = PetscFunctionListFind(f->user->flist_output,(const char*)functionname,(void(**)(void))&fp);CHKERRQ(ierr);
      if (!fp) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserOutput: Method with name \"%s\" has not been registered",functionname);
      f->output->userfunc[entry] = fp;
      
      jobj_k = cJSON_GetArrayItemNext(jobj_k);
    }
  }
  
  
  PetscFunctionReturn(0);
}

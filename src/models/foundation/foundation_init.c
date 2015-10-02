
#include "petsc.h"
#include "ptatin_utils.h"
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"
#include "foundation_output.h"


#undef __FUNCT__
#define __FUNCT__ "FoundationGetPtatinCtx"
PetscErrorCode FoundationGetPtatinCtx(Foundation f,pTatinCtx *ctx)
{
  if (ctx) {
    *ctx = f->ptatin_ctx;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationGetUserVariables"
PetscErrorCode FoundationGetUserVariables(Foundation f,FoundationUserVars *fv)
{
  if (fv) {
    *fv = f->user;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseJSONGetItemEssential"
PetscErrorCode FoundationParseJSONGetItemEssential(cJSON *jroot,const char itemname[],cJSON **jitem)
{
  PetscPrintf(PETSC_COMM_WORLD,"[foundation] Parsing contents of \"%s\"\n",itemname);
  *jitem = NULL;
  *jitem = cJSON_GetObjectItem(jroot,itemname);
  if (!(*jitem)) {
    SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"Foundation: Failed to locate essential struct \"%s\".\"%s\" ",jroot->string,itemname);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseJSONGetItemOptional"
PetscErrorCode FoundationParseJSONGetItemOptional(cJSON *jroot,const char itemname[],cJSON **jitem)
{
  PetscPrintf(PETSC_COMM_WORLD,"[foundation] Parsing contents of \"%s\" <optional...",itemname);
  *jitem = NULL;
  *jitem = cJSON_GetObjectItem(jroot,itemname);
  if (!*jitem) {
    PetscPrintf(PETSC_COMM_WORLD,"not found>\n",itemname);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"found>\n",itemname);
  }
  PetscFunctionReturn(0);
}

/* 
 Looks for particular keywords within jroot and prints the string
*/
#undef __FUNCT__
#define __FUNCT__ "FoundationParseJSONReportInfo"
PetscErrorCode FoundationParseJSONReportInfo(cJSON *jroot)
{
  cJSON *jitem;
  const char *comment_list[] = { "description", "comment", "doc", "note", "help", 0 };
  int k;
  
  for (k=0; k<5; k++) {
    jitem = NULL;
    FoundationParseJSONGetItemOptional(jroot,comment_list[k],&jitem);
    if (jitem) {
      PetscPrintf(PETSC_COMM_WORLD,"FoundationObject[%s].%s %s\n",jroot->string,comment_list[k],jitem->valuestring);
    }
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseJSONGetItemFromListEssential"
PetscErrorCode FoundationParseJSONGetItemFromListEssential(cJSON *jroot,const char itemname[],cJSON **jitem)
{
  int entry,nlist;
  cJSON *dummy,*jobj_k;
  PetscErrorCode ierr;
  
  *jitem = NULL;

  if (!cJSON_IsArrayItem(jroot)) {
    ierr = FoundationParseJSONGetItemEssential(jroot,itemname,jitem);CHKERRQ(ierr);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"[foundation] Parsing list contents within \"%s\"\n",jroot->string);

    nlist = cJSON_GetArraySize(jroot);
    jobj_k = cJSON_GetArrayItemRoot(jroot);
    for (entry=0; entry<nlist; entry++) {
      dummy = cJSON_GetObjectItem(jobj_k,itemname);
      if (dummy) { break; }
      jobj_k = cJSON_GetArrayItemNext(jobj_k);
    }
    
    *jitem = dummy;

    if (!(*jitem)) {
      SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_USER,"Foundation: Failed to locate essential struct \"%s\":[{\"%s\":{}}]",jroot->string,itemname);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseJSONGetItemFromListOptional"
PetscErrorCode FoundationParseJSONGetItemFromListOptional(cJSON *jroot,const char itemname[],cJSON **jitem)
{
  int entry,nlist;
  cJSON *dummy,*jobj_k;
  PetscErrorCode ierr;
  
  *jitem = NULL;
  
  if (!cJSON_IsArrayItem(jroot)) {
    ierr = FoundationParseJSONGetItemOptional(jroot,itemname,jitem);CHKERRQ(ierr);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"[foundation] Parsing list contents within \"%s\" <optional...",jroot->string);
    
    nlist = cJSON_GetArraySize(jroot);
    jobj_k = cJSON_GetArrayItemRoot(jroot);
    for (entry=0; entry<nlist; entry++) {
      dummy = cJSON_GetObjectItem(jobj_k,itemname);
      if (dummy) { break; }
      jobj_k = cJSON_GetArrayItemNext(jobj_k);
    }
    *jitem = dummy;

    if (!*jitem) {
      PetscPrintf(PETSC_COMM_WORLD,"not found>\n",itemname);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,"found>\n",itemname);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationJSONFindItemFromList"
PetscErrorCode FoundationJSONFindItemFromList(cJSON *jroot,const char *methodlist[],PetscInt *_index,cJSON **_jitem)
{
  cJSON *jitem;
  PetscInt k,index;
  char **method;
  
  index = 0;
  method = (char**)&methodlist[0];
  while (*method) {
    jitem = cJSON_GetObjectItem(jroot,*method);
    if (jitem) { break; }
    index++;
    method++;
  }

  if (!jitem) {
    PetscPrintf(PETSC_COMM_WORLD,"[foundation]:%s A valid constructor was not found. Valid constructors include\n",jroot->string);
    method = (char**)&methodlist[0];
    while (*method) {
      PetscPrintf(PETSC_COMM_WORLD,"  \"%s\"\n",*method);
      method++;
    }
    SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_SUP,"[foundation]:%s no support for method \"%s\"",jroot->string,jroot->child->string);
  }

  *_index = index;
  *_jitem = jitem;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationFindStringInList"
PetscErrorCode FoundationFindStringInList(const char *list[],const char itemname[],PetscInt *index)
{
  char **ii;
  PetscBool match;
  PetscInt idx;

  *index = -1;

  idx = 0;
  ii = (char**)&list[0];
  while (*ii) {
    PetscStrcmp(itemname,*ii,&match);
    if (match) break;
    ii++;
    idx++;
  }
  if (match) *index = idx;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Foundation"
PetscErrorCode ModelInitialize_Foundation(pTatinCtx c,void *ctx)
{
  Foundation     data = (Foundation)ctx;
  PetscErrorCode ierr;
  cJSON          *jfile;
  char           filename[PETSC_MAX_PATH_LEN];
  char           *pstr;
  PetscBool      flg;
  cJSON          *jobjectroot = NULL;
  int            found;
  
  /* defaults */
  data->ptatin_ctx           = c;
  data->mesh_geom_type       = FND_MeshGeomT_NULL;
  data->mesh_ref_type        = FND_MeshRefT_NULL;
  data->material_region_type = FND_MatRegT_NULL;
  
	flg = PETSC_FALSE;
	ierr = PetscOptionsGetString(NULL,"-fndin",filename,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Foundation: Must specify an input file via -fndin");
  }
  
  /* read json file - look for output path */
  cJSON_FileView(filename,&jfile);
  if (!jfile) {
    SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Failed to open JSON file: %s",filename);
  }
  
  /*
  jobjectroot = cJSON_GetObjectItem(jfile,"pTatinFoundation");
  if (!jobjectroot) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Foundation: Failed to locate struct \"pTatinFoundation\" in JSON file");
  }
  */
  ierr = FoundationParseJSONGetItemEssential(jfile,"pTatinFoundation",&jobjectroot);CHKERRQ(ierr);
  ierr = FoundationParseJSONReportInfo(jobjectroot);CHKERRQ(ierr);
  data->jfile = jfile;
  data->root = jobjectroot;
  
  PetscPrintf(PETSC_COMM_WORLD,"[foundation] Parsing contents of \"pTatinFoundation\"\n");
  
  cJSON_GetObjectValue_char(jobjectroot,"OutputPath",&found,&pstr);
  if (found) {
    PetscSNPrintf(c->outputpath,PETSC_MAX_PATH_LEN-1,pstr);
  } else {
    PetscSNPrintf(c->outputpath,PETSC_MAX_PATH_LEN-1,"./output");
  }
	ierr = pTatinCreateDirectory(c->outputpath);CHKERRQ(ierr);

  
  ierr = FoundationParseMaterialMetaData(c,data);CHKERRQ(ierr);
  {
    RheologyConstants *rheology;

    ierr = pTatinGetRheology(c,&rheology);CHKERRQ(ierr);
    rheology->rheology_type = RHEOLOGY_VP_STD;
    
    rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
    rheology->eta_upper_cutoff_global = 1.0e10;
    rheology->eta_lower_cutoff_global = 1.0e-10;
  }
  
  ierr = FoundationOutputParse(data,data->root);CHKERRQ(ierr);
  
  {
    cJSON *juservarroot;
    ierr = FoundationParseJSONGetItemOptional(jobjectroot,"FoundationUserVariables",&juservarroot);CHKERRQ(ierr);
    if (juservarroot) {
      ierr = FoundationUserVariablesParse(data,juservarroot);CHKERRQ(ierr);
    }
  }

  /*
  {
    cJSON *juservarroot;
    PetscReal coord[] = {1.01, 2.02, 3.3},value;
    
    ierr = FoundationUserApply_Evaluate(data,data->user,"UserEvaluator_Empty",coord,&value);CHKERRQ(ierr);
  }
  */
   
  PetscFunctionReturn(0);
}

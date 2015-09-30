
#include "petsc.h"
#include "ptatin_utils.h"
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"


#undef __FUNCT__
#define __FUNCT__ "FoundationParseJSONGetItemEssential"
PetscErrorCode FoundationParseJSONGetItemEssential(cJSON *jroot,const char itemname[],cJSON **jitem)
{
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
  *jitem = NULL;
  *jitem = cJSON_GetObjectItem(jroot,itemname);
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
  
  
  cJSON_GetObjectValue_char(jobjectroot,"OutputPath",&found,&pstr);
  if (found) {
    PetscSNPrintf(c->outputpath,PETSC_MAX_PATH_LEN-1,pstr);
  } else {
    PetscSNPrintf(c->outputpath,PETSC_MAX_PATH_LEN-1,"./output");
  }
	ierr = pTatinCreateDirectory(c->outputpath);CHKERRQ(ierr);

  
  ierr = FoundationParseMaterialMetaData(c,data);CHKERRQ(ierr);

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
    
    ierr = FoundationUserEvaluate(data,data->user,"UserEvaluator_Empty",coord,&value);CHKERRQ(ierr);
  }
  */
   
  PetscFunctionReturn(0);
}

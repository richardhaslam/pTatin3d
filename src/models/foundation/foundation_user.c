
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"


/*
 
{
"FoundationUserVariables": [
  {
    "datatype": "double",
    "a": 11.1,
    "b": 12.3,
    "v": 14.3,
    "gggg": 15.1
  },
  {
    "datatype": "double",
    "aaaa": 14.1
  },
  {
    "datatype": "int",
    "a": 11
  }
]}
 
*/


#undef __FUNCT__
#define __FUNCT__ "UserEvaluator_Empty"
PetscErrorCode UserEvaluator_Empty(pTatinCtx ctx,FoundationUserVars user,PetscReal coor[],PetscReal *value)
{
  PetscReal *userdata;
  PetscErrorCode ierr;
  
  PetscPrintf(PETSC_COMM_WORLD,"*** Evaluator: Empty\n");
  PetscPrintf(PETSC_COMM_WORLD,"***  ptatin    %p\n",ctx);
  PetscPrintf(PETSC_COMM_WORLD,"***  user vars %p\n",user);
  PetscPrintf(PETSC_COMM_WORLD,"***  coor      [%+1.4e, %+1.4e, %+1.4e]\n",coor[0],coor[1],coor[2]);
  
  ierr = UserVariablesGetDoubleFromName(user,"var1",&userdata);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"***  var1      %+1.4e\n",*userdata);
  
  *value = 0.0;
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FoundationUserVarsCreate"
PetscErrorCode FoundationUserVarsCreate(FoundationUserVars *uv)
{
  FoundationUserVars v;
  
  PetscMalloc(sizeof(struct _p_FoundationUserVars),&v);
  PetscMemzero(v,sizeof(struct _p_FoundationUserVars));
  *uv = v;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationUserVarsDestroy"
PetscErrorCode FoundationUserVarsDestroy(FoundationUserVars *uv)
{
  FoundationUserVars v;
  PetscErrorCode ierr;

  if (!uv) PetscFunctionReturn(0);
  v = *uv;
  ierr = PetscFunctionListDestroy(&v->flist_evaluator);CHKERRQ(ierr);
  PetscFree(v);
  *uv = NULL;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UserVariablesGetIndexFromName"
PetscErrorCode UserVariablesGetIndexFromName(FoundationUserVars uv,const char name[],UserVarDataTypes type,int *ii)
{
  int  nlist,i;
  char **list;
  
  *ii = -1;
  switch (type) {
    case UVDTYPE_INT:
      nlist = uv->vi_nreg;
      list = uv->var_names_i;
      break;
    case UVDTYPE_FLOAT:
      nlist = uv->vf_nreg;
      list = uv->var_names_f;
      break;
    case UVDTYPE_DOUBLE:
      nlist = uv->vd_nreg;
      list = uv->var_names_d;
      break;
  }

  for (i=0; i<nlist; i++) {
    if (strcmp(name,list[i]) == 0) { *ii = i; break; }
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UserVariablesAdd_Int"
PetscErrorCode UserVariablesAdd_Int(FoundationUserVars uv,const char name[],int value)
{
  int idx;
  PetscErrorCode ierr;

  if (uv->vi_space_used == FND_MAX_USER_SPACE) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable(int): Statically allocated space is exhausted");
  if (uv->vi_nreg == FND_MAX_USER_VARS) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable(int): Statically allocated registers are exhausted");

  ierr = UserVariablesGetIndexFromName(uv,name,UVDTYPE_INT,&idx);CHKERRQ(ierr);
  if (idx != -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable.\"%s\"(int) already registered\n",name);
  asprintf(&uv->var_names_i[uv->vi_nreg],"%s",name);
  uv->vars_i[uv->vi_nreg] = value;
  uv->vi_space_used++;
  uv->vi_nreg++;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UserVariablesAdd_Float"
PetscErrorCode UserVariablesAdd_Float(FoundationUserVars uv,const char name[],float value)
{
  int idx;
  PetscErrorCode ierr;
  
  if (uv->vf_space_used == FND_MAX_USER_SPACE) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable(float): Statically allocated space is exhausted");
  if (uv->vf_nreg == FND_MAX_USER_VARS) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable(float): Statically allocated registers are exhausted");
  
  ierr = UserVariablesGetIndexFromName(uv,name,UVDTYPE_FLOAT,&idx);CHKERRQ(ierr);
  if (idx != -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable.\"%s\"(float) already registered\n",name);
  asprintf(&uv->var_names_f[uv->vf_nreg],"%s",name);
  uv->vars_f[uv->vf_nreg] = value;
  uv->vf_space_used++;
  uv->vf_nreg++;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UserVariablesAdd_Double"
PetscErrorCode UserVariablesAdd_Double(FoundationUserVars uv,const char name[],double value)
{
  int idx;
  PetscErrorCode ierr;
  
  if (uv->vd_space_used == FND_MAX_USER_SPACE) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable(double): Statically allocated space is exhausted");
  if (uv->vd_nreg == FND_MAX_USER_VARS) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable(double): Statically allocated registers are exhausted");

  ierr = UserVariablesGetIndexFromName(uv,name,UVDTYPE_DOUBLE,&idx);CHKERRQ(ierr);
  if (idx != -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable.\"%s\"(double) already registered\n",name);
  asprintf(&uv->var_names_d[uv->vd_nreg],"%s",name);
  uv->vars_d[uv->vd_nreg] = value;
  uv->vd_space_used++;
  uv->vd_nreg++;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UserVariablesGetIntFromName"
PetscErrorCode UserVariablesGetIntFromName(FoundationUserVars uv,const char name[],int **value)
{
  int idx;
  PetscErrorCode ierr;
  
  ierr = UserVariablesGetIndexFromName(uv,name,UVDTYPE_INT,&idx);CHKERRQ(ierr);
  if (idx == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"UserVariable: Variable \"%s\"(int) was not registered",name);
  (*value) = &uv->vars_i[idx];
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UserVariablesGetFloatFromName"
PetscErrorCode UserVariablesGetFloatFromName(FoundationUserVars uv,const char name[],float **value)
{
  int idx;
  PetscErrorCode ierr;
  
  ierr = UserVariablesGetIndexFromName(uv,name,UVDTYPE_FLOAT,&idx);CHKERRQ(ierr);
  if (idx == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"UserVariable: Variable \"%s\"(float) was not registered",name);
  (*value) = &uv->vars_f[idx];
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UserVariablesGetDoubleFromName"
PetscErrorCode UserVariablesGetDoubleFromName(FoundationUserVars uv,const char name[],double **value)
{
  int idx;
  PetscErrorCode ierr;

  ierr = UserVariablesGetIndexFromName(uv,name,UVDTYPE_DOUBLE,&idx);CHKERRQ(ierr);
  if (idx == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"UserVariable: Variable \"%s\"(double) was not registered",name);
  (*value) = &uv->vars_d[idx];
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationUserEvaluate"
PetscErrorCode FoundationUserEvaluate(Foundation f,FoundationUserVars uv,const char functionname[],PetscReal *coor,PetscReal *value)
{
  FoundationUserEvaluator fp;
  pTatinCtx ptatin;
  PetscErrorCode ierr;
  
  ierr = PetscFunctionListFind(uv->flist_evaluator,(const char*)functionname,(void(**)(void))&fp);CHKERRQ(ierr);
  if (!fp) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Evaluator: Method with name \"%s\" has not been registered",functionname);
  ptatin = f->ptatin_ctx;
  ierr = fp(ptatin,uv,coor,value);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationUserRegisterEvaluator"
PetscErrorCode FoundationUserRegisterEvaluator(Foundation f,FoundationUserVars uv,const char functionname[],FoundationUserEvaluator fp)
{
  PetscErrorCode ierr;
  void (*eval)(void);
  
  ierr = PetscFunctionListFind(uv->flist_evaluator,(const char*)functionname,&eval);CHKERRQ(ierr);
  if (eval) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Evaluator: Method with name \"%s\" already registered",functionname);
  
  ierr = PetscFunctionListAdd(&uv->flist_evaluator,(const char*)functionname,(void(*)(void))fp);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FoundationUserVariablesParse"
PetscErrorCode FoundationUserVariablesParse(Foundation f,cJSON *root)
{
  cJSON *juservar,*jobj_k,*field;
  int   k,nvars,found;
  char *dt;
  UserVarDataTypes type;
  PetscErrorCode ierr;
  FoundationUserVars uv;
  
  uv = f->user;
  
  nvars = cJSON_GetArraySize(root);
  jobj_k = cJSON_GetArrayItemRoot(root);
  for (k=0; k<nvars; k++) {

    ierr = FoundationParseJSONGetItemEssential(jobj_k,"datatype",&juservar);CHKERRQ(ierr);
    cJSON_GetObjectValue_char(jobj_k,"datatype",&found,&dt);
    type = UVDTYPE_UNDEF;
    if (strcmp(dt,"int")    == 0) { type = UVDTYPE_INT; }
    if (strcmp(dt,"float")  == 0) { type = UVDTYPE_FLOAT; }
    if (strcmp(dt,"double") == 0) { type = UVDTYPE_DOUBLE; }
    if (type == UVDTYPE_UNDEF) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable: DataType \"%s\" is unknown",dt);
    
    field = jobj_k->child;
    if (field == NULL) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"FoundationUserVariable: Expected an array of variables");
    if (field == juservar) field = field->next;
    
    while (field != NULL) { /* if not datatype, check for next */
      if (field == juservar) {
        field = field->next;
        continue;
      } else {
        switch (type) {
          case UVDTYPE_INT:
            printf("name: %s : val %d\n",field->string,field->valueint);
            ierr = UserVariablesAdd_Int(uv,field->string,field->valueint);CHKERRQ(ierr);
            break;
          case UVDTYPE_FLOAT:
            printf("name: %s : val %f\n",field->string,(float)field->valuedouble);
            ierr = UserVariablesAdd_Float(uv,field->string,(float)field->valuedouble);CHKERRQ(ierr);
            break;
          case UVDTYPE_DOUBLE:
            printf("name: %s : val %1.4e\n",field->string,field->valuedouble);
            ierr = UserVariablesAdd_Double(uv,field->string,field->valuedouble);CHKERRQ(ierr);
            break;
        }
        field = field->next;
      }
    }
    
    jobj_k = cJSON_GetArrayItemNext(jobj_k);
  }
  
  /*
  PetscFunctionListAdd(&uv->plist_material_ic,"t1",Test1);
  PetscFunctionListAdd(&uv->plist_material_ic,"t2",Test2);
  
  {
    PetscErrorCode (*f1)(Foundation,FoundationUserVars);

    //PetscFunctionListFind(uv->plist_material_ic,"t1",&f1);
    PetscFunctionListFind(uv->plist_material_ic,"t2",&f1);
    ierr = f1(f,uv);CHKERRQ(ierr);
  }
  */
  ierr = FoundationUserRegisterEvaluator(f,uv,"UserEvaluator_Empty",UserEvaluator_Empty);CHKERRQ(ierr);
  

  PetscFunctionReturn(0);
}





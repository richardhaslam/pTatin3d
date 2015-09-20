
#ifndef __foundation_user_h__
#define __foundation_user_h__

#include "foundation.h"


#define FND_MAX_USER_SPACE 50
#define FND_MAX_USER_VARS 50


struct _p_FoundationUserVars {
  int    vi_nreg,vf_nreg,vd_nreg;
  int    vi_space_used,vf_space_used,vd_space_used;
  int    vars_i[FND_MAX_USER_SPACE];
  float  vars_f[FND_MAX_USER_SPACE];
  double vars_d[FND_MAX_USER_SPACE];
  char *var_names_i[FND_MAX_USER_VARS];
  char *var_names_f[FND_MAX_USER_VARS];
  char *var_names_d[FND_MAX_USER_VARS];
  PetscFunctionList flist_evaluator;
};

typedef PetscErrorCode (*FoundationUserEvaluator)(pTatinCtx,FoundationUserVars,PetscReal*,PetscReal*);

typedef enum { UVDTYPE_INT=0, UVDTYPE_FLOAT, UVDTYPE_DOUBLE, UVDTYPE_UNDEF } UserVarDataTypes;

PetscErrorCode FoundationUserVarsCreate(FoundationUserVars *uv);
PetscErrorCode FoundationUserVarsDestroy(FoundationUserVars *uv);

PetscErrorCode FoundationUserVariablesParse(Foundation f,cJSON *root);
PetscErrorCode FoundationUserRegisterEvaluator(FoundationUserVars uv,const char functioname[],FoundationUserEvaluator fp);
PetscErrorCode FoundationUserEvaluate(Foundation f,FoundationUserVars uv,const char functioname[],PetscReal *coor,PetscReal *value);

PetscErrorCode UserVariablesGetIntFromName(FoundationUserVars uv,const char name[],int **value);
PetscErrorCode UserVariablesGetFloatFromName(FoundationUserVars uv,const char name[],float **value);
PetscErrorCode UserVariablesGetDoubleFromName(FoundationUserVars uv,const char name[],double **value);

PetscErrorCode UserVariablesAdd_Int(FoundationUserVars uv,const char name[],int value);
PetscErrorCode UserVariablesAdd_Float(FoundationUserVars uv,const char name[],float value);
PetscErrorCode UserVariablesAdd_Double(FoundationUserVars uv,const char name[],double value);

#endif

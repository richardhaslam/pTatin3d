

#ifndef __foundation_output_h__
#define __foundation_output_h__

#include "foundation_user.h"

struct _p_FoundationOutput {
  PetscInt             n_types;
  FND_OutputType       type[100];
  PetscInt             n_userfuncs;
  FoundationUserOutput userfunc[100];
};

PetscErrorCode FoundationOutputCreate(FoundationOutput *fo);
PetscErrorCode FoundationOutputDestroy(FoundationOutput *fo);
PetscErrorCode FoundationOutputParse(Foundation f,cJSON *root);

#endif

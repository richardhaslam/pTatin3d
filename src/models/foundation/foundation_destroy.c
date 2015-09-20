
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Foundation"
PetscErrorCode ModelDestroy_Foundation(pTatinCtx c,void *ctx)
{
  Foundation data = (Foundation)ctx;
  PetscErrorCode ierr;
  
  cJSON_Delete(data->jfile);
  ierr = FoundationUserVarsDestroy(&data->user);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

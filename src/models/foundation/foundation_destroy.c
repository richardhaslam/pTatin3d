
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"
#include "foundation_output.h"

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Foundation"
PetscErrorCode ModelDestroy_Foundation(pTatinCtx c,void *ctx)
{
  Foundation data = (Foundation)ctx;
  PetscErrorCode ierr;
  
  cJSON_Delete(data->jfile);
  ierr = FoundationUserVarsDestroy(&data->user);CHKERRQ(ierr);
  ierr = FoundationOutputDestroy(&data->output);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

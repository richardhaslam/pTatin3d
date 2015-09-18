
#include "foundation.h"
#include "foundation_impl.h"

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Foundation"
PetscErrorCode ModelDestroy_Foundation(pTatinCtx c,void *ctx)
{
  Foundation data = (Foundation)ctx;
  
  cJSON_Delete(data->jfile);

  PetscFunctionReturn(0);
}


#include "petsc.h"
#include "ptatin3d.h"
#include "foundation.h"
#include "foundation_user.h"

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStandardViewer"
PetscErrorCode MaterialPointStandardViewer(Foundation f,Vec X,const char prefix[])
{
  pTatinCtx ptatinctx;
  PetscErrorCode ierr;
  
  ierr = FoundationGetPtatinCtx(f,&ptatinctx);CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_MPntStd(ptatinctx,prefix);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}
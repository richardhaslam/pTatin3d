
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"

/* jam methods in here */
#include "userevaluator_empty.c"



#undef __FUNCT__
#define __FUNCT__ "FoundationUserRegisterFunctions"
PetscErrorCode FoundationUserRegisterFunctions(Foundation f)
{
  PetscErrorCode ierr;

  /* evaluators */
  ierr = FoundationUserRegisterEvaluator(f->user,"UserEvaluator_Empty",UserEvaluator_Empty);CHKERRQ(ierr);

  /* vel dirichlet bc evaluators */

  /* stokes neumann bc evaluators */
  
  /* temperature dirichlet bc evaluators */

  /* vel-pressure ic */

  PetscFunctionReturn(0);
}

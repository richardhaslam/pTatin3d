/*
 Example of user registered function evaluator
 Illustrates usage of user variable dynamic registration / access via the function
 UserVariablesGetDoubleFromName()

*/

#include "petsc.h"
#include "ptatin3d.h"
#include "foundation.h"
#include "foundation_user.h"


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


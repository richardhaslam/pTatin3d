
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"
#include "foundation_output.h"

const char *FND_MeshGeomTypeNames[] = { "Cartesian", "Hex", "DMDAFile", "Lua", 0 };
const char *FND_MeshRefinementTypeName[] = { "A", "B", "Lua", 0 };
const char *FND_MaterialRegionTypeNames[] = { "JSON", "DBFile", "Lua", 0 };
const char *FND_BoundaryConditionTypeNames[] = { "Dirichlet", "Neumann", 0 };
const char *FND_BoundaryConditionModeNames[] = { "Scalar", "Component", "Vector", "Tensor", "User", "Lua", 0 };
const char *FND_OutputTypeNames[] = { "MeshFieldsDefault", "MeshFieldsLiteDefault", "MaterialPointsDefault", "User", 0 };

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Foundation"
PetscErrorCode pTatinModelRegister_Foundation(void)
{
  Foundation       data;
  pTatinModel      m;
  PetscErrorCode   ierr;
  
  PetscFunctionBegin;
  /* Allocate memory for the data structure for this model */
  ierr = PetscMalloc(sizeof(struct _p_Foundation),&data);CHKERRQ(ierr);
  ierr = PetscMemzero(data,sizeof(struct _p_Foundation));CHKERRQ(ierr);
  
  ierr = FoundationUserVarsCreate(&data->user);CHKERRQ(ierr);
  ierr = FoundationOutputCreate(&data->output);CHKERRQ(ierr);

  ierr = FoundationUserRegisterFunctions(data);CHKERRQ(ierr);
  
  /* register user model */
  ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
  
  /* Set name, model select via -ptatin_model foundation */
  ierr = pTatinModelSetName(m,"foundation");CHKERRQ(ierr);
  
  /* Set model data */
  ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
  
  /* Set function pointers */
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Foundation);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelInitialMeshGeometry_Foundation);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelInitialMaterialGeometry_Foundation);CHKERRQ(ierr);

	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,               (void (*)(void))ModelInitialConditionMeshFields_Foundation);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,(void (*)(void))ModelInitialConditonMarkerFields_Foundation);CHKERRQ(ierr);
  
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBC_Foundation);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBCMG_Foundation);CHKERRQ(ierr);

  //ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelUpdateMeshGeometry_Foundation);CHKERRQ(ierr);
  //ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelMaterialBoundaryCondition_Foundation);CHKERRQ(ierr);

  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Foundation);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Foundation);CHKERRQ(ierr);
  
  /* Insert model into list */
  ierr = pTatinModelRegister(m);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

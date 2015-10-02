
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"

#include "dmda_bcs.h"
#include "ptatin3d_energy.h"
#include "rheology.h"
#include "dmda_iterator.h"







#undef __FUNCT__
#define __FUNCT__ "FoundationMarkerFieldsIC_FromRheology"
PetscErrorCode FoundationMarkerFieldsIC_FromRheology(pTatinCtx c,Foundation f,Vec X,PetscBool ignore_pressure_dependence)
{
  PetscErrorCode  ierr;
	DM              dmstokes,dmu,dmp;
	PhysCompStokes  stokes;
	Vec             Uloc,Ploc;
	PetscScalar     *LA_Uloc,*LA_Ploc;
  DataField       PField;
  MaterialConst_MaterialType *mattypes;
  int                        nregions;
  DataBucket                 materialconstants;
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  DataBucketGetSizes(materialconstants,&nregions,0,0);
  DataBucketGetDataFieldByName(materialconstants,MaterialConst_MaterialType_classname,&PField);
  DataFieldGetEntries(PField,(void**)&mattypes);

	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
	
	ierr = DMCompositeGetEntries(dmstokes,&dmu,&dmp);CHKERRQ(ierr);
	ierr = DMCompositeGetLocalVectors(dmstokes,&Uloc,&Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(dmstokes,X,Uloc,Ploc);CHKERRQ(ierr);

  if (ignore_pressure_dependence) {
    ierr = VecZeroEntries(Ploc);CHKERRQ(ierr);
  }
  
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	ierr = pTatin_EvaluateRheologyNonlinearities(c,dmu,LA_Uloc,dmp,LA_Ploc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreLocalVectors(dmstokes,&Uloc,&Ploc);CHKERRQ(ierr);
  DataFieldRestoreEntries(PField,(void**)&mattypes);
  
  PetscFunctionReturn(0);
}

/* Leaving this out - users can do it themselves via a plugin */
#undef __FUNCT__
#define __FUNCT__ "Foundation_LinearInterpolateParse"
PetscErrorCode Foundation_LinearInterpolateParse(Foundation f,cJSON *root,DM dm,Vec X,PetscInt dof_idx)
{
  PetscErrorCode ierr;
	DMDAVecTraverse3d_InterpCtx IntpCtx;
  double value_0,value_1;
  PetscReal delta_v,delta_x;
  int dir;
	PetscReal MeshMin[3],MeshMax[3];
  int found;
  
	ierr = DMDAGetBoundingBox(dm,MeshMin,MeshMax);CHKERRQ(ierr);

  value_1 = 0.0;
  value_0 = 0.0;

  cJSON_GetObjectValue_double(root,"value0",&found,&value_0);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"LinearInterpolation.value0 missing");

  cJSON_GetObjectValue_double(root,"value1",&found,&value_1);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"LinearInterpolation.value1 missing");

  cJSON_GetObjectValue_int(root,"direction",&found,&dir);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"LinearInterpolation.dir missing");
  
  delta_v = (PetscReal)(value_1 - value_0);
  switch (dir) {
    case 0:
      delta_x = MeshMax[0] - MeshMin[0];
      ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,delta_v/delta_x,value_0,0.0);CHKERRQ(ierr);
      break;
    case 1:
      delta_x = MeshMax[1] - MeshMin[1];
      ierr = DMDAVecTraverse3d_InterpCtxSetUp_Y(&IntpCtx,delta_v/delta_x,value_0,0.0);CHKERRQ(ierr);
      break;
    case 2:
      delta_x = MeshMax[2] - MeshMin[2];
      ierr = DMDAVecTraverse3d_InterpCtxSetUp_Z(&IntpCtx,delta_v/delta_x,value_0,0.0);CHKERRQ(ierr);
      break;
  }
	ierr = DMDAVecTraverse3d(dm,X,dof_idx,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);

  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationVelocityICParse"
PetscErrorCode FoundationVelocityICParse(Foundation f,cJSON *root,DM dmu,Vec Xu)
{
  const char *methodlist[] = { "None", "LinearInterpolation", "User", 0 };
  cJSON *jitem;
  PetscInt type;
  PetscErrorCode ierr;
  
  /* Look for vx, vy, vz */
  /* Look for LinearInterpolation */
  ierr = FoundationJSONFindItemFromList(root,methodlist,&type,&jitem);CHKERRQ(ierr);
  
  switch (type) {
    case 0:
      break;

    case 1:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"LinearInterpolation is unsupported");
      break;

    case 2:
      break;

    default:
      break;
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationPressureICParse"
PetscErrorCode FoundationPressureICParse(Foundation f,cJSON *root,DM dmstokes,Vec X)
{
  
  
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationTemperatureICParse"
PetscErrorCode FoundationTemperatureICParse(Foundation f,cJSON *root,DM dmT,Vec X)
{
  
  
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialConditionMeshFields_Foundation"
PetscErrorCode ModelInitialConditionMeshFields_Foundation(pTatinCtx c,Vec X,void *ctx)
{
  Foundation     f = (Foundation)ctx;
  cJSON          *root,*jmeshic,*jitem[3];
  int            size;
  PetscBool      active_energy;
	DM             dmstokes,dmu,dmp;
	PhysCompStokes stokes;
	Vec            Xu,Xp;
  PetscErrorCode ierr;
  
  /* look for MeshFieldsIC */
  root = f->root;
  ierr = FoundationParseJSONGetItemEssential(root,"MeshFieldsIC",&jmeshic);CHKERRQ(ierr);
  size = cJSON_GetArraySize(jmeshic);
  PetscPrintf(PETSC_COMM_WORLD,"[foundation]--------- Found %d fields\n",size);

  jitem[0] = jitem[1] = jitem[2] = NULL;
  ierr = FoundationParseJSONGetItemFromListOptional(jmeshic,"Velocity",&jitem[0]);CHKERRQ(ierr);
  ierr = FoundationParseJSONGetItemFromListOptional(jmeshic,"Pressure",&jitem[1]);CHKERRQ(ierr);

  /* Throw errors if temperature loaded an no intial condition is specified  */
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
    ierr = FoundationParseJSONGetItemFromListEssential(jmeshic,"Temperature",&jitem[2]);CHKERRQ(ierr);
  }
  
  /* parse */
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = DMCompositeGetEntries(dmstokes,&dmu,&dmp);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(dmstokes,X,&Xu,&Xp);CHKERRQ(ierr);

  /* parse optional u,p */
  if (jitem[0]) ierr = FoundationVelocityICParse(f,jitem[0],dmu,Xu);CHKERRQ(ierr);
  if (jitem[1]) ierr = FoundationPressureICParse(f,jitem[1],dmp,Xp);CHKERRQ(ierr);
	
  if (active_energy) {
    PhysCompEnergy energy;
    Vec            temperature;
    DM             dmT;
    
    ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
    ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);
    ierr = pTatinPhysCompGetDM_Energy(c,&dmT);CHKERRQ(ierr);
    
    ierr = FoundationTemperatureICParse(f,jitem[2],dmT,temperature);CHKERRQ(ierr);
  }
  ierr = DMCompositeRestoreAccess(dmstokes,X,&Xu,&Xp);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialConditonMarkerFields_Foundation"
PetscErrorCode ModelInitialConditonMarkerFields_Foundation(pTatinCtx c,Vec X,void *ctx)
{
  Foundation      f = (Foundation)ctx;
  
  PetscFunctionReturn(0);
}


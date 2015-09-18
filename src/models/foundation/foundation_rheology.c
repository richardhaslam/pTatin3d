
#include "foundation.h"
#include "foundation_impl.h"
#include "data_bucket.h"
#include "material_constants.h"


#undef __FUNCT__
#define __FUNCT__ "FoundationParse_None"
PetscErrorCode FoundationParse_None(DataBucket mc,int ridx,cJSON *root)
{
  PetscFunctionReturn(0);
}

/* FLOW LAWS */
#undef __FUNCT__
#define __FUNCT__ "FoundationParseViscousFlowLaw_Const"
PetscErrorCode FoundationParseViscousFlowLaw_Const(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"eta0",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.Const.eta0 missing");
  MaterialConstantsSetValues_ViscosityConst(mc,ridx,val_double);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseViscousFlowLaw_FrankKamentskii"
PetscErrorCode FoundationParseViscousFlowLaw_FrankKamentskii(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"eta0",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.FrankKamentskii.eta0 missing");
  MaterialConstantsSetValues_ViscosityFK(mc,ridx,val_double,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"theta",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.FrankKamentskii.theta missing");
  MaterialConstantsSetValues_ViscosityFK(mc,ridx,PETSC_DEFAULT,val_double);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseViscousFlowLaw_Arrhenius"
PetscErrorCode FoundationParseViscousFlowLaw_Arrhenius(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"preexpA",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.Arrhenius.preexpA missing");
  MaterialConstantsSetValues_ViscosityArrh(mc,ridx,val_double,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"Ascale",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.Arrhenius.Ascale missing");
  MaterialConstantsSetValues_ViscosityArrh(mc,ridx,PETSC_DEFAULT,val_double,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"entalpy",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.Arrhenius.entalpy missing");
  MaterialConstantsSetValues_ViscosityArrh(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,val_double,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"Vmol",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.Arrhenius.Vmol missing");
  MaterialConstantsSetValues_ViscosityArrh(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,val_double,PETSC_DEFAULT,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"nexp",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.Arrhenius.nexp missing");
  MaterialConstantsSetValues_ViscosityArrh(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,val_double,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"Tref",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.Arrhenius.Tref missing");
  MaterialConstantsSetValues_ViscosityArrh(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,val_double);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseViscousFlowLaw_DepthDep"
PetscErrorCode FoundationParseViscousFlowLaw_DepthDep(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"eta0",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.ViscosityZ.eta0 missing");
  MaterialConstantsSetValues_ViscosityZ(mc,ridx,val_double,PETSC_DEFAULT,PETSC_DEFAULT);
  
  cJSON_GetObjectValue_double(root,"zeta",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.ViscosityZ.zeta missing");
  MaterialConstantsSetValues_ViscosityZ(mc,ridx,PETSC_DEFAULT,val_double,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"zref",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"ViscousFlowLaw.ViscosityZ.zref missing");
  MaterialConstantsSetValues_ViscosityZ(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,val_double);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseViscousFlowLaw"
PetscErrorCode FoundationParseViscousFlowLaw(pTatinCtx c,int ridx,cJSON *root)
{
  PetscErrorCode ierr;
  cJSON *jobj;
  DataBucket materialconstants;
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);

  ierr = FoundationParseJSONGetItemOptional(root,"Const",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ Const\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,VISCOUS_CONSTANT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseViscousFlowLaw_Const(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = FoundationParseJSONGetItemOptional(root,"FrankKamentskii",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ FrankKamentskii\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,VISCOUS_FRANKK,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseViscousFlowLaw_FrankKamentskii(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  
  ierr = FoundationParseJSONGetItemOptional(root,"Arrhenius",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ Arrhenius\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,VISCOUS_ARRHENIUS,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseViscousFlowLaw_Arrhenius(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = FoundationParseJSONGetItemOptional(root,"Arrhenius2",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ Arrhenius2\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,VISCOUS_ARRHENIUS_2,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseViscousFlowLaw_Arrhenius(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = FoundationParseJSONGetItemOptional(root,"DepthDependent",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ DepthDependent\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,VISCOUS_Z,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseViscousFlowLaw_DepthDep(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"%s.method missing: A valid constructor must be specified",root->string);
  PetscFunctionReturn(0);
}

/* YIELD FUNCTIONS */
#undef __FUNCT__
#define __FUNCT__ "FoundationParseYieldFunction_Mises"
PetscErrorCode FoundationParseYieldFunction_Mises(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"YieldStress",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"YieldFunction.Mises.YieldStress missing");
  MaterialConstantsSetValues_PlasticMises(mc,ridx,val_double,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"YieldStressInf",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"YieldFunction.Mises.YieldStressInf missing");
  MaterialConstantsSetValues_PlasticMises(mc,ridx,PETSC_DEFAULT,val_double);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseYieldFunction_DP"
PetscErrorCode FoundationParseYieldFunction_DP(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"FrictionCoefficient",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"YieldFunction.DP.FrictionCoefficient missing");
  MaterialConstantsSetValues_PlasticDP(mc,ridx,val_double,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"FrictionCoefficientInf",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"YieldFunction.DP.FrictionCoefficientInf missing");
  MaterialConstantsSetValues_PlasticDP(mc,ridx,PETSC_DEFAULT,val_double,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"Cohesion",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"YieldFunction.DP.Cohesion missing");
  MaterialConstantsSetValues_PlasticDP(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,val_double,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  
  cJSON_GetObjectValue_double(root,"CohesionInf",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"YieldFunction.DP.CohesionInf missing");
  MaterialConstantsSetValues_PlasticDP(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,val_double,PETSC_DEFAULT,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"MaxTensileStress",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"YieldFunction.DP.MaxTensileStress missing");
  MaterialConstantsSetValues_PlasticDP(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,val_double,PETSC_DEFAULT);
  
  cJSON_GetObjectValue_double(root,"MaxStrainRate",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"YieldFunction.DP.MaxStrainRate missing");
  MaterialConstantsSetValues_PlasticDP(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,val_double);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseYieldFunction"
PetscErrorCode FoundationParseYieldFunction(pTatinCtx c,int ridx,cJSON *root)
{
  PetscErrorCode ierr;
  cJSON *jobj;
  DataBucket materialconstants;
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  
  ierr = FoundationParseJSONGetItemOptional(root,"None",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ None\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PLASTIC_NONE,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParse_None(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  ierr = FoundationParseJSONGetItemOptional(root,"Mises",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ Mises\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PLASTIC_MISES,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseYieldFunction_Mises(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  ierr = FoundationParseJSONGetItemOptional(root,"DruckerPrager",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ DruckerPrager\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PLASTIC_DP,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseYieldFunction_DP(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = FoundationParseJSONGetItemOptional(root,"MisesH",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ MisesH\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PLASTIC_MISES_H,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseYieldFunction_Mises(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = FoundationParseJSONGetItemOptional(root,"DruckerPragerH",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ DruckerPragerH\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PLASTIC_DP_H,PETSC_DEFAULT,PETSC_DEFAULT);
    ierr = FoundationParseYieldFunction_DP(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"%s.method missing: A valid constructor must be specified",root->string);
  PetscFunctionReturn(0);
}

/* SOFTENING FUNCTIONS */
#undef __FUNCT__
#define __FUNCT__ "FoundationParseStrainSoftening_Linear"
PetscErrorCode FoundationParseStrainSoftening_Linear(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"StrainMin",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"StrainSoftening.Linear.StrainMin missing");
  MaterialConstantsSetValues_SoftLin(mc,ridx,val_double,PETSC_DEFAULT);
  
  cJSON_GetObjectValue_double(root,"StrainMax",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"StrainSoftening.Linear.StrainMax missing");
  MaterialConstantsSetValues_SoftLin(mc,ridx,PETSC_DEFAULT,val_double);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseStrainSoftening_Exponential"
PetscErrorCode FoundationParseStrainSoftening_Exponential(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"StrainMin",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"StrainSoftening.Exponential.StrainMin missing");
  MaterialConstantsSetValues_SoftExpo(mc,ridx,val_double,PETSC_DEFAULT);
  
  cJSON_GetObjectValue_double(root,"EFold",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"StrainSoftening.Exponential.EFold missing");
  MaterialConstantsSetValues_SoftExpo(mc,ridx,PETSC_DEFAULT,val_double);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseStrainSoftening"
PetscErrorCode FoundationParseStrainSoftening(pTatinCtx c,int ridx,cJSON *root)
{
  PetscErrorCode ierr;
  cJSON *jobj;
  DataBucket materialconstants;
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  
  ierr = FoundationParseJSONGetItemOptional(root,"None",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ None\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PETSC_DEFAULT,SOFTENING_NONE,PETSC_DEFAULT);
    ierr = FoundationParse_None(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  ierr = FoundationParseJSONGetItemOptional(root,"Linear",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ Linear\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PETSC_DEFAULT,SOFTENING_LINEAR,PETSC_DEFAULT);
    ierr = FoundationParseStrainSoftening_Linear(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  ierr = FoundationParseJSONGetItemOptional(root,"Exponential",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ Exponential\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PETSC_DEFAULT,SOFTENING_EXPONENTIAL,PETSC_DEFAULT);
    ierr = FoundationParseStrainSoftening_Exponential(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"%s.method missing: A valid constructor must be specified",root->string);
  PetscFunctionReturn(0);
}

/* DENSITY FUNCTIONS */
#undef __FUNCT__
#define __FUNCT__ "FoundationParseDensity_Constant"
PetscErrorCode FoundationParseDensity_Constant(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"rho0",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Density.Const.rho0 missing");
  MaterialConstantsSetValues_DensityConst(mc,ridx,val_double);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseDensity_Boussinesq"
PetscErrorCode FoundationParseDensity_Boussinesq(DataBucket mc,int ridx,cJSON *root)
{
  int found;
  double val_double;
  
  cJSON_GetObjectValue_double(root,"rho0",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Density.Boussinesq.rho0 missing");
  MaterialConstantsSetValues_DensityBoussinesq(mc,ridx,val_double,PETSC_DEFAULT,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"alpha",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Density.Boussinesq.alpha missing");
  MaterialConstantsSetValues_DensityBoussinesq(mc,ridx,PETSC_DEFAULT,val_double,PETSC_DEFAULT);

  cJSON_GetObjectValue_double(root,"beta",&found,&val_double);
  if (!found) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Density.Boussinesq.beta missing");
  MaterialConstantsSetValues_DensityBoussinesq(mc,ridx,PETSC_DEFAULT,PETSC_DEFAULT,val_double);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseDensity"
PetscErrorCode FoundationParseDensity(pTatinCtx c,int ridx,cJSON *root)
{
  PetscErrorCode ierr;
  cJSON *jobj;
  DataBucket materialconstants;
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  
  ierr = FoundationParseJSONGetItemOptional(root,"Const",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ Const\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,DENSITY_CONSTANT);
    ierr = FoundationParseDensity_Constant(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  ierr = FoundationParseJSONGetItemOptional(root,"Boussinesq",&jobj);CHKERRQ(ierr);
  if (jobj) {
    PetscPrintf(PETSC_COMM_WORLD,"    |_______\\ Boussinesq\n");
    MaterialConstantsSetValues_MaterialType(materialconstants,ridx,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,DENSITY_BOUSSINESQ);
    ierr = FoundationParseDensity_Boussinesq(materialconstants,ridx,jobj);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"%s.method missing: A valid constructor must be specified",root->string);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FoundationParseRegion"
PetscErrorCode FoundationParseRegion(pTatinCtx c,Foundation f,cJSON *region)
{
  PetscErrorCode ierr;
  cJSON *jobj;
  int found,ridx;
  DataBucket materialconstants;
  
  
  ierr = FoundationParseJSONGetItemEssential(region,"index",&jobj);CHKERRQ(ierr);
  cJSON_GetObjectValue_int(region,"index",&found,&ridx);
  PetscPrintf(PETSC_COMM_WORLD,"  Region[%d]\n",ridx);

  ierr = FoundationParseJSONGetItemEssential(region,"ViscousFlowLaw",&jobj);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"  ViscousFlowLaw\n");
  ierr = FoundationParseViscousFlowLaw(c,ridx,jobj);CHKERRQ(ierr);
  
  ierr = FoundationParseJSONGetItemEssential(region,"YieldFunction",&jobj);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"  YieldFunction\n");
  ierr = FoundationParseYieldFunction(c,ridx,jobj);CHKERRQ(ierr);

  ierr = FoundationParseJSONGetItemEssential(region,"StrainSoftening",&jobj);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"  StrainSoftening\n");
  ierr = FoundationParseStrainSoftening(c,ridx,jobj);CHKERRQ(ierr);

  ierr = FoundationParseJSONGetItemEssential(region,"Density",&jobj);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"  Density\n");
  ierr = FoundationParseDensity(c,ridx,jobj);CHKERRQ(ierr);
  
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  MaterialConstantsPrintAll(materialconstants,ridx);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FoundationParseMaterialMetaDataContents"
PetscErrorCode FoundationParseMaterialMetaDataContents(pTatinCtx c,Foundation f,cJSON *root)
{
  PetscErrorCode ierr;
  cJSON *jregionlist,*jobj_k;
  int k,nregions;
  
  /* Locate RegionList */
  ierr = FoundationParseJSONGetItemEssential(root,"RegionList",&jregionlist);CHKERRQ(ierr);

  nregions = cJSON_GetArraySize(jregionlist);
  jobj_k = cJSON_GetArrayItemRoot(jregionlist);

  PetscPrintf(PETSC_COMM_WORLD,"Parsing RegionList [%d members]\n",nregions);
  for (k=0; k<nregions; k++) {
    ierr = FoundationParseRegion(c,f,jobj_k);CHKERRQ(ierr);
    jobj_k = cJSON_GetArrayItemNext(jobj_k);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationParseMaterialMetaData"
PetscErrorCode FoundationParseMaterialMetaData(pTatinCtx c,Foundation f)
{
  PetscErrorCode ierr;
  cJSON *jitemfile,*jmetaroot;
  int found;
  DataBucket materialconstants;
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  MaterialConstantsSetDefaults(materialconstants);

  ierr = FoundationParseJSONGetItemEssential(f->root,"MaterialModelMetaData",&jmetaroot);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Parsing MaterialModelMetaData\n");

  ierr = FoundationParseJSONGetItemOptional(jmetaroot,"MaterialMetaDataFile",&jitemfile);CHKERRQ(ierr);
  if (jitemfile) {
    cJSON *jfile;
    char *filename;
    
    cJSON_GetObjectValue_char(jmetaroot,"MaterialMetaDataFile",&found,&filename);
    PetscPrintf(PETSC_COMM_WORLD,"Parsing MaterialMetaDataFile %s\n",filename);
    cJSON_FileView(filename,&jfile);
    if (!jfile) {
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Failed to open JSON file: %s",filename);
    }
    ierr = FoundationParseMaterialMetaDataContents(c,f,jfile);CHKERRQ(ierr);
    
    cJSON_Delete(jfile);
  }
  ierr = FoundationParseMaterialMetaDataContents(c,f,jmetaroot);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
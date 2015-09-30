
#include "foundation.h"
#include "foundation_impl.h"
#include "data_bucket.h"
#include "material_point_utils.h"

typedef enum {
  MPStokes_density=0, MPStokes_viscosity,
  MPPlastic_plstrain,
  MP_diffusivity, MP_heatsource,
  MP_NULL
} MaterialPointType;

const char *MaterialPointTypeNames[] = {
  "density", "viscosity",
  "plasticstrain",
  "diffusivity",
  "heatsource",
  0
};

#undef __FUNCT__
#define __FUNCT__ "MaterialMPEval_Const"
PetscErrorCode MaterialMPEval_Const(pTatinCtx c,Foundation f,MaterialPointType mptype,int nr,int regionlist[],cJSON *jitem,DataBucket matpoint)
{
  cJSON *dummy;
  int r,found,n_mp_points,p,ridx;
  PetscErrorCode  ierr;
  MPAccess helper;
  PetscErrorCode (*fp_MaterialPointSetter_double)(MPAccess,const int,double);
  PetscErrorCode (*fp_MaterialPointSetter_float)(MPAccess,const int,float);
  PetscErrorCode (*fp_MaterialPointSetter_int)(MPAccess,const int,int);
  int nv;
  double regionvals[100];

  ierr = FoundationParseJSONGetItemEssential(jitem,"value",&dummy);CHKERRQ(ierr);
  cJSON_GetObjectValue_doublearray(jitem,"value",&found,&nv,regionvals);
  if (nr != nv) {
    SETERRQ3(PETSC_COMM_WORLD,PETSC_ERR_USER,"%s.value: Number regions (%d) must equal number of values (%d)",jitem->string,nr,nv);
  }
  
  DataBucketGetSizes(matpoint,&n_mp_points,0,0);
  MaterialPointGetAccess(matpoint,&helper);
  
  fp_MaterialPointSetter_double = NULL;
  fp_MaterialPointSetter_float  = NULL;
  fp_MaterialPointSetter_int    = NULL;
  switch (mptype) {
    case MPStokes_density:
      fp_MaterialPointSetter_double = &MaterialPointSet_density;
      break;
    case MPStokes_viscosity:
      fp_MaterialPointSetter_double = &MaterialPointSet_viscosity;
      break;
    case MP_diffusivity:
      fp_MaterialPointSetter_double = &MaterialPointSet_diffusivity;
      break;
    case MP_heatsource:
      fp_MaterialPointSetter_double = &MaterialPointSet_heat_source;
      break;
  }
  switch (mptype) {
    case MPPlastic_plstrain:
      fp_MaterialPointSetter_float = &MaterialPointSet_plastic_strain;
      break;
  }
  
  if (fp_MaterialPointSetter_double) {
    for (p=0; p<n_mp_points; p++) {
      MaterialPointGet_phase_index(helper,p,&ridx);
      for (r=0; r<nr; r++) {
        if (regionlist[r] == ridx) fp_MaterialPointSetter_double(helper,p,regionvals[r]);
      }
    }
  }
  if (fp_MaterialPointSetter_float) {
    for (p=0; p<n_mp_points; p++) {
      MaterialPointGet_phase_index(helper,p,&ridx);
      for (r=0; r<nr; r++) {
        if (regionlist[r] == ridx) fp_MaterialPointSetter_float(helper,p,(float)regionvals[r]);
      }
    }
  }
  if (fp_MaterialPointSetter_int) {
    for (p=0; p<n_mp_points; p++) {
      MaterialPointGet_phase_index(helper,p,&ridx);
      for (r=0; r<nr; r++) {
        if (regionlist[r] == ridx) fp_MaterialPointSetter_int(helper,p,(int)regionvals[r]);
      }
    }
  }

  MaterialPointRestoreAccess(matpoint,&helper);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialMPEval_Rand"
PetscErrorCode MaterialMPEval_Rand(pTatinCtx c,Foundation f,MaterialPointType mptype,int nr,int regionlist[],cJSON *jitem,DataBucket matpoint)
{
  cJSON *dummy;
  int r,found,n_mp_points,p,ridx;
  PetscErrorCode  ierr;
  MPAccess helper;
  PetscErrorCode (*fp_MaterialPointSetter_double)(MPAccess,const int,double);
  PetscErrorCode (*fp_MaterialPointSetter_float)(MPAccess,const int,float);
  PetscErrorCode (*fp_MaterialPointSetter_int)(MPAccess,const int,int);
  int nv;
  int regionseed[100];
  double regionmin[100];
  double regionmax[100];
  
  ierr = FoundationParseJSONGetItemEssential(jitem,"min",&dummy);CHKERRQ(ierr);
  cJSON_GetObjectValue_doublearray(jitem,"min",&found,&nv,regionmin);
  if (nr != nv) SETERRQ3(PETSC_COMM_WORLD,PETSC_ERR_USER,"%s.min: Number regions (%d) must equal number of min values (%d)",jitem->string,nr,nv);

  ierr = FoundationParseJSONGetItemEssential(jitem,"max",&dummy);CHKERRQ(ierr);
  cJSON_GetObjectValue_doublearray(jitem,"max",&found,&nv,regionmax);
  if (nr != nv) SETERRQ3(PETSC_COMM_WORLD,PETSC_ERR_USER,"%s.max: Number regions (%d) must equal number of max values (%d)",jitem->string,nr,nv);

  ierr = FoundationParseJSONGetItemOptional(jitem,"seed",&dummy);CHKERRQ(ierr);
  if (dummy) {
    cJSON_GetObjectValue_intarray(jitem,"seed",&found,&nv,regionseed);
    if (nr != nv) SETERRQ3(PETSC_COMM_WORLD,PETSC_ERR_USER,"%s.seed: Number regions (%d) must equal number of seed values (%d)",jitem->string,nr,nv);
  } else {
    for (r=0; r<nr; r++) {
      regionseed[r] = 0;
    }
  }
  
  DataBucketGetSizes(matpoint,&n_mp_points,0,0);
  MaterialPointGetAccess(matpoint,&helper);
  
  fp_MaterialPointSetter_double = NULL;
  fp_MaterialPointSetter_float  = NULL;
  fp_MaterialPointSetter_int    = NULL;
  switch (mptype) {
    case MPStokes_density:
      fp_MaterialPointSetter_double = &MaterialPointSet_density;
      break;
    case MPStokes_viscosity:
      fp_MaterialPointSetter_double = &MaterialPointSet_viscosity;
      break;
    case MP_diffusivity:
      fp_MaterialPointSetter_double = &MaterialPointSet_diffusivity;
      break;
    case MP_heatsource:
      fp_MaterialPointSetter_double = &MaterialPointSet_heat_source;
      break;
  }
  switch (mptype) {
    case MPPlastic_plstrain:
      fp_MaterialPointSetter_float = &MaterialPointSet_plastic_strain;
      break;
  }
  
  if (fp_MaterialPointSetter_double) {
    for (r=0; r<nr; r++) {
      srand(regionseed[r]);
      for (p=0; p<n_mp_points; p++) {
        double zero2one = rand()/((double)RAND_MAX);
        double rval = (regionmax[r]-regionmin[r]) * zero2one + regionmin[r];
        
        MaterialPointGet_phase_index(helper,p,&ridx);
        if (regionlist[r] == ridx) fp_MaterialPointSetter_double(helper,p,rval);
      }
    }
  }
  
  if (fp_MaterialPointSetter_float) {
    for (r=0; r<nr; r++) {
      srand(regionseed[r]);
      for (p=0; p<n_mp_points; p++) {
        float zero2one = rand()/((float)RAND_MAX);
        float rval = (regionmax[r]-regionmin[r]) * zero2one + regionmin[r];
        
        MaterialPointGet_phase_index(helper,p,&ridx);
        if (regionlist[r] == ridx) fp_MaterialPointSetter_float(helper,p,rval);
      }
    }
  }
  
  
  MaterialPointRestoreAccess(matpoint,&helper);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialMaterialState_Foundation"
PetscErrorCode ModelInitialMaterialState_Foundation(pTatinCtx c,void *ctx)
{
  Foundation      data = (Foundation)ctx;
  cJSON           *jroot,*jitem,*jtype,*jdummy,*jobj_k;
  int             i,entry,nics,k;
  DataBucket      matpoint;
  MaterialPointType type;
  PetscErrorCode  ierr;
  
  /* look for MeshGeomICType */
  ierr = FoundationParseJSONGetItemEssential(data->root,"MaterialPointStateIC",&jroot);CHKERRQ(ierr);
  nics = cJSON_GetArraySize(jroot);
  PetscPrintf(PETSC_COMM_WORLD,"[foundation] --------- Found %d methods\n",nics);

  ierr = pTatinGetMaterialPoints(c,&matpoint,NULL);CHKERRQ(ierr);

  /* traverse list */
  jobj_k = cJSON_GetArrayItemRoot(jroot);
  for (entry=0; entry<nics; entry++) {
  
    /* parse field type */
    jitem = NULL;
    for (i=0; i<(PetscInt)MP_NULL-1; i++) {
      ierr = FoundationParseJSONGetItemOptional(jobj_k,MaterialPointTypeNames[i],&jitem);CHKERRQ(ierr);
      if (jitem) {
        int found,nr,regionlist[100];
        
        type = (MaterialPointType)i;

        /* parse regions */
        ierr = FoundationParseJSONGetItemEssential(jitem,"RegionIndex",&jdummy);CHKERRQ(ierr);
        cJSON_GetObjectValue_intarray(jitem,"RegionIndex",&found,&nr,regionlist);
        
        PetscPrintf(PETSC_COMM_WORLD,"MaterialPointStateIC.%s: RegionIndex [ ",MaterialPointTypeNames[i]);
        for (k=0; k<nr; k++) PetscPrintf(PETSC_COMM_WORLD,"%d ",regionlist[k]);
        PetscPrintf(PETSC_COMM_WORLD,"]\n");

        /* parse function eval type */
        ierr = FoundationParseJSONGetItemOptional(jitem,"Const",&jtype);CHKERRQ(ierr);
        if (jtype) {
          //PetscPrintf(PETSC_COMM_WORLD,"|_______\\ Const\n");
          ierr = MaterialMPEval_Const(c,data,type,nr,regionlist,jtype,matpoint);CHKERRQ(ierr);
        }
        
        ierr = FoundationParseJSONGetItemOptional(jitem,"RandomNumber",&jtype);CHKERRQ(ierr);
        if (jtype) {
          //PetscPrintf(PETSC_COMM_WORLD,"|_______\\ RandomNumber\n");
          ierr = MaterialMPEval_Rand(c,data,type,nr,regionlist,jtype,matpoint);CHKERRQ(ierr);
        }
        
      }
    }

    jobj_k = cJSON_GetArrayItemNext(jobj_k);
  }
  /* parse regions */
  
  
  /* apply method */
  
  PetscFunctionReturn(0);
}






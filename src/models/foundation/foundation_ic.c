
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"

#include "dmda_bcs.h"
#include "ptatin3d_energy.h"
#include "rheology.h"







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


#undef __FUNCT__
#define __FUNCT__ "ModelInitialConditionMeshFields_Foundation"
PetscErrorCode ModelInitialConditionMeshFields_Foundation(pTatinCtx c,Vec X,void *ctx)
{
  Foundation      f = (Foundation)ctx;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialConditonMarkerFields_Foundation"
PetscErrorCode ModelInitialConditonMarkerFields_Foundation(pTatinCtx c,Vec X,void *ctx)
{
  Foundation      f = (Foundation)ctx;
  
  PetscFunctionReturn(0);
}


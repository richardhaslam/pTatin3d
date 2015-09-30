
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"

#include "ptatin3d.h"
#include "material_point_utils.h"
#include "output_material_points.h"

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Foundation"
PetscErrorCode ModelOutput_Foundation(pTatinCtx ptatinctx,Vec X,const char prefix[],void *modelctx)
{
  Foundation       f = (Foundation)modelctx;
  PhysCompStokes   stokes;
  DM               dmstokes,dmv,dmp;
  DataBucket       materialpoint;
  
  PetscErrorCode   ierr;
  
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
  
  /* get the velocity mesh */
  ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = DMCompositeGetEntries(dmstokes,&dmv,&dmp);CHKERRQ(ierr);

  /* ---- Velocity-Pressure Mesh Output ---- */
  /* v and p are written */
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);

  /* Light weight viewer: Only v is written. v and coords are expressed as floats */
  ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);

  /* ---- Temperature Mesh Output ---- */

  /* ---- Material Point Output ---- */
  /* [1] Basic viewer: Only reports coords, regionid and other internal data */
  /*
  ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
  */
  
  /* [2] Customized viewer: User defines specific fields they want to view - NOTE not .pvd file will be created */
  {
    const int                 nf = 2;
    const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes };
    char                      mp_file_prefix[PETSC_MAX_PATH_LEN];
    
    ierr = pTatinGetMaterialPoints(ptatinctx,&materialpoint,NULL);CHKERRQ(ierr);
    PetscSNPrintf(mp_file_prefix,PETSC_MAX_PATH_LEN-1,"%s_mpoints",prefix);
    ierr = SwarmViewGeneric_ParaView(materialpoint,nf,mp_prop_list,ptatinctx->outputpath,mp_file_prefix);CHKERRQ(ierr);
  }
  
  /* [3] Customized marker->cell viewer: Marker data is projected onto the velocity mesh. User defines specific fields */
  {
    const int                    nf = 2;
    const MaterialPointVariable  mp_prop_list[] = { MPV_viscosity, MPV_density };
    
    ierr = pTatin3d_ModelOutput_MarkerCellFields(ptatinctx,nf,mp_prop_list,prefix);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

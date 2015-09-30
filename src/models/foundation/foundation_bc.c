
#include "foundation.h"
#include "foundation_impl.h"
#include "foundation_user.h"

#include "dmda_bcs.h"
#include "ptatin3d_energy.h"

struct _p_StokesBCSet {
  FND_BoundaryConditionType bctype_list[6][3]; /* six faces, 3 equation */
  FND_BoundaryConditionMode bcmode_list[6][3];
};

struct _p_EnergyBCSet {
  FND_BoundaryConditionType bctype_list[6][1]; /* six faces, 1 equation */
};


#undef __FUNCT__
#define __FUNCT__ "FoundationBC_ApplyNeumann_Stokes"
PetscErrorCode FoundationBC_ApplyNeumann_Stokes(pTatinCtx c,Foundation f,PetscReal time,DM dmstokes)
{
	PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationBC_ApplyDirichlet_Velocity"
PetscErrorCode FoundationBC_ApplyDirichlet_Velocity(pTatinCtx c,Foundation f,PetscReal time,DM dmv,BCList bclist)
{
  PetscErrorCode ierr;
  PetscReal      val;

	PetscFunctionBegin;
  val = 0.0;
  
	/* u.n = 0 at base */
	ierr = DMDABCListTraverse3d(bclist,dmv,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	
	/* sigma.n = 0 at surface */
	
	/* u.n = 0, tau.t = 0 on lateral faces */
	ierr = DMDABCListTraverse3d(bclist,dmv,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dmv,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	
	/* u.n = 0, tau.t = 0 on front/rear faces */
	ierr = DMDABCListTraverse3d(bclist,dmv,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dmv,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FoundationBC_ApplyDirichlet_Temperature"
PetscErrorCode FoundationBC_ApplyDirichlet_Temperature(pTatinCtx c,Foundation f,PetscReal time,DM dmT,BCList bclist)
{
  PetscErrorCode ierr;
  PetscReal      val;

	PetscFunctionBegin;
  val = 0.0;
  ierr = DMDABCListTraverse3d(bclist,dmT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
  
  val = 1.0;
  ierr = DMDABCListTraverse3d(bclist,dmT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBC_Foundation"
PetscErrorCode ModelApplyBC_Foundation(pTatinCtx c,void *ctx)
{
  Foundation      f = (Foundation)ctx;
	PhysCompStokes  stokes;
	DM              dmstokes,dmv,dmp;
  PetscReal       time;
  PetscBool       active_energy;
	PetscErrorCode  ierr;
  
	PetscFunctionBegin;
  ierr = pTatinGetTime(c,&time);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
	ierr = DMCompositeGetEntries(dmstokes,&dmv,&dmp);CHKERRQ(ierr);
  
	ierr = FoundationBC_ApplyNeumann_Stokes(c,f,time,dmstokes);CHKERRQ(ierr);
	ierr = FoundationBC_ApplyDirichlet_Velocity(c,f,time,dmv,stokes->u_bclist);CHKERRQ(ierr);
 
  ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PhysCompEnergy energy;
		BCList         bclist;
		DM             dmT;
		
		ierr   = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		dmT    = energy->daT;
		bclist = energy->T_bclist;
    
    ierr = FoundationBC_ApplyDirichlet_Temperature(c,f,time,dmT,bclist);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBCMG_Foundation"
PetscErrorCode ModelApplyBCMG_Foundation(PetscInt nl,BCList bclist[],DM dmv[],pTatinCtx c,void *ctx)
{
  Foundation      f = (Foundation)ctx;
  PetscInt        k;
  PetscReal       time;
	PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = pTatinGetTime(c,&time);CHKERRQ(ierr);
  for (k=0; k<nl; k++) {
    ierr = FoundationBC_ApplyDirichlet_Velocity(c,f,time,dmv[k],bclist[k]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "swarm_fields.h"

#include "rheology.h"
#include "stokes_rheology_viscous.h"


#undef __FUNCT__
#define __FUNCT__ "RheologyConstantsInitialise"
PetscErrorCode RheologyConstantsInitialise(RheologyConstants *R)
{	
  PetscInt p;
  PetscBool flg;
  PetscScalar vis;
  PetscErrorCode ierr; 
  
	PetscFunctionBegin;
	
	R->nphases_active = 0;
	
  /* Define defaults for the viscosity cut-offs */
  R->apply_viscosity_cutoff_global = PETSC_FALSE;
  R->eta_lower_cutoff_global = 1.0e-100;
  R->eta_upper_cutoff_global = 1.0e+100;
  
  
  flg = PETSC_FALSE;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-eta_lower_cutoff_global",&vis,&flg);CHKERRQ(ierr);
  
  if (flg == PETSC_TRUE) { 
    R->apply_viscosity_cutoff_global = PETSC_TRUE;
    R->eta_lower_cutoff_global       = vis;
  }
  
	flg = PETSC_FALSE;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-eta_upper_cutoff_global",&vis,&flg);CHKERRQ(ierr);  
  
  if (flg == PETSC_TRUE) { 
    R->apply_viscosity_cutoff_global = PETSC_TRUE;
    R->eta_upper_cutoff_global       = vis;
  }
  
  PetscPrintf(PETSC_COMM_WORLD," global viscosity cut-off, min= %1.6e, max = %1.6e  \n", R->eta_lower_cutoff_global, R->eta_upper_cutoff_global );
  
  
  /* phase cutoff is equal to global cutoff and for the moment there is no options to enforce it 
   I don't think it belongs here ... maybe to the model definition */ 
  R->apply_viscosity_cutoff        = PETSC_FALSE;
  
  for (p=0; p<MAX_PHASE; p++) {
    R->eta_lower_cutoff[p] = R->eta_lower_cutoff_global;
    R->eta_upper_cutoff[p] = R->eta_upper_cutoff_global;
  }
  
  /* Define defaults for the parameters for each rheology */
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin_EvaluateRheologyNonlinearitiesMarkers"
PetscErrorCode pTatin_EvaluateRheologyNonlinearitiesMarkers(pTatinCtx user,DM dau,PetscScalar u[],DM dap,PetscScalar p[])
{	
  RheologyConstants *rheo;
	int               npoints;
	DataField         PField_std;
	DataField         PField_stokes;
  MPntStd           *mp_std;
	MPntPStokes       *mp_stokes;
	PetscErrorCode    ierr;
  static int        been_here=0;
	PhysCompStokes    stokes;
	
	PetscFunctionBegin;
  rheo = &user->rheology_constants;
	switch (rheo->rheology_type) {

		case RHEOLOGY_VISCOUS:
			if (been_here==0) {
				PetscPrintf(PETSC_COMM_WORLD,"*** WARNING: Rheology update for RHEOLOGY_VISCOUS using markers is under development ***\n");
				been_here = 1;
			}
			/* update on markers */
			ierr = EvaluateRheologyNonlinearitiesMarkers_Viscous(user,dau,u,dap,p);CHKERRQ(ierr);
			break;
			
		case RHEOLOGY_VISCO_PLASTIC:
			if (been_here==0) {
				PetscPrintf(PETSC_COMM_WORLD,"*** WARNING: Rheology update for RHEOLOGY_VISCO_PLASTIC using markers is under development ***\n");
				been_here = 1;
			}
			/* update on markers */
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Rheology update for RHEOLOGY_VISCO_PLASTIC using markers not defined");
			//ierr = EvaluateRheologyNonlinearitiesMarkers_ViscoPlastic(user,dau,u,dap,p);CHKERRQ(ierr);
			break;
			
		case RHEOLOGY_VISCO_PLASTIC_STRAIN_WEAKENING:
			if (been_here==0) {
				PetscPrintf(PETSC_COMM_WORLD,"*** WARNING: Rheology update for RHEOLOGY_VISCO_PLASTIC_STRAIN_WEAKENING using markers is under development ***\n");
				been_here = 1;
			}
			/* update on markers */
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Rheology update for RHEOLOGY_VISCO_PLASTIC_STRAIN_WEAKENING using markers not defined");
			//ierr = EvaluateRheologyNonlinearitiesMarkers_ViscoPlasticStrainWeakening(user,dau,u,dap,p);CHKERRQ(ierr);
			break;
			
		case RHEOLOGY_VISCO_ELASTIC_PLASTIC:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Rheology update for RHEOLOGY_VISCO_ELASTIC_PLASTIC using markers not defined");
			break;
			
		default:
			
			break;
	}

  /* Marker -> quadrature point projection */
	DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
	DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
	
	DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
	mp_std    = PField_std->data; /* should write a function to do this */
	mp_stokes = PField_stokes->data; /* should write a function to do this */

	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	
	switch (user->coefficient_projection_type) {

		case 0:			/* Perform P0 projection over Q2 element directly onyo uadrature points */
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"P0 marker->quadrature projection not supported");
			break;
			
		case 1:			/* Perform Q1 projection over Q2 element and interpolate back to quadrature points */
			ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(npoints,mp_std,mp_stokes,stokes->dav,stokes->volQ);CHKERRQ(ierr);
			break;
			
		case 2: 			/* Perform Q2 projection and interpolate back to quadrature points */
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Q2 marker->quadrature projection not supported");
			break;
	}

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin_EvaluateRheologyNonlinearities"
PetscErrorCode pTatin_EvaluateRheologyNonlinearities(pTatinCtx user,DM dau,PetscScalar u[],DM dap,PetscScalar p[])
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = pTatin_EvaluateRheologyNonlinearitiesMarkers(user,dau,u,dap,p);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

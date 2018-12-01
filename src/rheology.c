/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   rheology.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "data_bucket.h"

#include "rheology.h"
#include "stokes_rheology_viscous.h"
#include "stokes_rheology_vp_std.h"
#include "stokes_rheology_lava.h"

PetscLogEvent PTATIN_CoefficientEvaluate;
PetscLogEvent PTATIN_CoefficientEvolve;

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
  ierr = PetscOptionsGetReal(NULL,NULL,"-eta_lower_cutoff_global",&vis,&flg);CHKERRQ(ierr);

  if (flg == PETSC_TRUE) {
    R->apply_viscosity_cutoff_global = PETSC_TRUE;
    R->eta_lower_cutoff_global       = vis;
  }

  flg = PETSC_FALSE;
  ierr = PetscOptionsGetReal(NULL,NULL,"-eta_upper_cutoff_global",&vis,&flg);CHKERRQ(ierr);

  if (flg == PETSC_TRUE) {
    R->apply_viscosity_cutoff_global = PETSC_TRUE;
    R->eta_upper_cutoff_global       = vis;
  }
  /*
  PetscPrintf(PETSC_COMM_WORLD,"RheologyConstantsInitialise: global viscosity cut-off, min= %1.6e, max = %1.6e  \n", R->eta_lower_cutoff_global, R->eta_upper_cutoff_global );
  */

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
      if (been_here == 0) {
        PetscPrintf(PETSC_COMM_WORLD,"*** Rheology update for RHEOLOGY_VISCOUS selected ***\n");
      }
      /* update on markers */
      ierr = EvaluateRheologyNonlinearitiesMarkers_Viscous(user,dau,u,dap,p);CHKERRQ(ierr);
      break;

    case RHEOLOGY_VP_STD:
      if (been_here == 0) {
        PetscPrintf(PETSC_COMM_WORLD,"*** Rheology update for RHEOLOGY_VP_STD selected ***\n");
      }
      ierr = EvaluateRheologyNonlinearitiesMarkers_VPSTD(user,dau,u,dap,p);CHKERRQ(ierr);
      ierr = ApplyViscosityCutOffMarkers_VPSTD(user);CHKERRQ(ierr);
      break;

    case RHEOLOGY_LAVA:
      if (been_here == 0) {
        PetscPrintf(PETSC_COMM_WORLD,"*** Rheology update for RHEOLOGY_LAVA selected ***\n");
      }
      ierr = EvaluateRheologyNonlinearitiesMarkers_LAVA(user,dau,u,dap,p);CHKERRQ(ierr);
      ierr = ApplyViscosityCutOffMarkers_VPSTD(user);CHKERRQ(ierr);
      break;

    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Rheology update is not defined");
      break;
  }

  /* Marker -> quadrature point projection */
  DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
  DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);

  DataBucketGetSizes(user->materialpoint_db,&npoints,NULL,NULL);
  mp_std    = PField_std->data; /* should write a function to do this */
  mp_stokes = PField_stokes->data; /* should write a function to do this */

  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);

  switch (user->coefficient_projection_type) {

    case -1:      /* Perform null projection use the values currently defined on the quadrature points */
      break;

    case 0:     /* Perform P0 projection over Q2 element directly onto quadrature points */
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"P0 [arithmetic avg] marker->quadrature projection not supported");
            ierr = MPntPStokesProj_P0(CoefAvgARITHMETIC,npoints,mp_std,mp_stokes,stokes->dav,stokes->volQ);CHKERRQ(ierr);
      break;
    case 10:
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"P0 [harmonic avg] marker->quadrature projection not supported");
            ierr = MPntPStokesProj_P0(CoefAvgHARMONIC,npoints,mp_std,mp_stokes,stokes->dav,stokes->volQ);CHKERRQ(ierr);
      break;
    case 20:
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"P0 [geometric avg] marker->quadrature projection not supported");
            ierr = MPntPStokesProj_P0(CoefAvgGEOMETRIC,npoints,mp_std,mp_stokes,stokes->dav,stokes->volQ);CHKERRQ(ierr);
      break;
    case 30:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"P0 [dominant phase] marker->quadrature projection not supported");
      break;

    case 1:     /* Perform Q1 projection over Q2 element and interpolate back to quadrature points */
      ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(npoints,mp_std,mp_stokes,stokes->dav,stokes->volQ);CHKERRQ(ierr);
      break;

    case 2:       /* Perform Q2 projection and interpolate back to quadrature points */
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Q2 marker->quadrature projection not supported");
      break;

    case 3:       /* Perform P1 projection and interpolate back to quadrature points */
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"P1 marker->quadrature projection not supported");
      break;
    case 4:
      ierr = SwarmUpdateGaussPropertiesOne2OneMap_MPntPStokes(npoints,mp_std,mp_stokes,stokes->volQ);CHKERRQ(ierr);
      break;

    case 5:
      ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_sorted(user->materialpoint_db,stokes->dav,stokes->volQ);CHKERRQ(ierr);
      break;
      
    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Viscosity projection type is not defined");
      break;
  }

  ierr = pTatin_ApplyStokesGravityModel(user);CHKERRQ(ierr);

  been_here = 1;
  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_EvaluateRheologyNonlinearities(pTatinCtx user,DM dau,PetscScalar u[],DM dap,PetscScalar p[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscLogEventBegin(PTATIN_CoefficientEvaluate,0,0,0,0);CHKERRQ(ierr);
  ierr = pTatin_EvaluateRheologyNonlinearitiesMarkers(user,dau,u,dap,p);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(PTATIN_CoefficientEvaluate,0,0,0,0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_EvaluateCoefficientNonlinearities_Stokes(pTatinCtx ptatin,Vec X)
{
  PetscErrorCode    ierr;
  DM                stokes_pack,dau,dap;
  Vec               Uloc,Ploc;
  PetscScalar       *LA_Uloc,*LA_Ploc;
  PhysCompStokes    stokes;

  PetscFunctionBegin;

  ierr = pTatinGetStokesContext(ptatin,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;

  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);

  /* get the local (ghosted) entries for each physics */
  ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
  ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);

  /* ======================================== */
  /*         UPDATE NON-LINEARITIES           */
  /* evaluate rheology and rhs using X        */
  /* map marker eta to quadrature points */
  /* map marker force to quadrature points */
  /* ======================================== */
  ierr = pTatin_EvaluateRheologyNonlinearities(ptatin,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);

  ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
  ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_StokesCoefficient_UpdateTimeDependentQuantities(pTatinCtx user,DM dau,PetscScalar u[],DM dap,PetscScalar p[])
{
  RheologyConstants *rheo;
  static int        been_here = 0;
  PetscErrorCode    ierr;

  PetscFunctionBegin;

  ierr = PetscLogEventBegin(PTATIN_CoefficientEvolve,0,0,0,0);CHKERRQ(ierr);
  ierr = pTatinGetRheology(user,&rheo);

  switch (rheo->rheology_type) {

    case RHEOLOGY_VISCOUS:
      if (been_here == 0) {
        PetscPrintf(PETSC_COMM_WORLD,"*** StokesCoefficientUpdate for RHEOLOGY_VISCOUS is NULL ***\n");
      }
      break;

    case RHEOLOGY_VP_STD:
    {
      DataBucket     db;
      static BTruth  found;

      if (been_here == 0) {
        /* access material point information */
        ierr = pTatinGetMaterialPoints(user,&db,NULL);CHKERRQ(ierr);
        /* check plastic marker type is loaded */
        DataBucketQueryDataFieldByName(db,MPntPStokesPl_classname,&found);
        if (found == BFALSE) {
          PetscPrintf(PETSC_COMM_WORLD,"*** WARNING: Update of plastic strain requires you register the marker field: %s. Ignoring update of plastic strain.\n", MPntPStokesPl_classname);
        }
      }

      /* call */
      if (found == BTRUE) {
        ierr = StokesCoefficient_UpdateTimeDependentQuantities_VPSTD(user,dau,u,dap,p);CHKERRQ(ierr);
      }
    }
      break;

    case RHEOLOGY_LAVA:
      if (been_here == 0) {
        PetscPrintf(PETSC_COMM_WORLD,"*** StokesCoefficientUpdate for RHEOLOGY_LAVA is NULL ***\n");
      }
      break;

    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"StokesCoefficientUpdate is not defined");
      break;
  }
  ierr = PetscLogEventEnd(PTATIN_CoefficientEvolve,0,0,0,0);CHKERRQ(ierr);


  been_here = 1;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_UpdateCoefficientTemporalDependence_Stokes(pTatinCtx ptatin,Vec X)
{
  PetscErrorCode    ierr;
  DM                stokes_pack,dau,dap;
  Vec               Uloc,Ploc;
  PetscScalar       *LA_Uloc,*LA_Ploc;
  PhysCompStokes    stokes;

  PetscFunctionBegin;

  ierr = pTatinGetStokesContext(ptatin,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;

  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);

  /* get the local (ghosted) entries for each physics */
  ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
  ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);

  ierr = pTatin_StokesCoefficient_UpdateTimeDependentQuantities(ptatin,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);

  ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
  ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin_ApplyStokesGravityModel(pTatinCtx ctx)
{
  PetscErrorCode    ierr;
  PhysCompStokes    stokes;
  PetscInt          e,nel,q,nqp;
  QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
  PetscFunctionBegin;

  ierr = pTatinGetStokesContext(ctx,&stokes);CHKERRQ(ierr);

  nel = stokes->volQ->n_elements;
  nqp = stokes->volQ->npoints;
  ierr = VolumeQuadratureGetAllCellData_Stokes(stokes->volQ,&all_gausspoints);CHKERRQ(ierr);
  for (e=0; e<nel; e++) {
    ierr = VolumeQuadratureGetCellData_Stokes(stokes->volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
    for (q=0; q<nqp; q++) {
      cell_gausspoints[q].Fu[0] = stokes->gravity_vector[0] * cell_gausspoints[q].rho;
      cell_gausspoints[q].Fu[1] = stokes->gravity_vector[1] * cell_gausspoints[q].rho;
      cell_gausspoints[q].Fu[2] = stokes->gravity_vector[2] * cell_gausspoints[q].rho;
    }
  }

  PetscFunctionReturn(0);
}

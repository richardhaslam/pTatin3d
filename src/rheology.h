/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      rheology.h
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

#ifndef __ptatin_rheology_h__
#define __ptatin_rheology_h__

#include "petsc.h"
#include "petscdm.h"
#include "ptatin3d.h"

typedef enum { 
	RHEOLOGY_VISCOUS=0, 
	RHEOLOGY_VP_STD
} RheologyType;





#define MAX_PHASE 100

struct _p_RheologyConstants {
  PetscBool apply_viscosity_cutoff_global;
  PetscReal eta_lower_cutoff_global, eta_upper_cutoff_global;
  PetscBool apply_viscosity_cutoff;
  PetscReal eta_lower_cutoff[MAX_PHASE], eta_upper_cutoff[MAX_PHASE];
  
  /* rheology; "constant" ==>> const (short name) */
  PetscReal const_eta0[MAX_PHASE];
  PetscReal const_rho0[MAX_PHASE];
	/* not in used in the marker version, elastic const,would require storing stress on marker... oh la la */
  PetscReal const_shearmod[MAX_PHASE];    
  /* rheology; "von mises" ==>> mises (short name) */
  PetscReal mises_tau_yield[MAX_PHASE]; /* would be nice to refactor it to DP_Co */ 
  PetscReal dp_pressure_dependance[MAX_PHASE]; /* would be nice to refactor it to DP_phi */
  PetscReal tens_cutoff[MAX_PHASE]; /* would be nice to refactor it to DP_tens_cutoff */
  PetscReal Hst_cutoff [MAX_PHASE]; /* would be nice to refactor it to DP_Hst_cutoff */
	
  PetscReal gamma_soft[MAX_PHASE]; /* not in used in the marker version */ 
  PetscReal mises_tau_yield_inf[MAX_PHASE]; /* not in used in the marker version */
	
  /* rheology; "strain_weakening" ==> soft (short name) */
  PetscReal soft_min_strain_cutoff [MAX_PHASE];
  PetscReal soft_max_strain_cutoff [MAX_PHASE];
  PetscReal soft_Co_inf [MAX_PHASE];
  PetscReal soft_phi_inf [MAX_PHASE];
	
	PetscInt     nphases_active;
	RheologyType rheology_type;
};


PetscErrorCode RheologyConstantsInitialise(RheologyConstants *R);
PetscErrorCode pTatin_EvaluateRheologyNonlinearities(pTatinCtx user,DM dau,PetscScalar u[],DM dap,PetscScalar p[]);
PetscErrorCode pTatin_StokesCoefficient_UpdateTimeDependentQuantities(pTatinCtx user,DM dau,PetscScalar u[],DM dap,PetscScalar p[]);

PetscErrorCode pTatin_EvaluateCoefficientNonlinearities_Stokes(pTatinCtx ptatin,Vec X);
PetscErrorCode pTatin_UpdateCoefficientTemporalDependence_Stokes(pTatinCtx ptatin,Vec X);

#endif


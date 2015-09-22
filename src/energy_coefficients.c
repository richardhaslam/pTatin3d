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
 **    filename:   energy_coefficients.c
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
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "phys_comp_energy.h"
#include "ptatin3d_energy.h"
#include "material_constants_energy.h"
#include "dmda_element_q1.h"
#include "element_utils_q1.h"
#include "MPntPEnergy_def.h"
#include "material_point_utils.h"


#undef __FUNCT__
#define __FUNCT__ "EnergyEvaluateCoefficients_MaterialPoints"
PetscErrorCode EnergyEvaluateCoefficients_MaterialPoints(pTatinCtx user,PetscReal time,DM dmT,PetscScalar LA_T[],PetscScalar LA_U[])
{
	PetscErrorCode ierr;
	DataBucket     material_constants,material_points;
	DataField      PField_MatConsts,PField_SourceConst,PField_SourceDecay,PField_SourceAdiAdv,PField_ConductivityConst,PField_ConductivityThreshold;
	DataField      PField_std,PField_energy;
	EnergyMaterialConstants        *mat_consts;
	EnergySourceConst              *source_const;
	EnergySourceDecay              *source_decay;
	EnergySourceAdiabaticAdvection *source_adi_adv;
  EnergyConductivityConst        *k_const;
  EnergyConductivityThreshold    *k_threshold;
	int       pidx,n_mp_points;
  PetscInt  k,nel,nen;
  const     PetscInt *elnidx;
	PetscReal el_T[Q1_NODES_PER_EL_3D],el_U[Q1_NODES_PER_EL_3D],NQ1[Q1_NODES_PER_EL_3D];
	PhysCompStokes stokes;
  PetscReal *grav_vec;

	PetscFunctionBegin;
	
	/* Get bucket of material constants */
	ierr = pTatinGetMaterialConstants(user,&material_constants);CHKERRQ(ierr);
  
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	grav_vec = stokes->gravity_vector;

	/* fetch array to data for material constants */
	DataBucketGetDataFieldByName(material_constants,EnergyMaterialConstants_classname,&PField_MatConsts);
	DataFieldGetEntries(PField_MatConsts,(void**)&mat_consts);

	/* fetch array to data for source method */
	DataBucketGetDataFieldByName(material_constants, EnergySourceConst_classname, &PField_SourceConst );
	DataFieldGetEntries(PField_SourceConst,(void**)&source_const);
	DataBucketGetDataFieldByName(material_constants, EnergySourceDecay_classname, &PField_SourceDecay );
	DataFieldGetEntries(PField_SourceDecay,(void**)&source_decay);
	DataBucketGetDataFieldByName(material_constants, EnergySourceAdiabaticAdvection_classname, &PField_SourceAdiAdv );
	DataFieldGetEntries(PField_SourceAdiAdv,(void**)&source_adi_adv);

	/* fetch array to data for conductivity method */
	DataBucketGetDataFieldByName(material_constants, EnergyConductivityConst_classname, &PField_ConductivityConst );
	DataFieldGetEntries(PField_ConductivityConst,(void**)&k_const);
	DataBucketGetDataFieldByName(material_constants, EnergyConductivityThreshold_classname, &PField_ConductivityThreshold );
	DataFieldGetEntries(PField_ConductivityThreshold,(void**)&k_threshold);
	
	/* Get bucket of material points */
	ierr = pTatinGetMaterialPoints(user,&material_points,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(material_points,&n_mp_points,0,0);
	
	DataBucketGetDataFieldByName(material_points,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataBucketGetDataFieldByName(material_points,MPntPEnergy_classname,&PField_energy);
	DataFieldGetAccess(PField_energy);

  ierr = DMDAGetElementsQ1(dmT,&nel,&nen,&elnidx);CHKERRQ(ierr);

	for (pidx=0; pidx<n_mp_points; pidx++) {
		MPntStd       *mp_std;
		MPntPEnergy   *mpp_energy;
		double        *xi_mp,T_mp,u_mp[3];
		int           t,eidx,region_idx;
		double        rho_mp,conductivity_mp,diffusivity_mp,H_mp,Cp;
		int           density_type,conductivity_type;
		int           *source_type;
    
		DataFieldAccessPoint(PField_std,    pidx,(void**)&mp_std);
		DataFieldAccessPoint(PField_energy, pidx,(void**)&mpp_energy);
		
		/* Get index of element containing this marker */
		eidx = mp_std->wil;
		/* Get marker local coordinate (for interpolation) */
		xi_mp = mp_std->xi;
		
		/* Get region index */
    region_idx = mp_std->phase;
				
    /* Get element temperature */
    ierr = DMDAEQ1_GetScalarElementField_3D(el_T,(PetscInt*)&elnidx[nen * eidx],LA_T);CHKERRQ(ierr);
    
		T_mp = 0.0;
    P3D_ConstructNi_Q1_3D(xi_mp,NQ1);
    
    for (k=0; k<Q1_NODES_PER_EL_3D; k++) {
      T_mp += NQ1[k] * el_T[k];
    }
    
		/* get velocity for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_U,(PetscInt*)&elnidx[nen * eidx],LA_U);CHKERRQ(ierr);

		u_mp[0] = u_mp[1] = u_mp[2] = 0.0;
    for (k=0; k<Q1_NODES_PER_EL_3D; k++) {
      u_mp[0] += NQ1[k] * el_U[3*k+0];      /* compute vx on the particle */
      u_mp[1] += NQ1[k] * el_U[3*k+1];      /* compute vy on the particle */
      u_mp[2] += NQ1[k] * el_U[3*k+2];      /* compute vz on the particle */
    }
		
		//density_type      = ematconsts[ region_idx ].density_type;
		density_type      = 0; /* Hardcoded to be ENERGYDENSITY_CONSTANT */
		conductivity_type = mat_consts[ region_idx ].conductivity_type;
		source_type       = mat_consts[ region_idx ].source_type;
		
		/* Fetch value for Cp */
		Cp = mat_consts[ region_idx ].Cp;

		/* Compute density */
		rho_mp = 1.0;
		switch (density_type) {
      case ENERGYDENSITY_NONE:
        rho_mp = 1.0;
        Cp = 1.0;
        break;
			
      case ENERGYDENSITY_CONSTANT:
				rho_mp = mat_consts[ region_idx ].rho_ref;
			break;

			case ENERGYDENSITY_BOUSSINESQ:
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"BOUSSINESQ is not available - sorry email GD for help");
				break;
		}

		/* Compute conductivity */
		conductivity_mp = 1.0;
		switch (conductivity_type) {
      case ENERGYCONDUCTIVITY_DEFAULT:
				conductivity_mp = mpp_energy->diffusivity;
        break;

			case ENERGYCONDUCTIVITY_CONSTANT:
        conductivity_mp = k_const[ region_idx ].k0;
				break;
			
      case ENERGYCONDUCTIVITY_TEMP_DEP_THRESHOLD:
        /*
        conductivity_mp = k_threshold[ region_idx ].k0;
        if (T_mp >= k_threshold[ region_idx ].T0) {
          conductivity_mp = k_threshold[ region_idx ].k1;
        }
        */
        conductivity_mp = k_threshold[ region_idx ].k0;
        if (k_threshold[ region_idx ].T_threshold - T_mp < k_threshold[ region_idx ].dT) {
          double shift_T = T_mp - (k_threshold[ region_idx ].T_threshold - k_threshold[ region_idx ].dT);
          double dk = k_threshold[ region_idx ].k1 - k_threshold[ region_idx ].k0;
          
          conductivity_mp = k_threshold[ region_idx ].k0 + (dk/k_threshold[ region_idx ].dT)*shift_T;
        } else if (T_mp >= k_threshold[ region_idx ].T_threshold) {
          conductivity_mp = k_threshold[ region_idx ].k1;
        }
				break;
		}

		/*
		  Compute heat sources
		    Note: We want to allow multiple heat sources to exists.
		    Presently 6 choices are available, we loop through all
		    possible cases and sum the resulting source
		*/
		H_mp = 0.0;
		for (t=0; t<6; t++) {
			switch (source_type[t]) {
				case ENERGYCONDUCTIVITY_DEFAULT:
					H_mp += mpp_energy->heat_source;
					break;
					
				case ENERGYSOURCE_CONSTANT:
					H_mp += source_const[ region_idx ].H;
					break;
					
				case ENERGYSOURCE_SHEAR_HEATING:
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"SHEAR-HEATING is not available");
					break;
					
				case ENERGYSOURCE_DECAY:
					H_mp += source_decay[ region_idx ].H0 * exp( -time * source_decay[ region_idx ].lambda );
					break;

				/*
						Taken from T. Gerya, "Introduction ot numerical geodynamic modelling"
				        page 156-157
        */
				case ENERGYSOURCE_ADIABATIC:
				{
					double g_dot_v; /* g_i * u_i */
          
					//g_dot_v = -(1.0)*u_mp[1]; /* todo - needs to be generalized to use gravity vector */
          
          g_dot_v = -( grav_vec[0]*u_mp[0] + grav_vec[1]*u_mp[1] + grav_vec[2]*u_mp[2] );
					
          H_mp += T_mp * mat_consts[ region_idx ].alpha * rho_mp * g_dot_v;
				}
					break;
					
        /*
         vector u point in the direction of gravity
        */
				case ENERGYSOURCE_ADIABATIC_ADVECTION:
				{
					double grav_nrm,u_vertical;
					     
					//u_vertical = u_mp[1]; /* todo - needs to be generalized to use gravity vector */

          grav_nrm = PetscSqrtReal( grav_vec[0]*grav_vec[0] + grav_vec[1]*grav_vec[1] + grav_vec[2]*grav_vec[2] );
          u_vertical = -(u_mp[0]*grav_vec[0] + u_mp[1]*grav_vec[1] + u_mp[2]*grav_vec[2])/grav_nrm;
          
          H_mp += rho_mp * Cp * u_vertical * ( source_adi_adv[ region_idx ].dTdy );
					
				}
					break;
			}
		}
		
		diffusivity_mp = conductivity_mp / (rho_mp * Cp);
		
		H_mp = H_mp / (rho_mp * Cp);
		
		MPntPEnergySetField_diffusivity(mpp_energy,diffusivity_mp);
		MPntPEnergySetField_heat_source(mpp_energy,H_mp);
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_energy);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EnergyEvaluateCoefficients"
PetscErrorCode EnergyEvaluateCoefficients(pTatinCtx user,PetscReal time,DM dmT,PetscScalar LA_T[],PetscScalar LA_U[])
{
	PetscErrorCode ierr;
  DataBucket     materialpoint;
  Quadrature     volQ;
	PhysCompEnergy energy;
  
	PetscFunctionBegin;
	
  /* Evaluate physics on material points */
	ierr = EnergyEvaluateCoefficients_MaterialPoints(user,time,dmT,LA_T,LA_U);CHKERRQ(ierr);
	
  /* Project effective diffusivity and source from material points to quadrature points */
  ierr = pTatinGetContext_Energy(user,&energy);CHKERRQ(ierr);
	volQ = energy->volQ;
	ierr = pTatinGetMaterialPoints(user,&materialpoint,NULL);CHKERRQ(ierr);

	ierr = MaterialPointQuadraturePointProjectionC0_Q2Energy(dmT,materialpoint,MPField_Energy,MPPEgy_diffusivity,volQ);CHKERRQ(ierr);
	ierr = MaterialPointQuadraturePointProjectionC0_Q2Energy(dmT,materialpoint,MPField_Energy,MPPEgy_heat_source,volQ);CHKERRQ(ierr);
  
	PetscFunctionReturn(0);
}


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

#undef __FUNCT__
#define __FUNCT__ "EnergyEvaluateCoefficients_MaterialPoints"
PetscErrorCode EnergyEvaluateCoefficients_MaterialPoints(pTatinCtx user,DM dmT,PetscScalar LA_T[],DM dmU,PetscScalar LA_U[])
{
	PetscErrorCode ierr;
	DataBucket     material_constants,material_points;
	DataField      PField_MatConsts,PField_SourceConst,PField_SourceDecay,PField_SourceAdiAdv;
	DataField      PField_std,PField_energy;
	EnergyMaterialConstants        *mat_consts;
	EnergySourceConst              *source_const;
	EnergySourceDecay              *source_decay;
	EnergySourceAdiabaticAdvection *source_adi_adv;
	int pidx,n_mp_points;
	PetscReal time;
	PetscFunctionBegin;

	ierr = pTatinGetTime(user,&time);CHKERRQ(ierr);
	
	/* Get bucket of material constants */
	ierr = pTatinGetMaterialConstants(user,&material_constants);CHKERRQ(ierr);
	
	/* fetch array to material constants */
	DataBucketGetDataFieldByName(material_constants,EnergyMaterialConstants_classname,&PField_MatConsts);
	DataFieldGetEntries(PField_MatConsts,(void**)&mat_consts);

	/* fetch array to source type */
	DataBucketGetDataFieldByName(material_constants, EnergySourceConst_classname, &PField_SourceConst );
	DataFieldGetEntries(PField_SourceConst,(void**)&source_const);
	DataBucketGetDataFieldByName(material_constants, EnergySourceDecay_classname, &PField_SourceDecay );
	DataFieldGetEntries(PField_SourceDecay,(void**)&source_decay);
	DataBucketGetDataFieldByName(material_constants, EnergySourceAdiabaticAdvection_classname, &PField_SourceAdiAdv );
	DataFieldGetEntries(PField_SourceAdiAdv,(void**)&source_adi_adv);

	
	/* Get bucket of material points */
	ierr = pTatinGetMaterialPoints(user,&material_points,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(material_points,&n_mp_points,0,0);
	
	DataBucketGetDataFieldByName(material_points,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataBucketGetDataFieldByName(material_points,MPntPEnergy_classname,&PField_energy);
	DataFieldGetAccess(PField_energy);

	
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
		
		
		T_mp = 0.0;
		u_mp[0] = u_mp[1] = u_mp[2] = 0.0;

		
		
		//density_type      = ematconsts[ region_idx ].density_type;
		density_type      = 0; /* Hardcoded to be ENERGYDENSITY_CONSTANT */
		conductivity_type = mat_consts[ region_idx ].conductivity_type;
		source_type       = mat_consts[ region_idx ].source_type;
		
		/* Fetch value for Cp */
		Cp = mat_consts[ region_idx ].Cp;

		/* Compute density */
		rho_mp = 1.0;
		switch (density_type) {
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
			case ENERGYCONDUCTIVITY_CONSTANT:
				conductivity_mp = mpp_energy->diffusivity;
				break;
			case ENERGYCONDUCTIVITY_TEMP_DEP_THRESHOLD:
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
				case ENERGYSOURCE_NONE:
					break;
					
				case ENERGYSOURCE_CONSTANT:
					H_mp += mpp_energy->heat_source;
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
					
					g_dot_v = -(1.0); /* todo */
					H_mp += T_mp * mat_consts[ region_idx ].alpha * rho_mp * g_dot_v;
				}
					break;
					
				case ENERGYSOURCE_ADIABATIC_ADVECTION:
				{
					double u_vertical;
					
					u_vertical = u_mp[1]; /* GD says this is fine */
					H_mp += rho_mp * Cp * u_vertical * (-source_adi_adv[ region_idx ].dTdy);
					
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
PetscErrorCode EnergyEvaluateCoefficients(pTatinCtx user,DM dmT,PetscScalar LA_T[],DM dmU,PetscScalar LA_U[])
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	ierr = EnergyEvaluateCoefficients_MaterialPoints(user,dmT,LA_T,dmU,LA_U);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


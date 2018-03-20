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
 **    filename:   material_constants_energy.c
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
#include "data_bucket.h"
#include "material_constants_energy.h"

PetscErrorCode MaterialConstantsEnergyInitialize(DataBucket db)
{
  DataBucketRegisterField(db,EnergyMaterialConstants_classname,       sizeof(EnergyMaterialConstants),NULL);
  DataBucketRegisterField(db,EnergyConductivityConst_classname,       sizeof(EnergyConductivityConst),NULL);
  DataBucketRegisterField(db,EnergyConductivityThreshold_classname,       sizeof(EnergyConductivityThreshold),NULL);
  DataBucketRegisterField(db,EnergySourceConst_classname,             sizeof(EnergySourceConst),NULL);
  DataBucketRegisterField(db,EnergySourceDecay_classname,             sizeof(EnergySourceDecay),NULL);
  DataBucketRegisterField(db,EnergySourceAdiabaticAdvection_classname,sizeof(EnergySourceAdiabaticAdvection),NULL);
  DataBucketFinalize(db);

  PetscFunctionReturn(0);
}

PetscErrorCode MaterialConstantsEnergySetDefaults(DataBucket db)
{
	int       nregions;
  DataField dfield;
  
	DataBucketGetSizes(db,&nregions,NULL,NULL);
  {
    EnergyMaterialConstants *data;

    DataBucketGetDataFieldByName(db,EnergyMaterialConstants_classname,&dfield);
    DataFieldGetEntries(dfield,(void**)&data);
    MaterialConstantsSetDefaultAll_EnergyMaterialConstants(nregions,data);
  }

  /* conductivity */
  {
    EnergyConductivityConst *data;
    
    DataBucketGetDataFieldByName(db,EnergyConductivityConst_classname,&dfield);
    DataFieldGetEntries(dfield,(void**)&data);
    MaterialConstantsSetDefaultAll_ConductivityConst(nregions,data);
  }

  {
    EnergyConductivityThreshold *data;
    
    DataBucketGetDataFieldByName(db,EnergyConductivityThreshold_classname,&dfield);
    DataFieldGetEntries(dfield,(void**)&data);
    MaterialConstantsSetDefaultAll_ConductivityThreshold(nregions,data);
  }
  
  /* sources */
  {
    EnergySourceConst *data;
    
    DataBucketGetDataFieldByName(db,EnergySourceConst_classname,&dfield);
    DataFieldGetEntries(dfield,(void**)&data);
    MaterialConstantsSetDefaultAll_SourceConst(nregions,data);
  }

  {
    EnergySourceDecay *data;
    
    DataBucketGetDataFieldByName(db,EnergySourceDecay_classname,&dfield);
    DataFieldGetEntries(dfield,(void**)&data);
    MaterialConstantsSetDefaultAll_SourceDecay(nregions,data);
  }

  {
    EnergySourceAdiabaticAdvection *data;
    
    DataBucketGetDataFieldByName(db,EnergySourceAdiabaticAdvection_classname,&dfield);
    DataFieldGetEntries(dfield,(void**)&data);
    MaterialConstantsSetDefaultAll_SourceAdiabaticAdv(nregions,data);
  }

  PetscFunctionReturn(0);
}

/* Assume temperature scale = 1 K */
PetscErrorCode MaterialConstantsEnergyScaleAll(DataBucket db,const int region_id,
                                               PetscReal length_scale,
                                               PetscReal time_scale,
                                               PetscReal pressure_scale)
{
  DataField                dfield;
  EnergyMaterialConstants *mdata;
  int st,type;
  PetscReal k_scale,H_scale,Cp_scale,density_scale;
  
  density_scale = pressure_scale * (time_scale*time_scale) / (length_scale*length_scale);
  
  k_scale = pressure_scale * length_scale * length_scale / time_scale;  /* W/(m.K) = kg.m.s^-3.K^-1 = (kg/m^3).m^4.s^-3 */
  H_scale = pressure_scale / time_scale; /* W/m^3 */
  Cp_scale = pressure_scale /(density_scale); /* J / (kg.K) = m^2.kg.s^-2.K^-1 = (kg/m^3).m^5.s^-2 . K^-1 */
  
	DataBucketGetDataFieldByName(db,EnergyMaterialConstants_classname,&dfield);
  DataFieldGetEntries(dfield,(void**)&mdata);
  
  /* NOTE: beta has units of 1/Pa, so scale by reciprocal */
  MaterialConstantsScaleValues_EnergyMaterialConstants(region_id,mdata,
                                                       1.0,                // alpha0
                                                       1.0/pressure_scale, // beta0,
                                                       density_scale,      // rho
                                                       Cp_scale);          // Cp0
 
  
  type = mdata[region_id].density_type;
  switch (type) {

    case ENERGYDENSITY_NONE:
      break;
    
    case ENERGYDENSITY_USE_MATERIALPOINT_VALUE:
      /* NOTE: I expect the user set the correct scaled values of the material point */
      break;
      
    case ENERGYDENSITY_CONSTANT:
      /* Nothing to do: rho is already scaled by MaterialConstantsScaleValues_EnergyMaterialConstants */
      break;
      
    case ENERGYDENSITY_BOUSSINESQ:
      break;
  }
  

  type = mdata[region_id].conductivity_type;
  switch (type) {

    case ENERGYCONDUCTIVITY_NONE:
      break;
    
    case ENERGYCONDUCTIVITY_USE_MATERIALPOINT_VALUE:
      /* NOTE: I expect the user set the correct scaled values of the material point */
      break;
      
    case ENERGYCONDUCTIVITY_CONSTANT: {
      EnergyConductivityConst *data;
      
      DataBucketGetDataFieldByName(db,EnergyConductivityConst_classname,&dfield);
      DataFieldGetEntries(dfield,(void**)&data);
      
      MaterialConstantsScaleValues_ConductivityConst(region_id,data,k_scale);
    }
      break;
      
    case ENERGYCONDUCTIVITY_TEMP_DEP_THRESHOLD: {
    	EnergyConductivityThreshold *data;

    	DataBucketGetDataFieldByName(db,EnergyConductivityThreshold_classname,&dfield);
    	DataFieldGetEntries(dfield,(void**)&data);

    MaterialConstantsScaleValues_ConductivityThreshold(region_id,data,k_scale,k_scale,1.0,1.0);
    }
      break;
  }
  
  for (st=0; st<7; st++) {
    type = mdata[region_id].source_type[st];

    switch (type) {
        
      case ENERGYSOURCE_NONE:
        break;

      case ENERGYSOURCE_USE_MATERIALPOINT_VALUE:
        /* NOTE: I expect the user set the correct scaled values of the material point */
        break;
        
      case ENERGYSOURCE_CONSTANT: {
        EnergySourceConst *data;
        
        DataBucketGetDataFieldByName(db,EnergySourceConst_classname,&dfield);
        DataFieldGetEntries(dfield,(void**)&data);
        
        MaterialConstantsScaleValues_SourceConst(region_id,data,H_scale);
      }
        
        break;
        
      case ENERGYSOURCE_SHEAR_HEATING:
        /* Nothing to scale - all interntal to EnergyEvaluateCoefficients_MaterialPoints */
        break;
        
      case ENERGYSOURCE_DECAY: {
        EnergySourceDecay *data;
        
        DataBucketGetDataFieldByName(db,EnergySourceDecay_classname,&dfield);
        DataFieldGetEntries(dfield,(void**)&data);
        
        /* NOTE: lambda has units of 1/s, so scale by reciprocal */
        MaterialConstantsScaleValues_SourceDecay(region_id,data,H_scale,1.0/time_scale);
      }        break;
        
      case ENERGYSOURCE_ADIABATIC:
        /* Nothing to scale - all interntal to EnergyEvaluateCoefficients_MaterialPoints */
        break;
        
      case ENERGYSOURCE_ADIABATIC_ADVECTION: {
        EnergySourceAdiabaticAdvection *data;
        
        DataBucketGetDataFieldByName(db,EnergySourceAdiabaticAdvection_classname,&dfield);
        DataFieldGetEntries(dfield,(void**)&data);
       
        /* NOTE: dT/dy has units of K/m scale by 1/m */
        MaterialConstantsScaleValues_SourceAdiabaticAdv(region_id,data,1.0/length_scale);
      }
        
        
        break;
    }
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode MaterialConstantsEnergyPrintAll(DataBucket db,const int region_id)
{
  DataField dfield;
  EnergyMaterialConstants *data;
  
	DataBucketGetDataFieldByName(db,EnergyMaterialConstants_classname,&dfield);
  DataFieldGetEntries(dfield,(void**)&data);
  
  MaterialConstantsPrintValues_EnergyMaterialConstants(NULL,region_id,data);
  
  PetscFunctionReturn(0);
}

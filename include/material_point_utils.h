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
 **    filename:   material_point_utils.h
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


#ifndef __ptatin_material_point_utils_h__
#define __ptatin_material_point_utils_h__

#include "data_bucket.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPStokesPl_def.h"
#include "MPntPEnergy_def.h"
#include "quadrature.h"

/* add material points into the list */
typedef enum { MPField_Std=0, MPField_Stokes, MPField_Energy, MPField_StokesPl } MaterialPointField;

typedef enum { CoefAvgNULL=0, CoefAvgARITHMETIC, CoefAvgHARMONIC, CoefAvgGEOMETRIC } CoefficientAveragingType;


typedef struct _p_MPAccess *MPAccess;
struct _p_MPAccess {
	DataBucket db;
	DataField  *PField;
	int        nfields;
	/* USER: add reference to all possible material point types here */
	int mp_std_field_idx;
	int mp_stokes_field_idx;
	int mp_stokespl_field_idx;
	int mp_energy_field_idx;
};


PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,DataBucket db,int *byte_offset,const int nfields,const MaterialPointField list[]);
PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,DataBucket db,const int nfields,const MaterialPointField list[]);
PetscErrorCode MaterialPointGeneric_PVTUWriteAllPPointDataFields(FILE *vtk_fp,const int nfields,const MaterialPointField list[]);

PetscErrorCode SwarmViewGeneric_VTUXML_binary_appended(DataBucket db,const int nfields,const MaterialPointField list[],const char name[]);
PetscErrorCode SwarmViewGeneric_PVTUXML(const int nfields,const MaterialPointField list[],const char prefix[],const char name[]);
PetscErrorCode SwarmViewGeneric_ParaView(DataBucket db,const int nfields,const MaterialPointField list[],const char path[],const char prefix[]);

/* projections [fine grid] */
PetscErrorCode MPntPStokesProj_P0(CoefficientAveragingType type,const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],DM da,Quadrature Q);
PetscErrorCode SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],DM da,Quadrature Q);
PetscErrorCode SwarmUpdateGaussPropertiesOne2OneMap_MPntPStokes(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],Quadrature Q);

/* projection [for levels in a hierarchy] */
PetscErrorCode SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(PetscInt coefficient_projection_type,const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,Mat R[],DM da[],Quadrature Q[]);

PetscErrorCode MProjection_P0Projection_onto_Q2_MPntPStokes_Level(CoefficientAveragingType eta_type,CoefficientAveragingType rho_type,const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,DM da[],PetscInt level,Quadrature Q_level);


/* depreciated */
PetscErrorCode MaterialPointQuadraturePointProjectionC0_Q2Stokes(DM da,DataBucket materialpoint_db,MaterialPointField field,const int member,Quadrature Q);

PetscErrorCode MProjection_Q1Projection_onto_Q2_MPntPStokes_Level(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,DM da[],PetscInt level,Quadrature Q_level);

PetscErrorCode _MaterialPointProjection_MapOntoQ2Mesh(
                                                      DM clone,Vec properties_A,Vec properties_B,CoefficientAveragingType avg_type,
                                                      const int npoints,MPntStd mp_std[],
                                                      size_t member_offset,size_t point_offset,void *point_data);
PetscErrorCode _MaterialPointProjection_MapOntoQ2Mesh_InterpolateToQuadraturePoint(
                                                                                   DM clone,Vec properties_A,
                                                                                   size_t member_offset,size_t qpoint_offset,void *qpoint_data,Quadrature Q) ;

PetscErrorCode DMDAEQ1_MaterialPointProjection_MapOntoQ2Mesh(
                                                             DM clone,Vec properties_A,Vec properties_B,CoefficientAveragingType avg_type,
                                                             const int npoints,MPntStd mp_std[],
                                                             size_t member_offset,size_t point_offset,void *point_data);
PetscErrorCode DMDAEQ1_MaterialPointProjection_MapOntoQ2Mesh_InterpolateToQuadraturePoint(
                                                                                          DM clone,Vec properties_A,
                                                                                          size_t member_offset,size_t qpoint_offset,void *qpoint_data,Quadrature Q) ;




PetscErrorCode MPntPStokesPlComputeMemberOffsets(size_t property_offsets[]);
PetscErrorCode MPntPEnergyComputeMemberOffsets(size_t property_offsets[]);
PetscErrorCode QPntVolCoefStokesComputeMemberOffsets(size_t property_offsets[]);
PetscErrorCode QPntVolCoefEnergyComputeMemberOffsets(size_t property_offsets[]);

PetscErrorCode MaterialPointGetAccess(DataBucket materialpoint_db,MPAccess *helper);
PetscErrorCode MaterialPointRestoreAccess(DataBucket matpoint_db,MPAccess *helper);

PetscErrorCode MaterialPointGet_global_coord(MPAccess X,const int p,double *var[]);
PetscErrorCode MaterialPointGet_local_coord(MPAccess X,const int p,double *var[]);
PetscErrorCode MaterialPointGet_local_element_index(MPAccess X,const int p,int *var);
PetscErrorCode MaterialPointGet_phase_index(MPAccess X,const int p,int *var);

PetscErrorCode MaterialPointGet_viscosity(MPAccess X,const int p,double *var);
PetscErrorCode MaterialPointGet_density(MPAccess X,const int p,double *var);
PetscErrorCode MaterialPointGet_plastic_strain(MPAccess X,const int p,float *var);
PetscErrorCode MaterialPointGet_yield_indicator(MPAccess X,const int p,short *var);
PetscErrorCode MaterialPointGet_diffusivity(MPAccess X,const int p,double *var);
PetscErrorCode MaterialPointGet_heat_source(MPAccess X,const int p,double *var);

PetscErrorCode MaterialPointSet_global_coord(MPAccess X,const int p,double var[]);
PetscErrorCode MaterialPointSet_local_coord(MPAccess X,const int p,double var[]);
PetscErrorCode MaterialPointSet_local_element_index(MPAccess X,const int p,int var);
PetscErrorCode MaterialPointSet_phase_index(MPAccess X,const int p,int var);
PetscErrorCode MaterialPointSet_viscosity(MPAccess X,const int p,double var);
PetscErrorCode MaterialPointSet_density(MPAccess X,const int p,double var);
PetscErrorCode MaterialPointSet_plastic_strain(MPAccess X,const int p,float var);
PetscErrorCode MaterialPointSet_yield_indicator(MPAccess X,const int p,short var);
PetscErrorCode MaterialPointSet_diffusivity(MPAccess X,const int p,double var);
PetscErrorCode MaterialPointSet_heat_source(MPAccess X,const int p,double var);

PetscErrorCode MaterialPointScale_global_coord(MPAccess X,double var);
PetscErrorCode MaterialPointScale_viscosity(MPAccess X,double var);
PetscErrorCode MaterialPointScale_density(MPAccess X,double var);
PetscErrorCode MaterialPointScale_plastic_strain(MPAccess X,double var);
PetscErrorCode MaterialPointScale_diffusivity(MPAccess X,double var);
PetscErrorCode MaterialPointScale_heat_source(MPAccess X,double var);

#endif


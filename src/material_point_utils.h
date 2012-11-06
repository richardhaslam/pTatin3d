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
 **    Filename:      material_point_utils.h
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


#ifndef __ptatin_material_point_utils_h__
#define __ptatin_material_point_utils_h__

#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "quadrature.h"

/* add material points into the list */
typedef enum { MPField_Std=0, MPField_Stokes } MaterialPointField;

PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,DataBucket db,int *byte_offset,const int nfields,const MaterialPointField list[]);
PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,DataBucket db,const int nfields,const MaterialPointField list[]);
PetscErrorCode MaterialPointGeneric_PVTUWriteAllPPointDataFields(FILE *vtk_fp,const int nfields,const MaterialPointField list[]);

PetscErrorCode SwarmViewGeneric_VTUXML_binary_appended(DataBucket db,const int nfields,const MaterialPointField list[],const char name[]);
PetscErrorCode SwarmViewGeneric_PVTUXML(const int nfields,const MaterialPointField list[],const char prefix[],const char name[]);
PetscErrorCode SwarmViewGeneric_ParaView(DataBucket db,const int nfields,const MaterialPointField list[],const char path[],const char prefix[]);

PetscErrorCode SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],DM da,Quadrature Q);
PetscErrorCode SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,Mat R[],DM da[],Quadrature Q[]);


#endif



#ifndef __ptatin_material_point_utils_h__
#define __ptatin_material_point_utils_h__

#include "swarm_fields.h"

/* add material points into the list */
typedef enum { MPField_Std=0, MPField_Stokes } MaterialPointField;

PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,DataBucket db,int *byte_offset,const int nfields,const MaterialPointField list[]);
PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,DataBucket db,const int nfields,const MaterialPointField list[]);
PetscErrorCode MaterialPointGeneric_PVTUWriteAllPPointDataFields(FILE *vtk_fp,const int nfields,const MaterialPointField list[]);

PetscErrorCode SwarmViewGeneric_VTUXML_binary_appended(DataBucket db,const int nfields,const MaterialPointField list[],const char name[]);
PetscErrorCode SwarmViewGeneric_PVTUXML(const int nfields,const MaterialPointField list[],const char prefix[],const char name[]);
PetscErrorCode SwarmViewGeneric_ParaView(DataBucket db,const int nfields,const MaterialPointField list[],const char path[],const char prefix[]);

#endif

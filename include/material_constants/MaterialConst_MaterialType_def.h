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
 **    filename:   MaterialConst_MaterialType_def.h
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
/*
  Auto generated by version 0.0 of swarm_class_generator.py
  on les-MacBook-Pro.local, at 2013-03-08 17:43:59.381411 by laetitia
*/


#ifndef __MaterialConst_MaterialType_DEF_H__
#define __MaterialConst_MaterialType_DEF_H__

typedef struct {
  int visc_type ;
  int plastic_type ;
  int softening_type ;
  int density_type ;
} MaterialConst_MaterialType ;


typedef enum {
  MaterialType_visc_type = 0,
  MaterialType_plastic_type,
  MaterialType_softening_type,
  MaterialType_density_type
} MaterialConst_MaterialTypeTypeName ;


extern const char MaterialConst_MaterialType_classname[];

extern const int MaterialConst_MaterialType_nmembers;

extern const size_t MaterialConst_MaterialType_member_sizes[];

extern const char *MaterialConst_MaterialType_member_names[];

/* prototypes */
void MaterialConst_MaterialTypeGetField_visc_type(MaterialConst_MaterialType *point,int *data);
void MaterialConst_MaterialTypeGetField_plastic_type(MaterialConst_MaterialType *point,int *data);
void MaterialConst_MaterialTypeGetField_softening_type(MaterialConst_MaterialType *point,int *data);
void MaterialConst_MaterialTypeGetField_density_type(MaterialConst_MaterialType *point,int *data);
void MaterialConst_MaterialTypeSetField_visc_type(MaterialConst_MaterialType *point,int data);
void MaterialConst_MaterialTypeSetField_plastic_type(MaterialConst_MaterialType *point,int data);
void MaterialConst_MaterialTypeSetField_softening_type(MaterialConst_MaterialType *point,int data);
void MaterialConst_MaterialTypeSetField_density_type(MaterialConst_MaterialType *point,int data);
void MaterialConst_MaterialTypeView(MaterialConst_MaterialType *point);
void MaterialConst_MaterialTypeVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const MaterialConst_MaterialType points[]);
void MaterialConst_MaterialTypePVTUWriteAllPPointDataFields(FILE *vtk_fp);
void MaterialConst_MaterialTypeVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const MaterialConst_MaterialType points[]);
void MaterialConst_MaterialTypeVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const MaterialConst_MaterialType points[]);

#endif

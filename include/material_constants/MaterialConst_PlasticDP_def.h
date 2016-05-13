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
 **    filename:   MaterialConst_PlasticDP_def.h
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
  on les-MacBook-Pro.local, at 2013-03-08 17:43:59.346273 by laetitia
*/


#ifndef __MaterialConst_PlasticDP_DEF_H__
#define __MaterialConst_PlasticDP_DEF_H__

typedef struct {
  double phi ;
  double Co ;
  double phi_inf ;
  double Co_inf ;
  double tens_cutoff ;
  double hst_cutoff ;
} MaterialConst_PlasticDP ;


typedef enum {
  PlasticDP_friction = 0,
  PlasticDP_cohesion,
  PlasticDP_friction_inf,
  PlasticDP_cohesion_inf,
  PlasticDP_tens_cutoff,
  PlasticDP_hst_cutoff
} MaterialConst_PlasticDPTypeName ;


extern const char MaterialConst_PlasticDP_classname[];

extern const int MaterialConst_PlasticDP_nmembers;

extern const size_t MaterialConst_PlasticDP_member_sizes[];

extern const char *MaterialConst_PlasticDP_member_names[];

/* prototypes */
void MaterialConst_PlasticDPGetField_friction(MaterialConst_PlasticDP *point,double *data);
void MaterialConst_PlasticDPGetField_cohesion(MaterialConst_PlasticDP *point,double *data);
void MaterialConst_PlasticDPGetField_friction_inf(MaterialConst_PlasticDP *point,double *data);
void MaterialConst_PlasticDPGetField_cohesion_inf(MaterialConst_PlasticDP *point,double *data);
void MaterialConst_PlasticDPGetField_tens_cutoff(MaterialConst_PlasticDP *point,double *data);
void MaterialConst_PlasticDPGetField_hst_cutoff(MaterialConst_PlasticDP *point,double *data);
void MaterialConst_PlasticDPSetField_friction(MaterialConst_PlasticDP *point,double data);
void MaterialConst_PlasticDPSetField_cohesion(MaterialConst_PlasticDP *point,double data);
void MaterialConst_PlasticDPSetField_friction_inf(MaterialConst_PlasticDP *point,double data);
void MaterialConst_PlasticDPSetField_cohesion_inf(MaterialConst_PlasticDP *point,double data);
void MaterialConst_PlasticDPSetField_tens_cutoff(MaterialConst_PlasticDP *point,double data);
void MaterialConst_PlasticDPSetField_hst_cutoff(MaterialConst_PlasticDP *point,double data);
void MaterialConst_PlasticDPView(MaterialConst_PlasticDP *point);
void MaterialConst_PlasticDPVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const MaterialConst_PlasticDP points[]);
void MaterialConst_PlasticDPPVTUWriteAllPPointDataFields(FILE *vtk_fp);
void MaterialConst_PlasticDPVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const MaterialConst_PlasticDP points[]);
void MaterialConst_PlasticDPVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const MaterialConst_PlasticDP points[]);

#endif
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
 **    filename:   MaterialConst_MaterialType_def.c
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
  on les-MacBook-Pro.local, at 2013-03-08 17:43:59.382972 by laetitia
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "material_constants/MaterialConst_MaterialType_def.h"


const char MaterialConst_MaterialType_classname[] = "MaterialConst_MaterialType";

const int MaterialConst_MaterialType_nmembers = 4;

const size_t MaterialConst_MaterialType_member_sizes[] = {
  1 * sizeof(int),
  1 * sizeof(int),
  1 * sizeof(int),
  1 * sizeof(int)
} ;

const char *MaterialConst_MaterialType_member_names[] = {
  "visc_type",
  "plastic_type",
  "softening_type",
  "density_type"
} ;


/* ===================================== */
/* Getters for MaterialConst_MaterialType */
/* ===================================== */
void MaterialConst_MaterialTypeGetField_visc_type(MaterialConst_MaterialType *point,int *data) 
{
  *data = point->visc_type;
}

void MaterialConst_MaterialTypeGetField_plastic_type(MaterialConst_MaterialType *point,int *data) 
{
  *data = point->plastic_type;
}

void MaterialConst_MaterialTypeGetField_softening_type(MaterialConst_MaterialType *point,int *data) 
{
  *data = point->softening_type;
}

void MaterialConst_MaterialTypeGetField_density_type(MaterialConst_MaterialType *point,int *data) 
{
  *data = point->density_type;
}


/* ===================================== */
/* Setters for MaterialConst_MaterialType */
/* ===================================== */
void MaterialConst_MaterialTypeSetField_visc_type(MaterialConst_MaterialType *point,int data) 
{
  point->visc_type = data;
}

void MaterialConst_MaterialTypeSetField_plastic_type(MaterialConst_MaterialType *point,int data) 
{
  point->plastic_type = data;
}

void MaterialConst_MaterialTypeSetField_softening_type(MaterialConst_MaterialType *point,int data) 
{
  point->softening_type = data;
}

void MaterialConst_MaterialTypeSetField_density_type(MaterialConst_MaterialType *point,int data) 
{
  point->density_type = data;
}


/* ===================================== */
/* C-viewer for MaterialConst_MaterialType */
/* ===================================== */
void MaterialConst_MaterialTypeView(MaterialConst_MaterialType *point)
{
  {
    int data;
    MaterialConst_MaterialTypeGetField_visc_type(point,&data);
    printf("field: visc_type = %d; [size %zu; type int; variable_name visc_type]\n",data, MaterialConst_MaterialType_member_sizes[0] );
  }
  {
    int data;
    MaterialConst_MaterialTypeGetField_plastic_type(point,&data);
    printf("field: plastic_type = %d; [size %zu; type int; variable_name plastic_type]\n",data, MaterialConst_MaterialType_member_sizes[1] );
  }
  {
    int data;
    MaterialConst_MaterialTypeGetField_softening_type(point,&data);
    printf("field: softening_type = %d; [size %zu; type int; variable_name softening_type]\n",data, MaterialConst_MaterialType_member_sizes[2] );
  }
  {
    int data;
    MaterialConst_MaterialTypeGetField_density_type(point,&data);
    printf("field: density_type = %d; [size %zu; type int; variable_name density_type]\n",data, MaterialConst_MaterialType_member_sizes[3] );
  }
}


/* ===================================== */
/* VTK viewer for MaterialConst_MaterialType */
/* ===================================== */
void MaterialConst_MaterialTypeVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const MaterialConst_MaterialType points[]) 
{
  int p;
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"visc_type\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%d\n",(int)points[p].visc_type);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"plastic_type\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%d\n",(int)points[p].plastic_type);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"softening_type\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%d\n",(int)points[p].softening_type);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"density_type\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%d\n",(int)points[p].density_type);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
}


/* ===================================== */
/* PVTU viewer for MaterialConst_MaterialType */
/* ===================================== */
void MaterialConst_MaterialTypePVTUWriteAllPPointDataFields(FILE *vtk_fp) 
{
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Int32\" Name=\"visc_type\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Int32\" Name=\"plastic_type\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Int32\" Name=\"softening_type\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Int32\" Name=\"density_type\" NumberOfComponents=\"1\"/>\n");
}


/* ===================================== */
/* VTK binary (appended header) viewer for MaterialConst_MaterialType */
/* ===================================== */
void MaterialConst_MaterialTypeVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const MaterialConst_MaterialType points[]) 
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"visc_type\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(int);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"plastic_type\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(int);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"softening_type\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(int);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"density_type\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(int);

}


/* ================================================== */
/* VTK binary (appended data) viewer for MaterialConst_MaterialType */
/* ==================================================== */
void MaterialConst_MaterialTypeVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const MaterialConst_MaterialType points[]) 
{
  int p,length;
  size_t atomic_size;

  atomic_size = sizeof(int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].visc_type,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].plastic_type,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].softening_type,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].density_type,atomic_size,1,vtk_fp);
  }

}


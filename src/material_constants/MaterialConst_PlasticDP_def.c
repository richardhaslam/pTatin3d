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
 **    filename:   MaterialConst_PlasticDP_def.c
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
  on les-MacBook-Pro.local, at 2013-03-08 17:43:59.346464 by laetitia
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "material_constants/MaterialConst_PlasticDP_def.h"


const char MaterialConst_PlasticDP_classname[] = "MaterialConst_PlasticDP";

const int MaterialConst_PlasticDP_nmembers = 6;

const size_t MaterialConst_PlasticDP_member_sizes[] = {
  1 * sizeof(double),
  1 * sizeof(double),
  1 * sizeof(double),
  1 * sizeof(double),
  1 * sizeof(double),
  1 * sizeof(double)
} ;

const char *MaterialConst_PlasticDP_member_names[] = {
  "friction",
  "cohesion",
  "friction_inf",
  "cohesion_inf",
  "tens_cutoff",
  "hst_cutoff"
} ;


/* ===================================== */
/* Getters for MaterialConst_PlasticDP */
/* ===================================== */
void MaterialConst_PlasticDPGetField_friction(MaterialConst_PlasticDP *point,double *data)
{
  *data = point->phi;
}

void MaterialConst_PlasticDPGetField_cohesion(MaterialConst_PlasticDP *point,double *data)
{
  *data = point->Co;
}

void MaterialConst_PlasticDPGetField_friction_inf(MaterialConst_PlasticDP *point,double *data)
{
  *data = point->phi_inf;
}

void MaterialConst_PlasticDPGetField_cohesion_inf(MaterialConst_PlasticDP *point,double *data)
{
  *data = point->Co_inf;
}

void MaterialConst_PlasticDPGetField_tens_cutoff(MaterialConst_PlasticDP *point,double *data)
{
  *data = point->tens_cutoff;
}

void MaterialConst_PlasticDPGetField_hst_cutoff(MaterialConst_PlasticDP *point,double *data)
{
  *data = point->hst_cutoff;
}


/* ===================================== */
/* Setters for MaterialConst_PlasticDP */
/* ===================================== */
void MaterialConst_PlasticDPSetField_friction(MaterialConst_PlasticDP *point,double data)
{
  point->phi = data;
}

void MaterialConst_PlasticDPSetField_cohesion(MaterialConst_PlasticDP *point,double data)
{
  point->Co = data;
}

void MaterialConst_PlasticDPSetField_friction_inf(MaterialConst_PlasticDP *point,double data)
{
  point->phi_inf = data;
}

void MaterialConst_PlasticDPSetField_cohesion_inf(MaterialConst_PlasticDP *point,double data)
{
  point->Co_inf = data;
}

void MaterialConst_PlasticDPSetField_tens_cutoff(MaterialConst_PlasticDP *point,double data)
{
  point->tens_cutoff = data;
}

void MaterialConst_PlasticDPSetField_hst_cutoff(MaterialConst_PlasticDP *point,double data)
{
  point->hst_cutoff = data;
}


/* ===================================== */
/* C-viewer for MaterialConst_PlasticDP */
/* ===================================== */
void MaterialConst_PlasticDPView(MaterialConst_PlasticDP *point)
{
  {
    double data;
    MaterialConst_PlasticDPGetField_friction(point,&data);
    printf("field: friction = %1.6e; [size %zu; type double; variable_name phi]\n",data, MaterialConst_PlasticDP_member_sizes[0] );
  }
  {
    double data;
    MaterialConst_PlasticDPGetField_cohesion(point,&data);
    printf("field: cohesion = %1.6e; [size %zu; type double; variable_name Co]\n",data, MaterialConst_PlasticDP_member_sizes[1] );
  }
  {
    double data;
    MaterialConst_PlasticDPGetField_friction_inf(point,&data);
    printf("field: friction_inf = %1.6e; [size %zu; type double; variable_name phi_inf]\n",data, MaterialConst_PlasticDP_member_sizes[2] );
  }
  {
    double data;
    MaterialConst_PlasticDPGetField_cohesion_inf(point,&data);
    printf("field: cohesion_inf = %1.6e; [size %zu; type double; variable_name Co_inf]\n",data, MaterialConst_PlasticDP_member_sizes[3] );
  }
  {
    double data;
    MaterialConst_PlasticDPGetField_tens_cutoff(point,&data);
    printf("field: tens_cutoff = %1.6e; [size %zu; type double; variable_name tens_cutoff]\n",data, MaterialConst_PlasticDP_member_sizes[4] );
  }
  {
    double data;
    MaterialConst_PlasticDPGetField_hst_cutoff(point,&data);
    printf("field: hst_cutoff = %1.6e; [size %zu; type double; variable_name hst_cutoff]\n",data, MaterialConst_PlasticDP_member_sizes[5] );
  }
}


/* ===================================== */
/* VTK viewer for MaterialConst_PlasticDP */
/* ===================================== */
void MaterialConst_PlasticDPVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const MaterialConst_PlasticDP points[])
{
  int p;
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"phi\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].phi);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Co\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].Co);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"phi_inf\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].phi_inf);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Co_inf\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].Co_inf);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"tens_cutoff\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].tens_cutoff);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"hst_cutoff\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].hst_cutoff);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
}


/* ===================================== */
/* PVTU viewer for MaterialConst_PlasticDP */
/* ===================================== */
void MaterialConst_PlasticDPPVTUWriteAllPPointDataFields(FILE *vtk_fp)
{
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"phi\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"Co\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"phi_inf\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"Co_inf\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"tens_cutoff\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"hst_cutoff\" NumberOfComponents=\"1\"/>\n");
}


/* ===================================== */
/* VTK binary (appended header) viewer for MaterialConst_PlasticDP */
/* ===================================== */
void MaterialConst_PlasticDPVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const MaterialConst_PlasticDP points[])
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"phi\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Co\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"phi_inf\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Co_inf\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"tens_cutoff\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"hst_cutoff\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

}


/* ================================================== */
/* VTK binary (appended data) viewer for MaterialConst_PlasticDP */
/* ==================================================== */
void MaterialConst_PlasticDPVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const MaterialConst_PlasticDP points[])
{
  int p,length;
  size_t atomic_size;

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].phi,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].Co,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].phi_inf,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].Co_inf,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].tens_cutoff,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].hst_cutoff,atomic_size,1,vtk_fp);
  }

}


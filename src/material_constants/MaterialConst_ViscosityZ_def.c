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
 **    filename:   MaterialConst_ViscosityZ_def.c
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
  on les-MacBook-Pro.local, at 2013-03-08 17:43:59.343799 by laetitia
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "material_constants/MaterialConst_ViscosityZ_def.h"


const char MaterialConst_ViscosityZ_classname[] = "MaterialConst_ViscosityZ";

const int MaterialConst_ViscosityZ_nmembers = 3;

const size_t MaterialConst_ViscosityZ_member_sizes[] = {
  1 * sizeof(double),
  1 * sizeof(double),
  1 * sizeof(double)
} ;

const char *MaterialConst_ViscosityZ_member_names[] = {
  "eta0",
  "zeta",
  "zref"
} ;


/* ===================================== */
/* Getters for MaterialConst_ViscosityZ */
/* ===================================== */
void MaterialConst_ViscosityZGetField_eta0(MaterialConst_ViscosityZ *point,double *data)
{
  *data = point->eta0;
}

void MaterialConst_ViscosityZGetField_zeta(MaterialConst_ViscosityZ *point,double *data)
{
  *data = point->zeta;
}

void MaterialConst_ViscosityZGetField_zref(MaterialConst_ViscosityZ *point,double *data)
{
  *data = point->zref;
}


/* ===================================== */
/* Setters for MaterialConst_ViscosityZ */
/* ===================================== */
void MaterialConst_ViscosityZSetField_eta0(MaterialConst_ViscosityZ *point,double data)
{
  point->eta0 = data;
}

void MaterialConst_ViscosityZSetField_zeta(MaterialConst_ViscosityZ *point,double data)
{
  point->zeta = data;
}

void MaterialConst_ViscosityZSetField_zref(MaterialConst_ViscosityZ *point,double data)
{
  point->zref = data;
}


/* ===================================== */
/* C-viewer for MaterialConst_ViscosityZ */
/* ===================================== */
void MaterialConst_ViscosityZView(MaterialConst_ViscosityZ *point)
{
  {
    double data;
    MaterialConst_ViscosityZGetField_eta0(point,&data);
    printf("field: eta0 = %1.6e; [size %zu; type double; variable_name eta0]\n",data, MaterialConst_ViscosityZ_member_sizes[0] );
  }
  {
    double data;
    MaterialConst_ViscosityZGetField_zeta(point,&data);
    printf("field: zeta = %1.6e; [size %zu; type double; variable_name zeta]\n",data, MaterialConst_ViscosityZ_member_sizes[1] );
  }
  {
    double data;
    MaterialConst_ViscosityZGetField_zref(point,&data);
    printf("field: zref = %1.6e; [size %zu; type double; variable_name zref]\n",data, MaterialConst_ViscosityZ_member_sizes[2] );
  }
}


/* ===================================== */
/* VTK viewer for MaterialConst_ViscosityZ */
/* ===================================== */
void MaterialConst_ViscosityZVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const MaterialConst_ViscosityZ points[])
{
  int p;
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"eta0\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].eta0);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"zeta\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].zeta);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"zref\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].zref);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
}


/* ===================================== */
/* PVTU viewer for MaterialConst_ViscosityZ */
/* ===================================== */
void MaterialConst_ViscosityZPVTUWriteAllPPointDataFields(FILE *vtk_fp)
{
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"eta0\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"zeta\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"zref\" NumberOfComponents=\"1\"/>\n");
}


/* ===================================== */
/* VTK binary (appended header) viewer for MaterialConst_ViscosityZ */
/* ===================================== */
void MaterialConst_ViscosityZVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const MaterialConst_ViscosityZ points[])
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"eta0\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"zeta\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"zref\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

}


/* ================================================== */
/* VTK binary (appended data) viewer for MaterialConst_ViscosityZ */
/* ==================================================== */
void MaterialConst_ViscosityZVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const MaterialConst_ViscosityZ points[])
{
  int p,length;
  size_t atomic_size;

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].eta0,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].zeta,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].zref,atomic_size,1,vtk_fp);
  }

}


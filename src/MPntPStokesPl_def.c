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
 **    filename:   MPntPStokesPl_def.c
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
  on musashi.ethz.ch, at 2013-02-27 18:40:06.481385 by dmay
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "MPntPStokesPl_def.h"


const char MPntPStokesPl_classname[] = "MPntPStokesPl";

const int MPntPStokesPl_nmembers = 2;

const size_t MPntPStokesPl_member_sizes[] = {
  1 * sizeof(float),
  1 * sizeof(short)
} ;

const char *MPntPStokesPl_member_names[] = {
  "plastic_strain",
  "yield_indicator"
} ;


/* ===================================== */
/* Getters for MPntPStokesPl */
/* ===================================== */
void MPntPStokesPlGetField_plastic_strain(MPntPStokesPl *point,float *data) 
{
  *data = point->e_plastic;
}

void MPntPStokesPlGetField_yield_indicator(MPntPStokesPl *point,short *data) 
{
  *data = point->is_yielding;
}


/* ===================================== */
/* Setters for MPntPStokesPl */
/* ===================================== */
void MPntPStokesPlSetField_plastic_strain(MPntPStokesPl *point,float data) 
{
  point->e_plastic = data;
}

void MPntPStokesPlSetField_yield_indicator(MPntPStokesPl *point,short data) 
{
  point->is_yielding = data;
}


/* ===================================== */
/* C-viewer for MPntPStokesPl */
/* ===================================== */
void MPntPStokesPlView(MPntPStokesPl *point)
{
  {
    float data;
    MPntPStokesPlGetField_plastic_strain(point,&data);
    printf("field: plastic_strain = %1.6e; [size %zu; type float; variable_name e_plastic]\n",data, MPntPStokesPl_member_sizes[0] );
  }
  {
    short data;
    MPntPStokesPlGetField_yield_indicator(point,&data);
    printf("field: yield_indicator = %d; [size %zu; type short; variable_name is_yielding]\n",data, MPntPStokesPl_member_sizes[1] );
  }
}


/* ===================================== */
/* VTK viewer for MPntPStokesPl */
/* ===================================== */
void MPntPStokesPlVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const MPntPStokesPl points[]) 
{
  int p;
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"e_plastic\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%f\n",(float)points[p].e_plastic);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int16\" Name=\"is_yielding\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%d\n",(short)points[p].is_yielding);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
}


/* ===================================== */
/* PVTU viewer for MPntPStokesPl */
/* ===================================== */
void MPntPStokesPlPVTUWriteAllPPointDataFields(FILE *vtk_fp) 
{
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float32\" Name=\"e_plastic\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Int16\" Name=\"is_yielding\" NumberOfComponents=\"1\"/>\n");
}


/* ===================================== */
/* VTK binary (appended header) viewer for MPntPStokesPl */
/* ===================================== */
void MPntPStokesPlVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const MPntPStokesPl points[]) 
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"e_plastic\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(float);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int16\" Name=\"is_yielding\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(short);

}


/* ================================================== */
/* VTK binary (appended data) viewer for MPntPStokesPl */
/* ==================================================== */
void MPntPStokesPlVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const MPntPStokesPl points[]) 
{
  int p,length;
  size_t atomic_size;

  atomic_size = sizeof(float);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].e_plastic,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(short);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].is_yielding,atomic_size,1,vtk_fp);
  }

}


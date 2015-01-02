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
 **    filename:   QPntSurfCoefStokes_def.c
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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "QPntSurfCoefStokes_def.h"


const char QPntSurfCoefStokes_classname[] = "QPntSurfCoefStokes";

const int QPntSurfCoefStokes_nmembers = 6;

const size_t QPntSurfCoefStokes_member_sizes[] = {
  3 * sizeof(double),
  3 * sizeof(double),
  3 * sizeof(double),
  3 * sizeof(double),
  1 * sizeof(double),
  1 * sizeof(double)
} ;

const char *QPntSurfCoefStokes_member_names[] = {
  "surface_normal",
  "surface_tangent1",
  "surface_tangent2",
  "surface_traction",
  "viscosity",
  "density"
} ;


/* ===================================== */
/* Getters for QPntSurfCoefStokes */
/* ===================================== */
void QPntSurfCoefStokesGetField_surface_normal(QPntSurfCoefStokes *point,double *data[]) 
{
  *data = point->normal;
}

void QPntSurfCoefStokesGetField_surface_tangent1(QPntSurfCoefStokes *point,double *data[]) 
{
  *data = point->tangent1;
}

void QPntSurfCoefStokesGetField_surface_tangent2(QPntSurfCoefStokes *point,double *data[]) 
{
  *data = point->tangent2;
}

void QPntSurfCoefStokesGetField_surface_traction(QPntSurfCoefStokes *point,double *data[]) 
{
  *data = point->traction;
}

void QPntSurfCoefStokesGetField_viscosity(QPntSurfCoefStokes *point,double *data) 
{
  *data = point->eta;
}

void QPntSurfCoefStokesGetField_density(QPntSurfCoefStokes *point,double *data) 
{
  *data = point->rho;
}


/* ===================================== */
/* Setters for QPntSurfCoefStokes */
/* ===================================== */
void QPntSurfCoefStokesSetField_surface_normal(QPntSurfCoefStokes *point,double data[]) 
{
  memcpy( &point->normal[0], data, sizeof(double)*3 );
}

void QPntSurfCoefStokesSetField_surface_tangent1(QPntSurfCoefStokes *point,double data[]) 
{
  memcpy( &point->tangent1[0], data, sizeof(double)*3 );
}

void QPntSurfCoefStokesSetField_surface_tangent2(QPntSurfCoefStokes *point,double data[]) 
{
  memcpy( &point->tangent2[0], data, sizeof(double)*3 );
}

void QPntSurfCoefStokesSetField_surface_traction(QPntSurfCoefStokes *point,double data[]) 
{
  memcpy( &point->traction[0], data, sizeof(double)*3 );
}

void QPntSurfCoefStokesSetField_viscosity(QPntSurfCoefStokes *point,double data) 
{
  point->eta = data;
}

void QPntSurfCoefStokesSetField_density(QPntSurfCoefStokes *point,double data) 
{
  point->rho = data;
}


/* ===================================== */
/* C-viewer for QPntSurfCoefStokes */
/* ===================================== */
void QPntSurfCoefStokesView(QPntSurfCoefStokes *point)
{
  {
    double *data;
    QPntSurfCoefStokesGetField_surface_normal(point,&data);
    printf("field: surface_normal[0] = %1.6e; [size %zu; type double; variable_name normal]\n",data[0], QPntSurfCoefStokes_member_sizes[0] );
    printf("field: surface_normal[1] = %1.6e; [size %zu; type double; variable_name normal]\n",data[1], QPntSurfCoefStokes_member_sizes[0] );
    printf("field: surface_normal[2] = %1.6e; [size %zu; type double; variable_name normal]\n",data[2], QPntSurfCoefStokes_member_sizes[0] );
  }
  {
    double *data;
    QPntSurfCoefStokesGetField_surface_tangent1(point,&data);
    printf("field: surface_tangent1[0] = %1.6e; [size %zu; type double; variable_name tangent1]\n",data[0], QPntSurfCoefStokes_member_sizes[1] );
    printf("field: surface_tangent1[1] = %1.6e; [size %zu; type double; variable_name tangent1]\n",data[1], QPntSurfCoefStokes_member_sizes[1] );
    printf("field: surface_tangent1[2] = %1.6e; [size %zu; type double; variable_name tangent1]\n",data[2], QPntSurfCoefStokes_member_sizes[1] );
  }
  {
    double *data;
    QPntSurfCoefStokesGetField_surface_tangent2(point,&data);
    printf("field: surface_tangent2[0] = %1.6e; [size %zu; type double; variable_name tangent2]\n",data[0], QPntSurfCoefStokes_member_sizes[2] );
    printf("field: surface_tangent2[1] = %1.6e; [size %zu; type double; variable_name tangent2]\n",data[1], QPntSurfCoefStokes_member_sizes[2] );
    printf("field: surface_tangent2[2] = %1.6e; [size %zu; type double; variable_name tangent2]\n",data[2], QPntSurfCoefStokes_member_sizes[2] );
  }
  {
    double *data;
    QPntSurfCoefStokesGetField_surface_traction(point,&data);
    printf("field: surface_traction[0] = %1.6e; [size %zu; type double; variable_name traction]\n",data[0], QPntSurfCoefStokes_member_sizes[3] );
    printf("field: surface_traction[1] = %1.6e; [size %zu; type double; variable_name traction]\n",data[1], QPntSurfCoefStokes_member_sizes[3] );
    printf("field: surface_traction[2] = %1.6e; [size %zu; type double; variable_name traction]\n",data[2], QPntSurfCoefStokes_member_sizes[3] );
  }
  {
    double data;
    QPntSurfCoefStokesGetField_viscosity(point,&data);
    printf("field: viscosity = %1.6e; [size %zu; type double; variable_name eta]\n",data, QPntSurfCoefStokes_member_sizes[4] );
  }
  {
    double data;
    QPntSurfCoefStokesGetField_density(point,&data);
    printf("field: density = %1.6e; [size %zu; type double; variable_name rho]\n",data, QPntSurfCoefStokes_member_sizes[5] );
  }
}


/* ===================================== */
/* VTK viewer for QPntSurfCoefStokes */
/* ===================================== */
void QPntSurfCoefStokesVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const QPntSurfCoefStokes points[]) 
{
  int p;
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"eta\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].eta);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"rho\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].rho);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
}


/* ===================================== */
/* PVTU viewer for QPntSurfCoefStokes */
/* ===================================== */
void QPntSurfCoefStokesPVTUWriteAllPPointDataFields(FILE *vtk_fp) 
{
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"eta\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"rho\" NumberOfComponents=\"1\"/>\n");
}


/* ===================================== */
/* VTK binary (appended header) viewer for QPntSurfCoefStokes */
/* ===================================== */
void QPntSurfCoefStokesVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const QPntSurfCoefStokes points[]) 
{
  /* Warning: swarm_class_generator.py is ignoring multi-component field normal[] */

  /* Warning: swarm_class_generator.py is ignoring multi-component field tangent1[] */

  /* Warning: swarm_class_generator.py is ignoring multi-component field tangent2[] */

  /* Warning: swarm_class_generator.py is ignoring multi-component field traction[] */

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"eta\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"rho\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

}


/* ================================================== */
/* VTK binary (appended data) viewer for QPntSurfCoefStokes */
/* ==================================================== */
void QPntSurfCoefStokesVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const QPntSurfCoefStokes points[]) 
{
  int p,length;
  size_t atomic_size;

  /* Warning: swarm_class_generator.py is ignoring multi-component field normal[] */

  /* Warning: swarm_class_generator.py is ignoring multi-component field tangent1[] */

  /* Warning: swarm_class_generator.py is ignoring multi-component field tangent2[] */

  /* Warning: swarm_class_generator.py is ignoring multi-component field traction[] */

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].eta,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].rho,atomic_size,1,vtk_fp);
  }

}


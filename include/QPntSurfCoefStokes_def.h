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
 **    filename:   QPntSurfCoefStokes_def.h
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

#ifndef __QPntSurfCoefStokes_DEF_H__
#define __QPntSurfCoefStokes_DEF_H__

typedef struct {
  double normal [ 3 ] ;
  double tangent1 [ 3 ] ;
  double tangent2 [ 3 ] ;
  double traction [ 3 ] ;
  double eta ;
  double rho ;
} QPntSurfCoefStokes ;


typedef enum {
  QPSCStk_surface_normal = 0,
  QPSCStk_surface_tangent1,
  QPSCStk_surface_tangent2,
  QPSCStk_surface_traction,
  QPSCStk_viscosity,
  QPSCStk_density
} QPntSurfCoefStokesTypeName ;


extern const char QPntSurfCoefStokes_classname[];

extern const int QPntSurfCoefStokes_nmembers;

extern const size_t QPntSurfCoefStokes_member_sizes[];

extern const char *QPntSurfCoefStokes_member_names[];

/* prototypes */
void QPntSurfCoefStokesGetField_surface_normal(QPntSurfCoefStokes *point,double *data[]);
void QPntSurfCoefStokesGetField_surface_tangent1(QPntSurfCoefStokes *point,double *data[]);
void QPntSurfCoefStokesGetField_surface_tangent2(QPntSurfCoefStokes *point,double *data[]);
void QPntSurfCoefStokesGetField_surface_traction(QPntSurfCoefStokes *point,double *data[]);
void QPntSurfCoefStokesGetField_viscosity(QPntSurfCoefStokes *point,double *data);
void QPntSurfCoefStokesGetField_density(QPntSurfCoefStokes *point,double *data);
void QPntSurfCoefStokesSetField_surface_normal(QPntSurfCoefStokes *point,double data[]);
void QPntSurfCoefStokesSetField_surface_tangent1(QPntSurfCoefStokes *point,double data[]);
void QPntSurfCoefStokesSetField_surface_tangent2(QPntSurfCoefStokes *point,double data[]);
void QPntSurfCoefStokesSetField_surface_traction(QPntSurfCoefStokes *point,double data[]);
void QPntSurfCoefStokesSetField_viscosity(QPntSurfCoefStokes *point,double data);
void QPntSurfCoefStokesSetField_density(QPntSurfCoefStokes *point,double data);
void QPntSurfCoefStokesView(QPntSurfCoefStokes *point);
void QPntSurfCoefStokesVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const QPntSurfCoefStokes points[]);
void QPntSurfCoefStokesPVTUWriteAllPPointDataFields(FILE *vtk_fp);
void QPntSurfCoefStokesVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const QPntSurfCoefStokes points[]);
void QPntSurfCoefStokesVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const QPntSurfCoefStokes points[]);

#endif

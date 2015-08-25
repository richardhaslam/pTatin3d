
/*
  Auto generated by version 0.0 of swarm_class_generator.py
  on geop-265.ethz.ch, at 2015-08-25 11:20:15.233652 by dmay
*/

#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "MPntPEVSS_def.h"


const char MPntPEVSS_classname[] = "MPntPEVSS";

const int MPntPEVSS_nmembers = 3;

const size_t MPntPEVSS_member_sizes[] = {
  6 * sizeof(double),
  1 * sizeof(double),
  1 * sizeof(double)
};

const char *MPntPEVSS_member_names[] = {
  "deviatoric_stress",
  "shear_modulus",
  "solid_viscosity"
};

MPI_Datatype MPI_MPNTPEVSS;


/* ===================================== */
/* Getters for MPntPEVSS */
/* ===================================== */
void MPntPEVSSGetField_deviatoric_stress(MPntPEVSS *point,double *data[]) 
{
  *data = point->tau;
}

void MPntPEVSSGetField_shear_modulus(MPntPEVSS *point,double *data) 
{
  *data = point->mu;
}

void MPntPEVSSGetField_solid_viscosity(MPntPEVSS *point,double *data) 
{
  *data = point->eta_e;
}


/* ===================================== */
/* Setters for MPntPEVSS */
/* ===================================== */
void MPntPEVSSSetField_deviatoric_stress(MPntPEVSS *point,double data[]) 
{
  memcpy( &point->tau[0], data, sizeof(double)*6 );
}

void MPntPEVSSSetField_shear_modulus(MPntPEVSS *point,double data) 
{
  point->mu = data;
}

void MPntPEVSSSetField_solid_viscosity(MPntPEVSS *point,double data) 
{
  point->eta_e = data;
}


/* ===================================== */
/* C-viewer for MPntPEVSS */
/* ===================================== */
void MPntPEVSSView(MPntPEVSS *point)
{
  {
    double *data;
    MPntPEVSSGetField_deviatoric_stress(point,&data);
    printf("field: deviatoric_stress[0] = %1.6e; [size %zu; type double; variable_name tau]\n",data[0], MPntPEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[1] = %1.6e; [size %zu; type double; variable_name tau]\n",data[1], MPntPEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[2] = %1.6e; [size %zu; type double; variable_name tau]\n",data[2], MPntPEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[3] = %1.6e; [size %zu; type double; variable_name tau]\n",data[3], MPntPEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[4] = %1.6e; [size %zu; type double; variable_name tau]\n",data[4], MPntPEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[5] = %1.6e; [size %zu; type double; variable_name tau]\n",data[5], MPntPEVSS_member_sizes[0] );
  }
  {
    double data;
    MPntPEVSSGetField_shear_modulus(point,&data);
    printf("field: shear_modulus = %1.6e; [size %zu; type double; variable_name mu]\n",data, MPntPEVSS_member_sizes[1] );
  }
  {
    double data;
    MPntPEVSSGetField_solid_viscosity(point,&data);
    printf("field: solid_viscosity = %1.6e; [size %zu; type double; variable_name eta_e]\n",data, MPntPEVSS_member_sizes[2] );
  }
}


/* ===================================== */
/* VTK viewer for MPntPEVSS */
/* ===================================== */
void MPntPEVSSVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const MPntPEVSS points[]) 
{
  int p;
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"mu\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].mu);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"eta_e\" format=\"ascii\">\n");
  for(p=0;p<N;p++) {
    fprintf( vtk_fp,"\t\t\t\t\t%lf\n",(double)points[p].eta_e);
  }
  fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
}


/* ===================================== */
/* PVTU viewer for MPntPEVSS */
/* ===================================== */
void MPntPEVSSPVTUWriteAllPPointDataFields(FILE *vtk_fp) 
{
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"mu\" NumberOfComponents=\"1\"/>\n");
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"eta_e\" NumberOfComponents=\"1\"/>\n");
}


/* ===================================== */
/* VTK binary (appended header) viewer for MPntPEVSS */
/* ===================================== */
void MPntPEVSSVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const MPntPEVSS points[]) 
{
  /* Warning: swarm_class_generator.py is ignoring multi-component field tau[] */

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"mu\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"eta_e\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);

}


/* ================================================== */
/* VTK binary (appended data) viewer for MPntPEVSS */
/* ==================================================== */
void MPntPEVSSVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const MPntPEVSS points[]) 
{
  int p,length;
  size_t atomic_size;

  /* Warning: swarm_class_generator.py is ignoring multi-component field tau[] */

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].mu,atomic_size,1,vtk_fp);
  }

  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].eta_e,atomic_size,1,vtk_fp);
  }

}


/* ===================================== */
/* MPI data type for MPntPEVSS */
/* ===================================== */
int MPntPEVSSCreateMPIDataType(MPI_Datatype *ptype)
{
  MPI_Datatype newtype;
  MPI_Datatype types[] = { MPI_DOUBLE , MPI_DOUBLE , MPI_DOUBLE };
  int blocklens[] = { 6 , 1 , 1 };
  MPI_Aint loc[4];
  MPI_Aint disp[3];
  MPntPEVSS dummy;
  int i,ierr;

  ierr = MPI_Get_address(&dummy,&loc[0]);
  ierr = MPI_Get_address(&dummy.tau,&loc[1]);
  ierr = MPI_Get_address(&dummy.mu,&loc[2]);
  ierr = MPI_Get_address(&dummy.eta_e,&loc[3]);

  for (i=0; i<3; i++) {
    disp[i] = loc[i+1] - loc[0];
  }

  ierr = MPI_Type_create_struct(3,blocklens,disp,types,&newtype);
  ierr = MPI_Type_commit(&newtype);
  *ptype = newtype;
  return 0;
}


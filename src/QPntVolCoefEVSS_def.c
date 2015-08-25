
/*
  Auto generated by version 0.0 of swarm_class_generator.py
  on geop-265.ethz.ch, at 2015-08-25 11:20:15.230518 by dmay
*/

#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "QPntVolCoefEVSS_def.h"


const char QPntVolCoefEVSS_classname[] = "QPntVolCoefEVSS";

const int QPntVolCoefEVSS_nmembers = 1;

const size_t QPntVolCoefEVSS_member_sizes[] = {
  6 * sizeof(double)
};

const char *QPntVolCoefEVSS_member_names[] = {
  "deviatoric_stress"
};

MPI_Datatype MPI_QPNTVOLCOEFEVSS;


/* ===================================== */
/* Getters for QPntVolCoefEVSS */
/* ===================================== */
void QPntVolCoefEVSSGetField_deviatoric_stress(QPntVolCoefEVSS *point,double *data[]) 
{
  *data = point->tau;
}


/* ===================================== */
/* Setters for QPntVolCoefEVSS */
/* ===================================== */
void QPntVolCoefEVSSSetField_deviatoric_stress(QPntVolCoefEVSS *point,double data[]) 
{
  memcpy( &point->tau[0], data, sizeof(double)*6 );
}


/* ===================================== */
/* C-viewer for QPntVolCoefEVSS */
/* ===================================== */
void QPntVolCoefEVSSView(QPntVolCoefEVSS *point)
{
  {
    double *data;
    QPntVolCoefEVSSGetField_deviatoric_stress(point,&data);
    printf("field: deviatoric_stress[0] = %1.6e; [size %zu; type double; variable_name tau]\n",data[0], QPntVolCoefEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[1] = %1.6e; [size %zu; type double; variable_name tau]\n",data[1], QPntVolCoefEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[2] = %1.6e; [size %zu; type double; variable_name tau]\n",data[2], QPntVolCoefEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[3] = %1.6e; [size %zu; type double; variable_name tau]\n",data[3], QPntVolCoefEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[4] = %1.6e; [size %zu; type double; variable_name tau]\n",data[4], QPntVolCoefEVSS_member_sizes[0] );
    printf("field: deviatoric_stress[5] = %1.6e; [size %zu; type double; variable_name tau]\n",data[5], QPntVolCoefEVSS_member_sizes[0] );
  }
}


/* ===================================== */
/* VTK viewer for QPntVolCoefEVSS */
/* ===================================== */
void QPntVolCoefEVSSVTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const QPntVolCoefEVSS points[]) 
{
  int p;
}


/* ===================================== */
/* PVTU viewer for QPntVolCoefEVSS */
/* ===================================== */
void QPntVolCoefEVSSPVTUWriteAllPPointDataFields(FILE *vtk_fp) 
{
}


/* ===================================== */
/* VTK binary (appended header) viewer for QPntVolCoefEVSS */
/* ===================================== */
void QPntVolCoefEVSSVTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const QPntVolCoefEVSS points[]) 
{
  /* Warning: swarm_class_generator.py is ignoring multi-component field tau[] */

}


/* ================================================== */
/* VTK binary (appended data) viewer for QPntVolCoefEVSS */
/* ==================================================== */
void QPntVolCoefEVSSVTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const QPntVolCoefEVSS points[]) 
{
  int p,length;
  size_t atomic_size;

  /* Warning: swarm_class_generator.py is ignoring multi-component field tau[] */

}


/* ===================================== */
/* MPI data type for QPntVolCoefEVSS */
/* ===================================== */
int QPntVolCoefEVSSCreateMPIDataType(MPI_Datatype *ptype)
{
  MPI_Datatype newtype;
  MPI_Datatype types[] = { MPI_DOUBLE };
  int blocklens[] = { 6 };
  MPI_Aint loc[2];
  MPI_Aint disp[1];
  QPntVolCoefEVSS dummy;
  int i,ierr;

  ierr = MPI_Get_address(&dummy,&loc[0]);
  ierr = MPI_Get_address(&dummy.tau,&loc[1]);

  for (i=0; i<1; i++) {
    disp[i] = loc[i+1] - loc[0];
  }

  ierr = MPI_Type_create_struct(1,blocklens,disp,types,&newtype);
  ierr = MPI_Type_commit(&newtype);
  *ptype = newtype;
  return 0;
}


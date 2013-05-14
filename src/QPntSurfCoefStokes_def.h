
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

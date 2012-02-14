
#ifndef __ptatin_quadrature_h__
#define __ptatin_quadrature_h__

#include "QPntVolCoefStokes_def.h"

PetscErrorCode QuadratureCreate(Quadrature *quadrature);
PetscErrorCode QuadratureDestroy(Quadrature *quadrature);
PetscErrorCode QuadratureView(Quadrature q);
PetscErrorCode VolumeQuadratureCreate_GaussLegendreStokes(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature);
PetscErrorCode VolumeQuadratureGetAllCellData_Stokes(Quadrature Q,QPntVolCoefStokes *coeffs[]);
PetscErrorCode VolumeQuadratureGetCellData_Stokes(Quadrature Q,QPntVolCoefStokes coeffs[],PetscInt cidx,QPntVolCoefStokes *cell[]);

#endif


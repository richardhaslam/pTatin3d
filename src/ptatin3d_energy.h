
#ifndef __ptatin3d_energy_h__
#define __ptatin3d_energy_h__

#include "QPntVolCoefEnergy_def.h"

PetscErrorCode AdvDiff_GetElementLocalIndices_Q1(PetscInt el_localIndices[],PetscInt elnid[]);

double AdvDiffResidualForceTerm_UpwindXiExact(double pecletNumber);
double AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption(double pecletNumber);
double AdvDiffResidualForceTerm_UpwindXiCriticalAssumption(double pecletNumber);

PetscErrorCode VolumeQuadratureCreate_GaussLegendreEnergy(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature);
PetscErrorCode VolumeQuadratureGetAllCellData_Energy(Quadrature Q,QPntVolCoefEnergy *coeffs[]);
PetscErrorCode VolumeQuadratureGetCellData_Energy(Quadrature Q,QPntVolCoefEnergy coeffs[],PetscInt cidx,QPntVolCoefEnergy *cell[]);

#endif

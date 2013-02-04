
#ifndef __phys_comp_energy_h__
#define __phys_comp_energy_h__

#include "QPntVolCoefEnergy_def.h"

PetscErrorCode PhysCompCreate_Energy(PhysCompEnergy *E);
PetscErrorCode PhysCompDestroy_Energy(PhysCompEnergy *E);

PetscErrorCode PhysCompNew_Energy(DM dav,PetscInt mx,PetscInt my, PetscInt mz,PetscInt mesh_generator_type,PhysCompEnergy *E);

PetscErrorCode PhysCompCreateMesh_Energy(PhysCompEnergy E,DM dav,PetscInt mx,PetscInt my, PetscInt mz,PetscInt mesh_generator_type);
PetscErrorCode PhysCompCreateBoundaryList_Energy(PhysCompEnergy E);
PetscErrorCode PhysCompCreateVolumeQuadrature_Energy(PhysCompEnergy E);

PetscErrorCode PhysCompLoad_Energy(void);
PetscErrorCode PhysCompSave_Energy(void);

PetscErrorCode VolumeQuadratureCreate_GaussLegendreEnergy(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature);
PetscErrorCode VolumeQuadratureGetAllCellData_Energy(Quadrature Q,QPntVolCoefEnergy *coeffs[]);
PetscErrorCode VolumeQuadratureGetCellData_Energy(Quadrature Q,QPntVolCoefEnergy coeffs[],PetscInt cidx,QPntVolCoefEnergy *cell[]);

#endif

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
 **    filename:   phys_comp_energy.h
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

PetscErrorCode PhysCompAddMaterialPointCoefficients_Energy(DataBucket db);
PetscErrorCode PhysCompCheckpointWrite_Energy(PhysCompEnergy e,PetscBool write_dmda,const char path[],const char prefix[]);
PetscErrorCode PhysCompLoad2_Energy(DM dav,const char jfilename[],PhysCompEnergy *E);

#endif

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
 **    filename:   ptatin3d_energy.h
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

#ifndef __ptatin3d_energy_h__
#define __ptatin3d_energy_h__

PetscErrorCode pTatinGetContext_Energy(pTatinCtx ctx,PhysCompEnergy *e);
PetscErrorCode pTatinContextValid_Energy(pTatinCtx ctx,PetscBool *exists);
PetscErrorCode pTatinPhysCompCreate_Energy(pTatinCtx user);
PetscErrorCode pTatinPhysCompActivate_Energy(pTatinCtx user,PetscBool load);

PetscErrorCode pTatinPhysCompAttachData_Energy(pTatinCtx user,Vec T,Mat A);
PetscErrorCode pTatinPhysCompGetData_Energy(pTatinCtx user,Vec *T,Mat *A);

PetscErrorCode pTatinPhysCompEnergy_UpdateALEVelocity(PhysCompStokes s,Vec X,PhysCompEnergy energy,PetscReal dt);
PetscErrorCode pTatinPhysCompEnergy_Update(PhysCompEnergy e,DM dav,Vec T);
PetscErrorCode pTatinPhysCompEnergy_Initialise(PhysCompEnergy e,Vec T);

PetscErrorCode MaterialPointQuadraturePointProjectionC0_Q2Energy(DM da,DataBucket materialpoint_db,MaterialPointField field,const int member,Quadrature Q);
PetscErrorCode pTatinPhysCompEnergy_MPProjectionQ1(pTatinCtx ctx);
PetscErrorCode pTatinPhysCompEnergy_ComputeTimestep(PhysCompEnergy energy,Vec X,PetscReal *timestep);

#endif

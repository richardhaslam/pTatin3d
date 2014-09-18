/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      dmda_redundant.h
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/
#ifndef __PTATIN3d_DMDA_REDUNDANT_H__
#define __PTATIN3d_DMDA_REDUNDANT_H__

#include "sub_comm.h"
#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

PetscErrorCode  x_DMDACreate3d(MPI_Comm comm,DMBoundaryType wrap[],DMDAStencilType stencil_type,
															 PetscInt M, PetscInt N,PetscInt P,PetscInt m,PetscInt n,PetscInt p,
															 PetscInt dof,PetscInt s,const PetscInt lx[],const PetscInt ly[],const PetscInt lz[],DM *da);

PetscErrorCode DMDACreate3dRedundant(DM da,PetscInt si, PetscInt ei, PetscInt sj, PetscInt ej, PetscInt sk, PetscInt ek, PetscInt n_dofs, DM *_seq_DM );
PetscErrorCode DMDACreate3dSemiRedundant(DM da,PetscInt nred,PetscMPISubComm *sub,DM *sda);

#endif


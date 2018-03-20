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
 **    filename:   ptatin3d_energy_impl.h
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

#ifndef __private_ptatin3d_energy_impl_h__
#define __private_ptatin3d_energy_impl_h__

#include <petscdm.h>
#include "dmda_bcs.h"

struct _p_PhysCompEnergy {
	PetscReal               time,dt;
	PetscInt                energy_mesh_type; /* 0-std dmda, 1-overlap, 2-nested */
	PetscInt                mx,my,mz; /* global mesh size */
	DM                      daT;
	BCList                  T_bclist;
	Vec                     u_minus_V;
	Quadrature              volQ;
	//	SurfaceQuadratureEnergy surfQ[QUAD_EDGES]; /* four edges */
	/* SUPG DATA */
	Vec                     Told; /* old temperature solution vector */
	Vec                     Xold; /* old coordinate vector */
};

#endif

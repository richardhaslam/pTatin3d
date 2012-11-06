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
 **    Filename:      viscous_sinker_ctx.h
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


#ifndef __ptatinmodel_viscous_sinker_ctx_h__
#define __ptatinmodel_viscous_sinker_ctx_h__

typedef enum { VSBC_FreeSlip=0, VSBC_NoSlip, VSBC_FreeSlipFreeSurface, VSBC_NoSlipFreeSurface, VSBC_Test } ViscousSinkerBC;

typedef struct {
	PetscInt  nmaterials;
	PetscReal eta0,eta1;
	PetscReal rho0,rho1;
	PetscBool is_sphere;
	PetscReal Lx,Ly,Lz;
	PetscReal origin[3];
	PetscReal length[3];
	ViscousSinkerBC  boundary_conditon_type; /* [ 0 free slip | 1 no slip | 2 free surface + free slip | 3 free surface + no slip ] */
} ModelViscousSinkerCtx;


#endif

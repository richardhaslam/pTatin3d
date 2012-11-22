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
 **    Filename:      ptatin_std_dirichlet_boundary_conditions.h
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

#ifndef _ptatin_std_dirichlet_boundary_conditions_h__
#define _ptatin_std_dirichlet_boundary_conditions_h__


typedef enum { FRONT_FACE=1, BACK_FACE, EAST_FACE, NORTH_FACE, SOUTH_FACE, WEST_FACE } BoundaryFaceType;

PetscErrorCode DirichletBC_FreeSlip(BCList list,DM dav,BoundaryFaceType face);
PetscErrorCode DirichletBC_ApplyNormalVelocity(BCList list,DM dav,BoundaryFaceType face,PetscReal v_normal);
PetscErrorCode DirichletBC_ApplyDirectStrainRate(BCList bclist,DM dav,PetscReal Evalue,PetscInt direction);
PetscErrorCode DirichletBC_ApplyStrainRateExx(BCList list,DM dav,PetscReal exx_bc);
PetscErrorCode DirichletBC_ApplyStrainRateExz(BCList bclist,DM dav,PetscReal exz_bc);
PetscErrorCode DirichletBC_ApplyStrainRateExz_b(BCList bclist,DM dav,PetscReal exz_bc);
PetscErrorCode DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ(BCList list,DM dav,PetscReal vx_bc);
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionX(BCList list,DM dav,PetscReal vx_bc);
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(BCList list,DM dav,PetscReal beta,PetscReal vx_bc);


#endif


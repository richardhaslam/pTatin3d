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
 **    Filename:      quadrature.h
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

#ifndef __ptatin_quadrature_h__
#define __ptatin_quadrature_h__

#include "element_type_Q2.h"


PetscErrorCode QuadratureCreate(Quadrature *quadrature);
PetscErrorCode QuadratureDestroy(Quadrature *quadrature);
PetscErrorCode QuadratureView(Quadrature q);

void QuadratureCreateGauss_2pnt_3D(PetscInt *ngp,PetscReal **_q_coor,PetscReal **_q_weight);
void QuadratureCreateGauss_3pnt_3D(PetscInt *ngp,PetscReal **_q_coor,PetscReal **_q_weight);

PetscErrorCode SurfaceQuadratureCreate(SurfaceQuadrature *quadrature);
PetscErrorCode SurfaceQuadratureDestroy(SurfaceQuadrature *quadrature);
PetscErrorCode _SurfaceQuadratureCreate(SurfaceQuadrature quadrature,HexElementFace index,PetscInt nfaces);
PetscErrorCode _SurfaceQuadratureCellIndexSetUp(SurfaceQuadrature Q,HexElementFace index,PetscInt nface_edge,DM da);

#endif


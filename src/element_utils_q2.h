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
 **    Filename:      element_utils_q2.h
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

#ifndef __ptatin_element_utils_q2_h__
#define __ptatin_element_utils_q2_h__

#include "petsc.h"
#include "ptatin3d_defs.h"

#define NQP 27
#define NPE Q2_NODES_PER_EL_3D

void P3D_ConstructNi_Q2_3D(PetscReal _xi[],PetscReal Ni[] );
void P3D_ConstructGNi_Q2_3D(PetscReal _xi[],PetscReal GNi[3][Q2_NODES_PER_EL_3D] );

void P3D_ConstructNi_P0_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);
void P3D_ConstructNi_P1L_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);
void P3D_ConstructNi_P1G_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);
void P3D_ConstructNi_P1GRel_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);

void P3D_prepare_elementQ2_3x3(PetscReal WEIGHT[NQP],PetscReal XI[NQP][3],PetscReal NI[NQP][NPE],PetscReal GNI[NQP][3][NPE]);
void P3D_prepare_elementQ2_2x2(PetscReal WEIGHT[NQP],PetscReal XI[NQP][3],PetscReal NI[NQP][NPE],PetscReal GNI[NQP][3][NPE]);
void P3D_prepare_elementQ2(PetscInt nqp,PetscReal WEIGHT[NQP],PetscReal XI[NQP][3],PetscReal NI[NQP][NPE],PetscReal GNI[NQP][3][NPE]);

void P3D_evaluate_geometry_elementQ2(PetscInt nqp,PetscReal el_coords[NPE*3],PetscReal GNI[][3][NPE],
																		 PetscReal detJ[],
																		 PetscReal dNudx[][NPE],
																		 PetscReal dNudy[][NPE],
																		 PetscReal dNudz[][NPE] );
void P3D_evaluate_geometry_elementQ2_1gp(
																				 PetscReal GNI_centre[3][NPE],
																				 PetscInt nqp,PetscReal el_coords[NPE*3],PetscReal GNI[][3][NPE],
																				 PetscReal detJ[],
																				 PetscReal dNudx[][NPE],
																				 PetscReal dNudy[][NPE],
																				 PetscReal dNudz[][NPE] );


#endif
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
 **    filename:   element_utils_q2.h
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

#ifndef __ptatin_element_utils_q2_h__
#define __ptatin_element_utils_q2_h__

#include "petsc.h"
#include "ptatin3d_defs.h"

#define NQP 27
#define NPE Q2_NODES_PER_EL_3D

void P3D_ConstructNi_Q2_2D(PetscReal _xi[],PetscReal Ni[]);
void P3D_ConstructGNi_Q2_2D(PetscReal _xi[],PetscReal GNi[2][Q2_NODES_PER_EL_2D]);

void P3D_ConstructNi_Q2_3D(PetscReal _xi[],PetscReal Ni[] );
void P3D_ConstructGNi_Q2_3D(PetscReal _xi[],PetscReal GNi[3][Q2_NODES_PER_EL_3D] );

void P3D_ConstructNi_P0_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);
void P3D_ConstructNi_P1L_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);
void P3D_ConstructNi_P1G_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);
void P3D_ConstructNi_P1GRel_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);
void P3D_ConstructNi_P1GRel_VertexBased_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[]);

void P3D_prepare_elementQ2_3x3(PetscReal WEIGHT[],PetscReal XI[][3],PetscReal NI[][NPE],PetscReal GNI[][3][NPE]);
void P3D_prepare_elementQ2_2x2(PetscReal WEIGHT[],PetscReal XI[][3],PetscReal NI[][NPE],PetscReal GNI[][3][NPE]);
void P3D_prepare_elementQ2(PetscInt nqp,PetscReal WEIGHT[],PetscReal XI[][3],PetscReal NI[][NPE],PetscReal GNI[][3][NPE]);

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
void P3D_evaluate_geometry_elementQ2_1gp_diagonal(
                                                  PetscReal GNI_centre[3][NPE],
                                                  PetscInt nqp,PetscReal el_coords[NPE*3],PetscReal GNI[][3][NPE],
                                                  PetscReal detJ[],
                                                  PetscReal dNudx[][NPE],
                                                  PetscReal dNudy[][NPE],
                                                  PetscReal dNudz[][NPE] );

void P3D_evaluate_global_derivatives_Q2(PetscReal el_coords[NPE*3],PetscReal GNI[3][NPE],
                                        PetscReal dNudx[NPE],
                                        PetscReal dNudy[NPE],
                                        PetscReal dNudz[NPE] );


#endif


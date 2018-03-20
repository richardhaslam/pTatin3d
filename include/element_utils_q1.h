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
 **    filename:   element_utils_q1.h
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

#ifndef __ptatin_element_utils_q1_h__
#define __ptatin_element_utils_q1_h__

#include "petsc.h"
#include "ptatin3d_defs.h"

#define NODES_PER_EL_Q1_3D 8

void P3D_ConstructNi_Q1_3D(PetscReal _xi[],PetscReal Ni[]);
void P3D_ConstructGNi_Q1_3D(PetscReal _xi[],PetscReal GNi[3][8]);
void P3D_evaluate_geometry_elementQ1(PetscInt nqp,PetscReal el_coords[8*3],PetscReal GNI[][3][8],
                                     PetscReal detJ[],
                                     PetscReal dNudx[][8],
                                     PetscReal dNudy[][8],
                                     PetscReal dNudz[][8] );
void P3D_evaluate_geometry_elementQ1_appliedQ2(PetscInt nqp,
                                               PetscReal detJ[],
                                               PetscReal GNIQ1[][3][8],
                                               PetscReal el_coords[27*3],
                                               PetscReal GNIQ2[][3][27],
                                               PetscReal dNudxQ2[][27],
                                               PetscReal dNudyQ2[][27],
                                               PetscReal dNudzQ2[][27] );
void P3D_evaluate_geometry_affine_appliedQ2(PetscInt nqp,
                                            PetscReal detJ[],
                                            PetscReal GNIQ1[][3][8],
                                            PetscReal el_coords[27*3],
                                            PetscReal GNIQ2[][3][27],
                                            PetscReal dNudxQ2[][27],
                                            PetscReal dNudyQ2[][27],
                                            PetscReal dNudzQ2[][27] );
void P3D_evaluate_geometry_affine2_appliedQ2(PetscInt nqp,
                                            PetscReal detJ[],
                                            PetscReal GNIQ1[][3][8],
                                            PetscReal el_coords[27*3],
                                            PetscReal GNIQ2[][3][27],
                                            PetscReal dNudxQ2[][27],
                                            PetscReal dNudyQ2[][27],
                                            PetscReal dNudzQ2[][27] );

void P3D_ConstructNi_Q1_2D(PetscReal _xi[],PetscReal Ni[]);

void P3D_ConstructGNi_Q1_2D(PetscReal _xi[],PetscReal GNix[],PetscReal GNiy[]);

void P3D_evaluate_geometry_elementQ1_2D(PetscReal el_coords[],PetscReal GNIx[],PetscReal GNIy[],
                                        PetscReal *detJ,PetscReal dNudx[],PetscReal dNudy[]);


#endif


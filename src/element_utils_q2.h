
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


#endif
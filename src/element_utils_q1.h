
#ifndef __ptatin_element_utils_q1_h__
#define __ptatin_element_utils_q1_h__

#include "petsc.h"
#include "ptatin3d_defs.h"

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

#endif
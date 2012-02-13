
#ifndef __ptatin3d_defs_h__
#define __ptatin3d_defs_h__


#define MAX_QUAD_PNTS_1D   10
#define MAX_QUAD_PNTS      MAX_QUAD_PNTS_1D * MAX_QUAD_PNTS_1D * MAX_QUAD_PNTS_1D
#define Q2_NODES_PER_EL_3D 27
#define Q2_NODES_PER_EL_2D 9
#define Q2_NODES_PER_EL_1D 3

#define NSD    3
#define U_DOFS NSD
#define P_DOFS 1
#define NSTRESS 6

#define U_BASIS_FUNCTIONS Q2_NODES_PER_EL_3D


#define PBASIS_P1_GLOBAL_RELATIVE
//#define PBASIS_P0

/* choose a default pressure space if nothing is already defined */
#ifndef PBASIS_P0
#ifndef PBASIS_P1_GLOBAL
#ifndef PBASIS_P1_LOCAL
#ifndef PBASIS_P1_GLOBAL_RELATIVE
#define PBASIS_P0
#endif
#endif
#endif
#endif

#ifdef PBASIS_P0
#define P_BASIS_FUNCTIONS 1
#define ConstructNi_pressure(_xi,coords,Ni) ConstructNi_P0_3D(_xi,coords,Ni)
#endif
#ifdef PBASIS_P1_GLOBAL
#define P_BASIS_FUNCTIONS 4
#define ConstructNi_pressure(_xi,coords,Ni) ConstructNi_P1G_3D(_xi,coords,Ni)
#endif
#ifdef PBASIS_P1_LOCAL
#define P_BASIS_FUNCTIONS 4
#define ConstructNi_pressure(_xi,coords,Ni) ConstructNi_P1L_3D(_xi,coords,Ni)
#endif
#ifdef PBASIS_P1_GLOBAL_RELATIVE
#define P_BASIS_FUNCTIONS 4
#define ConstructNi_pressure(_xi,coords,Ni) ConstructNi_P1GRel_3D(_xi,coords,Ni)
#endif

#endif

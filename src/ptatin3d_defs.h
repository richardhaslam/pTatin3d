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
 **    Filename:      ptatin3d_defs.h
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

#ifndef __ptatin3d_defs_h__
#define __ptatin3d_defs_h__


#define MAX_QUAD_PNTS_1D   10
#define MAX_QUAD_PNTS      MAX_QUAD_PNTS_1D * MAX_QUAD_PNTS_1D * MAX_QUAD_PNTS_1D
#define Q2_NODES_PER_EL_3D 27
#define Q2_NODES_PER_EL_2D 9
#define Q2_NODES_PER_EL_1D 3

#define Q1_NODES_PER_EL_3D 8
#define Q1_NODES_PER_EL_2D 4
#define Q1_NODES_PER_EL_1D 2

#define SURF_QUAD_FACES 8

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
#define ConstructNi_pressure(_xi,coords,Ni) P3D_ConstructNi_P0_3D(_xi,coords,Ni)
#endif
#ifdef PBASIS_P1_GLOBAL
#define P_BASIS_FUNCTIONS 4
#define ConstructNi_pressure(_xi,coords,Ni) P3D_ConstructNi_P1G_3D(_xi,coords,Ni)
#endif
#ifdef PBASIS_P1_LOCAL
#define P_BASIS_FUNCTIONS 4
#define ConstructNi_pressure(_xi,coords,Ni) P3D_ConstructNi_P1L_3D(_xi,coords,Ni)
#endif
#ifdef PBASIS_P1_GLOBAL_RELATIVE
#define P_BASIS_FUNCTIONS 4
#define ConstructNi_pressure(_xi,coords,Ni) P3D_ConstructNi_P1GRel_3D(_xi,coords,Ni)
#endif

#endif


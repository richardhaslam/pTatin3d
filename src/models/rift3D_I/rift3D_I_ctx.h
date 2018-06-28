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
 **    filename:   rift3D_T_ctx.h
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

#ifndef __ptatinmodel_rift3d_I_ctx_h__
#define __ptatinmodel_rift3d_I_ctx_h__

typedef struct {
  PetscInt  nmaterials;
  PetscBool runmises;
  PetscReal Lx,Ly,Lz,Ox,Oy,Oz,vx,vy,vz,vx_I,vz_I,rho0;
  PetscReal time_pause_start,time_pause,time_invert_start,time_invert;
  PetscReal Tbottom,Ttop,thermal_age0,Tlitho,h_prod,h_prod_rcp,y_prod,qm,ylab,k_ref;
  PetscReal thermal_age_anom,wx_anom,wz_anom,cx_anom,cz_anom;
  PetscBool dimensional;
  PetscReal density_bar;
  PetscReal length_bar;
  PetscReal viscosity_bar;
  PetscReal velocity_bar;
  PetscReal h_prod_bar;
  PetscReal k_bar;
  PetscReal time_bar;
  PetscReal pressure_bar;
  PetscBool use_semi_eulerian_mesh;
  PetscBool output_markers;
  GeometryObject G[100];
  PetscInt  ngo;
  PSwarm      pswarm;
} ModelRift3D_ICtx;

typedef struct _p_BC_SplitFace *BC_SplitFace;
struct _p_BC_SplitFace {
  PetscInt dim;
  PetscScalar x0,x1;
  PetscScalar v0,v1;
};
#endif
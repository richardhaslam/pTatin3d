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
 **    filename:   model_riftrh_ctx.h
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

#ifndef __ptatin3d_model_riftrh_ctx_h__
#define __ptatin3d_model_riftrh_ctx_h__

/* define user model */
typedef struct {
  PetscReal Lx,Ly,Lz,hc,hm,ha,dxs,dys,dzs,hvbx1,hvbx2;
  PetscReal vx_up,vx_down,rhoc,rhom,rhoa,etac,etam,etaa;

  PetscInt  nmaterials;
  PetscBool runmises;
  PetscBool dimensional;
  PetscReal density_bar;
  PetscReal length_bar;
  PetscReal viscosity_bar;
  PetscReal velocity_bar;
  PetscReal time_bar;
  PetscReal pressure_bar;
  PetscReal Ttop;
  PetscReal Tbottom;
  PetscBool use_semi_eulerian_mesh;
  PetscBool output_markers;
  PetscInt  param1,param2;
  DMDA_thermalfield_init_params thermalparams;
} ModelRiftrhCtx;

#endif

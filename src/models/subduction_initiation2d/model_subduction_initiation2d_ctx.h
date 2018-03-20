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
 **    filename:   model_subduction_initiation2d_ctx.h
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

#ifndef __ptatin3d_model_subduction_initiation2d_ctx_h__
#define __ptatin3d_model_subduction_initiation2d_ctx_h__

/* define user model */
typedef struct {
        PetscReal Lx,Ly,Lz,Ox,Oy,Oz;
  PetscReal velocity;
  PetscReal eta[1],rho[1];
  PetscReal C0[1],mu[1],C0_inf[1],mu_inf[1];
        PetscReal diffusivity[1],alpha[1],H[1];
  PetscReal Ttop,Tbot;
  PetscReal Thermal_age;
  PetscReal density_bar,length_bar,viscosity_bar,velocity_bar,time_bar,pressure_bar;
} ModelSubduction_Initiation2dCtx;

#endif

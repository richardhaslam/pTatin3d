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
 **    Filename:      rift3D_ctx.h
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
 **    $Id:$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#ifndef __ptatinmodel_rift3d_ctx_h__
#define __ptatinmodel_rift3d_ctx_h__


typedef struct {
	PetscInt  nmaterials;
    PetscBool runmises;
	PetscReal Lx,Ly,Lz,Ox,Oy,Oz,vx,vy,vz,rho0;
	PetscReal eta[20];
	PetscReal rho[20];
	PetscReal phi[20];
	PetscReal coh[20];
	PetscReal arrh_n[20];
	PetscReal arrh_E[20];
	PetscReal alpha[20];
	PetscBool dimensional;
	PetscReal density_bar;
	PetscReal length_bar;
	PetscReal viscosity_bar;
	PetscReal velocity_bar;
	PetscReal time_bar;
	PetscReal pressure_bar;
} ModelRift3DCtx;


#endif

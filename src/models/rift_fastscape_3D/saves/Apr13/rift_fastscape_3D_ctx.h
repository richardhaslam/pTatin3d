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
 **    Filename:      model_rift_oblique3d_ctx.h
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
 **    $Id: model_rift_oblique3d_ctx.h 3636 2012-11-06 19:14:24Z dmay $
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#ifndef __ptatin3d_model_rift_fastscape_3D_ctx_h__
#define __ptatin3d_model_rift_fastscape_3D_ctx_h__

/* define user model */
typedef struct {
	PetscReal Lx,Ly,Lz,hc,hm,ha,hvbx1,hvbx2;
	PetscReal notch_width,notch_height,notch_base,beta,fault_width;
	PetscReal vx_up,vx_down,vybottom,rhoc,rhom,rhowm,rhoa,etac,etam,etawm, etaa,eps1,eps2,phi1,phi2,coe1,coe2;
    
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
	PetscInt  param1,param2,notch_type;
    DMDA_thermalfield_init_params thermalparams;
	PSwarm			pswarm;
} Modelrift_fastscape_3DCtx;

#endif

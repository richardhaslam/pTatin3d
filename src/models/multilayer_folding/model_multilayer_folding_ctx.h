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
 **    filename:   model_multilayer_folding_ctx.h
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
#ifndef __ptatin3d_model_multilayer_folding_ctx_h__
#define __ptatin3d_model_multilayer_folding_ctx_h__

/* define user model */
typedef struct {
	PetscInt  max_layers,seed_layer_1;
	PetscInt  n_interfaces;
	PetscReal interface_heights[100];
    PetscInt  layer_res_j[99];
	PetscReal eta[100];
	PetscReal rho[100];
	PetscReal cohesion[100];
	PetscInt  bc_type, perturbation_type;
	PetscReal ezz;
	PetscReal exx;
	PetscReal vx_compression;
	PetscReal vz_compression;
	PetscReal Lx, Lz, Ly, L_char;
    PetscReal amp;
    PetscReal kx,kz;
    PetscReal A0;
	PetscBool visco_plastic,quasi2d;
} ModelMultilayerFoldingCtx;

#endif



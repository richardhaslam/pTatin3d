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
 **    filename:   model_gene3d_ctx.h
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


#ifndef __ptatinmodel_GENE3D_ctx_h__
#define __ptatinmodel_GENE3D_ctx_h__

typedef enum { GENEBC_FreeSlip=0, GENEBC_NoSlip, GENEBC_FreeSlipFreeSurface, GENEBC_NoSlipFreeSurface } GENE3DBC;
typedef enum { GENE_LayeredCake=0, GENE_ExtrudeFromMap, GENE_ReadFromCAD} GENE3DINIGEOM;
enum {LAYER_MAX = 100};

typedef struct {
	PetscInt  nmaterials;
	PetscReal Lx,Ly,Lz;
	PetscReal Ox,Oy,Oz;
	GENE3DBC  boundary_conditon_type; /* [ 0 free slip | 1 no slip | 2 free surface + free slip | 3 free surface + no slip ] */
	GENE3DINIGEOM  initial_geom;
} ModelGENE3DCtx;

PetscErrorCode ModelSetMarkerIndexLayeredCake_GENE3D(pTatinCtx c,void *ctx);
PetscErrorCode ModelSetMarkerIndexFromMap_GENE3D(pTatinCtx c,void *ctx);
PetscErrorCode ModelSetInitialStokesVariableOnMarker_GENE3D(pTatinCtx c,void *ctx);
		
#endif

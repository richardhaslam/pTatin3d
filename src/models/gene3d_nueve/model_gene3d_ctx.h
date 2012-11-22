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
 **    Filename:      model_gene3d_ctx.h
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


#ifndef __ptatinmodel_GENE3DNUEVE_ctx_h__
#define __ptatinmodel_GENE3DNUEVE_ctx_h__

typedef enum { GENEBC_FreeSlip=0, GENEBC_NoSlip, GENEBC_FreeSlipFreeSurface, GENEBC_NoSlipFreeSurface, GENEBC_ShearFreeSlipFreeSurface } GENE3DBC;
typedef enum { GENE_LayeredCake=0, GENE_ExtrudeFromMap, GENE_ReadFromCAD, GENE_ExtrudeByPartsFromMap} GENE3DINIGEOM;

typedef struct {
	PetscInt  nmaterials;
	PetscReal Lx,Ly,Lz;
	PetscReal Ox,Oy,Oz;
	GENE3DBC  boundary_conditon_type; /* [ 0 free slip | 1 no slip | 2 free surface + free slip | 3 free surface + no slip ] */
	GENE3DINIGEOM  initial_geom;
} ModelGene3DNueveCtx;

PetscErrorCode ModelSetMarkerIndexLayeredCake_Gene3DNueve(pTatinCtx c,void *ctx);
PetscErrorCode ModelSetMarkerIndexFromMap_Gene3DNueve(pTatinCtx c,void *ctx);
PetscErrorCode ModelSetInitialStokesVariableOnMarker_Gene3DNueve(pTatinCtx c,void *ctx);
		
#endif

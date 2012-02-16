

#ifndef __ptatinmodel_GENE3D_ctx_h__
#define __ptatinmodel_GENE3D_ctx_h__

typedef enum { GENEBC_FreeSlip=0, GENEBC_NoSlip, GENEBC_FreeSlipFreeSurface, GENEBC_NoSlipFreeSurface } GENE3DBC;
typedef enum { GENE_LayeredCake=0, GENE_ExtrudeFromMap, GENE_ReadFromCAD} GENE3DINIGEOM;
const PetscInt LAYER_MAX = 100; 

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

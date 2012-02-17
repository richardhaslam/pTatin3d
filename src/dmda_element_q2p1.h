
#ifndef __ptatin_dmda_element_q2p1_h__
#define __ptatin_dmda_element_q2p1_h__

#include "petsc.h"
#include "petscdm.h"
#include "private/daimpl.h" 

PetscErrorCode DMDAGetCornersElementQ2(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz);
PetscErrorCode DMDAGetSizeElementQ2(DM da,PetscInt *MX,PetscInt *MY,PetscInt *MZ);
PetscErrorCode DMDAGetLocalSizeElementQ2(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz);
PetscErrorCode DMDAGetCornersElementQ2(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz);
PetscErrorCode DMDAGetOwnershipRangesElementQ2(DM da,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt **si,PetscInt **sj,PetscInt **sk,PetscInt **_mx,PetscInt **_my,PetscInt **_mz);

PetscErrorCode DMDAGetElements_DA_Q2_3D(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx);
PetscErrorCode DMDAGetElements_DA_P0MD_3D(DM dm,PetscInt *nel,PetscInt *npe,const PetscInt **eidx);

PetscErrorCode DMDAGetElements_DA_Q2(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]);
PetscErrorCode DMDAGetElements_DA_P1(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]);

PetscErrorCode  DMDASetElementType_Q2(DM da);
PetscErrorCode  DMDASetElementType_P1(DM da);

PetscErrorCode DMDAGetElements_pTatinQ2P1(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]);

PetscErrorCode DMDAGetElementCoordinatesQ2_3D(PetscScalar elcoords[],PetscInt elnid[],PetscScalar LA_gcoords[]);
PetscErrorCode DMDAGetScalarElementFieldQ2_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDAGetVectorElementFieldQ2_3D(PetscScalar elfield[],PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDAGetScalarElementField(PetscScalar elfield[],PetscInt npe,PetscInt elnid[],PetscScalar LA_gfield[]);
PetscErrorCode DMDASetValuesLocalStencil_AddValues_DOF(PetscScalar *fields_F,PetscInt ndof,PetscInt eqn[],PetscScalar Fe[]);
PetscErrorCode Q2GetElementLocalIndicesDOF(PetscInt el_localIndices[],PetscInt ndof,PetscInt elnid[]);
PetscErrorCode DMDAGetScalarElementField(PetscScalar elfield[],PetscInt npe,PetscInt elnid[],PetscScalar LA_gfield[]);

#endif

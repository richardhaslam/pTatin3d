

#ifndef __ptatin3d_h__
#define __ptatin3d_h__

#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"



typedef struct _p_pTatinCtx *pTatinCtx;
typedef struct _p_pTatinModel *pTatinModel;
typedef struct _p_PhysCompStokes *PhysCompStokes;
typedef struct _p_Quadrature *Quadrature;

typedef enum { LINE_QUAD=0,SURFACE_QUAD,VOLUME_QUAD } QuadratureType;

typedef struct _p_RheologyConstants RheologyConstants;


/*
#include "dmda_bcs.h"
#include "dmda_checkpoint.h"
#include "dmda_compare.h"
#include "dmda_duplicate.h"
#include "dmda_element_q2p1.h"
#include "dmda_project_coords.h"
#include "dmda_redundant.h"
#include "dmda_remesh.h"
#include "dmda_update_coords.h"
#include "dmda_view_petscvtk.h"

#include "swarm_fields.h"
#include "data_exchanger.h"

#include "ptatin3d_defs.h"
#include "ptatin_utils.h"
#include "ptatin3d_stokes.h"
#include "ptatin_models.h"
#include "rheology.h"
*/

#include "swarm_fields.h"
#include "data_exchanger.h"
//#include "rheology.h"
#include "material_point_load.h"
#include "material_point_utils.h"


PetscErrorCode pTatin3d_PhysCompStokesCreate(pTatinCtx user);
PetscErrorCode pTatin3d_ModelOutput_VelocityPressure_Stokes(pTatinCtx ctx,Vec X,const char prefix[]);

PetscErrorCode pTatin3dCreateMaterialPoints(pTatinCtx ctx,DM dav);
PetscErrorCode MaterialPointCoordinateSetUp(pTatinCtx ctx,DM da);
PetscErrorCode pTatin3d_ModelOutput_MPntStd(pTatinCtx ctx,const char prefix[]);

PetscErrorCode pTatin3dCreateContext(pTatinCtx *ctx);
PetscErrorCode pTatin3dDestroyContext(pTatinCtx *ctx);
PetscErrorCode pTatin3dParseOptions(pTatinCtx ctx);
PetscErrorCode pTatinModelLoad(pTatinCtx ctx);

PetscErrorCode pTatinGetMaterialPoints(pTatinCtx ctx,DataBucket *db,DataEx *de);
PetscErrorCode pTatinGetModel(pTatinCtx ctx,pTatinModel *m);
PetscErrorCode pTatinGetRheology(pTatinCtx ctx,RheologyConstants **r);
PetscErrorCode pTatinGetStokesContext(pTatinCtx ctx,PhysCompStokes *s);
PetscErrorCode pTatinGetMaterialConstants(pTatinCtx ctx,DataBucket *db);

PetscErrorCode pTatin3dContextLoad(pTatinCtx *ctx,const char filename[]);
PetscErrorCode pTatin3dContextSave(pTatinCtx ctx,const char filename[]);
PetscErrorCode pTatin3dCheckpoint(pTatinCtx ctx,Vec X,const char prefix[]);

PetscErrorCode pTatin3d_PhysCompStokesLoad(pTatinCtx user,const char vname[],const char pname[]);
PetscErrorCode pTatin3dRestart(pTatinCtx ctx);

PetscErrorCode pTatinCtxGetModelData(pTatinCtx ctx,const char name[],void **data);
PetscErrorCode pTatinCtxAttachModelData(pTatinCtx ctx,const char name[],void *data);

PetscErrorCode pTatin3dCheckpointManager(pTatinCtx ctx,Vec X);

#endif

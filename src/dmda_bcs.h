
#ifndef __ptatin3d_dmda_bcs_h__
#define __ptatin3d_dmda_bcs_h__


typedef enum { 
	DMDABCList_INTERIOR_LOC=1, 
	DMDABCList_IMIN_LOC,
	DMDABCList_IMAX_LOC,
	DMDABCList_JMIN_LOC,
	DMDABCList_JMAX_LOC,
	DMDABCList_KMIN_LOC,
	DMDABCList_KMAX_LOC 
} DMDABCListConstraintLoc;

typedef enum { BCList_UNCONSTRAINED=-1, BCList_DIRICHLET=-2 } BCListConstraintType;

typedef struct _p_BCList *BCList;

struct _p_BCList {
	DM          dm;
	PetscInt    blocksize;
	PetscInt    N,L; /* L = N * blocksize */
	PetscInt    N_local,L_local;
	PetscScalar *vals_global;
	PetscInt    *dofidx_global;
	PetscInt    *dofidx_local;
	PetscScalar *vals_local;
	/* PetscScalar *scale_local; */
	PetscScalar *scale_global;
	PetscBool   allEmpty;
};

PetscErrorCode BCListIsDirichlet(PetscInt value,PetscBool *flg);
PetscErrorCode BCListInitialize(BCList list);
PetscErrorCode BCListCreate(BCList *list);
PetscErrorCode BCListDestroy(BCList *list);
PetscErrorCode BCListSetSizes(BCList list,PetscInt bs,PetscInt N,PetscInt N_local);
PetscErrorCode BCListUpdateCache(BCList list);
PetscErrorCode BCListInitGlobal(BCList list);
PetscErrorCode BCListGlobalToLocal(BCList list);
PetscErrorCode DMDABCListCreate(DM da,BCList *list);
PetscErrorCode BCListResidualDirichlet(BCList list,Vec X,Vec F);

PetscErrorCode BCListGetGlobalValues(BCList list,PetscInt *n,PetscScalar **vals);
PetscErrorCode BCListRestoreGlobalIndices(BCList list,PetscInt *n,PetscInt **idx);
PetscErrorCode BCListGetGlobalIndices(BCList list,PetscInt *n,PetscInt **idx);
PetscErrorCode BCListGetDofIdx(BCList list,PetscInt *Lg,PetscInt **dofidx_global,PetscInt *Ll,PetscInt **dofidx_local);

PetscErrorCode BCListInsert(BCList list,Vec y);
PetscErrorCode BCListInsertLocal(BCList list,Vec y);


/* 3d */
PetscErrorCode DMDABCListTraverse3d(BCList list,DM da,DMDABCListConstraintLoc doflocation,PetscInt dof_idx,PetscBool (*eval)(PetscScalar*,PetscScalar*,void*),void *ctx);
PetscBool BCListEvaluator_constant( PetscScalar position[], PetscScalar *value, void *ctx );

/* flattened (memory efficient versions) */
PetscErrorCode BCListFlattenedCreate(BCList std,BCList *flat);
PetscErrorCode BCListFlatInsert(BCList list,Vec y);
PetscErrorCode BCListFlatInsertLocal(BCList list,Vec y);
PetscErrorCode BCListFlatResidualDirichlet(BCList list,Vec X,Vec F);


#endif

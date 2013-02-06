
#ifndef __dmda_iterator_h__
#define __dmda_iterator_h__

PetscErrorCode DMDAVecTraverse3d(DM da,Vec X,PetscInt dof_idx,PetscBool (*eval)(PetscScalar*,PetscScalar*,void*),void *ctx);
PetscErrorCode DMDAVecTraverseIJK(DM da,Vec X,PetscInt dof_idx,PetscBool (*eval)(PetscScalar*,PetscInt*,PetscInt*,PetscScalar*,void*),void *ctx);

PetscBool DMDAVecTraverse3d_Constant(PetscScalar pos[],PetscScalar *val,void *ctx);
PetscBool DMDAVecTraverseIJK_Constant(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx);

/* 
 val0 = A[0].x + B[0] 
 val1 = A[1].y + B[1] 
 val2 = A[2].z + B[2] 
 val = val0 + val1 + val2
 */
typedef struct {
	PetscScalar Ox[3];
	PetscScalar A[3],B[3];	
} DMDAVecTraverse3d_InterpCtx;


PetscErrorCode DMDAVecTraverse3d_InterpCtxSetUp_X(DMDAVecTraverse3d_InterpCtx *c,PetscScalar a, PetscScalar b,PetscScalar ox);
PetscErrorCode DMDAVecTraverse3d_InterpCtxSetUp_Y(DMDAVecTraverse3d_InterpCtx *c,PetscScalar a, PetscScalar b,PetscScalar ox);
PetscErrorCode DMDAVecTraverse3d_InterpCtxSetUp_Z(DMDAVecTraverse3d_InterpCtx *c,PetscScalar a, PetscScalar b,PetscScalar ox);
PetscErrorCode DMDAVecTraverse3d_InterpCtxSetUp(DMDAVecTraverse3d_InterpCtx *c,PetscInt dir,PetscScalar a, PetscScalar b,PetscScalar ox);

PetscBool DMDAVecTraverse3d_Interp(PetscScalar pos[],PetscScalar *val,void *ctx);
PetscBool DMDAVecTraverseIJK_Interp(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx);


typedef struct {
	PetscScalar rho;
	PetscScalar grav;
	PetscScalar ref_height;
	PetscInt ref_N;
	PetscScalar surface_pressure;
	DM surface_da;
} DMDAVecTraverse3d_HydrostaticPressureCalcCtx;

PetscBool DMDAVecTraverseIJK_HydroStaticPressure_v1(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx);
PetscBool DMDAVecTraverseIJK_HydroStaticPressure_v2(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx);
PetscBool DMDAVecTraverseIJK_HydroStaticPressure(PetscScalar pos[],PetscInt global_index[],PetscInt local_index[],PetscScalar *val,void *ctx);



#endif

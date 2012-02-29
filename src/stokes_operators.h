

#ifndef __ptatin_stokes_operators_h__
#define __ptatin_stokes_operators_h__

typedef struct _p_MatStokesMF *MatStokesMF;
typedef struct _p_MatA11MF *MatA11MF;

struct _p_MatStokesMF {
	PetscInt    mu,mp,Mu,Mp;
	PetscInt    level;
	PetscInt    ii;
	DM          stokes_pack,daUVW,dap;
	BCList      u_bclist,p_bclist;
	Quadrature  volQ;
	DM          daU;
	IS          isUVW,isU,isV,isW,isP; /* Need the IS's for GetSubMatrix */
	PetscInt    refcnt;
};

struct _p_MatA11MF {
	PetscInt    mu,Mu;
	DM          daUVW;
	BCList      u_bclist;
	Quadrature  volQ;
	DM          daU; /* Optionally need this */
	IS          isUVW; /* Needed for full column space */
	IS          isU,isV,isW; /* Optionally: Need the IS's for GetSubMatrix */
	PetscInt    refcnt;
	/* Not sure I need this at all */
	PetscInt    level;
	PetscInt    ii;
};



PetscErrorCode StokesQ2P1CreateMatrix_Operator(PhysCompStokes user,Mat *B);

PetscErrorCode StokesQ2P1CreateMatrixNest_Operator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B);

PetscErrorCode StokesQ2P1CreateMatrixNest_PCOperator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B);

PetscErrorCode MatStokesMFCreate(MatStokesMF *B);
PetscErrorCode MatA11MFCreate(MatA11MF *B);
PetscErrorCode MatStokesMFSetup(MatStokesMF StkCtx,PhysCompStokes user);
PetscErrorCode MatA11MFSetup(MatA11MF A11Ctx,DM dav,Quadrature volQ,BCList u_bclist);
PetscErrorCode MatStokesMFDestroy(MatStokesMF *B);
PetscErrorCode MatA11MFDestroy(MatA11MF *B);
PetscErrorCode MatStokesMFCopy(MatStokesMF A,MatStokesMF *B);
PetscErrorCode MatA11MFCopy(MatA11MF A,MatA11MF *B);
PetscErrorCode MatCopy_StokesMF_A11MF(MatStokesMF A,MatA11MF *B);

PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11(MatA11MF A11,Mat *A);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A12(MatStokesMF Stk,Mat *A12);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A21(MatStokesMF Stk,Mat *A21);

/* non-optimal preallocation routines for A12,A21 */
PetscErrorCode StokesQ2P1CreateMatrix_A12(PhysCompStokes user,Mat *mat);
PetscErrorCode StokesQ2P1CreateMatrix_A21(PhysCompStokes user,Mat *mat);

PetscErrorCode MatShellGetMatStokesMF(Mat A,MatStokesMF *mf);
PetscErrorCode MatShellGetMatA11MF(Mat A,MatA11MF *mf);

#endif
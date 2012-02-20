

#ifndef __ptatin_stokes_operators_h__
#define __ptatin_stokes_operators_h__

typedef struct _p_MatStokesMF *MatStokesMF;
typedef struct _p_MatA11MF *MatA11MF;


PetscErrorCode StokesQ2P1CreateMatrix_Operator(PhysCompStokes user,Mat *B);

PetscErrorCode StokesQ2P1CreateMatrixNest_Operator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B);

PetscErrorCode StokesQ2P1CreateMatrixNest_PCOperator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B);

PetscErrorCode MatStokesMFCreate(MatStokesMF *B);
PetscErrorCode MatA11MFCreate(MatA11MF *B);
PetscErrorCode MatStokesMFSetup(MatStokesMF StkCtx,PhysCompStokes user);
PetscErrorCode MatA11MFSetup(MatA11MF A11Ctx,DM dav,Quadrature volQ,BCList u_bclist);
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11(MatA11MF A11,Mat *A);


#endif
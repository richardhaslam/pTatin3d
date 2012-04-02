
#ifndef __ptatin3d_stokes_operators_mf_h__
#define __ptatin3d_stokes_operators_mf_h__

/* MatMult wrappers */
PetscErrorCode MFStokesWrapper_A11(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);
PetscErrorCode MFStokesWrapper_A11PC(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);

PetscErrorCode MFStokesWrapper_A(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],PetscScalar Yu[],PetscScalar Yp[]);

PetscErrorCode MFStokesWrapper_A12(Quadrature volQ,DM dau,DM dap,PetscScalar Xp[],PetscScalar Yu[]);

PetscErrorCode MFStokesWrapper_A21(Quadrature volQ,DM dau,DM dap,PetscScalar Xu[],PetscScalar Yp[]);

/* MatGetDiagonal wrappers */
PetscErrorCode MFStokesWrapper_diagA11(Quadrature volQ,DM dau,PetscScalar Yu[]);

#endif

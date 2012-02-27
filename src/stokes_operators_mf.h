
#ifndef __ptatin3d_stokes_operators_mf_h__
#define __ptatin3d_stokes_operators_mf_h__

PetscErrorCode MFStokesWrapper_A11(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[]);

PetscErrorCode MFStokesWrapper_A(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],PetscScalar Yu[],PetscScalar Yp[]);

PetscErrorCode MFStokesWrapper_A12(Quadrature volQ,DM dau,DM dap,PetscScalar Xp[],PetscScalar Yu[]);

PetscErrorCode MFStokesWrapper_A21(Quadrature volQ,DM dau,DM dap,PetscScalar Xu[],PetscScalar Yp[]);

#endif

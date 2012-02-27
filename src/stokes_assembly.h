
#ifndef __ptatin3d_stokes_assembly_h__
#define __ptatin3d_stokes_assembly_h__

PetscErrorCode MatAssemble_StokesA_AUU(Mat A,DM dau,BCList u_dofs,Quadrature volQ);
PetscErrorCode MatAssemble_StokesPC_ScaledMassMatrix(Mat A,DM dau,DM dap,BCList p_bclist,Quadrature volQ);

PetscErrorCode MatAssemble_StokesA_A12(Mat A,DM dau,DM dap,BCList u_bclist,BCList p_bclist,Quadrature volQ);
PetscErrorCode MatAssemble_StokesA_A21(Mat A,DM dau,DM dap,BCList u_bclist,BCList p_bclist,Quadrature volQ);

#endif

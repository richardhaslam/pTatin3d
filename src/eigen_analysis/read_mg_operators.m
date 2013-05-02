%
% Read operators genrated from petsc
% Useful for debugging mg smoother spectrum and testing SLEPc
%

clear all

% Read matrices -----------------------------------
A_Lv0 = PetscBinaryRead('A_coarse.mat');
A_Lv1 = PetscBinaryRead('A_level_1.mat');
%spy(A_Lv0)
A_fine = A_Lv1;


% Check symmetry of coarse grid operators -----------------------------------
fprintf(1,'symm(A0)\n');
max(max(A_Lv0 - A_Lv0'))
fprintf(1,'symm(A1)\n');
max(max(A_Lv1 - A_Lv1'))


% Load ksp operators for each level -----------------------------------
ksp_A_Lv0 = PetscBinaryRead('A_ksp_coarse.mat');
ksp_A_Lv1 = PetscBinaryRead('A_ksp_level_1.mat');


% Load pc operator on fine level
pc_A_fine = PetscBinaryRead('A_pc_mg.mat');
fprintf(1,'symm(pc_fine)\n');
max(max(pc_A_fine - pc_A_fine'))
%tmp = pc_A_fine;
%tmp(tmp<1.0e-1) = 0.0;
%pc_A_fine = tmp;

% Compute eigenvalues -----------------------------------
% This should be compared with 
% right preconditioned operators which are fed into SLEPC
fprintf(1,'eigs(A_fine.pc_A_fine)\n');
eigs(A_fine * pc_A_fine,9,'LM')

% This should be compared with 
% LEFT preconditioned operators which are fed into SLEPC
%fprintf(1,'eigs(pc_A_fine.A_fine)\n');
%eigs(pc_A_fine * A_fine,9,'LM')


% Load ksp operator on fine level -----------------------------------
ksp_A_fine = PetscBinaryRead('A_ksp_mg.mat');



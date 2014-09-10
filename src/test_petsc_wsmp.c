/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      test_petsc_wsmp.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id:$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

static const char help[] = "Test program for PCWSMP\n\n";

#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "ptatin_init.h"

#undef __FUNCT__
#define __FUNCT__ "wssmp_ex1_serial_petsc_lu"
PetscErrorCode wssmp_ex1_serial_petsc_lu(void)
{
	PetscInt m = 9;
	Vec b,x;
	Mat A;
	KSP ksp;
	PC pc;
	PetscInt    rowidx,nz,colidx[9],_c[] = { 0,1,2,3,4,5,6,7,8 };
	PetscScalar vals[9],_v[] = { -3.0, -3.0, -4.0, -3.0, -3.0, -4.0, -4.0, 71.0, -4.0 };
	PetscErrorCode ierr;
	
    PetscPrintf(PETSC_COMM_WORLD,"%s:\n",__FUNCTION__);
	ierr = MatCreate(PETSC_COMM_SELF,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m);CHKERRQ(ierr);
    ierr = MatSetType(A,MATAIJ);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,9,NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A,9,NULL,9,NULL);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);

	/* ------------------Test from wssmp_ex1.f ------------------------- */
	rowidx = 0; nz = 4;
	colidx[0] = 0;     colidx[1] = 2;    colidx[2] = 6;    colidx[3] = 7;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -1.0; vals[3]   = -3.0;
	MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);

	rowidx = 1; nz = 4;
	colidx[0] = 1;     colidx[1] = 2;    colidx[2] = 7;    colidx[3] = 8;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -3.0; vals[3]   = -1.0;
	MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);

	rowidx = 2; nz = 6;
	colidx[0] = 0;     colidx[1] = 1;    colidx[2] = 2;    colidx[3] = 6;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -1.0; vals[2]   = 16.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = -2.0;
	MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);
	
	rowidx = 3; nz = 4;
	colidx[0] = 3;     colidx[1] = 5;    colidx[2] = 6;    colidx[3] = 7;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -1.0; vals[3]   = -3.0;
	MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);

	rowidx = 4; nz = 4;
	colidx[0] = 4;     colidx[1] = 5;    colidx[2] = 7;    colidx[3] = 8;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -3.0; vals[3]   = -1.0;
	MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);
	
	rowidx = 5; nz = 6;
	colidx[0] = 3;     colidx[1] = 4;    colidx[2] = 5;    colidx[3] = 6;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -1.0; vals[2]   = 16.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = -2.0;
	MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);

	rowidx = 6; nz = 6;
	colidx[0] = 0;     colidx[1] = 2;    colidx[2] = 3;    colidx[3] = 5;    colidx[4] = 6;    colidx[5] = 7;
	vals[0]   = -1.0;  vals[1]   = -2.0; vals[2]   = -1.0; vals[3]   = -2.0; vals[4]   = 16.0; vals[5]   = -4.0;
	MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);
	
	rowidx = 7; nz = 9;
	MatSetValues(A,1,&rowidx,nz,_c,_v,INSERT_VALUES);
	
	rowidx = 8; nz = 6;
	colidx[0] = 1;     colidx[1] = 2;    colidx[2] = 4;    colidx[3] = 5;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -2.0; vals[2]   = -1.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = 16.0;
	MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);
	
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = MatView(A,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
	
	ierr = MatGetVecs(A,&b,&x);CHKERRQ(ierr);
	ierr = VecSet(b,1.0);CHKERRQ(ierr);
	/* ------------------------------------------- */
	
	ierr = KSPCreate(PETSC_COMM_SELF,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	ierr = VecView(x,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
	
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "wssmp_ex1_serial_petsc_wsmp"
PetscErrorCode wssmp_ex1_serial_petsc_wsmp(void)
{
	PetscInt m = 9;
	Vec b,x;
	Mat A;
	KSP ksp;
	PC pc;
	PetscInt    rowidx,nz,colidx[9],_c[] = { 0,1,2,3,4,5,6,7,8 };
	PetscScalar vals[9],_v[] = { -3.0, -3.0, -4.0, -3.0, -3.0, -4.0, -4.0, 71.0, -4.0 };
	PetscErrorCode ierr;
	
	PetscPrintf(PETSC_COMM_WORLD,"%s:\n",__FUNCTION__);
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m);CHKERRQ(ierr);
    ierr = MatSetType(A,MATAIJ);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,9,NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A,9,NULL,9,NULL);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	
	/* ------------------Test from wssmp_ex1.f ------------------------- */
	rowidx = 0; nz = 4;
	colidx[0] = 0;     colidx[1] = 2;    colidx[2] = 6;    colidx[3] = 7;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -1.0; vals[3]   = -3.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 1; nz = 4;
	colidx[0] = 1;     colidx[1] = 2;    colidx[2] = 7;    colidx[3] = 8;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -3.0; vals[3]   = -1.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 2; nz = 6;
	colidx[0] = 0;     colidx[1] = 1;    colidx[2] = 2;    colidx[3] = 6;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -1.0; vals[2]   = 16.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = -2.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 3; nz = 4;
	colidx[0] = 3;     colidx[1] = 5;    colidx[2] = 6;    colidx[3] = 7;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -1.0; vals[3]   = -3.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 4; nz = 4;
	colidx[0] = 4;     colidx[1] = 5;    colidx[2] = 7;    colidx[3] = 8;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -3.0; vals[3]   = -1.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 5; nz = 6;
	colidx[0] = 3;     colidx[1] = 4;    colidx[2] = 5;    colidx[3] = 6;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -1.0; vals[2]   = 16.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = -2.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 6; nz = 6;
	colidx[0] = 0;     colidx[1] = 2;    colidx[2] = 3;    colidx[3] = 5;    colidx[4] = 6;    colidx[5] = 7;
	vals[0]   = -1.0;  vals[1]   = -2.0; vals[2]   = -1.0; vals[3]   = -2.0; vals[4]   = 16.0; vals[5]   = -4.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 7; nz = 9;
	ierr = MatSetValues(A,1,&rowidx,nz,_c,_v,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 8; nz = 6;
	colidx[0] = 1;     colidx[1] = 2;    colidx[2] = 4;    colidx[3] = 5;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -2.0; vals[2]   = -1.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = 16.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = MatGetVecs(A,&b,&x);CHKERRQ(ierr);
	ierr = VecSet(b,1.0);CHKERRQ(ierr);
	/* ------------------------------------------- */
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,"wsmp");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"first solve\n");
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"first solve* [identical]\n");
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"second solve [same nonzero pattern]\n");
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"third solve [different nonzero pattern]\n");
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "wssmp_ex1_mpi_petsc_wsmp"
PetscErrorCode wssmp_ex1_mpi_petsc_wsmp(void)
{
	PetscInt m = 9;
	Vec b,x;
	Mat A;
	KSP ksp;
	PC pc;
	PetscInt    rowidx,nz,colidx[9],_c[] = { 0,1,2,3,4,5,6,7,8 };
	PetscScalar vals[9],_v[] = { -3.0, -3.0, -4.0, -3.0, -3.0, -4.0, -4.0, 71.0, -4.0 };
	PetscMPIInt size,rank;
	PetscErrorCode ierr;
	
	PetscPrintf(PETSC_COMM_WORLD,"%s:\n",__FUNCTION__);
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    if (size == 3) { /* special case to match the example in the manual */
        if (rank == 0) {
            ierr = MatSetSizes(A,3,3,m,m);CHKERRQ(ierr);
        }
        if (rank == 1) {
            ierr = MatSetSizes(A,2,2,m,m);CHKERRQ(ierr);
        }
        if (rank == 2) {
            ierr = MatSetSizes(A,4,4,m,m);CHKERRQ(ierr);
        }
    } else {
        ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m);CHKERRQ(ierr);
    }
    ierr = MatSetType(A,MATAIJ);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,9,NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A,9,NULL,9,NULL);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	
	/* ------------------Test from wssmp_ex1.f ------------------------- */
	rowidx = 0; nz = 4;
	colidx[0] = 0;     colidx[1] = 2;    colidx[2] = 6;    colidx[3] = 7;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -1.0; vals[3]   = -3.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 1; nz = 4;
	colidx[0] = 1;     colidx[1] = 2;    colidx[2] = 7;    colidx[3] = 8;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -3.0; vals[3]   = -1.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 2; nz = 6;
	colidx[0] = 0;     colidx[1] = 1;    colidx[2] = 2;    colidx[3] = 6;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -1.0; vals[2]   = 16.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = -2.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 3; nz = 4;
	colidx[0] = 3;     colidx[1] = 5;    colidx[2] = 6;    colidx[3] = 7;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -1.0; vals[3]   = -3.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 4; nz = 4;
	colidx[0] = 4;     colidx[1] = 5;    colidx[2] = 7;    colidx[3] = 8;
	vals[0]   = 14.0;  vals[1]   = -1.0; vals[2]   = -3.0; vals[3]   = -1.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 5; nz = 6;
	colidx[0] = 3;     colidx[1] = 4;    colidx[2] = 5;    colidx[3] = 6;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -1.0; vals[2]   = 16.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = -2.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 6; nz = 6;
	colidx[0] = 0;     colidx[1] = 2;    colidx[2] = 3;    colidx[3] = 5;    colidx[4] = 6;    colidx[5] = 7;
	vals[0]   = -1.0;  vals[1]   = -2.0; vals[2]   = -1.0; vals[3]   = -2.0; vals[4]   = 16.0; vals[5]   = -4.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 7; nz = 9;
	ierr = MatSetValues(A,1,&rowidx,nz,_c,_v,INSERT_VALUES);CHKERRQ(ierr);
	
	rowidx = 8; nz = 6;
	colidx[0] = 1;     colidx[1] = 2;    colidx[2] = 4;    colidx[3] = 5;    colidx[4] = 7;    colidx[5] = 8;
	vals[0]   = -1.0;  vals[1]   = -2.0; vals[2]   = -1.0; vals[3]   = -2.0; vals[4]   = -4.0; vals[5]   = 16.0;
	ierr = MatSetValues(A,1,&rowidx,nz,colidx,vals,INSERT_VALUES);CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = MatGetVecs(A,&b,&x);CHKERRQ(ierr);
	ierr = VecSet(b,1.0);CHKERRQ(ierr);
	/* ------------------------------------------- */
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,"wsmp");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"first solve\n");
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    
	PetscPrintf(PETSC_COMM_WORLD,"first solve* [identical]\n");
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	PetscMPIInt size;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size == 1) {
		//ierr = wssmp_ex1_serial_petsc_lu();CHKERRQ(ierr);
		ierr = wssmp_ex1_serial_petsc_wsmp();CHKERRQ(ierr);
	} else {
		ierr = wssmp_ex1_mpi_petsc_wsmp();CHKERRQ(ierr);
	}
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

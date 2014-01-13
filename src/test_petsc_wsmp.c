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
#define __FUNCT__ "test_wsmp_serial_petsc"
PetscErrorCode test_wsmp_serial_petsc(void)
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
	MatCreate(PETSC_COMM_SELF,&A);
	MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m);
	MatSetFromOptions(A);

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
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	MatView(A,PETSC_VIEWER_STDOUT_SELF);
	
	MatGetVecs(A,&b,&x);
	VecSet(b,1.0);
	/* ------------------------------------------- */
	
	KSPCreate(PETSC_COMM_SELF,&ksp);
	KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
	KSPSetType(ksp,KSPPREONLY);
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCLU);
	KSPSetFromOptions(ksp);

	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	VecView(x,PETSC_VIEWER_STDOUT_SELF);
	
	KSPDestroy(&ksp);
	MatDestroy(&A);
	VecDestroy(&x);
	VecDestroy(&b);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "test_wsmp_pc_wsmp"
PetscErrorCode test_wsmp_pc_wsmp(void)
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
	MatCreate(PETSC_COMM_WORLD,&A);
	MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m);
	MatSetFromOptions(A);
	
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
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	MatGetVecs(A,&b,&x);
	VecSet(b,1.0);
	/* ------------------------------------------- */
	
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
	KSPSetType(ksp,KSPPREONLY);
	KSPGetPC(ksp,&pc);
	PCSetType(pc,"wsmp");
	KSPSetFromOptions(ksp);
	
	printf("first solve\n");
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	VecView(x,PETSC_COMM_WORLD);

	printf("first solve* [identical]\n");
	VecSet(x,0.0);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	VecView(x,PETSC_COMM_WORLD);
	
	printf("second solve [same nonzero pattern]\n");
	VecSet(x,0.0);
	KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	VecView(x,PETSC_COMM_WORLD);

	printf("third solve [different nonzero pattern]\n");
	VecSet(x,0.0);
	KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	VecView(x,PETSC_COMM_WORLD);
	
	KSPDestroy(&ksp);
	MatDestroy(&A);
	VecDestroy(&x);
	VecDestroy(&b);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	PetscMPIInt size;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	if (size == 1) {
		//ierr = test_wsmp_serial_petsc();CHKERRQ(ierr);
		ierr = test_wsmp_pc_wsmp();CHKERRQ(ierr);
	} else {
		ierr = test_wsmp_pc_wsmp();CHKERRQ(ierr);
	}
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

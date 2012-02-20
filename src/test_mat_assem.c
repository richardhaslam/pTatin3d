
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define _GNU_SOURCE

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>


/* ./test_mat_assem.app -A11_mat_type {sbaij,aij} */
PetscErrorCode test_1(void)
{
	PetscErrorCode ierr;
	PetscInt i,j,M,N,s,e;
	PetscScalar *LA_diag;
	Mat A;
	Vec diag;
	PetscBool same;
	
	
	PetscFunctionBegin;

	M = 20;
	N = 20;
	
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(A,"A11_");CHKERRQ(ierr);
	
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);

	ierr = PetscTypeCompare((PetscObject)A,MATSBAIJ,&same);CHKERRQ(ierr);
	if (same) {
		ierr = MatSetOption(A,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
	}
	
	ierr = MatGetOwnershipRange(A,&s,&e);CHKERRQ(ierr);
	for (i=s; i<e; i++) {
		if (i!=0) {   ierr = MatSetValue(A,i,i-1,-1.0,ADD_VALUES);CHKERRQ(ierr); }

								  ierr = MatSetValue(A,i,i,2.0,ADD_VALUES);CHKERRQ(ierr); 

		if (i!=M-1) { ierr = MatSetValue(A,i,i+1,-1.0,ADD_VALUES);CHKERRQ(ierr); }
	}
	
	ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (A,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);

	/* create an array of values to insert - only insert some of them */
	ierr = MatGetVecs(A,&diag,PETSC_NULL);CHKERRQ(ierr);
	ierr = VecSetRandom(diag,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = VecGetArray(diag,&LA_diag);CHKERRQ(ierr);
	for (i=s; i<e; i++) {
		PetscInt idx = i - s;
		
		if ( (i>10) && (i<15) ){
			ierr = MatSetValue(A,i,i,LA_diag[idx],INSERT_VALUES);CHKERRQ(ierr);
		}
	}
	ierr = VecRestoreArray(diag,&LA_diag);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	//ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE);CHKERRQ(ierr);	
	ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
				
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = VecDestroy(&diag);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


/* ./test_mat_assem.app -vel_da_mat_type {sbaij,aij} */
PetscErrorCode test_2(void)
{
	PetscErrorCode ierr;
	PetscInt i,j,M,N,m,s,e;
	PetscScalar *LA_diag;
	Mat B;
	Vec diag;
	DM da;
	PetscBool same;
	
	
	PetscFunctionBegin;

	M = 4;
	N = 5;
	
	ierr = DMDACreate2d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,M,N,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,&da);CHKERRQ(ierr);
	ierr = DMSetOptionsPrefix(da,"vel_");CHKERRQ(ierr);
	ierr = DMSetFromOptions(da);CHKERRQ(ierr);
	
	ierr = DMGetMatrix(da,MATAIJ,&B);CHKERRQ(ierr);

	ierr = PetscTypeCompare((PetscObject)B,MATSBAIJ,&same);CHKERRQ(ierr);
	if (same) {
		ierr = MatSetOption(B,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
	}
	
	ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = MatGetSize(B,&m,PETSC_NULL);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(B,&s,&e);CHKERRQ(ierr);
	for (i=s; i<e; i++) {
		if (i!=0) {   ierr = MatSetValue(B,i,i-1,-1.0,ADD_VALUES);CHKERRQ(ierr); }
		
		ierr = MatSetValue(B,i,i,2.0,ADD_VALUES);CHKERRQ(ierr); 
		
		if (i!=m-1) { ierr = MatSetValue(B,i,i+1,-1.0,ADD_VALUES);CHKERRQ(ierr); }
	}
	
	ierr = MatAssemblyBegin(B,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (B,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
	
	/* create an array of values to insert - only insert some of them */
	ierr = MatGetVecs(B,&diag,PETSC_NULL);CHKERRQ(ierr);
	ierr = VecSetRandom(diag,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = VecGetArray(diag,&LA_diag);CHKERRQ(ierr);
	for (i=s; i<e; i++) {
		PetscInt idx = i - s;
		
		if ( (i>10) && (i<15) ){
			ierr = MatSetValue(B,i,i,LA_diag[idx],INSERT_VALUES);CHKERRQ(ierr);
		}
	}
	ierr = VecRestoreArray(diag,&LA_diag);CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	ierr = VecDestroy(&diag);CHKERRQ(ierr);
	ierr = DMDestroy(&da);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


int main( int argc,char **argv )
{
	PetscErrorCode ierr;

	
	PetscInitialize(&argc,&argv,(char *)0,0);
	
	ierr = test_1();CHKERRQ(ierr);
	ierr = test_2();CHKERRQ(ierr);
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}

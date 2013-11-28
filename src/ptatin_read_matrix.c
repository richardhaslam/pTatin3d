
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_log.h"

#include "ptatin_init.h"

static const char help[] = "Driver to read petsc matrices";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	char filename[PETSC_MAX_PATH_LEN];
	PetscViewer viewer;
	PetscBool flg;
	Mat A;
	Vec b,x;
	KSP ksp;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	flg = PETSC_FALSE;
	PetscOptionsGetString(PETSC_NULL,"-matrix",filename,PETSC_MAX_PATH_LEN-1,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must provide a filename for the matrix to read using -matrix");
	}
	
	PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetType(A,MATAIJ);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	ierr = MatLoad(A,viewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	
	ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
	ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	
	ierr = MatGetVecs(A,&b,&x);CHKERRQ(ierr);
	ierr = VecSetRandom(b,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
		
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

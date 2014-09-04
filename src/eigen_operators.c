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
 **    Filename:      eigen_operators.c
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
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscksp.h"
#include "eigen_operators.h"


typedef struct _p_MatEigenOperator *MatEigenOperator;

struct _p_MatEigenOperator {
	Mat A;
	PC  pc;
	PCSide side;
	Vec t;
};

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MatEigenOperator"
PetscErrorCode MatMult_MatEigenOperator(Mat A,Vec X,Vec Y)
{
	MatEigenOperator  ctx;
  PetscErrorCode    ierr;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);

	switch (ctx->side) {
		case PC_LEFT:
			ierr = MatMult(ctx->A,X,ctx->t);CHKERRQ(ierr);
			ierr = PCApply(ctx->pc,ctx->t,Y);CHKERRQ(ierr);
			break;
		case PC_RIGHT:
			ierr = PCApply(ctx->pc,X,ctx->t);CHKERRQ(ierr);
			ierr = MatMult(ctx->A,ctx->t,Y);CHKERRQ(ierr);
			break;
		case PC_SYMMETRIC:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;
		case PC_SIDE_DEFAULT:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;			
	}
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatDestroy_MatEigenOperator"
PetscErrorCode MatDestroy_MatEigenOperator(Mat A)
{
	MatEigenOperator  ctx;
  PetscErrorCode    ierr;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	
	if (ctx->A) { ierr = MatDestroy(&ctx->A);CHKERRQ(ierr); }
	if (ctx->pc) { ierr = PCDestroy(&ctx->pc);CHKERRQ(ierr); }
	if (ctx->t) { ierr = VecDestroy(&ctx->t);CHKERRQ(ierr); }
	
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatCreateEigenOperatorFromKSPOperators"
PetscErrorCode MatCreateEigenOperatorFromKSPOperators(KSP ksp,Mat *A)
{
	Mat    Aop,Bop,B;
	PC     pc;
	PCSide side;
	PetscInt MA,NA,mA,nA,MB,NB,mB,nB;
	MatEigenOperator  ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	ierr = KSPGetOperators(ksp,&Aop,&Bop);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = KSPGetPCSide(ksp,&side);CHKERRQ(ierr);
	
	ierr = PetscMalloc(sizeof(struct _p_MatEigenOperator),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_MatEigenOperator));CHKERRQ(ierr);
	
	ctx->A   = Aop;     ierr = PetscObjectReference((PetscObject)Aop);CHKERRQ(ierr);
	ctx->pc   = pc;     ierr = PetscObjectReference((PetscObject)pc);CHKERRQ(ierr);
	ctx->side = side;
	
	switch (ctx->side) {
		case PC_LEFT:
			ierr = MatGetVecs(ctx->A,NULL,&ctx->t);CHKERRQ(ierr);
			break;
		case PC_RIGHT:
			ierr = MatGetVecs(ctx->A,NULL,&ctx->t);CHKERRQ(ierr);
			break;
		case PC_SYMMETRIC:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;
		case PC_SIDE_DEFAULT:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;
	}
	
	ierr = MatGetSize(Aop,&MA,&NA);CHKERRQ(ierr);
	ierr = MatGetLocalSize(Aop,&mA,&nA);CHKERRQ(ierr);
	ierr = MatGetSize(Bop,&MB,&NB);CHKERRQ(ierr);
	ierr = MatGetLocalSize(Bop,&mB,&nB);CHKERRQ(ierr);
	
	switch (ctx->side) {
		case PC_LEFT: /* PC.A MB x NB.MA x NA */
			ierr = MatCreateShell(PetscObjectComm((PetscObject)ksp),mB,nA,MB,NA,(void*)ctx,&B);CHKERRQ(ierr);
			break;
		case PC_RIGHT: /* A.PC */
			ierr = MatCreateShell(PetscObjectComm((PetscObject)ksp),mA,nB,MA,NB,(void*)ctx,&B);CHKERRQ(ierr);
			break;
		case PC_SYMMETRIC:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;
		case PC_SIDE_DEFAULT:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only PCSide PC_LEFT and PC_RIGHT are supported");
			break;
	}
	
	ierr = MatShellSetOperation(B,MATOP_MULT,         (void(*)(void))MatMult_MatEigenOperator);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,      (void(*)(void))MatDestroy_MatEigenOperator);CHKERRQ(ierr);
	
	*A = B;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MatEigenOperatorKSP"
PetscErrorCode MatMult_MatEigenOperatorKSP(Mat A,Vec X,Vec Y)
{
	KSP  ctx;
  PetscErrorCode    ierr;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	ierr = KSPSolve(ctx,X,Y);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatCreateEigenOperatorFromKSP"
PetscErrorCode MatCreateEigenOperatorFromKSP(KSP ksp,Mat *A)
{
	Mat    Aop,Bop,B;
	PetscInt MA,NA,mA,nA,MB,NB,mB,nB;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = KSPGetOperators(ksp,&Aop,&Bop);CHKERRQ(ierr);
	
	ierr = MatGetSize(Aop,&MA,&NA);CHKERRQ(ierr);
	ierr = MatGetLocalSize(Aop,&mA,&nA);CHKERRQ(ierr);
	ierr = MatGetSize(Bop,&MB,&NB);CHKERRQ(ierr);
	ierr = MatGetLocalSize(Bop,&mB,&nB);CHKERRQ(ierr);
	
	ierr = MatCreateShell(PetscObjectComm((PetscObject)ksp),mB,nA,MB,NA,(void*)ksp,&B);CHKERRQ(ierr);

	ierr = MatShellSetOperation(B,MATOP_MULT,         (void(*)(void))MatMult_MatEigenOperatorKSP);CHKERRQ(ierr);
	
	*A = B;
	
	PetscFunctionReturn(0);
}


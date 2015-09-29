/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   test_stokes_operators.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

static const char help[] = "Stokes solver using Q2-Pm1 mixed finite elements.\n"
"3D prototype of the (p)ragmatic version of Tatin. (pTatin3d_v0.0)\n\n";


#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "element_utils_q2.h"
#include "stokes_form_function.h"
#include "stokes_assembly.h"
#include "stokes_operators.h"


#undef __FUNCT__  
#define __FUNCT__ "_GenerateTestVector"
PetscErrorCode _GenerateTestVector(DM da,PetscInt dofs,PetscInt index,Vec x)
{
	PetscErrorCode ierr;
	Vec tmp;
	DMDACoor3d ***coors;
	PetscInt i,j,k,mstart,nstart,pstart,m,n,p;
	DM cda;
    ISLocalToGlobalMapping ltog;
	PetscInt NUM_GINDICES;
	const PetscInt *GINDICES;
	
	
	
    ierr = DMGetLocalToGlobalMapping(da, &ltog);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetSize(ltog, &NUM_GINDICES);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetIndices(ltog, &GINDICES);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,&tmp);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(cda,&mstart,&nstart,&pstart,&m,&n,&p);CHKERRQ(ierr);
	
	ierr = DMDAVecGetArray(cda,tmp,&coors);CHKERRQ(ierr);
	for (i=mstart; i<mstart+m; i++) {
		for (j=nstart; j<nstart+n; j++) {
			for (k=pstart; k<pstart+p; k++) {
				PetscInt LIDX = (i-mstart) + (j-nstart)*m + (k-pstart)*m*n;
				PetscInt GIDX = GINDICES[ dofs*LIDX + index ];
				PetscScalar f,XX,YY,ZZ;
				
				XX = coors[k][j][i].x;
				YY = coors[k][j][i].y;
				ZZ = coors[k][j][i].z;
				
				f = XX + YY + ZZ + (PetscScalar)index + 3.0;
				
				ierr = VecSetValue(x,GIDX,f,INSERT_VALUES);CHKERRQ(ierr);
			}}}
	ierr = DMDAVecRestoreArray(cda,tmp,&coors);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingRestoreIndices(ltog, &GINDICES);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "_GenerateTestVectorDAP"
PetscErrorCode _GenerateTestVectorDAP(DM da,PetscInt dofs,PetscInt index,Vec x)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,mstart,nstart,pstart,m,n,p;
    ISLocalToGlobalMapping ltog;
	PetscInt NUM_GINDICES;
	const PetscInt *GINDICES;
	PetscInt M,N,P;
	PetscReal dx,dy,dz;
	
	
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	dx = 2.0/(PetscReal)(M-1);
	dy = 2.0/(PetscReal)(N-1);
	dz = 2.0/(PetscReal)(P-1);
	
	
    ierr = DMGetLocalToGlobalMapping(da, &ltog);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetSize(ltog, &NUM_GINDICES);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetIndices(ltog, &GINDICES);CHKERRQ(ierr);
	
	ierr = DMDAGetGhostCorners(da,&mstart,&nstart,&pstart,&m,&n,&p);CHKERRQ(ierr);
	
	for (i=mstart; i<mstart+m; i++) {
		for (j=nstart; j<nstart+n; j++) {
			for (k=pstart; k<pstart+p; k++) {
				PetscInt LIDX = (i-mstart) + (j-nstart)*m + (k-pstart)*m*n;
				PetscInt GIDX = GINDICES[ dofs*LIDX + index ];
				PetscScalar f,XX,YY,ZZ;
				
				XX = -1.0 + i * dx;
				YY = -1.0 + j * dy;
				ZZ = -1.0 + k * dz;
				
				f = XX + YY + ZZ + (PetscScalar)index + 3.0;
				
				ierr = VecSetValue(x,GIDX,f,INSERT_VALUES);CHKERRQ(ierr);
			}}}
    ierr = ISLocalToGlobalMappingRestoreIndices(ltog, &GINDICES);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ass_A11"
PetscErrorCode ass_A11(PhysCompStokes stk)
{
	PetscErrorCode ierr;
	DM da;
	Mat B;
	PetscBool same;
	Vec x,y;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,stk->mx,stk->my,stk->mz );
	da = stk->dav;

	ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da,&B);CHKERRQ(ierr);
	
	ierr = PetscObjectTypeCompare((PetscObject)B,MATSBAIJ,&same);CHKERRQ(ierr);
	if (same) {
		ierr = MatSetOption(B,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
	}
	
	ierr = MatAssemble_StokesA_AUU(B,da,stk->u_bclist,stk->volQ);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
	
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,2,x);CHKERRQ(ierr);
	
	ierr = MatMult(B,x,y);CHKERRQ(ierr);
	
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ass_B22"
PetscErrorCode ass_B22(PhysCompStokes stk)
{
	PetscErrorCode ierr;
	DM dau,dap;
	Mat B;
	PetscBool same;
	Vec x,y;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,stk->mx,stk->my,stk->mz );
	dau = stk->dav;
	dap = stk->dap;
	
	ierr = DMSetMatType(dap,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(dap,&B);CHKERRQ(ierr);
	
	ierr = PetscObjectTypeCompare((PetscObject)B,MATSBAIJ,&same);CHKERRQ(ierr);
	if (same) {
		ierr = MatSetOption(B,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
	}
	
	ierr = MatAssemble_StokesPC_ScaledMassMatrix(B,dau,dap,stk->p_bclist,stk->volQ);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(dap,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);

	/*
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,2,x);CHKERRQ(ierr);
	
	ierr = MatMult(B,x,y);CHKERRQ(ierr);
	*/
	
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "compare_mf_A11"
PetscErrorCode compare_mf_A11(PhysCompStokes user)
{
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	Mat            Auu,B;
	Vec            x,y,y2;
	DM             da;
	PetscScalar    min,max;
	PetscReal      cmp;
	PetscLogDouble t0,t1;
	double tl,timeMIN,timeMAX;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );

	/* create the mf operators */
	da = user->dav;
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
	
	/* matrix free */
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
	
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,2,x);CHKERRQ(ierr);
	
	PetscTime(&t0);
	ierr = MatMult(Auu,x,y);CHKERRQ(ierr);
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(MF):      time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);


	/* assembled */
	ierr = VecDuplicate(x,&y2);CHKERRQ(ierr);

	ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da,&B);CHKERRQ(ierr);
	PetscTime(&t0);
	ierr = MatAssemble_StokesA_AUU(B,da,user->u_bclist,user->volQ);CHKERRQ(ierr);
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"MatAssemblyA11(ASM): time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
	
	PetscTime(&t0);
	ierr = MatMult(B,x,y2);CHKERRQ(ierr);
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(ASM):     time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);

	/*
	PetscPrintf(PETSC_COMM_WORLD,"y_mfo\n");
	ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"y_asm\n");
	ierr = VecView(y2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	*/
	
	/* compare result */
	ierr = VecDot(y,y,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y.y    = %+1.8e [mfo]\n", cmp );
	ierr = VecDot(y2,y2,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y2.y2  = %+1.8e [asm]\n", cmp );
		
	ierr = VecAXPY(y2,-1.0,y);CHKERRQ(ierr); /* y2 = y2 - y */
	ierr = VecMin(y2,NULL,&min);CHKERRQ(ierr);
	ierr = VecMax(y2,NULL,&max);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  min[A11_mfo.x-A11_asm.x]  = %+1.8e \n", min );
	PetscPrintf(PETSC_COMM_WORLD,"  max[A11_mfo.x-A11_asm.x]  = %+1.8e \n", max );

	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = VecDestroy(&y2);CHKERRQ(ierr);
	
	ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	ierr = MatDestroy(&Auu);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "compare_mf_A21"
PetscErrorCode compare_mf_A21(PhysCompStokes user)
{
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	Mat            A21,B;
	Vec            x,y,y2;
	DM             dav,dap;
	PetscScalar    min,max;
	PetscReal      cmp;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );
	
	/* create the mf operators */
	dav = user->dav;
	dap = user->dap;
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A21(StkCtx,&A21);CHKERRQ(ierr);
	
	/* matrix free */
	ierr = DMCreateGlobalVector(dav,&x);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(dap,&y);CHKERRQ(ierr);
	
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVector(dav,3,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(dav,3,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(dav,3,2,x);CHKERRQ(ierr);
	
	ierr = MatMult(A21,x,y);CHKERRQ(ierr);
	
	/* assembled */
	ierr = VecDuplicate(y,&y2);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrix_A21(user,&B);CHKERRQ(ierr);
	ierr = MatAssemble_StokesA_A21(B,dav,dap,user->u_bclist,user->p_bclist,user->volQ);CHKERRQ(ierr);

	ierr = MatMult(B,x,y2);CHKERRQ(ierr);


	/* compare result */
	ierr = VecDot(y,y,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y.y    = %+1.8e [mfo]\n", cmp );
	ierr = VecDot(y2,y2,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y2.y2  = %+1.8e [asm]\n", cmp );
	
	ierr = VecAXPY(y2,-1.0,y);CHKERRQ(ierr); /* y2 = y2 - y */
	ierr = VecMin(y2,NULL,&min);CHKERRQ(ierr);
	ierr = VecMax(y2,NULL,&max);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  min[A21_mfo.x-A21_asm.x]  = %+1.8e \n", min );
	PetscPrintf(PETSC_COMM_WORLD,"  max[A21_mfo.x-A21_asm.x]  = %+1.8e \n", max );
	
	//	ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = VecDestroy(&y2);CHKERRQ(ierr);
	
	ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	ierr = MatDestroy(&A21);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "compare_mf_A12"
PetscErrorCode compare_mf_A12(PhysCompStokes user)
{
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	Mat            A12,B;
	Vec            x,y,y2;
	DM             dav,dap;
	PetscScalar    min,max;
	PetscReal      cmp;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );
	
	/* create the mf operators */
	dav = user->dav;
	dap = user->dap;
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A12(StkCtx,&A12);CHKERRQ(ierr);
	
	/* matrix free */
	ierr = DMCreateGlobalVector(dap,&x);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(dav,&y);CHKERRQ(ierr);
	
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVectorDAP(dap,4,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVectorDAP(dap,4,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVectorDAP(dap,4,2,x);CHKERRQ(ierr);
	ierr = _GenerateTestVectorDAP(dap,4,3,x);CHKERRQ(ierr);
	
	ierr = MatMult(A12,x,y);CHKERRQ(ierr);
	
	/* assembled */
	ierr = VecDuplicate(y,&y2);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrix_A12(user,&B);CHKERRQ(ierr);
	ierr = MatAssemble_StokesA_A12(B,dav,dap,user->u_bclist,user->p_bclist,user->volQ);CHKERRQ(ierr);
	
	ierr = MatMult(B,x,y2);CHKERRQ(ierr);
	
	
	/* compare result */
	ierr = VecDot(y,y,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y.y    = %+1.8e [mfo]\n", cmp );
	ierr = VecDot(y2,y2,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y2.y2  = %+1.8e [asm]\n", cmp );
	
	ierr = VecAXPY(y2,-1.0,y);CHKERRQ(ierr); /* y2 = y2 - y */
	ierr = VecMin(y2,NULL,&min);CHKERRQ(ierr);
	ierr = VecMax(y2,NULL,&max);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  min[A12_mfo.x-A12_asm.x]  = %+1.8e \n", min );
	PetscPrintf(PETSC_COMM_WORLD,"  max[A12_mfo.x-A12_asm.x]  = %+1.8e \n", max );
	
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = VecDestroy(&y2);CHKERRQ(ierr);
	
	ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	ierr = MatDestroy(&A12);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "compare_mf_A"
PetscErrorCode compare_mf_A(PhysCompStokes user)
{
	Mat            A,B;
	Vec            x,xu,xp,y,y2;
	DM             pack;
	PetscScalar    min,max;
	PetscReal      cmp;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );
	
	/* create the mf operators */
	pack = user->stokes_pack;
	
	/* matrix free */
	ierr = DMCreateGlobalVector(pack,&x);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(pack,x,&xu,&xp);CHKERRQ(ierr);
		ierr = VecSetRandom(xu,NULL);CHKERRQ(ierr);
		ierr = VecSetRandom(xp,NULL);CHKERRQ(ierr);
//		ierr = VecZeroEntries(xp);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(pack,x,&xu,&xp);CHKERRQ(ierr);
	
	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrix_Operator(user,&A);CHKERRQ(ierr);
	ierr = MatMult(A,x,y);CHKERRQ(ierr);

//	PetscPrintf(PETSC_COMM_WORLD,"y_mfo\n");
//	ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	/* assembled */
	ierr = VecDuplicate(x,&y2);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrixNest_Operator(user,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE,&B);CHKERRQ(ierr);
	ierr = MatMult(B,x,y2);CHKERRQ(ierr);
	
//	 PetscPrintf(PETSC_COMM_WORLD,"y_asm\n");
//	 ierr = VecView(y2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	
	/* compare result */
	ierr = VecDot(y,y,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y.y    = %+1.8e [mfo]\n", cmp );
	ierr = VecDot(y2,y2,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y2.y2  = %+1.8e [asm]\n", cmp );
	
	ierr = VecAXPY(y2,-1.0,y);CHKERRQ(ierr); /* y2 = y2 - y */
	ierr = VecMin(y2,NULL,&min);CHKERRQ(ierr);
	ierr = VecMax(y2,NULL,&max);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  min[A_mfo.x-A_asm.x]  = %+1.8e \n", min );
	PetscPrintf(PETSC_COMM_WORLD,"  max[A_mfo.x-A_asm.x]  = %+1.8e \n", max );
	
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = VecDestroy(&y2);CHKERRQ(ierr);
	
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "compare_mf_diagA11"
PetscErrorCode compare_mf_diagA11(PhysCompStokes user)
{
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	Mat            Auu,B;
	Vec            y,y2;
	DM             da;
	PetscScalar    min,max;
	PetscReal      cmp;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );
	
	/* create the mf operators */
	da = user->dav;
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
	
	/* matrix free */
	ierr = DMCreateGlobalVector(da,&y);CHKERRQ(ierr);
	
	ierr = MatGetDiagonal(Auu,y);CHKERRQ(ierr);
	
	/* assembled */
	ierr = VecDuplicate(y,&y2);CHKERRQ(ierr);
	
	ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da,&B);CHKERRQ(ierr);
	ierr = MatAssemble_StokesA_AUU(B,da,user->u_bclist,user->volQ);CHKERRQ(ierr);
	
	ierr = MatGetDiagonal(B,y2);CHKERRQ(ierr);
	
	/* compare result */
	ierr = VecDot(y,y,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y.y    = %+1.8e [mfo]\n", cmp );
	ierr = VecDot(y2,y2,&cmp);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  y2.y2  = %+1.8e [asm]\n", cmp );
	
	ierr = VecAXPY(y2,-1.0,y);CHKERRQ(ierr); /* y2 = y2 - y */
	ierr = VecMin(y2,NULL,&min);CHKERRQ(ierr);
	ierr = VecMax(y2,NULL,&max);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  min[diagA11_mfo-diagA11_asm]  = %+1.8e \n", min );
	PetscPrintf(PETSC_COMM_WORLD,"  max[diagA11_mfo-diagA11_asm]  = %+1.8e \n", max );
	
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = VecDestroy(&y2);CHKERRQ(ierr);
	
	ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	ierr = MatDestroy(&Auu);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "apply_mf_A11"
PetscErrorCode apply_mf_A11(PhysCompStokes user)
{
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	Mat            Auu;
	Vec            x,y;
	DM             da;
	PetscLogDouble t0,t1;
	double         tl,timeMIN,timeMAX;
	PetscInt       ii,iterations;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );
	iterations = 5;
	ierr = PetscOptionsGetInt(NULL,"-iterations",&iterations,0);CHKERRQ(ierr);
	
	/* create the mf operators */
	da = user->dav;
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
	
	/* matrix free */
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
	
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,2,x);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (ii=0; ii<iterations; ii++) {
		ierr = MatMult(Auu,x,y);CHKERRQ(ierr);
	}
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(MF): iterations %.6d     time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",iterations,tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(MF): average               time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl/((double)iterations),100.0*(timeMIN/timeMAX),timeMIN/((double)iterations),timeMAX/((double)iterations));
	
	
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	
	ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	ierr = MatDestroy(&Auu);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "apply_asm_A11"
PetscErrorCode apply_asm_A11(PhysCompStokes user)
{
	Mat            B;
	Vec            x,y;
	DM             da;
	PetscLogDouble t0,t1;
	double         tl,timeMIN,timeMAX;
	PetscInt       ii,iterations;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );
	iterations = 5;
	ierr = PetscOptionsGetInt(NULL,"-iterations",&iterations,0);CHKERRQ(ierr);
	
	/* create the assembled operator */
	da = user->dav;
	
	/* assembled matrix */
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
	
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,2,x);CHKERRQ(ierr);
	
	ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da,&B);CHKERRQ(ierr);
	PetscTime(&t0);
	ierr = MatAssemble_StokesA_AUU(B,da,user->u_bclist,user->volQ);CHKERRQ(ierr);
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"MatAssemblyA11(ASM):                   time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
	
	PetscTime(&t0);
	for (ii=0; ii<iterations; ii++) {
		ierr = MatMult(B,x,y);CHKERRQ(ierr);
	}
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(ASM): iterations %.6d     time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",iterations,tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(ASM): average               time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl/((double)iterations),100.0*(timeMIN/timeMAX),timeMIN/((double)iterations),timeMAX/((double)iterations));
	
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "perform_viscous_solve"
PetscErrorCode perform_viscous_solve(PhysCompStokes user)
{
	Mat            A,B;
	Vec            x,y;
	DM             da;
	PetscLogDouble t0,t1;
	double         tl,timeMIN,timeMAX;
	PetscInt       its;
	PetscInt       ii,iterations;
	KSP            ksp;
	PetscViewer    monviewer;
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	PetscBool      use_mf_A;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );
	iterations = 5;
	ierr = PetscOptionsGetInt(NULL,"-iterations",&iterations,0);CHKERRQ(ierr);
	
	/* create the assembled operator */
	da = user->dav;
	
	
	/* assembled matrix */
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
	
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,2,x);CHKERRQ(ierr);
	
	ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da,&B);CHKERRQ(ierr);
	PetscTime(&t0);
	ierr = MatAssemble_StokesA_AUU(B,da,user->u_bclist,user->volQ);CHKERRQ(ierr);
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr); 
	PetscPrintf(PETSC_COMM_WORLD,"MatAssemblyA11(ASM):                time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);

	use_mf_A = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-use_mf_A11",&use_mf_A,NULL);CHKERRQ(ierr);
	if (use_mf_A) {
		ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
		ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
		ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);
		ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&A);CHKERRQ(ierr);
	} else {
		A = B;
		ierr = PetscObjectReference((PetscObject)B);CHKERRQ(ierr);
	}
	
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,B);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,1.0e-20,PETSC_DEFAULT,PETSC_DEFAULT,30);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    ierr = KSPSetDM(ksp,da);CHKERRQ(ierr);
    ierr = KSPSetDMActive(ksp,PETSC_FALSE);CHKERRQ(ierr);
    
	PetscTime(&t0);
	ierr = KSPSetUp(ksp);CHKERRQ(ierr);
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr); 

	PetscPrintf(PETSC_COMM_WORLD,"KSPSetUpA11:                        time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);

	
	ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)ksp),NULL,&monviewer);CHKERRQ(ierr);
	ierr = KSPMonitorSet(ksp,KSPMonitorDefaultShort,PETSC_VIEWER_STDOUT_WORLD,(PetscErrorCode (*)(void**))PetscViewerDestroy);CHKERRQ(ierr);

	PetscTime(&t0);
	ierr = KSPSolve(ksp,x,y);CHKERRQ(ierr);
	ierr = KSPMonitorCancel(ksp);CHKERRQ(ierr);

	for (ii=1; ii<iterations; ii++) {
		ierr = KSPSolve(ksp,x,y);CHKERRQ(ierr);
	}
	PetscTime(&t1);
	tl = (double)(t1 - t0);
	ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);

	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"KSPSolveA11(its = %d,cycles = %d)   time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",its,iterations,tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
	PetscPrintf(PETSC_COMM_WORLD,"KSPSolveA11: average                time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl/((double)iterations),100.0*(timeMIN/timeMAX),timeMIN/((double)iterations),timeMAX/((double)iterations));

	/*
	{
		PetscScalar min,max;
		PetscReal ydy,resnorm;
		Vec res;
		
		ierr = VecDot(y,y,&ydy);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  y.y     = %+1.8e \n", ydy );
		
		ierr = VecMin(y,NULL,&min);CHKERRQ(ierr);
		ierr = VecMax(y,NULL,&max);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  min[y]  = %+1.8e \n", min );
		PetscPrintf(PETSC_COMM_WORLD,"  max[y]  = %+1.8e \n", max );
		
		ierr = VecDuplicate(x,&res);CHKERRQ(ierr);
		ierr = MatMult(A,y,res);CHKERRQ(ierr);
		ierr = VecAXPY(res,-1.0,x);CHKERRQ(ierr);
		ierr = VecNorm(res,NORM_2,&resnorm);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  |res|     = %+1.8e \n", resnorm );
		
		ierr = VecDestroy(&res);CHKERRQ(ierr);
	}
	*/
	
	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	if (use_mf_A) {
		ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
		ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "perform_viscous_solve_warmup"
PetscErrorCode perform_viscous_solve_warmup(PhysCompStokes user)
{
  Mat            A,B;
  Vec            x,y;
  DM             da;
  PetscLogDouble t0,t1,t0all,t1all;
  double         tl,timeMIN,timeMAX,*time_,*timeMIN_,*timeMAX_;
  PetscInt       its;
  PetscInt       ii,iterations,iterations_warmup;
  KSP            ksp;
  PetscViewer    monviewer;
  MatStokesMF    StkCtx;
  MatA11MF       A11Ctx;
  PetscBool      use_mf_A;
  PetscLogStage stages[2];
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  
  PetscLogStageRegister("Warmup solve",&stages[0]);
  PetscLogStageRegister("Profiled solve",&stages[1]);

  
  PetscPrintf(PETSC_COMM_WORLD,"\n+  Test [%s]: Mesh %D x %D x %D \n", __FUNCT__,user->mx,user->my,user->mz );
  iterations = 5;
  iterations_warmup = 5;
  ierr = PetscOptionsGetInt(NULL,"-iterations",&iterations,0);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,"-iterations_warmup",&iterations_warmup,0);CHKERRQ(ierr);
  
  PetscMalloc(sizeof(double)*iterations,&time_);
  PetscMalloc(sizeof(double)*iterations,&timeMIN_);
  PetscMalloc(sizeof(double)*iterations,&timeMAX_);
  
  /* create the assembled operator */
  da = user->dav;
  
  /* assembled matrix */
  ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
  
  ierr = VecSet(x,0.0);CHKERRQ(ierr);
  ierr = _GenerateTestVector(da,3,0,x);CHKERRQ(ierr);
  ierr = _GenerateTestVector(da,3,1,x);CHKERRQ(ierr);
  ierr = _GenerateTestVector(da,3,2,x);CHKERRQ(ierr);
  
  /* set matrix type via -stk_velocity_dm_mat_type aij */
  ierr = DMCreateMatrix(da,&B);CHKERRQ(ierr);

  /* Assemble the preconditioned operator B11 */
  PetscTime(&t0);
  ierr = MatAssemble_StokesA_AUU(B,da,user->u_bclist,user->volQ);CHKERRQ(ierr);
  PetscTime(&t1);
  tl = (double)(t1 - t0);
  ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"MatAssemblyB11(ASM):                time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
  
  use_mf_A = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,"-use_mf_A11",&use_mf_A,NULL);CHKERRQ(ierr);
  if (use_mf_A) {
    ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
    ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
    ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);
    ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&A);CHKERRQ(ierr);
  } else {
    A = B;
    ierr = PetscObjectReference((PetscObject)B);CHKERRQ(ierr);
  }
  
  
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,B);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.0e-20,PETSC_DEFAULT,PETSC_DEFAULT,30);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  
  ierr = KSPSetDM(ksp,da);CHKERRQ(ierr);
  ierr = KSPSetDMActive(ksp,PETSC_FALSE);CHKERRQ(ierr);
  
  /* Force setup of the Krylov method and preconditioner */
  PetscTime(&t0);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  PetscTime(&t1);
  tl = (double)(t1 - t0);
  ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"KSPSetUpA11:                        time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
  
  
  ierr = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)ksp),NULL,&monviewer);CHKERRQ(ierr);
  ierr = KSPMonitorSet(ksp,KSPMonitorDefaultShort,PETSC_VIEWER_STDOUT_WORLD,(PetscErrorCode (*)(void**))PetscViewerDestroy);CHKERRQ(ierr);
  
  /* Warm up solve stage */
  PetscPrintf(PETSC_COMM_WORLD,"\n------ Warm up stage: Performing %D solves ------\n",iterations_warmup);
  PetscLogStagePush(stages[0]);
  
  PetscTime(&t0);
  ierr = KSPSolve(ksp,x,y);CHKERRQ(ierr);

  /* Cancel the KSP monitor for all subsequent solves */
  ierr = KSPMonitorCancel(ksp);CHKERRQ(ierr);

  /* Assume all solves produce identical iteration counts */
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

  /* Perform remaining warm up solves */
  for (ii=1; ii<iterations_warmup; ii++) {
    ierr = KSPSolve(ksp,x,y);CHKERRQ(ierr);
  }
  PetscTime(&t1);
  tl = (double)(t1 - t0);
  ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  
  PetscPrintf(PETSC_COMM_WORLD,"KSPSolveA11(kspits = %.4D,warmup cycles = %.4D)     time(0) %1.4e (sec) : min/max %1.4e %1.4e (sec) : ratio %1.4e%%\n",its,iterations_warmup,tl,timeMIN,timeMAX,100.0*(timeMIN/timeMAX));
  PetscPrintf(PETSC_COMM_WORLD,"KSPSolveA11: averages (over all runs - all ranks)   time(0) %1.4e (sec) : min/max %1.4e %1.4e (sec) : ratio %1.4e%%\n",tl/((double)iterations_warmup),timeMIN/((double)iterations_warmup),timeMAX/((double)iterations_warmup),100.0*(timeMIN/timeMAX));

  PetscLogStagePop();

  /* Profiled KSP solve */
  PetscPrintf(PETSC_COMM_WORLD,"\n------ Profiling stage: Performing %D solves ------\n",iterations);
  PetscLogStagePush(stages[1]);
  
  PetscTime(&t0all);
  for (ii=0; ii<iterations; ii++) {
    PetscTime(&t0);
  
    ierr = KSPSolve(ksp,x,y);CHKERRQ(ierr);
    
    PetscTime(&t1);
    tl = (double)(t1 - t0);
    time_[ii] = tl;
  }
  PetscTime(&t1all);
  
  /* Assume all solves produce identical iteration counts */
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

  /* Compute profile stats taken over all runs */
  tl = (double)(t1all - t0all);
  ierr = MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"KSPSolveA11(kspits = %.4D,profiled cycles = %.4D)     time(0) %1.4e (sec) : min/max %1.4e %1.4e (sec) : ratio %1.4e%%\n",its,iterations,tl,timeMIN,timeMAX,100.0*(timeMIN/timeMAX));
  PetscPrintf(PETSC_COMM_WORLD,"KSPSolveA11: averages (over all runs - all ranks)     time(0) %1.4e (sec) : min/max %1.4e %1.4e (sec) : ratio %1.4e%%\n",tl/((double)iterations),timeMIN/((double)iterations),timeMAX/((double)iterations),100.0*(timeMIN/timeMAX));

  /* Compute profile stats taken over individual runs */
  ierr = MPI_Allreduce(time_,timeMIN_,iterations,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(time_,timeMAX_,iterations,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"KSPSolveA11: averages (over each run - all ranks)\n");
  for (ii=0; ii<iterations; ii++) {
    PetscPrintf(PETSC_COMM_WORLD,"  [solve %.4D]: min/max %1.4e %1.4e (sec)\n",ii,timeMIN_[ii],timeMAX_[ii]);
  }
  
  PetscLogStagePop();
  
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  PetscFree(time_);
  PetscFree(timeMIN_);
  PetscFree(timeMAX_);
  
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  if (use_mf_A) {
    ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
    ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_assemble_stokes"
PetscErrorCode pTatin3d_assemble_stokes(int argc,char **argv)
{
	PetscErrorCode  ierr;
	DM              dav;
	pTatinCtx       user;
	PetscBool       found;

	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);

	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);

	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	user->pack = user->stokes_ctx->stokes_pack;

	/* fetch some local variables */
	dav           = user->stokes_ctx->dav;
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	/* interpolate point coordinates (needed if mesh was modified) */
	//ierr = QuadratureStokesCoordinateSetUp(user->stokes_ctx->Q,dav);CHKERRQ(ierr);
	//for (e=0; e<QUAD_EDGES; e++) {
	//	ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
	//}
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);


	/* update markers = >> gauss points */
	{
		int               npoints;
		DataField         PField_std;
		DataField         PField_stokes;
		MPntStd           *mp_std;
		MPntPStokes       *mp_stokes;
		
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
		
		DataBucketGetSizes(user->materialpoint_db,&npoints,NULL,NULL);
		mp_std    = PField_std->data; /* should write a function to do this */
		mp_stokes = PField_stokes->data; /* should write a function to do this */
		
		ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(npoints,mp_std,mp_stokes,user->stokes_ctx->dav,user->stokes_ctx->volQ);CHKERRQ(ierr);
	}
	
	/* perform tests */
	PetscPrintf(PETSC_COMM_WORLD,"\n\n\n====================================================================\n");
	
//	ierr = ass_A11(user->stokes_ctx);CHKERRQ(ierr);
//	ierr = ass_B22(user->stokes_ctx);CHKERRQ(ierr);


	found  = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-compare_operators",&found,0);CHKERRQ(ierr);
	if (found) {
		ierr = compare_mf_A11(user->stokes_ctx);CHKERRQ(ierr);

		ierr = compare_mf_A21(user->stokes_ctx);CHKERRQ(ierr);
		ierr = compare_mf_A12(user->stokes_ctx);CHKERRQ(ierr);
		
		ierr = compare_mf_A(user->stokes_ctx);CHKERRQ(ierr);

		ierr = compare_mf_diagA11(user->stokes_ctx);CHKERRQ(ierr);
	}
	
	found  = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-apply_A11mf_operator",&found,NULL);CHKERRQ(ierr);
	if (found) {
		ierr = apply_mf_A11(user->stokes_ctx);CHKERRQ(ierr);
	}

	found  = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-apply_A11asm_operator",&found,NULL);CHKERRQ(ierr);
	if (found) {
		ierr = apply_asm_A11(user->stokes_ctx);CHKERRQ(ierr);		
	}

	found  = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-perform_viscous_solve_A11asm_operator",&found,NULL);CHKERRQ(ierr);
	if (found) {
		ierr = perform_viscous_solve(user->stokes_ctx);CHKERRQ(ierr);		
	}
	
	found  = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,"-perform_viscous_solve_A11_warmup",&found,NULL);CHKERRQ(ierr);
	if (found) {
		ierr = perform_viscous_solve_warmup(user->stokes_ctx);CHKERRQ(ierr);
	}
	
	
	PetscPrintf(PETSC_COMM_WORLD,"\n\n\n====================================================================\n");


	ierr = pTatin3dDestroyContext(&user);

	PetscFunctionReturn(0);
}

extern PetscErrorCode PCCreate_SemiRedundant(PC pc);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	ierr = pTatin3d_assemble_stokes(argc,argv);CHKERRQ(ierr);
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

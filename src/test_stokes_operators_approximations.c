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
 **    Filename:      test_stokes_operators_approximations.c
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
#include "dmda_element_q2p1.h"


#undef __FUNCT__  
#define __FUNCT__ "_GenerateTestVector"
PetscErrorCode _GenerateTestVector(DM da,PetscInt dofs,PetscInt index,Vec x)
{
	PetscErrorCode ierr;
	Vec tmp;
	DMDACoor3d ***coors;
	PetscInt i,j,k,mstart,nstart,pstart,m,n,p;
	DM cda;
	PetscInt NUM_GINDICES, *GINDICES;
	
	
	
	ierr = DMDAGetGlobalIndices( da, &NUM_GINDICES, &GINDICES );CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(da,&tmp);CHKERRQ(ierr);
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
	
	ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "compare_mf_A11"
PetscErrorCode compare_mf_A11(PhysCompStokes user,Quadrature volQ_2x2x2)
{
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	MatStokesMF    StkCtx2x2x2;
	MatA11MF       A11Ctx2x2x2;
	Mat            Auu,Auu2x2x2,B;
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

	/* create the mf operators */
	da = user->dav;
	ierr = MatStokesMFCreate(&StkCtx2x2x2);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx2x2x2,user);CHKERRQ(ierr);
	ierr = MatCopy_StokesMF_A11MF(StkCtx2x2x2,&A11Ctx2x2x2);CHKERRQ(ierr);
	A11Ctx2x2x2->volQ = volQ_2x2x2;
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx2x2x2,&Auu2x2x2);CHKERRQ(ierr);
	
	
	
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
	
	ierr = VecSet(x,0.0);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,0,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,1,x);CHKERRQ(ierr);
	ierr = _GenerateTestVector(da,3,2,x);CHKERRQ(ierr);
	
	/* matrix free */
	PetscGetTime(&t0);
	ierr = MatMult(Auu,x,y);CHKERRQ(ierr);
	PetscGetTime(&t1);
	tl = (double)(t1 - t0);
	MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);	
	MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD); 
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(MF): time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);


	/* matrix free 2x2x2 */
	PetscGetTime(&t0);
	ierr = MatMult(Auu2x2x2,x,y);CHKERRQ(ierr);
	PetscGetTime(&t1);
	tl = (double)(t1 - t0);
	MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);	
	MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD); 
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11_2x2x2(MF): time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
	

	/* assembled */
	ierr = VecDuplicate(x,&y2);CHKERRQ(ierr);

	ierr = DMCreateMatrix(da,MATAIJ,&B);CHKERRQ(ierr);
	PetscGetTime(&t0);
	ierr = MatAssemble_StokesA_AUU(B,da,user->u_bclist,user->volQ);CHKERRQ(ierr);
	PetscGetTime(&t1);
	tl = (double)(t1 - t0);
	MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);	
	MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD); 
	PetscPrintf(PETSC_COMM_WORLD,"MatAssemblyA11(ASM): time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);
	
	PetscGetTime(&t0);
	ierr = MatMult(B,x,y2);CHKERRQ(ierr);
	PetscGetTime(&t1);
	tl = (double)(t1 - t0);
	MPI_Allreduce(&tl,&timeMIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);	
	MPI_Allreduce(&tl,&timeMAX,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD); 
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(ASM): time %1.4e (sec): ratio %1.4e%%: min/max %1.4e %1.4e (sec)\n",tl,100.0*(timeMIN/timeMAX),timeMIN,timeMAX);

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
	ierr = VecMin(y2,PETSC_NULL,&min);CHKERRQ(ierr);
	ierr = VecMax(y2,PETSC_NULL,&max);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"  min[A11_mfo.x-A11_asm.x]  = %+1.8e \n", min );
	PetscPrintf(PETSC_COMM_WORLD,"  max[A11_mfo.x-A11_asm.x]  = %+1.8e \n", max );

	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = VecDestroy(&y2);CHKERRQ(ierr);

	
	
	ierr = MatA11MFDestroy(&A11Ctx2x2x2);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx2x2x2);CHKERRQ(ierr);
	ierr = MatDestroy(&Auu2x2x2);CHKERRQ(ierr);
	
	ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	ierr = MatDestroy(&Auu);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_assemble_stokes"
PetscErrorCode pTatin3d_assemble_stokes(int argc,char **argv)
{
	PetscErrorCode ierr;
	DM              multipys_pack,dav,dap;
	pTatinCtx       user;
	Quadrature			volQ_2x2x2;
	
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
	multipys_pack = user->pack;
	dav           = user->stokes_ctx->dav;
	dap           = user->stokes_ctx->dap;
	
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
		
		DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
		mp_std    = PField_std->data; /* should write a function to do this */
		mp_stokes = PField_stokes->data; /* should write a function to do this */
		
		ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(npoints,mp_std,mp_stokes,user->stokes_ctx->dav,user->stokes_ctx->volQ);CHKERRQ(ierr);
	}


  {
		PetscInt ncells;
		PetscInt lmx,lmy,lmz;
		
		ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		ncells = lmx * lmy * lmz;
		ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,2,ncells,&volQ_2x2x2);CHKERRQ(ierr);
	}
	/* update markers = >> gauss points onto 2^3 quad points */
	{
		int               npoints;
		DataField         PField_std;
		DataField         PField_stokes;
		MPntStd           *mp_std;
		MPntPStokes       *mp_stokes;
		
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
		
		DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
		mp_std    = PField_std->data; /* should write a function to do this */
		mp_stokes = PField_stokes->data; /* should write a function to do this */
		
		ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(npoints,mp_std,mp_stokes,user->stokes_ctx->dav,volQ_2x2x2);CHKERRQ(ierr);
	}
	
	
	
	/* perform tests */
	PetscPrintf(PETSC_COMM_WORLD,"\n\n\n====================================================================\n");
	
	ierr = compare_mf_A11(user->stokes_ctx,volQ_2x2x2);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"\n\n\n====================================================================\n");


	ierr = pTatin3dDestroyContext(&user);

	PetscFunctionReturn(0);
}

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

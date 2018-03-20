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
 **    filename:   stokes_operators_mf.c
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

#include "petsc.h"
#include "petscvec.h" /*I   "petscmat.h"   I*/
#include "petscdm.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "dmda_element_q2p1.h"
#include "quadrature.h"

#include "element_utils_q2.h"
#include "element_utils_q1.h"

#include "stokes_operators.h"

#define PTAT3D_ELEMENT_MF_STANDARD
//#define PTAT3D_ELEMENT_MF_OPTIMIZED
#if defined(PTAT3D_ELEMENT_MF_STANDARD) && defined(PTAT3D_ELEMENT_MF_OPTIMIZED)
#error Only one of PTAT3D_ELEMENT_MF_STANDARD and PTAT3D_ELEMENT_MF_OPTIMIZED may be defined
#endif

#include "stokes_q2p1_mf_operators_def.c"
#include "stokes_q2p1_mf_operators_def_rolled.c"
#include "stokes_q2p1_mf_operators_diag_def.c"

//#define PTAT3D_USE_FORTRAN_MF_KERNELS
//#define PTAT3D_LOG_MF_OP

//#define NO_LOWORDER_OPERATORS
#define ONEPOINTQ_LOWORDER_OPERATORS
//#define ONEPOINTQ_DIAG_LOWORDER_OPERATORS
//#define Q1GEOM_LOWORDER_OPERATORS
//#define AFFINEGEOM_LOWORDER_OPERATORS

/* --- A11 --- */
PetscErrorCode MFStokesWrapper_diagA11(Quadrature volQ,DM dau,PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e;
	const PetscInt *elnidx_u;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	PetscLogDouble t0,t1;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3) );
		
		for (p=0; p<ngp; p++) {
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_diagB11(fac,el_eta[p],NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscTime(&t1);
#ifdef PTAT3D_LOG_MF_OP
	PetscPrintf(PETSC_COMM_WORLD,"MatGetDiagonalA11(MF): %1.4e (sec)\n",t1-t0);
#endif	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_diagA11LowOrder(Quadrature volQ,DM dau,PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e;
	const PetscInt *elnidx_u;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	PetscLogDouble t0,t1;
	PetscReal GNIQ1[NQP][3][8];
	PetscReal xiC[3],GNIC[3][NPE];
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);

	/* prepare Q1 basis */
	for (p=0; p<ngp; p++) {
		P3D_ConstructGNi_Q1_3D(XI[p],GNIQ1[p]);
	}

	xiC[0] = xiC[1] = xiC[2] = 0.0;
	P3D_ConstructGNi_Q2_3D(xiC,GNIC);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
#ifdef NO_LOWORDER_OPERATORS
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
#endif
#ifdef ONEPOINTQ_LOWORDER_OPERATORS
		P3D_evaluate_geometry_elementQ2_1gp(GNIC,ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
#endif
#ifdef ONEPOINTQ_DIAG_LOWORDER_OPERATORS
		P3D_evaluate_geometry_elementQ2_1gp_diagonal(GNIC,ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
#endif
#ifdef Q1GEOM_LOWORDER_OPERATORS
		P3D_evaluate_geometry_elementQ1_appliedQ2(ngp,detJ, GNIQ1, elcoords, GNI,dNudx,dNudy,dNudz );
#endif
#ifdef AFFINEGEOM_LOWORDER_OPERATORS
		/* this appears to be broken */
		//P3D_evaluate_geometry_affine_appliedQ2(ngp,detJ, NULL, elcoords, GNI,dNudx,dNudy,dNudz );
		
		/* double the iteration count */
		P3D_evaluate_geometry_affine2_appliedQ2(ngp,detJ, NULL, elcoords, GNI,dNudx,dNudy,dNudz );
#endif
		
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3) );
		
		for (p=0; p<ngp; p++) {
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_diagB11(fac,el_eta[p],NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscTime(&t1);
#ifdef PTAT3D_LOG_MF_OP	
	PetscPrintf(PETSC_COMM_WORLD,"MatGetDiagonalA11(MF): %1.4e (sec)\n",t1-t0);
#endif	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#ifdef PTAT3D_USE_FORTRAN_MF_KERNELS
extern void f_matmultmf_stokes_mixedfem3d_b11_(double*,double*, double*,double*,double*,double*, double*,double*,double*,double*, double*, double*);
#endif


PetscErrorCode MFStokesWrapper_A11(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e,k;
	const PetscInt *elnidx_u;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal elu[3*Q2_NODES_PER_EL_3D];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	PetscLogDouble t0,t1;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
//	if (ngp==27) { printf("ngp = 27\n");}
//	if (ngp==8) { printf("ngp = 8\n");}
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			PetscInt idx = 3*k;
			ux[k] = elu[  idx];
			uy[k] = elu[++idx];
			uz[k] = elu[++idx];
		}
		
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3) );
		
		
#ifdef PTAT3D_ELEMENT_MF_STANDARD
		for (p=0; p<ngp; p++) {
			//printf("[e=%4d,p=%d] %1.4e \n", e,p,cell_gausspoints[p].eta);
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			#ifdef PTAT3D_USE_FORTRAN_MF_KERNELS
			  f_matmultmf_stokes_mixedfem3d_b11_(&fac,&el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
			#else
			  MatMultMF_Stokes_MixedFEM3d_B11(fac,el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
      #endif
		}
#endif
#ifdef PTAT3D_ELEMENT_MF_OPTIMIZED
		for (p=0; p<ngp; p++) {
			fac = 1.0;
			el_eta[p] = cell_gausspoints[p].eta * WEIGHT[p] * detJ[p];
			
		#ifdef PTAT3D_USE_FORTRAN_MF_KERNELS
			f_matmultmf_stokes_mixedfem3d_b11_(&fac,&el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		#else
			MatMultMF_Stokes_MixedFEM3d_B11(fac,el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		#endif
		}
#endif
		
		
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscTime(&t1);
#ifdef PTAT3D_LOG_MF_OP
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(MF): %1.4e (sec)\n",t1-t0);
	{
		double ops_geom, ops_element_sum, ops_insert, ops, ops_total;
		double flops;
		
		ops_geom        = ( (18.0 + 15.0) * ngp * Q2_NODES_PER_EL_3D   +   (14.0 + 58.0) * ngp ) * nel;
		ops_element_sum = (double)( (1058 * ngp + ngp ) * nel ); /* extra ngp is for the weight x detJ done prior to MatMultMF_Stokes_MixedFEM3d_B11 */
		ops_insert      = (double)( (3 * Q2_NODES_PER_EL_3D ) * nel );
		ops            = ops_geom  +  ops_element_sum  +  ops_insert;
		ierr = MPI_Allreduce(&ops,&ops_total,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
		
		flops = ((double)ops_total)/(t1-t0);
		flops = flops * 1.0e-9;
		
		PetscPrintf(PETSC_COMM_WORLD,"MatMultA11(MF): %1.4e (GFLOPS)\n",flops);
	}	
#endif	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
 Supports using tri-linear element transformations
 
*/
PetscErrorCode MFStokesWrapper_A11PC(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e,k;
	const PetscInt *elnidx_u;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal elu[3*Q2_NODES_PER_EL_3D];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];

	PetscReal GNIQ1[NQP][3][8];
	PetscReal GNIC[3][Q2_NODES_PER_EL_3D],xiC[3];
	/* if you want to cache values - don't do it this way - store it on operator */
	/*
	static PetscReal *cell_detJ;
	static PetscReal *cell_J;
	static PetscBool beenhere;
	*/
	PetscReal fac;
	PetscLogDouble t0,t1;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* prepare Q1 basis */
	for (p=0; p<ngp; p++) {
		P3D_ConstructGNi_Q1_3D(XI[p],GNIQ1[p]);
	}

	xiC[0] = xiC[1] = xiC[2] = 0.0;
	P3D_ConstructGNi_Q2_3D(xiC,GNIC);
	
	/* 
	// update values //
	if (beenhere==PETSC_FALSE) {
		PetscInt ncells,lmx,lmy,lmz;
		
		ierr = DMDAGetSizeElementQ2(dau,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		ncells = lmx * lmy * lmz;
		
		PetscMalloc1(ncells,&cell_detJ);
		PetscMalloc1(3*ncells,&cell_J);
		
		
		beenhere = PETSC_TRUE;
	}
	 */	
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}

#ifdef NO_LOWORDER_OPERATORS
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
#endif
#ifdef ONEPOINTQ_LOWORDER_OPERATORS
		P3D_evaluate_geometry_elementQ2_1gp(GNIC,ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
#endif
#ifdef ONEPOINTQ_DIAG_LOWORDER_OPERATORS
		P3D_evaluate_geometry_elementQ2_1gp_diagonal(GNIC,ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
#endif
#ifdef Q1GEOM_LOWORDER_OPERATORS
		P3D_evaluate_geometry_elementQ1_appliedQ2(ngp,detJ, GNIQ1, elcoords, GNI,dNudx,dNudy,dNudz );
#endif
#ifdef AFFINEGEOM_LOWORDER_OPERATORS
		/* this appears to be broken */
		//P3D_evaluate_geometry_affine_appliedQ2(ngp,detJ, NULL, elcoords, GNI,dNudx,dNudy,dNudz );

		/* double the iteration count */
		P3D_evaluate_geometry_affine2_appliedQ2(ngp,detJ, NULL, elcoords, GNI,dNudx,dNudy,dNudz );
#endif
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3) );

		for (p=0; p<ngp; p++) {
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_B11(fac,el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		}
		/*
		 // this doesn't work //
		{
			p = 13;
			
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_B11(fac,el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		}
		*/
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscTime(&t1);
#ifdef PTAT3D_LOG_MF_OP
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11PC(MF): %1.4e (sec)\n",t1-t0);
#endif	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A11PC_2x2x2(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e,k;
	const PetscInt *elnidx_u;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal elu[3*Q2_NODES_PER_EL_3D];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	
	PetscReal GNIC[3][Q2_NODES_PER_EL_3D],xiC[3];
	PetscReal fac;
	PetscLogDouble t0,t1;
	PetscInt ii;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = 8;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	xiC[0] = xiC[1] = xiC[2] = 0.0;
	P3D_ConstructGNi_Q2_3D(xiC,GNIC);
	
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
		
		P3D_evaluate_geometry_elementQ2_1gp(GNIC,ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);


/*
		// local average of 27 viscosities
		for (ii=0; ii<2; ii++) {
			for (jj=0; jj<2; jj++) {
				for (kk=0; kk<2; kk++) {
					PetscInt si,sj,sk;

					idx_nqp_2x2 = ii + jj*2 + kk*4;
					el_eta[idx_nqp_2x2] = 0.0;
					
					for (si=0; si<2; si++) {
						for (sj=0; sj<2; sj++) {
							for (sk=0; sk<2; sk++) {
								idx_nqp_3x3 = (ii+si) + (jj+sj)*3 + (kk+sk)*9;
								//printf("idx_nqp_3x3 = %d \n", idx_nqp_3x3 );
								
								el_eta[idx_nqp_2x2] += cell_gausspoints[idx_nqp_3x3].eta;
							}
						}
					}
					el_eta[idx_nqp_2x2] = 0.125 * el_eta[idx_nqp_2x2];
				}
			}
			
		}
*/
		
		// average all 27 viscosities
		// converges at least //
		{
			double avg = 0.0;
			for (ii=0; ii<27; ii++) {
				avg += cell_gausspoints[ii].eta;
			}
			for (ii=0; ii<8; ii++) {
				el_eta[ii] = avg/27.0;
			}
		}
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3) );
		
		for (p=0; p<ngp; p++) {
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_B11(fac,el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		}

		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscTime(&t1);
#ifdef PTAT3D_LOG_MF_OP
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11PC(MF): %1.4e (sec)\n",t1-t0);
#endif	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A11PC_1x1x1(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e,k;
	const PetscInt *elnidx_u;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal elu[3*Q2_NODES_PER_EL_3D];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	
	PetscReal GNIC[3][Q2_NODES_PER_EL_3D],xiC[3];
	PetscReal GNI_1[1][3][Q2_NODES_PER_EL_3D];
	PetscReal fac;
	PetscLogDouble t0,t1;
	PetscInt ii;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = 1;
	
	xiC[0] = xiC[1] = xiC[2] = 0.0;
	P3D_ConstructGNi_Q2_3D(xiC,GNIC);

	P3D_ConstructGNi_Q2_3D(xiC,GNI_1[0]);
	
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
		
		P3D_evaluate_geometry_elementQ2_1gp(GNIC,ngp,elcoords,GNI_1, detJ,dNudx,dNudy,dNudz);
		
		{
			double avg = 0.0;
			for (ii=0; ii<27; ii++) {
				avg += cell_gausspoints[ii].eta;
			}

			el_eta[0] = avg/27.0;
		}
		//
		
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3) );
		
		for (p=0; p<ngp; p++) {
			fac       = 8.0 * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_B11(fac,el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscTime(&t1);
#ifdef PTAT3D_LOG_MF_OP
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA11PC(MF): %1.4e (sec)\n",t1-t0);
#endif	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


PetscErrorCode MFStokesWrapper_A(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],PetscScalar Yu[],PetscScalar Yp[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal elu[3*Q2_NODES_PER_EL_3D],elp[P_BASIS_FUNCTIONS];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D + P_BASIS_FUNCTIONS];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscInt  p_el_lidx[P_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = StokesPressure_GetElementLocalIndices(p_el_lidx,(PetscInt*)&elnidx_p[nen_p*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		ierr = DMDAGetScalarElementField(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*e],pfield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
		
		for (p=0; p<ngp; p++) {
			PetscScalar xip[] = { XI[p][0], XI[p][1], XI[p][2] };
			ConstructNi_pressure(xip,elcoords,NIp[p]);
		}
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3 + P_BASIS_FUNCTIONS ) );
		
		
		for (p=0; p<ngp; p++) {
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_B(fac,el_eta[p],ux,uy,uz,elp,NULL,dNudx[p],dNudy[p],dNudz[p],NIp[p],Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu,  vel_el_lidx,&Ye[0]);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(Yp,  p_el_lidx,  &Ye[81]);CHKERRQ(ierr);
	}
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A12(Quadrature volQ,DM dau,DM dap,PetscScalar Xp[],PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal elp[P_BASIS_FUNCTIONS];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscInt  p_el_lidx[P_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = StokesPressure_GetElementLocalIndices(p_el_lidx,(PetscInt*)&elnidx_p[nen_p*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetScalarElementField(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*e],Xp);CHKERRQ(ierr);
		
		for (p=0; p<ngp; p++) {
			PetscScalar xip[] = { XI[p][0], XI[p][1], XI[p][2] };
			ConstructNi_pressure(xip,elcoords,NIp[p]);
		}
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3 ) );
		
		
		for (p=0; p<ngp; p++) {
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_A12(fac,0,0,0,0,elp,NULL,dNudx[p],dNudy[p],dNudz[p],NIp[p],Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu,  vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A21(Quadrature volQ,DM dau,DM dap,PetscScalar Xu[],PetscScalar Yp[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal elu[3*Q2_NODES_PER_EL_3D];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[P_BASIS_FUNCTIONS];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscInt  p_el_lidx[P_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = StokesPressure_GetElementLocalIndices(p_el_lidx,(PetscInt*)&elnidx_p[nen_p*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],Xu);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
		
		for (p=0; p<ngp; p++) {
			PetscScalar xip[] = { XI[p][0], XI[p][1], XI[p][2] };
			ConstructNi_pressure(xip,elcoords,NIp[p]);
		}
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( P_BASIS_FUNCTIONS ) );
		
		
		for (p=0; p<ngp; p++) {
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_A21(fac,0,ux,uy,uz,0,0,dNudx[p],dNudy[p],dNudz[p],NIp[p],Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(Yp,  p_el_lidx,  Ye);CHKERRQ(ierr);
	}
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* --- A11_UPX --- */
PetscErrorCode MFStokesWrapper_diagA11_UPX(Quadrature volQ,DM dau,DM dax,PetscScalar xfield[],PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e;
	const PetscInt *elnidx_u;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	PetscLogDouble t0,t1;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	LA_gcoords = xfield;
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3) );
		
		for (p=0; p<ngp; p++) {
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_diagB11(fac,el_eta[p],NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscTime(&t1);
#ifdef PTAT3D_LOG_MF_OP
	PetscPrintf(PETSC_COMM_WORLD,"MatGetDiagonalA11_UPX(MF): %1.4e (sec)\n",t1-t0);
#endif	
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A11_UPX(Quadrature volQ,DM dau,PetscScalar ufield[],DM dax,PetscScalar xfield[],PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e,k;
	const PetscInt *elnidx_u;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal elu[3*Q2_NODES_PER_EL_3D];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	LA_gcoords = xfield;
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
		
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3) );
		
		for (p=0; p<ngp; p++) {
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_B11(fac,el_eta[p],ux,uy,uz,NULL,NULL,dNudx[p],dNudy[p],dNudz[p],NULL,Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A_UPX(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],DM dax,PetscScalar xfield[],PetscScalar Yu[],PetscScalar Yp[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal el_eta[MAX_QUAD_PNTS];
	PetscReal elu[3*Q2_NODES_PER_EL_3D],elp[P_BASIS_FUNCTIONS];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D + P_BASIS_FUNCTIONS];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscInt  p_el_lidx[P_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	LA_gcoords = xfield;
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = StokesPressure_GetElementLocalIndices(p_el_lidx,(PetscInt*)&elnidx_p[nen_p*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		ierr = DMDAGetScalarElementField(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*e],pfield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
		
		for (p=0; p<ngp; p++) {
			PetscScalar xip[] = { XI[p][0], XI[p][1], XI[p][2] };
			ConstructNi_pressure(xip,elcoords,NIp[p]);
		}
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3 + P_BASIS_FUNCTIONS ) );
		
		
		for (p=0; p<ngp; p++) {
			el_eta[p] = cell_gausspoints[p].eta;
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_B(fac,el_eta[p],ux,uy,uz,elp,NULL,dNudx[p],dNudy[p],dNudz[p],NIp[p],Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu,  vel_el_lidx,&Ye[0]);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(Yp,  p_el_lidx,  &Ye[81]);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A12_UPX(Quadrature volQ,DM dau,DM dap,PetscScalar pfield[],DM dax,PetscScalar xfield[],PetscScalar Yu[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal elp[P_BASIS_FUNCTIONS];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscInt  p_el_lidx[P_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	LA_gcoords = xfield;
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = StokesPressure_GetElementLocalIndices(p_el_lidx,(PetscInt*)&elnidx_p[nen_p*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetScalarElementField(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*e],pfield);CHKERRQ(ierr);
		
		for (p=0; p<ngp; p++) {
			PetscScalar xip[] = { XI[p][0], XI[p][1], XI[p][2] };
			ConstructNi_pressure(xip,elcoords,NIp[p]);
		}
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( Q2_NODES_PER_EL_3D*3 ) );
		
		for (p=0; p<ngp; p++) {
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_A12(fac,0,0,0,0,elp,NULL,dNudx[p],dNudy[p],dNudz[p],NIp[p],Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu,  vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A21_UPX(Quadrature volQ,DM dau,PetscScalar ufield[],DM dap,DM dax,PetscScalar xfield[],PetscScalar Yp[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal elu[3*Q2_NODES_PER_EL_3D];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Ye[P_BASIS_FUNCTIONS];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscInt  p_el_lidx[P_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	LA_gcoords = xfield;
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = StokesPressure_GetElementLocalIndices(p_el_lidx,(PetscInt*)&elnidx_p[nen_p*e]);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
		
		for (p=0; p<ngp; p++) {
			PetscScalar xip[] = { XI[p][0], XI[p][1], XI[p][2] };
			ConstructNi_pressure(xip,elcoords,NIp[p]);
		}
		P3D_evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* initialise element stiffness matrix */
		PetscMemzero( Ye, sizeof(PetscScalar)* ( P_BASIS_FUNCTIONS ) );
		
		
		for (p=0; p<ngp; p++) {
			fac       = WEIGHT[p] * detJ[p];
			
			MatMultMF_Stokes_MixedFEM3d_A21(fac,0,ux,uy,uz,0,0,dNudx[p],dNudy[p],dNudz[p],NIp[p],Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(Yp,  p_el_lidx,  Ye);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

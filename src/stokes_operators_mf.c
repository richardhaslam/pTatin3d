
#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
#include "private/matimpl.h" /*I   "petscmat.h"   I*/

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "dmda_element_q2p1.h"
#include "quadrature.h"

#include "element_utils_q2.h"

#include "stokes_operators.h"

#include "stokes_q2p1_mf_operators_def.c"
#include "stokes_q2p1_mf_operators_diag_def.c"



/* --- A11 --- */

#undef __FUNCT__
#define __FUNCT__ "MFStokesWrapper_A11"
PetscErrorCode MFStokesWrapper_A11(Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
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
	PetscInt  *gidx,elgidx[3*Q2_NODES_PER_EL_3D];
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
	ierr = DMDAGetCoordinateDA( dau, &cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
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
			
			MatMultMF_Stokes_MixedFEM3d_B11(fac,el_eta[p],ux,uy,uz,PETSC_NULL,PETSC_NULL,dNudx[p],dNudy[p],dNudz[p],PETSC_NULL,Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu, vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscGetTime(&t1);
	//PetscPrintf(PETSC_COMM_WORLD,"MatMultA11, = %1.4e (sec)\n",t1-t0);
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFStokesWrapper_A"
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
	PetscInt  *gidx,elgidx[3*Q2_NODES_PER_EL_3D];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	PetscLogDouble t0,t1;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA( dau, &cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
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
			
			MatMultMF_Stokes_MixedFEM3d_B(fac,el_eta[p],ux,uy,uz,elp,PETSC_NULL,dNudx[p],dNudy[p],dNudz[p],NIp[p],Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu,  vel_el_lidx,&Ye[0]);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(Yp,  p_el_lidx,  &Ye[81]);CHKERRQ(ierr);
	}
	
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA, = %1.4e (sec)\n",t1-t0);
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFStokesWrapper_A12"
PetscErrorCode MFStokesWrapper_A12(Quadrature volQ,DM dau,DM dap,PetscScalar Xp[],PetscScalar Yu[])
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
	PetscReal elp[P_BASIS_FUNCTIONS];
	PetscReal Ye[3*Q2_NODES_PER_EL_3D];
	PetscInt  vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscInt  p_el_lidx[P_BASIS_FUNCTIONS];
	PetscInt  *gidx,elgidx[3*Q2_NODES_PER_EL_3D];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	PetscLogDouble t0,t1;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA( dau, &cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
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
			
			MatMultMF_Stokes_MixedFEM3d_A12(fac,0,0,0,0,elp,PETSC_NULL,dNudx[p],dNudy[p],dNudz[p],NIp[p],Ye);
		}
		
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Yu,  vel_el_lidx,Ye);CHKERRQ(ierr);
	}
	
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA12, = %1.4e (sec)\n",t1-t0);
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFStokesWrapper_A21"
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
	PetscInt  *gidx,elgidx[3*Q2_NODES_PER_EL_3D];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscReal fac;
	PetscLogDouble t0,t1;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp = volQ->npoints;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA( dau, &cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
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
	
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"MatMultA21, = %1.4e (sec)\n",t1-t0);
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

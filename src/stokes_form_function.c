
#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
#include "petscsnes.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "swarm_fields.h"
#include "QPntVolCoefStokes_def.h"


#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal_P"
PetscErrorCode FormFunctionLocal_P(PhysCompStokes user,DM dau,PetscScalar u[],DM dap,PetscScalar p[],PetscScalar Rp[])
{	
	PetscErrorCode ierr;
	PetscInt ngp;
	PetscScalar *gp_xi,*gp_weight;
	DM cda;
	Vec gcoords;
	PetscScalar *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e,n,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscScalar elcoords[3*Q2_NODES_PER_EL_3D];
	PetscScalar elu[Q2_NODES_PER_EL_3D*3],elp[P_BASIS_FUNCTIONS];
	PetscScalar ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscScalar xc[Q2_NODES_PER_EL_3D],yc[Q2_NODES_PER_EL_3D],zc[Q2_NODES_PER_EL_3D];
	PetscScalar nx[Q2_NODES_PER_EL_3D],ny[Q2_NODES_PER_EL_3D],nz[Q2_NODES_PER_EL_3D];
	PetscScalar Fe[P_BASIS_FUNCTIONS];
	PetscScalar Be[P_BASIS_FUNCTIONS];
	PetscInt p_el_lidx[P_BASIS_FUNCTIONS];
	PetscInt *gidx,elgidx[3*Q2_NODES_PER_EL_3D];
	PetscInt *gidx_p,elgidx_p[P_BASIS_FUNCTIONS];
	PetscInt nbcs,i,j;
	const PetscInt *ubcidx;
	PetscLogDouble t0,t1;
	BCList u_bclist = user->u_bclist;
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	
	PetscFunctionBegin;
	/* quadrature */
	ngp       = user->volQ->npoints;
	gp_xi     = user->volQ->q_xi_coor;
	gp_weight = user->volQ->q_weight;
	/*
	for (k=0; k<ngp; k++) {
		PetscScalar *xip = &gp_xi[2*k];
		ConstructGNI_Q2_2D(xip,GNI[k]);
	}
	*/
	/* setup for coords */
	ierr = DMDAGetCoordinateDA( dau, &cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	ierr = DMDAGetGlobalIndices(dap,0,&gidx_p);CHKERRQ(ierr);
	/* zero entries */
	//	ierr = VecZeroEntries(Rp);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);

	ierr = VolumeQuadratureGetAllCellData_Stokes(user->volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
	for (e=0;e<nel;e++) {
		ierr = StokesPressure_GetElementLocalIndices(p_el_lidx,(PetscInt*)&elnidx_p[nen_p*e]);CHKERRQ(ierr);
		ierr = VolumeQuadratureGetCellData_Stokes(user->volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		/*
		ierr = DMDAGetElementCoordinatesQ2_2D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		ierr = DMDAGetVectorElementFieldQ2_2D(elu,(PetscInt*)&elnidx_u[nen_u*e],u);CHKERRQ(ierr);
		ierr = DMDAGetScalarElementField_2D(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*e],p);CHKERRQ(ierr);
		
		for( i=0; i<Q2_NODES_PER_EL_3D; i++ ) {
			ux[i] = elu[2*i  ];
			uy[i] = elu[2*i+1];
			
			xc[i] = elcoords[2*i  ];
			yc[i] = elcoords[2*i+1];
		}
		
		for (n=0; n<ngp; n++) {
			PetscScalar *xip = &gp_xi[2*n];
			ConstructNi_pressure(xip,elcoords,NIp[n]);
		}
		*/
		
		/* initialise element stiffness matrix */
		PetscMemzero( Fe, sizeof(PetscScalar)* P_BASIS_FUNCTIONS );
		PetscMemzero( Be, sizeof(PetscScalar)* P_BASIS_FUNCTIONS );
/*		
		for (n=0; n<ngp; n++) {
			PetscScalar  J[2][2], iJ[2][2];
			PetscScalar  div_u_gp;
			PetscScalar fac,J_p,ojp;
*/			
			/* coord transformation */
/*		
			for( i=0; i<2; i++ ) {
				for( j=0; j<2; j++ ) { J[i][j] = 0.0; }
			}
			for( k=0; k<Q2_NODES_PER_EL_2D; k++ ) {
				J[0][0] += GNI[n][0][k] * xc[k] ;
				J[0][1] += GNI[n][0][k] * yc[k] ;
				
				J[1][0] += GNI[n][1][k] * xc[k] ;
				J[1][1] += GNI[n][1][k] * yc[k] ;
			}
			J_p = (J[0][0]*J[1][1]) - (J[0][1]*J[1][0]);
			ojp = 1.0/J_p;
			iJ[0][0] =  J[1][1]*ojp;
			iJ[0][1] = -J[0][1]*ojp;
			iJ[1][0] = -J[1][0]*ojp;
			iJ[1][1] =  J[0][0]*ojp;
*/			
			/* global derivs */
/*		
			for( k=0; k<Q2_NODES_PER_EL_2D; k++ ) {
				nx[k] = iJ[0][0]*GNI[n][0][k] + iJ[0][1]*GNI[n][1][k];
				ny[k] = iJ[1][0]*GNI[n][0][k] + iJ[1][1]*GNI[n][1][k];
			}
			fac = gp_weight[n] * J_p;
*/			
			/* div(u) */
/*		
			div_u_gp = 0.0;
			for( i=0; i<U_BASIS_FUNCTIONS; i++ ) {
				div_u_gp += ( nx[i] * ux[i] + ny[i] * uy[i] );
			}
*/ 
//			div_u_gp = -div_u_gp * fac; /* note the -ve sign here */
/*		
			for( k=0; k<P_BASIS_FUNCTIONS; k++ ) { 
				Fe[k] += NIp[n][k] * div_u_gp; 
			}
*/			
			/* compute any body force terms here */
/*		
			for( k=0; k<P_BASIS_FUNCTIONS; k++ ) { 
				Be[k] = Be[k] + fac * NIp[n][k] * gausspoints[n].Fp;
			}
			
		}
*/		
		/* combine body force with A.x */
/*		
		for( k=0; k<P_BASIS_FUNCTIONS; k++ ) { 
			Fe[k] = Fe[k] - Be[k];
		}
*/		
//		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(Rp, p_el_lidx,Fe);CHKERRQ(ierr);
	}
	PetscGetTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"Assemble Rp, = %1.4e (sec)\n",t1-t0);
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

/*
 Computes r = Ax - b
 SNES will scale by -1, F = -r = b - Ax
 Thus, in OUR function, dirichlet slots become A_ii(x_i - phi)
 In SNES, these become A_ii(phi-x_i), and the updates on the dirichlet slots will be
 A_ii d_i = -F_i 
 = A_ii(phi-x_i)
 Then the update will be 
 x_i^new = x_i + d_i
 = x_i + inv(A_ii) A_ii(phi-x_i)
 = x_i + phi - x_i
 = phi
 */
#undef __FUNCT__  
#define __FUNCT__ "FormFunction_Stokes"
PetscErrorCode FormFunction_Stokes(SNES snes,Vec X,Vec F,void *ctx)
{
  PetscErrorCode    ierr;
  pTatinCtx         ptatin;
  DM                stokes_pack,dau,dap;
  DMDALocalInfo     infou,infop;
  Vec               Uloc,Ploc,FUloc,FPloc;
	Vec               u,p,Fu,Fp;
  PetscScalar       *LA_Uloc,*LA_Ploc;
  PetscScalar       *LA_FUloc,*LA_FPloc;
	PhysCompStokes    stokes;
	
  PetscFunctionBegin;
  
	ptatin      = (pTatinCtx)ctx;
//	stokes      = ptatin->stokes_ctx;
	ierr = pTatinGetStokesContext(ptatin,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dau,&infou);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dap,&infop);CHKERRQ(ierr);
	
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&FUloc,&FPloc);CHKERRQ(ierr);
	
	/* get the local (ghosted) entries for each physics */
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	/* insert boundary conditions into local vectors */
	ierr = BCListInsertLocal(stokes->u_bclist,Uloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(FUloc);CHKERRQ(ierr);
	ierr = VecZeroEntries(FPloc);CHKERRQ(ierr);
	ierr = VecGetArray(FUloc,&LA_FUloc);CHKERRQ(ierr);
	ierr = VecGetArray(FPloc,&LA_FPloc);CHKERRQ(ierr);
	
	/* ======================================== */
	/*         UPDATE NON-LINEARITIES           */
	/* evaluate rheology and rhs using X        */
	/* map marker eta to quadrature points */
	/* map marker force to quadrature points */
	/* ======================================== */
	ierr = pTatin_EvaluateRheologyNonlinearities(ptatin,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
	
	/* Form scaling for momentum */
	//ierr = FormScaling_U_etaMassMatrixDiagonal(user,dau,user->u_bclist);CHKERRQ(ierr);

	/* momentum */
	//ierr = FormFunctionLocal_U(stokes,dau,LA_Uloc,dap,LA_Ploc,LA_FUloc);CHKERRQ(ierr);
	//ierr = FormFunctionLocal_U_tractionBC(stokes,dau,LA_Uloc,dap,LA_Ploc,LA_FUloc);CHKERRQ(ierr);
	
	/* continuity */
	ierr = FormFunctionLocal_P(stokes,dau,LA_Uloc,dap,LA_Ploc,LA_FPloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(FPloc,&LA_FPloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(FUloc,&LA_FUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = VecZeroEntries(F);CHKERRQ(ierr);
	ierr = DMCompositeGather(stokes_pack,F,ADD_VALUES,FUloc,FPloc);CHKERRQ(ierr);
	
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&FUloc,&FPloc);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	/* modify F for the boundary conditions, F_k = scale_k(x_k - phi_k) */
	ierr = DMCompositeGetAccess(stokes_pack,F,&Fu,&Fp);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&u,&p);CHKERRQ(ierr);
	
	ierr = BCListResidualDirichlet(stokes->u_bclist,u,Fu);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&u,&p);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(stokes_pack,F,&Fu,&Fp);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

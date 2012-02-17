
#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
#include "petscsnes.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "swarm_fields.h"
#include "element_type_Q2.h"
#include "dmda_element_q2p1.h"
#include "QPntVolCoefStokes_def.h"
#include "quadrature.h"

#define NQP 27
#define NPE 27

void ConstructNi_P0_3D(const double _xi[],double coords[],double Ni[])
{
	//	printf("hey P0_2D\n");
	Ni[0] = 1.0;
}
void ConstructNi_P1L_3D(const double _xi[],double coords[],double Ni[])
{
	Ni[0] = 1.0;
	Ni[1] = _xi[0];
	Ni[2] = _xi[1];
	Ni[3] = _xi[2];
}
void ConstructNi_P1G_3D(const double _xi[],double coords[],double Ni[])
{
	double Ni_geom[27];
	double _xg[] = {0.0,0.0,0.0};
	int i;
	
	pTatin_ConstructNi_Q2_3D( _xi, Ni_geom );
	for( i=0; i<27; i++ ) {
		_xg[0] = _xg[0] + Ni_geom[i] * coords[3*i  ];
		_xg[1] = _xg[1] + Ni_geom[i] * coords[3*i+1];
		_xg[2] = _xg[2] + Ni_geom[i] * coords[3*i+2];
	}
	ConstructNi_P1L_3D( _xg, coords, Ni );
}
void ConstructNi_P1GRel_3D(const double _xi[],double coords[],double Ni[])
{
	double Ni_geom[27];
	double _xg[] = {0.0,0.0,0.0};
	double avg_x,avg_y,avg_z,Lx,Ly,Lz;
	int i;
	
	pTatin_ConstructNi_Q2_3D( _xi, Ni_geom );
	
	avg_x = avg_y = avg_z = 0.0;
	for( i=0; i<27; i++ ) {
		_xg[0] = _xg[0] + Ni_geom[i] * coords[3*i  ];
		_xg[1] = _xg[1] + Ni_geom[i] * coords[3*i+1];
		_xg[2] = _xg[2] + Ni_geom[i] * coords[3*i+2];
		
		avg_x = avg_x + coords[3*i  ];
		avg_y = avg_y + coords[3*i+1];
		avg_z = avg_z + coords[3*i+2];
	}
	
	avg_x = (1.0/27.0) * avg_x;
	avg_y = (1.0/27.0) * avg_y;
	avg_z = (1.0/27.0) * avg_z;
	/*
	 6--7--8
	 3--4--5
	 0--1--2
	 */
	Lx = coords[3*14  ] - coords[3*12  ];
	Ly = coords[3*16+1] - coords[3*10+1];
	Lz = coords[3*22+2] - coords[3*4+2];
	
	_xg[0] = ( _xg[0] - avg_x ) / Lx ;
	_xg[1] = ( _xg[1] - avg_y ) / Ly ;
	_xg[2] = ( _xg[2] - avg_z ) / Lz ;
	/*	
	 // -1 <= xi,eta <= 1.0 //
	 _xg[0] = 2.0*( _xg[0] - avg_x ) / Lx ;
	 _xg[1] = 2.0*( _xg[1] - avg_y ) / Ly ;
	 */
	ConstructNi_P1L_3D( _xg, coords, Ni );
}

void prepare_elementQ2_3x3(PetscReal WEIGHT[NQP],PetscReal XI[NQP][3],PetscReal NI[NQP][NPE],PetscReal GNI[NQP][3][NPE])
{
	const PetscReal  s   = 0.774596669241483; /* sqrt(3/5) */
	const PetscReal  w_1 = 0.555555555555556; /* 5/9 */
	const PetscReal  w_0 = 0.888888888888889; /* 8/9 */
	PetscReal od_xi[3], od_w[3];
	PetscInt i,j,k,cnt,p;
	
	/* setup gauss quadrature */
	od_xi[0] = -s;		od_xi[1] = 0.0;		od_xi[2] =  s;
	od_w[0] = w_1;		od_w[1] = w_0;		od_w[2] = w_1;
	
	cnt = 0;
	for( k=0; k<3; k++ ) {
		for( j=0; j<3; j++ ) {
			for( i=0; i<3; i++ ) {
				XI[cnt][0]  = od_xi[i];
				XI[cnt][1]  = od_xi[j];
				XI[cnt][2]  = od_xi[k];
				WEIGHT[cnt] = od_w[i] * od_w[j] * od_w[k];
				cnt++;
			}
		}
	} 
	
	/* evaluate derivatives at quadrature points */
	for( p=NQP; p--; ) {
		PetscInt   ix;
		PetscReal  xi = od_xi[p];
		PetscReal  eta = od_xi[p];
		PetscReal  zeta = od_xi[p];
		PetscReal  GNix[3][3], Nix[3][3];
		
		for (ix = 0; ix < 3; ix++) {
			Nix[ix][0] = 0.5 * XI[p][ix] * (XI[p][ix] - 1.0);
			Nix[ix][1] = (1.0 + XI[p][ix]) * (1.0 - XI[p][ix]);
			Nix[ix][2] = 0.5 * (1.0 + XI[p][ix]) * XI[p][ix];
			
			GNix[ix][0] = XI[p][ix] - 0.5;
			GNix[ix][1] = - 2.0 * XI[p][ix];
			GNix[ix][2] = 0.5 * ( 1.0 + 2.0 * XI[p][ix] );
		}
		
		cnt = 0;
		for( k=0; k<3; k++ ) {
			for( j=0; j<3; j++ ) {
				for( i=0; i<3; i++ ) {
					NI[p][cnt]    = Nix[0][i]  *  Nix[1][j]  * Nix[2][k];
					
					GNI[p][0][cnt] = GNix[0][i]  *  Nix[1][j]  *  Nix[2][k];
					GNI[p][1][cnt] =  Nix[0][i]  * GNix[1][j]  *  Nix[2][k];
					GNI[p][2][cnt] =  Nix[0][i]  *  Nix[1][j]  * GNix[2][k];
					
					cnt++;
				}
			}
		}
		
	}
	
}

void prepare_elementQ2_2x2(PetscReal WEIGHT[NQP],PetscReal XI[NQP][3],PetscReal NI[NQP][NPE],PetscReal GNI[NQP][3][NPE])
{
	const PetscReal  s   = 0.577350269189;
	const PetscReal  w0 = 1.0;
	PetscReal od_xi[]={-s,s}, od_w[]={w0,w0};
	PetscInt i,j,k,cnt,p;

	cnt = 0;
	for( k=0; k<2; k++ ) {
		for( j=0; j<2; j++ ) {
			for( i=0; i<2; i++ ) {
				XI[cnt][0]  = od_xi[i];
				XI[cnt][1]  = od_xi[j];
				XI[cnt][2]  = od_xi[k];
				WEIGHT[cnt] = od_w[i] * od_w[j] * od_w[k];
				cnt++;
			}
		}
	} 
	
	/* evaluate derivatives at quadrature points */
	for( p=NQP; p--; ) {
		PetscInt   ix;
		PetscReal  xi = od_xi[p];
		PetscReal  eta = od_xi[p];
		PetscReal  zeta = od_xi[p];
		PetscReal  GNix[3][3], Nix[3][3];
		
		for (ix = 0; ix < 3; ix++) {
			Nix[ix][0] = 0.5 * XI[p][ix] * (XI[p][ix] - 1.0);
			Nix[ix][1] = (1.0 + XI[p][ix]) * (1.0 - XI[p][ix]);
			Nix[ix][2] = 0.5 * (1.0 + XI[p][ix]) * XI[p][ix];
			
			GNix[ix][0] = XI[p][ix] - 0.5;
			GNix[ix][1] = - 2.0 * XI[p][ix];
			GNix[ix][2] = 0.5 * ( 1.0 + 2.0 * XI[p][ix] );
		}
		
		cnt = 0;
		for( k=0; k<3; k++ ) {
			for( j=0; j<3; j++ ) {
				for( i=0; i<3; i++ ) {
					NI[p][cnt]    = Nix[0][i]  *  Nix[1][j]  * Nix[2][k];
					
					GNI[p][0][cnt] = GNix[0][i]  *  Nix[1][j]  *  Nix[2][k];
					GNI[p][1][cnt] =  Nix[0][i]  * GNix[1][j]  *  Nix[2][k];
					GNI[p][2][cnt] =  Nix[0][i]  *  Nix[1][j]  * GNix[2][k];
					
					cnt++;
				}
			}
		}
		
	}

}

void prepare_elementQ2(PetscInt nqp,PetscReal WEIGHT[NQP],PetscReal XI[NQP][3],PetscReal NI[NQP][NPE],PetscReal GNI[NQP][3][NPE])
{
	if (nqp==27) {
		prepare_elementQ2_3x3(WEIGHT,XI,NI,GNI);
	}
	if (nqp==8) {
		prepare_elementQ2_2x2(WEIGHT,XI,NI,GNI);
	}
}

void v1_evaluate_geometry_elementQ2(PetscReal el_coords[NPE*3],PetscReal GNI[NQP][3][NPE],
																 PetscReal detJ[NQP],
																 PetscReal dNudx[NQP][NPE],
																 PetscReal dNudy[NQP][NPE],
																 PetscReal dNudz[NQP][NPE] )
{
	PetscInt k,p;
	PetscReal t4, t6, t8, t10, t12, t14, t17;
	PetscReal J[3][3],iJ[3][3];
	
	for (p=0; p<NQP; p++) {
//
		J[0][0] = J[0][1] = J[0][2] = 0.0;
		J[1][0] = J[1][1] = J[1][2] = 0.0;
		J[2][0] = J[2][1] = J[2][2] = 0.0;
// 
//		memset(J[0],0,sizeof(PetscReal)*3);
//		memset(J[1],0,sizeof(PetscReal)*3);
//		memset(J[2],0,sizeof(PetscReal)*3);
		for (k=0; k<NPE; k++) {
			PetscReal xc = el_coords[3*k+0];
			PetscReal yc = el_coords[3*k+1];
			PetscReal zc = el_coords[3*k+2];
			
			J[0][0] += GNI[p][0][k] * xc ;
			J[0][1] += GNI[p][0][k] * yc ;
			J[0][2] += GNI[p][0][k] * zc ;
			
			J[1][0] += GNI[p][1][k] * xc ;
			J[1][1] += GNI[p][1][k] * yc ;
			J[1][2] += GNI[p][1][k] * zc ;
			
			J[2][0] += GNI[p][2][k] * xc ;
			J[2][1] += GNI[p][2][k] * yc ;
			J[2][2] += GNI[p][2][k] * zc ;
		}
		/* flops = [NQP*NPE] * 18 */
		
		detJ[p] = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) // a
		- J[0][1]*(J[1][0]*J[2][2] + J[1][2]*J[2][0]) 
		+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]); // c
		/* flops = [NQP] * 14 */
		
		t4  = J[2][0] * J[0][1];
		t6  = J[2][0] * J[0][2];
		t8  = J[1][0] * J[0][1];
		t10 = J[1][0] * J[0][2];
		t12 = J[0][0] * J[1][1];
		t14 = J[0][0] * J[1][2]; // 6
		t17 = 0.1e1 / (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1]);  // 12
		
		iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17;  // 4
		iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17; // 5
		iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17;  // 4
		iJ[1][0] = -(-J[2][0] * J[1][2] + J[1][0] * J[2][2]) * t17;// 6
		iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17;                // 4
		iJ[1][2] = -(-t10 + t14) * t17;                            // 4
		iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17; // 5
		iJ[2][1] = -(-t4 + J[0][0] * J[2][1]) * t17;               // 5
		iJ[2][2] = (-t8 + t12) * t17;                              // 3
		/* flops = [NQP] * 58 */
		
		/* shape function derivatives */
		for (k=0; k<NPE; k++) {
			dNudx[p][k] = iJ[0][0]*GNI[p][0][k] + iJ[0][1]*GNI[p][1][k] + iJ[0][2]*GNI[p][2][k];
			
			dNudy[p][k] = iJ[1][0]*GNI[p][0][k] + iJ[1][1]*GNI[p][1][k] + iJ[1][2]*GNI[p][2][k];
			
			dNudz[p][k] = iJ[2][0]*GNI[p][0][k] + iJ[2][1]*GNI[p][1][k] + iJ[2][2]*GNI[p][2][k];
		}
	}
	/* flops = [NQP*NPE] * 15 */
	
	// TOTAL = 18 + 14 + 58 + 15 = 105
}

void evaluate_geometry_elementQ2(PetscInt nqp,PetscReal el_coords[NPE*3],PetscReal GNI[][3][NPE],
																 PetscReal detJ[],
																 PetscReal dNudx[][NPE],
																 PetscReal dNudy[][NPE],
																 PetscReal dNudz[][NPE] )
{
	PetscInt k,p;
	PetscReal t4, t6, t8, t10, t12, t14, t17;
	PetscReal J[3][3],iJ[3][3];
	
	for (p=0; p<nqp; p++) {
		//
		J[0][0] = J[0][1] = J[0][2] = 0.0;
		J[1][0] = J[1][1] = J[1][2] = 0.0;
		J[2][0] = J[2][1] = J[2][2] = 0.0;
		// 
		//		memset(J[0],0,sizeof(PetscReal)*3);
		//		memset(J[1],0,sizeof(PetscReal)*3);
		//		memset(J[2],0,sizeof(PetscReal)*3);
		for (k=0; k<NPE; k++) {
			PetscReal xc = el_coords[3*k+0];
			PetscReal yc = el_coords[3*k+1];
			PetscReal zc = el_coords[3*k+2];
			
			J[0][0] += GNI[p][0][k] * xc ;
			J[0][1] += GNI[p][0][k] * yc ;
			J[0][2] += GNI[p][0][k] * zc ;
			
			J[1][0] += GNI[p][1][k] * xc ;
			J[1][1] += GNI[p][1][k] * yc ;
			J[1][2] += GNI[p][1][k] * zc ;
			
			J[2][0] += GNI[p][2][k] * xc ;
			J[2][1] += GNI[p][2][k] * yc ;
			J[2][2] += GNI[p][2][k] * zc ;
		}
		/* flops = [NQP*NPE] * 18 */
		
		detJ[p] = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) // a
		- J[0][1]*(J[1][0]*J[2][2] + J[1][2]*J[2][0]) 
		+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]); // c
		/* flops = [NQP] * 14 */
		
		t4  = J[2][0] * J[0][1];
		t6  = J[2][0] * J[0][2];
		t8  = J[1][0] * J[0][1];
		t10 = J[1][0] * J[0][2];
		t12 = J[0][0] * J[1][1];
		t14 = J[0][0] * J[1][2]; // 6
		t17 = 0.1e1 / (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1]);  // 12
		
		iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17;  // 4
		iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17; // 5
		iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17;  // 4
		iJ[1][0] = -(-J[2][0] * J[1][2] + J[1][0] * J[2][2]) * t17;// 6
		iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17;                // 4
		iJ[1][2] = -(-t10 + t14) * t17;                            // 4
		iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17; // 5
		iJ[2][1] = -(-t4 + J[0][0] * J[2][1]) * t17;               // 5
		iJ[2][2] = (-t8 + t12) * t17;                              // 3
		/* flops = [NQP] * 58 */
		
		/* shape function derivatives */
		for (k=0; k<NPE; k++) {
			dNudx[p][k] = iJ[0][0]*GNI[p][0][k] + iJ[0][1]*GNI[p][1][k] + iJ[0][2]*GNI[p][2][k];
			
			dNudy[p][k] = iJ[1][0]*GNI[p][0][k] + iJ[1][1]*GNI[p][1][k] + iJ[1][2]*GNI[p][2][k];
			
			dNudz[p][k] = iJ[2][0]*GNI[p][0][k] + iJ[2][1]*GNI[p][1][k] + iJ[2][2]*GNI[p][2][k];
		}
	}
	/* flops = [NQP*NPE] * 15 */
	
	// TOTAL = 18 + 14 + 58 + 15 = 105
}

#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal_profile"
PetscErrorCode FormFunctionLocal_profile(PhysCompStokes user,DM dau,PetscScalar u[],DM dap,PetscScalar p[],PetscScalar Rp[])
{	
	PetscErrorCode ierr;
	DM cda;
	Vec gcoords;
	PetscScalar *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e,n,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscScalar elcoords[3*Q2_NODES_PER_EL_3D];
	PetscScalar nx[Q2_NODES_PER_EL_3D],ny[Q2_NODES_PER_EL_3D],nz[Q2_NODES_PER_EL_3D];
	PetscInt p_el_lidx[P_BASIS_FUNCTIONS];
	PetscInt *gidx,elgidx[3*Q2_NODES_PER_EL_3D];
	PetscInt *gidx_p,elgidx_p[P_BASIS_FUNCTIONS];
	PetscLogDouble t0,t1;
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	
	PetscFunctionBegin;

	/* quadrature */
	prepare_elementQ2(NQP,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA( dau, &cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	ierr = DMDAGetGlobalIndices(dap,0,&gidx_p);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(user->volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
	for (e=0;e<nel;e++) {
		ierr = VolumeQuadratureGetCellData_Stokes(user->volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);

		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);

		/* coord transformation */
		v1_evaluate_geometry_elementQ2(elcoords,GNI, detJ,dNudx,dNudy,dNudz);
//		evaluate_geometry_elementQ2(NQP,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
	
	}
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"Element geometry = %1.4e (sec)\n",t1-t0);
	{
		double ops = 105.0 * (double)( nel * NQP * NPE );
		PetscPrintf( PETSC_COMM_WORLD, " %1.4e (sec) : %1.4e (Mflop/s) : %d elements \n", t1-t0, 1.0e-6*ops/(t1-t0), nel );
	}	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal_U"
PetscErrorCode FormFunctionLocal_U(PhysCompStokes user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],PetscScalar Ru[])
{	
	PetscErrorCode ierr;
	PetscInt p,ngp;
	DM cda;
	Vec gcoords;
	PetscReal *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscReal elcoords[3*Q2_NODES_PER_EL_3D],el_eta[MAX_QUAD_PNTS];
	PetscReal elu[3*Q2_NODES_PER_EL_3D],elp[P_BASIS_FUNCTIONS];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal Fe[3*Q2_NODES_PER_EL_3D],Be[3*Q2_NODES_PER_EL_3D];
	PetscInt vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscInt *gidx,elgidx[3*Q2_NODES_PER_EL_3D];
	PetscLogDouble t0,t1;
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscFunctionBegin;
	
	/* quadrature */
	ngp = user->volQ->npoints;
	prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA( dau, &cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);

	ierr = VolumeQuadratureGetAllCellData_Stokes(user->volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
	for (e=0;e<nel;e++) {
		PetscScalar int_P, int_divu;
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);

		ierr = VolumeQuadratureGetCellData_Stokes(user->volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
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
		
		/* initialise element stiffness matrix */
		PetscMemzero( Fe, sizeof(PetscScalar)* Q2_NODES_PER_EL_3D*3 );
		PetscMemzero( Be, sizeof(PetscScalar)* Q2_NODES_PER_EL_3D*3 );
			
		evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		/* evaluate the viscosity */
		for (p=0; p<ngp; p++) {
			el_eta[p] = cell_gausspoints[p].eta;
			//printf("  [e=%d] eta = %1.4e \n", e, el_eta[p] );
		}
		
		int_P = int_divu = 0.0;
//#if 0
		for (p=0; p<ngp; p++) {
			PetscScalar  exx,eyy,ezz,exy,exz,eyz;
			PetscScalar  sxx,syy,szz,sxy,sxz,syz,pressure_gp;
			PetscScalar fac,d1,d2,J_p;
			
			fac = WEIGHT[p] * detJ[p];
			
			/* pressure */
			pressure_gp = 0.0;
			for (k=0; k<P_BASIS_FUNCTIONS; k++) {
				pressure_gp += NIp[p][k] * elp[k];
			}
			pressure_gp = pressure_gp * fac;
			
			/* strain rate, B u */
			exx=0.0;  eyy=0.0;  ezz=0.0;
			exy=0.0;  exz=0.0;  eyz=0.0;
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
				exx += (dNudx[p][k] * ux[k]);
				eyy += (dNudy[p][k] * uy[k]);
				ezz += (dNudz[p][k] * uz[k]);
				
				exy += dNudy[p][k] * ux[k] + dNudx[p][k] * uy[k];
				exz += dNudz[p][k] * ux[k] + dNudx[p][k] * uz[k];
				eyz += dNudz[p][k] * uy[k] + dNudy[p][k] * uz[k];
			}
			int_P    += fac * pressure_gp;
			int_divu += fac * (exx + eyy + ezz);
			
			/* constituive */
			d1 = 2.0 * el_eta[p] * fac;
			d2 =       el_eta[p] * fac;
			
			/* stress */
			sxx = d1 * exx;
			syy = d1 * eyy;
			szz = d1 * ezz;

			sxy = d2 * exy;
			sxz = d2 * exz;
			syz = d2 * eyz;
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) { 
				Fe[3*k]   += dNudx[p][k]*(sxx-pressure_gp)   + dNudy[p][k]*sxy   + dNudz[p][k]*sxz; 
				Fe[3*k+1] += dNudy[p][k]*(syy-pressure_gp)   + dNudx[p][k]*sxy   + dNudz[p][k]*syz; 
				Fe[3*k+2] += dNudz[p][k]*(szz-pressure_gp)   + dNudx[p][k]*sxz   + dNudy[p][k]*syz; 
			}

			/* compute any body force terms here */
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) { 
				Be[3*k  ] = Be[3*k  ] + fac * NI[p][k] * cell_gausspoints[p].Fu[0];
				Be[3*k+1] = Be[3*k+1] + fac * NI[p][k] * cell_gausspoints[p].Fu[1];
				Be[3*k+2] = Be[3*k+2] + fac * NI[p][k] * cell_gausspoints[p].Fu[2];
			}
		}
//#endif
#if 0		
		for (p=0; p<ngp; p++) {
			PetscReal  exx,eyy,ezz,exy,exz,eyz;
			PetscReal  sxx,syy,szz,sxy,sxz,syz,pressure_gp;
			PetscReal  fac,d1,d2,J_p;
			PetscReal  *dNudx_p = dNudx[p];
			PetscReal  *dNudy_p = dNudy[p];
			PetscReal  *dNudz_p = dNudz[p];
			PetscReal  *NI_p = NI[p];
			PetscReal  *NIp_p = NIp[p];
			
			
			fac = WEIGHT[p] * detJ[p];
			
			/* pressure */
			pressure_gp = 0.0;
			for (k=0; k<P_BASIS_FUNCTIONS; k++) {
				pressure_gp += NIp_p[k] * elp[k];
			}
			pressure_gp = pressure_gp * fac;
			
			/* strain rate, B u */
			exx=0.0;  eyy=0.0;  ezz=0.0;
			exy=0.0;  exz=0.0;  eyz=0.0;
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
				exx += (dNudx_p[k] * ux[k]);
				eyy += (dNudy_p[k] * uy[k]);
				ezz += (dNudz_p[k] * uz[k]);
				
				exy += dNudy_p[k] * ux[k] + dNudx_p[k] * uy[k];
				exz += dNudz_p[k] * ux[k] + dNudx_p[k] * uz[k];
				eyz += dNudz_p[k] * uy[k] + dNudy_p[k] * uz[k];
			}
			int_P    += fac * pressure_gp;
			int_divu += fac * (exx + eyy + ezz);
			
			/* constituive */
			d1 = 2.0 * el_eta[p] * fac;
			d2 =       el_eta[p] * fac;
			
			/* stress */
			sxx = d1 * exx;
			syy = d1 * eyy;
			szz = d1 * ezz;
			
			sxy = d2 * exy;
			sxz = d2 * exz;
			syz = d2 * eyz;
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) { 
				Fe[3*k]   += dNudx_p[k]*(sxx-pressure_gp)   + dNudy_p[k]*sxy   + dNudz_p[k]*sxz; 
				Fe[3*k+1] += dNudy_p[k]*(syy-pressure_gp)   + dNudx_p[k]*sxy   + dNudz_p[k]*syz; 
				Fe[3*k+2] += dNudz_p[k]*(szz-pressure_gp)   + dNudx_p[k]*sxz   + dNudy_p[k]*syz; 
			}
			
			/* compute any body force terms here */
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) { 
				Be[3*k  ] = Be[3*k  ] + fac * NI_p[k] * cell_gausspoints[p].Fu[0];
				Be[3*k+1] = Be[3*k+1] + fac * NI_p[k] * cell_gausspoints[p].Fu[1];
				Be[3*k+2] = Be[3*k+2] + fac * NI_p[k] * cell_gausspoints[p].Fu[2];
			}
		}
#endif		
		
		/* combine body force with A.x */
		for (k=0; k<Q2_NODES_PER_EL_3D; k++) { 
			Fe[3*k  ] = Fe[3*k  ] - Be[3*k  ];
			Fe[3*k+1] = Fe[3*k+1] - Be[3*k+1];
			Fe[3*k+2] = Fe[3*k+2] - Be[3*k+2];
		}
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(Ru, vel_el_lidx,Fe);CHKERRQ(ierr);
	}
	
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"Assemble Ru, = %1.4e (sec)\n",t1-t0);
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal_P"
PetscErrorCode FormFunctionLocal_P(PhysCompStokes user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],PetscScalar Rp[])
{	
	PetscErrorCode ierr;
	PetscInt ngp;
	PetscScalar *gp_xi,*gp_weight;
	DM cda;
	Vec gcoords;
	PetscScalar *LA_gcoords;
	PetscInt nel,nen_u,nen_p,e,p,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscScalar elcoords[3*Q2_NODES_PER_EL_3D];
	PetscScalar elu[3*Q2_NODES_PER_EL_3D],elp[P_BASIS_FUNCTIONS];
	PetscScalar ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscScalar xc[Q2_NODES_PER_EL_3D],yc[Q2_NODES_PER_EL_3D],zc[Q2_NODES_PER_EL_3D];
	PetscScalar nx[Q2_NODES_PER_EL_3D],ny[Q2_NODES_PER_EL_3D],nz[Q2_NODES_PER_EL_3D];
	PetscScalar Fe[P_BASIS_FUNCTIONS];
	PetscScalar Be[P_BASIS_FUNCTIONS];
	PetscInt p_el_lidx[P_BASIS_FUNCTIONS];
	PetscInt *gidx,elgidx[3*Q2_NODES_PER_EL_3D];
	PetscInt *gidx_p,elgidx_p[P_BASIS_FUNCTIONS];
	PetscLogDouble t0,t1;
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
	PetscReal detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscFunctionBegin;
	
	/* quadrature */
	ngp = user->volQ->npoints;
	prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);	
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA( dau, &cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates( dau,&gcoords );CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	ierr = DMDAGetGlobalIndices(dap,0,&gidx_p);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);

	ierr = VolumeQuadratureGetAllCellData_Stokes(user->volQ,&all_gausspoints);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
	for (e=0;e<nel;e++) {
		ierr = StokesPressure_GetElementLocalIndices(p_el_lidx,(PetscInt*)&elnidx_p[nen_p*e]);CHKERRQ(ierr);
		ierr = VolumeQuadratureGetCellData_Stokes(user->volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);

		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],ufield);CHKERRQ(ierr);
		ierr = DMDAGetScalarElementField(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*e],pfield);CHKERRQ(ierr);
		
		for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uy[k] = elu[3*k+2];
		}
		
		for (p=0; p<ngp; p++) {
			PetscScalar xip[] = { XI[p][0], XI[p][1], XI[p][2] };
			ConstructNi_pressure(xip,elcoords,NIp[p]);
		}
		
		/* initialise element stiffness matrix */
		PetscMemzero( Fe, sizeof(PetscScalar)* P_BASIS_FUNCTIONS );
		PetscMemzero( Be, sizeof(PetscScalar)* P_BASIS_FUNCTIONS );

		evaluate_geometry_elementQ2(ngp,elcoords,GNI, detJ,dNudx,dNudy,dNudz);
		
		for (p=0; p<ngp; p++) {
			PetscScalar  div_u_gp;
			PetscScalar  fac,J_p,ojp;

			fac = WEIGHT[p] * detJ[p];
			/* div(u) */

			div_u_gp = 0.0;
			for( k=0; k<U_BASIS_FUNCTIONS; k++) {
				div_u_gp += ( dNudx[p][k] * ux[k] + dNudy[p][k] * uy[k] + dNudz[p][k] * uz[k] );
			}
			div_u_gp = -div_u_gp * fac; /* note the -ve sign here */
			
			for( k=0; k<P_BASIS_FUNCTIONS; k++ ) { 
				Fe[k] += NIp[p][k] * div_u_gp; 
			}
			
			/* compute any body force terms here */
			for( k=0; k<P_BASIS_FUNCTIONS; k++ ) { 
				Be[k] = Be[k] + fac * NIp[p][k] * cell_gausspoints[p].Fp;
			}
			
		}
		/* combine body force with A.x */
		for( k=0; k<P_BASIS_FUNCTIONS; k++ ) { 
			Fe[k] = Fe[k] - Be[k];
		}
		ierr = DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(Rp, p_el_lidx,Fe);CHKERRQ(ierr);
	}
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"Assemble Rp, = %1.4e (sec)\n",t1-t0);
	
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
	ierr = FormFunctionLocal_U(stokes,dau,LA_Uloc,dap,LA_Ploc,LA_FUloc);CHKERRQ(ierr);
	//ierr = FormFunctionLocal_U_tractionBC(stokes,dau,LA_Uloc,dap,LA_Ploc,LA_FUloc);CHKERRQ(ierr);
	
	/* continuity */
	ierr = FormFunctionLocal_P(stokes,dau,LA_Uloc,dap,LA_Ploc,LA_FPloc);CHKERRQ(ierr);

	/* profiling */
	ierr = FormFunctionLocal_profile(stokes,dau,LA_Uloc,dap,LA_Ploc,LA_FPloc);CHKERRQ(ierr);
	
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

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
 **    Filename:      element_utils_q2.c
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

#include "element_utils_q2.h"

#define ELEMENT_OPERATION_STANDARD
//#define ELEMENT_OPERATION_OPTIMIZED


void P3D_ConstructNi_Q2_2D(PetscReal _xi[],PetscReal Ni[])
{
	PetscReal basis_NI[2][3];
	PetscInt  i,j,d,cnt;
	
	for (d=0; d<2; d++) {
		PetscReal xi = _xi[d];
		
		basis_NI[d][0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )
	}
	
	cnt = 0;
	for (j=0; j<3; j++) {
		for (i=0; i<3; i++) {
			Ni[cnt] = basis_NI[0][i] * basis_NI[1][j];
			cnt++;
		}
	}
}

/*
 ConstructGNi()
 + Defines the basis function derivatives wrt the local coordinates (xi,eta,zeta).
 + _1D implies the spatial degress of freedom expected in the _xi[] array
 */
void P3D_ConstructGNi_Q2_2D(PetscReal _xi[],PetscReal GNi[2][Q2_NODES_PER_EL_2D])
{
	PetscReal basis_NI[2][3];
	PetscReal basis_GNI[2][3];
	PetscInt  i,j,d,cnt;
	
	for (d=0; d<2; d++) {
		PetscReal xi = _xi[d];
		
		basis_NI[d][0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )
		
		basis_GNI[d][0] = 0.5 * ( 2.0*xi - 1.0 );
		basis_GNI[d][1] = - 2.0*xi;
		basis_GNI[d][2] = 0.5 * ( 2.0*xi + 1.0 );
	}
	
	cnt = 0;
	for (j=0; j<3; j++) {
		for (i=0; i<3; i++) {
			
			GNi[0][cnt] = basis_GNI[0][i]  *  basis_NI[1][j];
			GNi[1][cnt] =  basis_NI[0][i]  * basis_GNI[1][j];
			
			cnt++;
		}
	}
}

void P3D_ConstructNi_Q2_3D(PetscReal _xi[],PetscReal Ni[])
{
	PetscInt i,j,k,d,cnt;
	PetscReal basis_NI[3][3];
	
	for( d=0; d<3; d++ ) {
		double xi = _xi[d];
		
		basis_NI[d][0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )
	}
	
	cnt = 0;
	for( k=0; k<3; k++ ) {
		for( j=0; j<3; j++ ) {
			for( i=0; i<3; i++ ) {
				Ni[cnt] = basis_NI[0][i] * basis_NI[1][j] * basis_NI[2][k];
				cnt++;
			}
		}
	}
}

void P3D_ConstructGNi_Q2_3D(PetscReal _xi[],PetscReal GNi[3][Q2_NODES_PER_EL_3D])
{
	double basis_NI[3][3];
	double basis_GNI[3][3];
	int i,j,k,d,cnt;
	
	
	for( d=0; d<3; d++ ) {
		double xi = _xi[d];
		
		basis_NI[d][0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )
		
		basis_GNI[d][0] = 0.5 * ( 2.0*xi - 1.0 );
		basis_GNI[d][1] = - 2.0*xi;
		basis_GNI[d][2] = 0.5 * ( 2.0*xi + 1.0 );
	}
	
	cnt = 0;
	for( k=0; k<3; k++ ) {
		for( j=0; j<3; j++ ) {
			for( i=0; i<3; i++ ) {
				
				GNi[0][cnt] = basis_GNI[0][i]  *  basis_NI[1][j]  *  basis_NI[2][k];
				GNi[1][cnt] = basis_NI[0][i]   *  basis_GNI[1][j] *  basis_NI[2][k];
				GNi[2][cnt] = basis_NI[0][i]   *  basis_NI[1][j]  *  basis_GNI[2][k];
				
				cnt++;
			}
		}
	}
}

void P3D_ConstructNi_P0_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[])
{
	Ni[0] = 1.0;
}

void P3D_ConstructNi_P1L_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[])
{
	Ni[0] = 1.0;
	Ni[1] = _xi[0];
	Ni[2] = _xi[1];
	Ni[3] = _xi[2];
}

void P3D_ConstructNi_P1G_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[])
{
	PetscReal Ni_geom[27];
	PetscReal _xg[] = {0.0,0.0,0.0};
	PetscInt i;
	
	P3D_ConstructNi_Q2_3D( _xi, Ni_geom );
	for( i=0; i<27; i++ ) {
		_xg[0] = _xg[0] + Ni_geom[i] * coords[3*i  ];
		_xg[1] = _xg[1] + Ni_geom[i] * coords[3*i+1];
		_xg[2] = _xg[2] + Ni_geom[i] * coords[3*i+2];
	}
	P3D_ConstructNi_P1L_3D( _xg, coords, Ni );
}

void P3D_ConstructNi_P1GRel_3D(PetscReal _xi[],PetscReal coords[],PetscReal Ni[])
{
	PetscReal Ni_geom[27];
	PetscReal _xg[] = {0.0,0.0,0.0};
	PetscReal avg_x,avg_y,avg_z,Lx,Ly,Lz;
	PetscInt i;
	
	P3D_ConstructNi_Q2_3D( _xi, Ni_geom );
	
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
	P3D_ConstructNi_P1L_3D( _xg, coords, Ni );
}

void P3D_prepare_elementQ2_3x3(PetscReal WEIGHT[],PetscReal XI[][3],PetscReal NI[][NPE],PetscReal GNI[][3][NPE])
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

void P3D_prepare_elementQ2_2x2(PetscReal WEIGHT[],PetscReal XI[][3],PetscReal NI[][NPE],PetscReal GNI[][3][NPE])
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

void P3D_prepare_elementQ2(PetscInt nqp,PetscReal WEIGHT[],PetscReal XI[][3],PetscReal NI[][NPE],PetscReal GNI[][3][NPE])
{
	if (nqp==27) {
		P3D_prepare_elementQ2_3x3(WEIGHT,XI,NI,GNI);
	}
	if (nqp==8) {
		P3D_prepare_elementQ2_2x2(WEIGHT,XI,NI,GNI);
	}
}

#ifdef ELEMENT_OPERATION_STANDARD
void P3D_evaluate_geometry_elementQ2(PetscInt nqp,PetscReal el_coords[NPE*3],PetscReal GNI[][3][NPE],
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
	
	// TOTAL = [NQP*NPE]*18 + [NQP]*(14 + 58) + [NQP*NPE]*15
	PetscLogFlops(NQP*NPE*18 + NQP*(14+58) + NQP*NPE*15);
}
#endif

#ifdef ELEMENT_OPERATION_OPTIMIZED
void P3D_evaluate_geometry_elementQ2(PetscInt nqp,
                                     PetscReal el_coords[NPE*3],
                                     PetscReal GNI[][3][NPE],
                                     PetscReal detJ[],
                                     PetscReal dNudx[][NPE],
                                     PetscReal dNudy[][NPE],
                                     PetscReal dNudz[][NPE] )
{
	PetscInt k,p;
	PetscReal t4, t6, t8, t10, t12, t14, t17;
	PetscReal iJ[3][3];
	
	for (p=0; p<nqp; p++) {
		PetscReal J[][3] = { {0.0,0.0,0.0} , {0.0,0.0,0.0} , {0.0,0.0,0.0} };
		
		for (k=0; k<27; k++) {
			PetscInt  idx = 3*k;
			PetscReal xc = el_coords[  idx];
			PetscReal yc = el_coords[++idx];
			PetscReal zc = el_coords[++idx];
			
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

/* <option_A> */
		detJ[p] = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])  // a
		        - J[0][1]*(J[1][0]*J[2][2] + J[1][2]*J[2][0])  // b
		        + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]); // c
		/* flops = [NQP] * 14 */
		
		t4  = J[2][0] * J[0][1];
		t6  = J[2][0] * J[0][2];
		t8  = J[1][0] * J[0][1];
		t10 = J[1][0] * J[0][2];
		t12 = J[0][0] * J[1][1];
		t14 = J[0][0] * J[1][2]; // [NQP] * 6
		
		t17 = 1.0/detJ[p]; // [NQP] * 1
/* </option_A> */
#ifdef COMPUTE_DETERMINANT_OPTION_B
	 /* 
		 NOTE: Option A and option B are identical, however B uses 2 less operations.
		 We elected for option A as more unit tests fail with option B variant due to small floating point differences 
		*/
		t4  = J[2][0] * J[0][1];
		t6  = J[2][0] * J[0][2];
		t8  = J[1][0] * J[0][1];
		t10 = J[1][0] * J[0][2];
		t12 = J[0][0] * J[1][1];
		t14 = J[0][0] * J[1][2]; // [NQP] * 6
		detJ[p] = (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1]);  // [NQP] * 12
		t17 = 1.0/detJ[p]; // [NQP] * 1
#endif
		
		iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17;  // 4
		iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17; // 5
		iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17;  // 4
		iJ[1][0] = (J[2][0] * J[1][2] - J[1][0] * J[2][2]) * t17;  // 4
		iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17;                // 4
		iJ[1][2] = (t10 - t14) * t17;                              // 2
		iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17; // 5
		iJ[2][1] = (t4 - J[0][0] * J[2][1]) * t17;                 // 3
		iJ[2][2] = (-t8 + t12) * t17;                              // 3
		/* flops = [NQP] * 34 */
		
		/* shape function derivatives */
		for (k=0; k<27; k++) {
			dNudx[p][k] = iJ[0][0]*GNI[p][0][k] + iJ[0][1]*GNI[p][1][k] + iJ[0][2]*GNI[p][2][k];
			dNudy[p][k] = iJ[1][0]*GNI[p][0][k] + iJ[1][1]*GNI[p][1][k] + iJ[1][2]*GNI[p][2][k];
			dNudz[p][k] = iJ[2][0]*GNI[p][0][k] + iJ[2][1]*GNI[p][1][k] + iJ[2][2]*GNI[p][2][k];
		}
	}
	/* flops = [NQP*NPE] * 15 */
	
	// TOTAL = [NQP*NPE]*18 + [NQP]*(14 + 6 + 1 + 34) + [NQP*NPE]*15
	PetscLogFlops(NQP*NPE*18 + NQP*(14+6+1+34) + NQP*NPE*15);
}
#endif

void P3D_evaluate_geometry_elementQ2_1gp(PetscReal GNI_centre[3][NPE],
                                         PetscInt nqp,
                                         PetscReal el_coords[NPE*3],
                                         PetscReal GNI[][3][NPE],
                                         PetscReal detJ[],
                                         PetscReal dNudx[][NPE],
                                         PetscReal dNudy[][NPE],
                                         PetscReal dNudz[][NPE] )
{
	PetscInt k,p;
	PetscReal t4, t6, t8, t10, t12, t14, t17;
	PetscReal J[3][3],iJ[3][3],detJp;

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
		
		J[0][0] += GNI_centre[0][k] * xc ;
		J[0][1] += GNI_centre[0][k] * yc ;
		J[0][2] += GNI_centre[0][k] * zc ;
		
		J[1][0] += GNI_centre[1][k] * xc ;
		J[1][1] += GNI_centre[1][k] * yc ;
		J[1][2] += GNI_centre[1][k] * zc ;
		
		J[2][0] += GNI_centre[2][k] * xc ;
		J[2][1] += GNI_centre[2][k] * yc ;
		J[2][2] += GNI_centre[2][k] * zc ;
	}
	
	detJp = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) // a
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
		
	for (p=0; p<nqp; p++) {
		detJ[p] = detJp;
		
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

void P3D_evaluate_geometry_elementQ2_1gp_diagonal(PetscReal GNI_centre[3][NPE],
                                                  PetscInt nqp,
                                                  PetscReal el_coords[NPE*3],
                                                  PetscReal GNI[][3][NPE],
                                                  PetscReal detJ[],
                                                  PetscReal dNudx[][NPE],
                                                  PetscReal dNudy[][NPE],
                                                  PetscReal dNudz[][NPE] )
{
	PetscInt k,p;
	PetscReal J[3],iJ[3],detJp;
	
	J[0] = J[1] = J[2] = 0.0;

	for (k=0; k<NPE; k++) {
		PetscReal xc = el_coords[3*k+0];
		PetscReal yc = el_coords[3*k+1];
		PetscReal zc = el_coords[3*k+2];
		
		J[0] += GNI_centre[0][k] * xc ;
		J[1] += GNI_centre[1][k] * yc ;
		J[2] += GNI_centre[2][k] * zc ;
	}
	/* flops = NPE * 3 */
	
	detJp = J[0]*J[1]*J[2];
	/* flops = 2 */
	
	iJ[0] = 1.0/J[0];
	iJ[1] = 1.0/J[1];
	iJ[2] = 1.0/J[2];
	/* flops = 3 */
	
	for (p=0; p<nqp; p++) {
		detJ[p] = detJp;
		
		/* shape function derivatives */
		for (k=0; k<NPE; k++) {
			dNudx[p][k] = iJ[0]*GNI[p][0][k];
			
			dNudy[p][k] =                         iJ[1]*GNI[p][1][k];
			
			dNudz[p][k] =                                                 iJ[2]*GNI[p][2][k];
		}
	}
	/* flops = [NQP*NPE] * 3 */
	
	// TOTAL = NPE*3 + 2 + 3 + [NQP*NPE] * 3 = 2303
}

void P3D_evaluate_global_derivatives_Q2(PetscReal el_coords[NPE*3],PetscReal GNI[3][NPE],
                                        PetscReal dNudx[NPE],
                                        PetscReal dNudy[NPE],
                                        PetscReal dNudz[NPE] )
{
	PetscInt k;
	PetscReal t4, t6, t8, t10, t12, t14, t17;
	PetscReal J[3][3],iJ[3][3];
	
	J[0][0] = J[0][1] = J[0][2] = 0.0;
	J[1][0] = J[1][1] = J[1][2] = 0.0;
	J[2][0] = J[2][1] = J[2][2] = 0.0;
	
	for (k=0; k<NPE; k++) {
		PetscReal xc = el_coords[3*k+0];
		PetscReal yc = el_coords[3*k+1];
		PetscReal zc = el_coords[3*k+2];
		
		J[0][0] += GNI[0][k] * xc;
		J[0][1] += GNI[0][k] * yc;
		J[0][2] += GNI[0][k] * zc;
		
		J[1][0] += GNI[1][k] * xc;
		J[1][1] += GNI[1][k] * yc;
		J[1][2] += GNI[1][k] * zc;
		
		J[2][0] += GNI[2][k] * xc;
		J[2][1] += GNI[2][k] * yc;
		J[2][2] += GNI[2][k] * zc;
	}
	/* flops = [NPE] * 18 */
	
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
		dNudx[k] = iJ[0][0]*GNI[0][k] + iJ[0][1]*GNI[1][k] + iJ[0][2]*GNI[2][k];
		
		dNudy[k] = iJ[1][0]*GNI[0][k] + iJ[1][1]*GNI[1][k] + iJ[1][2]*GNI[2][k];
		
		dNudz[k] = iJ[2][0]*GNI[0][k] + iJ[2][1]*GNI[1][k] + iJ[2][2]*GNI[2][k];
	}
	/* flops = [NPE] * 15 */
	
	// TOTAL = 18 + 58 + 15 = 105
}

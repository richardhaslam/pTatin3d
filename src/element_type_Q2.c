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
 **    filename:   element_type_Q2.c
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
/*
 *  surface_integration_Q2_3D.c
 *
 *
 *  Created by Dave May on 4/12/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "element_type_Q2.h"
#include "ptatin3d_defs.h"

void pTatin_ConstructNi_Q2_3D( double _xi[], double Ni[] )
{
	int i,j,k,d,cnt;
	double basis_NI[3][3];

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

void pTatinConstructGNi_Q2_3D( double _xi[], double GNi[3][Q2_NODES_PER_EL_3D] )
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

void malloc2d( const int M, const int N, int ***A2d )
{
	int *a;
	int **a2d;
	int i;

	a = malloc( sizeof(int) * M*N );
	memset( a, 0, sizeof(int)*M*N );

	a2d = malloc( sizeof(int*) * M );
	for( i=0; i<M; i++ ) {
		a2d[i] = &a[i*N];
	}

	*A2d = a2d;
}

void free2d( int ***_a )
{
	int **a2 = *_a;
	free(a2[0]);
	free(a2);

	*_a = NULL;
}

/*
 generate_all_face_ids()
 + Defines the local node indices which define the face
*/
void generate_all_face_ids_Q2_1D( int **fid )
{
	fid[0][0] = 0;
	fid[1][0] = 2;
}

/*

 6--7--8
 |     |
 3  4  5
 |     |
 0--1--2

 // edge 2-5-8 //
 // edge 0-3-6 //
 // edge 6-7-8 //
 // edge 0-1-2 //

 order faces from
 +xi,-xi, +eta,-eta

 I DON'T UNDERSTAND THIS EDGE FLIPPING SHIT AT THE MOMENT....

 */
void generate_all_face_ids_Q2_2D( int **fid )
{
//	fid[0][0] = 8;		fid[0][1] = 5;		fid[0][2] = 2;
	fid[0][0] = 2;		fid[0][1] = 5;		fid[0][2] = 8; // Pxi needs to be reversed
	fid[1][0] = 0;		fid[1][1] = 3;		fid[1][2] = 6;
	fid[2][0] = 6;		fid[2][1] = 7;		fid[2][2] = 8;
//	fid[3][0] = 2;		fid[3][1] = 1;		fid[3][2] = 0;
	fid[3][0] = 0;		fid[3][1] = 1;		fid[3][2] = 2; // Neta needs to be reversed
}

/*

 faces:

 20-11-02 23-14-05 26-17-08 => +x
 00-09-18 03-12-21 06-15-24 => -x

 24-25-26 15-16-17 06-07-08 => +y
 00-01-02 09-10-11 18-19-20 => -y

 18-19-20 21-22-23 24-25-26 => +z
 02-01-00 05-04-03 08-07-06 => -z


 order faces from
 +xi,-xi, +eta,-eta, +zeta,-zeta

 */
void generate_all_face_ids_Q2_3D( int **fid )
{
	/*
	int fid_px[] = { 20,11,2, 23,14,5, 26,17,8 };
	int fid_mx[] = { 0,9,18,  3,12,21, 6,15,24 };

	int fid_pe[] = { 24,25,26, 15,16,17, 6,7,8 };
	int fid_me[] = { 0,1,2, 9,10,11, 18,19,20 };

	int fid_pz[] = { 18,19,20, 21,22,23, 24,25,26 };
	int fid_mz[] = { 2,1,0, 5,4,3, 8,7,6 };
	*/

	int fid_px[] = { 2,5,8, 11,14,17, 20,23,26 };
	int fid_mx[] = { 6,3,0, 15,12,9,  24,21,18 };

	int fid_pe[] = { 8,7,6, 17,16,15, 26,25,24 };
	int fid_me[] = { 0,1,2, 9,10,11, 18,19,20 };

	int fid_pz[] = { 18,19,20, 21,22,23, 24,25,26 };
	int fid_mz[] = { 2,1,0, 5,4,3, 8,7,6 };


	memcpy( fid[0], fid_px, sizeof(int)*9 );
	memcpy( fid[1], fid_mx, sizeof(int)*9 );

	memcpy( fid[2], fid_pe, sizeof(int)*9 );
	memcpy( fid[3], fid_me, sizeof(int)*9 );

	memcpy( fid[4], fid_pz, sizeof(int)*9 );
	memcpy( fid[5], fid_mz, sizeof(int)*9 );
}

/*
 ConstructNi()
 + Defines the basis function.
 + _1D implies the spatial degress of freedom expected in the _xi[] array
*/
void ConstructNi_Q2_1D( const QPoint1d *q, double Ni[] )
{
	double xi = q->xi;

	Ni[0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
	Ni[1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
	Ni[2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )
}

void ConstructGNi_Q2_1D( const QPoint1d *q, double **GNi )
{
	double xi = q->xi;

	GNi[0][0] = 0.5 * ( 2.0*xi - 1.0 );
	GNi[0][1] = - 2.0*xi;
	GNi[0][2] = 0.5 * ( 2.0*xi + 1.0 );
}

void ConstructNi_Q2_2D( const QPoint2d *q, double Ni[] )
{
	int i,j,d,cnt;
	double basis_NI[2][3],_xi[2];

	_xi[0] = q->xi;
	_xi[1] = q->eta;
	for( d=0; d<2; d++ ) {
		double xi = _xi[d];

		basis_NI[d][0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )
	}

	cnt = 0;
	for( j=0; j<3; j++ ) {
		for( i=0; i<3; i++ ) {
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
void ConstructGNi_Q2_2D( const QPoint2d *q, double **GNi )
{
	double basis_NI[2][3];
	double basis_GNI[2][3],_xi[2];
	int i,j,d,cnt;


	_xi[0] = q->xi;
	_xi[1] = q->eta;
	for( d=0; d<2; d++ ) {
		double xi = _xi[d];

		basis_NI[d][0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )

		basis_GNI[d][0] = 0.5 * ( 2.0*xi - 1.0 );
		basis_GNI[d][1] = - 2.0*xi;
		basis_GNI[d][2] = 0.5 * ( 2.0*xi + 1.0 );
	}

	cnt = 0;
	for( j=0; j<3; j++ ) {
		for( i=0; i<3; i++ ) {

			GNi[0][cnt] = basis_GNI[0][i]  *  basis_NI[1][j];
			GNi[1][cnt] =  basis_NI[0][i]  * basis_GNI[1][j];

			cnt++;
		}
	}
}

void ConstructNi_Q2_3D( const QPoint3d *q, double Ni[] )
{
	int i,j,k,d,cnt;
	double basis_NI[3][3],_xi[3];

	_xi[0] = q->xi;
	_xi[1] = q->eta;
	_xi[2] = q->zeta;
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

void ConstructGNi_Q2_3D( const QPoint3d *q, double **GNi )
{
	double basis_NI[3][3];
	double basis_GNI[3][3],_xi[3];
	int i,j,k,d,cnt;


	_xi[0] = q->xi;
	_xi[1] = q->eta;
	_xi[2] = q->zeta;
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

/*
 extract_point_field()
 + Give an element field, the value assocated with point_id will be returned in point_e_data[]
*/
void extract_point_field_Q2_1D( ConformingElementFamily e, int point_id, int dofs, double e_data[], double point_e_data[] )
{
	int n;
	int **pid_node_list;

	pid_node_list = e->point_node_list;

	for( n=0; n<e->n_nodes_0D; n++ ) {
		memcpy( &point_e_data[dofs*n] , &e_data[ dofs*pid_node_list[point_id][n] ], sizeof(double)*dofs );
	}
}

/*
 extract_point_field()
 + Give an element field, the value assocated with point_id will be returned in point_e_data[]
 */
void extract_edge_field_Q2_2D( ConformingElementFamily e, QuadElementEdge edge_id, int dofs, double e_data[], double edge_e_data[] )
{
	int n;
	int **eid_node_list;

	eid_node_list = e->edge_node_list;

	for( n=0; n<e->n_nodes_1D; n++ ) {
		memcpy( &edge_e_data[dofs*n] , &e_data[ dofs*eid_node_list[edge_id][n] ], sizeof(double)*dofs );
	}
}

void extract_surface_field_Q2_3D( ConformingElementFamily e, HexElementFace face_id, int dofs, double e_data[], double surf_e_data[] )
{
	int n;
	int **fid_node_list;

	fid_node_list = e->face_node_list;

	for( n=0; n<e->n_nodes_2D; n++ ) {
		memcpy( &surf_e_data[dofs*n] , &e_data[ dofs*fid_node_list[face_id][n] ], sizeof(double)*dofs );
	}
}

void extract_element_surface_field_Q2_3D( int face_id, int fid_node_list[][9], int dofs, double e_data[], double surf_e_data[] )
{
	int n;

	for( n=0; n<9; n++ ) {
		memcpy( &surf_e_data[dofs*n] , &e_data[ dofs*fid_node_list[face_id][n] ], sizeof(double)*dofs );
	}
}

/*
 gp_w[]    : the quad weights in 1-space (s integration space). They will sum to 2. Length 1 * 3.
 gp_xi_1[] : the local coordinates of the quad points in 1-space (s integration space). Length 1 * 3
 gp_xi_2[] : the local coordinates of the quad points in 2-space (xi-eta-zeta isoparametric space). Length 2 * 3
*/
void CLASS_get_face_quadrature_points_Q2_2D( ConformingElementFamily e, QuadElementEdge face_id, int *ngp_s, QPoint1d gp1[], QPoint2d gp2[] )
{
	int I;
	const double sqrt_15_on_5 = 0.774596669241483; /* sqrt(15)/5 */
	const double five_on_9 = 0.555555555555556;
	const double eight_on_9 = 0.888888888888889;
	const double w_1d[] = { five_on_9, eight_on_9, five_on_9 };
	const double xi_1d[] = { -sqrt_15_on_5, 0.0, sqrt_15_on_5 };


	/* standard 3 point quadrature */
	*ngp_s = 3;

	for( I=0; I<3; I++ ) {
		gp1[I].w  = w_1d[I];
		gp1[I].xi = xi_1d[I];

		gp2[I].w = gp1[I].w;
	}


	switch (face_id) {

		case QUAD_EDGE_Nxi:
			for( I=0; I<3; I++ ) {
				gp2[I].xi  = -1.0;
				gp2[I].eta = gp1[I].xi;
			}
			break;

		case QUAD_EDGE_Pxi:
			for( I=0; I<3; I++ ) {
				gp2[I].xi  = 1.0;
				gp2[I].eta = gp1[I].xi;
			}
			break;

		case QUAD_EDGE_Neta:
			for( I=0; I<3; I++ ) {
				gp2[I].xi  = gp1[I].xi;
				gp2[I].eta = -1.0;
			}
			break;

		case QUAD_EDGE_Peta:
			for( I=0; I<3; I++ ) {
				gp2[I].xi  = gp1[I].xi;
				gp2[I].eta = 1.0;
			}
			break;

		default:
			break;
	}
}

/*
 gp_w[]    : the quad weights in 2-space (s-t integration space). They will sum to 4. Length 2 * 9.
 gp_xi_2[] : the local coordinates of the quad points in 2-space (s-t integration space). Length 2 * 9
 gp_xi_3[] : the local coordinates of the quad points in 3-space (xi-eta-zeta isoparametric space). Length 3 * 9
 */
void CLASS_get_face_quadrature_points_Q2_3D( ConformingElementFamily e, HexElementFace face_id, int *ngp_s, QPoint2d gp2[], QPoint3d gp3[] )
{
	const double sqrt_15_on_5 = 0.774596669241483; /* sqrt(15)/5 */
	const double five_on_9 = 0.555555555555556;
	const double eight_on_9 = 0.888888888888889;
	const double w_1d[] = { five_on_9, eight_on_9, five_on_9 };
	const double xi_1d[] = { -sqrt_15_on_5, 0.0, sqrt_15_on_5 };
	double local_coords[] = {	-1.0,-1.0,-1.0,		0.0,-1.0,-1.0,		1.0,-1.0,-1.0,
		-1.0, 0.0,-1.0,		0.0, 0.0,-1.0,		1.0, 0.0,-1.0,
		-1.0, 1.0,-1.0,		0.0, 1.0,-1.0,		1.0, 1.0,-1.0,
		//
		-1.0,-1.0, 0.0,		0.0,-1.0, 0.0,		1.0,-1.0, 0.0,
		-1.0, 0.0, 0.0,		0.0, 0.0, 0.0,		1.0, 0.0, 0.0,
		-1.0, 1.0, 0.0,		0.0, 1.0, 0.0,		1.0, 1.0, 0.0,
		//
		-1.0,-1.0, 1.0,		0.0,-1.0, 1.0,		1.0,-1.0, 1.0,
		-1.0, 0.0, 1.0,		0.0, 0.0, 1.0,		1.0, 0.0, 1.0,
		-1.0, 1.0, 1.0,		0.0, 1.0, 1.0,		1.0, 1.0, 1.0 };
	double face_local_coord[9*3];
	double Ni_st[9];
	int I,J;


	e->extract_surface_field( e, face_id, 3, local_coords, face_local_coord );

	/* standard 3x3 point quadrature */
	*ngp_s = 9;

	for( I=0; I<3; I++ ) {
		for( J=0; J<3; J++ ) {
			int idx = I + J * 3;

			gp2[idx].w   = w_1d[I] * w_1d[J];
			gp2[idx].xi  = xi_1d[I];
			gp2[idx].eta = xi_1d[J];

			gp3[idx].w = gp2[idx].w;
		}
	}

	/* interpolate the quadrature points from s-t to xi-eta-zeta */
	for( I=0; I<3; I++ ) {
		for( J=0; J<3; J++ ) {
			double x_s,y_s,z_s;
			int k;
			int idx = I + J * 3;

			ConstructNi_Q2_2D( &gp2[idx], Ni_st );

			x_s = y_s = z_s = 0.0;
			for( k=0; k<9; k++ ) {
				x_s = x_s + Ni_st[k] * face_local_coord[3*k+0];
				y_s = y_s + Ni_st[k] * face_local_coord[3*k+1];
				z_s = z_s + Ni_st[k] * face_local_coord[3*k+2];
			}
			gp3[idx].xi = x_s;
			gp3[idx].eta = y_s;
			gp3[idx].zeta = z_s;
		}
	}
}

/*
 get_volume_quadrature_points()
 + Defines the weights and quadrature points for volume integration.
*/
void get_volume_quadrature_points_Q2_2D( int *ngp, QPoint2d gp2[] )
{
	const double sqrt_15_on_5 = 0.774596669241483; /* sqrt(15)/5 */
	const double five_on_9 = 0.555555555555556;
	const double eight_on_9 = 0.888888888888889;
	const double w_1d[] = { five_on_9, eight_on_9, five_on_9 };
	const double xi_1d[] = { -sqrt_15_on_5, 0.0, sqrt_15_on_5 };
	int I,J;


	/* standard 3x3 point quadrature */
	*ngp = 9;

	for( I=0; I<3; I++ ) {
		for( J=0; J<3; J++ ) {
			int idx = I + J*3;

			gp2[idx].w = w_1d[I] * w_1d[J];

			gp2[idx].xi = xi_1d[I];
			gp2[idx].eta = xi_1d[J];
		}
	}
}

void get_volume_quadrature_points_Q2_3D( int *ngp, QPoint3d gp3[] )
{
	const double sqrt_15_on_5 = 0.774596669241483; /* sqrt(15)/5 */
	const double five_on_9 = 0.555555555555556;
	const double eight_on_9 = 0.888888888888889;
	const double w_1d[] = { five_on_9, eight_on_9, five_on_9 };
	const double xi_1d[] = { -sqrt_15_on_5, 0.0, sqrt_15_on_5 };
	int I,J,K;


	/* standard 3x3x3 point quadrature */
	*ngp = 27;

	for( I=0; I<3; I++ ) {
		for( J=0; J<3; J++ ) {
			for( K=0; K<3; K++ ) {
				int idx = I + J*3 + K*3*3;

				gp3[idx].w = w_1d[I] * w_1d[J] * w_1d[K];

				gp3[idx].xi = xi_1d[I];
				gp3[idx].eta = xi_1d[J];
				gp3[idx].zeta = xi_1d[K];
			}
		}
	}
}

/*
 compute_surface_geometry()
 + Defines the Jacobian for surface integration.
 */
void compute_surface_geometry_Q2_2D(	ConformingElementFamily e,
										const double coords[], // should contain 4 points with dimension 2 (x,y) //
										QuadElementEdge edge_idx,       // face node ids (size = 2) //
										const QPoint1d *xi0,    // should contain 1 point with dimension 1 (xi)   //
										double n0[],double t0[], double *J )  // n0[] contains 1 point with dimension 2 (x,y) //
{
	double x_s, y_s;
	double mag;
	double __GNi_1d[3];
	double *GNi_1d[1];
	const int *fid = e->edge_node_list[edge_idx];
	int i;

	GNi_1d[0] = &__GNi_1d[0];

	// x_,s = dN_i/ds x_i //
	ConstructGNi_Q2_1D( xi0, GNi_1d );
//	x_s = GNi_1d[0][0] * coords[2*fid[0]  ] + GNi_1d[0][1] * coords[2*fid[1]  ] + GNi_1d[0][2] * coords[2*fid[2]  ];
//	y_s = GNi_1d[0][0] * coords[2*fid[0]+1] + GNi_1d[0][1] * coords[2*fid[1]+1] + GNi_1d[0][2] * coords[2*fid[2]+1];

	x_s = y_s = 0.0;
	for (i=0; i<3; i++) {
		x_s += GNi_1d[0][i] * coords[2*fid[i]  ];
		y_s += GNi_1d[0][i] * coords[2*fid[i]+1];
	}

	mag = sqrt( x_s*x_s + y_s*y_s );
	if ( (edge_idx==QUAD_EDGE_Nxi) || (edge_idx==QUAD_EDGE_Peta) ) {
		n0[0] = -y_s/mag;
		n0[1] =  x_s/mag;

		t0[0] = x_s/mag;
		t0[1] = y_s/mag;
	} else { /* reverse P_xi, N_eta */
		n0[0] =  y_s/mag;
		n0[1] = -x_s/mag;

		t0[0] = -x_s/mag;
		t0[1] = -y_s/mag;
	}

	*J = mag;
}

void compute_surface_geometry_Q2_3D(	ConformingElementFamily e,
										const double coords[],   // should contain 27 points with dimension 3 (x,y,z) //
										HexElementFace face_idx,         // face node ids (size = 9) //
										const QPoint2d *xi0,      // should contain 1 point with dimension 2 (xi,eta)   //
										double n0[], double t0[],double *J ) // n0[],t0[] contains 1 point with dimension 3 (x,y,z) //
{
	enum {SURFACE_DIM=2};
	enum {Q2_NPS=9};
	double x_s,y_s,z_s;
	double x_t,y_t,z_t;
	double n[3],t[3],mag,magt;
	int i;
	double __GNi_st[SURFACE_DIM*Q2_NPS];
	double *GNi_st[SURFACE_DIM];
	const int *fid = e->face_node_list[(int)face_idx];


	GNi_st[0] = &__GNi_st[0];
	GNi_st[1] = &__GNi_st[Q2_NPS];

/*
	 x_s \cross x_t  =	|  i   j   k  |
											| x_s y_s z_s |
						          | x_t y_t z_t |

	 n_i =  (y_s*z_t - z_s*y_t)
	 n_j = -(x_s*z_t - z_s*x_t)
	 n_k =  (x_s*y_t - y_s*x_t)
*/

	ConstructGNi_Q2_2D( xi0, GNi_st );

	// x_,s = dN_i/ds x_i //
	x_s = y_s = z_s = 0.0;
	x_t = y_t = z_t = 0.0;
	for( i=0; i<Q2_NPS; i++ ) {
		x_s = x_s + GNi_st[0][i] * coords[3*fid[i]  ];
		y_s = y_s + GNi_st[0][i] * coords[3*fid[i]+1];
		z_s = z_s + GNi_st[0][i] * coords[3*fid[i]+2];

		x_t = x_t + GNi_st[1][i] * coords[3*fid[i]  ];
		y_t = y_t + GNi_st[1][i] * coords[3*fid[i]+1];
		z_t = z_t + GNi_st[1][i] * coords[3*fid[i]+2];
	}
	n[0] =  (y_s*z_t - z_s*y_t);
	n[1] = -(x_s*z_t - z_s*x_t);
	n[2] =  (x_s*y_t - y_s*x_t);

	mag = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

	if (n0) {
		n0[0] = n[0]/mag;
		n0[1] = n[1]/mag;
		n0[2] = n[2]/mag;
	}

	/*
	 There are two tangent vectors, we return only the one in the "s" direction.

	 Tangents are constructured such that
	 t1 x t2 = normal
	 normal x t1 = t2
	 t2 x norma = t1
	*/
	if (t0) {
		t[0] = x_s;
		t[1] = y_s;
		t[2] = z_s;
		magt = sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
		t0[0] = t[0]/magt;
		t0[1] = t[1]/magt;
		t0[2] = t[2]/magt;
	}

	if (J) {
		*J = mag;
	}
}

void ElementHelper_matrix_inverse_2x2(double A[2][2],double B[2][2])
{
	double det;

	det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	det = 1.0/det;

	B[0][0] = det * A[1][1];
	B[0][1] = -det * A[0][1];
	B[1][0] = -det * A[1][0];
	B[1][1] = det * A[0][0];
}

void ElementHelper_matrix_inverse_3x3(double A[3][3],double B[3][3])
{
	double t4, t6, t8, t10, t12, t14, t17;

	t4 = A[2][0] * A[0][1];
	t6 = A[2][0] * A[0][2];
	t8 = A[1][0] * A[0][1];
	t10 = A[1][0] * A[0][2];
	t12 = A[0][0] * A[1][1];
	t14 = A[0][0] * A[1][2];
	t17 = 0.1e1 / (t4 * A[1][2] - t6 * A[1][1] - t8 * A[2][2] + t10 * A[2][1] + t12 * A[2][2] - t14 * A[2][1]);

	B[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * t17;
	B[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) * t17;
	B[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * t17;
	B[1][0] = -(-A[2][0] * A[1][2] + A[1][0] * A[2][2]) * t17;
	B[1][1] = (-t6 + A[0][0] * A[2][2]) * t17;
	B[1][2] = -(-t10 + t14) * t17;
	B[2][0] = (-A[2][0] * A[1][1] + A[1][0] * A[2][1]) * t17;
	B[2][1] = -(-t4 + A[0][0] * A[2][1]) * t17;
	B[2][2] = (-t8 + t12) * t17;
}

/*
 compute_volume_geometry()
 + Defines the Jacobian for volume integration.
 */
void compute_volume_geometry_Q2_3D( const double coords[], const double **GNi, double **GNx, double *det_J )
{
	const int npe = 27;
	const int nsd = 3;
	double detJ;
	double iJ[3][3],J[3][3];
	int i;
	double a,b,c,d,e,f,g,h,ii;


	J[0][0] = J[0][1] = J[0][2] = 0.0;
	J[1][0] = J[1][1] = J[1][2] = 0.0;
	J[2][0] = J[2][1] = J[2][2] = 0.0;
	for( i=0; i<npe; i++ ) {
		double cx = coords[ nsd*i + 0 ];
		double cy = coords[ nsd*i + 1 ];
		double cz = coords[ nsd*i + 2 ];
		// J_ij = d(x_j) / d(xi_i)  // J_ij = \sum _I GNi[j][I} * x_i
		J[0][0] += GNi[0][i] * cx;  // J_xx
		J[0][1] += GNi[0][i] * cy;  // J_xy = dx/deta
		J[0][2] += GNi[0][i] * cz;  // J_xz = dx/dzeta

		J[1][0] += GNi[1][i] * cx;  // J_yx = dy/dxi
		J[1][1] += GNi[1][i] * cy;  // J_yy
		J[1][2] += GNi[1][i] * cz;  // J_yz

		J[2][0] += GNi[2][i] * cx;  // J_zx
		J[2][1] += GNi[2][i] * cy;  // J_zy
		J[2][2] += GNi[2][i] * cz;  // J_zz
	}

	/* invert */
	/* Compute inverse and determinant (derived with MAPLE) */
	a = J[0][0]; b = J[0][1]; c = J[0][2];
	d = J[1][0]; e = J[1][1]; f = J[1][2];
	g = J[2][0]; h = J[2][1]; ii = J[2][2]; // nasty i => ii bug... //

	ElementHelper_matrix_inverse_3x3(J,iJ);

	detJ = a*e*ii - a*f*h - d*b*ii + d*c*h + g*b*f - g*c*e;

	if( detJ<1.0e-20 ) {
		printf( "***** element jacobian appears to be zero or negative ***** \n");
		printf( "*****           |J| = %1.6e  \n", detJ );
	}

	for( i=0; i<npe; i++ ) {
		// GNx[i][I] = \sum_j iJ[i][j] * GNi[j][I] // correct //
		GNx[0][i] = iJ[0][0]*GNi[0][i] + iJ[0][1]*GNi[1][i] + iJ[0][2]*GNi[2][i];
		GNx[1][i] = iJ[1][0]*GNi[0][i] + iJ[1][1]*GNi[1][i] + iJ[1][2]*GNi[2][i];
		GNx[2][i] = iJ[2][0]*GNi[0][i] + iJ[2][1]*GNi[1][i] + iJ[2][2]*GNi[2][i];
	}

	*det_J = detJ;
}

void compute_volume_geometry_Q2_2D( const double coords[], const double **GNi, double **GNx, double *det_J )
{
	const int npe = 9;
	const int nsd = 2;
	double J00,J01, J10,J11, detJ;
	double iJ00,iJ01, iJ10,iJ11;
	int i;


	J00 = J01 = 0.0;
	J10 = J11 = 0.0;
	for( i=0; i<npe; i++ ) {
		double cx = coords[ nsd*i + 0 ];
		double cy = coords[ nsd*i + 1 ];

		// J_ij = d(x_j) / d(xi_i)  // J_ij = \sum _I GNi[j][I} * x_i
		J00 = J00 + GNi[0][i] * cx;  // J_xx
		J01 = J01 + GNi[0][i] * cy;  // J_xy = dx/deta

		J10 = J10 + GNi[1][i] * cx;  // J_yx = dy/dxi
		J11 = J11 + GNi[1][i] * cy;  // J_yy
	}

	detJ = (J00*J11) - (J01*J10);

	iJ00 =  J11/detJ;
	iJ01 = -J01/detJ;
	iJ10 = -J10/detJ;
	iJ11 =  J00/detJ;

	for( i=0; i<npe; i++ ) {
		GNx[0][i] = GNi[0][i] * iJ00 + GNi[1][i] * iJ01;
		GNx[1][i] = GNi[0][i] * iJ10 + GNi[1][i] * iJ11;
	}

	if( detJ<1.0e-20 ) {
		printf( "***** element jacobian appears to be zero or negative ***** \n");
		printf( "*****           |J| = %1.6e  \n", detJ );
	}


	*det_J = detJ;
}

void compute_element_normal_Q2_2D(	ConformingElementFamily e,
																	const double coords[],   // should contain 9 points with dimension 2 (x,y) //
																	QuadElementEdge face_idx,  // face id //
																	const QPoint1d *xi0,      // should contain 1 point with dimension 1 (xi,)   //
																	double n0[]) // n0[] contains 1 point with dimension 2 (x,y) //
{
	enum {SURFACE_DIM=1};
	enum {Q2_NPS = 3};
	double x_s,y_s;
	double n[2],mag;
	int i;
	double __GNi_st[SURFACE_DIM*Q2_NPS];
	double *GNi_st[SURFACE_DIM];
	const int *fid = e->edge_node_list[(int)face_idx];


	GNi_st[0] = &__GNi_st[0];

	ConstructGNi_Q2_1D( xi0, GNi_st );

	// x_,s = dN_i/ds x_i //
	x_s = y_s = 0.0;
	for( i=0; i<Q2_NPS; i++ ) {
		x_s = x_s + GNi_st[0][i] * coords[2*fid[i]  ];
		y_s = y_s + GNi_st[0][i] * coords[2*fid[i]+1];
	}
	n[0] = -y_s;
	n[1] =  x_s;

	mag = sqrt( n[0]*n[0] + n[1]*n[1] );

	n0[0] = n[0]/mag;
	n0[1] = n[1]/mag;
}

void compute_element_normal_Q2_3D(	ConformingElementFamily e,
									const double coords[],   // should contain 27 points with dimension 3 (x,y,z) //
									HexElementFace face_idx,  // face id //
									const QPoint2d *xi0,      // should contain 1 point with dimension 2 (xi,eta)   //
									double n0[]) // n0[] contains 1 point with dimension 3 (x,y,z) //
{
	enum {SURFACE_DIM=2};
	enum {Q2_NPS=9};
	double x_s,y_s,z_s;
	double x_t,y_t,z_t;
	double n[3],mag;
	int i;
	double __GNi_st[SURFACE_DIM*Q2_NPS];
	double *GNi_st[SURFACE_DIM];
	const int *fid = e->face_node_list[(int)face_idx];


	GNi_st[0] = &__GNi_st[0];
	GNi_st[1] = &__GNi_st[Q2_NPS];

/*
   x_s \cross x_t  = |  i   j   k  |
										 | x_s y_s z_s |
										 | x_t y_t z_t |

    n_i =  (y_s*z_t - z_s*y_t)
    n_j = -(x_s*z_t - z_s*x_t)
    n_k =  (x_s*y_t - y_s*x_t)
*/

	ConstructGNi_Q2_2D( xi0, GNi_st );

	// x_,s = dN_i/ds x_i //
	x_s = y_s = z_s = 0.0;
	x_t = y_t = z_t = 0.0;
	for( i=0; i<Q2_NPS; i++ ) {
		x_s = x_s + GNi_st[0][i] * coords[3*fid[i]  ];
		y_s = y_s + GNi_st[0][i] * coords[3*fid[i]+1];
		z_s = z_s + GNi_st[0][i] * coords[3*fid[i]+2];

		x_t = x_t + GNi_st[1][i] * coords[3*fid[i]  ];
		y_t = y_t + GNi_st[1][i] * coords[3*fid[i]+1];
		z_t = z_t + GNi_st[1][i] * coords[3*fid[i]+2];
	}
	n[0] =  (y_s*z_t - z_s*y_t);
	n[1] = -(x_s*z_t - z_s*x_t);
	n[2] =  (x_s*y_t - y_s*x_t);

	mag = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

	n0[0] = n[0]/mag;
	n0[1] = n[1]/mag;
	n0[2] = n[2]/mag;

}

/*
 compute_surface_geometry()
 + Defines the Jacobian for surface integration.
 */
void compute_element_tangent_Q2_2D(	ConformingElementFamily e,
																		const double coords[], // should contain 4 points with dimension 2 (x,y) //
																		QuadElementEdge edge_idx,       // face node ids (size = 2) //
																		const QPoint1d *xi0,    // should contain 1 point with dimension 1 (xi)   //
																		double t0[])  // t0[] contains 1 point with dimension 2 (x,y) //
{
	double x_s, y_s;
	double mag;
	double __GNi_1d[3];
	double *GNi_1d[1];
	const int *fid = e->edge_node_list[edge_idx];
	int i;

	GNi_1d[0] = &__GNi_1d[0];

	// x_,s = dN_i/ds x_i //
	ConstructGNi_Q2_1D( xi0, GNi_1d );

	x_s = y_s = 0.0;
	for (i=0; i<3; i++) {
		x_s += GNi_1d[0][i] * coords[2*fid[i]  ];
		y_s += GNi_1d[0][i] * coords[2*fid[i]+1];
	}

	mag = sqrt( x_s*x_s + y_s*y_s );
	if ( (edge_idx==QUAD_EDGE_Nxi) || (edge_idx==QUAD_EDGE_Peta) ) {

		t0[0] = x_s/mag;
		t0[1] = y_s/mag;
	} else { /* reverse P_xi, N_eta */
		t0[0] = -x_s/mag;
		t0[1] = -y_s/mag;
	}

}

void compute_element_tangent_Q2_3D(	ConformingElementFamily e,
																	const double coords[],   // should contain 27 points with dimension 3 (x,y,z) //
																	HexElementFace face_idx,  // face id //
																	const QPoint2d *xi0,      // should contain 1 point with dimension 2 (xi,eta)   //
																	double t1[],double t2[] ) // t1[],t2[] contains 1 point with dimension 3 (x,y,z) //
{
	enum {SURFACE_DIM=2};
	enum {Q2_NPS=9};
	double x_s,y_s,z_s;
	double x_t,y_t,z_t;
	double mag;
	int i;
	double __GNi_st[SURFACE_DIM*Q2_NPS];
	double *GNi_st[SURFACE_DIM];
	const int *fid = e->face_node_list[(int)face_idx];


	GNi_st[0] = &__GNi_st[0];
	GNi_st[1] = &__GNi_st[Q2_NPS];

	ConstructGNi_Q2_2D( xi0, GNi_st );

	// x_,s = dN_i/ds x_i //
	x_s = y_s = z_s = 0.0;
	x_t = y_t = z_t = 0.0;
	for( i=0; i<Q2_NPS; i++ ) {
		x_s = x_s + GNi_st[0][i] * coords[3*fid[i]  ];
		y_s = y_s + GNi_st[0][i] * coords[3*fid[i]+1];
		z_s = z_s + GNi_st[0][i] * coords[3*fid[i]+2];

		x_t = x_t + GNi_st[1][i] * coords[3*fid[i]  ];
		y_t = y_t + GNi_st[1][i] * coords[3*fid[i]+1];
		z_t = z_t + GNi_st[1][i] * coords[3*fid[i]+2];
	}

	t1[0] = x_s;  t1[1] = y_s;  t1[2] = z_s;
	t2[0] = x_t;  t2[1] = y_t;  t2[2] = z_t;

	mag = sqrt( t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2] );
	t1[0] = t1[0]/mag;  t1[1] = t1[1]/mag;  t1[2] = t1[2]/mag;

	mag = sqrt( t2[0]*t2[0] + t2[1]*t2[1] + t2[2]*t2[2] );
	t2[0] = t2[0]/mag;  t2[1] = t2[1]/mag;  t2[2] = t2[2]/mag;
}

void ElementTypeDestroy_Q2(ConformingElementFamily *_e)
{
	ConformingElementFamily e;

	if (!_e) { return; }

	e = *_e;
	free(e->name);

	free2d(&e->point_node_list);
	free2d(&e->edge_node_list);
	free2d(&e->face_node_list);

	free(e);
	*_e = NULL;
}

void ElementTypeCreate_Q2(ConformingElementFamily *_e,const int dim)
{
	ConformingElementFamily e;

	e = calloc( 1, sizeof(struct _p_ConformingElementFamily) );

	if (asprintf( &e->name, "Q2(%d)", dim ) < 0) {printf("asprintf() error. Exiting ungracefully.\n"); exit(1);}
	e->nsd = dim;

	e->n_nodes_0D = 1;
	e->n_nodes_1D = 3;
	e->n_nodes_2D = 9;
	e->n_nodes_3D = 27;

	e->basis_NI_1D = ConstructNi_Q2_1D;
	e->basis_NI_2D = ConstructNi_Q2_2D;
	e->basis_NI_3D = ConstructNi_Q2_3D;

	e->basis_GNI_1D = ConstructGNi_Q2_1D;
	e->basis_GNI_2D = ConstructGNi_Q2_2D;
	e->basis_GNI_3D = ConstructGNi_Q2_3D;

	e->n_points = 2;
	e->generate_all_point_ids = generate_all_face_ids_Q2_1D;
	e->extract_point_field    = extract_point_field_Q2_1D;

	malloc2d( e->n_points, e->n_nodes_0D, &e->point_node_list );
	e->generate_all_point_ids( e->point_node_list );

	e->n_edges = 4;
	e->generate_all_edge_ids = generate_all_face_ids_Q2_2D;
	e->extract_edge_field    = extract_edge_field_Q2_2D;

	malloc2d( e->n_edges, e->n_nodes_1D, &e->edge_node_list );
	e->generate_all_edge_ids( e->edge_node_list );


	e->n_faces = 6;
	e->generate_all_face_ids = generate_all_face_ids_Q2_3D;
	e->extract_surface_field = extract_surface_field_Q2_3D; // Element,face,dofs,e_data[],surface_e_data[]

	malloc2d( e->n_faces, e->n_nodes_2D, &e->face_node_list );
	e->generate_all_face_ids( e->face_node_list );

	/* quadrature routines */
	e->generate_surface_quadrature_2D = CLASS_get_face_quadrature_points_Q2_2D;
	e->generate_surface_quadrature_3D = CLASS_get_face_quadrature_points_Q2_3D;

	e->generate_volume_quadrature_2D = get_volume_quadrature_points_Q2_2D;
	e->generate_volume_quadrature_3D = get_volume_quadrature_points_Q2_3D;

	/* integration - geometry routines */
	e->compute_surface_geometry_2D = compute_surface_geometry_Q2_2D;
	e->compute_surface_geometry_3D = compute_surface_geometry_Q2_3D;

	e->compute_volume_geometry_2D = compute_volume_geometry_Q2_2D;
	e->compute_volume_geometry_3D = compute_volume_geometry_Q2_3D;

	e->compute_surface_normal_2D   = compute_element_normal_Q2_2D;
	e->compute_surface_normal_3D   = compute_element_normal_Q2_3D;

	e->compute_surface_tangent_2D  = compute_element_tangent_Q2_2D;
	e->compute_surface_tangents_3D = compute_element_tangent_Q2_3D;

	*_e = e;
}

/* testing code */

void eval_cross_product(double x1[],double x2[],double r[])
{
	double x_s,y_s,z_s;
	double x_t,y_t,z_t;

	x_s = x1[0]; y_s = x1[1]; z_s = x1[2];
	x_t = x2[0]; y_t = x2[1]; z_t = x2[2];

	r[0] =  (y_s*z_t - z_s*y_t);
	r[1] = -(x_s*z_t - z_s*x_t);
	r[2] =  (x_s*y_t - y_s*x_t);
}

void insert_point_3D( const int i, const double x, const double y, const double z, double coords[] )
{
	coords[ 3*i    ] = x;
	coords[ 3*i +1 ] = y;
	coords[ 3*i +2 ] = z;
}

void insert_point_2D( const int i, const double x, const double y, double coords[] )
{
	coords[ 2*i    ] = x;
	coords[ 2*i +1 ] = y;
}

void testEQ2_SurfQuad3d_1( void )
{
	ConformingElementFamily e;
	int i,j,k;
	int cnt;
	double coords[27*3];
	double x,y,z;
	int ngp_s;
	QPoint2d gp2[9];
	QPoint3d gp3[9];
	int f,p;
	double sa,exact_sa;
	double dx, dy, dz;
	double xi,yi,zi,Lx,Ly,Lz,xf,yf,zf;


	printf("TEST(%s)\n", __func__);
	ElementTypeCreate_Q2(&e,3);

	xi = -10.0;		Lx = 24.0;
	yi = 1.0;		Ly = 4.0;
	zi = 3.0;		Lz = 122.2;

	xf = xi + Lx;
	yf = yi + Ly;
	zf = zi + Lz;

	dx = (xf-xi)/2.0;
	dy = (yf-yi)/2.0;
	dz = (zf-zi)/2.0;

	cnt = 0;
	for( k=0; k<3; k++ ) {
		for( j=0; j<3; j++ ) {
			for( i=0; i<3; i++ ) {
				x = xi + i * dx;
				y = yi + j * dy;
				z = zi + k * dz;

				insert_point_3D( cnt, x,y,z, coords );

				cnt++;
			}
		}
	}

	sa = 0.0;
	for( f=0; f<HEX_EDGES; f++ ) {
		e->generate_surface_quadrature_3D( e, f, &ngp_s, gp2, gp3 );

		for( p=0; p<ngp_s; p++ ) {
			double n0[3],t0[3],t1[3],t2[3],J;


			e->compute_surface_geometry_3D(	e, coords,	// should contain 27 points with dimension 3 (x,y,z) //
																		 f,			// face node ids (size = 9) //
																		 &gp2[p],    // should contain 1 point with dimension 2 (s,t)   //
																		 n0,t0, &J );	// n0[],t0 contains 1 point with dimension 3 (x,y,z) //

			compute_element_normal_Q2_3D(	e, coords,	// should contain 27 points with dimension 3 (x,y,z) //
											f,			// face node ids (size = 9) //
											&gp2[p],    // should contain 1 point with dimension 2 (xi,eta)   //
											n0 );	// n0[] contains 1 point with dimension 3 (x,y,z) //

			compute_element_tangent_Q2_3D(	e, coords,	// should contain 27 points with dimension 3 (x,y,z) //
																	 f,			// face node ids (size = 9) //
																	 &gp2[p],    // should contain 1 point with dimension 2 (xi,eta)   //
																	 t1,t2 );	// t1[],t2[] contains 1 point with dimension 3 (x,y,z) //

			if(p==0) {
				printf("f=%.2d: normal    = %+1.6lf, %+1.6lf, %+1.6lf \n", f, n0[0],n0[1],n0[2] );
				printf("f=%.2d: tangent1  = %+1.6lf, %+1.6lf, %+1.6lf \n", f, t1[0],t1[1],t1[2] );
				printf("f=%.2d: tangent2  = %+1.6lf, %+1.6lf, %+1.6lf \n", f, t2[0],t2[1],t2[2] );
			}

			sa = sa + gp2[p].w * 1.0 * J;
		}
		printf("\n");
	}
	exact_sa = 2.0*Lx*Ly + 2.0*Ly*Lz + 2.0*Lz*Lx;
	printf("sa = %+1.6f : exact = %+1.6f \n", sa, exact_sa );
}

void testEQ2_SurfQuad2d_1( void )
{
	ConformingElementFamily e;
	int i,j;
	int cnt;
	double coords[9*2];
	double x,y;
	int ngp_s;
	QPoint1d gp1[3];
	QPoint2d gp2[3];
	int f,p;
	double sa,exact_sa;
	double dx, dy;
	double xi,yi,Lx,Ly,xf,yf;


	printf("TEST(%s)\n", __func__);

	ElementTypeCreate_Q2(&e,2);

	xi = -10.0; Lx = 24.0;
	yi = 1.0;		Ly = 4.0;

	xf = xi + Lx;
	yf = yi + Ly;

	dx = (xf-xi)/2.0;
	dy = (yf-yi)/2.0;

	cnt = 0;
	for( j=0; j<3; j++ ) {
		for( i=0; i<3; i++ ) {
			x = xi + i * dx;
			y = yi + j * dy;

			insert_point_2D( cnt, x,y, coords );

			cnt++;
		}
	}

	sa = 0.0;
	for( f=0; f<QUAD_EDGES; f++ ) {
		e->generate_surface_quadrature_2D( e, f, &ngp_s, gp1, gp2 );

		for( p=0; p<ngp_s; p++ ) {
			double n0[2], t0[2], J;

			e->compute_surface_geometry_2D(	e, coords,	// should contain 27 points with dimension 3 (x,y,z) //
																		 f,			// face node ids (size = 9) //
																		 &gp1[p],    // should contain 1 point with dimension 1 (xi)   //
																		 n0,t0, &J );	// n0[],t0[] contains 1 point with dimension 2 (x,y) //
			if(p==0) {
				printf("f=%.2d: normal  = %+1.6lf, %+1.6lf \n", f, n0[0],n0[1] );
				printf("f=%.2d: tangent = %+1.6lf, %+1.6lf \n", f, t0[0],t0[1] );
			}

			sa = sa + gp1[p].w * 1.0 * J;
		}
	}
	exact_sa = 2.0*Lx + 2.0*Ly;
	printf("sa = %1.6f : exact = %1.6f \n", sa, exact_sa );
}

void testEQ2_FaceIds(void)
{
	ConformingElementFamily e;
	int f,i;

	printf("TEST(%s)\n", __func__);

	ElementTypeCreate_Q2(&e,3);

	printf("LINE \n");
	for (f=0; f<LINE_EDGES; f++) {
		for (i=0; i<e->n_nodes_0D; i++) {
			printf("  p:fid[%.2d] = %.2d \n", f,e->point_node_list[f][i]);
		}printf("\n");
	}

	printf("QUAD \n");
	for (f=0; f<QUAD_EDGES; f++) {
		for (i=0; i<e->n_nodes_1D; i++) {
			printf("  e:fid[%.2d] = %.2d \n", f,e->edge_node_list[f][i]);
		}printf("\n");
	}

	printf("HEX \n");
	for (f=0; f<HEX_EDGES; f++) {
		for (i=0; i<e->n_nodes_2D; i++) {
			printf("  f:fid[%.2d] = %.2d \n", f,e->face_node_list[f][i]);
		}printf("\n");
	}
}

void testEQ2_QuadratureStokesVol(void)
{
	ConformingElementFamily e;
	int p,ngp2,ngp3;
	QPoint2d gp2[9];
	QPoint3d gp3[27];
	double sum;

	printf("TEST(%s)\n", __func__);

	ElementTypeCreate_Q2(&e,3);

	printf("QUAD (volQ)\n");
	e->generate_volume_quadrature_2D(&ngp2,gp2);
	printf("ngp = %d \n", ngp2 );
	sum = 0.0;
	for (p=0; p<ngp2; p++) {
		printf("  [%.2d] %+1.6lf %+1.6lf | %+1.6lf \n",p,gp2[p].xi,gp2[p].eta,gp2[p].w );
		sum = sum + gp2[p].w;
	}
	printf("  sum w                      %+1.6lf \n",sum );


	printf("HEX (volQ)\n");
	e->generate_volume_quadrature_3D(&ngp3,gp3);
	printf("ngp = %d \n", ngp3 );
	sum = 0.0;
	for (p=0; p<ngp3; p++) {
		printf("  [%.2d] %+1.6lf %+1.6lf %+1.6lf | %+1.6lf \n",p,gp3[p].xi,gp3[p].eta,gp3[p].zeta,gp3[p].w );
		sum = sum + gp3[p].w;
	}
	printf("  sum w                                %+1.6lf \n",sum );
}

void testEQ2_QuadratureStokesSurf(void)
{
	ConformingElementFamily e;
	int p,ngp1,ngp2,f;
	QPoint1d gp1[3];
	QPoint2d gp2[9];
	QPoint3d gp3[27];
	double sum;

	printf("TEST(%s)\n", __func__);

	ElementTypeCreate_Q2(&e,3);

	printf("QUAD->LINE (surfQ)\n");
	for (f=0; f<QUAD_EDGES; f++) {
		e->generate_surface_quadrature_2D(e,f,&ngp1,gp1,gp2);
		printf("face[%d] \n", f);
		printf("  ngp = %.2d \n", ngp1 );
		sum = 0.0;
		for (p=0; p<ngp1; p++) {
			printf("    [%.2d] (s)-[%+1.6lf | %+1.6lf] ; (xi,eta) [%+1.6lf %+1.6lf | %+1.6lf] \n",
						 p,gp1[p].xi,gp1[p].w, gp2[p].xi,gp2[p].eta,gp2[p].w );
			sum = sum + gp1[p].w;
		}
		printf("    sum w                 %+1.6lf \n",sum );
	}

	printf("HEX->QUAD (surfQ)\n");
	for (f=0; f<HEX_EDGES; f++) {
		e->generate_surface_quadrature_3D(e,f,&ngp2,gp2,gp3);
		printf("face[%.2d] \n", f);
		printf("  ngp = %d \n", ngp2 );
		sum = 0.0;
		for (p=0; p<ngp2; p++) {
			printf("    [%.2d] (s,t)-[%+1.6lf %+1.6lf | %+1.6lf] ; (xi,eta,zeta) [%+1.6lf %+1.6lf %+1.6lf | %+1.6lf] \n",
						 p,gp2[p].xi,gp2[p].eta,gp2[p].w, gp3[p].xi,gp3[p].eta,gp3[p].zeta,gp3[p].w );
			sum = sum + gp2[p].w;
		}
		printf("    sum w                             %+1.6f \n",sum );
	}
}

void testEQ2_SurfQuad2d_int_rotated_2( void )
{
	ConformingElementFamily e;
	int i,j;
	int cnt;
	double coords[9*2];
	double x,y;
	int ngp_s;
	QPoint1d gp1[3];
	QPoint2d gp2[3];
	int f,p;
	double sa,exact_sa;
	double dx, dy;
	double xi,yi,Lx,Ly,xf,yf;


	printf("TEST(%s)\n", __func__);

	ElementTypeCreate_Q2(&e,2);

	xi = -10.0; Lx = 24.0;
	yi = 1.0;		Ly = 4.0;

	xf = xi + Lx;
	yf = yi + Ly;

	dx = (xf-xi)/2.0;
	dy = (yf-yi)/2.0;

	cnt = 0;
	for( j=0; j<3; j++ ) {
		for( i=0; i<3; i++ ) {
			x = xi + i * dx;
			y = yi + j * dy;

			insert_point_2D( cnt, x,y, coords );

			cnt++;
		}
	}

	/* rotate the coordinates */
	for( i=0; i<9; i++ ) {
		double R[2][2];
		double xold[2],xnew[2];
		double theta = 30.0*M_PI/180.0;

		R[0][0] = cos(theta);	R[0][1] = -sin(theta);
		R[1][0] = sin(theta);	R[1][1] = cos(theta);


		xold[0] = coords[2*i  ];
		xold[1] = coords[2*i+1];

		xnew[0] = R[0][0]*xold[0] + R[0][1]*xold[1];
		xnew[1] = R[1][0]*xold[0] + R[1][1]*xold[1];

		coords[2*i+0] = xnew[0];
		coords[2*i+1] = xnew[1];
	}

	sa = 0.0;
	for( f=0; f<QUAD_EDGES; f++ ) {
		e->generate_surface_quadrature_2D( e, f, &ngp_s, gp1, gp2 );

		for( p=0; p<ngp_s; p++ ) {
			double n0[2], t0[2], J;

			e->compute_surface_geometry_2D(	e, coords,	// should contain 27 points with dimension 3 (x,y,z) //
																		 f,			// face node ids (size = 9) //
																		 &gp1[p],    // should contain 1 point with dimension 1 (xi)   //
																		 n0,t0, &J );	// n0[],t0 contains 1 point with dimension 2 (x,y) //
			if(p==0) {
				printf("f=%.2d: normal  = %+1.6lf, %+1.6lf \n", f, n0[0],n0[1] );
				printf("f=%.2d: tangent = %+1.6lf, %+1.6lf \n", f, t0[0],t0[1] );
			}

			sa = sa + gp1[p].w * 1.0 * J;
		}
	}
	exact_sa = 2.0*Lx + 2.0*Ly;
	printf("sa = %1.6f : exact = %1.6f \n", sa, exact_sa );
}

void evaluate_divF_sin( double xp[], double *div )
{
	double dFx,dFy,dFz;
	double a = 0.9;
	double b = 0.6;
	double alpha = 3.0;
	double c = 0.4;

	/* F = ( sin(ax), cos(by), alpha.sin(cz) ) */
	/* dfx = a.cos(ax) */
	/* dfy = -b.sin(by) */
	/* dfz = alpha.c.cos(cz) */

	dFx =  a*cos(a*xp[0]);
	dFy = -b*sin(b*xp[1]);
	dFz =  alpha*c*cos(c*xp[2]);

	*div = dFx + dFy + dFz;
}

void evaluate_F_sin( double xp[], double *Fx,double *Fy,double *Fz )
{
	double a = 0.9;
	double b = 0.6;
	double alpha = 3.0;
	double c = 0.4;

	/* F = ( sin(ax), cos(by), alpha.sin(cz) ) */
	/* dfx = a.cos(ax) */
	/* dfy = -b.sin(by) */
	/* dfz = alpha.c.cos(cz) */

	if (Fx) { *Fx = sin(a*xp[0]); }
	if (Fy) { *Fy = cos(b*xp[1]); }
	if (Fz) { *Fz = alpha*sin(c*xp[2]); }
}

void evaluate_F_poly( double xp[], double *Fx,double *Fy,double *Fz )
{
	double x=xp[0];
	double y=xp[1];
	double z=xp[2];

	if (Fx) { *Fx = 3.0*x*x + 2.0*x + y*y; }
	if (Fy) { *Fy = 5.0*y + x + z*z; }
	if (Fz) { *Fz = 10.0*x*x + 5.0*y + 7.0*z*z; }
}

void evaluate_divF_poly( double xp[], double *div )
{
	double dFx,dFy,dFz;
	double x=xp[0];
	double z=xp[2];

	dFx = 6.0*x + 2.0;
	dFy = 5.0;
	dFz = 14.0*z;

	*div = dFx + dFy + dFz;
}

void evaluate_divF( const int type, double xp[], double *div )
{
	switch (type) {
		case 0:
			evaluate_divF_poly(xp,div);
			break;

		case 1:
			evaluate_divF_sin(xp,div);
			break;

		default:
			evaluate_divF_poly(xp,div);
			break;
	}
}

void evaluate_F( const int type, double xp[], double *Fx,double *Fy,double *Fz )
{
	switch (type) {
		case 0:
			evaluate_F_poly(xp,Fx,Fy,Fz);
			break;

		case 1:
			evaluate_F_sin(xp,Fx,Fy,Fz);
			break;

		default:
			evaluate_F_poly(xp,Fx,Fy,Fz);
			break;
	}
}

void testEQ2_SurfQuad3d_2( void )
{
	ConformingElementFamily e;
	int i,j,k;
	int cnt;
	double coords[27*3];
	double x,y,z,x_p,y_p,z_p,XP[3];
	int ngp2,ngp3;
	QPoint2d gp2[9];
	QPoint3d gp3[27],surf_gp3[9];
	int f,p;
	double div_F,F_dot_n, fx_p,fy_p,fz_p,Fn;
	double dx, dy, dz;
	double xi,yi,zi,Lx,Ly,Lz,xf,yf,zf;
	double Ni3D[27],__GNi[27*3], __GNx[27*3];
	double *GNi[3];
	double *GNx[3];
	double rx,ry,rz,fac;
	int func_type = 0;


	/* init memory */
	GNi[0] = &__GNi[0];
	GNi[1] = &__GNi[27*1];
	GNi[2] = &__GNi[27*2];

	GNx[0] = &__GNx[0];
	GNx[1] = &__GNx[27*1];
	GNx[2] = &__GNx[27*2];

	printf("TEST(%s)\n", __func__);
	ElementTypeCreate_Q2(&e,3);

	fac = 1.0e-2;
	rx = fac * 2.0*rand()/(RAND_MAX+1.0) - 1.0;
	ry = fac * 2.0*rand()/(RAND_MAX+1.0) - 1.0;
	rz = fac * 2.0*rand()/(RAND_MAX+1.0) - 1.0;

	xi = -2.0 + rx;		Lx = 4.0;
	yi = -1.0 + ry;		Ly = 5.0;
	zi = -3.0 + rz;		Lz = 7.0;

	rx = fac * 2.0*rand()/(RAND_MAX+1.0) - 1.0;
	ry = fac * 2.0*rand()/(RAND_MAX+1.0) - 1.0;
	rz = fac * 2.0*rand()/(RAND_MAX+1.0) - 1.0;

	xf = xi + Lx + rx;
	yf = yi + Ly + ry;
	zf = zi + Lz + rz;

	dx = (xf-xi)/2.0;
	dy = (yf-yi)/2.0;
	dz = (zf-zi)/2.0;

	cnt = 0;
	for( k=0; k<3; k++ ) {
		for( j=0; j<3; j++ ) {
			for( i=0; i<3; i++ ) {
				x = xi + i * dx;
				y = yi + j * dy;
				z = zi + k * dz;

				insert_point_3D( cnt, x,y,z, coords );

				cnt++;
			}
		}
	}

	// optional
	/* rotate the coordinates about z axis */
	for( i=0; i<27; i++ ) {
		double R[3][3];
		double xold[3],xnew[3];
		double theta = 15.0*M_PI/180.0;

		R[0][0] = cos(theta);	R[0][1] = -sin(theta);	R[0][2] = 0.0;
		R[1][0] = sin(theta);	R[1][1] = cos(theta);	R[1][2] = 0.0;
		R[2][0] = 0.0;			R[2][1] = 0.0;			R[2][2] = 1.0;


		xold[0] = coords[3*i  ];
		xold[1] = coords[3*i+1];
		xold[2] = coords[3*i+2];

		xnew[0] = R[0][0]*xold[0] + R[0][1]*xold[1] + R[0][2]*xold[2];
		xnew[1] = R[1][0]*xold[0] + R[1][1]*xold[1] + R[1][2]*xold[2];
		xnew[2] = R[2][0]*xold[0] + R[2][1]*xold[1] + R[2][2]*xold[2];

		coords[3*i+0] = xnew[0];
		coords[3*i+1] = xnew[1];
		coords[3*i+2] = xnew[2];
	}

	/* perform the volume integration of div(F) */
	e->generate_volume_quadrature_3D(&ngp3,gp3);
	div_F = 0.0;
	for( p=0; p<ngp3; p++ ) {
		double div, det_J;

		/* interpolate coordinates to the gauss point */
		e->basis_NI_3D( &gp3[p], Ni3D );
		x_p = y_p = z_p = 0.0;
		for( k=0; k<e->n_nodes_3D; k++ ) {
			x_p += Ni3D[k] * coords[3*k+0];
			y_p += Ni3D[k] * coords[3*k+1];
			z_p += Ni3D[k] * coords[3*k+2];
		}

		e->basis_GNI_3D( &gp3[p], GNi );
		e->compute_volume_geometry_3D( coords, (const double**)GNi, GNx, &det_J );

		XP[0] = x_p;
		XP[1] = y_p;
		XP[2] = z_p;

		evaluate_divF(func_type,XP,&div);

		div_F = div_F + gp3[p].w * ( div ) * det_J;
	}
	printf("int_vol( div(F) ) = %+1.6lf \n", div_F );

	/* perform the surface integration of F.n */
	F_dot_n = 0.0;
	for( f=0; f<HEX_EDGES; f++ ) {
		e->generate_surface_quadrature_3D(e,f,&ngp2,gp2,surf_gp3);

		for( p=0; p<ngp2; p++ ) {
			double normal[3],tangent[3],det_surfJ;

			/* interpolate coordinates to the gauss point */
			e->basis_NI_3D(&surf_gp3[p],Ni3D);
			x_p = y_p = z_p = 0.0;
			for( k=0; k<27; k++ ) {
				x_p += Ni3D[k] * coords[3*k+0];
				y_p += Ni3D[k] * coords[3*k+1];
				z_p += Ni3D[k] * coords[3*k+2];
			}
			XP[0] = x_p;
			XP[1] = y_p;
			XP[2] = z_p;
			evaluate_F( func_type, XP, &fx_p, &fy_p, &fz_p );

			e->compute_surface_geometry_3D(	e, coords,	// should contain 27 points with dimension 3 (x,y,z) //
																		 f,			// face node ids (size = 9) //
																		 &gp2[p],    // should contain 1 point with dimension 2 (s,t)   //
																		 normal,tangent, &det_surfJ );	// n0[],t0 contains 1 point with dimension 3 (x,y,z) //

			if (p==0) {
				double r[3],r1[3];

				printf("[f=%d:p=%d]: normal     %+1.2e %+1.2e %+1.2e \n", f,p,normal[0],normal[1],normal[2]);
				printf("[f=%d:p=%d]: tangent(s) %+1.2e %+1.2e %+1.2e \n", f,p,tangent[0],tangent[1],tangent[2]);

				eval_cross_product(normal,tangent,r); // t2
				eval_cross_product(r,normal,r1);
				printf("[f=%d:p=%d]: t2     %+1.2e %+1.2e %+1.2e \n", f,p,r[0],r[1],r[2]);
				printf("[f=%d:p=%d]: t1     %+1.2e %+1.2e %+1.2e \n", f,p,r1[0],r1[1],r1[2]);


			}
			Fn = fx_p*normal[0] + fy_p*normal[1] + fz_p*normal[2];
			F_dot_n = F_dot_n + gp2[p].w * Fn * det_surfJ;

		}
		printf("\n");
	}
	printf("int_surf( F.n )   = %+1.6lf \n", F_dot_n );
}

/*
void main( void )
{
	testEQ2_FaceIds();
	testEQ2_QuadratureStokesVol();
	testEQ2_QuadratureStokesSurf();

	testEQ2_SurfQuad2d_1();
	testEQ2_SurfQuad2d_int_rotated_2();

	testEQ2_SurfQuad3d_1();
	testEQ2_SurfQuad3d_2();
}
*/

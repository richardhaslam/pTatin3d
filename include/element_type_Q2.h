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
 **    filename:   element_type_Q2.h
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

#ifndef __pTatin_element_type_Q2_h__
#define __pTatin_element_type_Q2_h__

#include "ptatin3d_defs.h"


typedef struct {
	double xi,eta,zeta;
	double w;
} QPoint3d;
typedef struct {
	double xi,eta;
	double w;
} QPoint2d;
typedef struct {
	double xi;
	double w;
} QPoint1d;

typedef struct _p_ConformingElementFamily* ConformingElementFamily;

typedef enum { LINE_EDGES=1, QUAD_EDGES=4, HEX_EDGES=6 } HexElementFaceCount;

typedef enum { HEX_FACE_Pxi=0, HEX_FACE_Nxi, HEX_FACE_Peta, HEX_FACE_Neta, HEX_FACE_Pzeta, HEX_FACE_Nzeta } HexElementFace;
typedef enum { QUAD_EDGE_Pxi=0, QUAD_EDGE_Nxi, QUAD_EDGE_Peta, QUAD_EDGE_Neta } QuadElementEdge;
typedef enum { LINE_POINT_Pxi=0, LINE_POINT_Nxi } LineElementPoint;

struct _p_ConformingElementFamily {
	char *name;
	int nsd;

	int n_nodes_0D;
	int n_nodes_1D;
	int n_nodes_2D;
	int n_nodes_3D;

	void (*basis_NI_1D)(const QPoint1d*,double*); // xi[],Ni[]
	void (*basis_NI_2D)(const QPoint2d*,double*);
	void (*basis_NI_3D)(const QPoint3d*,double*);

	void (*basis_GNI_1D)(const QPoint1d*,double**); // xi[],GNi[]
	void (*basis_GNI_2D)(const QPoint2d*,double**);
	void (*basis_GNI_3D)(const QPoint3d*,double**);

	// 1D only
	int n_points;
	int **point_node_list;
	void (*generate_all_point_ids)(int**);
	void (*extract_point_field)(ConformingElementFamily,int,int,double*,double*); // Element,face,dofs,e_data[],edge_e_data[]

	// 2D only
	int n_edges;
	int **edge_node_list;
	void (*generate_all_edge_ids)(int**);
	void (*extract_edge_field)(ConformingElementFamily,QuadElementEdge,int,double*,double*); // Element,face,dofs,e_data[],edge_e_data[]

	// 3D only
	int n_faces;
	int **face_node_list;
	void (*generate_all_face_ids)(int**);
	void (*extract_surface_field)(ConformingElementFamily,HexElementFace,int,double*,double*); // Element,face,dofs,e_data[],surface_e_data[]

	void (*generate_surface_quadrature_2D)(ConformingElementFamily,QuadElementEdge,int*,QPoint1d*,QPoint2d*); // face_id,ngp_s,w[],xi_st[],xi_xyz[]
	void (*generate_volume_quadrature_2D)(int*,QPoint2d*); // ngp,w[],xi_xyz[]
	void (*compute_surface_geometry_2D)(ConformingElementFamily,const double*,QuadElementEdge,const QPoint1d*,double*,double*,double*); // coords[],edge_idx,,xi[],normal[],J
	void (*compute_volume_geometry_2D)(const double*,const double**,double**,double*); // coords[],GNi[][],GNx[][],J
	void (*compute_surface_normal_2D)(ConformingElementFamily,const double*,QuadElementEdge,const QPoint1d*,double*);
	void (*compute_surface_tangent_2D)(ConformingElementFamily,const double*,QuadElementEdge,const QPoint1d*,double*);

	void (*generate_surface_quadrature_3D)(ConformingElementFamily,HexElementFace,int*,QPoint2d*,QPoint3d*);
	void (*generate_volume_quadrature_3D)(int*,QPoint3d*);
	void (*compute_surface_geometry_3D)(ConformingElementFamily,const double*,HexElementFace,const QPoint2d*,double*,double*,double*);
	void (*compute_volume_geometry_3D)(const double*,const double**,double**,double*);
	void (*compute_surface_normal_3D)(ConformingElementFamily,const double*,HexElementFace,const QPoint2d*,double*);
	void (*compute_surface_tangents_3D)(ConformingElementFamily,const double*,HexElementFace,const QPoint2d*,double*,double*);

	/*
	 3D  2D  1D
	 volume   Y		X		X
	 face     Y		Y		X
	 edge     Y		Y		Y
	 point    Y		Y		Y
	 */

};


void ElementTypeCreate_Q2(ConformingElementFamily *_e,const int dim);
void ElementTypeDestroy_Q2(ConformingElementFamily *_e);

void extract_element_surface_field_Q2_3D( int face_id, int fid_node_list[][9], int dofs, double e_data[], double surf_e_data[] );

void pTatin_ConstructNi_Q2_3D( double _xi[], double Ni[] );
void pTatinConstructGNi_Q2_3D( double _xi[], double GNi[3][Q2_NODES_PER_EL_3D] );

void ElementHelper_matrix_inverse_2x2(double A[2][2],double B[2][2]);
void ElementHelper_matrix_inverse_3x3(double A[3][3],double B[3][3]);

#endif


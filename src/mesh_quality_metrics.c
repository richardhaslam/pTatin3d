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
 **    Filename:      mesh_quality_metrics.c
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

/*

 ## Mesh quality metric guidlines ##
 
 MESH_QUALITY_ASPECT_RATIO (AR)
 * Possible range of values:  1.0 < AR < \infinty 
 * Recommended bounds:        1.0 < AR < 10
 
 MESH_QUALITY_DISTORTION (D)
 * Possible range of values:  -\infinty < D < \infinty 
 * Recommended bounds:        0.7 < D < 1.0
 
 MESH_QUALITY_DIAGONAL_RATIO (DR)
 * Possible range of values:  0.0 < DR < 1.0
 * Recommended bounds:        0.5 < DR < 1.0

 MESH_QUALITY_VERTEX_ANGLE (VA)
 * Possible range of values:  -1.0 < VA < 1.0 
 * Recommended bounds:         0.5 < VA < 1.0
 
 MESH_QUALITY_FACE_AREA_RATIO (FAR)
 * Possible range of values:  0.0 < FAR < 1.0 
 * Recommended bounds:        0.5 < FAR < 1.0

*/ 
 
/*  
 Functions contributed by Thomas Philippe [phthomas@student.ethz.ch] 

 DMDAComputeMeshQualityMetric_FaceAreaRatio()
 DMDAComputeMeshQualityMetric_VertexAngle()
*/

#include "stdio.h"
#include "stdlib.h"
#include "petsc.h"
#include "petscdm.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_element_q2p1.h"
#include "quadrature.h"
#include "element_utils_q2.h"

#include "mesh_quality_metrics.h"


/*
#define Q2_FACE_WEST   2
#define Q2_FACE_EAST   3
#define Q2_FACE_SOUTH  5
#define Q2_FACE_NORTH  6
#define Q2_FACE_BACK   1
#define Q2_FACE_FRONT  4
*/

/* protoypes for helpers */
void get_node_coordinate(double el_coords[],int index,double pos[]);
void get_face_coordinates(double el_coords[],Q2Face face,double coords[][3]);
void get_vertex_vectors(double el_coords[],int vertex,double vectors[][3]);
double compute_quadrilateral_area_(double vertex[][3]);
void compute_cossin(double vectors[][3], double cossin[]);
double compute_distance3(double posA[],double posB[]);
double compute_norm3(double a[]);
double compute_dot_product3(double a[],double b[]);
void compute_quadrilateral_area_approximate(double quad_coords[],double *area);
void compute_cross_product3(double a[],double b[],double crossprod[]);


void get_node_coordinate(double el_coords[],int index,double pos[])
{
	pos[0] = el_coords[3*index+0];
	pos[1] = el_coords[3*index+1];
	pos[2] = el_coords[3*index+2];
}

/*
 Get the face node coordinates according to numbering:
 3----2
 |    |
 |    |
 0----1
 with the normal pointing outward.
 F_BACK:  0,  6,  8,  2
 F_WEST:  0,  18, 24, 6
 F_EAST:  2,  8,  26, 20
 F_FRONT: 18, 20, 26, 24
 F_SOUTH: 0,  2,  20, 18
 F_NORTH: 6,  24, 26, 8
 */
void get_face_coordinates(double el_coords[],Q2Face face,double coords[][3])
{
	int nodes[4], i;

	nodes[0] = nodes[1] = nodes[2] = nodes[3] = -1;	
	switch (face) {
			
		case Q2_FACE_WEST:
			nodes[0] = 0;  nodes[1] = 18;  nodes[2] = 24;  nodes[3] = 6;
			break;
			
		case Q2_FACE_EAST:
			nodes[0] = 2;  nodes[1] = 8;   nodes[2] = 26;  nodes[3] = 20;
			break;
			
		case Q2_FACE_SOUTH:
			nodes[0] = 0;  nodes[1] = 2;   nodes[2] = 20;  nodes[3] = 18;
			break;

		case Q2_FACE_NORTH:
			nodes[0] = 6;  nodes[1] = 24;  nodes[2] = 26;  nodes[3] = 8;
			break;
			
		case Q2_FACE_FRONT:
			nodes[0] = 18; nodes[1] = 20;  nodes[2] = 26;  nodes[3] = 24;
			break;
			
		case Q2_FACE_BACK:
			nodes[0] = 0;  nodes[1] = 6;   nodes[2] = 8;   nodes[3] = 2;
			break;			
	}
	
	for(i = 0; i < 4; i++){
		coords[i][0] = el_coords[3*nodes[i] + 0];
		coords[i][1] = el_coords[3*nodes[i] + 1];
		coords[i][2] = el_coords[3*nodes[i] + 2];
	}
}

/*
 The corner numbering is the corresponding node numbering. The vectors are right-handed oriented.
 
 C0: 0-18, 0-2, 0-6 
 C2: 2-0, 2-20, 2-18
 C8: 8-26, 8-6, 8-2
 C6: 6-8, 6-24, 6-0
 C18: 18-20, 18-0, 18-24
 C20: 20-2, 20-18, 20-26
 C24: 24-6, 24-26, 24-18
 C26: 26-24, 26-8, 26-20
*/
void get_vertex_vectors(double el_coords[],int vertex,double vectors[][3])
{
	int nodes[4], i;

	nodes[0] = nodes[1] = nodes[2] = nodes[3] = -1;	
	switch (vertex) {
			
		case Q2_VERTEX_0:
			nodes[0] = 0;   nodes[1] = 18;  nodes[2] = 2;   nodes[3] = 6;
			break;
			
		case Q2_VERTEX_1:
			nodes[0] = 2;   nodes[1] = 0;   nodes[2] = 20;  nodes[3] = 8;
			break;
			
		case Q2_VERTEX_2:
			nodes[0] = 8;   nodes[1] = 26;  nodes[2] = 6;   nodes[3] = 2;
			break;
	
		case Q2_VERTEX_3:
			nodes[0] = 6;   nodes[1] = 8;   nodes[2] = 24;  nodes[3] = 0;
			break;
			
		case Q2_VERTEX_4:
			nodes[0] = 18;  nodes[1] = 20;  nodes[2] = 0;   nodes[3] = 24;
			break;
			
		case Q2_VERTEX_5:
			nodes[0] = 20;  nodes[1] = 2;   nodes[2] = 18;  nodes[3] = 26;
			break;
		
		case Q2_VERTEX_6:
			nodes[0] = 24;  nodes[1] = 6;   nodes[2] = 26;  nodes[3] = 18;
			break;		
		
		case Q2_VERTEX_7:
			nodes[0] = 26;  nodes[1] = 24;  nodes[2] = 8;   nodes[3] = 20;
			break;					
	}
	for(i = 0; i < 3; i++){
		vectors[i][0] = el_coords[3*nodes[i+1] + 0] - el_coords[3*nodes[0] + 0];
		vectors[i][1] = el_coords[3*nodes[i+1] + 1] - el_coords[3*nodes[0] + 1];
		vectors[i][2] = el_coords[3*nodes[i+1] + 2] - el_coords[3*nodes[0] + 2];
	}  
}

/* 
 Compute a x b
*/
void compute_cross_product3(double a[],double b[],double crossprod[])
{
	crossprod[0] =    a[1]*b[2] - a[2]*b[1];
	crossprod[1] = -( a[0]*b[2] - a[2]*b[0] );
	crossprod[2] =    a[0]*b[1] - a[1]*b[0];
}

double compute_quadrilateral_area_(double vertex[][3])
{
	int i, c1, c2;
	double crossprods[4][3], scalarprod[4], normcrossprod[4], a[3], b[3];
	
	for(i = 0; i<4; i++){
		c1 = (i+1)%4;
		c2 = (i+3)%4;

        a[0] = vertex[c1][0] - vertex[i][0]; a[1] = vertex[c1][1] - vertex[i][1];  a[2] = vertex[c1][2] - vertex[i][2]; 
        b[0] = vertex[c2][0] - vertex[i][0]; b[1] = vertex[c2][1] - vertex[i][1];  b[2] = vertex[c2][2] - vertex[i][2]; 
        
        compute_cross_product3(a,b,crossprods[i]);
        scalarprod[i] = compute_dot_product3(crossprods[i],a);
		normcrossprod[i] = compute_norm3(crossprods[i]);
	}
	for(i = 0; i<4; i++){
		if(scalarprod[i] > 0.0 && scalarprod[(i+2)%4] > 0.0){
			break;
		}
	}
	return 0.5*(sqrt(normcrossprod[i]) + sqrt(normcrossprod[(i+2)%4]));
}

/*
 Given three vectors a,b,c, compute ||c.(a x b)|| = |a|*|b|*|c|*cos(alpha)*sin(theta),
 it returns cossin[0] = cos(alpha); cossin[1] = sin(theta)
 */
void compute_cossin(double vectors[][3], double cossin[])
{
	double crossprod[3], a[3], b[3], c[3], scalarprod, ncross, n0, n1, n2;
	
  a[0] = vectors[0][0]; a[1] = vectors[0][1]; a[2] = vectors[0][2];
  b[0] = vectors[1][0]; b[1] = vectors[1][1]; b[2] = vectors[1][2];
  c[0] = vectors[2][0]; c[1] = vectors[2][1]; c[2] = vectors[2][2];
    
	compute_cross_product3(a,b,crossprod);    

	ncross = compute_norm3(crossprod);
	n0 = compute_norm3(a);
	n1 = compute_norm3(b);
	n2 = compute_norm3(c);
    
	//scalarprod = |c|.|a|.|b| sin(theta)*cos(alpha)
	
	scalarprod = compute_dot_product3(crossprod,c);
	cossin[1] = ncross/(n0*n1);
	cossin[0] = scalarprod/(n0*n1*n2*cossin[1]); 
//    ierr = PetscPrintf(PETSC_COMM_SELF,"---> scalarprod: %f ncross:%f n0:%f n1:%f n2:%f cossin0:%f cossin1:%f   c(%f, %f, %f) cross(%f, %f, %f)\n", scalarprod, ncross, n0, n1, n2, cossin[0], cossin[1], c[0],c[1],c[2],  crossprod[0],crossprod[1],crossprod[2]);CHKERRQ(ierr);
//    ierr = PetscPrintf(PETSC_COMM_SELF,"---> scalarprod: %f c(%f, %f, %f) cross(%f, %f, %f)  a(%f, %f, %f)  b(%f, %f, %f)\n", scalarprod,c[0],c[1],c[2],  crossprod[0],crossprod[1],crossprod[2], a[0],a[1],a[2], b[0],b[1],b[2]);CHKERRQ(ierr);
}

double compute_distance3(double posA[],double posB[])
{
	double val;
	val = sqrt(  (posA[0]-posB[0])*(posA[0]-posB[0])  +  (posA[1]-posB[1])*(posA[1]-posB[1])  +  (posA[2]-posB[2])*(posA[2]-posB[2])  );
	return val;
}

double compute_norm3(double a[])
{
	double val;
	val = sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
	return val;
}

double compute_dot_product3(double a[],double b[])
{
	double val;
	val = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	return val;
}


/*
 Compute area of quadrilateral defined in 3 space. 
 
 3------ 2
 |       | 
 |       |
 |       |
 0 ----- 1
 
 area( trianlge [0,1,2] ) + area( triangle[0,2,3] )
 */
void compute_quadrilateral_area_approximate(double quad_coords[],double *area)
{
	double triangle_012[3][3];
	double triangle_023[3][3];
	double b,l,h,a_012,a_023,theta;
	
	get_node_coordinate( quad_coords, 0, triangle_012[0] );
	get_node_coordinate( quad_coords, 1, triangle_012[1] );
	get_node_coordinate( quad_coords, 2, triangle_012[2] );
	
	b = compute_distance3( triangle_012[1], triangle_012[0] ); /* length(01) */
	l = compute_distance3( triangle_012[2], triangle_012[0] ); /* length(02) */
	theta = acos(b/l);
	h = b * tan(theta);
	a_012 = 0.5 * b * h;
	
	get_node_coordinate( quad_coords, 0, triangle_023[0] );
	get_node_coordinate( quad_coords, 2, triangle_023[1] );
	get_node_coordinate( quad_coords, 3, triangle_023[2] );
	
	b = compute_distance3( triangle_023[1], triangle_023[0] ); /* length(02) */
	l = compute_distance3( triangle_023[2], triangle_023[0] ); /* length(03) */
	theta = acos(b/l);
	h = b * tan(theta);
	a_023 = 0.5 * b * h;
	
	*area = a_012 + a_023;
}

/*
  max_(over all elements) [  max_face_length_e / min_face_length_e ]
*/
#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetric_AspectRatio"
PetscErrorCode DMDAComputeMeshQualityMetric_AspectRatio(DM dm,PetscReal *value)
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e;
	const PetscInt  *el_nidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          dx,dy,dz,dl_min,dl_max,el_ar,max_ar,max_ar_g;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;

	/* setup for coords */
	ierr = DMGetCoordinateDM(dm,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	
	max_ar = -1.0e32;
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		
		dl_min = 1.0e32;
		dl_max = 1.0e-32;

		dx = fabs( el_coords[3*Q2_FACE_NODE_EAST +0] - el_coords[3*Q2_FACE_NODE_WEST +0]  );
		dy = fabs( el_coords[3*Q2_FACE_NODE_NORTH+1] - el_coords[3*Q2_FACE_NODE_SOUTH+1] );
		dz = fabs( el_coords[3*Q2_FACE_NODE_FRONT+2] - el_coords[3*Q2_FACE_NODE_BACK +2]  );
		//printf("e=%d: %1.4e %1.4e %1.4e \n", e, dx, dy, dz );
		
		if (dx < dl_min) { dl_min = dx; }
		if (dy < dl_min) { dl_min = dy; }
		if (dz < dl_min) { dl_min = dz; }

		if (dx > dl_max) { dl_max = dx; }
		if (dy > dl_max) { dl_max = dy; }
		if (dz > dl_max) { dl_max = dz; }
		
		//printf("dl_min %1.4e \n", dl_min );
		el_ar = dl_max / dl_min;
		if (el_ar > max_ar) { max_ar = el_ar; }

	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

	ierr = MPI_Allreduce(&max_ar,&max_ar_g,1,MPI_DOUBLE,MPI_MAX,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
											 
	*value = (PetscReal)max_ar_g;
		
	PetscFunctionReturn(0);
}

/*
 min_(over all elements) [  8 x min(|J|_e) / vol_e ]
*/
#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetric_Distortion"
PetscErrorCode DMDAComputeMeshQualityMetric_Distortion(DM dm,PetscReal *value)
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e,p;
	const PetscInt  *el_nidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          el_dist,min_dist,min_dist_g,el_vol;
	PetscInt        ngp;
	PetscReal       WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal       min_detJ,detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* setup quadrature */

	ngp = 27;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);

	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dm,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	
	min_dist = 1.0e32;
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);

		P3D_evaluate_geometry_elementQ2(ngp,el_coords,GNI, detJ,dNudx,dNudy,dNudz);

		el_vol = 0.0;
		min_detJ = 1.0e32;
		for (p=0; p<ngp; p++) {
			el_vol = el_vol + 1.0 * WEIGHT[p] * detJ[p];
			if (detJ[p] < min_detJ) {
				min_detJ = detJ[p];
			}
		}
		
		el_dist = 8.0 * min_detJ / el_vol;
		
		if (el_dist < min_dist) { min_dist = el_dist; }
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&min_dist,&min_dist_g,1,MPI_DOUBLE,MPI_MIN,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
	
	*value = (PetscReal)(min_dist_g);
	
	PetscFunctionReturn(0);
}

/*
 min_(over all elements) [  min_diagonal_length_e / max_diagonal_length_e ]
 */
#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetric_DiagonalRatio"
PetscErrorCode DMDAComputeMeshQualityMetric_DiagonalRatio(DM dm,PetscReal *value)
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e;
	const PetscInt  *el_nidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          posA[3],posB[3];
	double          diag,dl_min,dl_max,dr,dr_min,dr_min_g;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dm,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	
	dr_min = 1.0e32;
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		
		dl_min = 1.0e32;
		dl_max = -1.0e32;
		
		get_node_coordinate(el_coords,Q2_VERTEX_0,posA);
		get_node_coordinate(el_coords,Q2_VERTEX_7,posB);
		diag =	compute_distance3(posA,posB);
		if (diag < dl_min) { dl_min = diag; }
		if (diag > dl_max) { dl_max = diag; }

		get_node_coordinate(el_coords,Q2_VERTEX_1,posA);
		get_node_coordinate(el_coords,Q2_VERTEX_6,posB);
		diag =	compute_distance3(posA,posB);
		if (diag < dl_min) { dl_min = diag; }
		if (diag > dl_max) { dl_max = diag; }

		get_node_coordinate(el_coords,Q2_VERTEX_5,posA);
		get_node_coordinate(el_coords,Q2_VERTEX_2,posB);
		diag =	compute_distance3(posA,posB);
		if (diag < dl_min) { dl_min = diag; }
		if (diag > dl_max) { dl_max = diag; }
		
		get_node_coordinate(el_coords,Q2_VERTEX_4,posA);
		get_node_coordinate(el_coords,Q2_VERTEX_3,posB);
		diag =	compute_distance3(posA,posB);
		if (diag < dl_min) { dl_min = diag; }
		if (diag > dl_max) { dl_max = diag; }

		dr = dl_min / dl_max;
		if (dr < dr_min) { dr_min = dr; }
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&dr_min,&dr_min_g,1,MPI_DOUBLE,MPI_MIN,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
	
	*value = (PetscReal)dr_min_g;
	
	PetscFunctionReturn(0);
}

/*
 min_(over all elements) [  min_{each face pair} area_1/area_0  ]
 */
#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetric_FaceAreaRatio"
PetscErrorCode DMDAComputeMeshQualityMetric_FaceAreaRatio(DM dm,PetscReal *value)
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e;
	const PetscInt  *el_nidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          a1,a2,vertices[4][3];
	double          ratio,ratio_min,ratio_min_g;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dm,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	

	ratio_min = 1.0e32;
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
        //Loops through face pairs.
		get_face_coordinates(el_coords, Q2_FACE_WEST, vertices);
		a1 = compute_quadrilateral_area_(vertices);
    get_face_coordinates(el_coords, Q2_FACE_EAST, vertices);
		a2 = compute_quadrilateral_area_(vertices);
    ratio = (a1<a2)?a1/a2:a2/a1;
    ratio_min = (ratio<ratio_min)?ratio:ratio_min;

		get_face_coordinates(el_coords, Q2_FACE_SOUTH, vertices);
		a1 = compute_quadrilateral_area_(vertices);
		get_face_coordinates(el_coords, Q2_FACE_NORTH, vertices);
		a2 = compute_quadrilateral_area_(vertices);
		ratio = (a1<a2)?a1/a2:a2/a1;
		ratio_min = (ratio<ratio_min)?ratio:ratio_min;        
		
		get_face_coordinates(el_coords, Q2_FACE_BACK, vertices);
		a1 = compute_quadrilateral_area_(vertices);
		get_face_coordinates(el_coords, Q2_FACE_FRONT, vertices);
		a2 = compute_quadrilateral_area_(vertices);
		ratio = (a1<a2)?a1/a2:a2/a1;
		ratio_min = (ratio<ratio_min)?ratio:ratio_min;            
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&ratio_min,&ratio_min_g,1,MPI_DOUBLE,MPI_MIN,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
	
	*value = (PetscReal)ratio_min_g;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetric_VertexAngle"
PetscErrorCode DMDAComputeMeshQualityMetric_VertexAngle(DM dm,PetscReal *value)
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e;
	const PetscInt  *el_nidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          val, val_min, val_min_g, cossin[2], vectors[3][3];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;

	/* setup for coords */
	ierr = DMGetCoordinateDM(dm,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	val_min = 2.0; 
	for (e=0;e<nel;e++) {
	    ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
        //Loops through vertexes
        get_vertex_vectors(el_coords,Q2_VERTEX_0,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val;
        
        get_vertex_vectors(el_coords,Q2_VERTEX_1,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val;
        
        get_vertex_vectors(el_coords,Q2_VERTEX_2,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_vertex_vectors(el_coords,Q2_VERTEX_3,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_vertex_vectors(el_coords,Q2_VERTEX_4,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_vertex_vectors(el_coords,Q2_VERTEX_5,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_vertex_vectors(el_coords,Q2_VERTEX_6,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_vertex_vectors(el_coords,Q2_VERTEX_7,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
    
	ierr = MPI_Allreduce(&val_min,&val_min_g,1,MPI_DOUBLE,MPI_MIN,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
    
	*value = (PetscReal)val_min_g;
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetric"
PetscErrorCode DMDAComputeMeshQualityMetric(DM dm,MeshQualityMeasure measure,PetscReal *value)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	switch (measure) {

		case MESH_QUALITY_ASPECT_RATIO:
			ierr = DMDAComputeMeshQualityMetric_AspectRatio(dm,value);CHKERRQ(ierr);
			break;
		
		case MESH_QUALITY_DISTORTION:
			ierr = DMDAComputeMeshQualityMetric_Distortion(dm,value);CHKERRQ(ierr);
			break;
			
		case MESH_QUALITY_DIAGONAL_RATIO:
			ierr = DMDAComputeMeshQualityMetric_DiagonalRatio(dm,value);CHKERRQ(ierr);
			break;
        
		case MESH_QUALITY_VERTEX_ANGLE:
			ierr = DMDAComputeMeshQualityMetric_VertexAngle(dm,value);CHKERRQ(ierr);
			break;
    
    case MESH_QUALITY_FACE_AREA_RATIO:
			ierr = DMDAComputeMeshQualityMetric_FaceAreaRatio(dm,value);CHKERRQ(ierr);
			break;            
			
	}
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetricList"
PetscErrorCode DMDAComputeMeshQualityMetricList(DM dm,const PetscInt nmeasures,MeshQualityMeasure measure[],PetscReal value[])
{
	PetscInt       i;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	for (i=0; i<nmeasures; i++) {
		ierr = DMDAComputeMeshQualityMetric(dm,measure[i],&value[i]);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

/*
 Output:
 
 range[]:       contains lower and upper bounds for the given metric
 valid_range[]: contains recommended lower and upper bounds for the given metric
*/
#undef __FUNCT__
#define __FUNCT__ "DMDMeshQualityMetricGetInfo"
PetscErrorCode DMDMeshQualityMetricGetInfo(DM dm,MeshQualityMeasure measure,PetscReal range[],PetscReal valid_range[])
{	
	PetscFunctionBegin;
	switch (measure) {
			
		case MESH_QUALITY_ASPECT_RATIO:
			range[0] = 1.0;
			range[1] = FP_INFINITE;
			valid_range[0] = 1.0;
			valid_range[1] = 10.0;
			break;
			
		case MESH_QUALITY_DISTORTION:
			range[0] = -FP_INFINITE;
			range[1] =  FP_INFINITE;
			valid_range[0] = 0.7;
			valid_range[1] = 1.0;
			break;
			
		case MESH_QUALITY_DIAGONAL_RATIO:
			range[0] = 0.0;
			range[1] = 1.0;
			valid_range[0] = 0.5;
			valid_range[1] = 1.0;
			break;
			
		case MESH_QUALITY_VERTEX_ANGLE:
			range[0] = -1.0;
			range[1] =  1.0;
			valid_range[0] = 0.5;
			valid_range[1] = 1.0;
			break;
			
    case MESH_QUALITY_FACE_AREA_RATIO:
			range[0] = 0.0;
			range[1] = 1.0;
			valid_range[0] = 0.5;
			valid_range[1] = 1.0;
			break;            
			
	}
	
	PetscFunctionReturn(0);
}


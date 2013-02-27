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



#define Q2_CELL_NODE_CENTER 13

#define Q2_FACE_NODE_WEST   12
#define Q2_FACE_NODE_EAST   14
#define Q2_FACE_NODE_SOUTH  10
#define Q2_FACE_NODE_NORTH  16
#define Q2_FACE_NODE_BACK   4
#define Q2_FACE_NODE_FRONT  22

#define Q2_VERTEX_0   0
#define Q2_VERTEX_1   2
#define Q2_VERTEX_2   6
#define Q2_VERTEX_3   8
#define Q2_VERTEX_4   18
#define Q2_VERTEX_5   20
#define Q2_VERTEX_6   24
#define Q2_VERTEX_7   26

typedef enum {
	Q2_FACE_WEST = 1,
	Q2_FACE_EAST,
	Q2_FACE_SOUTH,
  Q2_FACE_NORTH,
  Q2_FACE_BACK,
  Q2_FACE_FRONT
} Q2Face;

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
void get_face_coordinates(double el_coords[],Q2Face face,double coords[]);
void get_corner_vectors(double el_coords[],int corner,double vectors[]);
double compute_quadrilateral_area_(double coners[]);
void compute_cossin(double vectors[], double cossin[]);
double compute_distance3(double posA[],double posB[]);
void compute_quadrilateral_area_approximate(double quad_coords[],double *area);


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
void get_face_coordinates(double el_coords[],Q2Face face,double coords[])
{
	int nodes[4], i;
	
	switch (face) {
			
		case Q2_FACE_WEST:
			nodes[0] = 0; nodes[1] = 18; nodes[2] = 24; nodes[3] = 6;
			break;
			
		case Q2_FACE_EAST:
			nodes[0] = 2; nodes[1] = 8; nodes[2] = 26; nodes[3] = 20;
			break;
			
		case Q2_FACE_SOUTH:
			nodes[0] = 0; nodes[1] = 2; nodes[2] = 20; nodes[3] = 18;
			break;
		case Q2_FACE_NORTH:
			nodes[0] = 6; nodes[1] = 24; nodes[2] = 26; nodes[3] = 8;
			break;
			
		case Q2_FACE_FRONT:
			nodes[0] = 18; nodes[1] = 20; nodes[2] = 26; nodes[3] = 24;
			break;
			
		case Q2_FACE_BACK:
			nodes[0] = 0; nodes[1] = 6; nodes[2] = 8; nodes[3] = 2;
			break;			
	}
	
	for(i = 0; i < 4; i++){
		coords[i*3 + 0] = el_coords[3*nodes[i] + 0];
		coords[i*3 + 1] = el_coords[3*nodes[i] + 1];
		coords[i*3 + 2] = el_coords[3*nodes[i] + 2];
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
void get_corner_vectors(double el_coords[],int corner,double vectors[])
{
	int nodes[4], i;
	
	switch (corner) {
			
		case 0:
			nodes[0] = 0; nodes[1] = 18; nodes[2] = 2; nodes[3] = 6;
			break;
			
		case 2:
			nodes[0] = 2; nodes[1] = 0; nodes[2] = 20; nodes[3] = 18;
			break;
			
		case 8:
			nodes[0] = 8; nodes[1] = 26; nodes[2] = 6; nodes[3] = 2;
			break;
		case 6:
			nodes[0] = 6; nodes[1] = 8; nodes[2] = 24; nodes[3] = 0;
			break;
			
		case 18:
			nodes[0] = 18; nodes[1] = 20; nodes[2] = 0; nodes[3] = 24;
			break;
			
		case 20:
			nodes[0] = 20; nodes[1] = 2; nodes[2] = 18; nodes[3] = 26;
			break;
		case 24:
			nodes[0] = 24; nodes[1] = 6; nodes[2] = 26; nodes[3] = 18;
			break;		
		case 26:
			nodes[0] = 26; nodes[1] = 24; nodes[2] = 8; nodes[3] = 20;
			break;					
	}
	for(i = 0; i < 3; i++){
		vectors[i*3 + 0] = el_coords[3*nodes[i+1] + 0] - el_coords[3*nodes[0] + 0];
		vectors[i*3 + 1] = el_coords[3*nodes[i+1] + 1] - el_coords[3*nodes[0] + 1];
		vectors[i*3 + 2] = el_coords[3*nodes[i+1] + 2] - el_coords[3*nodes[0] + 2];
	}  
}

double compute_quadrilateral_area_(double coners[])
{
	int i, c1, c2, idx;
	double crossprod[12], scalarprod[4], normcrossprod[4];
	
	for(i = 0; i<4; i++){
		c1 = (i+1)%4;
		c2 = (i+3)%4;
		idx = 0;
		crossprod[3*i+idx] = (coners[c1*3+(1+idx)%3] - coners[i*3+(1+idx)%3])*(coners[c2*3+(2+idx)%3] - coners[i*3+(2+idx)%3])
		- (coners[c2*3+(1+idx)%3] - coners[i*3+(1+idx)%3])*(coners[c1*3+(2+idx)%3] - coners[i*3+(1+idx)%3]); ++idx;
		crossprod[3*i+idx] = (coners[c1*3+(1+idx)%3] - coners[i*3+(1+idx)%3])*(coners[c2*3+(2+idx)%3] - coners[i*3+(2+idx)%3])
		- (coners[c2*3+(1+idx)%3] - coners[i*3+(1+idx)%3])*(coners[c1*3+(2+idx)%3] - coners[i*3+(1+idx)%3]); ++idx;
		crossprod[3*i+idx] = (coners[c1*3+(1+idx)%3] - coners[i*3+(1+idx)%3])*(coners[c2*3+(2+idx)%3] - coners[i*3+(2+idx)%3])
		- (coners[c2*3+(1+idx)%3] - coners[i*3+(1+idx)%3])*(coners[c1*3+(2+idx)%3] - coners[i*3+(1+idx)%3]);        
		idx = 0;
		scalarprod[i] = crossprod[3*i+idx]*(coners[c1*3+idx] - coners[i*3+idx]); ++idx;
		scalarprod[i] += crossprod[3*i+idx]*(coners[c1*3+idx] - coners[i*3+idx]); ++idx;
		scalarprod[i] += crossprod[3*i+idx]*(coners[c1*3+idx] - coners[i*3+idx]); 
		
		normcrossprod[i] = crossprod[3*i+0]*crossprod[3*i+0] + crossprod[3*i+1]*crossprod[3*i+1] +crossprod[3*i+2]*crossprod[3*i+2];
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
void compute_cossin(double vectors[], double cossin[])
{
	int idx;
	double crossprod[3], scalarprod, ncross, n0, n1, n2;
	
	idx = 0;
	crossprod[idx] = vectors[0*3+(1+idx)%3]*vectors[1*3+(2+idx)%3] - vectors[1*3+(1+idx)%3]*vectors[0*3+(2+idx)%3]; ++idx;
	crossprod[idx] = vectors[0*3+(1+idx)%3]*vectors[1*3+(2+idx)%3] - vectors[1*3+(1+idx)%3]*vectors[0*3+(2+idx)%3]; ++idx;
	crossprod[idx] = vectors[0*3+(1+idx)%3]*vectors[1*3+(2+idx)%3] - vectors[1*3+(1+idx)%3]*vectors[0*3+(2+idx)%3]; 
	
	
	ncross = sqrt(crossprod[0]*crossprod[0] + crossprod[1]*crossprod[1] + crossprod[2]*crossprod[2]);
	n0 = sqrt(vectors[0*3+0]*vectors[0*3+0] + vectors[0*3+1]*vectors[0*3+1] + vectors[0*3+2]*vectors[0*3+2]);
	n1 = sqrt(vectors[1*3+0]*vectors[1*3+0] + vectors[1*3+1]*vectors[1*3+1] + vectors[1*3+2]*vectors[1*3+2]);
	n2 = sqrt(vectors[2*3+0]*vectors[2*3+0] + vectors[2*3+1]*vectors[2*3+1] + vectors[2*3+2]*vectors[2*3+2]);
	//scalarprod = |c|.|a|.|b| sin(theta)*cos(alpha)
	scalarprod = crossprod[0]*vectors[2*3+0] + crossprod[0]*vectors[2*3+1] + crossprod[0]*vectors[2*3+2];
	cossin[1] = ncross/(n0*n1);
	cossin[0] = scalarprod/(n0*n1*n2*cossin[1]); 
}

double compute_distance3(double posA[],double posB[])
{
	double val;
	val = sqrt(  (posA[0]-posB[0])*(posA[0]-posB[0])  +  (posA[1]-posB[1])*(posA[1]-posB[1])  +  (posA[2]-posB[2])*(posA[2]-posB[2])  );
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
	PetscInt        *gidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          dx,dy,dz,dl_min,dl_max,el_ar,max_ar,max_ar_g;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;

	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dm,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dm,0,&gidx);CHKERRQ(ierr);
	
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

	ierr = MPI_Allreduce(&max_ar,&max_ar_g,1,MPI_DOUBLE,MPI_MAX,((PetscObject)dm)->comm);CHKERRQ(ierr);
											 
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
	PetscInt        *gidx;
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
	ierr = DMDAGetCoordinateDA(dm,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dm,0,&gidx);CHKERRQ(ierr);
	
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
	
	ierr = MPI_Allreduce(&min_dist,&min_dist_g,1,MPI_DOUBLE,MPI_MIN,((PetscObject)dm)->comm);CHKERRQ(ierr);
	
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
	PetscInt        *gidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          posA[3],posB[3];
	double          diag,dl_min,dl_min_g;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dm,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dm,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	
	dl_min = 1.0e32;
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		
		
		get_node_coordinate(el_coords,0,posA);
		get_node_coordinate(el_coords,26,posB);
		diag =	compute_distance3(posA,posB);
		if (diag < dl_min) { dl_min = diag; }

		get_node_coordinate(el_coords,2,posA);
		get_node_coordinate(el_coords,24,posB);
		diag =	compute_distance3(posA,posB);
		if (diag < dl_min) { dl_min = diag; }

		get_node_coordinate(el_coords,6,posA);
		get_node_coordinate(el_coords,20,posB);
		diag =	compute_distance3(posA,posB);
		if (diag < dl_min) { dl_min = diag; }
		
		get_node_coordinate(el_coords,8,posA);
		get_node_coordinate(el_coords,18,posB);
		diag =	compute_distance3(posA,posB);
		if (diag < dl_min) { dl_min = diag; }
				
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&dl_min,&dl_min_g,1,MPI_DOUBLE,MPI_MIN,((PetscObject)dm)->comm);CHKERRQ(ierr);
	
	*value = (PetscReal)dl_min_g;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetric_AreaRatio"
PetscErrorCode DMDAComputeMeshQualityMetric_FaceAreaRatio(DM dm,PetscReal *value)
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e;
	const PetscInt  *el_nidx;
	PetscInt        *gidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          a1,a2,corners[3*4];
	double          ratio,ratio_min,ratio_min_g;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dm,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dm,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	

	ratio_min = 1.0e32;
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
        //Loops through face pairs.
		get_face_coordinates(el_coords, Q2_FACE_WEST, corners);
        a1 = compute_quadrilateral_area_(corners);
        get_face_coordinates(el_coords, Q2_FACE_EAST, corners);
		a2 = compute_quadrilateral_area_(corners);
        ratio = (a1<a2)?a1/a2:a2/a1;
        ratio_min = (ratio<ratio_min)?ratio:ratio_min;

		get_face_coordinates(el_coords, Q2_FACE_SOUTH, corners);
        a1 = compute_quadrilateral_area_(corners);
        get_face_coordinates(el_coords, Q2_FACE_NORTH, corners);
		a2 = compute_quadrilateral_area_(corners);
        ratio = (a1<a2)?a1/a2:a2/a1;
        ratio_min = (ratio<ratio_min)?ratio:ratio_min;        
        
		get_face_coordinates(el_coords, Q2_FACE_BACK, corners);
        a1 = compute_quadrilateral_area_(corners);
        get_face_coordinates(el_coords, Q2_FACE_FRONT, corners);
		a2 = compute_quadrilateral_area_(corners);
        ratio = (a1<a2)?a1/a2:a2/a1;
        ratio_min = (ratio<ratio_min)?ratio:ratio_min;            
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&ratio_min,&ratio_min_g,1,MPI_DOUBLE,MPI_MIN,((PetscObject)dm)->comm);CHKERRQ(ierr);
	
	*value = (PetscReal)ratio_min_g;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshQualityMetric_CornerAngle"
PetscErrorCode DMDAComputeMeshQualityMetric_CornerAngle(DM dm,PetscReal *value)
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e;
	const PetscInt  *el_nidx;
	PetscInt        *gidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	double          val, val_min, val_min_g, cossin[2], vectors[3*3];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;

	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dm,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dm,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	
	val_min = 2.0; 
	for (e=0;e<nel;e++) {
	    ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
        //Loops through corners
        get_corner_vectors(el_coords,Q2_VERTEX_0,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val;
        
        get_corner_vectors(el_coords,Q2_VERTEX_1,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val;
        
        get_corner_vectors(el_coords,Q2_VERTEX_2,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_corner_vectors(el_coords,Q2_VERTEX_3,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_corner_vectors(el_coords,Q2_VERTEX_4,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_corner_vectors(el_coords,Q2_VERTEX_5,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_corner_vectors(el_coords,Q2_VERTEX_6,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
        
        get_corner_vectors(el_coords,Q2_VERTEX_7,vectors);
        compute_cossin(vectors,cossin); 
        val = (cossin[0]<cossin[1])?cossin[0]:cossin[1];
        val_min = (val_min<val)?val_min:val; 
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
    
	ierr = MPI_Allreduce(&val_min,&val_min_g,1,MPI_DOUBLE,MPI_MIN,((PetscObject)dm)->comm);CHKERRQ(ierr);
    
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
        
		case MESH_QUALITY_CORNER_ANGLE:
			ierr = DMDAComputeMeshQualityMetric_CornerAngle(dm,value);CHKERRQ(ierr);
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



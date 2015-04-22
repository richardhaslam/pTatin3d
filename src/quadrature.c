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
 **    filename:   quadrature.c
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

#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/quadrature_impl.h"
#include "data_bucket.h"
#include "dmda_element_q2p1.h"
#include "quadrature.h"

#undef __FUNCT__  
#define __FUNCT__ "QuadratureCreate"
PetscErrorCode QuadratureCreate(Quadrature *quadrature)
{
	Quadrature Q;
	PetscErrorCode  ierr;

	PetscFunctionBegin;

	ierr = PetscMalloc(sizeof(struct _p_Quadrature),&Q);CHKERRQ(ierr);
	ierr = PetscMemzero(Q,sizeof(struct _p_Quadrature));CHKERRQ(ierr);

	*quadrature = Q;

	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QuadratureDestroy"
PetscErrorCode QuadratureDestroy(Quadrature *quadrature)
{
	Quadrature Q;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;

	if (!quadrature) { PetscFunctionReturn(0); }
	Q = *quadrature;
	
	if (Q->q_xi_coor) { ierr = PetscFree(Q->q_xi_coor);CHKERRQ(ierr); }
	if (Q->q_weight) { ierr = PetscFree(Q->q_weight);CHKERRQ(ierr); }
	if (Q->properties_db) { DataBucketDestroy(&Q->properties_db); }
	
	ierr = PetscFree(Q);CHKERRQ(ierr);
	*quadrature = NULL;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QuadratureView"
PetscErrorCode QuadratureView(Quadrature q)
{
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"QuadratureView:\n");
	PetscPrintf(PETSC_COMM_WORLD,"  dim    %D\n", q->dim);
	PetscPrintf(PETSC_COMM_WORLD,"  type    %D\n", (PetscInt)q->type);
	PetscPrintf(PETSC_COMM_WORLD,"  npoints    %D\n", q->npoints);
	PetscPrintf(PETSC_COMM_WORLD,"  n_elements %D\n", q->n_elements);
	
	DataBucketView(PETSC_COMM_WORLD, q->properties_db,"GaussLegendre StokesCoefficients",DATABUCKET_VIEW_STDOUT);

	PetscFunctionReturn(0);
}

void QuadratureCreateGauss_2pnt_3D(PetscInt *ngp,PetscReal **_q_coor,PetscReal **_q_weight)
{
	const double s = 0.577350269189;
	const double w_1d[] = { 1.0, 1.0 };
	const double xi_1d[] = { -s, s};
	int nI,nJ,nK;
	PetscReal *q_coor,*q_weight;
	
	
	/* standard 2x2x2 point quadrature */
	*ngp = 8;
	PetscMalloc( sizeof(double)*(*ngp)*3, &q_coor );
	PetscMalloc( sizeof(double)*(*ngp)  , &q_weight );
	
	for( nI=0; nI<2; nI++ ) {
		for( nJ=0; nJ<2; nJ++ ) {
			for( nK=0; nK<2; nK++ ) {
				int idx = nI + nJ*2 + nK*2*2;
				
				q_weight[idx] = w_1d[nI] * w_1d[nJ] * w_1d[nK];
				
				q_coor[3*idx+0] = xi_1d[nI];
				q_coor[3*idx+1] = xi_1d[nJ];
				q_coor[3*idx+2] = xi_1d[nK];
			}
		}
	}
	*_q_coor = q_coor;
	*_q_weight = q_weight;
}

void QuadratureCreateGauss_3pnt_3D(PetscInt *ngp,PetscReal **_q_coor,PetscReal **_q_weight)
{
	const double sqrt_15_on_5 = 0.774596669241483; /* sqrt(15)/5 */
	const double five_on_9 = 0.555555555555556;
	const double eight_on_9 = 0.888888888888889;
	const double w_1d[] = { five_on_9, eight_on_9, five_on_9 };
	const double xi_1d[] = { -sqrt_15_on_5, 0.0, sqrt_15_on_5 };
	int nI,nJ,nK;
	PetscReal *q_coor,*q_weight;
	
	
	/* standard 3x3x3 point quadrature */
	*ngp = 27;
	PetscMalloc( sizeof(double)*(*ngp)*3, &q_coor );
	PetscMalloc( sizeof(double)*(*ngp)  , &q_weight );
	
	for( nI=0; nI<3; nI++ ) {
		for( nJ=0; nJ<3; nJ++ ) {
			for( nK=0; nK<3; nK++ ) {
				int idx = nI + nJ*3 + nK*3*3;
				
				q_weight[idx] = w_1d[nI] * w_1d[nJ] * w_1d[nK];
				
				q_coor[3*idx+0] = xi_1d[nI];
				q_coor[3*idx+1] = xi_1d[nJ];
				q_coor[3*idx+2] = xi_1d[nK];
			}
		}
	}
	*_q_coor = q_coor;
	*_q_weight = q_weight;
}

/* surface quadrature */
#undef __FUNCT__  
#define __FUNCT__ "SurfaceQuadratureCreate"
PetscErrorCode SurfaceQuadratureCreate(SurfaceQuadrature *quadrature)
{
	SurfaceQuadrature Q;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ierr = PetscMalloc(sizeof(struct _p_SurfaceQuadrature),&Q);CHKERRQ(ierr);
	ierr = PetscMemzero(Q,sizeof(struct _p_SurfaceQuadrature));CHKERRQ(ierr);
	
	*quadrature = Q;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_SurfaceQuadratureCreate"
PetscErrorCode _SurfaceQuadratureCreate(SurfaceQuadrature quadrature,HexElementFace index,PetscInt nfaces)
{
	ConformingElementFamily e;
	int ngp32;
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
	//PetscPrintf(PETSC_COMM_WORLD,"SurfaceQuadratureCreate:\n");
	
	ElementTypeCreate_Q2(&e,3);
	quadrature->e       = e;
	quadrature->face_id = index;
	e->generate_surface_quadrature_3D(e,index,&ngp32,quadrature->gp2,quadrature->gp3);	
	quadrature->ngp = (PetscInt)ngp32;
	
	quadrature->nfaces = nfaces;
	
	//PetscPrintf(PETSC_COMM_WORLD,"\t[SurfaceQuadrature]: attributing %D edge elements \n", nfaces );
	//PetscPrintf(PETSC_COMM_WORLD,"\t[SurfaceQPointCoefficient]: attributing %D surface quadrature points \n", nfaces*quadrature->ngp );
	if (nfaces != 0) {
		ierr = PetscMalloc( sizeof(PetscInt)*nfaces, &quadrature->element_list);CHKERRQ(ierr);
	} else {
		quadrature->element_list = NULL;
	}
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_SurfaceQuadratureCellIndexSetUp"
PetscErrorCode _SurfaceQuadratureCellIndexSetUp(SurfaceQuadrature Q,HexElementFace index,PetscInt nface_edge,DM da)
{
	PetscInt eli,elj,elk;
	PetscInt si,sj,sk,ni,nj,nk,M,N,P,lmx,lmy,lmz;
	PetscInt cnt,elidx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(da,0,&M,&N,&P, 0,0,0,0, 0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
	ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,&lmz);CHKERRQ(ierr);

	
	if (nface_edge == 0) { 
		PetscFunctionReturn(0); 
	}
  
	switch (index) {
			
		case HEX_FACE_Pxi:
			cnt = 0;
			eli = lmx - 1;
			for (elk=0; elk<lmz; elk++) {
				for (elj=0; elj<lmy; elj++) {
					elidx = eli + elj * lmx + elk * lmx*lmy;
					Q->element_list[ cnt ] = elidx;
					cnt++;
				}
			}
			break;
			
		case HEX_FACE_Nxi:
			cnt = 0;
			eli = 0;
			for (elk=0; elk<lmz; elk++) {
				for (elj=0; elj<lmy; elj++) {
					elidx = eli + elj * lmx + elk * lmx*lmy;
					Q->element_list[ cnt ] = elidx;
					cnt++;
				}
			}
			break;
			
		case HEX_FACE_Peta:
			cnt = 0;
			elj = lmy - 1;
			for (elk=0; elk<lmz; elk++) {
				for (eli=0; eli<lmx; eli++) {
					elidx = eli + elj * lmx + elk * lmx*lmy;
					Q->element_list[ cnt ] = elidx;
					cnt++;
				}
			}
			break;
			
		case HEX_FACE_Neta:
			cnt = 0;
			elj = 0;
			for (elk=0; elk<lmz; elk++) {
				for (eli=0; eli<lmx; eli++) {
					elidx = eli + elj * lmx + elk * lmx*lmy;
					Q->element_list[ cnt ] = elidx;
					cnt++;
				}
			}
			break;
			
		case HEX_FACE_Pzeta:
			cnt = 0;
			elk = lmz - 1;
			for (elj=0; elj<lmy; elj++) {
				for (eli=0; eli<lmx; eli++) {
					elidx = eli + elj * lmx + elk * lmx*lmy;
					Q->element_list[ cnt ] = elidx;
					cnt++;
				}
			}
			break;
			
			
		case HEX_FACE_Nzeta:
			cnt = 0;
			elk = 0;
			for (elj=0; elj<lmy; elj++) {
				for (eli=0; eli<lmx; eli++) {
					elidx = eli + elj * lmx + elk * lmx*lmy;
					Q->element_list[ cnt ] = elidx;
					cnt++;
				}
			}
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SurfaceQuadratureDestroy"
PetscErrorCode SurfaceQuadratureDestroy(SurfaceQuadrature *quadrature)
{
	SurfaceQuadrature Q;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	if (!quadrature) { PetscFunctionReturn(0); }
	
	Q = *quadrature;
	
	if (Q->element_list) { ierr = PetscFree(Q->element_list);CHKERRQ(ierr); }
	if (Q->properties_db) { DataBucketDestroy(&Q->properties_db); }
	ElementTypeDestroy_Q2(&Q->e);

	ierr = PetscFree(Q);CHKERRQ(ierr);
	
	*quadrature = NULL;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SurfaceQuadratureGetElementFamily"
PetscErrorCode SurfaceQuadratureGetElementFamily(SurfaceQuadrature q,ConformingElementFamily *e)
{
	if (e) { *e = q->e; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SurfaceQuadratureGetQuadratureInfo"
PetscErrorCode SurfaceQuadratureGetQuadratureInfo(SurfaceQuadrature q,PetscInt *nqp,QPoint2d **qp2,QPoint3d **qp3)
{
	if (nqp) { *nqp = q->ngp; }
	if (qp2) { *qp2 = q->gp2; }
	if (qp3) { *qp3 = q->gp3; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SurfaceQuadratureGetFaceInfo"
PetscErrorCode SurfaceQuadratureGetFaceInfo(SurfaceQuadrature q,HexElementFace *faceid,PetscInt *nfaces,PetscInt **ellist)
{
	if (faceid) { *faceid = q->face_id; }
	if (nfaces) { *nfaces = q->nfaces; }
	if (ellist) { *ellist = q->element_list; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SurfaceQuadratureInterpolate3D"
PetscErrorCode SurfaceQuadratureInterpolate3D(SurfaceQuadrature q,QPoint3d *qp3d,PetscInt ndof,PetscReal field[],PetscReal value[])
{
    int    k,d;
    double Ni[27];
    
    q->e->basis_NI_3D(qp3d,Ni);

    for (d=0; d<ndof; d++) {
        value[d] = 0.0;
        for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
            value[d] += Ni[k] * field[ndof*k + d];
        }
    }
    
	PetscFunctionReturn(0);
}

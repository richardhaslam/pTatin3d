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
 **    Filename:      mp_advection.c
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
 *  mp_advection.c
 *  
 *
 *  Created by Dave May on 3/4/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */


#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "petscdm.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "swarm_fields.h"
#include "data_exchanger.h"
#include "MPntStd_def.h"
#include "ptatin3d_stokes.h"
#include "element_utils_q2.h"
#include "dmda_element_q2p1.h"



#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_MPntStd_Euler"
PetscErrorCode SwarmUpdatePosition_MPntStd_Euler(DM da,Vec velocity,PetscReal step,int npoints,MPntStd marker[])
{
	Vec             Lvelocity;
	PetscScalar     *LA_velocity;
	PetscScalar     el_velocity[Q2_NODES_PER_EL_3D*NSD];
	PetscInt        e,i,p,wil;
	PetscScalar     Ni_p[Q2_NODES_PER_EL_3D],vel_p[NSD];
	PetscInt nel,nen_u;
	const PetscInt *elnidx_u;
	PetscInt vel_el_lidx[U_BASIS_FUNCTIONS*3];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* advect */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p = &marker[p];
		
		wil   = marker_p->wil;
		e     = wil;
		if (wil<0) { SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Point[%d] has wil_e < 0", wil ); }
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = DMDAGetVectorElementFieldQ2_3D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		P3D_ConstructNi_Q2_3D(marker_p->xi,Ni_p);
		
		vel_p[0] = vel_p[1] = vel_p[2] = 0.0;
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			vel_p[0] += Ni_p[i] * el_velocity[NSD*i+0];
			vel_p[1] += Ni_p[i] * el_velocity[NSD*i+1];
			vel_p[2] += Ni_p[i] * el_velocity[NSD*i+2];
		}
		
		marker_p->coor[0] = marker_p->coor[0] + step * vel_p[0];
		marker_p->coor[1] = marker_p->coor[1] + step * vel_p[1];
		marker_p->coor[2] = marker_p->coor[2] + step * vel_p[2];
	}
	
  ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

void LSFDeriv3dQ2(double _xi[],double **GNi)
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

void LSF3dQ2_CheckPartitionOfUnity(double xi[],double *val)
{
	const double tol = 1.0e-6;
	double Ni[27],sum;
	int i;
	
	P3D_ConstructNi_Q2_3D(xi,Ni);
	sum = 0.0;
	for (i=0; i<27; i++) {
		sum += Ni[i];
	}
	*val = sum;
	if (fabs(sum-1.0) > tol) {
		printf("**** sum( N_i(xi,eta) ) = %1.8e > 1.0 ==> partition of unity is not satisified, point is likely outside of the element ****\n",sum);
	}
}
void LSF3dQ2_CheckGlobalCoordinate(double element_coord[],double xi[],double xp[],double err[])
{
	const double tol = 1.0e-6;
	double Ni[27],xp_interp[3];
	int i;
	
	P3D_ConstructNi_Q2_3D(xi,Ni);
	xp_interp[0] = 0.0;
	xp_interp[1] = 0.0;
	xp_interp[2] = 0.0;
	
	for (i=0; i<27; i++) {
		xp_interp[0] += Ni[i] * element_coord[3*i + 0];
		xp_interp[1] += Ni[i] * element_coord[3*i + 1];
		xp_interp[2] += Ni[i] * element_coord[3*i + 2];
	}
	err[0] = (xp[0] - xp_interp[0]);
	err[1] = (xp[1] - xp_interp[1]);
	err[2] = (xp[2] - xp_interp[2]);
	
	if (fabs(err[0]) > tol) {
		printf("**** |xp - N_i(xi,eta).x_i| = %1.8e > %1.8e ==> x: coordinate interpolation error occurred ****\n",fabs(err[0]),tol);
	}
	if (fabs(err[1]) > tol) {
		printf("**** |yp - N_i(xi,eta).y_i| = %1.8e > %1.8e ==> y: coordinate interpolation error occurred ****\n",fabs(err[1]),tol);
	}
	if (fabs(err[2]) > tol) {
		printf("**** |zp - N_i(xi,eta).z_i| = %1.8e > %1.8e ==> z: coordinate interpolation error occurred ****\n",fabs(err[2]),tol);
	}
}

/* computes h = inv(J) f */
void _compute_deltaX( double A[3][3], double f[], double h[] )
{
	double B[3][3];
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
	
	h[0] = B[0][0]*f[0] + B[0][1]*f[1] + B[0][2]*f[2];
	h[1] = B[1][0]*f[0] + B[1][1]*f[1] + B[1][2]*f[2];
	h[2] = B[2][0]*f[0] + B[2][1]*f[1] + B[2][2]*f[2];
}

void _compute_J_3dQ2(double xi[],double vertex[],double J[3][3])
{
	int i,j;
	double GNi[3][27];
	
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			J[i][j] = 0.0;
		}
	}
	
	P3D_ConstructGNi_Q2_3D(xi,GNi);
	for (i=0; i<27; i++) {
		int i3 = i*3;
		double x = vertex[i3];
		double y = vertex[i3+1];
		double z = vertex[i3+2];
		
		J[0][0] += x * GNi[0][i];
		J[0][1] += x * GNi[1][i];
		J[0][2] += x * GNi[2][i];
		
		J[1][0] += y * GNi[0][i];
		J[1][1] += y * GNi[1][i];
		J[1][2] += y * GNi[2][i];
		
		J[2][0] += z * GNi[0][i];
		J[2][1] += z * GNi[1][i];
		J[2][2] += z * GNi[2][i];		
	}
}

void _compute_F_3dQ2(double xi[],double vertex[],double pos[],double f[])
{
	int i;
	double Ni[27];
	
	/* Update F for the next iteration */
	f[0] = f[1] = f[2] = 0.0;
	
	P3D_ConstructNi_Q2_3D(xi,Ni);
	for (i=0; i<27; i++) {
		int i3   = i*3;
		int i3p1 = i3+1;
		int i3p2 = i3+2;
		
		f[0] += vertex[i3  ] * Ni[i];
		f[1] += vertex[i3p1] * Ni[i];
		f[2] += vertex[i3p2] * Ni[i];
	}
	f[0] = - f[0] + pos[0];
	f[1] = - f[1] + pos[1];
	f[2] = - f[2] + pos[2];
}

void InverseMappingDomain_3dQ2( 
															 double tolerance, int max_its,
															 PetscBool use_nonzero_guess, 
															 PetscBool monitor, PetscBool log,
															 const double coords[], const int mx, const int my, const int mz,const int element[],
															 int np, MPntStd marker[] )
{
	const int nodesPerEl = Q2_NODES_PER_EL_3D; 
	double  h[NSD];
	double  Jacobian[NSD][NSD];
	double  f[NSD];
	int     i,p;
	int     its;
	double  residual2,tolerance2,F2;
	
	double    cxip[3],Lxip[3],Gxip[3];
	double    dxi,deta,dzeta,xi0,eta0,zeta0;
	int       I,J,K,wil_IJ,wil_2d,eid,k;
	double    vertex[NSD * Q2_NODES_PER_EL_3D];
	int       n0,n1,n2,n3;
	PetscBool point_found;
	
	tolerance2 = tolerance * tolerance; /* Eliminates the need to do a sqrt in the convergence test */
	
	if(log)printf("Domain: ncells = %d x %d x %d = %d \n", mx,my,mz,mx*my*mz );
	
	/* map domain to [-1,1]x[-1,1]x[-1,1] domain */
	dxi   = 2.0/((double)mx);
	deta  = 2.0/((double)my);
	dzeta = 2.0/((double)mz);
	if(log)printf("Domain: (dxi,eta,zeta) = (%1.8e,%1.8e,%1.8e)\n",dxi,deta,dzeta );
	
	for( p=0; p<np; p++ ) {
		MPntStd *marker_p = &marker[p];
		
		if(log)printf("POINT[%d]\n", p );

		/* copy these values */
		cxip[0] = marker_p->xi[0];
		cxip[1] = marker_p->xi[1];
		cxip[2] = marker_p->xi[2];
		
		/* Check for an initial guess initial guess */
		if( use_nonzero_guess == PETSC_FALSE ) {
			Gxip[0] = 0.0;
			Gxip[1] = 0.0;
			Gxip[2] = 0.0;
		}
		else {
			/* convert wil => IJ */
			wil_IJ = marker_p->wil;
			K = wil_IJ/(mx*my);
			wil_2d = wil_IJ - K*mx*my;
			J = wil_2d/mx;
			I = wil_2d - J*mx;
			if(log)printf("init I,J,K = %d %d %d [wil=%d]\n", I,J,K,wil_IJ );
			/* convert Lxip => Gxip */
			xi0   = -1.0 + I*dxi;
			eta0  = -1.0 + J*deta;
			zeta0 = -1.0 + K*dzeta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Gxip[0] = 2.0 * (cxip[0]-xi0)/dxi - 1.0;
			//			Gxip[1] = 2.0 * (cxip[1]-eta0)/deta - 1.0;
			// x*-x*0/dx = (x+1)/2
			Gxip[0] = dxi   * (cxip[0]+1.0)/2.0 + xi0;
			Gxip[1] = deta  * (cxip[1]+1.0)/2.0 + eta0;
			Gxip[2] = dzeta * (cxip[2]+1.0)/2.0 + zeta0;
			if(log)printf("[Lxi-init] = %1.8e %1.8e %1.8e\n", cxip[0], cxip[1], cxip[2] );
			if(log)printf("[Gxi-init] = %1.8e %1.8e %1.8e\n", Gxip[0], Gxip[1], Gxip[2] );
			
			/* check */
			if (log) {
				double err[3];

				for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
					int nid = element[wil_IJ*Q2_NODES_PER_EL_3D+k];
					
					vertex[3*k+0] = coords[3*nid+0];
					vertex[3*k+1] = coords[3*nid+1];
					vertex[3*k+2] = coords[3*nid+2];
				}
				LSF3dQ2_CheckGlobalCoordinate(vertex,cxip,marker_p->coor,err);
				printf("  interpolated coord err %1.4e %1.4e %1.4e \n",err[0],err[1],err[2]);
			}
			
		}
		if(monitor)printf("point[%d]: pos = ( %+1.8e, %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e, %+1.8e ) \n", p, marker_p->coor[0],marker_p->coor[1],marker_p->coor[2], marker_p->xi[0],marker_p->xi[1],marker_p->xi[2] );
		
		point_found = PETSC_FALSE;
		
		its = 0;
		do {
			if(log)printf("iteration: %d\n",its);
			/* convert Gxi to IJ */
			I = (Gxip[0]+1.0)/dxi;
			J = (Gxip[1]+1.0)/deta;
			K = (Gxip[2]+1.0)/dzeta;
			
			if( I==mx ) I--;
			if( J==my ) J--;
			if( K==mz ) K--;
			
			if( (I<0) || (J<0)|| (K<0) ) {
				if(log)printf("  I(%d),J(%d),K(%d) negative Gxip %1.8e,%1.8e,%1.8e \n",I,J,K,Gxip[0],Gxip[1],Gxip[2]);
				break;
			}
			if( I>=mx ) { 
				if(log)printf("  I too large \n");
				break;
			}
			if( J>=my ) {
				if(log)printf("  J too large \n");
				break;
			}
			if( K>=mz ) {
				if(log)printf("  K too large \n");
				break;
			}
			
			
			/* Get coords of cell IJ */
			wil_IJ = I + J*mx + K*mx*my;
			if(log)printf("  I,J,K=%d/%d/%d : wil_IJ %d : nid = ", I,J,K,wil_IJ);
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
				int nid = element[wil_IJ*Q2_NODES_PER_EL_3D+k];
				
				vertex[3*k+0] = coords[3*nid+0];
				vertex[3*k+1] = coords[3*nid+1];
				vertex[3*k+2] = coords[3*nid+2];
				if(log)printf("%d ", nid);
			}
			if(log)printf("\n");
			
			
			if(log) {
				printf("  [vertex] ");
				for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
					printf("(%1.8e , %1.8e , %1.8e) ",vertex[3*k+0],vertex[3*k+1],vertex[3*k+2] );
				}
				printf("\n");
			}
			
			/* convert global (domain) xi TO local (element) xi  */
			xi0   = -1.0 + I*dxi;
			eta0  = -1.0 + J*deta;
			zeta0 = -1.0 + K*dzeta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Lxip[0] = 0.5*(Gxip[0]+1.0)*dxi  + xi0;
			//			Lxip[1] = 0.5*(Gxip[1]+1.0)*deta + eta0;
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0  )/dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0 )/deta  - 1.0;
			Lxip[2] = 2.0*(Gxip[2]-zeta0)/dzeta - 1.0;
			
			if(log)printf("  Lxi,Lxeta,Lzeta = %1.8e, %1.8e, %1.8e (%d,%d,%d) \n", Lxip[0],Lxip[1],Lxip[2],I,J,K );
			
			_compute_F_3dQ2( Lxip, vertex, marker_p->coor, f );
			if( monitor ) printf("%4d InverseMapping : F = ( %+1.8e, %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e, %+1.8e ) \n", its, f[0],f[1],f[2], Lxip[0],Lxip[1],Lxip[2] );
			
			/* Check for convergence */
			F2 = (f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
			if( F2 < tolerance2 ) {
				if( monitor ) printf("%4d InverseMapping : converged : Norm of F %1.8e \n", its, sqrt(F2) );
				point_found = PETSC_TRUE;
				break;
			}
			
			_compute_J_3dQ2( Lxip, vertex, Jacobian );
			
			/* compute update */
			_compute_deltaX( Jacobian, f, h );
			if(log)printf("  [delta] = %1.8e %1.8e %1.8e \n", h[0],h[1],h[2] );
			
			/* update Lxip */
			Lxip[0] += 10.0e-1 *h[0];
			Lxip[1] += 10.0e-1 *h[1];
			Lxip[2] += 10.0e-1 *h[2];
			if(log)printf("  [corrected] Lxi,Leta,Lzeta = %1.8e, %1.8e, %1.8e \n", Lxip[0],Lxip[1],Lxip[2] );
			
			residual2 = ( h[0]*h[0] + h[1]*h[1] + h[2]*h[2] );
			if( residual2 < tolerance2 ) {
				if( monitor ) printf("%4d InverseMapping : converged : Norm of correction %1.8e \n", its, sqrt(residual2) );
				point_found = PETSC_TRUE;
				break;
			}
			
			/* convert Lxip => Gxip */
			//			xi0  = -1.0 + I*dxi;
			//			eta0 = -1.0 + J*deta;
			
			//			Gxip[0] = 2.0 * (Lxip[0]-xi0)/dxi   - 1.0;
			//			Gxip[1] = 2.0 * (Lxip[1]-eta0)/deta - 1.0;
			// x*-x*0/dx = (x+1)/2
			Gxip[0] = dxi   * (Lxip[0]+1.0)/2.0 + xi0;
			Gxip[1] = deta  * (Lxip[1]+1.0)/2.0 + eta0;
			Gxip[2] = dzeta * (Lxip[2]+1.0)/2.0 + zeta0;
			if(log)printf("  [Gxi] = %1.8e %1.8e %1.8e\n", Gxip[0], Gxip[1], Gxip[2] );
			
			if (Gxip[0]<-1.0) { 
				Gxip[0] = -1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			if (Gxip[1]<-1.0) {
				Gxip[1] = -1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			if (Gxip[2]<-1.0) {
				Gxip[2] = -1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			
			if (Gxip[0]>1.0) {
				Gxip[0] = 1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			if (Gxip[1]>1.0) {
				Gxip[1] = 1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			if (Gxip[2]>1.0) {
				Gxip[2] = 1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			
			its++;
		} while(its<max_its);
		
		if( monitor && point_found==PETSC_FALSE ){
			if( its>=max_its ) {
				printf("%4d %s : Reached maximum iterations (%d) without converging. \n", its, __FUNCTION__, max_its );
			}
			else {
				printf("%4d %s : Newton broke down, diverged or stagnated after (%d) iterations without converging. \n", its, __FUNCTION__, its );
			}
		}
		
		/* if at the end of the solve, it still looks like the point is outside the mapped domain, mark point as not being found */
		if( fabs(Gxip[0]) > 1.0 ) { point_found =PETSC_FALSE; }
		if( fabs(Gxip[1]) > 1.0 ) { point_found =PETSC_FALSE; }
		if( fabs(Gxip[2]) > 1.0 ) { point_found =PETSC_FALSE; }
		
		/* update local variables */
		if( point_found==PETSC_FALSE ) {
			Lxip[0] = NAN;
			Lxip[1] = NAN;
			Lxip[2] = NAN;
			wil_IJ  = -1;
		}
		else {
			/* convert Gxi to IJ */
			I = (Gxip[0]+1.0)/dxi;
			J = (Gxip[1]+1.0)/deta;
			K = (Gxip[2]+1.0)/dzeta;
			if( I==mx ) I--;
			if( J==my ) J--;
			if( K==mz ) K--;
			
			if( I>=mx ) {
				if(log)printf("  I too large \n");
				break;
			}
			if( J>=my ) {
				if(log)printf("  J too large \n");
				break;
			}
			if( K>=mz ) {
				if(log)printf("  K too large \n");
				break;
			}
			
			/* convert global (domain) xi TO local (element) xi  */
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			xi0   = -1.0 + I*dxi;
			eta0  = -1.0 + J*deta;
			zeta0 = -1.0 + K*dzeta;
			
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0  )/dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0 )/deta  - 1.0;
			Lxip[2] = 2.0*(Gxip[2]-zeta0)/dzeta - 1.0;
			
			wil_IJ = I + J*mx + K*mx*my;
		}
		
		/* set into vector */
		marker_p->xi[0] = Lxip[0];
		marker_p->xi[1] = Lxip[1];
		marker_p->xi[2] = Lxip[2];
		marker_p->wil   = wil_IJ;

		if(log)printf("  <<final>> xi,eta,zeta = %1.8e, %1.8e, %1.8e [wil=%d] \n", Lxip[0],Lxip[1],Lxip[2],wil_IJ);

	}
}



#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_ComputeCourantStep"
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep(DM da,Vec velocity,PetscReal *step)
{
	Vec             Lvelocity, gcoords;
	PetscScalar     *LA_velocity, *LA_coords;
	PetscScalar     el_coords[Q2_NODES_PER_EL_3D*NSD];
	PetscScalar     el_velocity[Q2_NODES_PER_EL_3D*NSD];
	PetscInt        e,i;
	PetscScalar     Ni_p[Q2_NODES_PER_EL_3D],xi_p[NSD];
	PetscInt        nel,nen_u,ii;
	const PetscInt  *elnidx_u;
	DM              cda;
	double          dt_min_local, dt_min;
	MPI_Comm        comm;
	PetscReal       hx,hy,hz,vavg[NSD],dtx,dty,dtz;
	PetscReal       xmin,xmax,ymin,ymax,zmin,zmax,coor;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
  /* setup for coords */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	/* setup velocity */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	dt_min_local = 1.0e32;
	for (e=0; e<nel; e++) {
		double vc[NSD];
		
		/* get coords for the element */
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx_u[nen_u*e],LA_coords);CHKERRQ(ierr);
		/* get velocity */
		ierr = DMDAGetVectorElementFieldQ2_3D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		/* find min,max x,y,z */
		xmin = ymin = zmin =  1.0e32;
		xmax = ymax = zmax = -1.0e32;
		vavg[0] = vavg[1] = vavg[2] = 0.0;
		for (ii=0; ii<Q2_NODES_PER_EL_3D; ii++) {
			coor = el_coords[3*ii+0];
			if (coor < xmin) { xmin = coor; }
			if (coor > xmax) { xmax = coor; }
			
			coor = el_coords[3*ii+1];
			if (coor < ymin) { ymin = coor; }
			if (coor > ymax) { ymax = coor; }
			
			coor = el_coords[3*ii+2];
			if (coor < zmin) { zmin = coor; }
			if (coor > zmax) { zmax = coor; }
			
			vavg[0] = vavg[0] + el_velocity[3*ii+0];
			vavg[1] = vavg[1] + el_velocity[3*ii+1];
			vavg[2] = vavg[2] + el_velocity[3*ii+2];
		}
		vavg[0] = vavg[0]/(PetscReal)(Q2_NODES_PER_EL_3D);
		vavg[1] = vavg[1]/(PetscReal)(Q2_NODES_PER_EL_3D);
		vavg[2] = vavg[2]/(PetscReal)(Q2_NODES_PER_EL_3D);
		
		hx = xmax - xmin;
		hy = ymax - ymin;
		hz = zmax - zmin;
		
		dtx = fabs( hx / vavg[0] );
		dty = fabs( hy / vavg[1] );
		dtz = fabs( hz / vavg[2] );
		
		if (dtx < dt_min_local) { dt_min_local = dtx; }
		if (dty < dt_min_local) { dt_min_local = dty; }
		if (dtz < dt_min_local) { dt_min_local = dtz; }
		
	}
	
  ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	
	ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&dt_min_local,&dt_min,1,MPI_DOUBLE,MPI_MIN,comm);CHKERRQ(ierr);
	
	{
		PetscScalar min,max;
		PetscInt lmin,lmax;
		ierr = VecMin(velocity,&lmin,&min);CHKERRQ(ierr);
		ierr = VecMax(velocity,&lmax,&max);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  vel: min %lf(%d): max %lf(%d)\n", min,lmin,max,lmax);
	}
	
	*step = dt_min;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdateProperties_MPntStd"
PetscErrorCode SwarmUpdateProperties_MPntStd(DataBucket db,pTatinCtx ctx,Vec X)
{
	BTruth         found;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	DataBucketQueryDataFieldByName(db,MPntStd_classname,&found);
	if(found==BFALSE) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cannot find DataField with name %s \n", MPntStd_classname );
	}
	
	
	PetscFunctionReturn(0);
}



/* ADVECT MARKERS */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateGlobalCoordinates"
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates(DataBucket materialpoints,DM dav,Vec velocity,PetscReal dt)
{
	PetscErrorCode ierr;
	int            npoints;
	MPntStd        *mp_std;
	DataField      PField;
	
	PetscFunctionBegin;
	DataBucketGetSizes(materialpoints,&npoints,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
	ierr = SwarmUpdatePosition_MPntStd_Euler(dav,velocity,dt,npoints,mp_std);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* UPDATE local coordinates */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateLocalCoordinates"
PetscErrorCode MaterialPointStd_UpdateLocalCoordinates(DataBucket materialpoints,DM dav)
{
	PetscErrorCode ierr;
	int            p,npoints;
	MPntStd        *mp_std;
	DataField      PField;
	double         tolerance;
	int            max_its;
	PetscBool      use_nonzero_guess, monitor, log;
	DM             cda;
	Vec            gcoords;
	PetscScalar    *LA_gcoords;
	const PetscInt *elnidx_u;
	PetscInt       nel,nen_u;
	PetscInt       *gidx;
	PetscInt       lmx,lmy,lmz;
	

	PetscFunctionBegin;
	/* get marker fields */
	DataBucketGetSizes(materialpoints,&npoints,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dav,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dav,0,&gidx);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	
	/* point location parameters */
	tolerance         = 1.0e-10;
	max_its           = 10;
	use_nonzero_guess = PETSC_TRUE;
	monitor           = PETSC_FALSE;
	log               = PETSC_FALSE;
	
	InverseMappingDomain_3dQ2( 		 tolerance, max_its,
														use_nonzero_guess, 
														monitor, log,
														(const double*)LA_gcoords, (const int)lmx,(const int)lmy,(const int)lmz, (const int*)elnidx_u,
														npoints, mp_std );
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* remove all points which didn't find a home */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_Removal"
PetscErrorCode MaterialPointStd_Removal(DataBucket materialpoints)
{
	int       p,npoints,escaped;
	MPntStd   *mp_std;
	DataField PField;
	
	
	PetscFunctionBegin;
	/* get marker fields */
	DataBucketGetSizes(materialpoints,&npoints,PETSC_NULL,PETSC_NULL);
	DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
	escaped = 0;
	for (p=0; p<npoints; p++) {
		if (mp_std[p].wil == -1) {
			escaped++;
		}
	}
	
	/* remove points which left processor */
	if (escaped!=0) {
		PetscPrintf(PETSC_COMM_SELF,"  *** MPntStd removal: Identified %d points which are not contained on subdomain (after communication) \n", escaped );
	}
	if (escaped!=0) {
		for (p=0; p<npoints; p++) {
			if (mp_std[p].wil == -1) {
				//printf("############################ Point %d not located in domain (%1.6e , %1.6e) ############################# \n",p,mp_std[p].coor[0],mp_std[p].coor[1]);
				
				/* kill point */
				DataBucketRemovePointAtIndex(materialpoints,p);
				DataBucketGetSizes(materialpoints,&npoints,0,0); /* you need to update npoints as the list size decreases! */
				p--; /* check replacement point */
			}
		}
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_Communication_Generic"
PetscErrorCode SwarmUpdatePosition_Communication_Generic(DataBucket db,DM da,DataEx de)
{
	DataField      PField_std;
	DataField      *PField_all;
	int            f,nfields;
	size_t         sizeof_marker_contents;
	int            p,npoints,npoints_global_init,npoints_global_fin;
	void           *recv_data;
	void           *data_p;
	PetscMPIInt    n,neighborcount, *neighborranks2;
	PetscInt       recv_length;
	int            npoints_accepted;
	PetscMPIInt    rank,size;
	MPntStd        *marker_std;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* communucate */
	MPI_Comm_size(((PetscObject)da)->comm,&size);
	if (size==1) {
		PetscFunctionReturn(0);
	}
	
	MPI_Comm_rank(((PetscObject)da)->comm,&rank);
	
	neighborcount  = de->n_neighbour_procs;
	neighborranks2 = de->neighbour_procs;
	
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataBucketGetDataFields(db,&nfields,&PField_all);
	for (f=1; f<nfields; f++) {
		DataFieldGetAccess(PField_all[f]);
	}
	
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	DataBucketGetSizes(db,&npoints,0,0);
	
	MPI_Allreduce(&npoints,&npoints_global_init,1,MPI_INT,MPI_SUM,de->comm);
	
	/* figure out how many points left processor */
	ierr = DataExInitializeSendCount(de);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		PetscBool onproc;
		MPntStd   *marker;
		
		DataFieldAccessPoint(PField_std,p,(void**)&marker);
		onproc = PETSC_TRUE;
		if (marker->wil==-1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc==PETSC_FALSE) {
			for (n=0; n<neighborcount; n++) {
				//	printf("  DataEx: rank %d sending %d to %d \n",rank,1,neighborranks2[n] );
				ierr = DataExAddToSendCount( de, neighborranks2[n], 1 );CHKERRQ(ierr);
			}
		}
	}
	ierr = DataExFinalizeSendCount(de);CHKERRQ(ierr);
	
	DataFieldRestoreAccess(PField_std);
	
	/* pack points which left processor */
	marker_std    = PField_std->data; /* should write a function to do this */
	
	sizeof_marker_contents = 0;
	sizeof_marker_contents = sizeof_marker_contents + sizeof(MPntStd);
	for (f=1; f<nfields; f++) {
		sizeof_marker_contents += PField_all[f]->atomic_size;
	}
	
	ierr = DataExPackInitialize(de,sizeof_marker_contents);CHKERRQ(ierr);
	
	/* allocate a temporary buffer whihc is large enough to store all marker fields from an individual point,p */
	ierr = PetscMalloc(sizeof_marker_contents,&data_p);CHKERRQ(ierr);
	
	for (p=0; p<npoints; p++) {
		PetscBool   onproc;
		MPntStd     *marker_p;
		void        *marker_p_prop;
		size_t      atomic_size,offset;
		
		/* access fields from the bucket */
		marker_p     = &marker_std[p];
		
		onproc = PETSC_TRUE;
		if (marker_p->wil==-1) {
			onproc = PETSC_FALSE;
			/* pack together */
			//(void*)((char*)field->data + index*field->atomic_size)
			
			/* copy all fields from the bucket into a temporary buffer */
			PetscMemcpy((void*)data_p,                         marker_p,    sizeof(MPntStd));
			
			offset = sizeof(MPntStd);
			for (f=1; f<nfields; f++) {
				atomic_size = PField_all[f]->atomic_size;
				DataFieldAccessPoint(PField_all[f],p,&marker_p_prop);
				
				PetscMemcpy((void*)((char*)data_p+offset),marker_p_prop,atomic_size);
				offset = offset + atomic_size;
			}
		}
		
		if (onproc==PETSC_FALSE) {
			for (n=0; n<neighborcount; n++) {
				ierr = DataExPackData( de, neighborranks2[n], 1,(void*)data_p );CHKERRQ(ierr);
			}
		}
	}		
	for (f=1; f<nfields; f++) {
		DataFieldRestoreAccess(PField_all[f]);
	}
	ierr = DataExPackFinalize(de);CHKERRQ(ierr);
	
	/* remove points which left processor */
	DataBucketGetSizes(db,&npoints,0,0);
	DataFieldGetAccess(PField_std);
	for (p=0; p<npoints; p++) {
		PetscBool onproc;
		MPntStd   *marker_p;
		
		DataFieldAccessPoint(PField_std,p,(void**)&marker_p);
		onproc = PETSC_TRUE;
		if (marker_p->wil==-1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc==PETSC_FALSE) { 
			/* kill point */
			DataBucketRemovePointAtIndex(db,p);
			DataBucketGetSizes(db,&npoints,0,0); /* you need to update npoints as the list size decreases! */
			p--; /* check replacement point */
		}
	}		
	DataFieldRestoreAccess(PField_std);
	
	// START communicate //
	ierr = DataExBegin(de);CHKERRQ(ierr);
	ierr = DataExEnd(de);CHKERRQ(ierr);
	// END communicate //
	
	// receive, if i own them, add new points to list //
	ierr = DataExGetRecvData( de, &recv_length, (void**)&recv_data );CHKERRQ(ierr);
	{
		PetscInt totalsent;
		MPI_Allreduce(&recv_length,&totalsent,1,MPIU_INT,MPI_SUM,de->comm);
		PetscPrintf(PETSC_COMM_WORLD,"  DataEx: total points sent = %d \n", totalsent);
	}
	
	
	
	/* update the local coordinates and cell owner for all recieved points */
	{
		DM cda;
		Vec gcoords;
		PetscScalar *LA_gcoords;
		double tolerance;
		int max_its;
		PetscBool use_nonzero_guess, monitor, log;
		PetscInt lmx,lmy,lmz;
		PetscInt nel,nen_u;
		MPntStd *marker_p;
		const PetscInt *elnidx_u;
		
		/* setup for coords */
		ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
		ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
		ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
		
		ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		
		/* point location parameters */
		tolerance         = 1.0e-10;
		max_its           = 10;
		use_nonzero_guess = PETSC_FALSE; /* for markers sent across processors, it is necessary to NOT use the last known values! */
		monitor           = PETSC_FALSE;
		log               = PETSC_FALSE;
		
		/* 
		 NOTE: My algorithm will break down if you try to use an non-zero initial guess
		 if the marker has been sent to another cpu. This could be fixed, however if would
		 require additional logic to be added into the function InverseMappingDomain_2dQ2()
		 
		 It seems more reasonable to simply assume that if the marker was sent to another cpu,
		 the value of wil currently stored on the marker is completely meaningless and thus we
		 should ALWAYS use a zero initial guess (i.e. assume we know nothing)
		 */
		
		for (p=0; p<recv_length; p++) {
			marker_p = (MPntStd*)( (char*)recv_data + p*(sizeof_marker_contents) );
			//printf("revc p=%d : (%lf,%lf) \n", p, marker_p->coor[0], marker_p->coor[1] );
			
			InverseMappingDomain_3dQ2( 		 tolerance, max_its,
																use_nonzero_guess, 
																monitor, log,
																(const double*)LA_gcoords, (const int)lmx,(const int)lmy,(const int)lmz, (const int*)elnidx_u,
																1, marker_p );
			//printf("++[%d] revc p=%d : wil=%d \n", rank,p, marker_p->wil );
		}
		
		ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	}
	
	/* accepte all points living locally */
	npoints_accepted = 0;
	for (p=0; p<recv_length; p++) {
		PetscBool   onproc;
		MPntStd     *marker_p;
		size_t      offset;
		
		offset = 0;
		marker_p     = (MPntStd*)(     (char*)recv_data + p*(sizeof_marker_contents) + offset);
		
		onproc = PETSC_TRUE;
		if (marker_p->wil==-1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc == PETSC_TRUE) {
			int end;
			
			DataBucketAddPoint(db);
			DataBucketGetSizes(db,&end,0,0);
			end = end - 1;
			
			offset = 0;
			for (f=0; f<nfields; f++) {
				void *data_p = (void*)( (char*)recv_data + p*(sizeof_marker_contents) + offset );
				
				DataFieldInsertPoint(PField_all[f], end, (void*)data_p );
				
				offset = offset + PField_all[f]->atomic_size;
			}
			
			
			npoints_accepted++;
		}
	}	
	//printf("  DataEx: rank %d accepted %d new points \n",rank,npoints_accepted );
	
	
	DataBucketGetSizes(db,&npoints,0,0);
	MPI_Allreduce(&npoints,&npoints_global_fin,1,MPI_INT,MPI_SUM,de->comm);
	PetscPrintf(PETSC_COMM_WORLD,"  SwarmUpdatePosition_GENERIC(Communication): num. points global ( init. = %d : final = %d )\n", npoints_global_init,npoints_global_fin);
	ierr = PetscFree(data_p);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateCoordinates"
PetscErrorCode MaterialPointStd_UpdateCoordinates(DataBucket materialpoints,DM dav,DataEx de)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	PetscLogDouble t0,t1,tl,tg,tgmn;
	
	PetscPrintf(PETSC_COMM_WORLD,"\n\n=========== MaterialPointStd_UpdateCoordinates ============= \n");
	PetscGetTime(&t0);
	ierr = MaterialPointStd_UpdateLocalCoordinates(materialpoints,dav);CHKERRQ(ierr);
	PetscGetTime(&t1);
	tl = t1 - t0;
	ierr = MPI_Allreduce(&tl,&tg,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&tgmn,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"    	MaterialPointStd_UpdateLocalCoordinates   %10.2e (sec) efficiency %1.1lf%%\n", tg, 100.0 - 100.0*(tg-tgmn)/tg );
	
	PetscGetTime(&t0);
	ierr = SwarmUpdatePosition_Communication_Generic(materialpoints,dav,de);CHKERRQ(ierr);
	PetscGetTime(&t1);
	ierr = MPI_Allreduce(&tl,&tg,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&tgmn,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"    	SwarmUpdatePosition_Communication_Generic %10.2e (sec) efficiency %1.1lf%%\n", tg, 100.0 - 100.0*(tg-tgmn)/tg );
	tl = t1 - t0;

	PetscGetTime(&t0);
	ierr = MaterialPointStd_Removal(materialpoints);CHKERRQ(ierr);
	PetscGetTime(&t1);
	tl = t1 - t0;
	ierr = MPI_Allreduce(&tl,&tg,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&tl,&tgmn,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"    	MaterialPointStd_Removal                  %10.2e (sec) efficiency %1.1lf%%\n", tg, 100.0 - 100.0*(tg-tgmn)/tg );
	
	PetscPrintf(PETSC_COMM_WORLD,"=========== ================================== ============= \n\n");

	
	PetscFunctionReturn(0);
}

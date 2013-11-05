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
 **    Filename:      material_point_point_location.c
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



#include "element_type_Q2.h"
#include "element_utils_q2.h"
#include "ptatin3d_defs.h"

#include "MPntStd_def.h"
#include "material_point_point_location.h"

/* debugging variable for point location routines */
//#define PNTLOC_LOG

#if 0
void _compute_deltaX_3d( double J[3][3], double f[], double h[] )
{
	int    i;
	double invJ[3][3];
	
	ElementHelper_matrix_inverse_3x3(J,invJ);

	h[0] = h[1] = h[2] = 0.0;	
	for (i=0; i<3; i++) {
		h[0] = h[0] + invJ[0][i] * f[i];
		h[1] = h[1] + invJ[1][i] * f[i];
		h[2] = h[2] + invJ[2][i] * f[i];
	}
}

void _compute_J_3dQ2(double xi[],double vertex[],double J[3][3])
{
	int    i;
	double GNi[3][Q2_NODES_PER_EL_3D];
	
	J[0][0] = J[0][1] = J[0][2] = 0.0;
	J[1][0] = J[1][1] = J[1][2] = 0.0;
	J[2][0] = J[2][1] = J[2][2] = 0.0;
	
	P3D_ConstructGNi_Q2_3D(xi,GNi);
	for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
		int i2 = i*3;
		double x = vertex[i2];
		double y = vertex[i2+1];
		double z = vertex[i2+2];
		
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
	int    i;
	double Ni[Q2_NODES_PER_EL_3D];
	
	/* Update F for the next iteration */
	f[0] = f[1] = f[2] = 0.0;
	
	P3D_ConstructNi_Q2_3D(xi,Ni);
	
	for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
		int i2   = i*3;
		int i2p1 = i2+1;
		int i2p2 = i2+2;
		
		f[0] += vertex[i2]   * Ni[i];
		f[1] += vertex[i2p1] * Ni[i];
		f[2] += vertex[i2p2] * Ni[i];
	}
	f[0] = - f[0] + pos[0];
	f[1] = - f[1] + pos[1];
	f[2] = - f[2] + pos[2];
}

void InverseMappingDomain_3dQ2( 
															 double tolerance, int max_its,
															 Truth use_nonzero_guess, 
															 Truth monitor, Truth log,
															 const double coords[], const int mx, const int my, const int mz, const int element[],
															 int np, MPntStd marker[] )
{
	const int dim = 3;
	const int nodesPerEl = Q2_NODES_PER_EL_3D; 
	double h[dim];
	double Jacobian[dim][dim];
	double f[dim];
	int    i,d;
	int    its;
	double residual2,tolerance2,F2;
	int    p;
	int    mx_origin, my_origin, wil_origin;
	double cxip[dim],Lxip[dim],Gxip[dim];
	double dxi,deta,dzeta,xi0,eta0,zeta0;
	int    I,J,K,wil_IJ,eid,k;
	double vertex[dim * nodesPerEl];
	int    n0,n1,n2,n3;
	Truth  point_found;
	
	tolerance2 = tolerance * tolerance; /* Eliminates the need to do a sqrt in the convergence test */
	
	if(log)printf("Domain: ncells = %d x %d x %d = %d \n", mx,my,mz,mx*my*mz );
	
	/* map domain to [-1,1]x[-1,1]x[-1,1] domain */
	dxi   = 2.0/((double)mx);
	deta  = 2.0/((double)my);
	dzeta = 2.0/((double)mz);
	if(log)printf("Domain: (dxi,eta,zeta) = (%1.8e,%1.8e,%1.8e)\n",dxi,deta,dzeta );
	
	for (p=0; p<np; p++) {
		MPntStd *marker_p = &marker[p];
		
		/* copy these values */
		cxip[0] = marker_p->xi[0];
		cxip[1] = marker_p->xi[1];
		cxip[2] = marker_p->xi[2];
		
		/* Check for an initial guess initial guess */
		if (use_nonzero_guess == _FALSE) {
			Gxip[0] = 0.0;
			Gxip[1] = 0.0;
			Gxip[2] = 0.0;
		} else {
			/* convert wil => IJ */
			wil_IJ = marker_p->wil;
			K = wil_IJ / (mx*my);
			J = wil_IJ - K*(mx*my);
			I = wil_IJ - J*mx;
			if(log)printf("init I,J,K = %d %d %d \n", I,J,K );
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
			if(log)printf("[Lxi-init] = %1.8e %1.8e %1.8e \n", cxip[0], cxip[1], cxip[2] );
			if(log)printf("[Gxi-init] = %1.8e %1.8e %1.8e \n", Gxip[0], Gxip[1], Gxip[2] );
			
		}
		if (monitor) {
			printf("point[%d]: pos = ( %+1.8e, %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e, %+1.8e ) \n", 
											p, marker_p->coor[0],marker_p->coor[1],marker_p->coor[2],
											marker_p->xi[0], marker_p->xi[1], marker_p->xi[2] );
		}
		
		point_found = _FALSE;
		
		its = 0;
		do {
			if(log)printf("iteration: %d\n",its);
			/* convert Gxi to IJ */
			I = (Gxip[0]+1.0)/dxi;
			J = (Gxip[1]+1.0)/deta;
			K = (Gxip[2]+1.0)/dzeta;
			
			if (I == mx) { I--; }
			if (J == my) { J--; }
			if (K == mz) { K--; }
			
			if ((I<0) || (J<0) || (K<0)) {
				if(log)printf("  I(%d),J(%d),K(%d) negative Gxip %1.8e,%1.8e,%1.8e \n",I,J,K,Gxip[0],Gxip[1],Gxip[2]);
				break;
			}
			if (I >= mx) { 
				if(log)printf("  I too large \n");
				break;
			}
			if (J >= my) {
				if(log)printf("  J too large \n");
				break;
			}
			if (K >= mz) {
				if(log)printf("  K too large \n");
				break;
			}
			
			
			/* Get coords of cell IJ */
			wil_IJ = I + J * mx + K * mx * my;
			if(log)printf("  I,J,K=%d/%d/%d : wil_IJ %d : nid = ", I,J,K,wil_IJ);
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
				int nid = element[wil_IJ*Q2_NODES_PER_EL_3D+k];
				
				vertex[dim*k+0] = coords[dim*nid+0];
				vertex[dim*k+1] = coords[dim*nid+1];
				vertex[dim*k+2] = coords[dim*nid+2];
				if(log)printf("%d ", nid);
			}
			if(log)printf("\n");
			
			
			if(log) {
				printf("  [vertex] ");
				for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
					printf("(%1.8e , %1.8e , %1.8e) ",vertex[dim*k+0],vertex[dim*k+1],vertex[dim*k+2] );
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
			
			if(log)printf("  Lxi,Lxeta,Lxzeta = %1.8e, %1.8e, %1.8e (%d,%d,%d) \n", Lxip[0],Lxip[1],Lxip[2],I,J,K );
			
			_compute_F_3dQ2( Lxip, vertex, marker_p->coor, f );
			if (monitor) {
				printf("%4d InverseMapping : F = ( %+1.8e, %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e, %+1.8e ) \n", 
							 its, f[0],f[1],f[2], Lxip[0],Lxip[1],Lxip[2] );
			}
			
			/* Check for convergence */
			F2 = (f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
			if (F2 < tolerance2) {
				if (monitor) printf("%4d InverseMapping : converged : Norm of F %1.8e \n", its, sqrt(F2) );
				point_found = _TRUE;
				break;
			}
			
			_compute_J_3dQ2( Lxip, vertex, Jacobian );
			
			/* compute update */
			_compute_deltaX_3d( Jacobian, f, h );
			if(log)printf("  [delta] = %1.8e %1.8e %1.8e \n", h[0],h[1],h[2] );
			
			/* update Lxip */
			Lxip[0] += 10.0e-1 *h[0];
			Lxip[1] += 10.0e-1 *h[1];
			Lxip[2] += 10.0e-1 *h[2];
			if(log)printf("  [corrected] Lxi,Lxeta,Lxzeta = %1.8e, %1.8e, %1.8e \n", Lxip[0],Lxip[1],Lxip[2] );
			
			residual2 = ( h[0]*h[0] + h[1]*h[1] + h[2]*h[2] );
			if (residual2 < tolerance2) {
				if (monitor) printf("%4d InverseMapping : converged : Norm of correction %1.8e \n", its, sqrt(residual2) );
				point_found = _TRUE;
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
			if(log)printf("  [Gxi] = %1.8e %1.8e %1.8e \n", Gxip[0], Gxip[1], Gxip[2] );
			
			for (d=0; d<dim; d++) {
				if (Gxip[d]<-1.0) { 
					Gxip[d] = -1.0;
					if(log)printf("  correction outside box: correcting \n");
				}
				
				if (Gxip[d]>1.0) {
					Gxip[d] = 1.0;
					if(log)printf("  correction outside box: correcting \n");
				}
			}
			
			
			its++;
		} while(its < max_its);
		
		if (monitor && point_found==_FALSE){
			if (its >= max_its) {
				printf("%4d %s : Reached maximum iterations (%d) without converging. \n", its, __FUNCTION__, max_its );
			}	else {
				printf("%4d %s : Newton broke down, diverged or stagnated after (%d) iterations without converging. \n", its, __FUNCTION__, its );
			}
		}
		
		/* if at the end of the solve, it still looks like the point is outside the mapped domain, mark point as not being found */
		if (fabs(Gxip[0]) > 1.0) { point_found =_FALSE; }
		if (fabs(Gxip[1]) > 1.0) { point_found =_FALSE; }
		if (fabs(Gxip[2]) > 1.0) { point_found =_FALSE; }
		
		/* update local variables */
		if (point_found == _FALSE) {
			Lxip[0] = NAN;
			Lxip[1] = NAN;
			Lxip[2] = NAN;
			wil_IJ  = -1;
		}	else {
			/* convert Gxi to IJ */
			I = (Gxip[0]+1.0)/dxi;
			J = (Gxip[1]+1.0)/deta;
			K = (Gxip[2]+1.0)/dzeta;
			if (I == mx) { I--; }
			if (J == my) { J--; }
			if (K == mz) { K--; }
			
			if (I >= mx) {
				if(log)printf("  I too large \n");
				break;
			}
			if (J >= my) {
				if(log)printf("  J too large \n");
				break;
			}
			if (K >= mz) {
				if(log)printf("  K too large \n");
				break;
			}
			
			/* convert global (domain) xi TO local (element) xi  */
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			xi0   = -1.0 + I*dxi;
			eta0  = -1.0 + J*deta;
			zeta0 = -1.0 + K*dzeta;
			
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0)  /dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0) /deta  - 1.0;
			Lxip[2] = 2.0*(Gxip[2]-zeta0)/dzeta - 1.0;
			
			wil_IJ = I + J * mx + K * mx * my;
		}
		
		/* set into vector */
		marker_p->xi[0] = Lxip[0];
		marker_p->xi[1] = Lxip[1];
		marker_p->xi[2] = Lxip[2];
		marker_p->wil   = wil_IJ;
	}
}
#endif

/* 2d implementation */
void _compute_deltaX_2d( double J[2][2], double f[], double h[] )
{
	int    i;
	double invJ[2][2];
	
	ElementHelper_matrix_inverse_2x2(J,invJ);
	
	h[0] = h[1] = 0.0;	
	for (i=0; i<2; i++) {
		h[0] = h[0] + invJ[0][i] * f[i];
		h[1] = h[1] + invJ[1][i] * f[i];
	}
}

void _compute_J_2dQ2(double xi[],double vertex[],double J[2][2])
{
	int    i;
	double GNi[2][Q2_NODES_PER_EL_2D];
	
	J[0][0] = J[0][1] = 0.0;
	J[1][0] = J[1][1] = 0.0;
	
	P3D_ConstructGNi_Q2_2D(xi,GNi);
	for (i=0; i<Q2_NODES_PER_EL_2D; i++) {
		int    i2 = i*2;
		double x = vertex[i2];
		double y = vertex[i2+1];
		
		J[0][0] += x * GNi[0][i];
		J[0][1] += x * GNi[1][i];
		
		J[1][0] += y * GNi[0][i];
		J[1][1] += y * GNi[1][i];
	}
}

void _compute_F_2dQ2(double xi[],double vertex[],double pos[],double f[])
{
	int    i;
	double Ni[Q2_NODES_PER_EL_2D];
	
	/* Update F for the next iteration */
	f[0] = f[1] = 0.0;
	
	P3D_ConstructNi_Q2_2D(xi,Ni);
	
	for (i=0; i<Q2_NODES_PER_EL_2D; i++) {
		int i2   = i*2;
		int i2p1 = i2+1;
		
		f[0] += vertex[i2]   * Ni[i];
		f[1] += vertex[i2p1] * Ni[i];
	}
	f[0] = - f[0] + pos[0];
	f[1] = - f[1] + pos[1];
}

void InverseMappingDomain_2dQ2( 
															 double tolerance, int max_its,
															 Truth use_nonzero_guess, 
															 Truth monitor,
															 const double coords[], const int mx, const int my, const int element[],
															 int np, MPntStd marker[] )
{
	const int dim = 2;
	const int nodesPerEl = Q2_NODES_PER_EL_2D; 
	double h[dim];
	double Jacobian[2][2];
	double f[2];
	int    i,d;
	int    its;
	double residual2,tolerance2,F2;
	int    p;
	int    mx_origin, my_origin, wil_origin;
	double cxip[2],Lxip[2],Gxip[2];
	double dxi,deta,xi0,eta0;
	int    I,J,wil_IJ,eid,k;
	double vertex[2 * Q2_NODES_PER_EL_2D];
	int    n0,n1,n2;
	Truth  point_found;
	
	tolerance2 = tolerance * tolerance; /* Eliminates the need to do a sqrt in the convergence test */

#ifdef PNTLOC_LOG
	printf("Domain: ncells = %d x %d = %d \n", mx,my,mx*my );
#endif
	
	/* map domain to [-1,1]x[-1,1] domain */
	dxi   = 2.0/((double)mx);
	deta  = 2.0/((double)my);
#ifdef PNTLOC_LOG
	printf("Domain: (dxi,eta) = (%1.8e,%1.8e)\n",dxi,deta );
#endif
	
	for (p=0; p<np; p++) {
		MPntStd *marker_p = &marker[p];
		
		/* copy these values */
		cxip[0] = marker_p->xi[0];
		cxip[1] = marker_p->xi[1];
		
		/* Check for an initial guess initial guess */
		if (use_nonzero_guess == _FALSE) {
			marker_p->xi[0] = 0.0;
			marker_p->xi[1] = 0.0;
			
			cxip[0] = marker_p->xi[0];
			cxip[1] = marker_p->xi[1];

			Gxip[0] = 0.0;
			Gxip[1] = 0.0;
		}	else {
			/* convert wil => IJ */
			wil_IJ = marker_p->wil;
			J = wil_IJ / mx;
			I = wil_IJ - J*mx;

			/* convert Lxip => Gxip */
			xi0   = -1.0 + I*dxi;
			eta0  = -1.0 + J*deta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Gxip[0] = 2.0 * (cxip[0]-xi0)/dxi - 1.0;
			//			Gxip[1] = 2.0 * (cxip[1]-eta0)/deta - 1.0;
			// x*-x*0/dx = (x+1)/2
			Gxip[0] = dxi   * (cxip[0]+1.0)/2.0 + xi0;
			Gxip[1] = deta  * (cxip[1]+1.0)/2.0 + eta0;
#ifdef PNTLOC_LOG
			printf("init I,J = %d %d \n", I,J );
			printf("[Lxi-init] = %1.8e %1.8e \n", cxip[0], cxip[1] );
			printf("[Gxi-init] = %1.8e %1.8e \n", Gxip[0], Gxip[1] );
#endif
		}
		if (monitor) {
			printf("point[%d]: pos = ( %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e ) \n", 
						 p, marker_p->coor[0],marker_p->coor[1],
						 marker_p->xi[0], marker_p->xi[1] );
		}
		
		point_found = _FALSE;
		
		its = 0;
		do {
#ifdef PNTLOC_LOG
			printf("iteration: %d\n",its);
#endif
			/* convert Gxi to IJ */
			I = (Gxip[0]+1.0)/dxi;
			J = (Gxip[1]+1.0)/deta;
			
			if (I == mx) { I--; }
			if (J == my) { J--; }
			
			if ((I<0) || (J<0)) {
#ifdef PNTLOC_LOG
				printf("  I(%d),J(%d) negative Gxip %1.8e,%1.8e \n",I,J,Gxip[0],Gxip[1]);
#endif
				break;
			}
			if (I >= mx) { 
#ifdef PNTLOC_LOG
				printf("  I too large \n");
#endif
				break;
			}
			if (J >= my) {
#ifdef PNTLOC_LOG
				printf("  J too large \n");
#endif
				break;
			}
			
			/* Get coords of cell IJ */
			wil_IJ = I + J * mx;
#ifdef PNTLOC_LOG
			printf("  I,J=%d/%d : wil_IJ %d : nid = ", I,J,wil_IJ);
#endif
			for (k=0; k<nodesPerEl; k++) {
				int nid = element[wil_IJ*nodesPerEl+k];
				
				vertex[dim*k+0] = coords[dim*nid+0];
				vertex[dim*k+1] = coords[dim*nid+1];
#ifdef PNTLOC_LOG
				printf("%d ", nid);
#endif
			}
#ifdef PNTLOC_LOG
			printf("\n");
#endif
			
#ifdef PNTLOC_LOG
			printf("  [vertex] ");
			for (k=0; k<nodesPerEl; k++) {
				printf("(%1.8e , %1.8e) ",vertex[dim*k+0],vertex[dim*k+1] );
			}
			printf("\n");
#endif
			
			/* convert global (domain) xi TO local (element) xi  */
			xi0   = -1.0 + I*dxi;
			eta0  = -1.0 + J*deta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Lxip[0] = 0.5*(Gxip[0]+1.0)*dxi  + xi0;
			//			Lxip[1] = 0.5*(Gxip[1]+1.0)*deta + eta0;
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0  )/dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0 )/deta  - 1.0;
			
#ifdef PNTLOC_LOG
			printf("  Lxi,Lxeta = %1.8e, %1.8e (%d,%d) \n", Lxip[0],Lxip[1],I,J );
#endif
			
			_compute_F_2dQ2( Lxip, vertex, marker_p->coor, f );
			if (monitor) {
				printf("%4d InverseMapping : F = ( %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e ) \n", 
							 its, f[0],f[1], Lxip[0],Lxip[1] );
			}
			
			/* Check for convergence */
			F2 = (f[0]*f[0]+f[1]*f[1]);
			if (F2 < tolerance2) {
				if (monitor) printf("%4d InverseMapping : converged : Norm of F %1.8e \n", its, sqrt(F2) );
				point_found = _TRUE;
				break;
			}
			
			_compute_J_2dQ2( Lxip, vertex, Jacobian );
			
			/* compute update */
			_compute_deltaX_2d( Jacobian, f, h );
#ifdef PNTLOC_LOG
			printf("  [delta] = %1.8e %1.8e \n", h[0],h[1] );
#endif
			
			/* update Lxip */
			Lxip[0] += 10.0e-1 *h[0];
			Lxip[1] += 10.0e-1 *h[1];
#ifdef PNTLOC_LOG
			printf("  [corrected] Lxi,Lxeta = %1.8e, %1.8e \n", Lxip[0],Lxip[1] );
#endif
			
			residual2 = ( h[0]*h[0] + h[1]*h[1] );
			if (residual2 < tolerance2) {
				if (monitor) printf("%4d InverseMapping : converged : Norm of correction %1.8e \n", its, sqrt(residual2) );
				point_found = _TRUE;
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
#ifdef PNTLOC_LOG
			printf("  [Gxi] = %1.8e %1.8e \n", Gxip[0], Gxip[1] );
#endif
			
			for (d=0; d<dim; d++) {
				if (Gxip[d] < -1.0) { 
					Gxip[d] = -1.0;
#ifdef PNTLOC_LOG
					printf("  correction outside box: correcting \n");
#endif
				}
				
				if (Gxip[d] > 1.0) {
					Gxip[d] = 1.0;
#ifdef PNTLOC_LOG
					printf("  correction outside box: correcting \n");
#endif
				}
			}
			
			its++;
		} while(its < max_its);
		
		if (monitor && point_found == _FALSE){
			if (its >= max_its) {
				printf("%4d %s : Reached maximum iterations (%d) without converging. \n", its, __FUNCTION__, max_its );
			} else {
				printf("%4d %s : Newton broke down, diverged or stagnated after (%d) iterations without converging. \n", its, __FUNCTION__, its );
			}
		}
		
		/* if at the end of the solve, it still looks like the point is outside the mapped domain, mark point as not being found */
		if (fabs(Gxip[0]) > 1.0) { point_found = _FALSE; }
		if (fabs(Gxip[1]) > 1.0) { point_found = _FALSE; }
		
		/* update local variables */
		if (point_found == _FALSE) {
			Lxip[0] = NAN;
			Lxip[1] = NAN;
			wil_IJ  = -1;
		}	else {
			/* convert Gxi to IJ */
			I = (Gxip[0]+1.0)/dxi;
			J = (Gxip[1]+1.0)/deta;
			if (I == mx){ I--; }
			if (J == my){ J--; }
			
			if (I >= mx) {
#ifdef PNTLOC_LOG
				printf("  I too large \n");
#endif
				break;
			}
			if (J >= my) {
#ifdef PNTLOC_LOG
				printf("  J too large \n");
#endif
				break;
			}
			
			/* convert global (domain) xi TO local (element) xi  */
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			xi0   = -1.0 + I*dxi;
			eta0  = -1.0 + J*deta;
			
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0)  /dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0) /deta  - 1.0;
			
			wil_IJ = I + J * mx;
		}
		
		/* set into vector */
		marker_p->xi[0] = Lxip[0];
		marker_p->xi[1] = Lxip[1];
		marker_p->wil   = wil_IJ;
	}
}


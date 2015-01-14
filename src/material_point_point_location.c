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
//#define PTAT3D_DBG_PointLocation

#define NSD2d 2
#define NSD3d 3

/* 2d implementation */
void _compute_deltaX_2d(PetscReal J[2][2],PetscReal f[],PetscReal h[])
{
	PetscInt  i;
	PetscReal invJ[NSD2d][NSD2d];
	
    
	ElementHelper_matrix_inverse_2x2(J,invJ);
	
	h[0] = h[1] = 0.0;	
	for (i=0; i<NSD2d; i++) {
		h[0] = h[0] + invJ[0][i] * f[i];
		h[1] = h[1] + invJ[1][i] * f[i];
	}
}

void _compute_J_2dQ2(PetscReal xi[],PetscReal vertex[],PetscReal J[2][2])
{
	PetscInt  i;
	PetscReal GNi[NSD2d][Q2_NODES_PER_EL_2D];
	
    
	J[0][0] = J[0][1] = 0.0;
	J[1][0] = J[1][1] = 0.0;
	
	P3D_ConstructGNi_Q2_2D(xi,GNi);
	for (i=0; i<Q2_NODES_PER_EL_2D; i++) {
		PetscInt  i2 = i * NSD2d;
		PetscReal x = vertex[i2];
		PetscReal y = vertex[i2+1];
		
		J[0][0] += x * GNi[0][i];
		J[0][1] += x * GNi[1][i];
		
		J[1][0] += y * GNi[0][i];
		J[1][1] += y * GNi[1][i];
	}
}

void _compute_F_2dQ2(PetscReal xi[],PetscReal vertex[],PetscReal pos[],PetscReal f[])
{
	PetscInt  i;
	PetscReal Ni[Q2_NODES_PER_EL_2D];
    
	
	/* Update F for the next iteration */
	f[0] = f[1] = 0.0;
	
	P3D_ConstructNi_Q2_2D(xi,Ni);
	
	for (i=0; i<Q2_NODES_PER_EL_2D; i++) {
		PetscInt i2   = i * NSD2d;
		PetscInt i2p1 = i2+1;
		
		f[0] += vertex[i2]   * Ni[i];
		f[1] += vertex[i2p1] * Ni[i];
	}
	f[0] = - f[0] + pos[0];
	f[1] = - f[1] + pos[1];
}

void InverseMappingDomain_2dQ2(PetscReal tolerance,PetscInt max_its,
                               PetscBool use_nonzero_guess,
                               PetscBool monitor,
                               const PetscReal coords[],const PetscInt mx,const PetscInt my,const PetscInt element[],
                               PetscInt np,MPntStd marker[] )
{
	PetscReal h[NSD2d];
	PetscReal Jacobian[NSD2d][NSD2d];
	PetscReal f[NSD2d];
	PetscInt  d;
	PetscInt  its;
	PetscReal residual2,tolerance2,F2;
	PetscInt  p;
	PetscReal cxip[NSD2d],Lxip[NSD2d],Gxip[NSD2d];
	PetscReal dxi,deta,xi0,eta0;
	PetscInt  nI,nJ,wil_nInJ,k;
	PetscReal vertex[NSD2d * Q2_NODES_PER_EL_2D];
	PetscBool point_found;
	
    
	tolerance2 = tolerance * tolerance; /* Eliminates the need to do a sqrt in the convergence test */

#ifdef PTAT3D_DBG_PointLocation
	PetscPrintf(PETSC_COMM_SELF,"Domain: ncells = %D x %D = %D \n", mx,my,mx*my );
#endif
	
	/* map domain to [-1,1]x[-1,1] domain */
	dxi   = 2.0/((PetscReal)mx);
	deta  = 2.0/((PetscReal)my);
#ifdef PTAT3D_DBG_PointLocation
	PetscPrintf(PETSC_COMM_SELF,"Domain: (dxi,eta) = (%1.8e,%1.8e)\n",dxi,deta );
#endif
	
	for (p=0; p<np; p++) {
		MPntStd *marker_p = &marker[p];
		
		/* copy these values */
		cxip[0] = marker_p->xi[0];
		cxip[1] = marker_p->xi[1];
		
		/* Check for an initial guess initial guess */
		if (use_nonzero_guess == PETSC_FALSE) {
			marker_p->xi[0] = 0.0;
			marker_p->xi[1] = 0.0;
			
			cxip[0] = marker_p->xi[0];
			cxip[1] = marker_p->xi[1];

			Gxip[0] = 0.0;
			Gxip[1] = 0.0;
		}	else {
			/* convert wil => nInJ */
			wil_nInJ = marker_p->wil;
			nJ = wil_nInJ / mx;
			nI = wil_nInJ - nJ*mx;

			/* convert Lxip => Gxip */
			xi0   = -1.0 + nI*dxi;
			eta0  = -1.0 + nJ*deta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Gxip[0] = 2.0 * (cxip[0]-xi0)/dxi - 1.0;
			//			Gxip[1] = 2.0 * (cxip[1]-eta0)/deta - 1.0;
			// x*-x*0/dx = (x+1)/2
			Gxip[0] = dxi   * (cxip[0]+1.0)/2.0 + xi0;
			Gxip[1] = deta  * (cxip[1]+1.0)/2.0 + eta0;
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"init nI,nJ = %D %D \n", nI,nJ );
			PetscPrintf(PETSC_COMM_SELF,"[Lxi-init] = %1.8e %1.8e \n", cxip[0], cxip[1] );
			PetscPrintf(PETSC_COMM_SELF,"[Gxi-init] = %1.8e %1.8e \n", Gxip[0], Gxip[1] );
#endif
		}
		if (monitor) {
			PetscPrintf(PETSC_COMM_SELF,"point[%d]: pos = ( %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e ) \n", 
						 p, marker_p->coor[0],marker_p->coor[1],
						 marker_p->xi[0], marker_p->xi[1] );
		}
		
		point_found = PETSC_FALSE;
		
		its = 0;
		do {
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"iteration: %D\n",its);
#endif
			/* convert Gxi to nInJ */
			nI = (Gxip[0]+1.0)/dxi;
			nJ = (Gxip[1]+1.0)/deta;
			
			if (nI == mx) { nI--; }
			if (nJ == my) { nJ--; }
			
			if ((nI<0) || (nJ<0)) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nI(%D),nJ(%D) negative Gxip %1.8e,%1.8e \n",nI,nJ,Gxip[0],Gxip[1]);
#endif
				break;
			}
			if (nI >= mx) { 
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nI too large \n");
#endif
				break;
			}
			if (nJ >= my) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nJ too large \n");
#endif
				break;
			}
			
			/* Get coords of cell nInJ */
			wil_nInJ = nI + nJ * mx;
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  nI,nJ=%D/%D : wil_nInJ %D : nid = ", nI,nJ,wil_nInJ);
#endif
			for (k=0; k<Q2_NODES_PER_EL_2D; k++) {
				PetscInt nid = element[wil_nInJ*Q2_NODES_PER_EL_2D+k];
				
				vertex[NSD2d*k+0] = coords[NSD2d*nid+0];
				vertex[NSD2d*k+1] = coords[NSD2d*nid+1];
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"%D ", nid);
#endif
			}
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"\n");
#endif
			
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  [vertex] ");
			for (k=0; k<Q2_NODES_PER_EL_2D; k++) {
				PetscPrintf(PETSC_COMM_SELF,"(%1.8e , %1.8e) ",vertex[NSD2d*k+0],vertex[NSD2d*k+1] );
			}
			PetscPrintf(PETSC_COMM_SELF,"\n");
#endif
			
			/* convert global (domain) xi TO local (element) xi  */
			xi0   = -1.0 + nI*dxi;
			eta0  = -1.0 + nJ*deta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Lxip[0] = 0.5*(Gxip[0]+1.0)*dxi  + xi0;
			//			Lxip[1] = 0.5*(Gxip[1]+1.0)*deta + eta0;
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0  )/dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0 )/deta  - 1.0;
			
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  Lxi,Lxeta = %1.8e, %1.8e (%D,%D) \n", Lxip[0],Lxip[1],nI,nJ );
#endif
			
			_compute_F_2dQ2( Lxip, vertex, marker_p->coor, f );
			if (monitor) {
				PetscPrintf(PETSC_COMM_SELF,"%4D InverseMapping : F = ( %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e ) \n", 
							 its, f[0],f[1], Lxip[0],Lxip[1] );
			}
			
			/* Check for convergence */
			F2 = (f[0]*f[0]+f[1]*f[1]);
			if (F2 < tolerance2) {
				if (monitor) PetscPrintf(PETSC_COMM_SELF,"%4D InverseMapping : converged : Norm of F %1.8e \n", its, sqrt(F2) );
				point_found = PETSC_TRUE;
				break;
			}
			
			_compute_J_2dQ2( Lxip, vertex, Jacobian );
			
			/* compute update */
			_compute_deltaX_2d( Jacobian, f, h );
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  [delta] = %1.8e %1.8e \n", h[0],h[1] );
#endif
			
			/* update Lxip */
			Lxip[0] += 10.0e-1 *h[0];
			Lxip[1] += 10.0e-1 *h[1];
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  [corrected] Lxi,Lxeta = %1.8e, %1.8e \n", Lxip[0],Lxip[1] );
#endif
			
			residual2 = ( h[0]*h[0] + h[1]*h[1] );
			if (residual2 < tolerance2) {
				if (monitor) PetscPrintf(PETSC_COMM_SELF,"%4D InverseMapping : converged : Norm of correction %1.8e \n", its, sqrt(residual2) );
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
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  [Gxi] = %1.8e %1.8e \n", Gxip[0], Gxip[1] );
#endif
			
			for (d=0; d<NSD2d; d++) {
				if (Gxip[d] < -1.0) { 
					Gxip[d] = -1.0;
#ifdef PTAT3D_DBG_PointLocation
					PetscPrintf(PETSC_COMM_SELF,"  correction outside box: correcting \n");
#endif
				}
				
				if (Gxip[d] > 1.0) {
					Gxip[d] = 1.0;
#ifdef PTAT3D_DBG_PointLocation
					PetscPrintf(PETSC_COMM_SELF,"  correction outside box: correcting \n");
#endif
				}
			}
			
			its++;
		} while (its < max_its);
		
		if (monitor && (point_found == PETSC_FALSE)){
			if (its >= max_its) {
				PetscPrintf(PETSC_COMM_SELF,"%4D %s : Reached maximum iterations (%D) without converging. \n", its, __FUNCTION__, max_its );
			} else {
				PetscPrintf(PETSC_COMM_SELF,"%4D %s : Newton broke down, diverged or stagnated after (%D) iterations without converging. \n", its, __FUNCTION__, its );
			}
		}
		
		/* if at the end of the solve, it still looks like the point is outside the mapped domain, mark point as not being found */
		if (fabs(Gxip[0]) > 1.0) { point_found = PETSC_FALSE; }
		if (fabs(Gxip[1]) > 1.0) { point_found = PETSC_FALSE; }
		
		/* update local variables */
		if (point_found == PETSC_FALSE) {
			Lxip[0] = NAN;
			Lxip[1] = NAN;
			wil_nInJ  = -1;
		}	else {
			/* convert Gxi to nInJ */
			nI = (Gxip[0]+1.0)/dxi;
			nJ = (Gxip[1]+1.0)/deta;
			if (nI == mx){ nI--; }
			if (nJ == my){ nJ--; }
			
			if (nI >= mx) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nI too large \n");
#endif
				break;
			}
			if (nJ >= my) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nJ too large \n");
#endif
				break;
			}
			
			/* convert global (domain) xi TO local (element) xi  */
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			xi0   = -1.0 + nI*dxi;
			eta0  = -1.0 + nJ*deta;
			
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0)  /dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0) /deta  - 1.0;
			
			wil_nInJ = nI + nJ * mx;
		}
		
		/* set into vector */
		marker_p->xi[0] = Lxip[0];
		marker_p->xi[1] = Lxip[1];
		marker_p->wil   = wil_nInJ;
	}
}

/* 3d implementation */
void LSFDeriv3dQ2(PetscReal _xi[],PetscReal **GNi)
{
	PetscReal basis_NI[NSD3d][NSD3d];
	PetscReal basis_GNI[NSD3d][NSD3d];
	PetscInt  i,j,k,d,cnt;
	
	
	for (d=0; d<NSD3d; d++) {
		PetscReal xi = _xi[d];
		
		basis_NI[d][0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )
		
		basis_GNI[d][0] = 0.5 * ( 2.0*xi - 1.0 );
		basis_GNI[d][1] = - 2.0*xi;
		basis_GNI[d][2] = 0.5 * ( 2.0*xi + 1.0 );
	}
	
	cnt = 0;
	for (k=0; k<NSD3d; k++) {
		for (j=0; j<NSD3d; j++) {
			for (i=0; i<NSD3d; i++) {
				
				GNi[0][cnt] = basis_GNI[0][i]  *  basis_NI[1][j]  *  basis_NI[2][k];
				GNi[1][cnt] = basis_NI[0][i]   *  basis_GNI[1][j] *  basis_NI[2][k];
				GNi[2][cnt] = basis_NI[0][i]   *  basis_NI[1][j]  *  basis_GNI[2][k];
				
				cnt++;
			}
		}
	}
}

void LSF3dQ2_CheckPartitionOfUnity(PetscReal xi[],PetscReal *val)
{
	const PetscReal tol = 1.0e-6;
	PetscReal       Ni[Q2_NODES_PER_EL_3D],sum;
	PetscInt        i;
	
    
	P3D_ConstructNi_Q2_3D(xi,Ni);
	sum = 0.0;
	for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
		sum += Ni[i];
	}
	*val = sum;
	if (fabs(sum-1.0) > tol) {
		PetscPrintf(PETSC_COMM_SELF,"**** sum( N_i(xi,eta) ) = %1.8e > 1.0 ==> partition of unity is not satisified, point is likely outside of the element ****\n",sum);
	}
}

void LSF3dQ2_CheckGlobalCoordinate(PetscReal element_coord[],PetscReal xi[],PetscReal xp[],PetscReal err[])
{
	const PetscReal tol = 1.0e-6;
	PetscReal       Ni[Q2_NODES_PER_EL_3D],xp_interp[NSD3d];
	PetscInt        i;
	
    
	P3D_ConstructNi_Q2_3D(xi,Ni);
	xp_interp[0] = 0.0;
	xp_interp[1] = 0.0;
	xp_interp[2] = 0.0;
	
	for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
		xp_interp[0] += Ni[i] * element_coord[NSD3d*i + 0];
		xp_interp[1] += Ni[i] * element_coord[NSD3d*i + 1];
		xp_interp[2] += Ni[i] * element_coord[NSD3d*i + 2];
	}
	err[0] = (xp[0] - xp_interp[0]);
	err[1] = (xp[1] - xp_interp[1]);
	err[2] = (xp[2] - xp_interp[2]);
	
	if (fabs(err[0]) > tol) {
		PetscPrintf(PETSC_COMM_SELF,"**** |xp - N_i(xi,eta).x_i| = %1.8e > %1.8e ==> x: coordinate interpolation error occurred ****\n",fabs(err[0]),tol);
	}
	if (fabs(err[1]) > tol) {
		PetscPrintf(PETSC_COMM_SELF,"**** |yp - N_i(xi,eta).y_i| = %1.8e > %1.8e ==> y: coordinate interpolation error occurred ****\n",fabs(err[1]),tol);
	}
	if (fabs(err[2]) > tol) {
		PetscPrintf(PETSC_COMM_SELF,"**** |zp - N_i(xi,eta).z_i| = %1.8e > %1.8e ==> z: coordinate interpolation error occurred ****\n",fabs(err[2]),tol);
	}
}

/* computes h = inv(J) f */
void _compute_deltaX(PetscReal A[NSD3d][NSD3d],PetscReal f[],PetscReal h[])
{
	PetscReal B[NSD3d][NSD3d];
	PetscReal t4, t6, t8, t10, t12, t14, t17;
	
    
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

void _compute_J_3dQ2(PetscReal xi[],PetscReal vertex[],PetscReal J[NSD3d][NSD3d])
{
	PetscInt  i,j;
	PetscReal GNi[NSD3d][Q2_NODES_PER_EL_3D];
	
	
    for (i=0; i<NSD3d; i++) {
		for (j=0; j<NSD3d; j++) {
			J[i][j] = 0.0;
		}
	}
	
	P3D_ConstructGNi_Q2_3D(xi,GNi);
	for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
		int i3 = i * NSD3d;
		PetscReal x = vertex[i3];
		PetscReal y = vertex[i3+1];
		PetscReal z = vertex[i3+2];
		
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

void _compute_F_3dQ2(PetscReal xi[],PetscReal vertex[],PetscReal pos[],PetscReal f[])
{
	PetscInt  i;
	PetscReal Ni[Q2_NODES_PER_EL_3D];
	
    
	/* Update F for the next iteration */
	f[0] = f[1] = f[2] = 0.0;
	
	P3D_ConstructNi_Q2_3D(xi,Ni);
	for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
		PetscInt i3   = i * NSD3d;
		PetscInt i3p1 = i3+1;
		PetscInt i3p2 = i3+2;
		
		f[0] += vertex[i3  ] * Ni[i];
		f[1] += vertex[i3p1] * Ni[i];
		f[2] += vertex[i3p2] * Ni[i];
	}
	f[0] = - f[0] + pos[0];
	f[1] = - f[1] + pos[1];
	f[2] = - f[2] + pos[2];
}

void InverseMappingDomain_3dQ2(PetscReal tolerance,PetscInt max_its,
                               PetscBool use_nonzero_guess,
                               PetscBool monitor,
                               const PetscReal coords[],const PetscInt mx,const PetscInt my,const PetscInt mz,const PetscInt element[],
                               int np,MPntStd marker[] )
{
	PetscReal h[NSD];
	PetscReal Jacobian[NSD][NSD];
	PetscReal f[NSD];
	int       p;
	PetscInt  its;
	PetscReal residual2,tolerance2,F2;
	PetscReal cxip[NSD3d],Lxip[NSD3d],Gxip[NSD3d];
	PetscReal dxi,deta,dzeta,xi0,eta0,zeta0;
	PetscInt  nI,nJ,nK,wil_nInJ,wil_2d,k;
	PetscReal vertex[NSD * Q2_NODES_PER_EL_3D];
	PetscBool point_found;
    
	tolerance2 = tolerance * tolerance; /* Eliminates the need to do a sqrt in the convergence test */
    
#ifdef PTAT3D_DBG_PointLocation
	PetscPrintf(PETSC_COMM_SELF,"Domain: ncells = %D x %D x %D = %D \n", mx,my,mz,mx*my*mz );
#endif
	
	/* map domain to [-1,1]x[-1,1]x[-1,1] domain */
	dxi   = 2.0/((PetscReal)mx);
	deta  = 2.0/((PetscReal)my);
	dzeta = 2.0/((PetscReal)mz);
#ifdef PTAT3D_DBG_PointLocation
	PetscPrintf(PETSC_COMM_SELF,"Domain: (dxi,eta,zeta) = (%1.8e,%1.8e,%1.8e)\n",dxi,deta,dzeta );
#endif
	
	for (p=0; p<np; p++) {
		MPntStd *marker_p = &marker[p];
		
#ifdef PTAT3D_DBG_PointLocation
		PetscPrintf(PETSC_COMM_SELF,"POINT[%d]\n", p );
#endif
		/* copy these values */
		cxip[0] = marker_p->xi[0];
		cxip[1] = marker_p->xi[1];
		cxip[2] = marker_p->xi[2];
		
		/* Check for an initial guess initial guess */
		if (use_nonzero_guess == PETSC_FALSE) {
			Gxip[0] = 0.0;
			Gxip[1] = 0.0;
			Gxip[2] = 0.0;
		} else {
			/* convert wil => nInJ */
			wil_nInJ = marker_p->wil;
			nK = wil_nInJ/(mx*my);
			wil_2d = wil_nInJ - nK*mx*my;
			nJ = wil_2d/mx;
			nI = wil_2d - nJ*mx;
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"init nI,nJ,nK = %D %D %D [wil=%D]\n", nI,nJ,nK,wil_nInJ );
#endif
			/* convert Lxip => Gxip */
			xi0   = -1.0 + nI*dxi;
			eta0  = -1.0 + nJ*deta;
			zeta0 = -1.0 + nK*dzeta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Gxip[0] = 2.0 * (cxip[0]-xi0)/dxi - 1.0;
			//			Gxip[1] = 2.0 * (cxip[1]-eta0)/deta - 1.0;
			// x*-x*0/dx = (x+1)/2
			Gxip[0] = dxi   * (cxip[0]+1.0)/2.0 + xi0;
			Gxip[1] = deta  * (cxip[1]+1.0)/2.0 + eta0;
			Gxip[2] = dzeta * (cxip[2]+1.0)/2.0 + zeta0;
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"[Lxi-init] = %1.8e %1.8e %1.8e\n", cxip[0], cxip[1], cxip[2] );
			PetscPrintf(PETSC_COMM_SELF,"[Gxi-init] = %1.8e %1.8e %1.8e\n", Gxip[0], Gxip[1], Gxip[2] );
#endif
			
			/* check */
#ifdef PTAT3D_DBG_PointLocation
			{
				PetscReal err[NSD3d];
                
				for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
					PetscInt nid = element[wil_nInJ*Q2_NODES_PER_EL_3D+k];
					
					vertex[NSD3d*k+0] = coords[NSD3d*nid+0];
					vertex[NSD3d*k+1] = coords[NSD3d*nid+1];
					vertex[NSD3d*k+2] = coords[NSD3d*nid+2];
				}
				LSF3dQ2_CheckGlobalCoordinate(vertex,cxip,marker_p->coor,err);
				PetscPrintf(PETSC_COMM_SELF,"  interpolated coord err %1.4e %1.4e %1.4e \n",err[0],err[1],err[2]);
			}
#endif
		}
		if (monitor) PetscPrintf(PETSC_COMM_SELF,"point[%d]: pos = ( %+1.8e, %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e, %+1.8e ) \n", p, marker_p->coor[0],marker_p->coor[1],marker_p->coor[2], marker_p->xi[0],marker_p->xi[1],marker_p->xi[2] );
		
		point_found = PETSC_FALSE;
		
		its = 0;
		do {
#ifdef PTAT3D_DBG_PointLocation
            PetscPrintf(PETSC_COMM_SELF,"iteration: %D\n",its);
#endif
			/* convert Gxi to nInJ */
			nI = (Gxip[0]+1.0)/dxi;
			nJ = (Gxip[1]+1.0)/deta;
			nK = (Gxip[2]+1.0)/dzeta;
			
			if (nI == mx) nI--;
			if (nJ == my) nJ--;
			if (nK == mz) nK--;
			
			if ( (nI<0) || (nJ<0)|| (nK<0) ) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nI(%D),nJ(%D),nK(%D) negative Gxip %1.8e,%1.8e,%1.8e \n",nI,nJ,nK,Gxip[0],Gxip[1],Gxip[2]);
#endif
				break;
			}
			if (nI >= mx) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nI too large \n");
#endif
				break;
			}
			if (nJ >= my) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nJ too large \n");
#endif
				break;
			}
			if (nK >= mz) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nK too large \n");
#endif
				break;
			}
			
			/* Get coords of cell nInJ */
			wil_nInJ = nI + nJ*mx + nK*mx*my;
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  nI,nJ,nK=%D/%D/%D : wil_nInJ %D : nid = ", nI,nJ,nK,wil_nInJ);
#endif
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
				PetscInt nid = element[wil_nInJ*Q2_NODES_PER_EL_3D+k];
				
				vertex[NSD3d*k+0] = coords[NSD3d*nid+0];
				vertex[NSD3d*k+1] = coords[NSD3d*nid+1];
				vertex[NSD3d*k+2] = coords[NSD3d*nid+2];
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"%D ", nid);
#endif
			}
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"\n");
#endif
#ifdef PTAT3D_DBG_PointLocation
			{
				PetscPrintf(PETSC_COMM_SELF,"  [vertex] ");
				for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
					PetscPrintf(PETSC_COMM_SELF,"(%1.8e , %1.8e , %1.8e) ",vertex[3*k+0],vertex[3*k+1],vertex[3*k+2] );
				}
				PetscPrintf(PETSC_COMM_SELF,"\n");
			}
#endif
			/* convert global (domain) xi TO local (element) xi  */
			xi0   = -1.0 + nI*dxi;
			eta0  = -1.0 + nJ*deta;
			zeta0 = -1.0 + nK*dzeta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Lxip[0] = 0.5*(Gxip[0]+1.0)*dxi  + xi0;
			//			Lxip[1] = 0.5*(Gxip[1]+1.0)*deta + eta0;
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0  )/dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0 )/deta  - 1.0;
			Lxip[2] = 2.0*(Gxip[2]-zeta0)/dzeta - 1.0;
			
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  Lxi,Lxeta,Lzeta = %1.8e, %1.8e, %1.8e (%D,%D,%D) \n", Lxip[0],Lxip[1],Lxip[2],nI,nJ,nK );
#endif
			_compute_F_3dQ2( Lxip, vertex, marker_p->coor, f );
#ifdef PTAT3D_DBG_PointLocation
			if (monitor) PetscPrintf(PETSC_COMM_SELF,"%4D InverseMapping : F = ( %+1.8e, %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e, %+1.8e ) \n", its, f[0],f[1],f[2], Lxip[0],Lxip[1],Lxip[2] );
#endif
			/* Check for convergence */
			F2 = (f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
			if (F2 < tolerance2) {
				if (monitor) PetscPrintf(PETSC_COMM_SELF,"%4D InverseMapping : converged : Norm of F %1.8e \n", its, sqrt(F2) );
				point_found = PETSC_TRUE;
				break;
			}
			
			_compute_J_3dQ2( Lxip, vertex, Jacobian );
			
			/* compute update */
			_compute_deltaX( Jacobian, f, h );
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  [delta] = %1.8e %1.8e %1.8e \n", h[0],h[1],h[2] );
#endif
			/* update Lxip */
			Lxip[0] += 10.0e-1 *h[0];
			Lxip[1] += 10.0e-1 *h[1];
			Lxip[2] += 10.0e-1 *h[2];
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  [corrected] Lxi,Leta,Lzeta = %1.8e, %1.8e, %1.8e \n", Lxip[0],Lxip[1],Lxip[2] );
#endif
			residual2 = ( h[0]*h[0] + h[1]*h[1] + h[2]*h[2] );
			if (residual2 < tolerance2) {
				if (monitor) PetscPrintf(PETSC_COMM_SELF,"%4D InverseMapping : converged : Norm of correction %1.8e \n", its, sqrt(residual2) );
				point_found = PETSC_TRUE;
				break;
			}
			
			/* convert Lxip => Gxip */
			//			xi0  = -1.0 + nI*dxi;
			//			eta0 = -1.0 + nJ*deta;
			
			//			Gxip[0] = 2.0 * (Lxip[0]-xi0)/dxi   - 1.0;
			//			Gxip[1] = 2.0 * (Lxip[1]-eta0)/deta - 1.0;
			// x*-x*0/dx = (x+1)/2
			Gxip[0] = dxi   * (Lxip[0]+1.0)/2.0 + xi0;
			Gxip[1] = deta  * (Lxip[1]+1.0)/2.0 + eta0;
			Gxip[2] = dzeta * (Lxip[2]+1.0)/2.0 + zeta0;
#ifdef PTAT3D_DBG_PointLocation
			PetscPrintf(PETSC_COMM_SELF,"  [Gxi] = %1.8e %1.8e %1.8e\n", Gxip[0], Gxip[1], Gxip[2] );
#endif
			if (Gxip[0] < -1.0) {
				Gxip[0] = -1.0;
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  correction outside box: correcting \n");
#endif
			}
			if (Gxip[1] < -1.0) {
				Gxip[1] = -1.0;
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  correction outside box: correcting \n");
#endif
			}
			if (Gxip[2] < -1.0) {
				Gxip[2] = -1.0;
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  correction outside box: correcting \n");
#endif
			}
			
			if (Gxip[0] > 1.0) {
				Gxip[0] = 1.0;
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  correction outside box: correcting \n");
#endif
			}
			if (Gxip[1] > 1.0) {
				Gxip[1] = 1.0;
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  correction outside box: correcting \n");
#endif
			}
			if (Gxip[2] > 1.0) {
				Gxip[2] = 1.0;
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  correction outside box: correcting \n");
#endif
			}
			
			its++;
		} while (its < max_its);
		
		if (monitor && (point_found == PETSC_FALSE)) {
			if (its >= max_its) {
				PetscPrintf(PETSC_COMM_SELF,"%4D %s : Reached maximum iterations (%D) without converging. \n", its, __FUNCTION__, max_its );
			} else {
				PetscPrintf(PETSC_COMM_SELF,"%4D %s : Newton broke down, diverged or stagnated after (%D) iterations without converging. \n", its, __FUNCTION__, its );
			}
		}
		
		/* if at the end of the solve, it still looks like the point is outside the mapped domain, mark point as not being found */
		if (fabs(Gxip[0]) > 1.0) { point_found =PETSC_FALSE; }
		if (fabs(Gxip[1]) > 1.0) { point_found =PETSC_FALSE; }
		if (fabs(Gxip[2]) > 1.0) { point_found =PETSC_FALSE; }
		
		/* update local variables */
		if (point_found == PETSC_FALSE) {
			Lxip[0] = NAN;
			Lxip[1] = NAN;
			Lxip[2] = NAN;
			wil_nInJ  = -1;
		} else {
			/* convert Gxi to nInJ */
			nI = (Gxip[0]+1.0)/dxi;
			nJ = (Gxip[1]+1.0)/deta;
			nK = (Gxip[2]+1.0)/dzeta;
			if (nI == mx) nI--;
			if (nJ == my) nJ--;
			if (nK == mz) nK--;
			
			if (nI >= mx) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nI too large \n");
#endif
				break;
			}
			if (nJ >= my) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nJ too large \n");
#endif
				break;
			}
			if (nK >= mz) {
#ifdef PTAT3D_DBG_PointLocation
				PetscPrintf(PETSC_COMM_SELF,"  nK too large \n");
#endif
				break;
			}

			/* convert global (domain) xi TO local (element) xi  */
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			xi0   = -1.0 + nI*dxi;
			eta0  = -1.0 + nJ*deta;
			zeta0 = -1.0 + nK*dzeta;
			
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0  )/dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0 )/deta  - 1.0;
			Lxip[2] = 2.0*(Gxip[2]-zeta0)/dzeta - 1.0;
			
			wil_nInJ = nI + nJ * mx + nK * mx * my;
		}
		
		/* set into vector */
		marker_p->xi[0] = Lxip[0];
		marker_p->xi[1] = Lxip[1];
		marker_p->xi[2] = Lxip[2];
		marker_p->wil   = wil_nInJ;
        
#ifdef PTAT3D_DBG_PointLocation
		PetscPrintf(PETSC_COMM_SELF,"  <<final>> xi,eta,zeta = %1.8e, %1.8e, %1.8e [wil=%D] \n", Lxip[0],Lxip[1],Lxip[2],wil_nInJ);
#endif
	}
}

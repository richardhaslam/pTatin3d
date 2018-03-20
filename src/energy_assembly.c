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
 **    filename:   energy_assembly.c
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

#include "petsc.h"
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_bcs.h"
#include "element_utils_q1.h"
#include "dmda_element_q1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"

#include "QPntVolCoefEnergy_def.h"
#include "phys_comp_energy.h"
#include "ptatin3d_energy.h"
#include "energy_assembly.h"
#include "energy_coefficients.h"

/*
 
 Advection-diffusion implementation for ALE mesh.
 
 We solve
 
 dS/dt + u.grad(S) = div( kappa grad(S) ) + f
 
 A stabilized formulation (SUPG) is used to discretise the above equation.
 We adopt a backward Euler time discretisation.
 Backward Euler is adopted as it doesn't oscillate, in contrast with Crank-Nicolson.
 
 (1/dt) ( M(x^k+1) T^k+1 - M(x^k) T^k ) + C(x^k+1) = K(x^k+1) T^k+1 + F(x^k+1)
 
 ( M(x^k+1) + dt C(x^k+1) - dt K(x^k+1) ) T^k+1 = dt F(x^k+1) + M(x^k) T^k
 
 We indicate the discrete operators are evaluated using particular coordinates
 which are time depenedent.
 x^k+1 are the updated coords, x^k are the previous coordinates.
 
 In non-linear residual form, we will have
 
 F := ( M(x^k+1) + dt C(x^k+1) - dt K(x^k+1) ) T^k+1 - dt F(x^k+1) - M(x^k) T^k
 
 
 
 LLP on SUPG:
 "Well it fucked up everythng like love..., always too fast or too slow never on time...,
 always trying to make tiny correction to make the result better...,
 but the correction are never in the right direction so it creates instabilities."
 - February 14, 2013 (Valentines day)
 
 */

#define SUPG_STAB 1
#define TAU_STAB  2

#define ADV_DIFF_STABILIZATION_TYPE SUPG_STAB
//#define ADV_DIFF_STABILIZATION_TYPE TAU_STAB

#define ADV_DIFF_STAB_EPS 1.0e-10

/* SUPG business */
/** AdvectionDiffusion_UpwindXiExact - Brooks, Hughes 1982 equation 2.4.2
 * \bar \xi = coth( \alpha ) - \frac{1}{\alpha} */
double AdvDiffResidualForceTerm_UpwindXiExact(double pecletNumber)
{
    if (fabs(pecletNumber) < 1.0e-8 ) {
        return 0.33333333333333 * pecletNumber;
    } else if (pecletNumber < -20.0) {
        return -1.0 - 1.0/pecletNumber;
    } else if (pecletNumber > 20.0) {
        return +1.0 - 1.0/pecletNumber;
    }
    return cosh( pecletNumber )/sinh( pecletNumber ) - 1.0/pecletNumber;
}

/** AdvectionDiffusion_UpwindXiDoublyAsymptoticAssumption - Brooks, Hughes 1982 equation 3.3.1
 * Simplification of \bar \xi = coth( \alpha ) - \frac{1}{\alpha} from Brooks, Hughes 1982 equation 2.4.2
 *
 *            { -1               for \alpha <= -3
 * \bar \xi ~ { \frac{/alpha}{3} for -3 < \alpha <= 3
 *            { +1               for \alpha > +3 */
double AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption(double pecletNumber)
{
    if (pecletNumber <= -3.0) {
        return -1;
    } else if (pecletNumber <= 3.0) {
        return 0.33333333333333 * pecletNumber;
    } else {
        return 1.0;
    }
}

/** AdvectionDiffusion_UpwindXiCriticalAssumption - Brooks, Hughes 1982 equation 3.3.2
 * Simplification of \bar \xi = coth( \alpha ) - \frac{1}{\alpha} from Brooks, Hughes 1982 equation 2.4.2
 *
 *            { -1 - \frac{1}{\alpha}   for \alpha <= -1
 * \bar \xi ~ { 0                       for -1 < \alpha <= +1
 *            { +1 - \frac{1}{\alpha}   for \alpha > +1             */
double AdvDiffResidualForceTerm_UpwindXiCriticalAssumption(double pecletNumber)
{
    if (pecletNumber <= -1.0) {
        return -1.0 - 1.0/pecletNumber;
    } else if (pecletNumber <= 1.0) {
        return 0.0;
    } else {
        return 1.0 - 1.0/pecletNumber;
    }
}

#undef __FUNCT__
#define __FUNCT__ "AdvDiff3dComputeAverageCellSize"
PetscErrorCode AdvDiff3dComputeAverageCellSize(PetscScalar el_coords[],PetscScalar DX[])
{
	PetscInt d,k;
	PetscReal min_x[NSD],max_x[NSD];
	PetscReal xc,yc,zc;
	
	PetscFunctionBegin;
	for (d=0; d<NSD; d++) {
		min_x[d] = 1.0e32;
		max_x[d] = -1.0e32;
	}
	
	for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		xc = el_coords[NSD*k+0];
		yc = el_coords[NSD*k+1];
		zc = el_coords[NSD*k+2];
		
		if (xc < min_x[0]) { min_x[0] = xc; }
		if (yc < min_x[1]) { min_x[1] = yc; }
		if (zc < min_x[2]) { min_x[2] = zc; }
		
		if (xc > max_x[0]) { max_x[0] = xc; }
		if (yc > max_x[1]) { max_x[1] = yc; }
		if (zc > max_x[2]) { max_x[2] = zc; }
	}
	
	for (d=0; d<NSD; d++) {
		DX[d] = max_x[d] - min_x[d];
	}
	PetscFunctionReturn(0);
}

/* Eqn 4.3.7 */
#undef __FUNCT__
#define __FUNCT__ "AdvDiff3dComputeElementPecletNumber_qp"
PetscErrorCode AdvDiff3dComputeElementPecletNumber_qp( PetscScalar el_coords[],PetscScalar u[],
                                                      PetscScalar kappa_el,
                                                      PetscScalar *alpha)
{
	PetscInt    d,k;
	PetscScalar DX[NSD];
    PetscScalar u_xi[NSD],one_dxi2[NSD];
    PetscScalar _alpha,u_norm,dxi[NSD];
	PetscErrorCode ierr;
	
    PetscFunctionBegin;
	
	ierr = AdvDiff3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
	
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
    for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] += 0.125 *  u[NSD*k + 0];
		u_xi[1] += 0.125 *  u[NSD*k + 1];
		u_xi[2] += 0.125 *  u[NSD*k + 2];
	}
	
	for (d=0; d<NSD; d++) {
		dxi[d]      = DX[d];
		one_dxi2[d] = 1.0/dxi[d]/dxi[d];
	}
	
    u_norm = sqrt( u_xi[0]*u_xi[0] + u_xi[1]*u_xi[1] + u_xi[2]*u_xi[2]);
	
	_alpha = 0.5 * u_norm;
	if (kappa_el > ADV_DIFF_STAB_EPS) {
		_alpha = _alpha / ( (kappa_el) * sqrt( one_dxi2[0] + one_dxi2[1] + one_dxi2[2]) );
	}
    *alpha  = _alpha;
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DASUPG3dComputeElementTimestep_qp"
PetscErrorCode DASUPG3dComputeElementTimestep_qp(PetscScalar el_coords[],PetscScalar u[],PetscReal kappa_el,PetscReal *dta,PetscReal *dtd)
{
	PetscInt    d,k;
	PetscScalar DX[NSD];
    PetscScalar u_xi[NSD],one_dxi2[NSD],dxi[NSD];
    PetscScalar CrFAC,dt_optimal,alpha;
    PetscScalar U,H,dt_diffusive;
	PetscErrorCode ierr;
	
    PetscFunctionBegin;
	
	ierr = AdvDiff3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
    
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
    for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] += 0.125 *  u[NSD*k + 0];
		u_xi[1] += 0.125 *  u[NSD*k + 1];
		u_xi[2] += 0.125 *  u[NSD*k + 2];
	}
	
	for (d=0; d<NSD; d++) {
		dxi[d]      = DX[d];
		one_dxi2[d] = 1.0/dxi[d]/dxi[d];
	}
	
    U = sqrt( u_xi[0]*u_xi[0] + u_xi[1]*u_xi[1] + u_xi[2]*u_xi[2] );
    H = 1.0 / sqrt( one_dxi2[0] + one_dxi2[1] + one_dxi2[2] );
	
    ierr = AdvDiff3dComputeElementPecletNumber_qp(el_coords,u,kappa_el,&alpha);CHKERRQ(ierr);
	
    if (U < ADV_DIFF_STAB_EPS) {
        dt_diffusive = 0.5 * H*H / kappa_el;
        *dtd = dt_diffusive;
        *dta = 1.0e10;
    }
	else if (kappa_el < ADV_DIFF_STAB_EPS) {
        if (alpha >= 100.0 ) {
            CrFAC = 0.4;
        } else {
            CrFAC = PetscMin(1.0,alpha);
        }
        dt_optimal = CrFAC * H / (U+ADV_DIFF_STAB_EPS);
        *dta = dt_optimal;
        *dtd = 1.0e10;
    } else {
        dt_diffusive = 0.5 * H*H / kappa_el;
        *dtd = dt_diffusive;
		
        if (alpha >= 100.0 ) {
            CrFAC = 0.4;
        } else {
            CrFAC = PetscMin(1.0,alpha);
        }
        dt_optimal = CrFAC * H / (U+ADV_DIFF_STAB_EPS);
        *dta = dt_optimal;
	}
	/*
	 printf("Element Pe:              %1.4f \n", alpha );
	 printf("        dx/dv:           %1.4f \n", H/(U+1.0e-32) );
	 printf("        0.5.dx.dx/kappa: %1.4f \n", 0.5*H*H/(kappa_el+1.0e-32) );
	 */
	
	//*dta = H/(U+ADV_DIFF_STAB_EPS);
	//*dtd = 0.5*H*H/(kappa_el+ADV_DIFF_STAB_EPS);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AdvDiff3dComputeElementTimestep_qp"
PetscErrorCode AdvDiff3dComputeElementTimestep_qp(PetscScalar el_coords[],PetscScalar u[],PetscReal kappa_el,PetscReal *dta,PetscReal *dtd)
{
	PetscInt    d,k;
	PetscScalar DX[NSD];
    PetscScalar u_xi[NSD],one_dxi2[NSD],dxi[NSD];
    PetscScalar U,H;
	PetscErrorCode ierr;
	
    PetscFunctionBegin;
	
	ierr = AdvDiff3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
    
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
    for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] += 0.125 *  u[NSD*k + 0];
		u_xi[1] += 0.125 *  u[NSD*k + 1];
		u_xi[2] += 0.125 *  u[NSD*k + 2];
	}
	
	for (d=0; d<NSD; d++) {
		dxi[d]      = DX[d];
		one_dxi2[d] = 1.0/dxi[d]/dxi[d];
	}
	
    U = sqrt( u_xi[0]*u_xi[0] + u_xi[1]*u_xi[1] + u_xi[2]*u_xi[2] );
    H = 1.0 / sqrt( one_dxi2[0] + one_dxi2[1] + one_dxi2[2] );
	
	*dta = H/(U+ADV_DIFF_STAB_EPS);
	*dtd = 0.5*H*H/(kappa_el+ADV_DIFF_STAB_EPS);
	
    PetscFunctionReturn(0);
}

/*
 Eqns 3.3.11 (transient) + 3.3.4, 3.3.5, 3.3.6
 */
#undef __FUNCT__
#define __FUNCT__ "DASUPG3dComputeElementStreamlineDiffusion_qp"
PetscErrorCode DASUPG3dComputeElementStreamlineDiffusion_qp(PetscScalar el_coords[],PetscScalar u[],
                                                            PetscInt nqp, PetscScalar qp_detJ[], PetscScalar qp_w[],
                                                            PetscScalar qp_kappa[],
                                                            PetscScalar *khat)
{
	PetscInt p,k,d;
	PetscScalar DX[NSD];
    PetscScalar u_xi[NSD];
    PetscScalar kappa_el,alpha_xi[NSD],_khat,dxi[NSD],xi[NSD],vol_el;
	PetscErrorCode ierr;
	
    PetscFunctionBegin;
	
	/* average velocity - cell centre velocity */
	ierr = AdvDiff3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
	
	/* average velocity - cell centre velocity */
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
    for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] += 0.125 *  u[NSD*k + 0];
		u_xi[1] += 0.125 *  u[NSD*k + 1];
		u_xi[2] += 0.125 *  u[NSD*k + 2];
	}
	
	/* integral average of diffusivity */
	kappa_el = 0.0;
	vol_el   = 0.0;
	for (p=0; p<nqp; p++) {
		kappa_el += qp_w[p] * qp_kappa[p] * qp_detJ[p];
		vol_el   += qp_w[p] * qp_detJ[p];
	}
    kappa_el = kappa_el / vol_el;
	
	/* could replace with / (1.0e-30 + kappa_el) */
	for (d=0; d<NSD; d++) {
		dxi[d]      = DX[d];
		if (kappa_el < ADV_DIFF_STAB_EPS) {
			alpha_xi[d] = 1.0/ADV_DIFF_STAB_EPS;
		} else {
			alpha_xi[d] = 0.5 * u_xi[d]  * dxi[d]  / kappa_el;
		}
		xi[d]       = AdvDiffResidualForceTerm_UpwindXiExact(alpha_xi[d]);
	}
	
	//  khat = 0.5*(xi*u_xi*dxi + eta*v_eta*deta); /* steady state */
    _khat = 1.0 * 0.258198889747161 * ( xi[0]*u_xi[0]*dxi[0] + xi[1]*u_xi[1]*dxi[1] + xi[2]*u_xi[2]*dxi[2] ); /* transient case, sqrt(1/15) */
    *khat = _khat;
	
    PetscFunctionReturn(0);
}

/* SUPG IMPLEMENTATION */
void ConstructNiSUPG_Q1_3D(PetscScalar Up[],PetscScalar kappa_hat,PetscScalar Ni[],PetscScalar GNx[NSD][NODES_PER_EL_Q1_3D],PetscScalar Ni_supg[])
{
    PetscScalar uhat[NSD],unorm;
    PetscInt i;
	
	uhat[0] = Up[0];
	uhat[1] = Up[1];
	uhat[2] = Up[2];
    unorm = PetscSqrtScalar(Up[0]*Up[0] + Up[1]*Up[1] + Up[2]*Up[2]);
	if (unorm > ADV_DIFF_STAB_EPS) {
		uhat[0] = Up[0]/unorm;
		uhat[1] = Up[1]/unorm;
		uhat[2] = Up[2]/unorm;
	}
	
    if (kappa_hat < ADV_DIFF_STAB_EPS) {
        for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
            Ni_supg[i] = Ni[i];
        }
    } else {
        for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
            Ni_supg[i] = Ni[i] + kappa_hat * ( uhat[0] * GNx[0][i] + uhat[1] * GNx[1][i] + uhat[2] * GNx[2][i] );
        }
    }
}

#undef __FUNCT__
#define __FUNCT__ "AElement_FormJacobian_T_supg"
PetscErrorCode AElement_FormJacobian_T_supg( PetscScalar Re[],PetscReal dt,PetscScalar el_coords[],
                                            PetscScalar gp_kappa[],
                                            PetscScalar el_V[],
                                            PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] )
{
    PetscInt    p,i,j;
    PetscReal   Ni_p[NODES_PER_EL_Q1_3D],Ni_supg_p[NODES_PER_EL_Q1_3D];
    PetscReal   GNi_p[NSD][NODES_PER_EL_Q1_3D],GNx_p[NSD][NODES_PER_EL_Q1_3D],gp_detJ[27];
    PetscScalar J_p,fac;
    PetscScalar kappa_p,v_p[NSD];
    PetscScalar kappa_hat;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
    /* compute constants for the element */
    for (p = 0; p < ngp; p++) {
        P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
        P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&gp_detJ[p],&GNx_p[0],&GNx_p[1],&GNx_p[2]);
	}
	
	ierr = DASUPG3dComputeElementStreamlineDiffusion_qp(el_coords,el_V,
                                                        ngp,gp_detJ,gp_weight,gp_kappa,
                                                        &kappa_hat);CHKERRQ(ierr);
    /* evaluate integral */
    for (p = 0; p < ngp; p++) {
        P3D_ConstructNi_Q1_3D(&gp_xi[NSD*p],Ni_p);
        P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
        P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&J_p,&GNx_p[0],&GNx_p[1],&GNx_p[2]);
		
        fac = gp_weight[p]*J_p;
		
		kappa_p = gp_kappa[p];
        v_p[0] = 0.0;
		v_p[1] = 0.0;
		v_p[2] = 0.0;
        for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
            v_p[0]     += Ni_p[j] * el_V[NSD*j+0];      /* compute vx on the particle */
            v_p[1]     += Ni_p[j] * el_V[NSD*j+1];      /* compute vy on the particle */
            v_p[2]     += Ni_p[j] * el_V[NSD*j+2];      /* compute vy on the particle */
        }
        ConstructNiSUPG_Q1_3D(v_p,kappa_hat,Ni_p,GNx_p,Ni_supg_p);
		
		// F := ( M(x^k+1) + dt C(x^k+1) - dt K(x^k+1) ) T^k+1 - dt F(x^k+1) - M(x^k) T^k
		// J = dF/d(T^k+1) = M(x^k+1) + dt C(x^k+1) - dt K(x^k+1)
		for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
			for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
				Re[j+i*NODES_PER_EL_Q1_3D] += fac * (
                                                     //Ni_p[i] * Ni_p[j]
                                                     Ni_supg_p[i] * Ni_p[j]
                                                     + dt * Ni_supg_p[i] * ( v_p[0]*GNx_p[0][j] + v_p[1]*GNx_p[1][j] + v_p[2]*GNx_p[2][j] )
                                                     + dt * kappa_p * ( GNx_p[0][i]*GNx_p[0][j] + GNx_p[1][i]*GNx_p[1][j] + GNx_p[2][i]*GNx_p[2][j] )
                                                     );
			}
		}
    }
    
	/*
	 printf("e=\n");
	 for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
	 for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
	 printf("%lf ", Re[j+i*NODES_PER_EL_Q1_3D]);
	 }printf("\n");
	 }
	 */
    
	PetscFunctionReturn(0);
}

/*
 Computes M + dt.(L + A)
 */
#undef __FUNCT__
#define __FUNCT__ "TS_FormJacobianEnergy"
PetscErrorCode TS_FormJacobianEnergy(PetscReal time,Vec X,PetscReal dt,Mat A,Mat B,void *ctx)
{
  pTatinCtx      ptatin = (pTatinCtx)ctx;
  PhysCompEnergy data;
	PetscInt          nqp;
	PetscScalar       *qp_xi,*qp_weight;
	Quadrature        volQ;
	QPntVolCoefEnergy *all_quadpoints,*cell_quadpoints;
	PetscScalar       qp_kappa[27];
    DM            da,cda;
	Vec           gcoords;
	PetscScalar   *LA_gcoords;
	PetscScalar   ADe[NODES_PER_EL_Q1_3D*NODES_PER_EL_Q1_3D];
	PetscScalar   el_coords[NSD*NODES_PER_EL_Q1_3D];
	PetscScalar   el_V[NSD*NODES_PER_EL_Q1_3D];
	/**/
	PetscInt       nel,nen,e,n,ii;
	const PetscInt *elnidx;
	BCList         bclist;
    ISLocalToGlobalMapping ltog;
	PetscInt       NUM_GINDICES,T_el_lidx[Q1_NODES_PER_EL_3D],ge_eqnums[Q1_NODES_PER_EL_3D];
	const PetscInt *GINDICES;
	Vec            V;
    Vec            local_V,local_X;
    PetscScalar    *LA_V,*LA_X;
	PetscBool      mat_mffd;
	PetscErrorCode ierr;
	
	
    PetscFunctionBegin;
	ierr = pTatinGetContext_Energy(ptatin,&data);CHKERRQ(ierr);
	da     = data->daT;
	V      = data->u_minus_V;
	bclist = data->T_bclist;
	volQ   = data->volQ;
	
    
	/* setup for coords */
  ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

	/* get solution */
  ierr = DMGetLocalVector(da,&local_X);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_X);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,local_X);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (da,X,INSERT_VALUES,local_X);CHKERRQ(ierr);
  ierr = VecGetArray(local_X,&LA_X);CHKERRQ(ierr);
  
  /* get acces to the vector V */
  ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
  ierr = DMGetLocalVector(cda,&local_V);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(cda,V,INSERT_VALUES,local_V);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(  cda,V,INSERT_VALUES,local_V);CHKERRQ(ierr);
  ierr = VecGetArray(local_V,&LA_V);CHKERRQ(ierr);

  /* ===================== */
  /* Evaluate coefficients */
  ierr = EnergyEvaluateCoefficients(ptatin,time,da,LA_X,LA_V);CHKERRQ(ierr);
  
  /* trash old entries */
	ierr = PetscObjectTypeCompare((PetscObject)A,MATMFFD,&mat_mffd);CHKERRQ(ierr);
  
	if (!mat_mffd) {
		ierr = MatZeroEntries(A);CHKERRQ(ierr);
	} else {
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
  
	if (B) { ierr = MatZeroEntries(B);CHKERRQ(ierr); }
  
	if (B == NULL) { PetscFunctionReturn(0); }
	
	

	/* quadrature */
	volQ      = data->volQ;
	nqp       = volQ->npoints;
	qp_xi     = volQ->q_xi_coor;
	qp_weight = volQ->q_weight;
	
	ierr = VolumeQuadratureGetAllCellData_Energy(volQ,&all_quadpoints);CHKERRQ(ierr);
	
  
	/* stuff for eqnums */
    ierr = DMGetLocalToGlobalMapping(da, &ltog);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetSize(ltog, &NUM_GINDICES);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetIndices(ltog, &GINDICES);CHKERRQ(ierr);
	ierr = BCListApplyDirichletMask(NUM_GINDICES,(PetscInt*)GINDICES,bclist);CHKERRQ(ierr);
	
	ierr = DMDAGetElementsQ1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		/* get coords for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		
		/* get velocity for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_V,(PetscInt*)&elnidx[nen*e],LA_V);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Energy(volQ,all_quadpoints,e,&cell_quadpoints);CHKERRQ(ierr);
		
		/* copy the diffusivity */
		for (n=0; n<nqp; n++) {
			qp_kappa[n] = cell_quadpoints[n].diffusivity;
		}
		
		/* initialise element stiffness matrix */
		ierr = PetscMemzero(ADe,sizeof(PetscScalar)*NODES_PER_EL_Q1_3D*NODES_PER_EL_Q1_3D);CHKERRQ(ierr);
		
		/* form element stiffness matrix */
#if (ADV_DIFF_STABILIZATION_TYPE == SUPG_STAB)
		ierr = AElement_FormJacobian_T_supg( ADe,dt,el_coords, qp_kappa, el_V, nqp,qp_xi,qp_weight );CHKERRQ(ierr);
#elif (ADV_DIFF_STABILIZATION_TYPE == TAU_STAB)
		{
			PetscReal tau;
			PetscReal theta = 1.0; /* 0.5 crank nicoloson; 1.0 backward difference */
			PetscReal kappa_cell;
			
			kappa_cell = 0.0;
			for (n=0; n<nqp; n++) {
				kappa_cell = kappa_cell + qp_kappa[n];
			}
			kappa_cell = kappa_cell / ((PetscReal)nqp);
			
			//ierr = AdvDiffComputeTau_BrooksHughes(el_coords,el_V,kappa_cell,&tau);CHKERRQ(ierr);
			//printf("AdvDiffComputeTau_BrooksHughes : tau = %1.4e\n",tau);
			
			ierr = AdvDiffComputeTau_TezduyarOsawa(el_coords,el_V,kappa_cell,theta,dt,&tau);CHKERRQ(ierr);
			//printf("AdvDiffComputeTau_TezduyarOsawa : tau = %1.4e\n",tau);
			
			//ierr = AdvDiffComputeTau_UserDefinedConstant(NULL,&tau);CHKERRQ(ierr);
			//printf("AdvDiffComputeTau_UserDefinedConstant : tau = %1.4e\n",tau);
			
			ierr = AElement_FormJacobian_T_tau( ADe,dt,tau,el_coords, qp_kappa, el_V, nqp,qp_xi,qp_weight );CHKERRQ(ierr);
		}
#endif
		
		/* get indices */
		ierr = DMDAEQ1_GetElementLocalIndicesDOF(T_el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		for (ii=0; ii<NODES_PER_EL_Q1_3D; ii++) {
			const PetscInt NID = elnidx[ NODES_PER_EL_Q1_3D * e + ii ];
            
			ge_eqnums[ii] = GINDICES[ NID ];
		}
		
		/* insert element matrix into global matrix */
		ierr = MatSetValues(B,NODES_PER_EL_Q1_3D,ge_eqnums, NODES_PER_EL_Q1_3D,ge_eqnums, ADe, ADD_VALUES );CHKERRQ(ierr);
    }
	/* tidy up */
  ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  ierr = VecRestoreArray(local_X,&LA_X);CHKERRQ(ierr);
  ierr = VecRestoreArray(local_V,&LA_V);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(cda,&local_V);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&local_X);CHKERRQ(ierr);
	
	/* partial assembly */
	ierr = MatAssemblyBegin(B, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(B, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = BCListRemoveDirichletMask(NUM_GINDICES,(PetscInt*)GINDICES,bclist);CHKERRQ(ierr);
	ierr = BCListInsertScaling(B,NUM_GINDICES,(PetscInt*)GINDICES,bclist);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingRestoreIndices(ltog, &GINDICES);CHKERRQ(ierr);
	
	/* assemble */
	ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    if (A != B) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNES_FormJacobianEnergy"
PetscErrorCode SNES_FormJacobianEnergy(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  pTatinCtx      ptatin = (pTatinCtx)ctx;
  PhysCompEnergy energy  = (PhysCompEnergy)ctx;
  PetscErrorCode ierr;

	PetscFunctionBegin;
	ierr = pTatinGetContext_Energy(ptatin,&energy);CHKERRQ(ierr);
	ierr = TS_FormJacobianEnergy(energy->time,X,energy->dt,A,B,ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AElement_FormFunction_T_supg"
PetscErrorCode AElement_FormFunction_T_supg(
                                            PetscScalar Re[],
                                            PetscReal dt,
                                            PetscScalar el_coords[],PetscScalar el_coords_old[],
                                            PetscScalar el_V[],
                                            PetscScalar el_phi[],PetscScalar el_phi_last[],
                                            PetscScalar gp_kappa[],PetscScalar gp_Q[],
                                            PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] )
{
    PetscInt    p,i,j;
    PetscScalar Ni_p[NODES_PER_EL_Q1_3D],Ni_supg_p[NODES_PER_EL_Q1_3D];
    PetscScalar GNi_p[NSD][NODES_PER_EL_Q1_3D];
	PetscScalar GNx_p[NSD][NODES_PER_EL_Q1_3D], GNx_p_old[NSD][NODES_PER_EL_Q1_3D];
    PetscScalar phi_p,phi_last_p,f_p,v_p[NSD],kappa_p, gradphi_p[NSD],gradphiold_p[NSD],M_dotT_p;
    PetscScalar J_p,J_p_old,fac,gp_detJ[27];
    PetscScalar kappa_hat;
	
	PetscFunctionBegin;
    /* compute constants for the element */
    for (p = 0; p < ngp; p++) {
		P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
		P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&gp_detJ[p],&GNx_p[0],&GNx_p[1],&GNx_p[2]);
	}
	
	DASUPG3dComputeElementStreamlineDiffusion_qp(el_coords,el_V,
                                                 ngp,gp_detJ,gp_weight,gp_kappa,
                                                 &kappa_hat);
	
    /* evaluate integral */
    for (p = 0; p < ngp; p++) {
        P3D_ConstructNi_Q1_3D(&gp_xi[NSD*p],Ni_p);
		P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
        
		P3D_evaluate_geometry_elementQ1(1,el_coords,    &GNi_p,&J_p,    &GNx_p[0],    &GNx_p[1],    &GNx_p[2]);
		P3D_evaluate_geometry_elementQ1(1,el_coords_old,&GNi_p,&J_p_old,&GNx_p_old[0],&GNx_p_old[1],&GNx_p_old[2]);
		
        fac     = gp_weight[p] * J_p;
        //fac_old = gp_weight[p] * J_p_old;
		//printf("p=%d: J_p %1.4e ; j_p_old %1.4e ; fac %1.4e ; fac_old %1.4e \n", p,J_p,J_p_old,fac,fac_old);
		
		kappa_p      = gp_kappa[p];
		f_p          = gp_Q[p];
        phi_p        = 0.0;
        phi_last_p   = 0.0;
        
        gradphi_p[0] = 0.0;
		gradphi_p[1] = 0.0;
		gradphi_p[2] = 0.0;
        
		gradphiold_p[0] = 0.0;
		gradphiold_p[1] = 0.0;
		gradphiold_p[2] = 0.0;
        
		v_p[0] = 0.0;
		v_p[1] = 0.0;
		v_p[2] = 0.0;
        
		for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
            phi_p      += Ni_p[j] * el_phi[j];      /* compute phi on the particle */
            phi_last_p += Ni_p[j] * el_phi_last[j];  /* compute phi_dot on the particle */
            
            gradphiold_p[0] += GNx_p[0][j] * el_phi_last[j];
            gradphiold_p[1] += GNx_p[1][j] * el_phi_last[j];
            gradphiold_p[2] += GNx_p[2][j] * el_phi_last[j];
			
            gradphi_p[0] += GNx_p[0][j] * el_phi[j];
            gradphi_p[1] += GNx_p[1][j] * el_phi[j];
            gradphi_p[2] += GNx_p[2][j] * el_phi[j];
			
            v_p[0] += Ni_p[j] * el_V[NSD*j+0];      /* compute vx on the particle */
            v_p[1] += Ni_p[j] * el_V[NSD*j+1];      /* compute vy on the particle */
            v_p[2] += Ni_p[j] * el_V[NSD*j+2];      /* compute vz on the particle */
        }
        ConstructNiSUPG_Q1_3D(v_p,kappa_hat,Ni_p,GNx_p,Ni_supg_p);
        
		/*
         printf("qp=%d : v_p = %1.4e %1.4e %1.4e : kappa = %1.4e : kappa_hat = %1.4e : phi = %1.4e : phi_last = %1.4e \n",
         p,v_p[0],v_p[1],v_p[2],kappa_p,kappa_hat,phi_p,phi_last_p);
         */
		/*
         printf("Ni/Ni_supg 0(%1.4e %1.4e) 3(%1.4e %1.4e) 7(%1.4e %1.4e) \n",
         Ni_p[0],Ni_supg_p[0],Ni_p[3],Ni_supg_p[3],Ni_p[7],Ni_supg_p[7]);
         */
		
		M_dotT_p = 0.0;
		for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
			M_dotT_p += Ni_supg_p[j] * el_phi_last[j];
		}
		
		// F := ( M(x^k+1) + dt C(x^k+1) - dt K(x^k+1) ) T^k+1 - dt F(x^k+1) - M(x^k) T^k
        for (i = 0; i < NODES_PER_EL_Q1_3D; i++) {
            Re[ i ] += fac * (
                              // -dt.F
                              - dt * Ni_p[i] * f_p
                              // + dt.(C + L) T^k+1
                              + dt * kappa_p * ( GNx_p[0][i]*gradphi_p[0] + GNx_p[1][i]*gradphi_p[1] + GNx_p[2][i]*gradphi_p[2] )
                              + dt * Ni_supg_p[i] * ( v_p[0]*gradphi_p[0] + v_p[1]*gradphi_p[1] + v_p[2]*gradphi_p[2] )
                              // + M T^k+1
                              + Ni_supg_p[i] * phi_p
                              //+ Ni_p[i] * phi_p
                              // + Ni_p[i] T^k+1
                              )
            + fac * (
                     // - M(x^k) T^k
                     // - Ni_p[i] T^k
                     - Ni_supg_p[i] * phi_last_p
                     //- Ni_p[i] * phi_last_p
                     );
        }
    }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal_T"
PetscErrorCode FormFunctionLocal_T(
                                   PhysCompEnergy data,
                                   PetscReal dt,
                                   DM da,
                                   PetscScalar *LA_V,
                                   PetscScalar *LA_phi,
                                   PetscScalar *LA_philast,
                                   PetscScalar *LA_R)
{
    DM                     cda;
    Vec                    gcoords,gcoords_old;
    PetscScalar            *LA_gcoords, *LA_gcoords_old;
    PetscScalar            Re[NODES_PER_EL_Q1_3D];
    PetscScalar            el_coords[NSD*NODES_PER_EL_Q1_3D];
    PetscScalar            el_coords_old[NSD*NODES_PER_EL_Q1_3D];
	PetscScalar            el_V[NSD*NODES_PER_EL_Q1_3D];
	PetscScalar            el_phi[NODES_PER_EL_Q1_3D];
	PetscScalar            el_philast[NODES_PER_EL_Q1_3D];
	
	PetscInt       nel,nen,e,n;
	const PetscInt *elnidx;
	PetscInt       ge_eqnums[NODES_PER_EL_Q1_3D];
	PetscInt          nqp;
	PetscScalar       *qp_xi,*qp_weight;
	Quadrature        volQ;
	QPntVolCoefEnergy *all_quadpoints,*cell_quadpoints;
	PetscScalar       qp_kappa[27],qp_Q[27];
    PetscErrorCode   ierr;
	
    PetscFunctionBegin;
	//da     = data->daT;
	//V      = data->u_minus_V;
	//bclist = data->T_bclist;
	volQ   = data->volQ;
	
	/* quadrature */
	nqp       = volQ->npoints;
	qp_xi     = volQ->q_xi_coor;
	qp_weight = volQ->q_weight;
	
	ierr = VolumeQuadratureGetAllCellData_Energy(volQ,&all_quadpoints);CHKERRQ(ierr);
	
	/* setup for coords */
    ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
    ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
    
	/* setup for old coordinates */
    ierr = DMGetLocalVector(cda,&gcoords_old);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(cda,data->Xold,INSERT_VALUES,gcoords_old);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (cda,data->Xold,INSERT_VALUES,gcoords_old);CHKERRQ(ierr);
    ierr = VecGetArray(gcoords_old,&LA_gcoords_old);CHKERRQ(ierr);
	
	/* stuff for eqnums */
	ierr = DMDAGetElementsQ1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = DMDAEQ1_GetElementLocalIndicesDOF(ge_eqnums,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		/* get coords for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		ierr = DMDAEQ1_GetVectorElementField_3D(el_coords_old,(PetscInt*)&elnidx[nen*e],LA_gcoords_old);CHKERRQ(ierr);
		
		ierr = DMDAEQ1_GetScalarElementField_3D(el_phi,(PetscInt*)&elnidx[nen*e],LA_phi);CHKERRQ(ierr);
		ierr = DMDAEQ1_GetScalarElementField_3D(el_philast,(PetscInt*)&elnidx[nen*e],LA_philast);CHKERRQ(ierr);
		
		/* get velocity for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_V,(PetscInt*)&elnidx[nen*e],LA_V);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Energy(volQ,all_quadpoints,e,&cell_quadpoints);CHKERRQ(ierr);
		
		/* copy the diffusivity and force */
		for (n=0; n<nqp; n++) {
			qp_kappa[n] = cell_quadpoints[n].diffusivity;
			qp_Q[n]     = cell_quadpoints[n].heat_source;
		}
		
		/* initialise element stiffness matrix */
		ierr = PetscMemzero(Re,sizeof(PetscScalar)*NODES_PER_EL_Q1_3D);CHKERRQ(ierr);
		
		/* form element stiffness matrix */
#if (ADV_DIFF_STABILIZATION_TYPE == SUPG_STAB)
		ierr = AElement_FormFunction_T_supg(Re,dt,el_coords,el_coords_old,el_V,el_phi,el_philast,qp_kappa,qp_Q,nqp,qp_xi,qp_weight);CHKERRQ(ierr);
#elif (ADV_DIFF_STABILIZATION_TYPE == TAU_STAB)
		{
			PetscReal tau;
			PetscReal theta = 1.0; /* 0.5 crank nicoloson; 1.0 backward difference */
			PetscReal kappa_cell;
			
			kappa_cell = 0.0;
			for (n=0; n<nqp; n++) {
				kappa_cell = kappa_cell + qp_kappa[n];
			}
			kappa_cell = kappa_cell / ((PetscReal)nqp);
			
			//ierr = AdvDiffComputeTau_BrooksHughes(el_coords,el_V,kappa_cell,&tau);CHKERRQ(ierr);
			//printf("AdvDiffComputeTau_BrooksHughes : tau = %1.4e\n",tau);
			
			ierr = AdvDiffComputeTau_TezduyarOsawa(el_coords,el_V,kappa_cell,theta,dt,&tau);CHKERRQ(ierr);
			//printf("AdvDiffComputeTau_TezduyarOsawa : tau = %1.4e\n",tau);
			
			//ierr = AdvDiffComputeTau_UserDefinedConstant(NULL,&tau);CHKERRQ(ierr);
			//printf("AdvDiffComputeTau_UserDefinedConstant : tau = %1.4e\n",tau);
			
			ierr = AElement_FormFunction_T_tau(Re,dt,tau,el_coords,el_coords_old,el_V,el_phi,el_philast,qp_kappa,qp_Q,nqp,qp_xi,qp_weight);CHKERRQ(ierr);
		}
#endif
		
		ierr = DMDAEQ1_SetValuesLocalStencil_AddValues_DOF(LA_R,1,ge_eqnums,Re);CHKERRQ(ierr);
		
		/* diagnostics */
		/*
         for (p=0; p<nqp; p++) {
         P3D_ConstructNi_Q1_3D(&qp_xi[NSD*p],Ni_p);
         P3D_ConstructGNi_Q1_3D(&qp_xi[NSD*p],GNi_p);
         P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&J_p,&GNx_p[0],&GNx_p[1],&GNx_p[2]);
         
         phi_p = 0.0;
         for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
         phi_p = phi_p + Ni_p[i] * el_philast[i];
         }
         
         fac = qp_weight[p]*J_p;
         c = c + phi_p * fac;
         }
         */
	}
    //ierr = MPI_Allreduce(&c,&cg,1,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)da));CHKERRQ(ierr);
	//PetscPrintf(PETSC_COMM_WORLD,"\\int \\phi dV = %1.12e \n", cg );
	
    /* tidy up local arrays (input) */
    ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = VecRestoreArray(gcoords_old,&LA_gcoords_old);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(cda,&gcoords_old);CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TS_FormFunctionEnergy"
PetscErrorCode TS_FormFunctionEnergy(PetscReal time,Vec X,PetscReal dt,Vec F,void *ctx)
{
  pTatinCtx      ptatin = (pTatinCtx)ctx;
  PhysCompEnergy data;
    DM             da,cda;
	Vec            philoc, philastloc, Fphiloc;
	Vec            Vloc;
	PetscScalar    *LA_philoc, *LA_philastloc, *LA_Fphiloc;
	PetscScalar    *LA_V;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;

	ierr = pTatinGetContext_Energy(ptatin,&data);CHKERRQ(ierr);
	da = data->daT;
	
	
	ierr = VecZeroEntries(F);CHKERRQ(ierr);
	
    ierr = DMGetLocalVector(da,&philoc);CHKERRQ(ierr);
    ierr = DMGetLocalVector(da,&philastloc);CHKERRQ(ierr);
	
    ierr = DMGetLocalVector(da,&Fphiloc);CHKERRQ(ierr);
	
	/* get local solution and time derivative */
	ierr = VecZeroEntries(philoc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,philoc);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd  (da,X,INSERT_VALUES,philoc);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(philastloc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,data->Told,INSERT_VALUES,philastloc);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd  (da,data->Told,INSERT_VALUES,philastloc);CHKERRQ(ierr);
	
	/* insert boundary conditions into local vectors */
	ierr = BCListInsertLocal(data->T_bclist,philoc);CHKERRQ(ierr);
	
	/* init residual */
	ierr = VecZeroEntries(Fphiloc);CHKERRQ(ierr);
	
	/* get arrays */
	ierr = VecGetArray(philoc,    &LA_philoc);CHKERRQ(ierr);
	ierr = VecGetArray(philastloc,&LA_philastloc);CHKERRQ(ierr);
	ierr = VecGetArray(Fphiloc,   &LA_Fphiloc);CHKERRQ(ierr);
	
    /* get acces to the vector V */
    ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
    ierr = DMGetLocalVector(cda,&Vloc);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(cda,data->u_minus_V,INSERT_VALUES,Vloc);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(  cda,data->u_minus_V,INSERT_VALUES,Vloc);CHKERRQ(ierr);
    ierr = VecGetArray(Vloc,&LA_V);CHKERRQ(ierr);
	
	/* ===================== */
	/* Evaluate coefficients */
  ierr = EnergyEvaluateCoefficients(ptatin,time,da,LA_philoc,LA_V);CHKERRQ(ierr);
  
	/* ============= */
	/* FORM_FUNCTION */
	ierr = FormFunctionLocal_T(data,dt,da,LA_V, LA_philoc,LA_philastloc,LA_Fphiloc);CHKERRQ(ierr);
	
    ierr = VecRestoreArray(Vloc,&LA_V);CHKERRQ(ierr);
	
    ierr = DMRestoreLocalVector(cda,&Vloc);CHKERRQ(ierr);
	/* ============= */
	
	ierr = VecRestoreArray(Fphiloc,   &LA_Fphiloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(philastloc,&LA_philastloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(philoc,    &LA_philoc);CHKERRQ(ierr);
	
	
	/* do global fem summation */
	ierr = VecZeroEntries(F);CHKERRQ(ierr);
	
    ierr = DMLocalToGlobalBegin(da,Fphiloc,ADD_VALUES,F);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd  (da,Fphiloc,ADD_VALUES,F);CHKERRQ(ierr);
	
    ierr = DMRestoreLocalVector(da,&Fphiloc);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da,&philoc);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da,&philastloc);CHKERRQ(ierr);
	
	/* modify F for the boundary conditions, F_k = scale_k(x_k - phi_k) */
	ierr = BCListResidualDirichlet(data->T_bclist,X,F);CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNES_FormFunctionEnergy"
PetscErrorCode SNES_FormFunctionEnergy(SNES snes,Vec X,Vec F,void *ctx)
{
  pTatinCtx      ptatin = (pTatinCtx)ctx;
  PhysCompEnergy energy;
  PetscErrorCode ierr;

	PetscFunctionBegin;
	ierr = pTatinGetContext_Energy(ptatin,&energy);CHKERRQ(ierr);
	ierr = TS_FormFunctionEnergy(energy->time,X,energy->dt,F,ctx);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

/*
 Stabilized adv-diff written in terms of tau
 
 */

/* compute tau */
#undef __FUNCT__
#define __FUNCT__ "AdvDiffComputeTau_BrooksHughes"
PetscErrorCode AdvDiffComputeTau_BrooksHughes(PetscScalar el_coords[],PetscScalar el_vel[],PetscScalar kappa_el,PetscScalar *tau)
{
	PetscInt k,d;
	PetscScalar DX[NSD];
    PetscScalar u_xi[NSD];
    PetscScalar alpha_xi[NSD],_khat,dxi[NSD],xi[NSD];
	PetscErrorCode ierr;
	
    PetscFunctionBegin;
	
	/* average velocity - cell centre velocity */
	ierr = AdvDiff3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
	
	/* average velocity - cell centre velocity */
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
    for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] += 0.125 *  el_vel[NSD*k + 0];
		u_xi[1] += 0.125 *  el_vel[NSD*k + 1];
		u_xi[2] += 0.125 *  el_vel[NSD*k + 2];
	}
	
	/* could replace with / (1.0e-30 + kappa_el) */
	for (d=0; d<NSD; d++) {
		dxi[d]      = DX[d];
		if (kappa_el < ADV_DIFF_STAB_EPS) {
			alpha_xi[d] = 1.0/ADV_DIFF_STAB_EPS;
		} else {
			alpha_xi[d] = 0.5 * u_xi[d]  * dxi[d]  / kappa_el;
		}
		xi[d] = AdvDiffResidualForceTerm_UpwindXiExact(alpha_xi[d]);
	}
	
	//  khat = 0.5*(xi*u_xi*dxi + eta*v_eta*deta); /* steady state */
    _khat = 1.0 * 0.258198889747161 * ( xi[0]*u_xi[0]*dxi[0] + xi[1]*u_xi[1]*dxi[1] + xi[2]*u_xi[2]*dxi[2] ); /* transient case, sqrt(1/15) */
    *tau = _khat;
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AdvDiffComputeTau_TezduyarOsawa"
PetscErrorCode AdvDiffComputeTau_TezduyarOsawa(PetscScalar el_coords[],PetscScalar el_vel[],PetscScalar kappa_cell,PetscScalar theta,PetscScalar dt,PetscScalar *tau)
{
	PetscErrorCode ierr;
	PetscScalar h_cell,v_cell;
	PetscScalar tau1,tau2,tau3;
	PetscScalar t,DX[NSD],u_xi[NSD];
	PetscInt k;
	const PetscScalar r = 1.0; /* or 0.5 */
    
	PetscFunctionBegin;
	
	ierr = AdvDiff3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
	h_cell = sqrt( DX[0]*DX[0] + DX[1]*DX[1] + DX[2]*DX[2] );
    
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
    for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] += 0.125 *  el_vel[NSD*k + 0];
		u_xi[1] += 0.125 *  el_vel[NSD*k + 1];
		u_xi[2] += 0.125 *  el_vel[NSD*k + 2];
	}
	v_cell = sqrt( u_xi[0]*u_xi[0] + u_xi[1]*u_xi[1] + u_xi[2]*u_xi[2] );
	
	t = 0.0;
    
	if (v_cell > ADV_DIFF_STAB_EPS) {
		tau1 = 0.5 * h_cell / v_cell;
		t = t + 1.0/pow(tau1,r);
	}
	
	tau2 = theta * dt;
	t = t + 1.0/pow(tau2,r);
    
	if (kappa_cell > ADV_DIFF_STAB_EPS) {
		tau3 = (h_cell * h_cell ) / kappa_cell;
		t = t + 1.0/pow(tau3,r);
	}
	
	t = pow( t, -1.0/r );
	*tau = t;
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AdvDiffComputeTau_UserDefinedConstant"
PetscErrorCode AdvDiffComputeTau_UserDefinedConstant(PetscScalar const_t,PetscScalar *tau)
{
	static int been_here;
	static PetscReal user_tau=0.0;
	PetscBool flg;
	
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
    
	if (been_here==0) {
		flg = PETSC_FALSE;
		ierr = PetscOptionsGetReal(NULL,NULL,"-adv_diff_stab_tau",&user_tau,&flg);CHKERRQ(ierr);
		been_here = 1;
	}
	
	*tau = user_tau;
    PetscFunctionReturn(0);
}


/*
 
 Alternative implementation of SUPG in terms of stabilization parameter
 
 R(S) = dS/dt + div[ u , grad(S) ] - div( kappa grad(S) ) - f
 
 \int R(S) ( W + \tau div[ u, grad(S)) ] = 0
 
 M'(S,W) = \int S ( W + \tau div[ u, grad(S) ] )
 G(S,W)  = \int ( div[u,grad(S)] ) W
 K'(S,W) = \int grad(S) . ( I + \tau u \cross u ) . grad(W)
 F(S,W)    = \int ( W + \tau div[ u, grad(S) ] ) f
 
 R(S,W) := M'(S,W) ( dS/dt ) + (G(S,W) + K'(S,W) )S - F(S,W)
 
 R'(S) := M'(x^k+1) S^k+1 - M'(x^k) S^k + dt.G(x^k+1)S^k+1 + dt.K'(x^k+1)S^k+1 - dt.F(x^k+1)
 = [ M'(x^k+1) + dt.G(x^k+1) + dt.K'(x^k+1) ] S^k+1 - M'(x^k) S^k - dt.F(x^k+1)
 
 
 */


/*
 
 R'(S) := [ M'(x^k+1) + dt.G(x^k+1) + dt.K'(x^k+1) ] S^k+1 - M'(x^k) S^k - dt.F(x^k+1)
 J := M'(x^k+1) + dt.G(x^k+1) + dt.K'(x^k+1)
 
 */
#undef __FUNCT__
#define __FUNCT__ "AElement_FormJacobian_T_tau"
PetscErrorCode AElement_FormJacobian_T_tau(
                                           PetscScalar Re[],
                                           PetscReal dt,
                                           PetscScalar tau,
                                           PetscScalar el_coords[],
                                           PetscScalar gp_kappa[],
                                           PetscScalar el_V[],
                                           PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] )
{
    PetscInt    p,i,j,di,dj;
    PetscReal   Ni_p[NODES_PER_EL_Q1_3D];
    PetscReal   GNi_p[NSD][NODES_PER_EL_Q1_3D];
	PetscReal   GNx_p[NSD][NODES_PER_EL_Q1_3D];
	PetscReal   v_dot_gradW_i,v_dot_gradW_j, S_i;
    PetscScalar J_p,fac;
    PetscScalar kappa_p,v_p[NSD];
    PetscScalar kappa_3[NSD][NSD],u_cross_u[NSD][NSD],d_dx_kappa_d_dx,grad_i[NSD],grad_j[NSD];
	
	PetscFunctionBegin;
    
    /* evaluate integral */
    for (p=0; p<ngp; p++) {
        P3D_ConstructNi_Q1_3D(&gp_xi[NSD*p],Ni_p);
        P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
        P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&J_p,&GNx_p[0],&GNx_p[1],&GNx_p[2]);
		
        fac = gp_weight[p] * J_p;
		
		kappa_p = gp_kappa[p];
        v_p[0] = 0.0;
		v_p[1] = 0.0;
		v_p[2] = 0.0;
        for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
            v_p[0]     += Ni_p[j] * el_V[NSD*j+0];      /* compute vx on the particle */
            v_p[1]     += Ni_p[j] * el_V[NSD*j+1];      /* compute vy on the particle */
            v_p[2]     += Ni_p[j] * el_V[NSD*j+2];      /* compute vy on the particle */
        }
		//printf("%1.4e %1.4e %1.4e\n",v_p[0],v_p[1],v_p[2]);
		for (di=0; di<NSD; di++) {
			for (dj=0; dj<NSD; dj++) {
				u_cross_u[di][dj] = v_p[di] * v_p[dj];
			}
		}
		
		
		kappa_3[0][0] = kappa_p + tau * u_cross_u[0][0];
		kappa_3[0][1] =           tau * u_cross_u[0][1];
		kappa_3[0][2] =           tau * u_cross_u[0][2];
        
		kappa_3[1][0] =           tau * u_cross_u[1][0];
		kappa_3[1][1] = kappa_p + tau * u_cross_u[1][1];
		kappa_3[1][2] =           tau * u_cross_u[1][2];
        
		kappa_3[2][0] =           tau * u_cross_u[2][0];
		kappa_3[2][1] =           tau * u_cross_u[2][1];
		kappa_3[2][2] = kappa_p + tau * u_cross_u[2][2];
		
		
		for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
			grad_i[0] = GNx_p[0][i];
			grad_i[1] = GNx_p[1][i];
			grad_i[2] = GNx_p[2][i];
            
			v_dot_gradW_i = v_p[0]*grad_i[0] + v_p[1]*grad_i[1] + v_p[2]*grad_i[2];
			S_i = tau * v_dot_gradW_i;
            
			for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
				grad_j[0] = GNx_p[0][j];
				grad_j[1] = GNx_p[1][j];
				grad_j[2] = GNx_p[2][j];
				
				v_dot_gradW_j = v_p[0]*grad_j[0] + v_p[1]*grad_j[1] + v_p[2]*grad_j[2];
				
				d_dx_kappa_d_dx = 0.0;
				for (di=0; di<NSD; di++) {
					for (dj=0; dj<NSD; dj++) {
						d_dx_kappa_d_dx = d_dx_kappa_d_dx + grad_i[di] * kappa_3[di][dj] * grad_j[dj];
					}
				}
                
				Re[j+i*NODES_PER_EL_Q1_3D] += fac * (
                                                     ( Ni_p[i] + S_i ) * Ni_p[j]
                                                     + dt * Ni_p[i] * v_dot_gradW_j
                                                     + dt * d_dx_kappa_d_dx
                                                     );
                
                /*
                 Re[j+i*NODES_PER_EL_Q1_3D] += fac * (
                 Ni_p[i] * ( Ni_p[j] + S_j )
                 + dt * ( v_p[0]*GNx_p[0][i] + v_p[1]*GNx_p[1][i] + v_p[2]*GNx_p[2][i] ) * Ni_p[j]
                 + dt * d_dx_kappa_d_dx
                 );
                 */
                
			}
		}
    }
	
	/*
	 printf("e=\n");
	 for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
	 for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
	 printf("%lf ", Re[j+i*NODES_PER_EL_Q1_3D]);
	 }printf("\n");
	 }
	 */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AElement_FormFunction_T_tau"
PetscErrorCode AElement_FormFunction_T_tau(
                                           PetscScalar Re[],
                                           PetscReal dt,
                                           PetscReal tau,
                                           PetscScalar el_coords[],
                                           PetscScalar el_coords_old[],
                                           PetscScalar el_V[],
                                           PetscScalar el_phi[],PetscScalar el_phi_old[],
                                           PetscScalar gp_kappa[],PetscScalar gp_Q[],
                                           PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] )
{
    PetscInt    p,i,j,di,dj;
    PetscScalar Ni_p[NODES_PER_EL_Q1_3D],Ni_tau_p[NODES_PER_EL_Q1_3D];
    PetscScalar GNi_p[NSD][NODES_PER_EL_Q1_3D];
	PetscScalar GNx_p[NSD][NODES_PER_EL_Q1_3D], GNx_p_old[NSD][NODES_PER_EL_Q1_3D];
    PetscScalar phi_p,phi_p_old,f_p,v_p[NSD],kappa_p,gradphi_p[NSD],gradphiold_p[NSD],k_gradphi_p[NSD];
    PetscScalar J_p,J_p_old,fac;
    PetscScalar kappa_3[NSD][NSD],u_cross_u[NSD][NSD],grad_i[NSD],v_dot_gradphi;
	PetscScalar phi_p_old_tau,phi_p_tau;
	
	PetscFunctionBegin;
    /* evaluate integral */
    for (p=0; p<ngp; p++) {
        P3D_ConstructNi_Q1_3D(&gp_xi[NSD*p],Ni_p);
		P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
		
		P3D_evaluate_geometry_elementQ1(1,el_coords,    &GNi_p,&J_p,    &GNx_p[0],    &GNx_p[1],    &GNx_p[2]);
		P3D_evaluate_geometry_elementQ1(1,el_coords_old,&GNi_p,&J_p_old,&GNx_p_old[0],&GNx_p_old[1],&GNx_p_old[2]);
		
        fac     = gp_weight[p] * J_p;
        //fac_old = gp_weight[p] * J_p_old;
		
		kappa_p   = gp_kappa[p];
		f_p       = gp_Q[p];
        phi_p     = 0.0;
        phi_p_old = 0.0;
		for (di=0; di<NSD; di++) {
			v_p[di]          = 0.0;
			gradphi_p[di]    = 0.0;
			gradphiold_p[di] = 0.0;
		}
        for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
            phi_p     += Ni_p[j] * el_phi[j];      /* compute phi on the particle */
            phi_p_old += Ni_p[j] * el_phi_old[j];  /* compute phi_dot on the particle */
            
            gradphi_p[0] += GNx_p[0][j] * el_phi[j];
            gradphi_p[1] += GNx_p[1][j] * el_phi[j];
            gradphi_p[2] += GNx_p[2][j] * el_phi[j];
            
            gradphiold_p[0] += GNx_p_old[0][j] * el_phi_old[j];
            gradphiold_p[1] += GNx_p_old[1][j] * el_phi_old[j];
            gradphiold_p[2] += GNx_p_old[2][j] * el_phi_old[j];
			
            v_p[0] += Ni_p[j] * el_V[NSD*j+0];      /* compute vx on the particle */
            v_p[1] += Ni_p[j] * el_V[NSD*j+1];      /* compute vy on the particle */
            v_p[2] += Ni_p[j] * el_V[NSD*j+2];      /* compute vz on the particle */
        }
        
		for (di=0; di<NSD; di++) {
			for (dj=0; dj<NSD; dj++) {
				u_cross_u[di][dj] = v_p[di] * v_p[dj];
			}
		}
		
		kappa_3[0][0] = kappa_p + tau * u_cross_u[0][0];
		kappa_3[0][1] =           tau * u_cross_u[0][1];
		kappa_3[0][2] =           tau * u_cross_u[0][2];
		
		kappa_3[1][0] =           tau * u_cross_u[1][0];
		kappa_3[1][1] = kappa_p + tau * u_cross_u[1][1];
		kappa_3[1][2] =           tau * u_cross_u[1][2];
		
		kappa_3[2][0] =           tau * u_cross_u[2][0];
		kappa_3[2][1] =           tau * u_cross_u[2][1];
		kappa_3[2][2] = kappa_p + tau * u_cross_u[2][2];
		
        k_gradphi_p[0] = kappa_3[0][0]*gradphi_p[0] + kappa_3[0][1]*gradphi_p[1] + kappa_3[0][2]*gradphi_p[2];
        k_gradphi_p[1] = kappa_3[1][0]*gradphi_p[0] + kappa_3[1][1]*gradphi_p[1] + kappa_3[1][2]*gradphi_p[2];
        k_gradphi_p[2] = kappa_3[2][0]*gradphi_p[0] + kappa_3[2][1]*gradphi_p[1] + kappa_3[2][2]*gradphi_p[2];
		
		v_dot_gradphi = v_p[0]*gradphi_p[0] + v_p[1]*gradphi_p[1] + v_p[2]*gradphi_p[2];
		//S_j = tau * v_dot_gradphi;
        
		//v_dot_gradphi_old = v_p[0]*gradphiold_p[0] + v_p[1]*gradphiold_p[1] + v_p[2]*gradphiold_p[2];
		//S_j_old = tau * v_dot_gradphi_old;
		
        for (i = 0; i < NODES_PER_EL_Q1_3D; i++) {
			Ni_tau_p[i] = Ni_p[i] + tau * ( v_p[0]*GNx_p[0][i] + v_p[1]*GNx_p[1][i] + v_p[2]*GNx_p[2][i] );
		}
		phi_p_tau = phi_p_old_tau = 0.0;
        for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
            phi_p_tau     += Ni_tau_p[j] * el_phi[j];
            phi_p_old_tau += Ni_tau_p[j] * el_phi_old[j];
		}			
		
		// R'(S) = [ M'(x^k+1) + dt.G(x^k+1) + dt.K'(x^k+1) ] S^k+1 - M'(x^k) S^k - dt.F(x^k+1)
        for (i = 0; i < NODES_PER_EL_Q1_3D; i++) {
			
			grad_i[0] = GNx_p[0][i];
			grad_i[1] = GNx_p[1][i];
			grad_i[2] = GNx_p[2][i];
            
            
            Re[i] += fac * (
                            // - dt.F(x^k+1) 
                            - dt * Ni_tau_p[i] * f_p
                            // [ M'(x^k+1) + dt.G(x^k+1) + dt.K'(x^k+1) ] S^k+1
                            + Ni_tau_p[i] * phi_p
                            + dt * Ni_p[i] * v_dot_gradphi
                            + dt * ( grad_i[0]*k_gradphi_p[0] + grad_i[1]*k_gradphi_p[1] + grad_i[2]*k_gradphi_p[2] ) 
                            );
            Re[i] += fac * (
                            //- M'(x^k) S^k
                            - Ni_tau_p[i] * phi_p_old 
                            );
            
            /*			
             Re[i] += fac * (
             - dt * Ni_tau_p[i] * f_p
             + dt * ( v_p[0]*grad_i[0] + v_p[1]*grad_i[1] + v_p[2]*grad_i[2] ) * phi_p
             + dt * ( grad_i[0]*k_gradphi_p[0] + grad_i[1]*k_gradphi_p[1] + grad_i[2]*k_gradphi_p[2] ) 
             + Ni_p[i] * phi_p_tau
             - Ni_p[i] * phi_p_old_tau 
             );
             */			
        }
    }
	
	PetscFunctionReturn(0);	
}

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
 **    filename:   stokes_rheology_evss.c
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
 EVSS formulation

 Defines the deviatoric stress (tau)
 
 tau = 2 eta_v E[u] + tau_e
     = tau_v + tau_e

 Conservation of momentum then becomes
   div(2 eta_v E[u]) + div(tau_e) - grad(p) = f

 Conservation of mass is
   div(u) = 0
 
 The constitutive relationship for the extra stress (tau_e) is
 
   tau_e + lambda D(tau_e)/Dt + lambda B[u,tau_e] = 2 eta_e E[u]
 
 where 
   lambda = eta_e / mu (i.e. lambda is a relaxation time)
   B[] is an objective derivative (e.g. W tau - tau W)
   E[] is the strain-rate tensor
   D()/Dt is the material derivative
 
 Discretizing in time (semi-implicit)
 
 tau_e^{k+1} + lambda (tau_e^{k+1} - tau_e^{k})/dt 
             + lambda B[u^k,tau_e^k] = 2 eta_e E[u^{k+1}]
 
 tau_e^{k+1} = 2 eta_e E[u^{k+1}] (dt/(dt+lambda))
             - B[u^k,tau_e^k] (dt lambda/(dt + lambda))
             + lambda/(dt + lambda) tau_e^{k}
 
 Operating splitting with Lagranian particles storing tau_e
 
 [1] Initialize tau_e
 
 [2] Solve
   div(2 eta' E[u]) - grad(p) = f - div(  tau_e (lambda/(dt + lambda)) )
   eta' = eta_v + eta_e (dt/(dt+lambda))
 
 [3] Compute new stress at current position, e.g. local rate of change d(tau_e)/dt
   tau_e^{k+1} = 2 eta_e E[u^{k+1}] (dt/(dt+lambda))
               - B[u^k,tau_e^k] (dt lambda/(dt + lambda))
               + lambda/(dt + lambda) tau_e^{k}

 [4] Advect markers, thereby accounting for material derivative D(tau_e)/Dt
 
 Non-linear residual is this the weak form of
   div(2 eta' E[u]) - grad(p) - f + div(  tau_e (lambda/(dt + lambda)) )
 which is
   \int_v symgrad[v]: eta' symgrad[u]
   - \int_v div(v) p
   - \int_s v . (2 eta' E[u] - p I) . n
   + \int v . f
   + \int_v symgrad[v]: tau_e (lambda/(dt + lambda)
   - \int_s v . (tau_e) . n
 
 Note: The sign changes for \int v . f arise due to integrating by parts and multiplying by -1
 
*/

#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPEnergy_def.h"
#include "MPntPStokesPl_def.h"
#include "MPntPEVSS_def.h"

#include "QPntVolCoefStokes_def.h"
#include "QPntVolCoefEnergy_def.h"
#include "QPntVolCoefSymTens_def.h"

#include "material_constants.h"

#include "dmda_duplicate.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1.h"
#include "data_bucket.h"
#include "output_paraview.h"
#include "quadrature.h"
#include "element_type_Q2.h"
#include "material_point_utils.h"
#include "element_utils_q2.h"
#include "element_utils_q1.h"


#undef __FUNCT__
#define __FUNCT__ "EVSSComputeViscosities"
PetscErrorCode EVSSComputeViscosities(PetscReal alpha,PetscReal eta_ref,PetscReal *eta_v,PetscReal *eta_e)
{
  if (eta_v) *eta_v = alpha * eta_ref;
  if (eta_e) *eta_e = (1.0-alpha) * eta_ref;
  PetscFunctionReturn(0);
}

/*
 tau = 2 eta_v E[u] + tau_e
*/
#undef __FUNCT__
#define __FUNCT__ "EVSSComputeDeviatoricStress"
PetscErrorCode EVSSComputeDeviatoricStress(PetscReal eta_v,PetscReal symE[],MPntPEVSS *mp,PetscReal tau[])
{
  PetscInt k;
  for (k=0; k<6; k++) {
    tau[k] = 2.0 * eta_v * symE[k] + mp->tau[k];
  }
  PetscFunctionReturn(0);
}

/*
 tau_e^{k+1} = 2 eta_e E[u^{k+1}] (dt/(dt+lambda))
             - B[u^k,tau_e^k] (dt lambda/(dt + lambda))
             + lambda/(dt + lambda) tau_e^{k}
*/
#undef __FUNCT__
#define __FUNCT__ "EVSSUpdateExtraStress_Jaumman"
PetscErrorCode EVSSUpdateExtraStress_Jaumman(PetscReal dt,PetscReal L[],MPntPEVSS *mp)
{
  PetscInt k,ii,jj,kk,inc,nsubsteps;
  PetscReal E[6],W[3][3],T[3][3],B[3][3],tau[6],vT[6],eta_e,lambda_1,dt_s;

  
  eta_e = mp->eta_e;
  lambda_1 = eta_e/mp->mu;
  
  /* compute strain-rate */
  E[voigt_xx] = L[0];
  E[voigt_yy] = L[4];
  E[voigt_zz] = L[8];

  E[voigt_xy] = 0.5*(L[1] + L[3]);
  E[voigt_xz] = 0.5*(L[2] + L[6]);
  E[voigt_yz] = 0.5*(L[5] + L[7]);
 
  /* compute vorticity */
  for (ii=0; ii<3; ii++) {
    for (jj=0; jj<3; jj++) {
      W[ii][jj] = 0.5 * ( L[ii+jj*3] - L[jj+ii*3] );
    }
  }

  nsubsteps = 10;
  dt_s = dt/((PetscReal)(nsubsteps-1));

  /* Update spin over maxinc time sub-divisions */
  /* A better approach would be to use RK4 */
  for (inc=0; inc<nsubsteps; inc++) {

    /* copy voigt vector */
    tau[voigt_xx] = mp->tau[voigt_xx];
    tau[voigt_xy] = mp->tau[voigt_xy];
    tau[voigt_xz] = mp->tau[voigt_xz];
    
    tau[voigt_yy] = mp->tau[voigt_yy];
    tau[voigt_yz] = mp->tau[voigt_yz];
    
    tau[voigt_zz] = mp->tau[voigt_zz];

    /* convert voigt vector to symmetric tensor */
    T[0][0] = mp->tau[voigt_xx];
    T[0][1] = mp->tau[voigt_xy];
    T[0][2] = mp->tau[voigt_xz];

    T[1][0] = mp->tau[voigt_xy];
    T[1][1] = mp->tau[voigt_yy];
    T[1][2] = mp->tau[voigt_yz];

    T[2][0] = mp->tau[voigt_xz];
    T[2][1] = mp->tau[voigt_yz];
    T[2][2] = mp->tau[voigt_zz];

    /* B = WT - TW */
    for (ii=0; ii<3; ii++) {
      for (jj=0; jj<3; jj++) {
        B[ii][jj] = 0.0;
        for (kk=0; kk<3; kk++) {
          B[ii][jj] += W[ii][kk] * T[kk][jj];
          B[ii][jj] -= T[ii][kk] * W[kk][jj];
        }
      }
    }

    /* compute viscous part of extra stress */
    for (k=0; k<3; k++) {
      vT[k] = 2.0 * eta_e * E[k];
      vT[k+3] =     eta_e * E[k+3];
    }

    for (k=0; k<6; k++) {
      vT[k] = ( dt_s * vT[k] + lambda_1 * tau[k] ) / (dt_s + lambda_1);
    }
    
    vT[voigt_xx] -= ( dt_s * lambda_1/(dt_s + lambda_1) ) * B[0][0];
    vT[voigt_yy] -= ( dt_s * lambda_1/(dt_s + lambda_1) ) * B[1][1];
    vT[voigt_zz] -= ( dt_s * lambda_1/(dt_s + lambda_1) ) * B[2][2];

    vT[voigt_xy] -= ( dt_s * lambda_1/(dt_s + lambda_1) ) * B[0][1];
    vT[voigt_xz] -= ( dt_s * lambda_1/(dt_s + lambda_1) ) * B[0][2];
    vT[voigt_yz] -= ( dt_s * lambda_1/(dt_s + lambda_1) ) * B[1][2];

    for (k=0; k<6; k++) {
      mp->tau[k] = vT[k];
    }
  }
  
  PetscFunctionReturn(0);
}

/*
 Performs a P0 (constant) projection of the quantity
 
   phi tau_{ij}, phi = lambda/(lambda + dt), lambda = eta_e/mu
 
 from the material points onto each Q2 element
*/
#undef __FUNCT__
#define __FUNCT__ "MPntPEVSSProjection_P0"
PetscErrorCode MPntPEVSSProjection_P0(PetscReal dt,const int npoints,MPntStd mp_std[],MPntPEVSS mp_evss[],DM dmu,Quadrature Q)
{
	PetscErrorCode     ierr;
  DataField          PField;
  PetscInt           nqp,q,p,e,nel,k;
  QPntVolCoefSymTens *quadrature_T,*cell_quadrature_T;
	
	PetscFunctionBegin;
  
  /* Get quadrature point data */
  nqp = Q->npoints;
  
  if (nqp < 2)  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"P0 projection requires at least two quadrature points per cell <lame/lazy design issue>");
  
  DataBucketGetDataFieldByName(Q->properties_db,QPntVolCoefSymTens_classname,&PField);
  DataFieldGetEntries(PField,(void**)&quadrature_T);
  
  ierr = DMDAGetElements_pTatinQ2P1(dmu,&nel,NULL,NULL);CHKERRQ(ierr);
  
	/* traverse elements and initialize zeroth component of stress on each quadrature point */
	for (e=0; e<nel; e++) {
    cell_quadrature_T = &quadrature_T[e*nqp];
    
    cell_quadrature_T[0].T[voigt_xx] = 0.0;
    cell_quadrature_T[0].T[voigt_yy] = 0.0;
    cell_quadrature_T[0].T[voigt_zz] = 0.0;
    
    cell_quadrature_T[0].T[voigt_xy] = 0.0;
    cell_quadrature_T[0].T[voigt_xz] = 0.0;
    cell_quadrature_T[0].T[voigt_yz] = 0.0;
    
    cell_quadrature_T[1].T[voigt_xx] = 0.0;
	}
  
  /* traverse points and sum */
  for (p=0; p<npoints; p++) {
    PetscReal lambda_1,factor;
    
    e = mp_std[p].wil;
    cell_quadrature_T = &quadrature_T[e*nqp];
    
    lambda_1 = mp_evss[p].eta_e / mp_evss[p].mu;
    factor = lambda_1/(lambda_1 + dt);
    
    for (k=0; k<6; k++) {
      cell_quadrature_T[0].T[k] += factor * mp_evss[p].tau[k];
    }
    cell_quadrature_T[1].T[voigt_xx] += 1.0;
  }
  
  /* assign constant values to each cell */
  for (e=0; e<nel; e++) {
    PetscReal T0[6],np_per_el;
    
    cell_quadrature_T = &quadrature_T[e*nqp];
    
    np_per_el = cell_quadrature_T[1].T[voigt_xx];
    for (k=0; k<6; k++) {
      T0[k] = cell_quadrature_T[0].T[k] / np_per_el;
    }
    
    for (q=0; q<nqp; q++) {
      cell_quadrature_T[q].T[voigt_xx] = T0[voigt_xx];
      cell_quadrature_T[q].T[voigt_yy] = T0[voigt_yy];
      cell_quadrature_T[q].T[voigt_zz] = T0[voigt_zz];
      
      cell_quadrature_T[q].T[voigt_xy] = T0[voigt_xy];
      cell_quadrature_T[q].T[voigt_xz] = T0[voigt_xz];
      cell_quadrature_T[q].T[voigt_yz] = T0[voigt_yz];
    }
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateRheologyNonlinearitiesMarkers_ViscousEVSS"
PetscErrorCode EvaluateRheologyNonlinearitiesMarkers_ViscousEVSS(pTatinCtx user,DM dau,PetscScalar u[],DM dap,PetscScalar p[])
{
  PetscErrorCode ierr;
  
  int            pidx,n_mp_points;
  DataBucket     db;
  DataField      PField_std,PField_stokes,PField_evss;
  double         eta_ref,eta_v,eta_e,mu,lambda,factor,eta0;
  PetscReal      dt,alpha;
  int            region_idx;
	DataBucket                   material_constants;
	MaterialConst_ViscosityConst *ViscConst_data;
	DataField                    PField_ViscConst;
  
  PetscFunctionBegin;
  
  ierr = pTatinGetTimestep(user,&dt);CHKERRQ(ierr);
  
  ierr = pTatinGetMaterialPoints(user,&db,NULL);CHKERRQ(ierr);
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
  
  DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
  DataFieldGetAccess(PField_stokes);
  DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));
  
  DataBucketGetDataFieldByName(db,MPntPEVSS_classname,&PField_evss);
  DataFieldGetAccess(PField_evss);
  DataFieldVerifyAccess(PField_evss,sizeof(MPntPEVSS));
  
  DataBucketGetSizes(db,&n_mp_points,0,0);
  
	ierr = pTatinGetMaterialConstants(user,&material_constants);CHKERRQ(ierr);

	DataBucketGetDataFieldByName(material_constants,MaterialConst_ViscosityConst_classname,&PField_ViscConst);
  DataFieldGetEntries(PField_ViscConst,(void**)&ViscConst_data);
  
  alpha = 0.5;
  PetscOptionsGetReal(NULL,"-evss_alpha",&alpha,NULL);
  if (alpha < 1.0e-10) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-evss_alpha > 0");
  if (alpha > 1.0) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-evss_alpha <= 1 ");
  
  for (pidx=0; pidx<n_mp_points; pidx++) {
    MPntStd     *mp;
    MPntPStokes *mpprop_stokes;
    MPntPEVSS   *mpprop_evss;
    
    DataFieldAccessPoint(PField_std,   pidx,(void**)&mp);
    DataFieldAccessPoint(PField_stokes,pidx,(void**)&mpprop_stokes);
    DataFieldAccessPoint(PField_evss,  pidx,(void**)&mpprop_evss);
    
    MPntStdGetField_phase_index(mp,&region_idx);
    
    eta_ref = ViscConst_data[region_idx].eta0;

    MPntPEVSSGetField_shear_modulus(mpprop_evss,&mu);

    ierr = EVSSComputeViscosities(alpha,eta_ref,&eta_v,&eta_e);CHKERRQ(ierr);
    
    lambda = eta_e/mu;
    factor = dt / (dt + lambda);
    
    eta0 = eta_v + factor * eta_e;
    
    MPntPEVSSSetField_solid_viscosity(mpprop_evss,eta_e);
    MPntPStokesSetField_eta_effective(mpprop_stokes,eta0);
  }
  DataFieldRestoreAccess(PField_evss);
  DataFieldRestoreAccess(PField_stokes);
  DataFieldRestoreAccess(PField_std);
  
  PetscFunctionReturn(0);
}

static inline void ComputeGradU3d(PetscReal ux[],PetscReal uy[],PetscReal uz[],
                                  PetscReal dNudx[],PetscReal dNudy[],PetscReal dNudz[],PetscReal L[])
{
	PetscInt k;

	for (k=0; k<9; k++) {
    L[k] = 0.0;
  }
	
	for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
		L[0] += dNudx[k] * ux[k]; // du/dx_j
		L[1] += dNudy[k] * ux[k];
		L[2] += dNudz[k] * ux[k];

		L[3] += dNudx[k] * uy[k]; // dv/dx_j
		L[4] += dNudy[k] * uy[k];
		L[5] += dNudz[k] * uy[k];

    L[6] += dNudx[k] * uz[k]; // dw/dx_j
		L[7] += dNudy[k] * uz[k];
		L[8] += dNudz[k] * uz[k];
  }
}

#undef __FUNCT__
#define __FUNCT__ "StokesCoefficient_UpdateTimeDependentQuantities_EVSS"
PetscErrorCode StokesCoefficient_UpdateTimeDependentQuantities_EVSS(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[])
{
  PetscErrorCode ierr;
  DM             cda;
  Vec            gcoords;
  PetscScalar    *LA_gcoords;
  int            pidx,n_mp_points;
  DataBucket     db;
  DataField      PField_std,PField;
  PetscInt       nel,nen_u,k;
  const PetscInt *elnidx_u;
  PetscReal      elcoords[3*Q2_NODES_PER_EL_3D];
  PetscReal      elu[3*Q2_NODES_PER_EL_3D];
  PetscReal      ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
  PetscReal      GNI[3][Q2_NODES_PER_EL_3D];
  PetscReal      dNudx[Q2_NODES_PER_EL_3D],dNudy[Q2_NODES_PER_EL_3D],dNudz[Q2_NODES_PER_EL_3D];
  int            eidx_mp;
  double         *xi_mp;
  PetscReal      L_mp[9];
  PetscReal      dt;
  PetscFunctionBegin;
  
  /* access current time step */
  ierr = pTatinGetTimestep(user,&dt);CHKERRQ(ierr);
  
  /* access material point information */
  ierr = pTatinGetMaterialPoints(user,&db,NULL);CHKERRQ(ierr);
  
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  
  DataBucketGetDataFieldByName(db,MPntPEVSS_classname,&PField);
  DataFieldGetAccess(PField);
  
  DataBucketGetSizes(db,&n_mp_points,0,0);
  
  /* setup for coords */
  ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  
  /* setup for elements */
  ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
  
  /* marker loop */
  for (pidx=0; pidx<n_mp_points; pidx++) {
    MPntStd   *mpprop_std;
    MPntPEVSS *mpprop;
    
    DataFieldAccessPoint(PField_std, pidx,(void**)&mpprop_std);
    DataFieldAccessPoint(PField,     pidx,(void**)&mpprop);
    
    eidx_mp = mpprop_std->wil;
    
    /* Get element coordinates */
    ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*eidx_mp],LA_gcoords);CHKERRQ(ierr);
    
    /* Get element velocity */
    ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*eidx_mp],ufield);CHKERRQ(ierr);

    /* get velocity components */
    for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
      ux[k] = elu[3*k  ];
      uy[k] = elu[3*k+1];
      uz[k] = elu[3*k+2];
    }
    
    /* Get local coordinate of marker */
    xi_mp = mpprop_std->xi;
    
    /* Prepare basis functions derivatives */
    P3D_ConstructGNi_Q2_3D(xi_mp,GNI);
    P3D_evaluate_global_derivatives_Q2(elcoords,GNI,dNudx,dNudy,dNudz);
    
    /* strain rate */
    ComputeGradU3d(ux,uy,uz,dNudx,dNudy,dNudz,L_mp);
    
    ierr = EVSSUpdateExtraStress_Jaumman(dt,L_mp,mpprop);CHKERRQ(ierr);
  }
  DataFieldRestoreAccess(PField);
  DataFieldRestoreAccess(PField_std);
  ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

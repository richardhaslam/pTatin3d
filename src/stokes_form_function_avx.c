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
 **    filename:   stokes_form_function.c
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

#if defined(__AVX__)

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>
#include <petscsnes.h>

#include <ptatin3d.h>
#include <ptatin3d_defs.h>
#include <private/ptatin_impl.h>
#include <data_bucket.h>
#include <element_type_Q2.h>
#include <dmda_element_q2p1.h>
#include <element_utils_q2.h>
#include <element_utils_q1.h>
#include <QPntVolCoefStokes_def.h>
#include <quadrature.h>
#include <stokes_operators.h>
#include <immintrin.h>


#ifndef __FMA__
#  define _mm256_fmadd_pd(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b),c)
#endif

#define ALIGN32 __attribute__((aligned(32))) /* AVX packed instructions need 32-byte alignment */

/*
 * Performs three tensor contractions: y[l,a,b,c] += T[a,k] S[b,j] R[c,i] x[l,k,j,i]
 */
PetscErrorCode formf_TensorContractNEV_AVX(PetscReal Rf[][3],PetscReal Sf[][3],PetscReal Tf[][3],GradMode gmode,PetscReal x[][NQP][NEV],PetscReal y[][NQP][NEV])
{
  PetscReal R[3][3],S[3][3],T[3][3];
  PetscReal u[3][NQP][NEV] ALIGN32,v[3][NQP][NEV] ALIGN32;
  PetscInt i,j,k,l,kj,ji,a,b,c;
  
  for (j=0; j<3; j++) {
    for (i=0; i<3; i++) {
      R[i][j] = i<3 ? (gmode == GRAD ? Rf[i][j] : Rf[j][i]) : 0.;
      S[i][j] = i<3 ? (gmode == GRAD ? Sf[i][j] : Sf[j][i]) : 0.;
      T[i][j] = i<3 ? (gmode == GRAD ? Tf[i][j] : Tf[j][i]) : 0.;
    }
  }
  
  // u[l,k,j,c] = R[c,i] x[l,k,j,i]
  for (l=0; l<3; l++) {
    for (c=0; c<3; c++) {
      __m256d r[3] = {_mm256_set1_pd(R[c][0]),_mm256_set1_pd(R[c][1]),_mm256_set1_pd(R[c][2])};
      for (kj=0; kj<9; kj++) {
        __m256d u_lkjc = _mm256_setzero_pd();
        for (i=0; i<3; i++) {
          __m256d x_lkji = _mm256_load_pd(x[l][kj*3+i]);
          u_lkjc = _mm256_fmadd_pd(r[i],x_lkji,u_lkjc);
        }
        _mm256_store_pd(u[l][kj*3+c],u_lkjc);
      }
    }
  }
  
  // v[l,k,b,c] = S[b,j] u[l,k,j,c]
  for (l=0; l<3; l++) {
    for (k=0; k<3; k++) {
      for (b=0; b<3; b++) {
        __m256d s[3] = {_mm256_set1_pd(S[b][0]),_mm256_set1_pd(S[b][1]),_mm256_set1_pd(S[b][2])};
        for (c=0; c<3; c++) {
          __m256d v_lkbc = _mm256_setzero_pd();
          for (j=0; j<3; j++) {
            __m256d u_lkjc = _mm256_load_pd(u[l][(k*3+j)*3+c]);
            v_lkbc = _mm256_fmadd_pd(s[j],u_lkjc,v_lkbc);
          }
          _mm256_store_pd(v[l][(k*3+b)*3+c],v_lkbc);
        }
      }
    }
  }
  
  // y[l,a,b,c] = T[a,k] v[l,k,b,c]
  for (a=0; a<3; a++) {
    __m256d t[3] = {_mm256_set1_pd(T[a][0]),_mm256_set1_pd(T[a][1]),_mm256_set1_pd(T[a][2])};
    for (ji=0; ji<9; ji++) {
      for (l=0; l<3; l++) {
        __m256d y_laji = _mm256_load_pd(y[l][a*9+ji]);
        for (k=0; k<3; k++) {
          __m256d v_lkji = _mm256_load_pd(v[l][k*9+ji]);
          y_laji = _mm256_fmadd_pd(v_lkji,t[k],y_laji);
        }
        _mm256_store_pd(y[l][a*9+ji],y_laji);
      }
    }
  }
  return 0;
}

PetscErrorCode formf_JacobianInvertNEV_AVX(PetscScalar dx[3][3][NQP][NEV],PetscScalar dxdet[NQP][NEV])
{
  PetscInt i,j,k,e;
  
  for (i=0; i<NQP; i++) {
    PetscScalar a[3][3][NEV] ALIGN32;
    for (e=0; e<NEV; e++) {
      PetscScalar b0,b3,b6,det,idet;
      for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
          a[j][k][e] = dx[j][k][i][e];
        }
      }
      b0 =  (a[1][1][e]*a[2][2][e] - a[2][1][e]*a[1][2][e]);
      b3 = -(a[1][0][e]*a[2][2][e] - a[2][0][e]*a[1][2][e]);
      b6 =  (a[1][0][e]*a[2][1][e] - a[2][0][e]*a[1][1][e]);
      det = a[0][0][e]*b0 + a[0][1][e]*b3 + a[0][2][e]*b6;
      idet = 1.0 / det;
      dx[0][0][i][e] =  idet*b0;
      dx[0][1][i][e] = -idet*(a[0][1][e]*a[2][2][e] - a[2][1][e]*a[0][2][e]);
      dx[0][2][i][e] =  idet*(a[0][1][e]*a[1][2][e] - a[1][1][e]*a[0][2][e]);
      dx[1][0][i][e] =  idet*b3;
      dx[1][1][i][e] =  idet*(a[0][0][e]*a[2][2][e] - a[2][0][e]*a[0][2][e]);
      dx[1][2][i][e] = -idet*(a[0][0][e]*a[1][2][e] - a[1][0][e]*a[0][2][e]);
      dx[2][0][i][e] =  idet*b6;
      dx[2][1][i][e] = -idet*(a[0][0][e]*a[2][1][e] - a[2][0][e]*a[0][1][e]);
      dx[2][2][i][e] =  idet*(a[0][0][e]*a[1][1][e] - a[1][0][e]*a[0][1][e]);
      dxdet[i][e] =  det;
    }
  }
  return 0;
}

__attribute__((noinline))
static PetscErrorCode formf_QuadratureAction_A_AVX(const QPntVolCoefStokes *gausspt[],
                                             PetscScalar dx[3][3][NQP][NEV],
                                             PetscScalar dxdet[NQP][NEV],
                                             PetscReal w[NQP],
                                             PetscScalar du[3][3][NQP][NEV],
                                             PetscScalar dp[NQP][NEV],
                                             PetscScalar dv[3][3][NQP][NEV])
{
  PetscInt i,l,k,e;
  
  for (i=0; i<NQP; i++) {
    PetscScalar Du[6][NEV] ALIGN32,Dv[6][NEV] ALIGN32; /* Symmetric gradient with respect to physical coordinates, xx, yy, zz, xy+yx, xz+zx, yz+zy */
    __m256d dux[3][3],mhalf = _mm256_set1_pd(0.5),dvx[3][3];
    __m256d mweight = _mm256_mul_pd(_mm256_set1_pd(w[i]),_mm256_load_pd(dxdet[i]));
    
    for (k=0; k<3; k++) { // directions
      __m256d dxk[3] = {_mm256_load_pd(dx[k][0][i]),_mm256_load_pd(dx[k][1][i]),_mm256_load_pd(dx[k][2][i])};
      for (l=0; l<3; l++) { // fields
        dux[k][l] = _mm256_mul_pd(_mm256_load_pd(du[0][l][i]),dxk[0]);
        dux[k][l] = _mm256_fmadd_pd(_mm256_load_pd(du[1][l][i]),dxk[1],dux[k][l]);
        dux[k][l] = _mm256_fmadd_pd(_mm256_load_pd(du[2][l][i]),dxk[2],dux[k][l]);
      }
    }
    _mm256_store_pd(Du[0],dux[0][0]);
    _mm256_store_pd(Du[1],dux[1][1]);
    _mm256_store_pd(Du[2],dux[2][2]);
    _mm256_store_pd(Du[3],_mm256_mul_pd(mhalf,_mm256_add_pd(dux[0][1],dux[1][0])));
    _mm256_store_pd(Du[4],_mm256_mul_pd(mhalf,_mm256_add_pd(dux[0][2],dux[2][0])));
    _mm256_store_pd(Du[5],_mm256_mul_pd(mhalf,_mm256_add_pd(dux[1][2],dux[2][1])));
    
    for (e=0; e<NEV; e++) {
      for (k=0; k<6; k++) { /* Stress is coefficient of test function */
        Dv[k][e] = 2 * gausspt[e][i].eta * Du[k][e];
      }
      for (k=0; k<3; k++) { /* shift diagonal components of stress by -p */
        Dv[k][e] -= dp[i][e];
      }
    }
    
    dvx[0][0] = _mm256_load_pd(Dv[0]);
    dvx[0][1] = _mm256_load_pd(Dv[3]);
    dvx[0][2] = _mm256_load_pd(Dv[4]);
    dvx[1][0] = _mm256_load_pd(Dv[3]);
    dvx[1][1] = _mm256_load_pd(Dv[1]);
    dvx[1][2] = _mm256_load_pd(Dv[5]);
    dvx[2][0] = _mm256_load_pd(Dv[4]);
    dvx[2][1] = _mm256_load_pd(Dv[5]);
    dvx[2][2] = _mm256_load_pd(Dv[2]);
    
    for (l=0; l<3; l++) { // fields
      for (k=0; k<3; k++) { // directions
        __m256d sum = _mm256_mul_pd(dvx[0][l],_mm256_load_pd(dx[0][k][i]));
        sum = _mm256_fmadd_pd(dvx[1][l],_mm256_load_pd(dx[1][k][i]),sum);
        sum = _mm256_fmadd_pd(dvx[2][l],_mm256_load_pd(dx[2][k][i]),sum);
        _mm256_store_pd(dv[k][l][i],_mm256_mul_pd(mweight,sum));
      }
    }
  }
  //PetscLogFlops(NQP*NEV*(5*9+6+6+3+6*9));
  return 0;
}


PetscErrorCode FormFunctionLocal_UP(Quadrature volQ,
                                     DM dau,PetscScalar ufield[],
                                     DM dap,PetscScalar pfield[],
                                     PetscScalar Yu[],
                                     PetscScalar Yp[])
{
  PetscErrorCode ierr;
  DM cda;
  Vec gcoords;
  const PetscReal *LA_gcoords;
  PetscInt e,i,j,k,p;
  const PetscInt *elnidx_u;
  const PetscInt *elnidx_p;
  PetscInt nel,nen_u,nen_p;
  QPntVolCoefStokes *all_gausspoints;
  const QPntVolCoefStokes *cell_gausspoints[NEV];
  PetscReal x1[3],w1[3],B[3][3],D[3][3],w[NQP],xi[NQP][3];
  PetscReal NI_geom[NQP][Q1_NODES_PER_EL_3D],NIp[NEV][NQP][4],NIu[NQP][Q2_NODES_PER_EL_3D];
  
  PetscFunctionBegin;
  ierr = PetscDTGaussQuadrature(3,-1,1,x1,w1);CHKERRQ(ierr);
  for (i=0; i<3; i++) {
    B[i][0] = .5*(PetscSqr(x1[i]) - x1[i]);
    B[i][1] = 1 - PetscSqr(x1[i]);
    B[i][2] = .5*(PetscSqr(x1[i]) + x1[i]);
    D[i][0] = x1[i] - .5;
    D[i][1] = -2*x1[i];
    D[i][2] = x1[i] + .5;
  }
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        w[(i*3+j)*3+k] = w1[i] * w1[j] * w1[k];
      }}}
  
  p = 0;
  for (k=0; k<3; k++) {
    for (j=0; j<3; j++) {
      for (i=0; i<3; i++) {
        PetscReal XI[3];
        
        XI[0] = x1[i];
        XI[1] = x1[j];
        XI[2] = x1[k];
        xi[p][0] = XI[0];
        xi[p][1] = XI[1];
        xi[p][2] = XI[2];
        P3D_ConstructNi_Q1_3D(XI,NI_geom[p]);
        //_w[p] = w1[i] * w1[j] * w1[k];
        p++;
      }}}
  
  for (p=0; p<NQP; p++) {
    P3D_ConstructNi_Q2_3D(xi[p],NIu[p]);
  }
  
  /* setup for coords */
  ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);
  
  ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
  ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
  
  ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);
  
  for (e=0;e<nel;e+=NEV) {
    PetscInt     ee,l;
    PetscScalar  elu[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32;
    PetscScalar  elx[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32;
    PetscScalar  elv[3][Q2_NODES_PER_EL_3D][NEV] ALIGN32;
    PetscScalar  dx[3][3][NQP][NEV] ALIGN32;
    PetscScalar  dxdet[NQP][NEV];
    PetscScalar  du[3][3][NQP][NEV] ALIGN32;
    PetscScalar  dv[3][3][NQP][NEV] ALIGN32;
    PetscReal    _elx[NEV][3*Q2_NODES_PER_EL_3D] ALIGN32;
    PetscReal    elp[4][NEV] ALIGN32;
    PetscReal    elq[4][NEV] ALIGN32;
    PetscScalar  dp[NQP][NEV] ALIGN32;
    PetscReal    Be_u[3][NQP][NEV] ALIGN32;
    PetscReal    Be_p[NQP][NEV] ALIGN32;

    
    for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
      for (ee=0; ee<NEV; ee++) {
        PetscInt E = elnidx_u[nen_u*PetscMin(e+ee,nel-1)+i]; // Pad up to length NEV by duplicating last element
        for (l=0; l<3; l++) {
          elx[l][i][ee] = LA_gcoords[3*E+l];
          _elx[ee][3*i+l] = LA_gcoords[3*E+l];
          
          elu[l][i][ee] = ufield[3*E+l];
        }
      }
    }
    
    /* pressure */
    for (i=0; i<4; i++) {
      for (ee=0; ee<NEV; ee++) {
        PetscInt E = elnidx_p[nen_p*PetscMin(e+ee,nel-1)+i]; // Pad up to length NEV by duplicating last element
        
        elp[i][ee] = pfield[E];
      }
    }
    
    for (p=0; p<NQP; p++) {
      for (ee=0; ee<NEV; ee++) {
        ConstructNi_pressure(xi[p],_elx[ee],NIp[ee][p]);
      }
    }
    
    for (ee=0; ee<NEV; ee++) {
      ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,PetscMin(e+ee,nel-1),(QPntVolCoefStokes**)&cell_gausspoints[ee]);CHKERRQ(ierr);
    }
    
    ierr = PetscMemzero(dx,sizeof dx);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(D,B,B,GRAD,elx,dx[0]);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(B,D,B,GRAD,elx,dx[1]);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(B,B,D,GRAD,elx,dx[2]);CHKERRQ(ierr);
    
    ierr = formf_JacobianInvertNEV_AVX(dx,dxdet);CHKERRQ(ierr);
    
    // copy source terms
    for (p=0; p<NQP; p++) {
      for (ee=0; ee<NEV; ee++) {
        for (i=0; i<3; i++) {
          Be_u[i][p][ee] = cell_gausspoints[ee][p].Fu[i];
        }
      }
    }
    
    for (p=0; p<NQP; p++) {
      for (ee=0; ee<NEV; ee++) {
        Be_p[p][ee] = cell_gausspoints[ee][p].Fp;
      }
    }

    // u - momentum
    ierr = PetscMemzero(du,sizeof du);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(D,B,B,GRAD,elu,du[0]);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(B,D,B,GRAD,elu,du[1]);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(B,B,D,GRAD,elu,du[2]);CHKERRQ(ierr);
    
    // evaluate p at each quadrature point in each element
    ierr = PetscMemzero(dp,sizeof dp);CHKERRQ(ierr);
    for (p=0; p<NQP; p++) {
      for (i=0; i<4; i++) {
        for (ee=0; ee<NEV; ee++) {
          dp[p][ee] += NIp[ee][p][i] * elp[i][ee];
        }
      }
    }

    ierr = formf_QuadratureAction_A_AVX(cell_gausspoints,dx,dxdet,w,du,dp,dv);CHKERRQ(ierr);
    
    ierr = PetscMemzero(elv,sizeof elv);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(D,B,B,GRAD_TRANSPOSE,dv[0],elv);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(B,D,B,GRAD_TRANSPOSE,dv[1],elv);CHKERRQ(ierr);
    ierr = formf_TensorContractNEV_AVX(B,B,D,GRAD_TRANSPOSE,dv[2],elv);CHKERRQ(ierr);

    // this could also be a tensor contraction
    // elv[3][Q2_NODES_PER_EL_3D][NEV] Be_u[3][NQP][NEV]
    for (p=0; p<NQP; p++) {
      for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
        for (ee=0; ee<NEV; ee++) {
          for (l=0; l<3; l++) {
            elv[l][i][ee] -= w[p] * ( NIu[p][i] * Be_u[l][p][ee] ) * dxdet[p][ee];
          }
        }
      }
    }
    
    // p - continuity
    // evaluate -div(u) at each quadrature point in each element
    ierr = PetscMemzero(dp,sizeof dp);CHKERRQ(ierr);
    for (p=0; p<NQP; p++) {
      for (ee=0; ee<NEV; ee++) {
        PetscReal gradx[3];
        
        //dp[p][ee] = -(du[0][0][p][ee] + du[1][1][p][ee] + du[2][2][p][ee]);
        
        for (k=0; k<3; k++) { // directions
          gradx[k] = du[0][k][p][ee] * dx[k][0][p][ee]
                   + du[1][k][p][ee] * dx[k][1][p][ee]
                   + du[2][k][p][ee] * dx[k][2][p][ee];
        }
        dp[p][ee] = -(gradx[0] + gradx[1] + gradx[2]);
      }
    }
    
    ierr = PetscMemzero(elq,sizeof elq);CHKERRQ(ierr);
    for (p=0; p<NQP; p++) {
      for (i=0; i<4; i++) {
        for (ee=0; ee<NEV; ee++) {//NIp[NEV][NQP][4]
          elq[i][ee] += w[p] * ( NIp[ee][p][i] * (dp[p][ee] - Be_p[p][ee]) ) * dxdet[p][ee];
        }
      }
    }
    
    /* add contibutions */
    for (ee=0; ee<PetscMin(NEV,nel-e); ee++) {
      
      for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
        PetscInt E = elnidx_u[nen_u*(e+ee)+i];
        
        for (l=0; l<3; l++) {
          Yu[3*E+l] += elv[l][i][ee];
        }
      }
      
      for (i=0; i<4; i++) {
        PetscInt E = elnidx_p[nen_p*(e+ee)+i];
        
        Yp[E] += elq[i][ee];
      }
    }
    
  }
  
  ierr = VecRestoreArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);
  
  {
    PetscLogDouble flops = 0.0;
    
    flops += 3 + NQP*(Q1_NODES_PER_EL_3D*(2*3) + 2*3); /* pressure basis evaluation */
    flops += 9 * 3*NQP*(6+6+6);           /* 9 tensor contractions per element */
    flops += NQP*(14 + 1/* division */ + 27); /* 3x3 matrix inversion + determinant per element */
    flops += NQP*(5*9+6+6+3+6*9);             /* quadrature action per element */
    // haven't logged flops associated with -Be_u and -Be_p
    flops = flops * ((PetscLogDouble)nel);
    PetscLogFlops(flops);
  }
  
  PetscFunctionReturn(0);
}


PetscErrorCode FormFunction_Stokes_AVX(SNES snes,Vec X,Vec F,void *ctx)
{
  PetscErrorCode    ierr;
  pTatinCtx         ptatin;
  DM                stokes_pack,dau,dap;
  DMDALocalInfo     infou,infop;
  Vec               Uloc,Ploc,FUloc,FPloc;
  Vec               u,p,Fu,Fp;
  PetscScalar       *LA_Uloc,*LA_Ploc;
  PetscScalar       *LA_FUloc,*LA_FPloc;
  PhysCompStokes    stokes;
  
  PetscFunctionBegin;
  
  ptatin = (pTatinCtx)ctx;
  ierr = pTatinGetStokesContext(ptatin,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;
  
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dau,&infou);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dap,&infop);CHKERRQ(ierr);
  
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&FUloc,&FPloc);CHKERRQ(ierr);
  
  /* get the local (ghosted) entries for each physics */
  ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
  /* insert boundary conditions into local vectors */
  ierr = BCListInsertLocal(stokes->u_bclist,Uloc);CHKERRQ(ierr);
  
  ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
  
  /* compute Ax - b */
  ierr = VecZeroEntries(FUloc);CHKERRQ(ierr);
  ierr = VecZeroEntries(FPloc);CHKERRQ(ierr);
  ierr = VecGetArray(FUloc,&LA_FUloc);CHKERRQ(ierr);
  ierr = VecGetArray(FPloc,&LA_FPloc);CHKERRQ(ierr);
  
  /* ======================================== */
  /*         UPDATE NON-LINEARITIES           */
  /* evaluate rheology and rhs using X        */
  /* map marker eta to quadrature points */
  /* map marker force to quadrature points */
  /* ======================================== */
  ierr = pTatin_EvaluateRheologyNonlinearities(ptatin,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
  
  /* momentum + continuity + forces */
  ierr = FormFunctionLocal_UP(stokes->volQ,dau,LA_Uloc,dap,LA_Ploc,LA_FUloc,LA_FPloc);CHKERRQ(ierr);
  //ierr = FormFunctionLocal_U_tractionBC(stokes,dau,LA_Uloc,dap,LA_Ploc,LA_FUloc);CHKERRQ(ierr);

  ierr = VecRestoreArray(FPloc,&LA_FPloc);CHKERRQ(ierr);
  ierr = VecRestoreArray(FUloc,&LA_FUloc);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
  ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  
  /* do global fem summation */
  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  ierr = DMCompositeGather(stokes_pack,ADD_VALUES,F,FUloc,FPloc);CHKERRQ(ierr);
  
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&FUloc,&FPloc);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
  
  /* modify F for the boundary conditions, F_k = scale_k(x_k - phi_k) */
  ierr = DMCompositeGetAccess(stokes_pack,F,&Fu,&Fp);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,X,&u,&p);CHKERRQ(ierr);
  
  ierr = BCListResidualDirichlet(stokes->u_bclist,u,Fu);CHKERRQ(ierr);
  
  ierr = DMCompositeRestoreAccess(stokes_pack,X,&u,&p);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(stokes_pack,F,&Fu,&Fp);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#endif /* define(__AVX__) */

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
 **    filename:   stokes_operators_tensor.c
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
// -*- indent-tabs-mode:t c-basic-offset:8 -*-

#include <petscfe.h>
#include <ptatin3d.h>
#include <ptatin3d_stokes.h>
#include <dmda_element_q2p1.h>
#include <stokes_operators.h>

/*
 * Performs three tensor contractions: y[l,a,b,c] += T[a,k] S[b,j] R[c,i] x[l,k,j,i]
 */
static PetscErrorCode TensorContractNEV(PetscReal Rf[][3],PetscReal Sf[][3],PetscReal Tf[][3],GradMode gmode,PetscReal x[][NQP][NEV],PetscReal y[][NQP][NEV])
{
  PetscReal R[3][3],S[3][3],T[3][3];
  PetscReal u[3][NQP][NEV],v[3][NQP][NEV];
    PetscInt i,j,k,l,kj,ji,a,b,c,e;

  PetscFunctionBegin;
  for (j=0; j<3; j++) {
    for (i=0; i<3; i++) {
      R[i][j] = i<3 ? (gmode == GRAD ? Rf[i][j] : Rf[j][i]) : 0.;
      S[i][j] = i<3 ? (gmode == GRAD ? Sf[i][j] : Sf[j][i]) : 0.;
      T[i][j] = i<3 ? (gmode == GRAD ? Tf[i][j] : Tf[j][i]) : 0.;
    }
  }

  // u[l,k,j,c] = R[c,i] x[l,k,j,i]
  PetscMemzero(u,sizeof u);
  for (l=0; l<3; l++) {
    for (kj=0; kj<9; kj++) {
      for (i=0; i<3; i++) {
        for (c=0; c<3; c++) {
          for (e=0; e<NEV; e++) u[l][kj*3+c][e] += R[c][i] * x[l][kj*3+i][e];
        }
      }
    }
  }

  // v[l,k,b,c] = S[b,j] u[l,k,j,c]
  PetscMemzero(v,sizeof v);
  for (l=0; l<3; l++) {
    for (k=0; k<3; k++) {
      for (j=0; j<3; j++) {
        for (c=0; c<3; c++) {
          for (b=0; b<3; b++) {
            for (e=0; e<NEV; e++) v[l][(k*3+b)*3+c][e] += S[b][j] * u[l][(k*3+j)*3+c][e];
          }
        }
      }
    }
  }

  // y[l,a,b,c] = T[a,k] v[l,k,b,c]
  for (k=0; k<3; k++) {
    for (l=0; l<3; l++) {
      for (a=0; a<3; a++) {
        for (ji=0; ji<9; ji++) {
          for (e=0; e<NEV; e++) y[l][a*9+ji][e] += T[a][k] * v[l][k*9+ji][e];
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode JacobianInvertNEV(PetscScalar dx[3][3][NQP][NEV],PetscScalar dxdet[NQP][NEV])
{
  PetscInt i,j,k,e;

  for (i=0; i<NQP; i++) {
    PetscScalar a[3][3][NEV];
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

static PetscErrorCode QuadratureAction(const QPntVolCoefStokes *gausspt[],
               PetscScalar dx[3][3][Q2_NODES_PER_EL_3D][NEV],
               PetscScalar dxdet[Q2_NODES_PER_EL_3D][NEV],
               PetscReal w[Q2_NODES_PER_EL_3D],
               PetscScalar du[3][3][Q2_NODES_PER_EL_3D][NEV],
               PetscScalar dv[3][3][Q2_NODES_PER_EL_3D][NEV])
{
  PetscInt i,l,k,e;

  for (i=0; i<NQP; i++) {
    PetscScalar Du[6][NEV],Dv[6][NEV]; /* Symmetric gradient with respect to physical coordinates, xx, yy, zz, xy+yx, xz+zx, yz+zy */
    for (e=0; e<NEV; e++) {
      PetscScalar dux[3][3];
      for (l=0; l<3; l++) { // fields
        for (k=0; k<3; k++) { // directions
          dux[k][l] = du[0][l][i][e] * dx[k][0][i][e] + du[1][l][i][e] * dx[k][1][i][e] + du[2][l][i][e] * dx[k][2][i][e];
        }
      }
      Du[0][e] = dux[0][0];
      Du[1][e] = dux[1][1];
      Du[2][e] = dux[2][2];
      Du[3][e] = 0.5*(dux[0][1] + dux[1][0]);
      Du[4][e] = 0.5*(dux[0][2] + dux[2][0]);
      Du[5][e] = 0.5*(dux[1][2] + dux[2][1]);
    }

    for (e=0; e<NEV; e++) {
      for (k=0; k<6; k++) { /* Stress is coefficient of test function */
        Dv[k][e] = 2 * gausspt[e][i].eta * Du[k][e];
      }
    }

    for (e=0; e<NEV; e++) {
      PetscScalar dvx[3][3];
      dvx[0][0] = Dv[0][e];
      dvx[0][1] = Dv[3][e];
      dvx[0][2] = Dv[4][e];
      dvx[1][0] = Dv[3][e];
      dvx[1][1] = Dv[1][e];
      dvx[1][2] = Dv[5][e];
      dvx[2][0] = Dv[4][e];
      dvx[2][1] = Dv[5][e];
      dvx[2][2] = Dv[2][e];

      for (l=0; l<3; l++) { // fields
        for (k=0; k<3; k++) { // directions
          dv[k][l][i][e] = w[i] * dxdet[i][e] * (dvx[0][l] * dx[0][k][i][e] + dvx[1][l] * dx[1][k][i][e] + dvx[2][l] * dx[2][k][i][e]);
        }
      }
    }
  }
  return 0;
}

PetscErrorCode MFStokesWrapper_A11_Tensor(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{
  PetscErrorCode ierr;
  DM cda;
  Vec gcoords;
  const PetscReal *LA_gcoords;
  PetscInt nel,nen_u;
  const PetscInt *elnidx_u;
  QPntVolCoefStokes *all_gausspoints;
  PetscReal x1[3],w1[3],B[3][3],D[3][3],w[NQP];
  PetscInt i,j,k,e;

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
        w[(i*3+j)*3+k] = w1[i] * w1[j] * w1[k];}}}

  /* setup for coords */
  ierr = DMGetCoordinateDM( dau, &cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal( dau,&gcoords );CHKERRQ(ierr);
  ierr = VecGetArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);

  ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);

  ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);

#if defined(_OPENMP)
  #define OPENMP_CHKERRQ(x)
#else
  #define OPENMP_CHKERRQ(x)   CHKERRQ(x)
#endif

#if defined(_OPENMP)
    #pragma omp parallel for private(i)
#endif
  for (e=0;e<nel;e+=NEV) {
    PetscScalar elu[3][Q2_NODES_PER_EL_3D][NEV]={},elx[3][Q2_NODES_PER_EL_3D][NEV]={},elv[3][Q2_NODES_PER_EL_3D][NEV];
    PetscScalar dx[3][3][NQP][NEV],dxdet[NQP][NEV],du[3][3][NQP][NEV],dv[3][3][NQP][NEV];
    const QPntVolCoefStokes *cell_gausspoints[NEV];
    PetscInt ee,l;

    for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
      for (ee=0; ee<NEV; ee++) {
        PetscInt E = elnidx_u[nen_u*PetscMin(e+ee,nel-1)+i]; // Pad up to length NEV by duplicating last element
        for (l=0; l<3; l++) {
          elx[l][i][ee] = LA_gcoords[3*E+l];
          elu[l][i][ee] = ufield[3*E+l];
        }
      }
    }
    for (ee=0; ee<NEV; ee++) {
      ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,PetscMin(e+ee,nel-1),(QPntVolCoefStokes**)&cell_gausspoints[ee]);OPENMP_CHKERRQ(ierr);
    }

    ierr = PetscMemzero(dx,sizeof dx);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(D,B,B,GRAD,elx,dx[0]);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(B,D,B,GRAD,elx,dx[1]);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(B,B,D,GRAD,elx,dx[2]);OPENMP_CHKERRQ(ierr);

    ierr = JacobianInvertNEV(dx,dxdet);OPENMP_CHKERRQ(ierr);

    ierr = PetscMemzero(du,sizeof du);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(D,B,B,GRAD,elu,du[0]);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(B,D,B,GRAD,elu,du[1]);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(B,B,D,GRAD,elu,du[2]);OPENMP_CHKERRQ(ierr);

    ierr = QuadratureAction(cell_gausspoints,dx,dxdet,w,du,dv);OPENMP_CHKERRQ(ierr);

    ierr = PetscMemzero(elv,sizeof elv);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(D,B,B,GRAD_TRANSPOSE,dv[0],elv);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(B,D,B,GRAD_TRANSPOSE,dv[1],elv);OPENMP_CHKERRQ(ierr);
    ierr = TensorContractNEV(B,B,D,GRAD_TRANSPOSE,dv[2],elv);OPENMP_CHKERRQ(ierr);

#if defined(_OPENMP)
    #pragma omp critical
#endif
    for (ee=0; ee<PetscMin(NEV,nel-e); ee++) {
      for (i=0; i<NQP; i++) {
        PetscInt E = elnidx_u[nen_u*(e+ee)+i];
        for (l=0; l<3; l++) {
          Yu[3*E+l] += elv[l][i][ee];
        }
      }
    }

  }

#undef OPENMP_CHKERRQ
    PetscLogFlops((nel * 9) * 3*NQP*(6+6+6));           /* 9 tensor contractions per element */
    PetscLogFlops(nel*NQP*(14 + 1/* division */ + 27)); /* 1 Jacobi inversion per element */
    PetscLogFlops(nel*NQP*(5*9+6+6+6*9));               /* 1 quadrature action per element */

  ierr = VecRestoreArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

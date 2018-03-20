

#include <petsc.h>
#include <petscvec.h>
#include <ptatin3d.h>
#include <ptatin3d_defs.h>
#include <private/ptatin_impl.h>
#include <element_type_Q2.h>
#include <element_utils_q2.h>
#include <dmda_element_q2p1.h>


PetscErrorCode Evaluate_uL2_uH1_pL2(DM dav,DM dap,Vec sv,Vec sp,Quadrature volQ,
                                    PetscErrorCode (*exact_v)(PetscReal*,PetscReal*,void*),
                                    PetscErrorCode (*exact_L)(PetscReal*,PetscReal*,void*),
                                    PetscErrorCode (*exact_p)(PetscReal*,PetscReal*,void*),
                                    void *data,
                                    PetscReal *ev_L2,PetscReal *ev_H1s,PetscReal *ev_H1,PetscReal *ep_L2)
{
  DM              cda;
  Vec             gcoords,sv_local,sp_local;
  PetscScalar     *LA_gcoords,*LA_sv,*LA_sp;
  PetscInt        nel,nen_u,nen_p,e,k;
  const PetscInt  *elnidx_u;
  const PetscInt  *elnidx_p;
  PetscReal       el_coords[3*Q2_NODES_PER_EL_3D],el_v[3*Q2_NODES_PER_EL_3D],el_p[P_BASIS_FUNCTIONS];
  PetscReal       ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
  PetscInt        q,nqp;
  PetscReal       WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
  PetscReal       detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
  PetscReal         norms_l[4],norms_g[4];
  PetscErrorCode    ierr;


  PetscFunctionBegin;

  /* setup quadrature */
  nqp = volQ->npoints;
  P3D_prepare_elementQ2(nqp,WEIGHT,XI,NI,GNI);

  /* setup local coords */
  ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dav,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

  /* setup local velocity */
  ierr = DMGetLocalVector(dav,&sv_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dav,sv,INSERT_VALUES,sv_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(  dav,sv,INSERT_VALUES,sv_local);CHKERRQ(ierr);
  ierr = VecGetArray(sv_local,&LA_sv);CHKERRQ(ierr);

  /* setup local pressure */
  ierr = DMGetLocalVector(dap,&sp_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dap,sp,INSERT_VALUES,sp_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(  dap,sp,INSERT_VALUES,sp_local);CHKERRQ(ierr);
  ierr = VecGetArray(sp_local,&LA_sp);CHKERRQ(ierr);

  ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
  ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);

  norms_l[0] = norms_l[1] = norms_l[2] = norms_l[3] = 0.0;
  for (e=0; e<nel; e++) {

    ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);

    ierr = DMDAGetVectorElementFieldQ2_3D(el_v,(PetscInt*)&elnidx_u[nen_u*e],LA_sv);CHKERRQ(ierr);
    for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
      ux[k] = el_v[3*k  ];
      uy[k] = el_v[3*k+1];
      uz[k] = el_v[3*k+2];
    }

    ierr = DMDAGetScalarElementField(el_p,nen_p,(PetscInt*)&elnidx_p[nen_p*e],LA_sp);CHKERRQ(ierr);

    for (q=0; q<nqp; q++) {
      PetscScalar xip[] = { XI[q][0], XI[q][1], XI[q][2] };
      ConstructNi_pressure(xip,el_coords,NIp[q]);
    }

    P3D_evaluate_geometry_elementQ2(nqp,el_coords,GNI,detJ,dNudx,dNudy,dNudz);

    for (q=0; q<nqp; q++) {
      PetscReal p_qp,v_qp[3],L_qp[3*3],integrand[4],x_qp[3];
      PetscReal p_ex,v_ex[3],L_ex[3*3];
      PetscInt ii,jj;

      /* coordinate */
      x_qp[0] = x_qp[1] = x_qp[2] = 0.0;
      for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
        x_qp[0] += NI[q][k] * el_coords[3*k+0];
        x_qp[1] += NI[q][k] * el_coords[3*k+1];
        x_qp[2] += NI[q][k] * el_coords[3*k+2];
      }

      /* velocity */
      v_qp[0] = v_qp[1] = v_qp[2] = 0.0;
      for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
        v_qp[0] += NI[q][k] * ux[k];
        v_qp[1] += NI[q][k] * uy[k];
        v_qp[2] += NI[q][k] * uz[k];
      }

      /* pressure */
      p_qp = 0.0;
      for (k=0; k<P_BASIS_FUNCTIONS; k++) {
        p_qp += NIp[q][k] * el_p[k];
      }

      /* velocity gradient, L_ij = u_{i,j} */
      for (ii=0; ii<3; ii++) {
        for (jj=0; jj<3; jj++) {
          L_qp[ii*3+jj] = 0.0;
        }
      }

      for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
        L_qp[0*3+0] += dNudx[q][k] * ux[k];
        L_qp[0*3+1] += dNudy[q][k] * ux[k];
        L_qp[0*3+2] += dNudz[q][k] * ux[k];

        L_qp[1*3+0] += dNudx[q][k] * uy[k];
        L_qp[1*3+1] += dNudy[q][k] * uy[k];
        L_qp[1*3+2] += dNudz[q][k] * uy[k];

        L_qp[2*3+0] += dNudx[q][k] * uz[k];
        L_qp[2*3+1] += dNudy[q][k] * uz[k];
        L_qp[2*3+2] += dNudz[q][k] * uz[k];
      }

      /* compute analytic funcitons at x_qp */
      v_ex[0] = v_ex[1] = v_ex[2] = 0.0;
      if (exact_v) {
        ierr = exact_v(x_qp,v_ex,data);CHKERRQ(ierr);
      }

      for (ii=0; ii<3; ii++) {
        for (jj=0; jj<3; jj++) {
          L_ex[ii*3+jj] = 0.0;
        }
      }
      if (exact_L) {
        ierr = exact_L(x_qp,L_ex,data);CHKERRQ(ierr);
      }

      p_ex = 0.0;
      if (exact_p) {
        ierr = exact_p(x_qp,&p_ex,data);CHKERRQ(ierr);
      }

      /* compute errors */
      integrand[0] = 0.0;
      for (ii=0; ii<3; ii++) {
        integrand[0] += (v_qp[ii] - v_ex[ii])*(v_qp[ii] - v_ex[ii]);
      }

      integrand[1] = 0.0;
      for (ii=0; ii<3; ii++) {
        for (jj=0; jj<3; jj++) {
          integrand[1] += (L_qp[ii*3+jj] - L_ex[ii*3+jj])*(L_qp[ii*3+jj] - L_ex[ii*3+jj]);
        }
      }
      integrand[2] = integrand[0] + integrand[1];
      integrand[3] = (p_qp - p_ex)*(p_qp - p_ex);

      norms_l[0] += WEIGHT[q] * integrand[0] * detJ[q];
      norms_l[1] += WEIGHT[q] * integrand[1] * detJ[q];
      norms_l[2] += WEIGHT[q] * integrand[2] * detJ[q];
      norms_l[3] += WEIGHT[q] * integrand[3] * detJ[q];
    }
  }

  ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  ierr = VecRestoreArray(sv_local,&LA_sv);CHKERRQ(ierr);
  ierr = VecRestoreArray(sp_local,&LA_sp);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dav,&sv_local);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dap,&sp_local);CHKERRQ(ierr);

  ierr = MPI_Allreduce(norms_l,norms_g,4,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)dav));CHKERRQ(ierr);

  *ev_L2  = PetscSqrtReal(norms_g[0]);
  *ev_H1s = PetscSqrtReal(norms_g[1]);
  *ev_H1  = PetscSqrtReal(norms_g[2]);
  *ep_L2  = PetscSqrtReal(norms_g[3]);

  PetscFunctionReturn(0);
}

PetscErrorCode Evaluate_uL2_symuH1_pL2(DM dav,DM dap,Vec sv,Vec sp,Quadrature volQ,
                                    PetscErrorCode (*exact_v)(PetscReal*,PetscReal*,void*),
                                    PetscErrorCode (*exact_E)(PetscReal*,PetscReal*,void*),
                                    PetscErrorCode (*exact_p)(PetscReal*,PetscReal*,void*),
                                    void *data,
                                    PetscReal *ev_L2,PetscReal *ev_H1s,PetscReal *ev_H1,PetscReal *ep_L2)
{
  DM              cda;
  Vec             gcoords,sv_local,sp_local;
  PetscScalar     *LA_gcoords,*LA_sv,*LA_sp;
  PetscInt        nel,nen_u,nen_p,e,k;
  const PetscInt  *elnidx_u;
  const PetscInt  *elnidx_p;
  PetscReal       el_coords[3*Q2_NODES_PER_EL_3D],el_v[3*Q2_NODES_PER_EL_3D],el_p[P_BASIS_FUNCTIONS];
  PetscReal       ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
  PetscInt        q,nqp;
  PetscReal       WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE],NIp[NQP][P_BASIS_FUNCTIONS];
  PetscReal       detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
  PetscReal         norms_l[4],norms_g[4];
  PetscErrorCode    ierr;


  PetscFunctionBegin;

  /* setup quadrature */
  nqp = volQ->npoints;
  P3D_prepare_elementQ2(nqp,WEIGHT,XI,NI,GNI);

  /* setup local coords */
  ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dav,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

  /* setup local velocity */
  ierr = DMGetLocalVector(dav,&sv_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dav,sv,INSERT_VALUES,sv_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(  dav,sv,INSERT_VALUES,sv_local);CHKERRQ(ierr);
  ierr = VecGetArray(sv_local,&LA_sv);CHKERRQ(ierr);

  /* setup local pressure */
  ierr = DMGetLocalVector(dap,&sp_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dap,sp,INSERT_VALUES,sp_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(  dap,sp,INSERT_VALUES,sp_local);CHKERRQ(ierr);
  ierr = VecGetArray(sp_local,&LA_sp);CHKERRQ(ierr);

  ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
  ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);

  norms_l[0] = norms_l[1] = norms_l[2] = norms_l[3] = 0.0;
  for (e=0; e<nel; e++) {

    ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);

    ierr = DMDAGetVectorElementFieldQ2_3D(el_v,(PetscInt*)&elnidx_u[nen_u*e],LA_sv);CHKERRQ(ierr);
    for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
      ux[k] = el_v[3*k  ];
      uy[k] = el_v[3*k+1];
      uz[k] = el_v[3*k+2];
    }

    ierr = DMDAGetScalarElementField(el_p,nen_p,(PetscInt*)&elnidx_p[nen_p*e],LA_sp);CHKERRQ(ierr);

    for (q=0; q<nqp; q++) {
      PetscScalar xip[] = { XI[q][0], XI[q][1], XI[q][2] };
      ConstructNi_pressure(xip,el_coords,NIp[q]);
    }

    P3D_evaluate_geometry_elementQ2(nqp,el_coords,GNI,detJ,dNudx,dNudy,dNudz);

    for (q=0; q<nqp; q++) {
      PetscReal p_qp,v_qp[3],E_qp[3*3],integrand[4],x_qp[3];
      PetscReal p_ex,v_ex[3],E_ex[3*3];
      PetscInt ii,jj;

      /* coordinate */
      x_qp[0] = x_qp[1] = x_qp[2] = 0.0;
      for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
        x_qp[0] += NI[q][k] * el_coords[3*k+0];
        x_qp[1] += NI[q][k] * el_coords[3*k+1];
        x_qp[2] += NI[q][k] * el_coords[3*k+2];
      }

      /* velocity */
      v_qp[0] = v_qp[1] = v_qp[2] = 0.0;
      for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
        v_qp[0] += NI[q][k] * ux[k];
        v_qp[1] += NI[q][k] * uy[k];
        v_qp[2] += NI[q][k] * uz[k];
      }

      /* pressure */
      p_qp = 0.0;
      for (k=0; k<P_BASIS_FUNCTIONS; k++) {
        p_qp += NIp[q][k] * el_p[k];
      }

      /* symmetric velocity gradient, E_ij = u_{i,j} */
      for (ii=0; ii<3; ii++) {
        for (jj=0; jj<3; jj++) {
          E_qp[ii*3+jj] = 0.0;
        }
      }

      for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
        E_qp[0*3+0] += dNudx[q][k] * ux[k];
        E_qp[0*3+1] += 0.5*(dNudy[q][k] * ux[k] + dNudx[q][k] * uy[k]);
        E_qp[0*3+2] += 0.5*(dNudz[q][k] * ux[k] + dNudx[q][k] * uz[k]);

        E_qp[1*3+0] += 0.5*(dNudx[q][k] * uy[k] + dNudy[q][k] * ux[k]);
        E_qp[1*3+1] += dNudy[q][k] * uy[k];
        E_qp[1*3+2] += 0.5*(dNudz[q][k] * uy[k] + dNudy[q][k] * uz[k]);

        E_qp[2*3+0] += 0.5*(dNudx[q][k] * uz[k] + dNudz[q][k] * ux[k]);
        E_qp[2*3+1] += 0.5*(dNudy[q][k] * uz[k] + dNudz[q][k] * uy[k]);
        E_qp[2*3+2] += dNudz[q][k] * uz[k];
      }

      /* compute analytic funcitons at x_qp */
      v_ex[0] = v_ex[1] = v_ex[2] = 0.0;
      if (exact_v) {
        ierr = exact_v(x_qp,v_ex,data);CHKERRQ(ierr);
        //printf("n (%+1.4e,%+1.4e,%+1.4e) : e (%+1.4e,%+1.4e,%+1.4e)\n",v_qp[0],v_qp[1],v_qp[2],v_ex[0],v_ex[1],v_ex[2]);
        //printf("%+1.4e %+1.4e %+1.4e   %+1.4e %+1.4e %+1.4e    %+1.4e %+1.4e %+1.4e\n",x_qp[0],x_qp[1],x_qp[2],v_qp[0],v_qp[1],v_qp[2],v_ex[0],v_ex[1],v_ex[2]);
      }

      for (ii=0; ii<3; ii++) {
        for (jj=0; jj<3; jj++) {
          E_ex[ii*3+jj] = 0.0;
        }
      }
      if (exact_E) {
        ierr = exact_E(x_qp,E_ex,data);CHKERRQ(ierr);
      }

      p_ex = 0.0;
      if (exact_p) {
        ierr = exact_p(x_qp,&p_ex,data);CHKERRQ(ierr);
      }

      /* compute errors */
      integrand[0] = 0.0;
      for (ii=0; ii<3; ii++) {
        integrand[0] += (v_qp[ii] - v_ex[ii])*(v_qp[ii] - v_ex[ii]);
      }

      integrand[1] = 0.0;
      for (ii=0; ii<3; ii++) {
        for (jj=0; jj<3; jj++) {
          integrand[1] += (E_qp[ii*3+jj] - E_ex[ii*3+jj])*(E_qp[ii*3+jj] - E_ex[ii*3+jj]);
        }
      }
      integrand[2] = integrand[0] + integrand[1];
      integrand[3] = (p_qp - p_ex)*(p_qp - p_ex);

      norms_l[0] += WEIGHT[q] * integrand[0] * detJ[q];
      norms_l[1] += WEIGHT[q] * integrand[1] * detJ[q];
      norms_l[2] += WEIGHT[q] * integrand[2] * detJ[q];
      norms_l[3] += WEIGHT[q] * integrand[3] * detJ[q];
    }
  }

  ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  ierr = VecRestoreArray(sv_local,&LA_sv);CHKERRQ(ierr);
  ierr = VecRestoreArray(sp_local,&LA_sp);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dav,&sv_local);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dap,&sp_local);CHKERRQ(ierr);

  ierr = MPI_Allreduce(norms_l,norms_g,4,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)dav));CHKERRQ(ierr);

  *ev_L2  = PetscSqrtReal(norms_g[0]);
  *ev_H1s = PetscSqrtReal(norms_g[1]);
  *ev_H1  = PetscSqrtReal(norms_g[2]);
  *ep_L2  = PetscSqrtReal(norms_g[3]);

  PetscFunctionReturn(0);
}

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
 **    filename:   ptatin3d_stokes.c
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
#include "dmda_element_q2p1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"
#include "element_type_Q2.h"
#include "data_bucket.h"

#include "petsc/private/snesimpl.h" /* for snes->ttol */

#include "QPntVolCoefStokes_def.h"
#include "QPntSurfCoefStokes_def.h"


PetscErrorCode StokesVelocity_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[])
{
  PetscInt n;

  PetscFunctionBegin;
  for (n=0; n<27; n++) {
    el_localIndices[3*n  ] = 3*elnid[n]  ;
    el_localIndices[3*n+1] = 3*elnid[n]+1;
    el_localIndices[3*n+2] = 3*elnid[n]+2;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode StokesPressure_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[])
{
  PetscInt n;

  PetscFunctionBegin;
  for (n=0; n<P_BASIS_FUNCTIONS; n++) {
    el_localIndices[n] = elnid[n];
  }
  PetscFunctionReturn(0);
}
PetscErrorCode StokesVelocityScalar_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[])
{
  PetscInt n;

  PetscFunctionBegin;
  for (n=0; n<27; n++) {
    el_localIndices[n] = elnid[n] ;
  }
  PetscFunctionReturn(0);
}


/* physics component loader */
PetscErrorCode PhysCompCreate_Stokes(PhysCompStokes *ctx)
{
  PetscErrorCode ierr;
  PhysCompStokes stokes;

  PetscFunctionBegin;
  ierr = PetscMalloc(sizeof(struct _p_PhysCompStokes),&stokes);CHKERRQ(ierr);
  ierr = PetscMemzero(stokes,sizeof(struct _p_PhysCompStokes));CHKERRQ(ierr);
  *ctx = stokes;
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompDestroy_Stokes(PhysCompStokes *ctx)
{
  PetscErrorCode ierr;
  PhysCompStokes user;
  PetscInt e;

  PetscFunctionBegin;

  if (!ctx) {PetscFunctionReturn(0);}
  user = *ctx;

  if (user->surfQ) {
    for (e=0; e<HEX_EDGES; e++) {
      if (user->surfQ[e]) {
        ierr = SurfaceQuadratureDestroy(&user->surfQ[e]);CHKERRQ(ierr);
        user->surfQ[e] = NULL;
      }
    }
    ierr = PetscFree(user->surfQ);CHKERRQ(ierr);
  }

  if (user->volQ) { ierr = QuadratureDestroy(&user->volQ);CHKERRQ(ierr); }
  if (user->p_bclist) { ierr = BCListDestroy(&user->p_bclist);CHKERRQ(ierr); }
  if (user->u_bclist) { ierr = BCListDestroy(&user->u_bclist);CHKERRQ(ierr); }
  if (user->stokes_pack) { ierr = DMDestroy(&user->stokes_pack);CHKERRQ(ierr); }
  if (user->dap) { ierr = DMDestroy(&user->dap);CHKERRQ(ierr); }
  if (user->dav) { ierr = DMDestroy(&user->dav);CHKERRQ(ierr); }
  if (user) { ierr = PetscFree(user);CHKERRQ(ierr); }

  *ctx = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompCreateMesh_Stokes3d(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx)
{
  DM dav,dap,multipys_pack;
  PetscInt vbasis_dofs;
  PetscInt pbasis_dofs;
  const PetscInt *lxp,*lyp,*lzp;
  PetscInt MX,MY,MZ,Mp,Np,Pp,*lxv,*lyv,*lzv;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  MX = mx;
  MY = my;
  MZ = mz;

  /* velocity */
  vbasis_dofs = 3;
  ierr = DMDACreate3d( PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, 2*MX+1,2*MY+1,2*MZ+1, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, vbasis_dofs,2, NULL,NULL,NULL, &dav );CHKERRQ(ierr);
  ierr = DMSetUp(dav);CHKERRQ(ierr);
  ierr = DMDASetElementType_Q2(dav);CHKERRQ(ierr);
  ierr = DMSetMatType(dav,MATSBAIJ);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dav,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRangesElementQ2(dav, 0,0,0, 0,0,0, &lxv,&lyv,&lzv);CHKERRQ(ierr);

  /* pressure */
  pbasis_dofs = P_BASIS_FUNCTIONS;
  ierr = DMDACreate3d( PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX,MY,MZ, Mp,Np,Pp, pbasis_dofs,0, lxv,lyv,lzv, &dap );CHKERRQ(ierr);
  ierr = DMSetUp(dap);CHKERRQ(ierr);
  ierr = DMDASetElementType_P1(dap);CHKERRQ(ierr);
  ierr = DMSetMatType(dap,MATSBAIJ);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRanges(dap,&lxp,&lyp,&lzp);CHKERRQ(ierr);

  /* set an initial geometry */
  ierr = DMDASetUniformCoordinates(dav,0.0,1.0, 0.0,1.0, 0.0,1.0);CHKERRQ(ierr);

  /* stokes */
  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&multipys_pack);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(multipys_pack,dav);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(multipys_pack,dap);CHKERRQ(ierr);

  ierr = DMDASetFieldName(dav,0,"ux");CHKERRQ(ierr);
  ierr = DMDASetFieldName(dav,1,"uy");CHKERRQ(ierr);
  ierr = DMDASetFieldName(dav,2,"uz");CHKERRQ(ierr);
  switch (P_BASIS_FUNCTIONS) {
    case 1:
      ierr = DMDASetFieldName(dap,0,"P0_p");CHKERRQ(ierr);
      break;
    case 4:
      ierr = DMDASetFieldName(dap,0,"P1_p");CHKERRQ(ierr);
      ierr = DMDASetFieldName(dap,1,"P1_dpdx");CHKERRQ(ierr);
      ierr = DMDASetFieldName(dap,2,"P1_dpdy");CHKERRQ(ierr);
      ierr = DMDASetFieldName(dap,3,"P1_dpdz");CHKERRQ(ierr);
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Pressure space may ONLY contain 1 or 4 basis functions");
      break;
  }
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dav,"stk_velocity_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dap,"stk_pressure_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)multipys_pack,"stk_pack_");CHKERRQ(ierr);

  ierr = DMSetFromOptions(dav);CHKERRQ(ierr);
  ierr = DMSetFromOptions(dap);CHKERRQ(ierr);
  ctx->dav  = dav;
  ctx->dap  = dap;
  ctx->stokes_pack = multipys_pack;

  ierr = PetscFree(lxv);CHKERRQ(ierr);
  ierr = PetscFree(lyv);CHKERRQ(ierr);
  ierr = PetscFree(lzv);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompCreateBoundaryList_Stokes(PhysCompStokes ctx)
{
  DM dav;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* vel bc's */
  dav = ctx->dav;
  if (!dav) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"dav must be set"); }

  ierr = DMDABCListCreate(dav,&ctx->u_bclist);CHKERRQ(ierr);

  /* pressure bc's */
  ctx->p_bclist = NULL;

  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompCreateVolumeQuadrature_Stokes(PhysCompStokes ctx)
{
  DM dav;
  PetscInt lmx,lmy,lmz;
  PetscInt np_per_dim,ncells;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  dav = ctx->dav;

  np_per_dim = 3;
    ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
  ncells = lmx * lmy * lmz;
  ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&ctx->volQ);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesSetGravityUnitVector(PhysCompStokes ctx,PetscReal grav[])
{
  PetscReal      norm_g;

  PetscFunctionBegin;
  norm_g = PetscSqrtScalar(grav[0]*grav[0] + grav[1]*grav[1] + grav[2]*grav[2]);
  ctx->gravity_vector[0] = grav[0]/norm_g;
  ctx->gravity_vector[1] = grav[1]/norm_g;
  ctx->gravity_vector[2] = grav[2]/norm_g;
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesScaleGravityVector(PhysCompStokes ctx,PetscReal fac)
{
  PetscFunctionBegin;
  ctx->gravity_vector[0] = ctx->gravity_vector[0]*fac;
  ctx->gravity_vector[1] = ctx->gravity_vector[1]*fac;
  ctx->gravity_vector[2] = ctx->gravity_vector[2]*fac;
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesSetGravityVector(PhysCompStokes ctx,PetscReal grav[])
{
  PetscFunctionBegin;
  ctx->gravity_vector[0] = grav[0];
  ctx->gravity_vector[1] = grav[1];
  ctx->gravity_vector[2] = grav[2];
  PetscFunctionReturn(0);
}

PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[])
{
  PetscInt n,idx;

  PetscFunctionBegin;
  for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
    idx = u_eqn[3*n  ];
    fields_F[idx] += Fe_u[NSD*n  ];

    idx = u_eqn[3*n+1];
    fields_F[idx] += Fe_u[NSD*n+1];

    idx = u_eqn[3*n+2];
    fields_F[idx] += Fe_u[NSD*n+2];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMDASetValuesLocalStencil_InsertValues_Stokes_Velocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[])
{
  PetscInt n,idx;

  PetscFunctionBegin;
  for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
    idx = u_eqn[3*n  ];
    fields_F[idx] = Fe_u[NSD*n  ];

    idx = u_eqn[3*n+1];
    fields_F[idx] = Fe_u[NSD*n+1];

    idx = u_eqn[3*n+2];
    fields_F[idx] = Fe_u[NSD*n+2];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(PetscScalar *fields_F,PetscInt p_eqn[],PetscScalar Fe_p[])
{
  PetscInt n,idx;

  PetscFunctionBegin;
  for (n = 0; n<P_BASIS_FUNCTIONS; n++) {
    idx = p_eqn[n];
    fields_F[idx] += Fe_p[n];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_ScalarVelocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[])
{
  PetscInt n,idx;

  PetscFunctionBegin;
  for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
    idx = u_eqn[n];
    fields_F[idx] += Fe_u[n];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompSetup_Stokes(PhysCompStokes ctx,DM dav)
{
  DM dap,multipys_pack;
  PetscInt pbasis_dofs;
  PetscInt Mp,Np,Pp,*lxv,*lyv,*lzv,MX,MY,MZ;
  const PetscInt *lxp,*lyp,*lzp;
  MPI_Comm comm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)dav,&comm);CHKERRQ(ierr);

  /* velocity */
  ierr = DMDASetElementType_Q2(dav);CHKERRQ(ierr);
  ierr = DMSetMatType(dav,MATSBAIJ);CHKERRQ(ierr);

  ierr = DMDAGetInfo(dav,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRangesElementQ2(dav, 0,0,0, 0,0,0, &lxv,&lyv,&lzv);CHKERRQ(ierr);

  ierr = DMDAGetSizeElementQ2(dav,&ctx->mx,&ctx->my,&ctx->mz);CHKERRQ(ierr);

  /* pressure */
  MX = ctx->mx;
  MY = ctx->my;
  MZ = ctx->mz;

  pbasis_dofs = P_BASIS_FUNCTIONS;
  ierr = DMDACreate3d( comm, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX,MY,MZ, Mp,Np,Pp, pbasis_dofs,0, lxv,lyv,lzv, &dap );CHKERRQ(ierr);
  ierr = DMSetUp(dap);CHKERRQ(ierr);
  ierr = DMDASetElementType_P1(dap);CHKERRQ(ierr);
  ierr = DMSetMatType(dap,MATSBAIJ);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRanges(dap,&lxp,&lyp,&lzp);CHKERRQ(ierr);

  ierr = PetscFree(lxv);CHKERRQ(ierr);
  ierr = PetscFree(lyv);CHKERRQ(ierr);
  ierr = PetscFree(lzv);CHKERRQ(ierr);


  /* stokes */
  ierr = DMCompositeCreate(comm,&multipys_pack);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(multipys_pack,dav);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(multipys_pack,dap);CHKERRQ(ierr);

  ierr = DMDASetFieldName(dav,0,"ux");CHKERRQ(ierr);
  ierr = DMDASetFieldName(dav,1,"uy");CHKERRQ(ierr);
  ierr = DMDASetFieldName(dav,2,"uz");CHKERRQ(ierr);
  switch (P_BASIS_FUNCTIONS) {
    case 1:
      ierr = DMDASetFieldName(dap,0,"P0_p");CHKERRQ(ierr);
      break;
    case 4:
      ierr = DMDASetFieldName(dap,0,"P1_p");CHKERRQ(ierr);
      ierr = DMDASetFieldName(dap,1,"P1_dpdx");CHKERRQ(ierr);
      ierr = DMDASetFieldName(dap,2,"P1_dpdy");CHKERRQ(ierr);
      ierr = DMDASetFieldName(dap,3,"P1_dpdz");CHKERRQ(ierr);
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Pressure space may ONLY contain 1 or 4 basis functions");
      break;
  }
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dav,"stk_velocity_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dap,"stk_pressure_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)multipys_pack,"stk_pack_");CHKERRQ(ierr);

  ierr = DMSetFromOptions(dav);CHKERRQ(ierr);
  ierr = DMSetFromOptions(dap);CHKERRQ(ierr);

  ctx->dav  = dav;
  ctx->dap  = dap;
  ctx->stokes_pack = multipys_pack;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinStokesKSPMonitorBlocks(KSP ksp,PetscInt n,PetscReal rnorm,void *data)
{
  PetscErrorCode ierr;
  pTatinCtx ctx;
  PetscReal norms[4];
  Vec X,Xu,Xp,v,w;
  Mat A;

  PetscFunctionBegin;
  ctx = (pTatinCtx)data;
  ierr = KSPGetOperators(ksp,&A,0);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&w,&v);CHKERRQ(ierr);

  ierr = KSPBuildResidual(ksp,v,w,&X);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);

  ierr = VecStrideNorm(Xu,0,NORM_2,&norms[0]);CHKERRQ(ierr);
  ierr = VecStrideNorm(Xu,1,NORM_2,&norms[1]);CHKERRQ(ierr);
  ierr = VecStrideNorm(Xu,2,NORM_2,&norms[2]);CHKERRQ(ierr);
  ierr = VecNorm(Xp,NORM_2,&norms[3]);CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"%3D KSP Component U,V,W,P residual norm [ %1.12e, %1.12e, %1.12e, %1.12e ]\n",n,norms[0],norms[1],norms[2],norms[3]);

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeQuadratureCreate_GaussLegendreStokes(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature)
{
  Quadrature Q;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = QuadratureCreate(&Q);CHKERRQ(ierr);
  Q->dim  = nsd;
  Q->type = VOLUME_QUAD;

  /*PetscPrintf(PETSC_COMM_WORLD,"VolumeQuadratureCreate_GaussLegendreStokes:\n");*/
  switch (np_per_dim) {
    case 1:
      /*PetscPrintf(PETSC_COMM_WORLD,"\tUsing 1 pnt Gauss Legendre quadrature\n");*/
      //QuadratureCreateGauss_1pnt_3D(&ngp,gp_xi,gp_weight);
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"This will result in a rank-deficient operator");
      break;

    case 2:
      /*PetscPrintf(PETSC_COMM_WORLD,"\tUsing 2x2 pnt Gauss Legendre quadrature\n");*/
      QuadratureCreateGauss_2pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
      break;

    case 3:
      /*PetscPrintf(PETSC_COMM_WORLD,"\tUsing 3x3 pnt Gauss Legendre quadrature\n");*/
      QuadratureCreateGauss_3pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
      break;

    default:
      /*PetscPrintf(PETSC_COMM_WORLD,"\tUsing 3x3 pnt Gauss Legendre quadrature\n");*/
      QuadratureCreateGauss_3pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
      break;
  }

  Q->n_elements = ncells;
  if (ncells != 0) {

    DataBucketCreate(&Q->properties_db);
    DataBucketRegisterField(Q->properties_db,QPntVolCoefStokes_classname, sizeof(QPntVolCoefStokes),NULL);
    DataBucketFinalize(Q->properties_db);

    DataBucketSetInitialSizes(Q->properties_db,Q->npoints*ncells,1);

    /*
    // Note: This call will hang if any rank contained zero elements
    DataBucketView(PETSC_COMM_WORLD, Q->properties_db,"GaussLegendre StokesCoefficients",DATABUCKET_VIEW_STDOUT);
    */
  }

  *quadrature = Q;
  PetscFunctionReturn(0);
}

PetscErrorCode VolumeQuadratureGetAllCellData_Stokes(Quadrature Q,QPntVolCoefStokes *coeffs[])
{
  QPntVolCoefStokes *quadraturepoint_data;
  DataField          PField;
  PetscFunctionBegin;

  DataBucketGetDataFieldByName(Q->properties_db, QPntVolCoefStokes_classname ,&PField);
  quadraturepoint_data = PField->data;
  *coeffs = quadraturepoint_data;

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeQuadratureGetCellData_Stokes(Quadrature Q,QPntVolCoefStokes coeffs[],PetscInt cidx,QPntVolCoefStokes *cell[])
{
  PetscFunctionBegin;
  if (cidx>=Q->n_elements) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_SIZ,"cidx > max cells");
  }

  *cell = &coeffs[cidx*Q->npoints];
  PetscFunctionReturn(0);
}

/* surface quadrature */
PetscErrorCode SurfaceQuadratureCreate_GaussLegendreStokes(DM da,HexElementFace index,SurfaceQuadrature *quadrature)
{
  SurfaceQuadrature Q;
  PetscInt nface_list[HEX_EDGES];
  PetscInt lmx,lmy,lmz,nfaces,M,N,P,si,sj,sk,ni,nj,nk;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if ((int)index > HEX_EDGES) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"face index for 3d hex must be in the range [0,5]");
  }

  ierr = SurfaceQuadratureCreate(&Q);CHKERRQ(ierr);
  //Q->dim  = 3;
  //Q->type = SURFACE_QUAD;

  ierr = DMDAGetInfo(da,0,&M,&N,&P, 0,0,0,0, 0,0,0,0,0);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,&lmz);CHKERRQ(ierr);

  nface_list[0] = nface_list[1] = 0;
  nface_list[2] = nface_list[3] = 0;
  nface_list[4] = nface_list[5] = 0;
  if (si+ni == M) { nface_list[HEX_FACE_Pxi]   = lmy*lmz; }
  if (si == 0)    { nface_list[HEX_FACE_Nxi]   = lmy*lmz; }
  if (sj+nj == N) { nface_list[HEX_FACE_Peta]  = lmx*lmz; }
  if (sj == 0)    { nface_list[HEX_FACE_Neta]  = lmx*lmz; }
  if (sk+nk == P) { nface_list[HEX_FACE_Pzeta] = lmx*lmy; }
  if (sk == 0)    { nface_list[HEX_FACE_Nzeta] = lmx*lmy; }

  nfaces = nface_list[ index ];
  //Q->ncells = nfaces;

  /* setup surface gauss points (weights and local coordinates) */
  ierr = _SurfaceQuadratureCreate(Q,index,nfaces);CHKERRQ(ierr);

  /* setup element lists */
  ierr = _SurfaceQuadratureCellIndexSetUp(Q,index,nfaces,da);CHKERRQ(ierr);

  /* setup properties */
  /*
  // note: the parallel viewer wont work if db passed in is null //
  if (nfaces != 0) {
    DataBucketCreate(&Q->properties_db);
    DataBucketRegisterField(Q->properties_db,QPntSurfCoefStokes_classname, sizeof(QPntSurfCoefStokes),NULL);
    DataBucketFinalize(Q->properties_db);

    DataBucketSetInitialSizes(Q->properties_db,Q->ngp*nfaces,1);

    DataBucketView(PETSC_COMM_WORLD, Q->properties_db,"SurfaceGaussLegendre StokesCoefficients",DATABUCKET_VIEW_STDOUT);
  } else {
    Q->properties_db = NULL;
  }
  */

  DataBucketCreate(&Q->properties_db);
  DataBucketRegisterField(Q->properties_db,QPntSurfCoefStokes_classname, sizeof(QPntSurfCoefStokes),NULL);
  DataBucketFinalize(Q->properties_db);

  if (nfaces != 0) {
    DataBucketSetInitialSizes(Q->properties_db,Q->ngp*nfaces,1);
  } else {
    DataBucketSetInitialSizes(Q->properties_db,1,1);
        DataBucketSetSizes(Q->properties_db,0,-1);
  }
//  DataBucketView(PetscObjectComm((PetscObject)da), Q->properties_db,"SurfaceGaussLegendre StokesCoefficients",DATABUCKET_VIEW_STDOUT);

  *quadrature = Q;
  PetscFunctionReturn(0);
}



PetscErrorCode SurfaceQuadratureGeometrySetUpStokes(SurfaceQuadrature Q,DM da)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  //ierr = SurfaceQuadratureStokesCoordinatesSetUp(Q,da);CHKERRQ(ierr); // not storing coordinates //
  ierr = SurfaceQuadratureOrientationSetUpStokes(Q,da);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SurfaceQuadratureOrientationSetUpStokes(SurfaceQuadrature Q,DM da)
{
  PetscErrorCode ierr;
  DM             cda;
  Vec            gcoords;
  PetscScalar    *LA_gcoords;
  PetscInt       nel,nen,fe,e,k,gp;
  const PetscInt *elnidx;
  ConformingElementFamily element;
  double         elcoords[3*Q2_NODES_PER_EL_3D];
  double         Ni[27];
  QPntSurfCoefStokes *all_qpoint;
  QPntSurfCoefStokes *cell_qpoint;

  PetscFunctionBegin;

  /* setup for coords */
  ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

  ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
  element = Q->e;

  ierr = SurfaceQuadratureGetAllCellData_Stokes(Q,&all_qpoint);CHKERRQ(ierr);
  for (fe=0; fe<Q->nfaces; fe++) {

    e = Q->element_list[fe];
    ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);

    ierr =  SurfaceQuadratureGetCellData_Stokes(Q,all_qpoint,fe,&cell_qpoint);CHKERRQ(ierr);

    for (gp=0; gp<Q->ngp; gp++) {
      //double normal[3],tangent1[3],tangent2[3],xp,yp,zp;
      double *normal,*tangent1,*tangent2,xp,yp,zp;
      QPntSurfCoefStokes *qpoint = &cell_qpoint[gp];


      QPntSurfCoefStokesGetField_surface_normal(qpoint,&normal);
      QPntSurfCoefStokesGetField_surface_tangent1(qpoint,&tangent1);
      QPntSurfCoefStokesGetField_surface_tangent2(qpoint,&tangent2);

      element->compute_surface_normal_3D(
                                                 element,
                                                 elcoords,    // should contain 27 points with dimension 3 (x,y,z) //
                                                 Q->face_id,   // edge index 0,1,2,3,4,5,6,7 //
                                                 &Q->gp2[gp], // should contain 1 point with dimension 2 (xi,eta)   //
                                                 normal ); // normal[] contains 1 point with dimension 3 (x,y,z) //
      element->compute_surface_tangents_3D(
                                                 element,
                                                 elcoords,    // should contain 27 points with dimension 3 (x,y,z) //
                                                 Q->face_id,
                                                 &Q->gp2[gp], // should contain 1 point with dimension 2 (xi,eta)   //
                                                 tangent1,tangent2 ); // t1[],t2[] contains 1 point with dimension 3 (x,y,z) //

      /* interpolate global coords */
      element->basis_NI_3D(&Q->gp3[gp],Ni);
      xp = yp = zp = 0.0;
      for (k=0; k<element->n_nodes_3D; k++) {
        xp += Ni[k] * elcoords[3*k  ];
        yp += Ni[k] * elcoords[3*k+1];
        zp += Ni[k] * elcoords[3*k+2];
      }

      /*
      printf("[face=%d] fe=%d p=%d (s,t) %1.4e %1.4e: (xi,eta,zeta) %1.4e %1.4e %1.4e: (x,y,z) %1.4e %1.4e %1.4e \n",
                   Q->face_id,fe,gp, Q->gp2[gp].xi,Q->gp2[gp].eta, Q->gp3[gp].xi,Q->gp3[gp].eta,Q->gp3[gp].zeta,
                   xp,yp,zp);
      */
      //printf("%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n",xp,yp,zp,0.1*normal[0],0.1*normal[1],0.1*normal[2]);

    }
  }
  ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SurfaceQuadratureOrientationViewGnuplotStokes(SurfaceQuadrature Q,DM da,const char name[])
{
  PetscErrorCode ierr;
  DM             cda;
  Vec            gcoords;
  PetscScalar    *LA_gcoords;
  PetscInt       nel,nen,fe,e,k,gp;
  const PetscInt *elnidx;
  ConformingElementFamily element;
  double         elcoords[3*Q2_NODES_PER_EL_3D];
  double         Ni[27];
  QPntSurfCoefStokes *all_qpoint;
  QPntSurfCoefStokes *cell_qpoint;
  char fname[PETSC_MAX_PATH_LEN];
  PetscMPIInt rank;
  FILE *file;

  PetscFunctionBegin;


  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  if (name) {
    PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"%s-surfquadrature-face%D-r%D.gp",name,Q->face_id,rank);
  } else {
    PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"surfquadrature-face%D-r%D.gp",Q->face_id,rank);
  }
  file = fopen(fname,"w");
  if (!file) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",fname);
  }

  PetscFPrintf(PETSC_COMM_SELF,file,"# Surface quadrature data (n,t1,t1,traction) for face %D \n",Q->face_id);
  PetscFPrintf(PETSC_COMM_SELF,file,"# nfaces = %D \n",Q->nfaces);

  if (Q->nfaces == 0) { PetscFunctionReturn(0); }

  /* setup for coords */
  ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

  ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
  element = Q->e;

  ierr = SurfaceQuadratureGetAllCellData_Stokes(Q,&all_qpoint);CHKERRQ(ierr);
  for (fe=0; fe<Q->nfaces; fe++) {

    e = Q->element_list[fe];
    ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);

    ierr =  SurfaceQuadratureGetCellData_Stokes(Q,all_qpoint,fe,&cell_qpoint);CHKERRQ(ierr);

    for (gp=0; gp<Q->ngp; gp++) {
      //double normal[3],tangent1[3],tangent2[3],xp,yp,zp;
      double *normal,*tangent1,*tangent2,*traction,xp,yp,zp;
      QPntSurfCoefStokes *qpoint = &cell_qpoint[gp];


      QPntSurfCoefStokesGetField_surface_normal(qpoint,&normal);
      QPntSurfCoefStokesGetField_surface_tangent1(qpoint,&tangent1);
      QPntSurfCoefStokesGetField_surface_tangent2(qpoint,&tangent2);

      QPntSurfCoefStokesGetField_surface_traction(qpoint,&traction);

      element->compute_surface_normal_3D(
                                         element,
                                         elcoords,    // should contain 27 points with dimension 3 (x,y,z) //
                                         Q->face_id,   // edge index 0,1,2,3,4,5,6,7 //
                                         &Q->gp2[gp], // should contain 1 point with dimension 2 (xi,eta)   //
                                         normal ); // normal[] contains 1 point with dimension 3 (x,y,z) //
      element->compute_surface_tangents_3D(
                                           element,
                                           elcoords,    // should contain 27 points with dimension 3 (x,y,z) //
                                           Q->face_id,
                                           &Q->gp2[gp], // should contain 1 point with dimension 2 (xi,eta)   //
                                           tangent1,tangent2 ); // t1[],t2[] contains 1 point with dimension 3 (x,y,z) //

      /* interpolate global coords */
      element->basis_NI_3D(&Q->gp3[gp],Ni);
      xp = yp = zp = 0.0;
      for (k=0; k<element->n_nodes_3D; k++) {
        xp += Ni[k] * elcoords[3*k  ];
        yp += Ni[k] * elcoords[3*k+1];
        zp += Ni[k] * elcoords[3*k+2];
      }

      fprintf(file,"%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e  %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e  %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e  %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n",
              xp,yp,zp,0.1*normal[0],0.1*normal[1],0.1*normal[2],
              xp,yp,zp,0.1*tangent1[0],0.1*tangent1[1],0.1*tangent1[2],
              xp,yp,zp,0.1*tangent2[0],0.1*tangent2[1],0.1*tangent2[2],
              xp,yp,zp,0.1*traction[0],0.1*traction[1],0.1*traction[2]
              );

    }
  }
  ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  fclose(file);

  PetscFunctionReturn(0);
}

PetscErrorCode SurfaceQuadratureGetAllCellData_Stokes(SurfaceQuadrature Q,QPntSurfCoefStokes *coeffs[])
{
  QPntSurfCoefStokes *quadraturepoint_data;
  DataField          PField;
  PetscFunctionBegin;

  if (Q->nfaces) {
    DataBucketGetDataFieldByName(Q->properties_db, QPntSurfCoefStokes_classname ,&PField);
    quadraturepoint_data = PField->data;
    *coeffs = quadraturepoint_data;
  } else {
    *coeffs = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SurfaceQuadratureGetCellData_Stokes(SurfaceQuadrature Q,QPntSurfCoefStokes coeffs[],PetscInt cidx,QPntSurfCoefStokes *cell[])
{
  PetscFunctionBegin;
  *cell = NULL;
  if (cidx >= Q->nfaces) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_SIZ,"cidx > max cells");
  }
  if (Q->nfaces) {
    *cell = &coeffs[cidx*Q->ngp];
  } else {
    *cell = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompCreateSurfaceQuadrature_Stokes(PhysCompStokes ctx)
{
  DM       dav;
  PetscInt face_index;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  dav = ctx->dav;

  ierr = PetscMalloc(sizeof(SurfaceQuadrature)*HEX_EDGES,&ctx->surfQ);CHKERRQ(ierr);
  for (face_index=0; face_index<HEX_EDGES; face_index++) {
    ctx->surfQ[face_index] = NULL;
  }

  for (face_index=0; face_index<HEX_EDGES; face_index++) {
    SurfaceQuadrature surfQ;

    ierr = SurfaceQuadratureCreate_GaussLegendreStokes(dav,face_index,&surfQ);CHKERRQ(ierr);
    ierr = SurfaceQuadratureGeometrySetUpStokes(surfQ,dav);CHKERRQ(ierr);

    // debugging
    //ierr = SurfaceQuadratureOrientationViewGnuplotStokes(surfQ,dav,"init");CHKERRQ(ierr);

    ctx->surfQ[face_index] = surfQ;
  }

  /*
  for (face_index=0; face_index<HEX_EDGES; face_index++) {
    char name[PETSC_MAX_PATH_LEN];

    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"SurfaceGaussLegendre StokesCoefficients[face %D]",face_index);
    DataBucketView(PetscObjectComm((PetscObject)dav), ctx->surfQ[face_index]->properties_db,name,DATABUCKET_VIEW_STDOUT);
  }
  */

  PetscFunctionReturn(0);
}


PetscErrorCode SNESStokes_ConvergenceTest_UPstol(SNES snes,PetscInt it,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason,void *ctx)
{
  Vec X,dX,Xu,Xp,dXu,dXp;
  PetscReal atol,rtol,stol;
  PetscInt maxit,maxf;
  PetscReal xnormUP[2],snormUP[2];
  PetscReal alpha[2];
  pTatinCtx       user;
  PhysCompStokes  stokes;
  DM              stokes_pack;
  PetscErrorCode  ierr;

  PetscFunctionBegin;

  *reason = SNES_CONVERGED_ITERATING;

  user = (pTatinCtx)ctx;
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;

  ierr = SNESGetTolerances(snes,&atol,&rtol,&stol,&maxit,&maxf);CHKERRQ(ierr);

  ierr = SNESGetSolution(snes,&X);CHKERRQ(ierr);
  ierr = SNESGetSolutionUpdate(snes,&dX);CHKERRQ(ierr);

  ierr = DMCompositeGetAccess(stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,dX,&dXu,&dXp);CHKERRQ(ierr);

  ierr = VecNorm(dXu,NORM_2,&snormUP[0]);CHKERRQ(ierr);
  ierr = VecNorm(dXp,NORM_2,&snormUP[1]);CHKERRQ(ierr);

  ierr = VecNorm(Xu,NORM_2,&xnormUP[0]);CHKERRQ(ierr);
  ierr = VecNorm(Xp,NORM_2,&xnormUP[1]);CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(stokes_pack,dX,&dXu,&dXp);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);

  if (it==0) {
    /* set parameter for default relative tolerance convergence test */
        snes->ttol = fnorm*rtol;
  }

  if (fnorm < atol) {
    *reason = SNES_CONVERGED_FNORM_ABS;
        ierr = PetscInfo2(snes,"Converged due to function norm %14.12e < %14.12e\n",(double)fnorm,(double)atol);CHKERRQ(ierr);
  }

  if (it>0 && !*reason) {
        ierr = PetscInfo2(snes,"ConvergenceTest : function norm %14.12e ?<? %14.12e\n",(double)fnorm,(double)atol);CHKERRQ(ierr);

    ierr = PetscInfo2(snes,"ConvergenceTest : small update length (U): %14.12e ?<? %14.12e \n",(double)snormUP[0]/(double)xnormUP[0],(double)stol);CHKERRQ(ierr);
    ierr = PetscInfo2(snes,"ConvergenceTest : small update length (P): %14.12e ?<? %14.12e \n",(double)snormUP[1]/(double)xnormUP[1],(double)stol);CHKERRQ(ierr);

    ierr = PetscInfo2(snes,"ConvergenceTest : function norm %14.12e ?<? %14.12e (relative tolerance)\n",(double)fnorm,(double)snes->ttol);CHKERRQ(ierr);

    // ||dX|| < eps ||X||
    alpha[0] = 1.0;
    alpha[1] = 0.0;
    if ( snormUP[0] < alpha[0] * stol * xnormUP[0] ) {
      *reason = SNES_CONVERGED_SNORM_RELATIVE;
      ierr = PetscInfo3(snes,"Converged due to small update length (U): %14.12e < %14.12e * %14.12e\n",(double)snormUP[0],(double)stol,(double)xnormUP[0]);CHKERRQ(ierr);
    }
    if ( snormUP[1] < alpha[1] * stol * xnormUP[1] ) {
      *reason = SNES_CONVERGED_SNORM_RELATIVE;
      ierr = PetscInfo3(snes,"Converged due to small update length (P): %14.12e < %14.12e * %14.12e\n",(double)snormUP[1],(double)stol,(double)xnormUP[1]);CHKERRQ(ierr);
    }

    if (fnorm <= snes->ttol) {
      *reason = SNES_CONVERGED_FNORM_RELATIVE;
      ierr = PetscInfo2(snes,"Converged due to function norm %14.12e < %14.12e (relative tolerance)\n",(double)fnorm,(double)snes->ttol);CHKERRQ(ierr);
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode SNESStokes_SetConvergenceTest_UPstol(SNES snes,pTatinCtx user)
{
  const char *prefix;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Activating \"ConvergenceTest_UPstol\" on SNES (%s)\n",prefix);

  //ierr = SNESSetApplicationContext(snes,(void*)user);CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(snes,SNESStokes_ConvergenceTest_UPstol,(void**)user,NULL);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode PetscOptionsInsertPrefixString(const char prefix[],const char option[])
{
  PetscErrorCode ierr;
  char opt[PETSC_MAX_PATH_LEN];

  if (option[0] != '-') {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Option must start with a -");
  }

  if (!prefix) {
    ierr = PetscOptionsInsertString(NULL,option);CHKERRQ(ierr);
  } else {
    sprintf(opt,"-%s%s",prefix,&option[1]);
    ierr = PetscOptionsInsertString(NULL,opt);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SNESStokesPCSetOptions_A(SNES snes)
{
  const char *prefix;

  PetscErrorCode ierr;

  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-ksp_type fgmres");CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-pc_type fieldsplit");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-pc_fieldsplit_type schur");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-pc_fieldsplit_schur_factorization_type upper");CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_p_ksp_type preonly");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_p_pc_type jacobi");CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_ksp_type fgmres");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_ksp_rtol 1.0e-2");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_ksp_max_it 40");CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_ksp_type chebychev");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_ksp_norm_type NONE");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_ksp_max_it 6");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_pc_type jacobi");CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_esteig_ksp_norm_type NONE");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_ksp_chebychev_estimate_eigenvalues 0,0.2,0,1.1");CHKERRQ(ierr);

  //ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_type gmres");CHKERRQ(ierr);
  //ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_pc_type bjacobi");CHKERRQ(ierr);
  //ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_pc_type lu");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SNESStokesPCMGSetOptions(SNES snes,PetscInt maxits,PetscBool mglog)
{
  const char *prefix;
  char opt[PETSC_MAX_PATH_LEN];

  PetscErrorCode ierr;

  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_pc_type mg");CHKERRQ(ierr);

  {
    KSP ksp,*split_ksp;
    PC pc,split_pc;
    PetscInt nsplits,nlevels;

    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&split_ksp);CHKERRQ(ierr);

    ierr = KSPGetPC(split_ksp[0],&split_pc);CHKERRQ(ierr);
    ierr = PCMGGetLevels(split_pc,&nlevels);CHKERRQ(ierr);

    PetscSNPrintf(opt,PETSC_MAX_PATH_LEN-1,"-fieldsplit_u_pc_mg_levels %D",nlevels);
    ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);
  }

  if (mglog) {
    ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_pc_mg_log");CHKERRQ(ierr);
    ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_converged_reason");CHKERRQ(ierr);
  }

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_ksp_type chebychev");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_ksp_norm_type NONE");CHKERRQ(ierr);

  PetscSNPrintf(opt,PETSC_MAX_PATH_LEN-1,"-fieldsplit_u_mg_levels_ksp_max_it %D",maxits);
  ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_esteig_ksp_norm_type NONE");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_ksp_chebychev_estimate_eigenvalues 0,0.2,0,1.1");CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_levels_pc_type jacobi");CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

PetscErrorCode SNESStokesPCMGCoarseSetOptions_IterativeASM(SNES snes,PetscReal rtol,PetscInt maxits,PetscInt overlap)
{
  const char *prefix;
  char opt[PETSC_MAX_PATH_LEN];

  PetscErrorCode ierr;

  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_type gmres");CHKERRQ(ierr);

  PetscSNPrintf(opt,PETSC_MAX_PATH_LEN-1,"-fieldsplit_u_mg_coarse_ksp_max_it %D",maxits);
  ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);

  PetscSNPrintf(opt,PETSC_MAX_PATH_LEN-1,"-fieldsplit_u_mg_coarse_ksp_rtol %1.4e",rtol);
  ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_pc_type asm");CHKERRQ(ierr);

  PetscSNPrintf(opt,PETSC_MAX_PATH_LEN-1,"-fieldsplit_u_mg_coarse_pc_asm_overlap %D",overlap);
  ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_sub_pc_type ilu");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SNESStokesPCMGCoarseSetOptions_NestedIterativeASM(SNES snes,PetscReal rtol,PetscInt maxits,PetscInt maxitsnested,PetscInt overlap)
{
  const char *prefix;
  char opt[PETSC_MAX_PATH_LEN];

  PetscErrorCode ierr;

  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_type fgmres");CHKERRQ(ierr);

  PetscSNPrintf(opt,PETSC_MAX_PATH_LEN-1,"-fieldsplit_u_mg_coarse_ksp_max_it %D",maxits);
  ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);

  sprintf(opt,"-fieldsplit_u_mg_coarse_ksp_rtol %1.4e",rtol);
  ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);

  /* defined nested ksp solver on coarse grid */
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_pc_type ksp");CHKERRQ(ierr);

  PetscSNPrintf(opt,PETSC_MAX_PATH_LEN-1,"-fieldsplit_u_mg_coarse_ksp_ksp_max_it %D",maxitsnested);
  ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_ksp_type chebychev");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_ksp_norm_type NONE");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_esteig_ksp_norm_type NONE");CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_ksp_chebychev_estimate_eigenvalues 0,0.2,0,1.1");CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_pc_type asm");CHKERRQ(ierr);
  PetscSNPrintf(opt,PETSC_MAX_PATH_LEN-1,"-fieldsplit_u_coarse_ksp_pc_asm_overlap %D",overlap);
  ierr = PetscOptionsInsertPrefixString(prefix,opt);CHKERRQ(ierr);
  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_ksp_sub_pc_type ilu");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SNESStokesPCMGCoarseSetOptions_SparseDirect(SNES snes)
{
  const char *prefix;
  MPI_Comm comm;
  PetscErrorCode ierr;
  PetscMPIInt nproc;

  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);

  ierr = PetscObjectGetComm((PetscObject)snes,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);

  ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_pc_type lu");CHKERRQ(ierr);

  if (nproc == 1) {
    ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_pc_factor_mat_solver_package umfpack");CHKERRQ(ierr);
  } else {
    ierr = PetscOptionsInsertPrefixString(prefix,"-fieldsplit_u_mg_coarse_pc_factor_mat_solver_package superlu_dist");CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesGetDMComposite(PhysCompStokes stokes,DM *dmc)
{
    if (dmc) { *dmc = stokes->stokes_pack; }
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesGetDMs(PhysCompStokes stokes,DM *dmv,DM *dmp)
{
    if (dmv) { *dmv = stokes->dav; }
    if (dmp) { *dmp = stokes->dap; }
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesGetBCList(PhysCompStokes stokes,BCList *ulist,BCList *plist)
{
    if (ulist) { *ulist = stokes->u_bclist; }
    if (plist) { *plist = stokes->p_bclist; }
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesGetVolumeQuadrature(PhysCompStokes stokes,Quadrature *q)
{
    if (q) { *q = stokes->volQ; }
  PetscFunctionReturn(0);
}
PetscErrorCode PhysCompStokesGetVolumeQuadratureAllCellData(PhysCompStokes stokes,QPntVolCoefStokes *coeffs[])
{
    PetscErrorCode ierr;
    if (coeffs) {
        ierr = VolumeQuadratureGetAllCellData_Stokes(stokes->volQ,coeffs);CHKERRQ(ierr);
    }
  PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesGetSurfaceQuadrature(PhysCompStokes stokes,HexElementFace fid,SurfaceQuadrature *sq)
{
    if (sq) { *sq = stokes->surfQ[fid]; };
    PetscFunctionReturn(0);
}

PetscErrorCode PhysCompStokesGetSurfaceQuadratureAllCellData(PhysCompStokes stokes,HexElementFace fid,QPntSurfCoefStokes *coeffs[])
{
  PetscErrorCode ierr;
    if (coeffs) {
        ierr = SurfaceQuadratureGetAllCellData_Stokes(stokes->surfQ[fid],coeffs);CHKERRQ(ierr);
    }
  PetscFunctionReturn(0);
}

PetscErrorCode Stokes_KSPConvergenceTest_ScaledResiduals(KSP ksp,PetscInt it,PetscReal rnorm,KSPConvergedReason *reason,void *data)
{
  PetscErrorCode  ierr;
  pTatinCtx       ctx;
    const char      *names[] = { "U", "V", "W", "P" };
    PetscInt        s,component_cnt;
  PetscReal       norms[4],max[4],maxX[4],minX[4],fabs_minX[4];
  Vec             X,Xu,Xp,F,Fu,Fp,v,w;
  Mat             A;
  PetscReal       atol,rtol;
    static PetscInt initial_norms[4];
    PetscBool       relative_component_norm_rtol;

  PetscFunctionBegin;

    if (it < 2) {
        PetscPrintf(PETSC_COMM_WORLD,"  *** Stokes_KSPConvergenceTest_ScaledResiduals: always do two iterations ***\n");
        *reason = KSP_CONVERGED_ITERATING;
        PetscFunctionReturn(0);
    }


  ctx = (pTatinCtx)data;

    ierr = KSPGetTolerances(ksp,&rtol,&atol,NULL,NULL);CHKERRQ(ierr);

  ierr = KSPGetOperators(ksp,&A,0);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&w,&v);CHKERRQ(ierr);

  ierr = KSPBuildResidual(ksp,v,w,&F);CHKERRQ(ierr);
    ierr = VecAbs(F);CHKERRQ(ierr);
    ierr = VecChop(F,1.0e-12);CHKERRQ(ierr);

  ierr = DMCompositeGetAccess(ctx->stokes_ctx->stokes_pack,F,&Fu,&Fp);CHKERRQ(ierr);

  ierr = VecStrideNorm(Fu,0,NORM_2,&norms[0]);CHKERRQ(ierr);
  ierr = VecStrideNorm(Fu,1,NORM_2,&norms[1]);CHKERRQ(ierr);
  ierr = VecStrideNorm(Fu,2,NORM_2,&norms[2]);CHKERRQ(ierr);
  ierr = VecNorm(Fp,NORM_2,&norms[3]);CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"  *** Stokes_KSPConvergenceTest_ScaledResiduals:        ResidualNorms(%s,%s,%s,%s) = [ %1.4e %1.4e %1.4e %1.4e ] ***\n",
                names[0],names[1],names[2],names[3],norms[0],norms[1],norms[2],norms[3]);

    if (it == 2) {
        initial_norms[0] = norms[0];
        initial_norms[1] = norms[1];
        initial_norms[2] = norms[2];
        initial_norms[3] = norms[3];
        PetscPrintf(PETSC_COMM_WORLD,"  *** Stokes_KSPConvergenceTest_ScaledResiduals: InitialResidualNorms(%s,%s,%s,%s) = [ %1.4e %1.4e %1.4e %1.4e ] ***\n",
                    names[0],names[1],names[2],names[3],initial_norms[0],initial_norms[1],initial_norms[2],initial_norms[3]);
    }

  ierr = VecStrideMax(Fu,0,NULL,&max[0]);CHKERRQ(ierr);
  ierr = VecStrideMax(Fu,1,NULL,&max[1]);CHKERRQ(ierr);
  ierr = VecStrideMax(Fu,2,NULL,&max[2]);CHKERRQ(ierr);
  ierr = VecMax(Fp,NULL,&max[3]);CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"  *** Stokes_KSPConvergenceTest_ScaledResiduals:              max(r%s,r%s,r%s,r%s) = [ %1.4e %1.4e %1.4e %1.4e ] ***\n",
                names[0],names[1],names[2],names[3],max[0],max[1],max[2],max[3]);


    *reason = KSP_CONVERGED_ITERATING;

    relative_component_norm_rtol = PETSC_TRUE;
    for (s=0; s<4; s++) {
        if (max[s] > atol) {
            relative_component_norm_rtol = PETSC_FALSE;
        }
    }

    /* residuals have been reduced sufficiently, compare with the solution we have */
  ierr = KSPBuildSolution(ksp,v,&X);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);

  /* aggresive point wise checking */
  /*
  {
    PetscScalar *LA_Fi,*LA_Xi;
    PetscInt fails = 0,fail_sum,i,m;

    ierr = VecGetLocalSize(F,&m);CHKERRQ(ierr);
    ierr = VecGetArray(F,&LA_Fi);CHKERRQ(ierr);
    ierr = VecGetArray(X,&LA_Xi);CHKERRQ(ierr);

    for (i=0; i<m; i++) {
      if (PetscAbsReal(PetscRealPart(LA_Xi[i])) > 1.0e-16) {
        if (PetscAbsReal(PetscRealPart(LA_Fi[i])) > rtol * PetscAbsReal(PetscRealPart(LA_Xi[i]))) {
          fails++;
        }
      }
    }
    ierr = VecRestoreArray(X,&LA_Xi);CHKERRQ(ierr);
    ierr = VecRestoreArray(F,&LA_Fi);CHKERRQ(ierr);

    ierr = MPI_Allreduce(&fails,&fail_sum,1,MPIU_INT,MPI_SUM,PetscObjectComm((PetscObject)ksp));CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"  *** Stokes_KSPConvergenceTest_ScaledResiduals: failed sum %D ***\n",fail_sum);

    *reason = KSP_CONVERGED_ITERATING;
    if (fail_sum == 0) {
      *reason = KSP_CONVERGED_RTOL;
    }
  }
  */

    /* more relaxed vector-wise hecking */
    ierr = VecStrideMin(Xu,0,NULL,&minX[0]);CHKERRQ(ierr);
    ierr = VecStrideMin(Xu,1,NULL,&minX[1]);CHKERRQ(ierr);
    ierr = VecStrideMin(Xu,2,NULL,&minX[2]);CHKERRQ(ierr);
    ierr = VecMin(Xp,NULL,&minX[3]);CHKERRQ(ierr);

    ierr = VecStrideMax(Xu,0,NULL,&maxX[0]);CHKERRQ(ierr);
    ierr = VecStrideMax(Xu,1,NULL,&maxX[1]);CHKERRQ(ierr);
    ierr = VecStrideMax(Xu,2,NULL,&maxX[2]);CHKERRQ(ierr);
    ierr = VecMax(Xp,NULL,&maxX[3]);CHKERRQ(ierr);

    for (s=0; s<4; s++) {
        fabs_minX[s] = PetscMin(PetscAbsReal(minX[s]),PetscAbsReal(maxX[s]));
    }
    PetscPrintf(PETSC_COMM_WORLD,"  *** Stokes_KSPConvergenceTest_ScaledResiduals:             fmin(X%s,X%s,X%s,X%s) = [ %1.4e %1.4e %1.4e %1.4e ] ***\n",
                names[0],names[1],names[2],names[3],fabs_minX[0],fabs_minX[1],fabs_minX[2],fabs_minX[3]);

    if (relative_component_norm_rtol) {
        component_cnt = 0;
        for (s=0; s<4; s++) {
            if (max[s] < rtol * fabs_minX[s]) {
                PetscPrintf(PETSC_COMM_WORLD,"  *** Stokes_KSPConvergenceTest_ScaledResiduals:                    ++ max(r%s) < %1.1e min(X%s) ***\n",names[s],rtol,names[s]);
                component_cnt++;
            } else {
                PetscPrintf(PETSC_COMM_WORLD,"  *** Stokes_KSPConvergenceTest_ScaledResiduals:                    -- max(r%s) > %1.1e min(X%s) ***\n",names[s],rtol,names[s]);
            }
        }
        // all residuals pass //
        *reason = KSP_CONVERGED_ITERATING;
        if (component_cnt == 4) {
            *reason = KSP_CONVERGED_RTOL;
        }
    }

  ierr = DMCompositeRestoreAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(ctx->stokes_ctx->stokes_pack,F,&Fu,&Fp);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SNESStokes_KSPSetConvergenceTest_ScaledResiduals(SNES snes,pTatinCtx user)
{
  KSP            ksp;
  const char     *prefix;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Activating \"KSPSetConvergenceTest_ScaledResiduals\" on SNES->KSP (%s)\n",prefix);

  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPSetConvergenceTest(ksp,Stokes_KSPConvergenceTest_ScaledResiduals,(void**)user,NULL);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

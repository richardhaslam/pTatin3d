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
 **    filename:   analytics_vv.c
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


#include <petsc.h>
#include <petscvec.h>
#include <ptatin3d.h>
#include <private/ptatin_impl.h>
#include <dmda_bcs.h>
#include <data_bucket.h>
#include <MPntStd_def.h>
#include <MPntPStokes_def.h>
#include <ptatin_std_dirichlet_boundary_conditions.h>
#include <dmda_element_q2p1.h>
#include <material_point_std_utils.h>
#include <output_material_points.h>
#include <pswarm.h>

#include <inorms.h>
#include <models/analytics_vv/SolKxSolution.h>
#include <models/analytics_vv/SolCxSolution.h>

PSwarm pswarm;

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_AnlVV"
PetscErrorCode ModelInitialize_AnlVV(pTatinCtx c,void *ctx)
{
  RheologyConstants *rheology;
  PetscErrorCode    ierr;
  
  PetscFunctionBegin;
  ierr = pTatinGetRheology(c,&rheology);CHKERRQ(ierr);
  rheology->rheology_type = RHEOLOGY_VISCOUS;
  
  ierr = PSwarmCreate(PETSC_COMM_WORLD,&pswarm);CHKERRQ(ierr);
  ierr = PSwarmSetOptionsPrefix(pswarm,"passive_");CHKERRQ(ierr);
  ierr = PSwarmSetPtatinCtx(pswarm,c);CHKERRQ(ierr);
  ierr = PSwarmSetTransportModeType(pswarm,PSWARM_TM_EULERIAN);CHKERRQ(ierr);
  ierr = PSwarmSetFromOptions(pswarm);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_AnlVV"
PetscErrorCode ModelApplyBoundaryCondition_AnlVV(pTatinCtx c,void *ctx)
{
  PetscScalar    zero = 0.0;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_AnlVV"
PetscErrorCode ModelApplyBoundaryConditionMG_AnlVV(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
  PetscScalar    zero = 0.0;
  PetscInt       n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  for (n=0; n<nl; n++) {
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_AnlVV"
PetscErrorCode ModelApplyInitialMeshGeometry_AnlVV(pTatinCtx c,void *ctx)
{
  PetscReal      Lx,Ly,Lz;
  PetscReal      gvec[] = { 0.0, -1.0, 0.0 };
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PhysCompStokesSetGravityVector(c->stokes_ctx,gvec);CHKERRQ(ierr);
  Lx = 1.0;
  Ly = 1.0;
  Lz = 1.0;
  ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,0.0,Lx, 0.0,Ly, 0.0,Lz);CHKERRQ(ierr);
  ierr = PSwarmSetUp(pswarm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_AnlVV_solkx"
PetscErrorCode ModelApplyInitialMaterialGeometry_AnlVV_solkx(pTatinCtx c,void *ctx)
{
  int            p,n_mp_points;
  DataBucket     db;
  DataField      PField_std,PField_stokes;
  PetscReal      B,m,km,kn;
  PetscInt       n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  /* define properties on material points */
  ierr = pTatinGetMaterialPoints(c,&db,NULL);CHKERRQ(ierr);
  DataBucketGetSizes(db,&n_mp_points,0,0);
  
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  
  DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
  DataFieldGetAccess(PField_stokes);

  m = 1.2;
  n = 2;
  B = 1.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-solkx_m",&m,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-solkx_B",&B,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-solkx_n",&n,NULL);CHKERRQ(ierr);
  
  {
    PetscReal pos[2],vel[2],pressure;
    
    pos[0] = 0.2;
    pos[1] = 0.2;
    SolKxSolution(pos,m,n,B,vel,&pressure, NULL,NULL,NULL);
    printf("%+1.12e %+1.12e vx %+1.12e\n",pos[0],pos[1],vel[0]);
    pos[0] = 0.2;
    pos[1] = 0.7;
    SolKxSolution(pos,m,n,B,vel,&pressure, NULL,NULL,NULL);
    printf("%+1.12e %+1.12e vx %+1.12e\n",pos[0],pos[1],vel[0]);

    pos[0] = 0.0;
    pos[1] = 0.5;
    SolKxSolution(pos,m,n,B,vel,&pressure, NULL,NULL,NULL);
    printf("%+1.12e %+1.12e vy %+1.12e\n",pos[0],pos[1],vel[1]);
    pos[0] = 0.5;
    pos[1] = 0.5;
    SolKxSolution(pos,m,n,B,vel,&pressure, NULL,NULL,NULL);
    printf("%+1.12e %+1.12e vy %+1.12e\n",pos[0],pos[1],vel[1]);

    pos[0] = 0.25;
    pos[1] = 1.0;
    SolKxSolution(pos,m,n,B,vel,&pressure, NULL,NULL,NULL);
    printf("%+1.12e %+1.12e vy %+1.12e <top>\n",pos[0],pos[1],vel[1]);
    pos[0] = 0.75;
    pos[1] = 1.0;
    SolKxSolution(pos,m,n,B,vel,&pressure, NULL,NULL,NULL);
    printf("%+1.12e %+1.12e vy %+1.12e <top>\n",pos[0],pos[1],vel[1]);
    
    
    pos[0] = 0.1;
    pos[1] = 0.1;
    SolKxSolution(pos,m,n,B,vel,&pressure, NULL,NULL,NULL);
    printf("%+1.12e %+1.12e p %+1.12e\n",pos[0],pos[1],pressure);
    pos[0] = 0.1;
    pos[1] = 0.8;
    SolKxSolution(pos,m,n,B,vel,&pressure, NULL,NULL,NULL);
    printf("%+1.12e %+1.12e p %+1.12e\n",pos[0],pos[1],pressure);
  }
  
  km = m * PETSC_PI; /* solution valid for km not zero -- should get trivial solution if km=0 */
  kn = (PetscReal)n * PETSC_PI;
  
  for (p=0; p<n_mp_points; p++) {
    MPntStd     *material_point;
    MPntPStokes *material_point_properties_stokes;
    double      *position;
    double      eta,rho;
    int         phase;
    
    DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
    DataFieldAccessPoint(PField_stokes,p,(void**)&material_point_properties_stokes);
    
    MPntStdGetField_global_coord(material_point,&position);
    
    phase = 0;
    eta = exp(2.0*B*position[0]);
    rho = -PetscSinReal(km*position[1])*PetscCosReal(kn*position[0]);
    
    MPntStdSetField_phase_index(material_point,phase);
    MPntPStokesSetField_eta_effective(material_point_properties_stokes,eta);
    MPntPStokesSetField_density(material_point_properties_stokes,rho);
  }
  
  DataFieldRestoreAccess(PField_std);
  DataFieldRestoreAccess(PField_stokes);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_AnlVV_solcx"
PetscErrorCode ModelApplyInitialMaterialGeometry_AnlVV_solcx(pTatinCtx c,void *ctx)
{
  int            p,n_mp_points;
  DataBucket     db;
  DataField      PField_std,PField_stokes;
  PetscReal      eta0,eta1,xc,m,km,kn;
  PetscInt       n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  /* define properties on material points */
  ierr = pTatinGetMaterialPoints(c,&db,NULL);CHKERRQ(ierr);
  DataBucketGetSizes(db,&n_mp_points,0,0);
  
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  
  DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
  DataFieldGetAccess(PField_stokes);
  
  eta0 = 1.0;
  eta1 = 1.0;
  xc = 0.5;
  n = 2;
  m = 1.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-solcx_eta0",&eta0,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-solcx_eta1",&eta1,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-solcx_xc",&xc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-solcx_n",&n,NULL);CHKERRQ(ierr);
  
  km = m * PETSC_PI; /* solution valid for km not zero -- should get trivial solution if km=0 */
  kn = (PetscReal)n * PETSC_PI;
  
  for (p=0; p<n_mp_points; p++) {
    MPntStd     *material_point;
    MPntPStokes *material_point_properties_stokes;
    double      *position;
    double      eta,rho;
    int         phase;
    
    DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
    DataFieldAccessPoint(PField_stokes,p,(void**)&material_point_properties_stokes);
    
    MPntStdGetField_global_coord(material_point,&position);
    
    phase = 0;
    if (position[0] <= xc) {
      eta = eta0;
    } else {
      eta = eta1;
    }
    rho = -PetscCosReal(km*position[0])*PetscSinReal(kn*position[1]);
    
    MPntStdSetField_phase_index(material_point,phase);
    MPntPStokesSetField_eta_effective(material_point_properties_stokes,eta);
    MPntPStokesSetField_density(material_point_properties_stokes,rho);
  }
  
  DataFieldRestoreAccess(PField_std);
  DataFieldRestoreAccess(PField_stokes);
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelOutput_AnlVV"
PetscErrorCode ModelOutput_AnlVV(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
  PetscErrorCode ierr;
  PhysCompStokes stokes_ctx;
  DM dms,dmv,dmp;
  Vec Xv,Xp;
  Quadrature volQ;
  PetscReal ds,ev_L2,eE_H1s,eE_H1,ep_L2;
  
  PetscFunctionBegin;
  ierr = PSwarmFieldUpdateAll(pswarm);CHKERRQ(ierr);
  /*ierr = PSwarmViewInfo(pswarm);CHKERRQ(ierr);*/
  
  /* write paraview files - non-essential */
  ierr = PSwarmView(pswarm,PSW_VT_SINGLETON);CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);

  ierr = pTatinGetStokesContext(c,&stokes_ctx);CHKERRQ(ierr);
  ierr = PhysCompStokesGetVolumeQuadrature(stokes_ctx,&volQ);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes_ctx,&dms);CHKERRQ(ierr);
  ierr = DMCompositeGetEntries(dms,&dmv,&dmp);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(dms,X,&Xv,&Xp);CHKERRQ(ierr);


  /*
  {
    PetscReal m,B;
    PetscInt n;
    ParamsSolKx solkx;

    m = 1.2;
    n = 2;
    B = 1.0;
    ierr = PetscOptionsGetReal(NULL,NULL,"-solkx_m",&m,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-solkx_B",&B,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-solkx_n",&n,NULL);CHKERRQ(ierr);
    ierr = SolKxParamsSetValues(&solkx,m,n,B);CHKERRQ(ierr);
    
    ierr = Evaluate_uL2_symuH1_pL2(dmv,dmp,Xv,Xp,volQ,
                                   EvaluateV_SolKx,EvaluateE_SolKx,EvaluateP_SolKx,(void*)&solkx,
                                   &ev_L2,&eE_H1s,&eE_H1,&ep_L2);CHKERRQ(ierr);
  }
  */

  {
    PetscReal eta0,eta1,xc;
    PetscInt n;
    ParamsSolCx solcx;
    
    eta0 = 1.0;
    eta1 = 1.0;
    xc = 0.5;
    n = 2;
    ierr = PetscOptionsGetReal(NULL,NULL,"-solcx_eta0",&eta0,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-solcx_eta1",&eta1,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-solcx_xc",&xc,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-solcx_n",&n,NULL);CHKERRQ(ierr);
    ierr = SolCxParamsSetValues(&solcx,eta0,eta1,xc,n);CHKERRQ(ierr);
    
    ierr = Evaluate_uL2_symuH1_pL2(dmv,dmp,Xv,Xp,volQ,
                                   EvaluateV_SolCx,EvaluateE_SolCx,EvaluateP_SolCx,(void*)&solcx,
                                   &ev_L2,&eE_H1s,&eE_H1,&ep_L2);CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccess(dms,X,&Xv,&Xp);CHKERRQ(ierr);

  ds = 1.0 / (PetscReal)(c->mx);         /* ds     v_L2   H1_semi     H1      p_L2 */
  PetscPrintf(PETSC_COMM_WORLD,"[inorms] %+1.12e,%+1.12e,%+1.12e,%+1.12e,%+1.12e\n",ds,ev_L2,eE_H1s,eE_H1,ep_L2);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialCondition_AnlVV"
PetscErrorCode ModelInitialCondition_AnlVV(pTatinCtx c,Vec X,void *ctx)
{
  PhysCompStokes stokes_ctx;
  DM             stokes_pack,dau,dap;
  Vec            velocity,pressure;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = pTatinGetStokesContext(c,&stokes_ctx);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes_ctx,&stokes_pack);CHKERRQ(ierr);
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
  ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  
  ierr = PSwarmAttachStateVecVelocityPressure(pswarm,X);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_AnlVV"
PetscErrorCode ModelDestroy_AnlVV(pTatinCtx c,void *ctx)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[ModelDestroy_AnlVV]]\n");
  ierr = PSwarmDestroy(&pswarm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelCreate_AnlVV"
PetscErrorCode pTatinModelCreate_AnlVV(pTatinModel m)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  /* Allocate memory for the data structure for this model */
  
  /* Set model data */
  ierr = pTatinModelSetUserData(m,NULL);CHKERRQ(ierr);
  
  /* Set function pointers */
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_AnlVV);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelInitialCondition_AnlVV);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_AnlVV);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_AnlVV);CHKERRQ(ierr);
  /*ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_AnlVV);CHKERRQ(ierr);*/
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_AnlVV);CHKERRQ(ierr);
  
  //ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_AnlVV_solkx);CHKERRQ(ierr);
  
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_AnlVV_solcx);CHKERRQ(ierr);
  
  /*ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_AnlVV);CHKERRQ(ierr);*/
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_AnlVV);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_AnlVV);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

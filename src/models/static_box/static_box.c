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
 **    filename:   static_box.c
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


PSwarm pswarm;

PetscErrorCode ModelInitialize_StaticBox(pTatinCtx c,void *ctx)
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

PetscErrorCode ModelApplyBoundaryCondition_StaticBox(pTatinCtx c,void *ctx)
{
  PetscScalar    zero = 0.0;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyBoundaryConditionMG_StaticBox(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
  PetscScalar    zero = 0.0;
  PetscInt       n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  for (n=0; n<nl; n++) {
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
    ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMeshGeometry_StaticBox(pTatinCtx c,void *ctx)
{
  PetscReal      Lx,Ly,Lz;
  PetscReal      gvec[] = { 0.0, -10.0, 0.0 };
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PhysCompStokesSetGravityVector(c->stokes_ctx,gvec);CHKERRQ(ierr);
  Lx = 6.0;
  Ly = 6.0;
  Lz = 6.0;
  ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,0.0,Lx, 0.0,Ly, 0.0,Lz);CHKERRQ(ierr);
  ierr = PSwarmSetUp(pswarm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_StaticBox(pTatinCtx c,void *ctx)
{
  int            p,n_mp_points;
  DataBucket     db;
  DataField      PField_std,PField_stokes;
  PetscBool      variable_density = PETSC_FALSE;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscOptionsGetBool(NULL,NULL,"-variable_density",&variable_density,NULL);CHKERRQ(ierr);
  
  /* define properties on material points */
  ierr = pTatinGetMaterialPoints(c,&db,NULL);CHKERRQ(ierr);
  DataBucketGetSizes(db,&n_mp_points,0,0);
  
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  
  DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
  DataFieldGetAccess(PField_stokes);
  
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
    eta = 10.0;
    rho = 2.0;
    if (variable_density) {
      if (position[1] <= 3.0) {
        phase = 1;
        eta = 10.0;
        rho = 5.0;
      }
    }
    
    MPntStdSetField_phase_index(material_point,phase);
    MPntPStokesSetField_eta_effective(material_point_properties_stokes,eta);
    MPntPStokesSetField_density(material_point_properties_stokes,rho);
  }
  
  DataFieldRestoreAccess(PField_std);
  DataFieldRestoreAccess(PField_stokes);
  
  PetscFunctionReturn(0);
}

/*
 depth  = 0
 rho0 = 2.0
 depth_interface = 3.0
 rho1 = 5.0
 */
static PetscErrorCode ComputeHydrostaticPressure(PetscReal depth,PetscReal rho0,PetscReal rho1,PetscReal depth_interface,PetscReal *pressure)
{
  PetscReal p0;
  
  PetscFunctionBegin;
  p0 = 0.0;
  if (depth <= depth_interface) {
    *pressure = rho0 * 10.0 * depth;
  } else {
    p0 = rho0 * 10.0 * depth_interface;
    
    *pressure = p0 + rho1 * 10.0 * (depth - depth_interface);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ModelOutput_StaticBox(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
  PetscReal      nrm[3];
  Vec            Xu,Xp;
  PhysCompStokes stokes_ctx;
  DM             stokes;
  PetscInt       k;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PSwarmFieldUpdateAll(pswarm);CHKERRQ(ierr);
  /*ierr = PSwarmViewInfo(pswarm);CHKERRQ(ierr);*/
  
  /* write paraview files - non-essential */
  /*
  ierr = PSwarmView(pswarm,PSW_VT_SINGLETON);CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
  */

  if (c->step != c->nsteps) { PetscFunctionReturn(0); }

  /* examine max magnitude of vx,vy,vz */
  ierr = pTatinGetStokesContext(c,&stokes_ctx);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes_ctx,&stokes);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes,X,&Xu,&Xp);CHKERRQ(ierr);
  ierr = VecStrideNorm(Xu,0,NORM_INFINITY,&nrm[0]);CHKERRQ(ierr);
  ierr = VecStrideNorm(Xu,1,NORM_INFINITY,&nrm[1]);CHKERRQ(ierr);
  ierr = VecStrideNorm(Xu,1,NORM_INFINITY,&nrm[2]);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(stokes,X,&Xu,&Xp);CHKERRQ(ierr);
  for (k=0; k<3; k++) {
    if (nrm[k] > 1.0e-12) {
      PetscPrintf(PETSC_COMM_WORLD,"[staticBox] ||v%D - v%D_exact||_inf = %+1.12e <fail>\n",k,k,nrm[k]);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,"[staticBox] ||v%D - v%D_exact||_inf <= 1.0e-12 <pass>\n",k,k);
    }
  }
  
  /* Compare computed pressure with rho.g.z */
  {
    DataBucket pdb;
    int        p,npoints;
    DataField  datafield_tracers,datafield;
    MPntStd    *tracer;
    double     *tracer_pressure;
    PetscReal  rho0,rho1,depth_interface;
    PetscBool  variable_density = PETSC_FALSE;
    PetscReal  max_abs_delta_p,max_abs_delta_p_g;
    
    ierr = PSwarmGetDataBucket(pswarm,&pdb);CHKERRQ(ierr);
    DataBucketGetSizes(pdb,&npoints,0,0);
    
    DataBucketGetDataFieldByName(pdb,MPntStd_classname,&datafield_tracers); DataFieldGetEntries(datafield_tracers,(void**)&tracer);
    DataBucketGetDataFieldByName(pdb,"pressure",&datafield);       DataFieldGetEntries(datafield,(void**)&tracer_pressure);
    
    max_abs_delta_p = 0.0;
    rho0 = 2.0;
    rho1 = 5.0;
    depth_interface = 1.0e32;
    ierr = PetscOptionsGetBool(NULL,NULL,"-variable_density",&variable_density,NULL);CHKERRQ(ierr);
    if (variable_density) {
      depth_interface = 3.0;
    }
    
    for (p=0; p<npoints; p++) {
      PetscReal depth,p_hydrostatic,delta;
      
      depth = 6.0 - (PetscReal)tracer[p].coor[1];
      ierr = ComputeHydrostaticPressure(depth,rho0,rho1,depth_interface,&p_hydrostatic);CHKERRQ(ierr);
      
      delta = PetscAbsReal(tracer_pressure[p] - p_hydrostatic);
      /*PetscPrintf(PETSC_COMM_SELF,"delta %+1.12e : numeric %+1.12e : exact %+1.12e\n",delta,tracer_pressure[p],p_hydrostatic);*/
      if (delta > max_abs_delta_p) {
        max_abs_delta_p = delta;
      }
    }
    
    DataFieldRestoreEntries(datafield_tracers,(void**)&tracer);
    DataFieldRestoreEntries(datafield,(void**)&tracer_pressure);
    
    ierr = MPI_Allreduce(&max_abs_delta_p,&max_abs_delta_p_g,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
    if (max_abs_delta_p_g > 1.0e-12) {
      PetscPrintf(PETSC_COMM_WORLD,"[staticBox] ||p - p_exact||_inf = %+1.12e <fail>\n",max_abs_delta_p_g);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,"[staticBox] ||p - p_exact||_inf <= 1.0e-12 <pass>\n");
    }
    if (!variable_density) PetscPrintf(PETSC_COMM_WORLD,"[staticBox] rho = const.\n");
    else PetscPrintf(PETSC_COMM_WORLD,"[staticBox] rho = variable\n");
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode ModelInitialCondition_StaticBox(pTatinCtx c,Vec X,void *ctx)
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

PetscErrorCode ModelDestroy_StaticBox(pTatinCtx c,void *ctx)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[ModelDestroy_StaticBox]]\n");
  ierr = PSwarmDestroy(&pswarm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinModelCreate_StaticBox(pTatinModel m)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  /* Allocate memory for the data structure for this model */
  
  /* Set model data */
  ierr = pTatinModelSetUserData(m,NULL);CHKERRQ(ierr);
  
  /* Set function pointers */
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_StaticBox);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelInitialCondition_StaticBox);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_StaticBox);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_StaticBox);CHKERRQ(ierr);
  /*ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_StaticBox);CHKERRQ(ierr);*/
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_StaticBox);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_StaticBox);CHKERRQ(ierr);
  /*ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_StaticBox);CHKERRQ(ierr);*/
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_StaticBox);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_StaticBox);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

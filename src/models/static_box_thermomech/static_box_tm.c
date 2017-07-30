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
#include <ptatin_std_dirichlet_boundary_conditions.h>
#include <dmda_element_q2p1.h>
#include <material_point_utils.h>
#include <output_material_points.h>
#include <pswarm.h>
#include <ptatin3d_energy.h>
#include <energy_output.h>
#include <material_constants_energy.h>

PSwarm pswarm;

/*
 Coefficients for the energy uses a combination of region-wise (rho,Cp) and material-point-wised defined properties (kappa,H)
 Coefficients for the mechanics defined strictly via users definition of material-point properties
*/
#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_StaticBoxTM_variant1"
PetscErrorCode ModelInitialize_StaticBoxTM_variant1(pTatinCtx c,void *ctx)
{
  RheologyConstants       *rheology;
  DataBucket              materialconstants;
  DataField               PField;
  EnergyMaterialConstants *matconstants_e;
  int                     regionidx;
  double                  rho_ref,Cp;
  PetscErrorCode          ierr;
  
  PetscFunctionBegin;
  ierr = pTatinGetRheology(c,&rheology);CHKERRQ(ierr);
  rheology->rheology_type = RHEOLOGY_VISCOUS;
  
  ierr = PetscOptionsInsertString(NULL,"-activate_energy true");CHKERRQ(ierr);
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  ierr = MaterialConstantsSetDefaults(materialconstants);CHKERRQ(ierr);
  
  DataBucketGetDataFieldByName(materialconstants,EnergyMaterialConstants_classname,&PField); DataFieldGetEntries(PField,(void**)&matconstants_e);
  
  regionidx = 0;
  rho_ref = 1.0;
  Cp = 1.0;
  ierr = MaterialConstantsSetValues_EnergyMaterialConstants(regionidx,matconstants_e,0.0,0.0,rho_ref,Cp,ENERGYDENSITY_CONSTANT,ENERGYCONDUCTIVITY_USE_MATERIALPOINT_VALUE,NULL);CHKERRQ(ierr);
  EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[regionidx],ENERGYSOURCE_NONE);
  EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[regionidx],0,ENERGYSOURCE_USE_MATERIALPOINT_VALUE);
  
  ierr = PSwarmCreate(PETSC_COMM_WORLD,&pswarm);CHKERRQ(ierr);
  ierr = PSwarmSetOptionsPrefix(pswarm,"passive_");CHKERRQ(ierr);
  ierr = PSwarmSetPtatinCtx(pswarm,c);CHKERRQ(ierr);
  ierr = PSwarmSetTransportModeType(pswarm,PSWARM_TM_EULERIAN);CHKERRQ(ierr);
  ierr = PSwarmSetFromOptions(pswarm);CHKERRQ(ierr);
  
  if (c->mx != 4) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Model only valid for -mx 4");
  if (c->my != 4) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Model only valid for -my 4");
  if (c->mz != 4) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Model only valid for -mz 4");
  
  PetscFunctionReturn(0);
}

/*
 Coefficients for the energy prescribed region-wise
 Coefficients for the mechanics defined region-wise
*/
#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_StaticBoxTM_variant2"
PetscErrorCode ModelInitialize_StaticBoxTM_variant2(pTatinCtx c,void *ctx)
{
  RheologyConstants       *rheology;
  DataBucket              materialconstants;
  DataField               PField;
  EnergyMaterialConstants *matconstants_e;
  int                     regionidx;
  double                  rho,eta,k,Cp,Q;
  PetscBool               use_source = PETSC_FALSE;
  PetscErrorCode          ierr;
  
  PetscFunctionBegin;
  ierr = pTatinGetRheology(c,&rheology);CHKERRQ(ierr);
  rheology->rheology_type  = RHEOLOGY_VP_STD;
  rheology->nphases_active = 1;
  
  ierr = PetscOptionsInsertString(NULL,"-activate_energy true");CHKERRQ(ierr);
  
  ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
  ierr = MaterialConstantsSetDefaults(materialconstants);CHKERRQ(ierr);
  
  DataBucketGetDataFieldByName(materialconstants,EnergyMaterialConstants_classname,&PField); DataFieldGetEntries(PField,(void**)&matconstants_e);
  
  regionidx = 0;
  rho = 2.0;
  eta = 10.0;
  k   = 2.0e1;
  Cp  = 1.0;
  Q   = 0.0;
  ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_source",&use_source,NULL);CHKERRQ(ierr);
  if (use_source) {
    Q = 1000.0;
  }
  
  ierr = MaterialConstantsSetValues_EnergyMaterialConstants(regionidx,matconstants_e,0.0,0.0,rho,Cp,ENERGYDENSITY_CONSTANT,ENERGYCONDUCTIVITY_CONSTANT,NULL);CHKERRQ(ierr);
  EnergyMaterialConstantsSetFieldAll_SourceMethod(&matconstants_e[regionidx],ENERGYSOURCE_NONE);
  EnergyMaterialConstantsSetFieldByIndex_SourceMethod(&matconstants_e[regionidx],0,ENERGYSOURCE_CONSTANT);
  
  {
    EnergyConductivityConst *data;
    DataField               PField_;

    DataBucketGetDataFieldByName(materialconstants,EnergyConductivityConst_classname,&PField_);
    DataFieldGetEntries(PField_,(void**)&data);
    EnergyConductivityConstSetField_k0(&data[regionidx],k);
  }
  {
    EnergySourceConst *data;
    DataField         PField_;

    DataBucketGetDataFieldByName(materialconstants,EnergySourceConst_classname,&PField_);
    DataFieldGetEntries(PField_,(void**)&data);
    EnergySourceConstSetField_HeatSource(&data[regionidx],Q);
  }
  
  ierr = MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_CONSTANT,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);CHKERRQ(ierr);
  ierr = MaterialConstantsSetValues_ViscosityConst(materialconstants,regionidx,eta);CHKERRQ(ierr);
  ierr = MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,rho);CHKERRQ(ierr);

  
  ierr = PSwarmCreate(PETSC_COMM_WORLD,&pswarm);CHKERRQ(ierr);
  ierr = PSwarmSetOptionsPrefix(pswarm,"passive_");CHKERRQ(ierr);
  ierr = PSwarmSetPtatinCtx(pswarm,c);CHKERRQ(ierr);
  ierr = PSwarmSetTransportModeType(pswarm,PSWARM_TM_EULERIAN);CHKERRQ(ierr);
  ierr = PSwarmSetFromOptions(pswarm);CHKERRQ(ierr);
  
  if (c->mx != 4) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Model only valid for -mx 4");
  if (c->my != 4) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Model only valid for -my 4");
  if (c->mz != 4) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Model only valid for -mz 4");
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_StaticBoxTM"
PetscErrorCode ModelApplyBoundaryCondition_StaticBoxTM(pTatinCtx c,void *ctx)
{
  PetscScalar    zero = 0.0;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  /* Impose zero Neumann conditions (sigma.n = 0, sigma.t = 0) on top face */
  /* Impose freeslip conditions (u.n = 0, tau.t = 0) on left/right/front/back/bottom faces */
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  ierr = DMDABCListTraverse3d(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
  
  /* Impose Dirichlet conditions on bottom/top faces */
  {
    PetscReal      val_T;
    PhysCompEnergy energy;
    BCList         bclist;
    DM             daT;
    
    ierr   = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
    daT    = energy->daT;
    bclist = energy->T_bclist;
    
    val_T = 273.0;
    ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);
    
    val_T = 1000.0 + 273.0;
    ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_StaticBoxTM"
PetscErrorCode ModelApplyBoundaryConditionMG_StaticBoxTM(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
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

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_StaticBoxTM"
PetscErrorCode ModelApplyInitialMeshGeometry_StaticBoxTM(pTatinCtx c,void *ctx)
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

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_StaticBoxTM_variant1"
PetscErrorCode ModelApplyInitialMaterialGeometry_StaticBoxTM_variant1(pTatinCtx c,void *ctx)
{
  int            p,npoints;
  DataBucket     db;
  MPAccess       mpX;
  PetscReal      Q;
  PetscBool      use_source = PETSC_FALSE;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  Q = 0.0;
  ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_source",&use_source,NULL);CHKERRQ(ierr);
  if (use_source) {
    Q = 1000.0;
  }
  
  /* define properties on material points */
  ierr = pTatinGetMaterialPoints(c,&db,NULL);CHKERRQ(ierr);
  DataBucketGetSizes(db,&npoints,0,0);
  ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
  
  for (p=0; p<npoints; p++) {
    double *position;
    double eta,rho,kappa,H;
    int    phase;
    
    /* Access using the getter function provided for you (recommeneded for beginner user) */
    ierr = MaterialPointGet_global_coord(mpX,p,&position);CHKERRQ(ierr);
    
    phase = 0;
    
    eta = 10.0;
    rho = 2.0;
    
    kappa = 2.0e1/(rho * 1.0);
    H     = Q/(rho * 1.0);
    
    ierr = MaterialPointSet_phase_index(mpX,p,phase);CHKERRQ(ierr);
    
    ierr = MaterialPointSet_viscosity(mpX,p,eta);CHKERRQ(ierr);
    ierr = MaterialPointSet_density(mpX,p,rho);CHKERRQ(ierr);
    
    ierr = MaterialPointSet_diffusivity(mpX,p,kappa);CHKERRQ(ierr);
    ierr = MaterialPointSet_heat_source(mpX,p,H);CHKERRQ(ierr);
  }
  ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_StaticBoxTM_variant2"
PetscErrorCode ModelApplyInitialMaterialGeometry_StaticBoxTM_variant2(pTatinCtx c,void *ctx)
{
  int            p,npoints;
  DataBucket     db;
  MPAccess       mpX;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  /* define properties on material points */
  ierr = pTatinGetMaterialPoints(c,&db,NULL);CHKERRQ(ierr);
  DataBucketGetSizes(db,&npoints,0,0);
  ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
  
  for (p=0; p<npoints; p++) {
    int    phase;
    
    phase = 0;
    ierr = MaterialPointSet_phase_index(mpX,p,phase);CHKERRQ(ierr);
  }
  ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*
 depth  = 0
 rho0 = 2.0
 depth_interface = 3.0
 rho1 = 5.0
 */
#undef __FUNCT__
#define __FUNCT__ "ComputeHydrostaticPressure"
static PetscErrorCode ComputeHydrostaticPressure(PetscReal depth,PetscReal rho0,PetscReal *pressure)
{
  PetscFunctionBegin;
  *pressure = rho0 * 10.0 * depth;
  PetscFunctionReturn(0);
}

/* 
 Evaluates analytic solution of
   0 = k d^T/dt^2 + Q
 or
   T'' = -gamma, where gamma = Q/k
 over the domain
   y0=0 <= y <= y1=6
 with Dirichlet conditions
   T0 = T(y0)
   T1 = T(y1)
*/
#undef __FUNCT__
#define __FUNCT__ "Compute1DTemperature"
static PetscErrorCode Compute1DTemperature(PetscReal y,PetscReal Q,PetscReal *t,PetscReal *dtdy)
{
  PetscReal k,T0,T1,y1,gamma,A,B;
  
  PetscFunctionBegin;
  y1 = 6.0;
  T0 = 1273.0;
  T1 = 273.0;
  k = 2.0e1;
  
  gamma = Q/k;
  B = T0;
  A = (T1 - B + 0.5 * (gamma) * y1*y1)/y1;
  *t = -0.5 * (gamma) * y*y + A*y + B;
  *dtdy = -(gamma) * y + A;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_StaticBoxTM"
PetscErrorCode ModelOutput_StaticBoxTM(pTatinCtx c,Vec X,const char prefix[],void *ctx)
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
  {
    PhysCompEnergy energy;
    Vec            temperature;
    
    ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
    ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);
    ierr = pTatin3d_ModelOutput_Temperature_Energy(c,temperature,prefix);CHKERRQ(ierr);
  }
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
    PetscReal  rho0;
    PetscReal  max_abs_delta_p,max_abs_delta_p_g;
    
    ierr = PSwarmGetDataBucket(pswarm,&pdb);CHKERRQ(ierr);
    DataBucketGetSizes(pdb,&npoints,0,0);
    
    DataBucketGetDataFieldByName(pdb,MPntStd_classname,&datafield_tracers); DataFieldGetEntries(datafield_tracers,(void**)&tracer);
    DataBucketGetDataFieldByName(pdb,"pressure",&datafield);       DataFieldGetEntries(datafield,(void**)&tracer_pressure);
    
    max_abs_delta_p = 0.0;
    rho0 = 2.0;
    
    for (p=0; p<npoints; p++) {
      PetscReal depth,p_hydrostatic,delta;
      
      depth = 6.0 - (PetscReal)tracer[p].coor[1];
      ierr = ComputeHydrostaticPressure(depth,rho0,&p_hydrostatic);CHKERRQ(ierr);
      
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
  }
  
  /* Compare computed temperature with rho.g.z */
  {
    PhysCompEnergy energy;
    Vec            temperature,coor;
    DM             daT,daC;
    PetscScalar    ***T,****nodecoor;
    PetscInt       i,j,k,si,sj,sk,mi,mj,mk;
    PetscReal      max_abs_delta_t,max_abs_delta_t_g,delta,Q;
    PetscBool      use_source = PETSC_FALSE;

    Q = 0.0;
    ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_source",&use_source,NULL);CHKERRQ(ierr);
    if (use_source) {
      Q = 1000.0;
    }

    ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
    daT = energy->daT;
    ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);
    ierr = DMGetCoordinates(daT,&coor);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(daT,&daC);CHKERRQ(ierr);
    ierr = DMDAGetCorners(daT,&si,&sj,&sk,&mi,&mj,&mk);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(daT,temperature,&T);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(daC,coor,&nodecoor);CHKERRQ(ierr);

    max_abs_delta_t = 0.0;
    for (k=sk; k<(sk+mk); k++) {
      for (j=sj; j<(sj+mj); j++) {
        for (i=si; i<(si+mi); i++) {
          PetscReal y,T_numeric,T_exact,gradT_exact;
          
          y = nodecoor[k][j][i][1];
          ierr = Compute1DTemperature(y,Q,&T_exact,&gradT_exact);CHKERRQ(ierr);
          T_numeric = T[k][j][i];
          
          /*PetscPrintf(PETSC_COMM_WORLD,"%D %+1.12e %+1.12e %+1.12e\n",j,y,T_numeric,T_exact);*/
          delta = PetscAbsReal(T_numeric - T_exact);
          if (delta > max_abs_delta_t) {
            max_abs_delta_t = delta;
          }
        }
      }
    }
    ierr = DMDAVecRestoreArrayDOF(daC,coor,&nodecoor);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(daT,temperature,&T);CHKERRQ(ierr);
    
    ierr = MPI_Allreduce(&max_abs_delta_t,&max_abs_delta_t_g,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
    if (max_abs_delta_t_g > 1.0e-12) {
      PetscPrintf(PETSC_COMM_WORLD,"[staticBox] ||T - T_exact||_inf = %+1.12e <fail>\n",max_abs_delta_t_g);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,"[staticBox] ||T - T_exact||_inf <= 1.0e-12 <pass>\n");
    }
    PetscPrintf(PETSC_COMM_WORLD,"[staticBox] Q = %+1.2e\n",Q);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialCondition_StaticBoxTM"
PetscErrorCode ModelInitialCondition_StaticBoxTM(pTatinCtx c,Vec X,void *ctx)
{
  PhysCompStokes stokes_ctx;
  DM             stokes_pack,dau,dap;
  Vec            velocity,pressure;
  PetscBool      phys_component_energy_valid;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = pTatinGetStokesContext(c,&stokes_ctx);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes_ctx,&stokes_pack);CHKERRQ(ierr);
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
  ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  
  ierr = pTatinContextValid_Energy(c,&phys_component_energy_valid);CHKERRQ(ierr);
  if (phys_component_energy_valid) {
    PhysCompEnergy energy;
    Vec            temperature;
    
    ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
    ierr = pTatinPhysCompGetData_Energy(c,&temperature,NULL);CHKERRQ(ierr);
    ierr = VecSet(temperature,273.0);CHKERRQ(ierr);
  } else SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Model only valid if option -activate_energy true is used");

  ierr = PSwarmAttachStateVecVelocityPressure(pswarm,X);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_StaticBoxTM"
PetscErrorCode ModelDestroy_StaticBoxTM(pTatinCtx c,void *ctx)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[ModelDestroy_StaticBoxTM]]\n");
  ierr = PSwarmDestroy(&pswarm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelCreate_StaticBoxTM"
PetscErrorCode pTatinModelCreate_StaticBoxTM(pTatinModel m)
{
  PetscBool use_v1 = PETSC_FALSE;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  /* Allocate memory for the data structure for this model */
  
  /* Set model data */
  ierr = pTatinModelSetUserData(m,NULL);CHKERRQ(ierr);
  
  /* Set function pointers */
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelInitialCondition_StaticBoxTM);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_StaticBoxTM);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_StaticBoxTM);CHKERRQ(ierr);
  /*ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_StaticBox);CHKERRQ(ierr);*/
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_StaticBoxTM);CHKERRQ(ierr);
  /*ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_StaticBox);CHKERRQ(ierr);*/
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_StaticBoxTM);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_StaticBoxTM);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,NULL,"-use_model_variant1",&use_v1,NULL);CHKERRQ(ierr);
  if (use_v1) {
    ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_StaticBoxTM_variant1);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_StaticBoxTM_variant1);CHKERRQ(ierr);
  } else {
    ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_StaticBoxTM_variant2);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_StaticBoxTM_variant2);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

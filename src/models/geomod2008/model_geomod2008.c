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
 **    filename:   model_geomod2008.c
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
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_models.h"
#include "ptatin3d_stokes.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "material_point_utils.h"
#include "dmda_update_coords.h"
#include "material_point_std_utils.h"
#include "material_point_popcontrol.h"
#include "mesh_update.h"

typedef struct _p_ModelCtxGeoMod2008 *ModelCtxGeoMod2008;
struct _p_ModelCtxGeoMod2008 {
  PetscInt experiment;
  PetscBool three_dimensional;
  PetscBool use_free_surface;
  PetscReal Lx,Ly,Lz;
  PetscReal vx_bc;
  PetscReal frictional_boundary_layer_delta;

  PetscReal quartz_eta,quartz_rho;
  PetscReal quartz_C0,quartz_mu;
  PetscReal quartz_C0_inf,quartz_mu_inf;

  PetscReal quartz_bdy_C0,quartz_bdy_mu;
  PetscReal quartz_bdy_C0_inf,quartz_bdy_mu_inf;

  PetscReal corundum_eta,corundum_rho;
  PetscReal corundum_C0,corundum_mu;
  PetscReal corundum_C0_inf,corundum_mu_inf;

  PetscReal corundum_bdy_C0,corundum_bdy_mu;
  PetscReal corundum_bdy_C0_inf,corundum_bdy_mu_inf;

  PetscReal air_eta,air_rho;
};

const enum { RegionAir=0, RegionQuartz, RegionQuartzBdy, RegionCorundum, RegionCorundumBdy, RegionQuartzBdySkin, RegionCorundumBdySkin } RegionTags;


PetscReal cm2m    = 1.0e-2;
PetscReal h2s     = 60.0 * 60.0;
PetscReal deg2rad = M_PI/180.0;


PetscErrorCode ModelInitialize_GeoMod2008(pTatinCtx c,void *ctx)
{
  ModelCtxGeoMod2008  data = (ModelCtxGeoMod2008)ctx;
  RheologyConstants   *rheology;
  DataBucket          materialconstants = c->material_constants;
  PetscBool flg;
  PetscReal fac,length_scale,velocity_scale,time_scale,viscosity_scale,density_scale,pressure_scale;
  const PetscReal one_third = 0.333333333333333;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* experiment info */
  data->experiment = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-gm08_experiment",&data->experiment,&flg);CHKERRQ(ierr);
  data->three_dimensional = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-gm08_three_dimensional",&data->three_dimensional,&flg);CHKERRQ(ierr);
  data->use_free_surface = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-gm08_use_free_surface",&data->use_free_surface,&flg);CHKERRQ(ierr);
  data->frictional_boundary_layer_delta = 0.2 * cm2m;
  ierr = PetscOptionsGetReal(NULL,NULL,"-gm08_skin_thickness",&data->frictional_boundary_layer_delta,&flg);CHKERRQ(ierr);

  /* materials */
  fac = 1.0;
  data->quartz_eta    = 1.0e12;  /* Pa s */
  data->quartz_rho    = 1560.0; /* kg/m3 */
  data->quartz_C0     = 30.0*fac;    /* Pa */
  data->quartz_C0_inf = 30.0;
  data->quartz_mu     = 36.0 * deg2rad;    /* deg */
  data->quartz_mu_inf = 31.0 * deg2rad;

  data->quartz_bdy_C0     = 30.0*fac;
  data->quartz_bdy_C0_inf = 30.0;
  data->quartz_bdy_mu     = 16.0 * deg2rad;
  data->quartz_bdy_mu_inf = 14.0 * deg2rad;


  data->corundum_eta    = 1.0e12; /* Pa s */
  data->corundum_rho    = 1890.0; /* kg/m3 */
  data->corundum_C0     = 30.0*fac;   /* Pa */
  data->corundum_C0_inf = 30.0;
  data->corundum_mu     = 36.0 * deg2rad;   /* deg */
  data->corundum_mu_inf = 31.0 * deg2rad;

  data->corundum_bdy_C0     = 30.0*fac;
  data->corundum_bdy_C0_inf = 30.0;
  data->corundum_bdy_mu     = 24.0 * deg2rad;
  data->corundum_bdy_mu_inf = 23.0 * deg2rad;


  data->air_eta = 1.0e4;
  data->air_rho = 0.0;

  /* Material constant */
  MaterialConstantsSetDefaults(materialconstants);

  /* air */
  MaterialConstantsSetValues_MaterialType(materialconstants,  RegionAir,VISCOUS_CONSTANT,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);
  MaterialConstantsSetValues_ViscosityConst(materialconstants,RegionAir,data->air_eta);
  MaterialConstantsSetValues_DensityConst(materialconstants,  RegionAir,data->air_rho);

  /* Quartz */
  MaterialConstantsSetValues_MaterialType(materialconstants,  RegionQuartz, VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
  MaterialConstantsSetValues_ViscosityConst(materialconstants,RegionQuartz, data->quartz_eta);
  MaterialConstantsSetValues_DensityConst(materialconstants,  RegionQuartz, data->quartz_rho);
  MaterialConstantsSetValues_PlasticDP(materialconstants,     RegionQuartz, data->quartz_mu,data->quartz_mu_inf, data->quartz_C0,data->quartz_C0_inf, one_third*data->quartz_C0,1.0e20);
  MaterialConstantsSetValues_SoftLin(materialconstants,       RegionQuartz, 0.5,1.0);

  MaterialConstantsSetValues_MaterialType(materialconstants,  RegionQuartzBdy, VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
  MaterialConstantsSetValues_ViscosityConst(materialconstants,RegionQuartzBdy, data->quartz_eta);
  MaterialConstantsSetValues_DensityConst(materialconstants,  RegionQuartzBdy, data->quartz_rho);
  MaterialConstantsSetValues_PlasticDP(materialconstants,     RegionQuartzBdy, data->quartz_bdy_mu,data->quartz_bdy_mu_inf, data->quartz_bdy_C0,data->quartz_bdy_C0_inf, one_third*data->quartz_bdy_C0,1.0e20);
  MaterialConstantsSetValues_SoftLin(materialconstants,       RegionQuartzBdy, 0.5,1.0);

  MaterialConstantsSetValues_MaterialType(materialconstants,  RegionQuartzBdySkin, VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
  MaterialConstantsSetValues_ViscosityConst(materialconstants,RegionQuartzBdySkin, data->quartz_eta);
  MaterialConstantsSetValues_DensityConst(materialconstants,  RegionQuartzBdySkin, data->quartz_rho);
  MaterialConstantsSetValues_PlasticDP(materialconstants,     RegionQuartzBdySkin, data->quartz_bdy_mu,data->quartz_bdy_mu_inf, data->quartz_bdy_C0,data->quartz_bdy_C0_inf, one_third*data->quartz_bdy_C0,1.0e20);
  MaterialConstantsSetValues_SoftLin(materialconstants,       RegionQuartzBdySkin, 0.5,1.0);

  /* Corundum */
  MaterialConstantsSetValues_MaterialType(materialconstants,  RegionCorundum, VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
  MaterialConstantsSetValues_ViscosityConst(materialconstants,RegionCorundum, data->corundum_eta);
  MaterialConstantsSetValues_DensityConst(materialconstants,  RegionCorundum, data->corundum_rho);
  MaterialConstantsSetValues_PlasticDP(materialconstants,     RegionCorundum, data->corundum_mu,data->corundum_mu_inf, data->corundum_C0,data->corundum_C0_inf, one_third*data->corundum_C0,1.0e20);
  MaterialConstantsSetValues_SoftLin(materialconstants,       RegionCorundum, 0.5,1.0);

  MaterialConstantsSetValues_MaterialType(materialconstants,  RegionCorundumBdy, VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
  MaterialConstantsSetValues_ViscosityConst(materialconstants,RegionCorundumBdy, data->corundum_eta);
  MaterialConstantsSetValues_DensityConst(materialconstants,  RegionCorundumBdy, data->corundum_rho);
  MaterialConstantsSetValues_PlasticDP(materialconstants,     RegionCorundumBdy, data->corundum_bdy_mu,data->corundum_bdy_mu_inf, data->corundum_bdy_C0,data->corundum_bdy_C0_inf, one_third*data->corundum_bdy_C0,1.0e20);
  MaterialConstantsSetValues_SoftLin(materialconstants,       RegionCorundumBdy, 0.5,1.0);

  MaterialConstantsSetValues_MaterialType(materialconstants,  RegionCorundumBdySkin, VISCOUS_CONSTANT,PLASTIC_DP,SOFTENING_LINEAR,DENSITY_CONSTANT);
  MaterialConstantsSetValues_ViscosityConst(materialconstants,RegionCorundumBdySkin, data->corundum_eta);
  MaterialConstantsSetValues_DensityConst(materialconstants,  RegionCorundumBdySkin, data->corundum_rho);
  MaterialConstantsSetValues_PlasticDP(materialconstants,     RegionCorundumBdySkin, data->corundum_bdy_mu,data->corundum_bdy_mu_inf, data->corundum_bdy_C0,data->corundum_bdy_C0_inf, one_third*data->corundum_bdy_C0,1.0e20);
  MaterialConstantsSetValues_SoftLin(materialconstants,       RegionCorundumBdySkin, 0.5,1.0);

  rheology                = &c->rheology_constants;
  rheology->rheology_type = RHEOLOGY_VP_STD;
  rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
  rheology->eta_upper_cutoff_global = 1.0e+25;
  rheology->eta_lower_cutoff_global = 1.0e-20;

  /*
  switch (data->experiment) {
    case 1:
      rheology->nphases_active = 3;
      break;
    case 2:
      rheology->nphases_active = 5;
      break;
    case 3:
      rheology->nphases_active = 5;
      break;
  }
  */
  rheology->nphases_active = 7;

  length_scale    = 1.0;
  velocity_scale  = 1.0e-5;
  viscosity_scale = 1.0e12;

  time_scale      = length_scale / velocity_scale;
  pressure_scale  = viscosity_scale / time_scale;
  density_scale   = pressure_scale / length_scale;

  PetscPrintf(PETSC_COMM_WORLD,"geomod2008: length scale    = %1.4e \n", length_scale);
  PetscPrintf(PETSC_COMM_WORLD,"geomod2008: velocity scale  = %1.4e \n", velocity_scale);
  PetscPrintf(PETSC_COMM_WORLD,"geomod2008: viscosity scale = %1.4e \n", viscosity_scale);
  PetscPrintf(PETSC_COMM_WORLD,"geomod2008: time scale      = %1.4e \n", time_scale);
  PetscPrintf(PETSC_COMM_WORLD,"geomod2008: pressure scale  = %1.4e \n", pressure_scale);
  PetscPrintf(PETSC_COMM_WORLD,"geomod2008: density scale   = %1.4e \n", density_scale);

  // scale material properties
  {
    PetscInt                   regionidx;

    for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
      MaterialConstantsScaleAll(materialconstants,regionidx,length_scale,velocity_scale,time_scale,viscosity_scale,density_scale,pressure_scale);
    }
  }

  /* experiment info */
  data->vx_bc = - 2.5 * cm2m / h2s;
  //data->vx_bc = 0.0 * cm2m / h2s;
  data->vx_bc = data->vx_bc / velocity_scale;

  switch (data->experiment) {

    case 1:
      data->Lx = (5.0 + 8.24) * cm2m;
      data->Ly = (3.0 + 2.0) * cm2m; /* direction of g */
      data->Lz = (10.0) * cm2m;
      break;

    case 2:
      data->Lx = (35.0) * cm2m;
      data->Ly = (3.0 + 2.0) * cm2m; /* direction of g */
      data->Lz = (10.0) * cm2m;
      break;

    case 3:
      data->Lx = (35.0) * cm2m;
      data->Ly = (3.0 + 2.0) * cm2m; /* direction of g */
      data->Lz = (10.0) * cm2m;
      break;

    default:
      break;
  }

  if (data->use_free_surface) {
    data->Ly = (3.0) * cm2m; /* direction of g */
  }

  if (!data->three_dimensional) {
    PetscReal dx;
    PhysCompStokes stokes;

    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);

    ierr = PetscOptionsInsertString(NULL,"-stk_velocity_da_refine_hierarchy_z 1,1,1,1,1,1,1,1");CHKERRQ(ierr);
    ierr = PetscOptionsInsertString(NULL,"-stk_velocity_da_refine_z 1");CHKERRQ(ierr);
    //ierr = DMDASetRefinementFactor(stokes->dav,PETSC_IGNORE,PETSC_IGNORE,1);CHKERRQ(ierr);

    c->mz = 1;
    dx = data->Lx / ( (PetscReal)c->mx );
    data->Lz = dx;
  }








  PetscFunctionReturn(0);
}


PetscErrorCode ApplyStokesVelocityBC_GeoMod2008(pTatinCtx c,void *ctx,BCList bclist,DM dav)
{
  ModelCtxGeoMod2008 data = (ModelCtxGeoMod2008)ctx;
  PetscReal vx,vy,vz;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);


  vx = data->vx_bc;

  switch (data->experiment) {

    case 1:
      // Experiment 1
      //vx = 0.0;
      //ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&vx);CHKERRQ(ierr);
      vy = 0.0;
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,vy);CHKERRQ(ierr);

      vx = data->vx_bc;
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,EAST_FACE,vx);CHKERRQ(ierr);
      //vy = 0.0;
      //ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&vy);CHKERRQ(ierr);

      vz = 0.0;
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,FRONT_FACE,vz);CHKERRQ(ierr);
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,BACK_FACE,vz);CHKERRQ(ierr);

      /* sticky air sides - option to leave open? */
      vx = vy = 0.0;
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,WEST_FACE,vx);CHKERRQ(ierr);
      //ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,NORTH_FACE,vy);CHKERRQ(ierr);



      break;

    case 2:
      // Experiment 1
      vx = data->vx_bc;
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,EAST_FACE,vx);CHKERRQ(ierr);

      vy = 0.0;
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,vy);CHKERRQ(ierr);

      vz = 0.0;
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,FRONT_FACE,vz);CHKERRQ(ierr);
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,BACK_FACE,vz);CHKERRQ(ierr);

      /* sticky air sides - option to leave open? */
      vx = vy = 0.0;
      ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,WEST_FACE,vx);CHKERRQ(ierr);
      //ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,NORTH_FACE,vy);CHKERRQ(ierr);
      break;

    case 3:
      break;

    default:
      break;
  }



  PetscFunctionReturn(0);
}


/* wrapper function for applying boundary condition on fine grid */
PetscErrorCode ModelApplyBoundaryCondition_GeoMod2008(pTatinCtx c,void *ctx)
{
  PhysCompStokes stokes;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* stokes dirichlet boundary conditions on u */
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  ierr = ApplyStokesVelocityBC_GeoMod2008(c,ctx,stokes->u_bclist,stokes->dav);CHKERRQ(ierr);

  /* any other physiscs? - apply bcs here */

  PetscFunctionReturn(0);
}

/* wrapper function for applying boundary condition on all multi-grid levels */
PetscErrorCode ModelApplyStokesVelocityBoundaryConditionMG_GeoMod2008(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
  PetscInt       n;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  for (n=0; n<nl; n++) {
    /* Define boundary conditions for each level in the MG hierarchy */
    ierr = ApplyStokesVelocityBC_GeoMod2008(user,ctx,bclist[n],dav[n]);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode GeomMod2008ApplyFrictionalBoundarySkin(ModelCtxGeoMod2008 data,DM dav,DataBucket material_points)
{
  PetscInt         M,N,P;
  PetscReal        gmin[3],gmax[3];
  PetscReal        dx,dy,dz,delta,dl;
  MPAccess         mpX;
  int              p,n_mp_points;
  PetscBool        apply_skin;
  PetscErrorCode   ierr;


  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* reset all skin particles */
  DataBucketGetSizes(material_points,&n_mp_points,0,0);
  ierr = MaterialPointGetAccess(material_points,&mpX);CHKERRQ(ierr);
  for (p=0; p<n_mp_points; p++) {
    int region_index;

    ierr = MaterialPointGet_phase_index(mpX,p,&region_index);CHKERRQ(ierr);

    if (region_index == RegionQuartzBdySkin) {
      ierr = MaterialPointSet_phase_index(mpX,p,RegionQuartz);CHKERRQ(ierr);
    } else if (region_index == RegionCorundumBdySkin) {
      ierr = MaterialPointSet_phase_index(mpX,p,RegionCorundum);CHKERRQ(ierr);
    }

  }
  ierr = MaterialPointRestoreAccess(material_points,&mpX);CHKERRQ(ierr);


  /* compute dx */
  ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
  dx = (gmax[0]-gmin[0])/((PetscReal)((M-1)/2));
  dy = (gmax[1]-gmin[1])/((PetscReal)((N-1)/2));
  dz = (gmax[2]-gmin[2])/((PetscReal)((P-1)/2));

  delta = dx; if (dy > delta) { delta = dy; } if (dz > delta) { delta = dz; }
  dl = data->frictional_boundary_layer_delta;
  if ( data->frictional_boundary_layer_delta < 2.0 * delta) {
    PetscPrintf(PETSC_COMM_WORLD,"WARNING (geomod2008): skin thickness is smaller than 2 elements\n");
    PetscPrintf(PETSC_COMM_WORLD,"WARNING (geomod2008): suggest you re-run with -gm08_skin_thickness %1.4e\n", 2.0 * delta);
  }

  PetscPrintf(PETSC_COMM_WORLD,"geomod2008: dx,dy,dz = (%1.4e , %1.4e , %1.4e); skin thickness %1.4e \n", dx,dy,dz,dl );

  /* check for delta epsilon distance from boundary */
  /* reset all skin particles */
  DataBucketGetSizes(material_points,&n_mp_points,0,0);
  ierr = MaterialPointGetAccess(material_points,&mpX);CHKERRQ(ierr);
  for (p=0; p<n_mp_points; p++) {
    int region_index;
    double *xp;

    ierr = MaterialPointGet_phase_index(mpX,p,&region_index);CHKERRQ(ierr);
    ierr = MaterialPointGet_global_coord(mpX,p,&xp);CHKERRQ(ierr);

    apply_skin = PETSC_FALSE;
    if ( fabs(xp[0] - gmin[0]) < dl ) { apply_skin = PETSC_TRUE; }
    if ( fabs(xp[0] - gmax[0]) < dl ) { apply_skin = PETSC_TRUE; }
    if ( fabs(xp[1] - gmin[1]) < dl ) { apply_skin = PETSC_TRUE; }

    if (data->three_dimensional) {
      if ( fabs(xp[2] - gmin[2]) < dl ) { apply_skin = PETSC_TRUE; }
      if ( fabs(xp[2] - gmax[2]) < dl ) { apply_skin = PETSC_TRUE; }
    }


    if (apply_skin) {
      if (region_index == RegionQuartz) {
        ierr = MaterialPointSet_phase_index(mpX,p,RegionQuartzBdySkin);CHKERRQ(ierr);
      } else if (region_index == RegionCorundum) {
        ierr = MaterialPointSet_phase_index(mpX,p,RegionCorundumBdySkin);CHKERRQ(ierr);
      }
    }

  }
  ierr = MaterialPointRestoreAccess(material_points,&mpX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelAdaptMaterialPointResolution_GeoMod2008(pTatinCtx c,void *ctx)
{
  ModelCtxGeoMod2008 data = (ModelCtxGeoMod2008)ctx;
  DataBucket         material_points;
  PhysCompStokes     stokes;

  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = pTatinGetMaterialPoints(c,&material_points,NULL);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);

  ierr = GeomMod2008ApplyFrictionalBoundarySkin(data,stokes->dav,material_points);CHKERRQ(ierr);
  /* Perform injection and cleanup of markers */
  ierr = MaterialPointPopulationControl_v1(c);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMeshGeometry_GeoMod2008(pTatinCtx c,void *ctx)
{
  ModelCtxGeoMod2008 data = (ModelCtxGeoMod2008)ctx;
  PhysCompStokes stokes;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(stokes->dav,0.0,data->Lx,0.0,data->Ly,0.0,data->Lz);CHKERRQ(ierr);
  {
    PetscReal gvec[] = { 0.0, -9.81, 0.0 };
    ierr = PhysCompStokesSetGravityVector(stokes,gvec);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_GeoMod2008_exp1(pTatinCtx c,void *ctx)
{
  ModelCtxGeoMod2008 data = (ModelCtxGeoMod2008)ctx;
  DataBucket     material_points;
  MPAccess       mpX;
  int            p,n_mp_points;
  PhysCompStokes     stokes;
  DM                 stokes_pack,dav,dap;
  DataBucket     material_points_face;
  PetscInt       face_index;
  PetscInt       Nxp[2];
  PetscReal      perturb;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = pTatinGetMaterialPoints(c,&material_points,NULL);CHKERRQ(ierr);
  DataBucketGetSizes(material_points,&n_mp_points,0,0);

  ierr = MaterialPointGetAccess(material_points,&mpX);CHKERRQ(ierr);
  for (p=0; p<n_mp_points; p++) {
    int region_index;
    double viscosity,density,*xp;

    region_index = RegionAir;
    viscosity    = data->air_eta;
    density      = data->air_rho;

    ierr = MaterialPointGet_global_coord(mpX,p,&xp);CHKERRQ(ierr);

    // Experiment 1
    if (xp[0] > 5.0 * cm2m) {
      double angle = 20.0 * M_PI/180.0;
      double h,shift_x;

      shift_x = xp[0] - 5.0 * cm2m;
      h = shift_x * tan( angle );

      if (xp[1] <= h) {
        region_index = RegionQuartz;
        viscosity    = data->quartz_eta;
        density      = data->quartz_rho;
      }
    }

    viscosity = viscosity / 1.0e10;
    density   = density / 1.0e10;

    ierr = MaterialPointSet_phase_index(mpX,p,region_index);CHKERRQ(ierr);
    //ierr = MaterialPointSet_viscosity(mpX,p,viscosity);CHKERRQ(ierr);
    //ierr = MaterialPointSet_density(mpX,p,density);CHKERRQ(ierr);
  }
  ierr = MaterialPointRestoreAccess(material_points,&mpX);CHKERRQ(ierr);


  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

  /* add boundary points along base */
  face_index = 3;
  DataBucketDuplicateFields(material_points,&material_points_face);
  /* reset size */
  DataBucketSetSizes(material_points_face,0,-1);
  /* assign coords */
  Nxp[0]  = 8;
  Nxp[1]  = 8;
  perturb = 0.1;
  ierr = SwarmMPntStd_CoordAssignment_FaceLatticeLayout3d(dav,Nxp,perturb,face_index,material_points_face);CHKERRQ(ierr);

  /* assign values */
  DataBucketGetSizes(material_points_face,&n_mp_points,0,0);
  ierr = MaterialPointGetAccess(material_points_face,&mpX);CHKERRQ(ierr);
  for (p=0; p<n_mp_points; p++) {
    int region_index;
    double viscosity,density,*xp;

    ierr = MaterialPointGet_global_coord(mpX,p,&xp);CHKERRQ(ierr);

    region_index = RegionAir;
    viscosity    = data->air_eta;
    density      = data->air_rho;

    if (xp[0] >= 5.0 * cm2m) {
      region_index = RegionQuartzBdy;
      viscosity    = data->quartz_eta;
      density      = data->quartz_rho;
    }

    density = -9.81 * density;

    viscosity = viscosity / 1.0e10;
    density   = density / 1.0e10;

    ierr = MaterialPointSet_phase_index(mpX,p,region_index);CHKERRQ(ierr);
    //ierr = MaterialPointSet_viscosity(mpX,p,viscosity);CHKERRQ(ierr);
    //ierr = MaterialPointSet_density(mpX,p,density);CHKERRQ(ierr);
  }
  ierr = MaterialPointRestoreAccess(material_points_face,&mpX);CHKERRQ(ierr);
  /* insert into volume bucket */
  DataBucketInsertValues(material_points,material_points_face);
  /* delete */
  DataBucketDestroy(&material_points_face);




  /* add boundary points along right */
  face_index = 0;
  DataBucketDuplicateFields(material_points,&material_points_face);
  /* reset size */
  DataBucketSetSizes(material_points_face,0,-1);
  /* assign coords */
  Nxp[0]  = 8;
  Nxp[1]  = 8;
  perturb = 0.1;
  ierr = SwarmMPntStd_CoordAssignment_FaceLatticeLayout3d(dav,Nxp,perturb,face_index,material_points_face);CHKERRQ(ierr);

  /* assign values */
  DataBucketGetSizes(material_points_face,&n_mp_points,0,0);
  ierr = MaterialPointGetAccess(material_points_face,&mpX);CHKERRQ(ierr);
  for (p=0; p<n_mp_points; p++) {
    int region_index;
    double viscosity,density,*xp;

    ierr = MaterialPointGet_global_coord(mpX,p,&xp);CHKERRQ(ierr);

    region_index = RegionAir;
    viscosity    = data->air_eta;
    density      = data->air_rho;

    if (xp[1] <= 3.0 * cm2m) {
      region_index = RegionQuartzBdy;
      viscosity    = data->quartz_eta;
      density      = data->quartz_rho;
    }

    viscosity = viscosity / 1.0e10;
    density   = density / 1.0e10;

    ierr = MaterialPointSet_phase_index(mpX,p,region_index);CHKERRQ(ierr);
    //ierr = MaterialPointSet_viscosity(mpX,p,viscosity);CHKERRQ(ierr);
    //ierr = MaterialPointSet_density(mpX,p,density);CHKERRQ(ierr);
  }
  ierr = MaterialPointRestoreAccess(material_points_face,&mpX);CHKERRQ(ierr);
  /* insert into volume bucket */
  DataBucketInsertValues(material_points,material_points_face);
  /* delete */
  DataBucketDestroy(&material_points_face);


  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_GeoMod2008_exp2(pTatinCtx c,void *ctx)
{
  DataBucket     material_points;
  MPAccess       mpX;
  int            p,n_mp_points;
  PhysCompStokes     stokes;
  DM                 stokes_pack,dav,dap;
  DataBucket     material_points_face;
  PetscInt       face,face_index[3];
  PetscInt       Nxp[2];
  PetscReal      perturb;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = pTatinGetMaterialPoints(c,&material_points,NULL);CHKERRQ(ierr);
  DataBucketGetSizes(material_points,&n_mp_points,0,0);

  ierr = MaterialPointGetAccess(material_points,&mpX);CHKERRQ(ierr);
  for (p=0; p<n_mp_points; p++) {
    int region_index;
    double *xp,y_pos;

    ierr = MaterialPointGet_global_coord(mpX,p,&xp);CHKERRQ(ierr);
    y_pos = xp[1] / cm2m;

    // Experiment 2
    region_index = RegionAir;
    if (y_pos < 1.0) {
      region_index = RegionQuartz;
    }
    if ( (y_pos >= 1.0) && (y_pos <= 3.0) ) {
      region_index = RegionCorundum;
    }
    if ( (y_pos > 2.0) && (y_pos <= 3.0) ) {
      region_index = RegionQuartz;
    }


    if (y_pos < 0.1) {
    //  region_index = RegionQuartzBdy;
    }

    ierr = MaterialPointSet_phase_index(mpX,p,region_index);CHKERRQ(ierr);
  }
  ierr = MaterialPointRestoreAccess(material_points,&mpX);CHKERRQ(ierr);

  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

  /* add boundary points along base, left and right walls */

  face_index[0] = 0;
  face_index[1] = 1;
  face_index[2] = 3;

  for (face=0; face<3; face++) {
    DataBucketDuplicateFields(material_points,&material_points_face);
    /* reset size */
    DataBucketSetSizes(material_points_face,0,-1);
    /* assign coords */
    Nxp[0]  = 8;
    Nxp[1]  = 8;
    perturb = 0.1;
    ierr = SwarmMPntStd_CoordAssignment_FaceLatticeLayout3d(dav,Nxp,perturb,face_index[face],material_points_face);CHKERRQ(ierr);

    /* assign values */
    DataBucketGetSizes(material_points_face,&n_mp_points,0,0);
    ierr = MaterialPointGetAccess(material_points_face,&mpX);CHKERRQ(ierr);
    for (p=0; p<n_mp_points; p++) {
      int region_index;
      double *xp,y_pos;

      ierr = MaterialPointGet_global_coord(mpX,p,&xp);CHKERRQ(ierr);
      y_pos = xp[1] / cm2m;

      region_index = RegionAir;
      if (y_pos < 1.0) {
        region_index = RegionQuartzBdy;
      }
      if ( (y_pos >= 1.0) && (y_pos <= 3.0) ) {
        region_index = RegionCorundumBdy;
      }
      if ( (y_pos > 2.0) && (y_pos <= 3.0) ) {
        region_index = RegionQuartzBdy;
      }
      ierr = MaterialPointSet_phase_index(mpX,p,region_index);CHKERRQ(ierr);
    }
    ierr = MaterialPointRestoreAccess(material_points_face,&mpX);CHKERRQ(ierr);
    /* insert into volume bucket */
    DataBucketInsertValues(material_points,material_points_face);
    /* delete */
    DataBucketDestroy(&material_points_face);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_GeoMod2008(pTatinCtx c,void *ctx)
{
  ModelCtxGeoMod2008 data = (ModelCtxGeoMod2008)ctx;
  DataBucket         material_points;
  PhysCompStokes     stokes;
  PetscErrorCode     ierr;

  PetscFunctionBegin;

  switch (data->experiment) {

    case 1:
      ierr = ModelApplyInitialMaterialGeometry_GeoMod2008_exp1(c,ctx);CHKERRQ(ierr);
      break;

    case 2:
      ierr = ModelApplyInitialMaterialGeometry_GeoMod2008_exp2(c,ctx);CHKERRQ(ierr);
      break;
  }

  ierr = pTatinGetMaterialPoints(c,&material_points,NULL);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  ierr = GeomMod2008ApplyFrictionalBoundarySkin(data,stokes->dav,material_points);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialSolution_GeoMod2008(pTatinCtx c,Vec X,void *ctx)
{
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyUpdateMeshGeometry_GeoMod2008(pTatinCtx c,Vec X,void *ctx)
{
  ModelCtxGeoMod2008 data = (ModelCtxGeoMod2008)ctx;
  PetscReal        step;
  PhysCompStokes   stokes;
  DM               stokes_pack,dav,dap;
  Vec              velocity,pressure;
  PetscReal        gmin[3],gmax[3];
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* fully lagrangian update */
  ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);

  stokes_pack = stokes->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

  if (!data->use_free_surface) {
    ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

    ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates(dav,0.0,gmax[0],0.0,data->Ly,0.0,gmax[2]);CHKERRQ(ierr);
    ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
  } else {
    ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

    ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialStokesVariableMarkers_GeoMod2008(pTatinCtx user,Vec X,void *ctx)
{

  PhysCompStokes    stokes;
  DM                stokes_pack,dau,dap;
  Vec               Uloc,Ploc;
  PetscScalar       *LA_Uloc,*LA_Ploc;
  PetscErrorCode    ierr;
  PetscFunctionBegin;

  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;

  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);

  ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
  ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
  ierr = pTatin_EvaluateRheologyNonlinearities(user,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);

  ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);

  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelOutput_GeoMod2008(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
  DataBucket     materialpoint_db;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* mesh */
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

  /* markers */
  ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
  {
    const int nf = 2;
    const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
    char mp_file_prefix[256];

    sprintf(mp_file_prefix,"%s_mpoints",prefix);
    ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelDestroy_GeoMod2008(pTatinCtx c,void *ctx)
{
  ModelCtxGeoMod2008 data = (ModelCtxGeoMod2008)ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* Free contents of structure */

  /* Free structure */
  ierr = PetscFree(data);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinModelRegister_GeoMod2008(void)
{
  ModelCtxGeoMod2008 data;
  pTatinModel m;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* register user model */
  ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

  /* Set name, model select via -ptatin_model NAME */
  ierr = pTatinModelSetName(m,"geomod2008");CHKERRQ(ierr);

  /* Set model data */
  ierr = PetscMalloc(sizeof(struct _p_ModelCtxGeoMod2008),&data);CHKERRQ(ierr);
  ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);

  /* Set function pointers */
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyStokesVelocityBoundaryConditionMG_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_ADAPT_MP_RESOLUTION,   (void (*)(void))ModelAdaptMaterialPointResolution_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_GeoMod2008);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_GeoMod2008);CHKERRQ(ierr);

  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,   (void (*)(void))ModelApplyInitialStokesVariableMarkers_GeoMod2008);CHKERRQ(ierr);

  /* Insert model into list */
  ierr = pTatinModelRegister(m);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

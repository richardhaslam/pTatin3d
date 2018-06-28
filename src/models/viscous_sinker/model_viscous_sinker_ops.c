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
 **    filename:   model_viscous_sinker_ops.c
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

#include "dmda_bcs.h"
#include "data_bucket.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_iterator.h"
#include "mesh_deformation.h"
#include "mesh_update.h"
#include "dmda_remesh.h"
#include "dmda_element_q2p1.h"
#include "material_point_std_utils.h"
#include "material_point_popcontrol.h"
#include "output_material_points.h"
#include "xdmf_writer.h"
#include "output_material_points_p0.h"
#include "private/quadrature_impl.h"
#include "quadrature.h"
#include "QPntSurfCoefStokes_def.h"

#include "viscous_sinker_ctx.h"

PetscErrorCode ModelInitialize_ViscousSinker(pTatinCtx c,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  RheologyConstants      *rheology;
  PetscBool flg;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  rheology                = &c->rheology_constants;
  rheology->rheology_type = RHEOLOGY_VISCOUS;

  /* box geometry */
  data->Lx = 1.0;
  data->Ly = 1.0;
  data->Lz = 1.0;

  /* inclusion geometry */
  data->is_sphere = PETSC_TRUE;

  data->origin[0] = 0.5 * data->Lx;
  data->origin[1] = 0.5 * data->Ly;
  data->origin[2] = 0.5 * data->Lz;

  /* spheriod diameter or box length */
  data->length[0] = 0.5 * data->Lx;
  data->length[1] = 0.5 * data->Ly;
  data->length[2] = 0.5 * data->Lz;

  /* bc type */
  data->boundary_conditon_type = VSBC_FreeSlipFreeSurface;

  /* parse from command line */
  rheology->const_eta0[0] = 1.0;
  rheology->const_eta0[1] = 1.0;

  rheology->const_rho0[0] = 0.0;
  rheology->const_rho0[1] = 1.0;

  rheology->nphases_active = 2;
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_eta0",&rheology->const_eta0[0],&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_eta1",&rheology->const_eta0[1],&flg);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_rho0",&rheology->const_rho0[0],&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_rho1",&rheology->const_rho0[1],&flg);CHKERRQ(ierr);

  /* set initial values for model parameters */
  /* material properties */
  data->nmaterials = rheology->nphases_active;
  data->eta0 = rheology->const_eta0[0];
  data->eta1 = rheology->const_eta0[1];
  data->rho0 = rheology->const_rho0[0];
  data->rho1 = rheology->const_rho0[1];

  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_Lx",&data->Lx,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_Ly",&data->Ly,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_Lz",&data->Lz,&flg);CHKERRQ(ierr);

  flg = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-model_viscous_sinker_cube",&data->is_sphere,&flg);CHKERRQ(ierr);
  if (flg == PETSC_TRUE) { data->is_sphere = PETSC_FALSE; }

  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_Ox",&data->origin[0],&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_Oy",&data->origin[1],&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_Oz",&data->origin[2],&flg);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_dx",&data->length[0],&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_dy",&data->length[1],&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_dz",&data->length[2],&flg);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(NULL,NULL,"-model_viscous_sinker_bc_type",(PetscInt*)&data->boundary_conditon_type,&flg);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscBool BCEvaluator_SpinXZ_u(PetscScalar coor[],PetscScalar *value,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscBool             impose_dirichlet = PETSC_TRUE;
  PetscReal             vx,vz;
  
  vx =   coor[2] - 0.5 * data->Lx;
  vz = -(coor[0] - 0.5 * data->Lz);
  vx = vx * 0.5 * (1.0 - PetscCosReal(2.0 * PETSC_PI * coor[0]/data->Lx));
  vz = vz * 0.5 * (1.0 - PetscCosReal(2.0 * PETSC_PI * coor[2]/data->Lz));
  *value = vx;
  return impose_dirichlet;
}

PetscBool BCEvaluator_SpinXZ_w(PetscScalar coor[],PetscScalar *value,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscBool             impose_dirichlet = PETSC_TRUE;
  PetscReal             vx,vz;
  
  vx =   coor[2] - 0.5 * data->Lx;
  vz = -(coor[0] - 0.5 * data->Lz);
  vx = vx * 0.5 * (1.0 - PetscCosReal(2.0 * PETSC_PI * coor[0]/data->Lx));
  vz = vz * 0.5 * (1.0 - PetscCosReal(2.0 * PETSC_PI * coor[2]/data->Lz));
  *value = vz;
  return impose_dirichlet;
}

PetscBool BCEvaluator_SpinYZ_v(PetscScalar coor[],PetscScalar *value,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscBool             impose_dirichlet = PETSC_TRUE;
  PetscReal             vy,vz;
  
  vy =   coor[2] - 0.5 * data->Ly;
  vz = -(coor[1] - 0.5 * data->Lz);
  vy = vy * 0.5 * (1.0 - PetscCosReal(2.0 * PETSC_PI * coor[1]/data->Ly));
  vz = vz * 0.5 * (1.0 - PetscCosReal(2.0 * PETSC_PI * coor[2]/data->Lz));
  *value = vy;
  return impose_dirichlet;
}

PetscBool BCEvaluator_SpinYZ_w(PetscScalar coor[],PetscScalar *value,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscBool             impose_dirichlet = PETSC_TRUE;
  PetscReal             vy,vz;
  
  vy =   coor[2] - 0.5 * data->Ly;
  vz = -(coor[1] - 0.5 * data->Lz);
  vy = vy * 0.5 * (1.0 - PetscCosReal(2.0 * PETSC_PI * coor[1]/data->Ly));
  vz = vz * 0.5 * (1.0 - PetscCosReal(2.0 * PETSC_PI * coor[2]/data->Lz));
  *value = vz;
  return impose_dirichlet;
}

PetscErrorCode ModelApplyBoundaryCondition_ViscousSinker(pTatinCtx user,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscScalar zero = 0.0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  switch (data->boundary_conditon_type) {
    case VSBC_FreeSlip:
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      break;

    case VSBC_NoSlip:
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      break;

    case VSBC_FreeSlipFreeSurface:
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      break;

    case VSBC_NoSlipFreeSurface:
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);


      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      break;

    case VSBC_TimeDependent:
    {
      BCList bclist_u;
      DM     dm_u;
      
      dm_u     = user->stokes_ctx->dav;
      bclist_u = user->stokes_ctx->u_bclist;
      ierr = BCListReset(bclist_u);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      
      ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      if (user->step < 20) {
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,0,BCEvaluator_SpinXZ_u,(void*)data);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,2,BCEvaluator_SpinXZ_w,(void*)data);CHKERRQ(ierr);

      } else {
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,1,BCEvaluator_SpinYZ_v,(void*)data);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,2,BCEvaluator_SpinYZ_w,(void*)data);CHKERRQ(ierr);
      }
    }
    break;
    
    case VSBC_Test:
      {
        PhysCompStokes stokes;
        BCList         bclist;
        DM             dav;

        ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
        bclist = stokes->u_bclist;
        ierr = PhysCompStokesGetDMs(stokes,&dav,NULL);CHKERRQ(ierr);

        /* free slip */
        /*
        // passed
        ierr = DirichletBC_FreeSlip(bclist,dav,FRONT_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist,dav,BACK_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist,dav,EAST_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist,dav,NORTH_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist,dav,SOUTH_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist,dav,WEST_FACE);CHKERRQ(ierr);
        */

        /* strain rate xx, fixed z boundaries, free slip base and free surface */
        /*
        // passed : -model_viscous_sinker_bc_type 4 -model_viscous_sinker_rho0 0.0 -model_viscous_sinker_rho1 0.0
        ierr = DirichletBC_ApplyStrainRateExx(bclist,dav,2.2);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        */

        /* strain rate zz, fixed x boundaries, free slip base and free surface */
        /*
        // passed
        ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,3.3,2);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        */

        /* shear along north/south faces + free slip base, free surface */
        /*
        // passed
        ierr = DirichletBC_ApplyStrainRateExz(bclist,dav,4.4);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        */

        /* extension */
        /*
        // passed
        ierr = DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(bclist,dav,1.0,10.0);CHKERRQ(ierr);
        */

        /* extension in x / compression in y to conserve volume + free slip base, free surface */
        // passed
        //ierr = DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ(bclist,dav,3.5);CHKERRQ(ierr);
        //ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        /* free slip, free surface, normal stress on IMAX */
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        //
        {
          SurfaceQuadrature surfQ_east;
          QPntSurfCoefStokes *surfQ_coeff,*surfQ_cell_coeff;
          PetscInt c,q,nqp,nfaces,*element_list;
          QPoint3d *qp3d;

          DM cda;
          Vec gcoords;
          PetscReal *LA_gcoords;
          PetscInt nel,nen_u;
          const PetscInt *elnidx_u;
          PetscReal elcoords[3*Q2_NODES_PER_EL_3D];

          ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
          ierr = PhysCompStokesGetDMs(stokes,&dav,NULL);CHKERRQ(ierr);

          ierr = PhysCompStokesGetSurfaceQuadrature(stokes,HEX_FACE_Pxi,&surfQ_east);CHKERRQ(ierr);
          ierr = SurfaceQuadratureGetQuadratureInfo(surfQ_east,&nqp,NULL,&qp3d);CHKERRQ(ierr);
          ierr = SurfaceQuadratureGetFaceInfo(surfQ_east,NULL,&nfaces,&element_list);CHKERRQ(ierr);

          ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
          /* setup for coords */
          ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
          ierr = DMGetCoordinatesLocal(dav,&gcoords );CHKERRQ(ierr);
          ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

          ierr = SurfaceQuadratureGetAllCellData_Stokes(surfQ_east,&surfQ_coeff);CHKERRQ(ierr);

          for (c=0; c<nfaces; c++) {
            PetscInt eidx;

            eidx = element_list[c];

            ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*eidx],LA_gcoords);CHKERRQ(ierr);

            ierr = SurfaceQuadratureGetCellData_Stokes(surfQ_east,surfQ_coeff,c,&surfQ_cell_coeff);CHKERRQ(ierr);

            for (q=0; q<nqp; q++) {
              PetscReal E[3][3],exx,exy,eyy,eta;
              PetscReal pos_qp[3],hydro_pressure_qp,y_qp;
              double *normal,*traction;

              QPntSurfCoefStokesGetField_surface_normal(&surfQ_cell_coeff[q],&normal);
              QPntSurfCoefStokesGetField_surface_traction(&surfQ_cell_coeff[q],&traction);

              //printf("normal %1.4e %1.4e %1.4e \n",normal[0],normal[1],normal[2]);

              ierr = SurfaceQuadratureInterpolate3D(surfQ_east,&qp3d[q],3,elcoords,pos_qp);CHKERRQ(ierr);
              y_qp = pos_qp[1];

              hydro_pressure_qp = (1.0 - y_qp) * 1.000 * 1.0;

              /*
                 surfQ_cell_coeff[q].traction[0] = -1.0;
                 surfQ_cell_coeff[q].traction[1] = 0.0;
                 surfQ_cell_coeff[q].traction[2] = 0.0;
              */
              /*
                 surfQ_cell_coeff[q].traction[0] = -(1.0 - y_qp) * 1.000 * 1.0;
                 surfQ_cell_coeff[q].traction[1] = 0.0;
                 surfQ_cell_coeff[q].traction[2] = 0.0;
              */

              exx = 0.0;
              exy = 0.0;
              eyy = 0.0;
              E[0][0] = exx; E[0][1] = exy; E[0][2] = 0.0;
              E[1][0] = exy; E[1][1] = eyy; E[1][2] = 0.0;
              E[2][0] = 0.0; E[2][1] = 0.0; E[2][2] = 0.0;

              eta = 1.0;
              traction[0] = (2.0 * eta * E[0][0] - hydro_pressure_qp) * (normal[0]);
              traction[1] = (2.0 * eta * E[1][0]                    ) * (normal[0]);
              traction[2] = (2.0 * eta * E[2][0]                    ) * (normal[0]);
            }
          }

          ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
        }
        //
        //ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

      }

      break;
  }

  /*
     {
     BCList flat;

     ierr = BCListFlattenedCreate(user->stokes_ctx->u_bclist,&flat);CHKERRQ(ierr);
     ierr = BCListDestroy(&user->stokes_ctx->u_bclist);CHKERRQ(ierr);
     user->stokes_ctx->u_bclist = flat;
     }
  */

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyBoundaryConditionMG_ViscousSinker(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscScalar zero = 0.0;
  PetscInt n;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  for (n=0; n<nl; n++) {

    switch (data->boundary_conditon_type) {
      case VSBC_FreeSlip:
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        break;

      case VSBC_NoSlip:
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        break;

      case VSBC_FreeSlipFreeSurface:
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        break;

      case VSBC_NoSlipFreeSurface:
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        break;

      case VSBC_TimeDependent:
      {
        BCList bclist_u;
        DM     dm_u;
        
        dm_u     = dav[n];
        bclist_u = bclist[n];
        ierr = BCListReset(bclist_u);CHKERRQ(ierr);
        /*
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,0,BCEvaluator_SpinXZ_u,(void*)data);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,2,BCEvaluator_SpinXZ_w,(void*)data);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        */
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        
        if (user->step < 20) {
          ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
          
          ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,0,BCEvaluator_SpinXZ_u,(void*)data);CHKERRQ(ierr);
          ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
          ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,2,BCEvaluator_SpinXZ_w,(void*)data);CHKERRQ(ierr);
          
        } else {
          ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
          
          ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
          ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,1,BCEvaluator_SpinYZ_v,(void*)data);CHKERRQ(ierr);
          ierr = DMDABCListTraverse3d(bclist_u,dm_u,DMDABCList_IMAX_LOC,2,BCEvaluator_SpinYZ_w,(void*)data);CHKERRQ(ierr);
        }
        
      }
      break;

      case VSBC_Test:
        /* free slip */
        /*
        // passed
        ierr = DirichletBC_FreeSlip(bclist[n],dav[n],FRONT_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist[n],dav[n],BACK_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist[n],dav[n],EAST_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist[n],dav[n],NORTH_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist[n],dav[n],SOUTH_FACE);CHKERRQ(ierr);
        ierr = DirichletBC_FreeSlip(bclist[n],dav[n],WEST_FACE);CHKERRQ(ierr);
        */

        /* strain rate xx */
        /*
        // passed
        ierr = DirichletBC_ApplyStrainRateExx(bclist[n],dav[n],2.2);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        */

        /* strain rate zz, fixed x boundaries, free slip base and free surface */
        /*
        // passed
        ierr = DirichletBC_ApplyDirectStrainRate(bclist[n],dav[n],3.3,2);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        */

        /* shear along north/south faces + free slip base, free surface */
        /*
        // passed
        ierr = DirichletBC_ApplyStrainRateExz(bclist[n],dav[n],4.4);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        */

        /* extension */
        /*
        // passed
        ierr = DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(bclist[n],dav[n],1.0,10.0);CHKERRQ(ierr);
        */

        /* extension in x / compression in y to conserve volume + free slip base, free surface */
        // passed
        //ierr = DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ(bclist[n],dav[n],3.5);CHKERRQ(ierr);
        //ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        //ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist[n],dav[n],DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        break;
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyMaterialBoundaryCondition_ViscousSinker(pTatinCtx c,void *ctx)
{
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  PetscPrintf(PETSC_COMM_WORLD,"  NOT IMPLEMENTED \n", PETSC_FUNCTION_NAME);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMeshGeometry_ViscousSinker(pTatinCtx c,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav,0.0,data->Lx, 0.0,data->Ly, 0.0,data->Lz);CHKERRQ(ierr);

  //
  PetscPrintf(PETSC_COMM_WORLD,"RUNNING DEFORMED MESH EXAMPLE \n");
  ierr = MeshDeformation_GaussianBump_YMAX(c->stokes_ctx->dav,-0.3,-5.6);CHKERRQ(ierr);
  //ierr = DMDASetGraduatedCoordinates1D(c->stokes_ctx->dav,2,1,2.0);CHKERRQ(ierr);
  // [test c] remesh interp
  //ierr = DMDASetCoordinatesCentralSqueeze1D(c->stokes_ctx->dav,0,4.0,0.0,0.4,0.8,1.0);CHKERRQ(ierr);
  //ierr = DMDASetCoordinatesColumnRefinement(c->stokes_ctx->dav,1,4.0,0.66,1.0);CHKERRQ(ierr);

  /* diffusion example */
  /*
     ierr = UpdateMeshGeometry_ApplyDiffusionJMAX(c->stokes_ctx->dav,1.0e-2,0.44,
     PETSC_TRUE,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE, PETSC_FALSE);CHKERRQ(ierr);
  */
  // [test d] remesh vertically and preserve topography
  /*
     {
     PetscInt npoints,dir;
     PetscReal xref[10],xnat[10];

     npoints = 6;
     xref[0] = 0.0;
     xref[1] = 0.2;
     xref[2] = 0.4;
     xref[3] = 0.6;
     xref[4] = 0.8;
     xref[5] = 1.0;

     xnat[0] = 0.0;
     xnat[1] = 0.67;
     xnat[2] = 0.92;
     xnat[3] = 0.97;
     xnat[4] = 0.985;
     xnat[5] = 1.0;

     dir = 1;
     ierr = DMDACoordinateRefinementTransferFunction(c->stokes_ctx->dav,dir,PETSC_TRUE,npoints,xref,xnat);CHKERRQ(ierr);
     ierr = DMDABilinearizeQ2Elements(c->stokes_ctx->dav);CHKERRQ(ierr);
     }
  */

  PetscFunctionReturn(0);
}

PetscErrorCode ViscousSinker_ApplyInitialMaterialGeometry_SingleInclusion(pTatinCtx c,ModelViscousSinkerCtx *data)
{
  int                    p,n_mp_points;
  DataBucket             db;
  DataField              PField_std,PField_stokes;
  int                    phase;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* define properties on material points */
  db = c->materialpoint_db;
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  DataFieldVerifyAccess(PField_std,sizeof(MPntStd));

  DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
  DataFieldGetAccess(PField_stokes);
  DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));

  DataBucketGetSizes(db,&n_mp_points,0,0);

  for (p=0; p<n_mp_points; p++) {
    MPntStd     *material_point;
    MPntPStokes *mpprop_stokes;
    double      *position;
    double      eta,rho;

    DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
    DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);

    /* Access using the getter function provided for you (recommeneded for beginner user) */
    MPntStdGetField_global_coord(material_point,&position);

    phase = 0;
    eta =  data->eta0;
    rho = -data->rho0;

    if (data->is_sphere) {
      double rx = (position[0]-data->origin[0])/(0.5*data->length[0]);
      double ry = (position[1]-data->origin[1])/(0.5*data->length[1]);
      double rz = (position[2]-data->origin[2])/(0.5*data->length[2]);

      if ( rx*rx + ry*ry + rz*rz < 1.0 ) {
        phase = 1;
        eta =  data->eta1;
        rho = -data->rho1;
      }

    } else { /* box */
      /*
         printf("length = %lf %lf %lf \n",data->length[0],data->length[1],data->length[2]);
         printf("origin = %lf %lf %lf \n",data->origin[0],data->origin[1],data->origin[2]);
         printf("x: range [%lf %lf] xp = %lf \n",data->origin[0] - 0.5*data->length[0],data->origin[0] + 0.5*data->length[0],position[0]);
         printf("y: range [%lf %lf] yp = %lf \n",data->origin[1] - 0.5*data->length[1],data->origin[1] + 0.5*data->length[1],position[1]);
         printf("z: range [%lf %lf] zp = %lf \n",data->origin[2] - 0.5*data->length[2],data->origin[2] + 0.5*data->length[2],position[2]);
         */
      if ( (position[0]>data->origin[0] - 0.5*data->length[0]) && (position[0]<data->origin[0] + 0.5*data->length[0]) ) {
        if ( (position[1]>data->origin[1] - 0.5*data->length[1]) && (position[1]<data->origin[1] + 0.5*data->length[1]) ) {
          if ( (position[2]>data->origin[2] - 0.5*data->length[2]) && (position[2]<data->origin[2] + 0.5*data->length[2]) ) {
            phase = 1;
            eta =  data->eta1;
            rho = -data->rho1;
          }
        }
      }
    }

    /* user the setters provided for you */
    MPntStdSetField_phase_index(material_point,phase);

    MPntPStokesSetField_eta_effective(mpprop_stokes,eta);
    MPntPStokesSetField_density(mpprop_stokes,rho);
  }

  DataFieldRestoreAccess(PField_std);
  DataFieldRestoreAccess(PField_stokes);

  PetscFunctionReturn(0);
}

PetscErrorCode compute_inclusion_origins(PetscInt ninclusions,PetscReal rmax,PetscReal Lx,PetscReal Ly,PetscReal Lz,
    PetscReal **_pos)
{
  PetscReal      *pos;
  PetscInt       p,found=0,overlap,attempt,loops=0;

  PetscFunctionBegin;

  PetscMalloc(sizeof(PetscReal)*3*ninclusions,&pos);

  srand(0);

  loops = 0;
  PetscPrintf(PETSC_COMM_WORLD,"  Commencing inclusion generation \n");
START_INCLUSION:

  loops++;
  found = 0;
  attempt = 0;
  while (found < ninclusions) {
    PetscReal xp = rand()/( (PetscReal)RAND_MAX );
    PetscReal yp = rand()/( (PetscReal)RAND_MAX );
    PetscReal zp = rand()/( (PetscReal)RAND_MAX );
    PetscReal dx,dy,dz,range[3],cp[3];

    if (attempt == 50000) { goto START_INCLUSION; }

    //xp = 2.1*rmax + xp * (Lx-2.1*rmax);
    //yp = 2.1*rmax + yp * (Ly-2.1*rmax);
    //zp = 2.1*rmax + zp * (Lz-2.1*rmax);

    xp = xp * (Lx);
    yp = yp * (Ly);
    zp = zp * (Lz);
    attempt++;
    //
    dx = 1.5*rmax;
    range[0] = xp - dx;
    if (range[0] < 0.0) { continue; }
    range[0] = xp + dx;
    if (range[0] > Lx) { continue; }

    dy = 1.5*rmax;
    range[1] = yp - dy;
    if (range[1] < 0.0) { continue; }
    range[1] = yp + dy;
    if (range[1] > Ly) { continue; }

    dz = 1.5*rmax;
    range[2] = zp - dz;
    if (range[2] < 0.0) { continue; }
    range[2] = zp + dz;
    if (range[2] > Lz) { continue; }
    //
    /* check others for overlap */
    cp[0] = xp;
    cp[1] = yp;
    cp[2] = zp;
    overlap = 0;
    for (p=0; p<found; p++) {
      PetscScalar sep;

      sep = PetscSqrtReal(
          (pos[3*p+0]-cp[0])*(pos[3*p+0]-cp[0])
          + (pos[3*p+1]-cp[1])*(pos[3*p+1]-cp[1])
          + (pos[3*p+2]-cp[2])*(pos[3*p+2]-cp[2]) );

      if (sep < 2.1*rmax) {
        overlap = 1;
        break;
      }
    }
    if (overlap == 1) { continue; }

    pos[3*found+0] = xp;
    pos[3*found+1] = yp;
    pos[3*found+2] = zp;
    found++;
  }
  PetscPrintf(PETSC_COMM_WORLD,"  inclusion generation performed %D loops: Made %D attempts and correctly defined %D of %D inclusions\n",loops,attempt,found,ninclusions);

  *_pos   = pos;

  PetscFunctionReturn(0);
}


PetscErrorCode ViscousSinker_ApplyInitialMaterialGeometry_MultipleInclusions(pTatinCtx c,ModelViscousSinkerCtx *data,PetscInt ninclusions)
{
  int            p,n_mp_points;
  DataBucket     db;
  DataField      PField_std,PField_stokes;
  int            phase;
  PetscInt       cc;
  PetscReal      max_radius,*inclusion_pos;
  PetscBool      inclusion_view = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  PetscOptionsGetBool(NULL,NULL,"-model_viscous_sinker_inclusion_view",&inclusion_view,NULL);

  /* define properties on material points */
  ierr = pTatinGetMaterialPoints(c,&db,NULL);CHKERRQ(ierr);
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  DataFieldVerifyAccess(PField_std,sizeof(MPntStd));

  DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
  DataFieldGetAccess(PField_stokes);
  DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));

  DataBucketGetSizes(db,&n_mp_points,0,0);

  //  max_radius = 0.5*data->length[0];
  //  if ( 0.5*data->length[1] > max_radius ) { max_radius = 0.5*data->length[0]; }
  //  if ( 0.5*data->length[2] > max_radius ) { max_radius = 0.5*data->length[2]; }
  max_radius = 0.25*data->length[0]*data->length[0]
    + 0.25*data->length[1]*data->length[1]
    + 0.25*data->length[2]*data->length[2];
  max_radius = PetscSqrtReal(max_radius);

  max_radius = PetscMax(0.5*data->length[0],0.5*data->length[1]);
  max_radius = PetscMax(max_radius,0.5*data->length[2]);

  ierr = compute_inclusion_origins(ninclusions,max_radius,data->Lx,data->Ly,data->Lz,&inclusion_pos);CHKERRQ(ierr);
  if (inclusion_view) {
    for (cc=0; cc<ninclusions; cc++) {
      PetscReal *cp = &inclusion_pos[3*cc];
      PetscPrintf(PETSC_COMM_WORLD," inclusion[%d]: Ox = %1.4e %1.4e %1.4e \n",cc,cp[0],cp[1],cp[2]);
    }
  }

  for (p=0; p<n_mp_points; p++) {
    MPntStd     *material_point;
    MPntPStokes *mpprop_stokes;
    double      *position;
    double      eta,rho;
    PetscBool   inside_inclusion;

    DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
    DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);

    /* Access using the getter function provided for you (recommeneded for beginner user) */
    MPntStdGetField_global_coord(material_point,&position);

    phase = 0;
    eta =  data->eta0;
    rho = -data->rho0;

    inside_inclusion = PETSC_FALSE;

    if (data->is_sphere) {
      for (cc=0; cc<ninclusions; cc++) {
        PetscReal *cp = &inclusion_pos[3*cc];
        PetscReal sep,rx,ry,rz;

        rx = (position[0]-cp[0])/(0.5*data->length[0]);
        ry = (position[1]-cp[1])/(0.5*data->length[1]);
        rz = (position[2]-cp[2])/(0.5*data->length[2]);

        sep = rx*rx + ry*ry + rz*rz;
        if (sep < 1.0) {
          inside_inclusion = PETSC_TRUE;
          break;
        }

      }
    } else { /* box */

      for (cc=0; cc<ninclusions; cc++) {
        PetscReal *cp = &inclusion_pos[3*cc];

        if ( (position[0]>cp[0] - 0.5*data->length[0]) && (position[0]<cp[0] + 0.5*data->length[0]) ) {
          if ( (position[1]>cp[1] - 0.5*data->length[1]) && (position[1]<cp[1] + 0.5*data->length[1]) ) {
            if ( (position[2]>cp[2] - 0.5*data->length[2]) && (position[2]<cp[2] + 0.5*data->length[2]) ) {
              inside_inclusion = PETSC_TRUE;
              break;
            }
          }
        }

      }
    }

    if (inside_inclusion) {
      phase = 1;
      eta =  data->eta1;
      rho = -data->rho1;
    }

    /* user the setters provided for you */
    MPntStdSetField_phase_index(material_point,phase);

    MPntPStokesSetField_eta_effective(mpprop_stokes,eta);
    MPntPStokesSetField_density(mpprop_stokes,rho);
  }

  DataFieldRestoreAccess(PField_std);
  DataFieldRestoreAccess(PField_stokes);

  PetscFree(inclusion_pos);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_ViscousSinker(pTatinCtx c,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscInt              ninclusions = 1;
  PetscErrorCode        ierr;

  PetscOptionsGetInt(NULL,NULL,"-model_viscous_sinker_ninclusions",&ninclusions,0);
  if (ninclusions == 1) {
    ierr = ViscousSinker_ApplyInitialMaterialGeometry_SingleInclusion(c,data);CHKERRQ(ierr);
  } else {
    ierr = ViscousSinker_ApplyInitialMaterialGeometry_MultipleInclusions(c,data,ninclusions);CHKERRQ(ierr);
  }

#if 0
  /* Test: Face marker insertion based on tagged elements which live on the north face of the DMDA */
  {
    PetscInt Nxp[] = { 20, 20 };
    PetscInt start_pidx,n_pidx,face_idx,p;
    DataBucket db;
    DM dav;
    PetscInt ncells_list,*cell_list;
    PetscInt e,si,sj,sk,lnx,lny,lnz,M,N,P,nel,nen,elnidx,lmx,lmy,lmz,ei,ej,ek;
    PetscBool contains_north_face;
    MPAccess mpX;
    PetscInt ncc,npp,*pcell_list;
    PSortCtx *plist;

    face_idx = 2;
    dav = c->stokes_ctx->dav;
    db = c->materialpoint_db;

    /* Identify if the sub-domain of this DMDA is connected to the north face of the mesh */
    ierr = DMDAGetCorners(dav,&si,&sj,&sk,&lnx,&lny,&lnz);CHKERRQ(ierr);
    ierr = DMDAGetInfo(dav,0, &M,&N,&P, 0,0,0, 0,0, 0,0,0, 0);CHKERRQ(ierr);
    contains_north_face = PETSC_FALSE; if (sj+lny == N) { contains_north_face = PETSC_TRUE; }

    ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
    ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen,&elnidx);CHKERRQ(ierr);

    //
    ncells_list = 0;
    cell_list = NULL;

    /* generate cell list, identifying only those living on the surface */
    if (contains_north_face) {
      PetscInt ecnt;

      ncells_list = lmx*lmy*lmz;
      PetscMalloc(sizeof(PetscInt)*lmx*lmy*lmz,&cell_list);

      // init //
      for (e=0; e<lmx*lmy*lmz; e++) { cell_list[e] = -1; }

      // label each cell on the north surface //
      ej = lmy-1;
      for (ek=0; ek<lmz; ek++) {
        for (ei=0; ei<lmx; ei++) {
          ecnt = ei + ej * lmx + ek * lmx * lmy;

          cell_list[ecnt] = ecnt;
        }
      }
    }
    //

    /* mask out empty cells */
    ierr = MPPCCreateSortedCtx(db,dav,&npp,&ncc,&plist,&pcell_list);CHKERRQ(ierr);

    if (contains_north_face) {
      ej = lmy-1;
      for (ek=0; ek<lmz; ek++) {
        for (ei=0; ei<lmx; ei++) {
          PetscInt nppercell,p;
          PetscInt eidx;
          MPntStd *point;

          eidx = ei + ej * lmx + ek * lmx * lmy;
          ierr = MPPCSortedCtxGetNumberOfPointsPerCell(db,eidx,pcell_list,&nppercell);CHKERRQ(ierr);

          //PetscPrintf(PETSC_COMM_SELF,"cell %d : nppercell %d \n",eidx,nppercell);
          for (p=0; p<nppercell; p++) {
            ierr = MPPCSortedCtxGetPointByCell(db,eidx,p,plist,pcell_list,&point);CHKERRQ(ierr);

            if (point->phase != 1) {
              cell_list[eidx] = -1;
            }
          }

        }
      }
    }

    ierr = MPPCDestroySortedCtx(db,dav,&plist,&pcell_list);CHKERRQ(ierr);

    /* insert new markers */
    ierr = SwarmMPntStd_CoordAssignmentFromElementList_FaceLatticeLayout3d(dav,Nxp,0.0,face_idx,ncells_list,cell_list,&start_pidx,&n_pidx,db);CHKERRQ(ierr);

    //PetscPrintf(PETSC_COMM_SELF,"start of new points %d \n",start_pidx);
    //PetscPrintf(PETSC_COMM_SELF,"num   of new points %d \n",n_pidx);

    /* set phase on all newly added face markers */
    ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
    for (p=start_pidx; p<start_pidx+n_pidx; p++) {
      ierr = MaterialPointSet_phase_index(mpX,p,5);CHKERRQ(ierr);
    }
    ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);

    if (contains_north_face) {
      PetscFree(cell_list);
    }
  }
#endif

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyUpdateMeshGeometry_ViscousSinker(pTatinCtx c,Vec X,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  Vec velocity,pressure;
  DM stokes_pack,dav,dap;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  stokes_pack = c->stokes_ctx->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  if (data->boundary_conditon_type != (PetscInt)VSBC_TimeDependent) {
    //ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,c->dt);CHKERRQ(ierr);
    ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,c->dt);CHKERRQ(ierr);
  }

  // [test c] remesh interp
  //ierr = UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX(dav,velocity,NULL,c->dt);CHKERRQ(ierr);
  //ierr = DMDASetCoordinatesColumnRefinement(dav,1,4.0,0.75,1.0);CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelOutput_ViscousSinker(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
  // testing
  //ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(c,X,prefix);CHKERRQ(ierr);
  //ierr = pTatinOutputLiteMeshVelocitySlicedPVTS(c->stokes_ctx->stokes_pack,c->outputpath,prefix);CHKERRQ(ierr);
  //ierr = ptatin3d_StokesOutput_VelocityXDMF(c,X,prefix);CHKERRQ(ierr);
  // testing
  //ierr = pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
  // tests for alternate (output/load)ing of "single file marker" formats
  //ierr = SwarmDataWriteToPetscVec(c->materialpoint_db,prefix);CHKERRQ(ierr);
  //ierr = SwarmDataLoadFromPetscVec(c->materialpoint_db,prefix);CHKERRQ(ierr);
  /*
     {
     MaterialPointVariable vars[] = { MPV_viscosity, MPV_density, MPV_region };
  //ierr = pTatin3dModelOutput_MarkerCellFieldsP0_ParaView(c,sizeof(vars)/sizeof(MaterialPointVariable),vars,PETSC_TRUE,prefix);CHKERRQ(ierr);
  ierr = pTatin3dModelOutput_MarkerCellFieldsP0_PetscVec(c,PETSC_TRUE,sizeof(vars)/sizeof(MaterialPointVariable),vars,prefix);CHKERRQ(ierr);
  }
  */
  PetscFunctionReturn(0);
}

PetscErrorCode ModelInitialCondition_ViscousSinker(pTatinCtx c,Vec X,void *ctx)
{
  DM stokes_pack,dau,dap;
  Vec velocity,pressure;
  //DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
  //DMDAVecTraverse3d_InterpCtx IntpCtx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  stokes_pack = c->stokes_ctx->stokes_pack;

  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
  /* apply -5 < vx 5 across the domain x \in [0,1] */
  //ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,10.0,-5.0,0.0);CHKERRQ(ierr);
  //ierr = DMDAVecTraverse3d(dau,velocity,0,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);

  ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  /*
     ierr = PetscOptionsGetReal(NULL,NULL,"-model_viscous_sinker_rho0",&rho0,0);CHKERRQ(ierr);

     HPctx.surface_pressure = 0.0;
     HPctx.ref_height = data->Ly;
     HPctx.ref_N = c->stokes_ctx->my-1;
     HPctx.grav = 1.0;
     HPctx.rho = rho0;

     ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
     ierr = DMDAVecTraverseIJK(dap,pressure,0,DMDAVecTraverseIJK_HydroStaticPressure_v2,(void*)&HPctx);
     ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

     ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,"testHP");CHKERRQ(ierr);
     */

  PetscFunctionReturn(0);
}

PetscErrorCode ModelDestroy_ViscousSinker(pTatinCtx c,void *ctx)
{
  ModelViscousSinkerCtx *data = (ModelViscousSinkerCtx*)ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  /* Free contents of structure */

  /* Free structure */
  ierr = PetscFree(data);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinModelRegister_ViscousSinker(void)
{
  ModelViscousSinkerCtx *data;
  pTatinModel m;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* Allocate memory for the data structure for this model */
  ierr = PetscMalloc(sizeof(ModelViscousSinkerCtx),&data);CHKERRQ(ierr);
  ierr = PetscMemzero(data,sizeof(ModelViscousSinkerCtx));CHKERRQ(ierr);

  /* register user model */
  ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

  /* Set name, model select via -ptatin_model NAME */
  ierr = pTatinModelSetName(m,"viscous_sinker");CHKERRQ(ierr);

  /* Set model data */
  ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);

  /* Set function pointers */
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelInitialCondition_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_ViscousSinker);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_ViscousSinker);CHKERRQ(ierr);

  /* Insert model into list */
  ierr = pTatinModelRegister(m);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

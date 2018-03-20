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
 **    filename:   ptatin3d_energy.c
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

#include "ptatin_utils.h"
#include "dmda_bcs.h"
#include "dmda_duplicate.h"
#include "dmdae.h"
#include "element_utils_q1.h"
#include "dmda_element_q1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"
#include "material_point_utils.h"

#include "MPntPEnergy_def.h"
#include "QPntVolCoefEnergy_def.h"
#include "phys_comp_energy.h"
#include "energy_assembly.h"
#include "ptatin3d_energy.h"


PetscErrorCode pTatinGetContext_Energy(pTatinCtx ctx,PhysCompEnergy *e)
{
  PetscFunctionBegin;
  if (e) { *e = ctx->energy_ctx; }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinContextValid_Energy(pTatinCtx ctx,PetscBool *exists)
{
  PetscFunctionBegin;
  *exists = PETSC_FALSE;
  if (ctx->energy_ctx) {
    *exists = PETSC_TRUE;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompCreate_Energy(pTatinCtx user)
{
  PetscErrorCode ierr;
  PhysCompStokes stokes_ctx;
  PetscInt energy_mesh_type;

  PetscFunctionBegin;
  stokes_ctx = user->stokes_ctx;
  /* create from data */

  energy_mesh_type = 1; /* default is Q1 overlapping Q2 */
  ierr = PetscOptionsGetInt(NULL,NULL,"-energy_mesh_type",&energy_mesh_type,0);CHKERRQ(ierr);
  ierr = PhysCompNew_Energy(stokes_ctx->dav,-1,-1,-1,energy_mesh_type,&user->energy_ctx);CHKERRQ(ierr);

  if (user->restart_from_file) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatinPhysCompCreate_Energy should not be called during restart");
  } else {
    ierr = PhysCompAddMaterialPointCoefficients_Energy(user->materialpoint_db);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompActivate_Energy(pTatinCtx user,PetscBool load)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (load && (user->energy_ctx == NULL)) {
    ierr = pTatinPhysCompCreate_Energy(user);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompAttachData_Energy(pTatinCtx user,Vec T,Mat A)
{
  PhysCompEnergy e;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = pTatinGetContext_Energy(user,&e);CHKERRQ(ierr);

  if (T) {
    ierr = pTatinCtxAttachModelData(user,"PhysCompEnergy_T",(void*)T);CHKERRQ(ierr);
  }
  if (A) {
    ierr = pTatinCtxAttachModelData(user,"PhysCompEnergy_JE",(void*)A);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompGetData_Energy(pTatinCtx user,Vec *T,Mat *A)
{
  PhysCompEnergy e;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = pTatinGetContext_Energy(user,&e);CHKERRQ(ierr);

  if (T) {
    ierr = pTatinCtxGetModelData(user,"PhysCompEnergy_T",(void**)T);CHKERRQ(ierr);
  }
  if (A) {
    ierr = pTatinCtxGetModelData(user,"PhysCompEnergy_JE",(void**)A);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode MaterialPointQuadraturePointProjectionC0_Q2Energy(DM da,DataBucket materialpoint_db,MaterialPointField field,const int member,Quadrature Q)
{
  DMDAE          dae,dae_clone;
  PetscInt       dof;
  DM             clone;
  Vec            properties_A,properties_B;
  int            npoints;
  DataField      PField_std;
  DataField      PField_material_point_property;
  MPntStd        *mp_std;
  void           *material_point_property;
  size_t         mp_field_offset, mp_offset, qp_field_offset, qp_offset;
  size_t         mp_property_offsets[MPntPEnergy_nmembers];
  size_t         qp_property_offsets[QPntVolCoefEnergy_nmembers];
  QPntVolCoefEnergy *all_quadpoints;
  PetscBool      view;
  PetscInt       nel,nen;
  const PetscInt *els;
  PetscErrorCode ierr;

  PetscFunctionBegin;


  if (field != MPField_Energy) {
    /* error - these is only valid for energy fields defined on Q2 */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"User must choose either properties from MPntPEnergy which are to be projected onto a Q2 space");
  }

  DataBucketGetDataFieldByName(materialpoint_db, MPntStd_classname,&PField_std);
  DataBucketGetSizes(materialpoint_db,&npoints,NULL,NULL);
  mp_std  = PField_std->data;


  ierr = MPntPEnergyComputeMemberOffsets(mp_property_offsets);CHKERRQ(ierr);
  ierr = QPntVolCoefEnergyComputeMemberOffsets(qp_property_offsets);CHKERRQ(ierr);

  /* setup */
  dof = 1;
  ierr = DMDADuplicateLayout(da,dof,1,DMDA_STENCIL_BOX,&clone);CHKERRQ(ierr);
  ierr = DMGetDMDAE(da,&dae);CHKERRQ(ierr);

  ierr = DMAttachDMDAE(clone);CHKERRQ(ierr);
  ierr = DMGetDMDAE(clone,&dae_clone);CHKERRQ(ierr);
  {
    PetscInt NP[3];

    ierr = DMDAGetInfo(da,0,0,0,0,&NP[0],&NP[1],&NP[2],0,0, 0,0,0, 0);CHKERRQ(ierr);
    ierr = DMDAEDeepCopy(dae,NP,dae_clone);CHKERRQ(ierr);
  }
  //ierr = DMDAECopy(dae,dae_clone);CHKERRQ(ierr);

  ierr = DMDASetElementType_Q1(clone);CHKERRQ(ierr);
  ierr = DMDAGetElements_DA_Q1_3D(clone,&nel,&nen,&els);CHKERRQ(ierr);


  ierr = DMGetGlobalVector(clone,&properties_A);CHKERRQ(ierr);
  ierr = DMGetGlobalVector(clone,&properties_B);CHKERRQ(ierr);

  ierr = VecZeroEntries(properties_A);CHKERRQ(ierr);
  ierr = VecZeroEntries(properties_B);CHKERRQ(ierr);


  switch (field) {

    case MPField_Energy:
    {
      MPntPEnergyTypeName member_name = (MPntPEnergyTypeName)member;

      mp_offset = sizeof(MPntPEnergy);
      qp_offset = sizeof(QPntVolCoefEnergy);

      DataBucketGetDataFieldByName(materialpoint_db, MPntPEnergy_classname,&PField_material_point_property);
      material_point_property = PField_material_point_property->data;

      switch (member_name) {
        case MPPEgy_diffusivity:
          ierr = PetscObjectSetName( (PetscObject)properties_A, "kappa");CHKERRQ(ierr);
          mp_field_offset = mp_property_offsets[ MPPEgy_diffusivity ];
          qp_field_offset = qp_property_offsets[ QPVCEgy_diffusivity ];
          break;
          /* ----------------------------------- */
        case MPPEgy_heat_source:
          ierr = PetscObjectSetName( (PetscObject)properties_A, "H");CHKERRQ(ierr);
          mp_field_offset = mp_property_offsets[ MPPEgy_heat_source ];
          qp_field_offset = qp_property_offsets[ QPVCEgy_heat_source ];
          break;
          /* ----------------------------------- */
        default:
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"User must choose either {MPPEgy_diffusivity, MPPEgy_heat_source}");
          break;
      }
    }
      break;

    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"User must choose either {MPntPEnergy}");
      break;
  }

  /* compute */
  //
  ierr = DMDAEQ1_MaterialPointProjection_MapOntoQ2Mesh(
                                                clone,properties_A,properties_B,
                                                //CoefAvgHARMONIC,
                                                CoefAvgARITHMETIC,
                                                npoints,mp_std,
                                                mp_field_offset,mp_offset,material_point_property);CHKERRQ(ierr);
  //

  /*
   ierr = _MaterialPointProjection_MapOntoNestedQ1Mesh(
   clone,properties_A,properties_B,
   //CoefAvgHARMONIC,
   CoefAvgARITHMETIC,
   npoints,mp_std,
   mp_field_offset,mp_offset,material_point_property);CHKERRQ(ierr);
   */

  /* interpolate to quad points */
  ierr = VolumeQuadratureGetAllCellData_Energy(Q,&all_quadpoints);CHKERRQ(ierr);
  ierr = DMDAEQ1_MaterialPointProjection_MapOntoQ2Mesh_InterpolateToQuadraturePoint(
                        clone,properties_A,
                        qp_field_offset,qp_offset,(void*)all_quadpoints,Q);CHKERRQ(ierr);


  /* view */
  view = PETSC_FALSE;
  PetscOptionsGetBool(NULL,NULL,"-view_projected_marker_fields",&view,NULL);
  if (view) {
    char filename[256];
    PetscViewer viewer;

    sprintf(filename,"MaterialPointProjection_energy_member_%d.vtk",(int)member );
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
    ierr = DMView(clone, viewer);CHKERRQ(ierr);
    ierr = VecView(properties_A, viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* destroy */
  ierr = DMRestoreGlobalVector(clone,&properties_B);CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(clone,&properties_A);CHKERRQ(ierr);

  ierr = DMDestroyDMDAE(clone);CHKERRQ(ierr);
  ierr = DMDestroy(&clone);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompEnergy_MPProjectionQ1(pTatinCtx ctx)
{
  PhysCompEnergy energy;
  DM             daT;
  DataBucket     materialpoint_db;
  Quadrature     volQ;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = pTatinGetContext_Energy(ctx,&energy);CHKERRQ(ierr);
  daT  = energy->daT;
  volQ = energy->volQ;
  ierr = pTatinGetMaterialPoints(ctx,&materialpoint_db,NULL);CHKERRQ(ierr);

  ierr = MaterialPointQuadraturePointProjectionC0_Q2Energy(daT,materialpoint_db,MPField_Energy,MPPEgy_diffusivity,volQ);CHKERRQ(ierr);
  ierr = MaterialPointQuadraturePointProjectionC0_Q2Energy(daT,materialpoint_db,MPField_Energy,MPPEgy_heat_source,volQ);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/* update V */
/*
 u - V = u - (X_current - X_old)/dt
       = (dt.u - X_current + X_old)/dt
*/
PetscErrorCode _pTatinPhysCompEnergy_UpdateALEVelocity(PhysCompEnergy energy,PetscReal dt)
{
  Vec            coordinates;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = VecScale(energy->u_minus_V,dt);CHKERRQ(ierr);
  ierr = DMGetCoordinates(energy->daT,&coordinates);CHKERRQ(ierr);
  ierr = VecAXPY(energy->u_minus_V,-1.0,coordinates);CHKERRQ(ierr);
  ierr = VecAXPY(energy->u_minus_V, 1.0,energy->Xold);CHKERRQ(ierr);
  ierr = VecScale(energy->u_minus_V,1.0/dt);CHKERRQ(ierr);

  /*
  {
    PetscReal      min,max;

    ierr = VecMin(energy->u_minus_V,0,&min);CHKERRQ(ierr);
    ierr = VecMax(energy->u_minus_V,0,&max);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"ALE(vel) min = %1.4e : max = %1.4e \n", min,max);
  }
  */
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompEnergy_UpdateALEVelocity(PhysCompStokes s,Vec X,PhysCompEnergy energy,PetscReal dt)
{
  DM             cdaT;
  Vec            velocity,pressure;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = DMGetCoordinateDM(energy->daT,&cdaT);CHKERRQ(ierr);

  /* Project fluid velocity from Q2 space into Q1 space */
  ierr = DMCompositeGetAccess(s->stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
  ierr = DMDAProjectVectorQ2toQ1(s->dav,velocity,cdaT,energy->u_minus_V,energy->energy_mesh_type);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(s->stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

  /* Compute ALE velocity in Q1 space */
  ierr = _pTatinPhysCompEnergy_UpdateALEVelocity(energy,dt);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode pTatinPhysCompEnergy_Update(PhysCompEnergy e,DM dav,Vec T)
{
  Vec            coords;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* save current coords before advecting */
  ierr = DMGetCoordinates(e->daT,&coords);CHKERRQ(ierr);
  ierr = VecCopy(coords,e->Xold);CHKERRQ(ierr);
  //ierr = DMDAUpdateGhostedCoordinates(daq1);CHKERRQ(ierr);

  /* update solution */
  ierr = VecCopy(T,e->Told);CHKERRQ(ierr);

  /* update coords */
  ierr = DMDAProjectCoordinatesQ2toQ1(dav,e->daT,e->energy_mesh_type);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompEnergy_Initialise(PhysCompEnergy e,Vec T)
{
  Vec            coords;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* update solution */
  ierr = VecCopy(T,e->Told);CHKERRQ(ierr);

  /* update coordinates */
  ierr = DMGetCoordinates(e->daT,&coords);CHKERRQ(ierr);
  ierr = VecCopy(coords,e->Xold);CHKERRQ(ierr);

  //ierr = VecZeroEntries(T);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompEnergy_ComputeTimestep(PhysCompEnergy energy,Vec X,PetscReal *timestep)
{
  Vec         philoc,Vloc,coordsloc;
  DM          da,cda;
  PetscScalar *LA_philoc,*LA_Vloc,*LA_coordsloc;
  PetscReal   el_volume,el_kappa_const,dt_adv,dt_diff,g_dt_adv,g_dt_diff,min_dt_adv,min_dt_diff;
  PetscReal   el_coords[NSD*NODES_PER_EL_Q1_3D],el_V[NSD*NODES_PER_EL_Q1_3D],el_phi[NODES_PER_EL_Q1_3D];
  PetscInt    ge_eqnums[NODES_PER_EL_Q1_3D];
  Quadrature        volQ;
  QPntVolCoefEnergy *all_quadpoints,*cell_quadpoints;
  PetscInt          nqp;
  PetscScalar       *qp_xi,*qp_weight;
  PetscScalar       qp_kappa[27];//,qp_Q[27];
  PetscScalar       Ni_p[NODES_PER_EL_Q1_3D];
  PetscScalar       GNi_p[NSD][NODES_PER_EL_Q1_3D],GNx_p[NSD][NODES_PER_EL_Q1_3D];
  PetscScalar       J_p,fac;
  PetscInt          i,p,e,nel,nen;
  const PetscInt    *elnidx;
  PetscReal         phi_p,cg,c = 0.0;
  PetscErrorCode    ierr;

  PetscFunctionBegin;

  da   = energy->daT;
  volQ = energy->volQ;

  /* quadrature */
  nqp       = volQ->npoints;
  qp_xi     = volQ->q_xi_coor;
  qp_weight = volQ->q_weight;
  ierr = VolumeQuadratureGetAllCellData_Energy(volQ,&all_quadpoints);CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&philoc);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
  ierr = DMGetLocalVector(cda,&Vloc);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,philoc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (da,X,INSERT_VALUES,philoc);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(cda,energy->u_minus_V,INSERT_VALUES,Vloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(  cda,energy->u_minus_V,INSERT_VALUES,Vloc);CHKERRQ(ierr);

  ierr = DMGetCoordinatesLocal(da,&coordsloc);CHKERRQ(ierr);

  ierr = VecGetArray(philoc,   &LA_philoc);CHKERRQ(ierr);
  ierr = VecGetArray(Vloc,     &LA_Vloc);CHKERRQ(ierr);
  ierr = VecGetArray(coordsloc,&LA_coordsloc);CHKERRQ(ierr);

  ierr = DMDAGetElementsQ1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);

  c = 0.0;
  min_dt_adv  = 1.0e32;
  min_dt_diff = 1.0e32;
  for (e=0; e<nel; e++) {

    ierr = DMDAEQ1_GetElementLocalIndicesDOF(ge_eqnums,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);

    /* get coords for the element */
    ierr = DMDAEQ1_GetVectorElementField_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_coordsloc);CHKERRQ(ierr);
    /* get current temperature */
    ierr = DMDAEQ1_GetScalarElementField_3D(el_phi,(PetscInt*)&elnidx[nen*e],LA_philoc);CHKERRQ(ierr);
    /* get velocity for the element */
    ierr = DMDAEQ1_GetVectorElementField_3D(el_V,(PetscInt*)&elnidx[nen*e],LA_Vloc);CHKERRQ(ierr);
    //printf("el_v %1.4e %1.4e %1.4e \n", el_V[0],el_V[1],el_V[2]);

    ierr = VolumeQuadratureGetCellData_Energy(volQ,all_quadpoints,e,&cell_quadpoints);CHKERRQ(ierr);

    /* copy the diffusivity and force */
    for (p=0; p<nqp; p++) {
      qp_kappa[p] = cell_quadpoints[p].diffusivity;
      //qp_Q[p]     = cell_quadpoints[p].heat_source;
    }

    /*
    el_kappa_const = 0.0;
    for (n=0; n<nqp; n++) {
      el_kappa_const += qp_kappa[n];
    }
    el_kappa_const = el_kappa_const / ((PetscReal)nqp);
    */

    /* diagnostics */

    el_volume = 0.0;
    el_kappa_const = 0.0;
    for (p=0; p<nqp; p++) {
      P3D_ConstructNi_Q1_3D(&qp_xi[NSD*p],Ni_p);
      P3D_ConstructGNi_Q1_3D(&qp_xi[NSD*p],GNi_p);
      P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&J_p,&GNx_p[0],&GNx_p[1],&GNx_p[2]);

      fac = qp_weight[p]*J_p;

      phi_p = 0.0;
      for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
        phi_p = phi_p + Ni_p[i] * el_phi[i];
      }

      el_volume      = el_volume      + 1.0 * fac;
      el_kappa_const = el_kappa_const + qp_kappa[p] * fac;
      c              = c              + phi_p * fac;
    }
    el_kappa_const = el_kappa_const / el_volume;

    //ierr = DASUPG3dComputeElementTimestep_qp(el_coords,el_V,el_kappa_const,&dt_adv,&dt_diff);CHKERRQ(ierr);
    ierr = AdvDiff3dComputeElementTimestep_qp(el_coords,el_V,el_kappa_const,&dt_adv,&dt_diff);CHKERRQ(ierr);
    if (dt_adv < min_dt_adv) {
      min_dt_adv = dt_adv;
    }
    if (dt_diff < min_dt_diff) {
      min_dt_diff = dt_diff;
    }
  }
  ierr = MPI_Allreduce(&min_dt_adv, &g_dt_adv, 1,MPIU_REAL,MPIU_MIN,PetscObjectComm((PetscObject)da));CHKERRQ(ierr);
  ierr = MPI_Allreduce(&min_dt_diff,&g_dt_diff,1,MPIU_REAL,MPIU_MIN,PetscObjectComm((PetscObject)da));CHKERRQ(ierr);

  ierr = MPI_Allreduce(&c,&cg,1,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)da));CHKERRQ(ierr);
  /*
  PetscPrintf(PETSC_COMM_WORLD,"Diffusion dt = %1.12e \n", g_dt_diff );
  PetscPrintf(PETSC_COMM_WORLD,"Advective dt = %1.12e \n", g_dt_adv );
  PetscPrintf(PETSC_COMM_WORLD,"\\int \\phi dV = %1.12e \n", cg );
  */
  *timestep = g_dt_adv;
  if (g_dt_diff < g_dt_adv) {
    *timestep = g_dt_diff;
  }

  ierr = VecRestoreArray(Vloc,     &LA_coordsloc);CHKERRQ(ierr);
  ierr = VecRestoreArray(coordsloc,&LA_Vloc);CHKERRQ(ierr);
  ierr = VecRestoreArray(philoc,   &LA_philoc);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(cda,&Vloc);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da, &philoc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


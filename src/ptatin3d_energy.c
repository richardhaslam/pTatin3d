/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      ptatin3d_energy.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


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
#include "ptatin3d_energy.h"


#undef __FUNCT__
#define __FUNCT__ "pTatinGetContext_Energy"
PetscErrorCode pTatinGetContext_Energy(pTatinCtx ctx,PhysCompEnergy *e)
{
	PetscFunctionBegin;
	if (e) { *e = ctx->energy_ctx; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinContextValid_Energy"
PetscErrorCode pTatinContextValid_Energy(pTatinCtx ctx,PetscBool *exists)
{
	PetscFunctionBegin;
	*exists = PETSC_FALSE;
	if (ctx->energy_ctx) {
		*exists = PETSC_TRUE;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinPhysCompCreate_Energy"
PetscErrorCode pTatinPhysCompCreate_Energy(pTatinCtx user)
{
	PetscErrorCode ierr;
	PhysCompStokes stokes_ctx;
	
	PetscFunctionBegin;
	
	stokes_ctx = user->stokes_ctx;
	
	if (user->restart_from_file) {
		/* load from file */
		char t_name[PETSC_MAX_PATH_LEN];
		char told_name[PETSC_MAX_PATH_LEN];
		char xold_name[PETSC_MAX_PATH_LEN];
		char v_name[PETSC_MAX_PATH_LEN];
		
		/* dav,dap */
		if (!StringEmpty(user->restart_prefix)) {
			sprintf(t_name,   "%s/ptat3dcpf.dmda-t_%s",user->restart_dir,user->restart_prefix);
			sprintf(told_name,"%s/ptat3dcpf.dmda-t_old_%s",user->restart_dir,user->restart_prefix);
			sprintf(xold_name,"%s/ptat3dcpf.dmda-x_old_%s",user->restart_dir,user->restart_prefix);
			sprintf(v_name,   "%s/ptat3dcpf.dmda-u_minus_v_%s",user->restart_dir,user->restart_prefix);
		} else {
			sprintf(t_name,   "%s/ptat3dcpf.dmda-t",user->restart_dir);
			sprintf(told_name,"%s/ptat3dcpf.dmda-t_old",user->restart_dir);
			sprintf(xold_name,"%s/ptat3dcpf.dmda-x_old",user->restart_dir);
			sprintf(v_name,   "%s/ptat3dcpf.dmda-u_minus_v",user->restart_dir);
		}
		PetscPrintf(PETSC_COMM_WORLD,"  reading %s \n", t_name );
		PetscPrintf(PETSC_COMM_WORLD,"  reading %s \n", told_name );
		PetscPrintf(PETSC_COMM_WORLD,"  reading %s \n", xold_name );
		PetscPrintf(PETSC_COMM_WORLD,"  reading %s \n", v_name );
		
		ierr = PhysCompLoad_Energy();CHKERRQ(ierr);
	} else {
		/* create from data */
		PetscInt energy_mesh_type;
		
		energy_mesh_type = 1; /* default is Q1 overlapping Q2 */
		ierr = PetscOptionsGetInt(PETSC_NULL,"-energy_mesh_type",&energy_mesh_type,0);CHKERRQ(ierr);
		ierr = PhysCompNew_Energy(stokes_ctx->dav,-1,-1,-1,energy_mesh_type,&user->energy_ctx);CHKERRQ(ierr);
	}	
	
	
	if (user->restart_from_file) {
		
	} else {
		ierr = PhysCompAddMaterialPointCoefficients_Energy(user->materialpoint_db);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinPhysCompActivate_Energy"
PetscErrorCode pTatinPhysCompActivate_Energy(pTatinCtx user,PetscBool load)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	if (load && (user->energy_ctx == PETSC_NULL)) {
		ierr = pTatinPhysCompCreate_Energy(user);CHKERRQ(ierr);
	}	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinPhysCompAttachData_Energy"
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

#undef __FUNCT__  
#define __FUNCT__ "pTatinPhysCompGetData_Energy"
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

#undef __FUNCT__
#define __FUNCT__ "MaterialPointQuadraturePointProjectionC0_Q2Energy"
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
	DataBucketGetSizes(materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
	mp_std  = PField_std->data;
	
	
	ierr = MPntPEnergyComputeMemberOffsets(mp_property_offsets);CHKERRQ(ierr);
	ierr = QPntVolCoefEnergyComputeMemberOffsets(qp_property_offsets);CHKERRQ(ierr);
	
	/* setup */
	dof = 1;
	ierr = DMDADuplicateLayout(da,dof,1,DMDA_STENCIL_BOX,&clone);CHKERRQ(ierr);
	ierr = DMGetDMDAE(da,&dae);CHKERRQ(ierr);
	
	ierr = DMAttachDMDAE(clone);CHKERRQ(ierr);
	ierr = DMGetDMDAE(clone,&dae_clone);CHKERRQ(ierr);
	/*
	{
		PetscInt NP[3];

		ierr = DMDAGetInfo(da,0,0,0,0,&NP[0],&NP[1],&NP[2],0,0, 0,0,0, 0);CHKERRQ(ierr);		
		ierr = DMDAEDeepCopy(dae,NP,dae_clone);CHKERRQ(ierr);
	}*/
	ierr = DMDAECopy(dae,dae_clone);CHKERRQ(ierr);

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
	PetscOptionsGetBool(PETSC_NULL,"-view_projected_marker_fields",&view,PETSC_NULL);
	if (view) {
		char filename[256];
		PetscViewer viewer;
		
		sprintf(filename,"MaterialPointProjection_energy_member_%d.vtk",(int)member );
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = DMView(clone, viewer);CHKERRQ(ierr);
		ierr = VecView(properties_A, viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	}
	
	/* destroy */
	ierr = DMRestoreGlobalVector(clone,&properties_B);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(clone,&properties_A);CHKERRQ(ierr);
	
	ierr = DMDestroyDMDAE(clone);CHKERRQ(ierr);
	ierr = DMDestroy(&clone);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinPhysCompEnergy_MPProjectionQ1"
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
	ierr = pTatinGetMaterialPoints(ctx,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = MaterialPointQuadraturePointProjectionC0_Q2Energy(daT,materialpoint_db,MPField_Energy,MPPEgy_diffusivity,volQ);
	ierr = MaterialPointQuadraturePointProjectionC0_Q2Energy(daT,materialpoint_db,MPField_Energy,MPPEgy_heat_source,volQ);
	
	PetscFunctionReturn(0);
}


/* update V */
/*
 u - V = u - (X_current - X_old)/dt
			 = (dt.u - X_current + X_old)/dt
*/
#undef __FUNCT__  
#define __FUNCT__ "_pTatinPhysCompEnergy_UpdateALEVelocity"
PetscErrorCode _pTatinPhysCompEnergy_UpdateALEVelocity(PhysCompEnergy energy,PetscReal dt)
{
	Vec            coordinates;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = VecScale(energy->u_minus_V,dt);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(energy->daT,&coordinates);CHKERRQ(ierr);
	ierr = VecAXPY(energy->u_minus_V,-1.0,coordinates);CHKERRQ(ierr);
	ierr = VecAXPY(energy->u_minus_V, 1.0,energy->Xold);CHKERRQ(ierr);
	ierr = VecScale(energy->u_minus_V,1.0/dt);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinPhysCompEnergy_UpdateALEVelocity"
PetscErrorCode pTatinPhysCompEnergy_UpdateALEVelocity(PhysCompStokes s,Vec X,PhysCompEnergy energy,PetscReal dt)
{
	DM             cdaT;
	Vec            velocity,pressure;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetCoordinateDA(energy->daT,&cdaT);CHKERRQ(ierr);   
	
	/* Project fluid velocity from Q2 space into Q1 space */
	ierr = DMCompositeGetAccess(s->stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	ierr = DMDAProjectVectorQ2toQ1(s->dav,velocity,cdaT,energy->u_minus_V,energy->energy_mesh_type);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(s->stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	/* Compute ALE velocity in Q1 space */
	ierr = _pTatinPhysCompEnergy_UpdateALEVelocity(energy,dt);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "pTatinPhysCompEnergy_Update"
PetscErrorCode pTatinPhysCompEnergy_Update(PhysCompEnergy e,DM dav,Vec T)
{
	Vec            coords;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* update coords */
	ierr = DMDAProjectCoordinatesQ2toQ1(dav,e->daT,e->energy_mesh_type);CHKERRQ(ierr);
	
	/* update solution */
	ierr = VecCopy(T,e->Told);CHKERRQ(ierr);
	
	/* update coordinates */
	ierr = DMDAGetCoordinates(e->daT,&coords);CHKERRQ(ierr);
	ierr = VecCopy(coords,e->Xold);CHKERRQ(ierr);
	//ierr = DMDAUpdateGhostedCoordinates(daq1);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinPhysCompEnergy_Initialise"
PetscErrorCode pTatinPhysCompEnergy_Initialise(PhysCompEnergy e,Vec T)
{
	Vec            coords;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* update solution */
	ierr = VecCopy(T,e->Told);CHKERRQ(ierr);
	
	/* update coordinates */
	ierr = DMDAGetCoordinates(e->daT,&coords);CHKERRQ(ierr);
	ierr = VecCopy(coords,e->Xold);CHKERRQ(ierr);

	ierr = VecZeroEntries(T);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


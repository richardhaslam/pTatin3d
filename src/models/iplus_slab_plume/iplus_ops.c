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
 **    Filename:      iplus_ops.c
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
 **    $Id:$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


/*
 Lx,Ly,Lz
 
 eta {mantle,plume,slab}
   83, 5, ?

 rho {mantle,plume,slab}
   1413, 1373, ?
 
 typical A0,r0 we should compare with
 
*/

#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "dmda_element_q2p1.h"
#include "mesh_deformation.h"
#include "dmda_remesh.h"
#include "mesh_update.h"
#include "ptatin_models.h"
#include "model_utils.h"
#include "ptatin_utils.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "geometry_object.h"

#include "iplus_ctx.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_iPLUS"
PetscErrorCode ModelInitialize_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscBool        flg;
	char             logfile[PETSC_MAX_PATH_LEN];
	PetscInt         modeltype;
	PetscErrorCode   ierr;
	
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscOptionsGetInt(PETSC_NULL,"-iplus_modeltype",&modeltype,&flg);
	if (flg) {
		data->modeltype = (iPLUSModelType)modeltype;
	}
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: Using model type %D \n",(PetscInt)modeltype);
	
	data->mantle_eta = 85.0;    data->mantle_rho = 1413.0;
	data->plume_eta  = 5.0;     data->plume_rho  = 1373.0;
	data->slab_eta   = 64000.0; data->slab_rho   = 1495.0;
	
	PetscOptionsGetReal(PETSC_NULL,"-iplus_slab_eta",&data->slab_eta,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-iplus_slab_rho",&data->slab_rho,&flg);

	PetscOptionsGetReal(PETSC_NULL,"-iplus_mantle_eta",&data->mantle_eta,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-iplus_mantle_rho",&data->mantle_rho,&flg);

	PetscOptionsGetReal(PETSC_NULL,"-iplus_plume_eta",&data->plume_eta,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-iplus_plume_rho",&data->plume_rho,&flg);
	
	data->length_scale = 1.0;
	data->vel_scale    = 1.0;
	data->time_scale   = data->length_scale / data->vel_scale;
	data->eta_scale    = 1.0e3;
	if (data->modeltype == iPLUsModelPlume) {
		data->eta_scale = 1.0;
	}
	PetscOptionsGetReal(PETSC_NULL,"-iplus_eta_scale",&data->eta_scale,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-iplus_velocity_scale",&data->vel_scale,&flg);
	data->time_scale   = data->length_scale / data->vel_scale;
	
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: Using viscosity scale, eta* = %1.4e Pa s \n", data->eta_scale);
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: Using velocity scale,    v* = %1.4e m/s \n", data->vel_scale);
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: Using time scale,        t* = %1.4e sec \n", data->time_scale);
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: dimensional quantities \n");
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [mantle] eta = %1.4e Pa s ; rho = %1.4e kg/m^3 \n",data->mantle_eta,data->mantle_rho);
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [plume]  eta = %1.4e Pa s ; rho = %1.4e kg/m^3 \n",data->plume_eta,data->plume_rho);
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [slab]   eta = %1.4e Pa s ; rho = %1.4e kg/m^3 \n",data->slab_eta,data->slab_rho);
	
	data->mantle_eta = data->mantle_eta / data->eta_scale;
	data->plume_eta  = data->plume_eta / data->eta_scale;
	data->slab_eta   = data->slab_eta / data->eta_scale;

	data->mantle_rho = data->mantle_rho * data->length_scale * data->length_scale / (data->eta_scale * data->vel_scale);
	data->plume_rho  = data->plume_rho * data->length_scale * data->length_scale  / (data->eta_scale * data->vel_scale);
	data->slab_rho   = data->slab_rho * data->length_scale * data->length_scale   / (data->eta_scale * data->vel_scale);
	
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [mantle] eta* = %1.4e ; rho* = %1.4e \n",data->mantle_eta,data->mantle_rho);
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [plume]  eta* = %1.4e ; rho* = %1.4e \n",data->plume_eta,data->plume_rho);
	PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [slab]   eta* = %1.4e ; rho* = %1.4e \n",data->slab_eta,data->slab_rho);
	
	if (data->modeltype == iPLUsModelPlume) {
		data->plume_pos[0] = 0.15;
		data->plume_pos[2] = 0.15;
		data->plume_pos[1] = 0.0;
	} else {
		data->plume_pos[0] = 0.15;
		data->plume_pos[2] = 0.15;

		PetscOptionsGetReal(PETSC_NULL,"-iplus_plume_x",&data->plume_pos[0],&flg);
		PetscOptionsGetReal(PETSC_NULL,"-iplus_plume_z",&data->plume_pos[2],&flg);
		
		data->plume_pos[1] = 0.0;

	}
	
	data->plume_radius = 0.015625;
	data->plume_A0     = 1.0681;
	PetscOptionsGetReal(PETSC_NULL,"-iplus_plume_r",&data->plume_radius,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-iplus_plume_A0",&data->plume_A0,&flg);

	data->refinement_type = 0;
	PetscOptionsGetInt(PETSC_NULL,"-iplus_mesh_refinement_type",&data->refinement_type,&flg);
	
	data->np_plume_x = 10;
	data->np_plume_z = 10;
	PetscOptionsGetInt(PETSC_NULL,"-iplus_plume_npx",&data->np_plume_x,&flg);
	PetscOptionsGetInt(PETSC_NULL,"-iplus_plume_npz",&data->np_plume_z,&flg);

	if ( (data->modeltype == iPLUsModelPlume) || (data->modeltype == iPLUsModelSlabPlume) ) {
		PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [plume origin] %1.4e m ; %1.4e m ; %1.4e m\n",data->plume_pos[0],data->plume_pos[1],data->plume_pos[2]);
		PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [plume radius] %1.4e m \n",data->plume_radius);
		PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [plume A0]     %1.4e 1/(ms) \n",data->plume_A0);
		PetscPrintf(PETSC_COMM_WORLD, "  iPLUS: [plume points] %D x %D \n",data->np_plume_x,data->np_plume_z);
	}
	
	/* slab geometry */
	ierr = iPLUS_CreateSlabGeometry(data);CHKERRQ(ierr);
	
	sprintf(logfile,"%s/iplus.logfile",c->outputpath);
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,logfile,&data->logviewer);CHKERRQ(ierr);
	
	data->iplus_output_frequency = 1;
	PetscOptionsGetInt(PETSC_NULL,"-iplus_output_frequency",&data->iplus_output_frequency,&flg);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_iPLUS"
PetscErrorCode ModelApplyInitialMeshGeometry_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscBool        ribe_slab,wouter_slab,shallow_mantle;
	PetscReal        Gmin[3],Gmax[3],center_x,center_z,Lx,Ly,Lz;
	PetscErrorCode   ierr;
	
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

	if (data->modeltype == iPLUsModelPlume) {
		/* small tank */
		ierr = DMDASetUniformCoordinates(dav,0.0,0.3,0.0,0.3,0.0,0.3);CHKERRQ(ierr);
	} else {
		/* bigger tank */
		wouter_slab = PETSC_FALSE;
		PetscOptionsGetBool(NULL,"-iplus_slab_type_schellart_g3_2008",&wouter_slab,NULL);

		ribe_slab = PETSC_FALSE;
		PetscOptionsGetBool(NULL,"-iplus_slab_type_liribe_jgr_2012",&ribe_slab,NULL);
		
		shallow_mantle = PETSC_FALSE;
		PetscOptionsGetBool(NULL,"-iplus_slab_domain_shallow_mantle",&shallow_mantle,NULL);
		
		if(wouter_slab || ribe_slab) {

			if (shallow_mantle) {
				ierr = DMDASetUniformCoordinates(dav,0.0,1.0, 0.0,0.12, 0.0,0.6);CHKERRQ(ierr);
			} else {
				ierr = DMDASetUniformCoordinates(dav,0.0,1.0, 0.0,0.38, 0.0,0.6);CHKERRQ(ierr);
			}
			
		} else {
			ierr = DMDASetUniformCoordinates(dav,0.0,1.0,0.0,0.4,0.0,0.6);CHKERRQ(ierr);
		}

		
	}
	/* refine? */
	ierr = DMDAGetBoundingBox(dav,Gmin,Gmax);CHKERRQ(ierr);
	Lx = (Gmax[0] - Gmin[0]);
	Ly = (Gmax[1] - Gmin[1]);
	Lz = (Gmax[2] - Gmin[2]);
	center_x = (Gmax[0] + Gmin[0]) * 0.5;
	center_z = (Gmax[2] + Gmin[2]) * 0.5;
	switch (data->refinement_type) {
			
		case 1:
			PetscPrintf(PETSC_COMM_WORLD,"  iPLUS: [Mesh refinement] Type 1 (Subtle)\n");
			ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 0, 4.0, Gmin[0], center_x - 0.2*Lx*0.5, center_x + 0.2*Lx*0.5, Gmax[0]);CHKERRQ(ierr);
			ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 2, 4.0, Gmin[2], center_z - 0.2*Lz*0.5, center_z + 0.2*Lz*0.5, Gmax[2]);CHKERRQ(ierr);
			break;
		case 2:
			PetscPrintf(PETSC_COMM_WORLD,"  iPLUS: [Mesh refinement] Type 2 (Aggressive)\n");
			ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 0, 8.0, Gmin[0], center_x - 0.2*Lx*0.5, center_x + 0.2*Lx*0.5, Gmax[0]);CHKERRQ(ierr);
			ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 2, 8.0, Gmin[2], center_z - 0.2*Lz*0.5, center_z + 0.2*Lz*0.5, Gmax[2]);CHKERRQ(ierr);
			break;
		case 3:
		PetscPrintf(PETSC_COMM_WORLD,"  iPLUS: [Mesh refinement] Type 3 (x-z refinement)\n");
			ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 0, 4.0, Gmin[0], 0.25, 0.4, Gmax[0]);CHKERRQ(ierr);
			//ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 1, 4.0, Gmin[1], 0.75*Ly, Gmax[1], Gmax[1]);CHKERRQ(ierr);
			ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 2, 3.0, Gmin[2], center_z - 0.35*Lz*0.5, center_z + 0.35*Lz*0.5, Gmax[2]);CHKERRQ(ierr);

			break;
		case 4:
			PetscPrintf(PETSC_COMM_WORLD,"  iPLUS: [Mesh refinement] Type 4 (x-y-z refinement)\n");
			ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 0, 4.0, Gmin[0], 0.25, 0.4, Gmax[0]);CHKERRQ(ierr);
			//ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 1, 4.0, Gmin[1], 0.75*Ly, Gmax[1], Gmax[1]);CHKERRQ(ierr);
			ierr = DMDASetCoordinatesCentralSqueeze1D(dav, 2, 3.0, Gmin[2], center_z - 0.35*Lz*0.5, center_z + 0.35*Lz*0.5, Gmax[2]);CHKERRQ(ierr);
			
			ierr = DMDASetCoordinatesColumnRefinement(dav,1,4.0, 0.75,1.0);CHKERRQ(ierr);
			break;
	}
	
	/* determine elements located within the plume */
	ierr = iPLUS_DetermineElementsContainingPlumeInlet(dav,data);CHKERRQ(ierr);

	/* compute initial volume of the domain */
	ierr = DMDAComputeMeshVolume(dav,&data->intial_domain_volume);CHKERRQ(ierr);	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_iPLUS"
PetscErrorCode ModelApplyInitialMaterialGeometry_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	DataBucket       materialpoint_db;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;

	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointSet_phase_index(mpX,p,iPLUSMatMantle);CHKERRQ(ierr);
		ierr = MaterialPointSet_viscosity(mpX,p,data->mantle_eta);CHKERRQ(ierr);
		ierr = MaterialPointSet_density(mpX,p,-GRAVITY*data->mantle_rho);CHKERRQ(ierr);
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	/* insert slab? */
	if ( (data->modeltype == iPLUsModelSlab) || (data->modeltype == iPLUsModelSlabPlume) ) {
		ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
		stokes_pack = stokes->stokes_pack;
		ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
		ierr = iPLUS_DefineSlabMaterial(dav,materialpoint_db,data);CHKERRQ(ierr);
	}
	
	/* insert plume? */
	if ( (data->modeltype == iPLUsModelPlume) || (data->modeltype == iPLUsModelSlabPlume) ) {
		ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
		stokes_pack = stokes->stokes_pack;
		ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
		ierr = iPLUS_InsertPlumeMaterial(dav,materialpoint_db,data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_iPLUS"
PetscErrorCode ModelApplyInitialSolution_iPLUS(pTatinCtx c,Vec X,void *ctx)
{
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUS_VelocityBC"
PetscErrorCode iPLUS_VelocityBC(BCList bclist,DM dav,pTatinCtx c,iPLUSCtx *data)
{
	PetscScalar    zero;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	zero = 0.0;

	/* south face */
	if ( (data->modeltype == iPLUsModelPlume) || (data->modeltype == iPLUsModelSlabPlume) ) {
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,iPLUS_Inflow_BCListEvaluator,(void*)data);CHKERRQ(ierr);
	} else {
		/* no normal flow */
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}
	
	/* no slip bottom */
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

	/* no slip */
	/*
	for (d=0; d<3; d++) {
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,d,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,d,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,d,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,d,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}
	*/
	ierr = DirichletBC_FreeSlip(bclist,dav,FRONT_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,EAST_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,BACK_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,WEST_FACE);CHKERRQ(ierr);

	
	/* north face is free surface */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_iPLUS"
PetscErrorCode ModelApplyBoundaryCondition_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;

	
	PetscFunctionBegin;

	/* Define velocity boundary conditions */
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = iPLUS_VelocityBC(stokes->u_bclist,dav,c,data);CHKERRQ(ierr);

	/* Define boundary conditions for any other physics */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_iPLUS"
PetscErrorCode ModelApplyBoundaryConditionMG_iPLUS(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscInt         n;
	PetscErrorCode   ierr;
	
	
	PetscFunctionBegin;
	
	/* Define velocity boundary conditions on each level within the MG hierarchy */
	for (n=0; n<nl; n++) {
		ierr = iPLUS_VelocityBC(bclist[n],dav[n],c,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_iPLUS"
PetscErrorCode ModelApplyMaterialBoundaryCondition_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	if ( (data->modeltype == iPLUsModelPlume) || (data->modeltype == iPLUsModelSlabPlume) ) {
		ierr = iPLUS_ApplyMaterialBoundaryCondition_Plume(c,data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_iPLUS"
PetscErrorCode ModelApplyUpdateMeshGeometry_iPLUS(pTatinCtx c,Vec X,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscReal        step;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	PetscErrorCode   ierr;

	
	PetscFunctionBegin;
	
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	//ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);

	ierr = UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX(dav,velocity,PETSC_NULL,c->dt);CHKERRQ(ierr);
	/* perform vertical remeshing column wise */
	switch (data->refinement_type) {
		case 4:
			ierr = DMDASetCoordinatesColumnRefinement(dav,1,4.0, 0.75,1.0);CHKERRQ(ierr);
			break;
	}
	
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUSOutput_ComputeVerticalRangeOfRegion"
PetscErrorCode iPLUSOutput_ComputeVerticalRangeOfRegion(DataBucket materialpoint_db,PetscInt region_idx,PetscReal range_yp[])
{
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	PetscReal        _range_yp[2];
	double           *pos_p;
	int              region_p;
	PetscInt         found_region = 0,found_region_g;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;

	/* initialize */
	_range_yp[0] =  1.0e32;
	_range_yp[1] = -1.0e32;
	
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointGet_phase_index(mpX,p,&region_p);CHKERRQ(ierr);
		ierr = MaterialPointGet_global_coord(mpX,p,&pos_p);CHKERRQ(ierr);
		
		if (region_p == region_idx) {
			if (pos_p[1] < _range_yp[0]) { _range_yp[0] = pos_p[1]; }
			if (pos_p[1] > _range_yp[1]) { _range_yp[1] = pos_p[1]; }
			found_region = 1;
		}
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&_range_yp[0],&range_yp[0],1,MPIU_REAL,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&_range_yp[1],&range_yp[1],1,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);

	ierr = MPI_Allreduce(&found_region,&found_region_g,1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);

	/* no particles of the desired region were found on any processors, set min/max accordingly */
	if (found_region_g == 0) {
		range_yp[0] = -1.0e32;
		range_yp[1] =  1.0e32;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUSOutput_ComputeDomainVolume"
PetscErrorCode iPLUSOutput_ComputeDomainVolume(DM dav,PetscReal *volume)
{
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = DMDAComputeMeshVolume(dav,volume);CHKERRQ(ierr);	

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_iPLUS"
PetscErrorCode ModelOutput_iPLUS(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	static int       beenhere = 0;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	
	if (c->step%data->iplus_output_frequency == 0) {
		/* ---- Velocity-Pressure Mesh Output ---- */
		/* [1] Standard viewer: v,p written out as binary in double */
		//ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

		/* [2] Light weight viewer: Only v is written out. v and coords are expressed as floats */
		ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(c,X,prefix);CHKERRQ(ierr);
		
		/* [3] Write out v,p into PETSc Vec. These can be used to restart pTatin */
		/*
		ierr = pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
		*/

		
		/* ---- Material Point Output ---- */
		/* [1] Basic viewer: Only reports coords, regionid and other internal data */
		ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
		
		/* [2] Customized viewer: User defines specific fields they want to view - NOTE not .pvd file will be created */
		/*
		{
			DataBucket                materialpoint_db;
			const int                 nf = 4;
			const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
			char                      mp_file_prefix[256];
			
			ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
			sprintf(mp_file_prefix,"%s_mpoints",prefix);
			ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
		}
		*/
		/* [3] Customized marker->cell viewer: Marker data is projected onto the velocity mesh. User defines specific fields */
		/*
		{
			const int                    nf = 3;
			const MaterialPointVariable  mp_prop_list[] = { MPV_viscosity, MPV_density, MPV_plastic_strain }; 
			
			ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
		}	
		*/
	}
	
	/* iPlus specific output */
	{
		PhysCompStokes   stokes;
		DM               stokes_pack,dav,dap;
		DataBucket       materialpoint_db;
		PetscReal        volume,plume_range_yp[2],slab_range_yp[2];
		
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
		
		ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
		stokes_pack = stokes->stokes_pack;
		ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
		
		ierr = iPLUSOutput_ComputeDomainVolume(dav,&volume);CHKERRQ(ierr);
		ierr = iPLUSOutput_ComputeVerticalRangeOfRegion(materialpoint_db,iPLUSMatPlume,plume_range_yp);
		ierr = iPLUSOutput_ComputeVerticalRangeOfRegion(materialpoint_db,iPLUSMatSlab, slab_range_yp);
		
		if (beenhere == 0) {
			PetscViewerASCIIPrintf(data->logviewer,"# iPLUS logfile\n");
			PetscViewerASCIIPrintf(data->logviewer,"# Note: If the slab (or plume) is not present, the min/max y coordinatate reported will be -1.0e32/+1.0e32 \n");
			PetscViewerASCIIPrintf(data->logviewer,"# ----------------------------------------------------------------------------------------------------------------- \n");
			PetscViewerASCIIPrintf(data->logviewer,"# step  time (ND)       time (sec)      Omega(t=0)      Omega(t)        plume_{y_min,y_max}             slab_{y_min,y_max} \n");
			beenhere = 1;
		}
		PetscViewerASCIIPrintf(data->logviewer,"%D\t%1.4e\t%1.4e\t%1.6e\t%1.6e\t%+1.4e\t%+1.4e\t%+1.4e\t%+1.4e\n",
								c->step, c->time,c->time*data->time_scale, 
								data->intial_domain_volume, volume,
								plume_range_yp[0], plume_range_yp[1] ,slab_range_yp[0], slab_range_yp[1]);
		
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_iPLUS"
PetscErrorCode ModelDestroy_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	if (data->plume_element) {
		ierr = PetscFree(data->plume_element);CHKERRQ(ierr);
	}
	if (data->logviewer) {
		ierr = PetscViewerDestroy(&data->logviewer);CHKERRQ(ierr);
	}
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_iPLUS"
PetscErrorCode pTatinModelRegister_iPLUS(void)
{
	iPLUSCtx       *data;
	pTatinModel    m;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(iPLUSCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(iPLUSCtx));CHKERRQ(ierr);
	
	/* set initial values for model parameters */
	data->modeltype = iPLUsModelPlume;
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"iplus");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_iPLUS);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

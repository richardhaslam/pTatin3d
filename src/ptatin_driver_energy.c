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
 **    Filename:      ptatin_driver_energy.c
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

static const char help[] = "Stokes solver using Q2-Pm1 mixed finite elements.\n"
"3D prototype of the (p)ragmatic version of Tatin with energy solver. (pTatin3d_v0.0)\n\n";


#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "stokes_form_function.h"
#include "stokes_operators.h"
#include "dmda_redundant.h"
#include "dmda_update_coords.h"
#include "dmda_element_q1.h"
#include "dmda_view_petscvtk.h"
#include "ptatin3d_energy.h"

#include "MPntPEnergy_def.h"
#include "QPntVolCoefEnergy_def.h"


#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_energy_tester"
PetscErrorCode pTatin3d_energy_tester(int argc,char **argv)
{
	DM              multipys_pack,dav,dap;
	PetscErrorCode  ierr;
	pTatinCtx       user;
	Vec             X,F;
	PetscBool       active_energy;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dParseOptions(user);CHKERRQ(ierr);

	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);

	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	user->pack = user->stokes_ctx->stokes_pack;

	/* fetch some local variables */
	multipys_pack = user->pack;
	dav           = user->stokes_ctx->dav;
	dap           = user->stokes_ctx->dap;
	
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);	
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);

	/* generate energy solver */
	/* NOTE - Generating the thermal solver here will ensure that the initial geom on the mechanical model is copied */
	/* NOTE - Calling pTatinPhysCompActivate_Energy() after pTatin3dCreateMaterialPoints() is essential */
	{
		PetscBool activate_energy;
		
		PetscOptionsGetBool(PETSC_NULL,"-activate_energy",&activate_energy,0);
		ierr = pTatinPhysCompActivate_Energy(user,activate_energy);CHKERRQ(ierr);
	}
	ierr = pTatinContextValid_Energy(user,&active_energy);CHKERRQ(ierr);
	
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);

	{
		Vec Xu,Xp;
		ierr = DMCompositeGetAccess(multipys_pack,X,&Xu,&Xp);CHKERRQ(ierr);
		ierr = BCListInsert(user->stokes_ctx->u_bclist,Xu);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&Xu,&Xp);CHKERRQ(ierr);
	}
	
	/* test form function */
	{
		SNES snes;
		
		
		ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
		ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
		ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
		ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	}
	
	/* test data bucket viewer */
	DataBucketView(((PetscObject)multipys_pack)->comm, user->materialpoint_db,"materialpoint_stokes",DATABUCKET_VIEW_STDOUT);
	DataBucketView(PETSC_COMM_SELF, user->material_constants,"material_constants",DATABUCKET_VIEW_STDOUT);
	
	/* write out the initial condition */
	ierr = pTatinModel_Output(user->model,user,X,"icbc");CHKERRQ(ierr);
	
	/* test generic viewer */
	if (active_energy) {
		Vec T,FE,Told;
		Mat JE;
		MatStructure mstr;
		PhysCompEnergy energy;
		BCList bclist;
		DM     daT;
		PetscScalar zero;
		KSP kspT;
		
		
		/* mash in some diffusivity */
		{
			MPntPEnergy *material_point;
			DataField PField_energy;
			int p;
			int npoints;

			DataBucketGetSizes(user->materialpoint_db,&npoints,0,0);
			DataBucketGetDataFieldByName(user->materialpoint_db,MPntPEnergy_classname,&PField_energy);
			DataFieldGetAccess(PField_energy);
			
			for (p=0; p<npoints; p++) {
				DataFieldAccessPoint(PField_energy,p,   (void**)&material_point);
				material_point->diffusivity = 1.0;
				material_point->heat_source = 0.0;
			}
			DataFieldRestoreAccess(PField_energy);
		}
		{
			QPntVolCoefEnergy *quad_point;
			DataField PField_energy;
			int p;
			int npoints;

			DataBucketGetSizes(user->energy_ctx->volQ->properties_db,&npoints,0,0);
			DataBucketGetDataFieldByName(user->energy_ctx->volQ->properties_db,QPntVolCoefEnergy_classname,&PField_energy);
			DataFieldGetAccess(PField_energy);
			
			for (p=0; p<npoints; p++) {
				DataFieldAccessPoint(PField_energy,p,   (void**)&quad_point);
				quad_point->diffusivity = 1.0;
				quad_point->heat_source = 0.0;
			}
			DataFieldRestoreAccess(PField_energy);
		}
		
		
		/*  THERMAL ENERGY SOLVE  */
		energy = user->energy_ctx;
		Told   = energy->Told;
		bclist = energy->T_bclist;
		daT    = energy->daT;
		
		ierr = DMCreateGlobalVector(energy->daT,&T);CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(energy->daT,&FE);CHKERRQ(ierr);
		ierr = DMGetMatrix(energy->daT,MATAIJ,&JE);CHKERRQ(ierr);
		
		/* map velocity vector from Q2 space onto Q1 space */
		//ierr = DMCompositeGetAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
		//ierr = DMDAProjectVelocityQ2toQ1_2d(user->dav,velocity,energy->daT,energy->u_minus_V);CHKERRQ(ierr);
		//ierr = DMCompositeRestoreAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
		
		/* update V */
		/*
		 u - V = u - (X_current - X_old)/dt
		 = (dt.u - X_current + X_old)/dt
		 */
		
		user->dt = 1.0e-3;
		user->time = 1.0;

		zero = 0.0;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_KMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_KMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		
		ierr = FormJacobianEnergy(user->time,T,user->dt,&JE,&JE,&mstr,(void*)energy);CHKERRQ(ierr);
		ierr = FormFunctionEnergy(user->time,T,user->dt,FE,(void*)energy);CHKERRQ(ierr);

		ierr = VecSetRandom(FE,0);CHKERRQ(ierr);
		
		
		ierr = KSPCreate(PETSC_COMM_WORLD,&kspT);CHKERRQ(ierr);
		ierr = KSPSetOptionsPrefix(kspT,"T_");CHKERRQ(ierr);
		ierr = KSPSetOperators(kspT,JE,JE,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = KSPSetFromOptions(kspT);CHKERRQ(ierr);
		
		ierr = KSPSolve(kspT,FE,T);CHKERRQ(ierr);

		ierr = KSPDestroy(&kspT);CHKERRQ(ierr);
		
		{
			const int nf = 2;
			const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Energy }; 
			ierr = SwarmViewGeneric_ParaView(user->materialpoint_db,nf,mp_prop_list,user->outputpath,"test_MPStd_MPEnergy");CHKERRQ(ierr);
		}

		/* output energy mesh */
		ierr = DMDAViewPetscVTK(user->energy_ctx->daT,user->energy_ctx->Told,"phiOld_overlapping_q1.vtk");CHKERRQ(ierr);

	
		ierr = VecDestroy(&T);CHKERRQ(ierr);
		ierr = MatDestroy(&JE);CHKERRQ(ierr);
		ierr = VecDestroy(&FE);CHKERRQ(ierr);
	}	
	
	
	
	
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = pTatin3dDestroyContext(&user);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	
	ierr = PetscInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	ierr = pTatin3d_energy_tester(argc,argv);CHKERRQ(ierr);
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}

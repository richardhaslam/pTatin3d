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
 **    filename:   ptatin_driver_energy.c
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

static const char help[] = "Stokes solver using Q2-Pm1 mixed finite elements.\n"
"3D prototype of the (p)ragmatic version of Tatin with energy solver. (pTatin3d_v0.0)\n\n";


#include "ptatin3d.h"
#include "ptatin3d_defs.h"
#include "ptatin_init.h"
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
#include "dmda_iterator.h"
#include "dmda_view_petscvtk.h"
#include "ptatin3d_energy.h"
#include "energy_assembly.h"

#include "MPntPEnergy_def.h"
#include "QPntVolCoefEnergy_def.h"


#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_energy_tester"
PetscErrorCode pTatin3d_energy_tester(int argc,char **argv)
{
	pTatinModel     model;
	PhysCompStokes  stokes;
	PhysCompEnergy  energy;
	DataBucket      materialpoint_db,materialconstants_db;
	DM              multipys_pack,dav;
	PetscErrorCode  ierr;
	pTatinCtx       user;
	Vec             X,F;
	Mat             JE;
	Vec             T,f;
	PetscBool       active_energy;
	PetscBool       use_JFNK_T = PETSC_FALSE;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);

	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(model,user);CHKERRQ(ierr);
	ierr = pTatinGetMaterialConstants(user,&materialconstants_db);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	
	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)stokes->stokes_pack);CHKERRQ(ierr);
	user->pack = stokes->stokes_pack;

	/* fetch some local variables */
	multipys_pack = user->pack;
	dav           = stokes->dav;
	
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);	
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	ierr = pTatinGetMaterialPoints(user,&materialpoint_db,NULL);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(model,user);CHKERRQ(ierr);

	/* generate energy solver */
	/* NOTE - Generating the thermal solver here will ensure that the initial geometry on the mechanical model is copied */
	/* NOTE - Calling pTatinPhysCompActivate_Energy() after pTatin3dCreateMaterialPoints() is essential */
	{
		PetscBool load_energy = PETSC_FALSE;
		
		PetscOptionsGetBool(NULL,NULL,"-activate_energy",&load_energy,NULL);
		ierr = pTatinPhysCompActivate_Energy(user,load_energy);CHKERRQ(ierr);
		ierr = pTatinContextValid_Energy(user,&active_energy);CHKERRQ(ierr);
	}
	if (active_energy) {
		ierr = pTatinGetContext_Energy(user,&energy);CHKERRQ(ierr);

		ierr = PetscOptionsGetBool(NULL,NULL,"-use_jfnk_energy",&use_JFNK_T,NULL);CHKERRQ(ierr);
		
		ierr = DMCreateGlobalVector(energy->daT,&T);CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(energy->daT,&f);CHKERRQ(ierr);

		JE = NULL;
		if (!use_JFNK_T) {
			ierr = DMSetMatType(energy->daT,MATAIJ);CHKERRQ(ierr);
			ierr = DMCreateMatrix(energy->daT,&JE);CHKERRQ(ierr);
			ierr = MatSetFromOptions(JE);CHKERRQ(ierr);
		}

		ierr = pTatinPhysCompAttachData_Energy(user,T,JE);CHKERRQ(ierr);
	}
	
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(model,user);CHKERRQ(ierr);
	if (active_energy) {
		ierr = pTatinPhysCompEnergy_MPProjectionQ1(user);CHKERRQ(ierr);
	}

	/* test data bucket viewer */
	DataBucketView(PetscObjectComm((PetscObject)multipys_pack), materialpoint_db,"materialpoints",DATABUCKET_VIEW_STDOUT);
	DataBucketView(PetscObjectComm((PetscObject)multipys_pack), materialconstants_db,"materialconstants",DATABUCKET_VIEW_STDOUT);
	
	/* initial condition */
	ierr = pTatinModel_ApplyInitialSolution(model,user,X);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(model,user);CHKERRQ(ierr);

	/* insert boundary conditions */
	{
		Vec Xu,Xp;

		ierr = DMCompositeGetAccess(multipys_pack,X,&Xu,&Xp);CHKERRQ(ierr);
		ierr = BCListInsert(stokes->u_bclist,Xu);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&Xu,&Xp);CHKERRQ(ierr);

		if (active_energy) {
			ierr = BCListInsert(energy->T_bclist,T);CHKERRQ(ierr);
		}
	}

	/* write out the initial condition */
	ierr = pTatinModel_Output(model,user,X,"icbc");CHKERRQ(ierr);
	
	if (active_energy) {
		PetscReal    dx;
		SNES         snesT;
		PetscInt     tk;
		
		/*  THERMAL ENERGY SOLVE  */
		
		/* MAP V into adv_diff_v TODO */
		/* map velocity vector from Q2 space onto Q1 space */
		//ierr = DMGetCoordinateDM(daT,&cdaT);CHKERRQ(ierr);   
		//ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		//ierr = DMDAProjectVectorQ2toQ1(dav,velocity,cdaT,energy->u_minus_V,energy->energy_mesh_type);CHKERRQ(ierr);
		//ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		//ierr = DMDAProjectCoordinatesQ2toQ1(dav,daT,energy->energy_mesh_type);CHKERRQ(ierr);
		//ierr = DMDAProjectCoordinatesQ2toOverlappingQ1_3d(dav,daT);CHKERRQ(ierr);
		
		/* update V */
		/*
		 u - V = u - (X_current - X_old)/dt
		 = (dt.u - X_current + X_old)/dt
		 */
		


		ierr = pTatinPhysCompEnergy_Initialise(energy,T);CHKERRQ(ierr);
		energy->dt = 1.0;
		ierr = pTatinPhysCompEnergy_UpdateALEVelocity(stokes,X,energy,energy->dt);CHKERRQ(ierr);

		dx = 1.0/((PetscReal)(user->mx));
		user->dt   = 0.8 * (dx * dx) / 1.0;
		user->dt   = 0.1 * (dx) / 1.0;
		{ 
			PetscReal timestep;
			ierr = pTatinPhysCompEnergy_ComputeTimestep(energy,energy->Told,&timestep);CHKERRQ(ierr);
			ierr = pTatin_SetTimestep(user,"EnergyCFL",timestep);CHKERRQ(ierr);
		}
		energy->dt   = user->dt;
		
		user->time = 0.0;
		energy->time = user->time;
		
		for (tk=1; tk<=user->nsteps; tk++) {
			char stepname[PETSC_MAX_PATH_LEN];
			
			/* MAP Tin into Told */
			//ierr = VecCopy(T,Told);CHKERRQ(ierr);
			//ierr = VecZeroEntries(T);CHKERRQ(ierr);
			
			ierr = pTatinPhysCompEnergy_UpdateALEVelocity(stokes,X,energy,energy->dt);CHKERRQ(ierr);
			
			// crappy way - make it non-linear
	#if 0		
			ierr = TS_FormJacobianEnergy(user->time,T,user->dt,JE,JE,(void*)energy);CHKERRQ(ierr);
			ierr = TS_FormFunctionEnergy(user->time,T,user->dt,f,(void*)energy);CHKERRQ(ierr);

			//ierr = VecSetRandom(f,0);CHKERRQ(ierr);
			//ierr = VecSet(f,12.1);CHKERRQ(ierr);
					
			ierr = KSPCreate(PETSC_COMM_WORLD,&kspT);CHKERRQ(ierr);
			ierr = KSPSetOptionsPrefix(kspT,"T_");CHKERRQ(ierr);
			ierr = KSPSetOperators(kspT,JE,JE);CHKERRQ(ierr);
			ierr = KSPSetFromOptions(kspT);CHKERRQ(ierr);
			
			ierr = KSPSolve(kspT,f,T);CHKERRQ(ierr);

			ierr = KSPDestroy(&kspT);CHKERRQ(ierr);
	#endif
			
	#if 1
			ierr = SNESCreate(PETSC_COMM_WORLD,&snesT);CHKERRQ(ierr);
			ierr = SNESSetOptionsPrefix(snesT,"T_");CHKERRQ(ierr);

			ierr = SNESSetFunction(snesT,f,    SNES_FormFunctionEnergy,(void*)energy);CHKERRQ(ierr);
			if (use_JFNK_T) {
				ierr = SNESSetJacobian(snesT,NULL,NULL,SNES_FormJacobianEnergy,(void*)energy);CHKERRQ(ierr);
			} else {
				ierr = SNESSetJacobian(snesT,JE,JE,SNES_FormJacobianEnergy,(void*)energy);CHKERRQ(ierr);
			}
					
			ierr = SNESSetType(snesT,SNESKSPONLY);
			ierr = SNESSetFromOptions(snesT);CHKERRQ(ierr);

			ierr = SNESSolve(snesT,NULL,T);CHKERRQ(ierr);
			
			ierr = SNESDestroy(&snesT);CHKERRQ(ierr);
	#endif
			ierr = pTatinPhysCompEnergy_Update(energy,dav,T);CHKERRQ(ierr);
			
			user->time = user->time + user->dt;
			
			PetscSNPrintf(stepname,PETSC_MAX_PATH_LEN-1,"step%.4D",tk);
			ierr = pTatinModel_Output(model,user,X,stepname);CHKERRQ(ierr);
		
		}
	}	
	
	
	ierr = VecDestroy(&T);CHKERRQ(ierr);
	if (JE) { ierr = MatDestroy(&JE);CHKERRQ(ierr); }
	ierr = VecDestroy(&f);CHKERRQ(ierr);
	
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
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	ierr = pTatin3d_energy_tester(argc,argv);CHKERRQ(ierr);
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

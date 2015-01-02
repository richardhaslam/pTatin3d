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
 **    filename:   phys_comp_energy.c
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

/*
 
 Physics componenet: Energy (advection-diffusion equation)
 
 Require functions:
 
 PhysCompCreateMesh_Energy()
 PhysCompCreateBoundaryList_Energy()
 PhysCompCreateVolumeQuadrature_Energy()
 
 PhysCompCreate_Energy() : 
 PhysCompDestroy_Energy() : 
 
 PhysCompsNew_Energy() : Allocates new structure
 PhysCompLoad_Energy() : Loads from checkpoint file 
 PhysCompSave_Energy() : Saves to checkpoint file

 pTatin3d_PhysCompCreate_Energy() : Reads checkpoint data or calls PhysCompsNew_Energy()
 
*/


#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "data_bucket.h"
#include "dmda_bcs.h"
#include "element_utils_q1.h"
#include "dmdae.h"
#include "dmda_element_q1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"
#include "QPntVolCoefEnergy_def.h"
#include "phys_comp_energy.h"


#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreate_Energy"
PetscErrorCode PhysCompCreate_Energy(PhysCompEnergy *E)
{
	PetscErrorCode ierr;
	PhysCompEnergy energy;
	
	PetscFunctionBegin;
	ierr = PetscMalloc(sizeof(struct _p_PhysCompEnergy),&energy);CHKERRQ(ierr);
	ierr = PetscMemzero(energy,sizeof(struct _p_PhysCompEnergy));CHKERRQ(ierr);
	*E = energy;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompDestroy_Energy"
PetscErrorCode PhysCompDestroy_Energy(PhysCompEnergy *E)
{
	PetscErrorCode ierr;
	PhysCompEnergy ctx;
	
	PetscFunctionBegin;
	
	if (!E) {PetscFunctionReturn(0);}
	ctx = *E;
	
	//	for (e=0; e<HEX_FACES; e++) {
	//		if (ctx->surfQ[e]) { ierr = SurfaceQuadratureDestroy(&ctx->surfQ[e]);CHKERRQ(ierr); }
	//	}
	if (ctx->volQ) { ierr = QuadratureDestroy(&ctx->volQ);CHKERRQ(ierr); }
	if (ctx->T_bclist) { ierr = BCListDestroy(&ctx->T_bclist);CHKERRQ(ierr); }
  if (ctx->daT) { ierr = DMDestroy(&ctx->daT);CHKERRQ(ierr); }
  if (ctx->u_minus_V) { ierr = VecDestroy(&ctx->u_minus_V);CHKERRQ(ierr); }
  if (ctx->Told) {      ierr = VecDestroy(&ctx->Told);CHKERRQ(ierr); }
  if (ctx->Xold) {      ierr = VecDestroy(&ctx->Xold);CHKERRQ(ierr); }
	if (ctx) { ierr = PetscFree(ctx);CHKERRQ(ierr); }
	
	*E = NULL;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateMesh_Energy"
PetscErrorCode PhysCompCreateMesh_Energy(PhysCompEnergy E,DM dav,PetscInt mx,PetscInt my, PetscInt mz,PetscInt mesh_generator_type)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	E->energy_mesh_type = mesh_generator_type;

	switch (mesh_generator_type) {
		DMDAE dae;
			
		case 0:
			PetscPrintf(PETSC_COMM_WORLD,"PhysCompCreateMesh_Energy: Generating standard Q1 DMDA\n");
			E->mx = mx;
			E->my = my;
			E->mz = mz;
			
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only overlapping and nested supported {1,2} ");
			break;
		
		case 1:
			PetscPrintf(PETSC_COMM_WORLD,"PhysCompCreateMesh_Energy: Generating overlapping Q1 DMDA\n");
			if (!dav) {
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Require valid DM dav");
			}
			ierr = DMDACreateOverlappingQ1FromQ2(dav,1,&E->daT);CHKERRQ(ierr);
			ierr = DMGetDMDAE(E->daT,&dae);CHKERRQ(ierr);
			E->mx = dae->mx;
			E->my = dae->my;
			E->mz = dae->mz;
			break;
		
		case 2:
			PetscPrintf(PETSC_COMM_WORLD,"PhysCompCreateMesh_Energy: Generating nested Q1 DMDA\n");
			if (!dav) {
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Require valid DM dav");
			}
			ierr = DMDACreateNestedQ1FromQ2(dav,1,&E->daT);CHKERRQ(ierr);
			ierr = DMGetDMDAE(E->daT,&dae);CHKERRQ(ierr);
			E->mx = dae->mx;
			E->my = dae->my;
			E->mz = dae->mz;
			break;

		default:
			E->energy_mesh_type = 1;

			PetscPrintf(PETSC_COMM_WORLD,"PhysCompCreateMesh_Energy: Generating overlapping Q1 DMDA\n");
			if (!dav) {
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Require valid DM dav");
			}
			ierr = DMDACreateOverlappingQ1FromQ2(dav,1,&E->daT);CHKERRQ(ierr);
			ierr = DMGetDMDAE(E->daT,&dae);CHKERRQ(ierr);
			E->mx = dae->mx;
			E->my = dae->my;
			E->mz = dae->mz;
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateBoundaryList_Energy"
PetscErrorCode PhysCompCreateBoundaryList_Energy(PhysCompEnergy E)
{
	DM daT;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	daT = E->daT;
	if (!daT) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"daT must be set"); }
	ierr = DMDABCListCreate(daT,&E->T_bclist);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateVolumeQuadrature_Energy"
PetscErrorCode PhysCompCreateVolumeQuadrature_Energy(PhysCompEnergy E)
{
	PetscInt dim, np_per_dim, ncells;
	DMDAE dae;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	np_per_dim = 2;
	dim = 3;

	ierr = DMGetDMDAE(E->daT,&dae);CHKERRQ(ierr);
	ncells = dae->lmx * dae->lmy * dae->lmz;
	
	ierr = VolumeQuadratureCreate_GaussLegendreEnergy(dim,np_per_dim,ncells,&E->volQ);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompNew_Energy"
PetscErrorCode PhysCompNew_Energy(DM dav,PetscInt mx,PetscInt my, PetscInt mz,PetscInt mesh_generator_type,PhysCompEnergy *E)
{
	PetscErrorCode  ierr;
	PhysCompEnergy  energy;
	DM              cda;
	
	PetscFunctionBegin;
	
	ierr = PhysCompCreate_Energy(&energy);CHKERRQ(ierr);

	ierr = PhysCompCreateMesh_Energy(energy,dav,mx,my,mz,mesh_generator_type);CHKERRQ(ierr);
	ierr = PhysCompCreateBoundaryList_Energy(energy);CHKERRQ(ierr);
	ierr = PhysCompCreateVolumeQuadrature_Energy(energy);CHKERRQ(ierr);
	
	/* create aux vectors */
	ierr = DMCreateGlobalVector(energy->daT,&energy->Told);CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(energy->daT,&cda);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(cda,&energy->u_minus_V);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(cda,&energy->Xold);CHKERRQ(ierr);
	
	*E = energy;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompLoad_Energy"
PetscErrorCode PhysCompLoad_Energy(void)
{
	PetscFunctionBegin;
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Currently unavailable");
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompSave_Energy"
PetscErrorCode PhysCompSave_Energy(void)
{
	PetscFunctionBegin;
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Currently unavailable");
	PetscFunctionReturn(0);
}

/* quadrature */
#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureCreate_GaussLegendreEnergy"
PetscErrorCode VolumeQuadratureCreate_GaussLegendreEnergy(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature)
{
	Quadrature Q;
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
	
	ierr = QuadratureCreate(&Q);CHKERRQ(ierr);
	Q->dim  = nsd;
	Q->type = VOLUME_QUAD;
	
	PetscPrintf(PETSC_COMM_WORLD,"VolumeQuadratureCreate_GaussLegendreEnergy:\n");
	switch (np_per_dim) {
		case 1:
			PetscPrintf(PETSC_COMM_WORLD,"\tUsing 1 pnt Gauss Legendre quadrature\n");
			//QuadratureCreateGauss_1pnt_3D(&ngp,gp_xi,gp_weight);
			break;
			
		case 2:
			PetscPrintf(PETSC_COMM_WORLD,"\tUsing 2x2 pnt Gauss Legendre quadrature\n");
			QuadratureCreateGauss_2pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
			
		case 3:
			PetscPrintf(PETSC_COMM_WORLD,"\tUsing 3x3 pnt Gauss Legendre quadrature\n");
			QuadratureCreateGauss_3pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
			
		default:
			PetscPrintf(PETSC_COMM_WORLD,"\tUsing 3x3 pnt Gauss Legendre quadrature\n");
			QuadratureCreateGauss_3pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
	}
	
	Q->n_elements = ncells;
	if (ncells!=0) {
		
		DataBucketCreate(&Q->properties_db);
		DataBucketRegisterField(Q->properties_db,QPntVolCoefEnergy_classname, sizeof(QPntVolCoefEnergy),NULL);
		DataBucketFinalize(Q->properties_db);
		
		DataBucketSetInitialSizes(Q->properties_db,Q->npoints*ncells,1);
		
		DataBucketView(PETSC_COMM_WORLD, Q->properties_db,"GaussLegendre EnergyCoefficients",DATABUCKET_VIEW_STDOUT);
	}
	
	*quadrature = Q;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureGetAllCellData_Energy"
PetscErrorCode VolumeQuadratureGetAllCellData_Energy(Quadrature Q,QPntVolCoefEnergy *coeffs[])
{
	QPntVolCoefEnergy *quadraturepoint_data;
  DataField          PField;
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(Q->properties_db, QPntVolCoefEnergy_classname ,&PField);
	quadraturepoint_data = PField->data;
	*coeffs = quadraturepoint_data;
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureGetCellData_Energy"
PetscErrorCode VolumeQuadratureGetCellData_Energy(Quadrature Q,QPntVolCoefEnergy coeffs[],PetscInt cidx,QPntVolCoefEnergy *cell[])
{
  PetscFunctionBegin;
	if (cidx>=Q->n_elements) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_SIZ,"cidx > max cells");
	}
	
	*cell = &coeffs[cidx*Q->npoints];
  PetscFunctionReturn(0);
}

/* material points */
#undef __FUNCT__  
#define __FUNCT__ "PhysCompAddMaterialPointCoefficients_Energy"
PetscErrorCode PhysCompAddMaterialPointCoefficients_Energy(DataBucket db)
{
	PetscFunctionBegin;
	/* register marker structures here */
	DataBucketRegisterField(db,MPntPEnergy_classname,sizeof(MPntPEnergy),NULL);

	PetscFunctionReturn(0);
}

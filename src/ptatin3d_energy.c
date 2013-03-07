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
#include "element_utils_q1.h"
#include "dmda_element_q1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"

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
	if (load) {
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


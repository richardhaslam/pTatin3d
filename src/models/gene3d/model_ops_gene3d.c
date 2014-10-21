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
 **    Filename:      model_ops_gene3d.c
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

/*
 Developed by Laetitia Le Pourhiet [laetitia.le_pourhiet@upmc.fr]
 */

#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "data_bucket.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "phase_map.h"

#include "model_gene3d_ctx.h"

const char MODEL_NAME[] = "model_GENE3D_";

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Gene3D"
PetscErrorCode ModelInitialize_Gene3D(pTatinCtx c,void *ctx)
{
    ModelGENE3DCtx    *data = (ModelGENE3DCtx*)ctx;
    RheologyConstants *rheology;
    PetscBool         flg, found;
    char              *option_name;
    PetscInt          i,nphase;
    PetscErrorCode    ierr;
	
    PetscFunctionBegin;
	
	
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    rheology = &c->rheology_constants;
	
	
    /* model geometry */
    PetscPrintf(PETSC_COMM_WORLD,"reading model initial geometry from options\n");
    ierr = PetscOptionsGetInt(NULL, "-initial_geom",(PetscInt *) & data->initial_geom, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide a type of material index initialisation \n");
	}
	
    /* box geometry */
    PetscPrintf(PETSC_COMM_WORLD, "reading box geometry from options\n");
	
    ierr = PetscOptionsGetReal(NULL, "-Lx", &data->Lx, &flg);CHKERRQ(ierr);
    if (found == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide model length Lx \n");
	}
	
    ierr = PetscOptionsGetReal(NULL, "-Ly", &data->Ly, &flg);CHKERRQ(ierr);
    if (found == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Expected user to provide model length Ly \n");
	}
	
    ierr = PetscOptionsGetReal(NULL, "-Lz", &data->Lz, &flg);CHKERRQ(ierr);
    if (found == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide model length Lz \n");
	}
	
	
    ierr = PetscOptionsGetReal(NULL, "-Ox", &data->Ox, &flg);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide model Origin Ox \n");
	}
	
    ierr = PetscOptionsGetReal(NULL, "-Oy", &data->Oy, &flg);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide model Origin Oy \n");
	}
	
    ierr = PetscOptionsGetReal(NULL, "-Oz", &data->Oz, &flg);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide model Origin Oz \n");
	}
	
	
	
    /* rheology type */
    PetscPrintf(PETSC_COMM_WORLD, "reading rheology type from options\n");
    ierr = PetscOptionsGetInt(MODEL_NAME, "-rheol",(PetscInt*) & rheology->rheology_type, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Expected user to provide value for rheology type \n");
	}
	
    /* material properties */
    PetscPrintf(PETSC_COMM_WORLD,"reading material properties from options\n");
    ierr = PetscOptionsGetInt(MODEL_NAME, "-nphase", &nphase, &found);CHKERRQ(ierr);
	
	/* set the active phases on the rheo structure */
	rheology->nphases_active = nphase;
	
    if (found == PETSC_FALSE)	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Expected user to provide value for number of materials \n");
	}
    for (i = 0; i < nphase; i++){
		
		switch (rheology->rheology_type){
			case 0:
			{
				asprintf(&option_name, "-eta0_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->const_eta0[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
					free (option_name);
				}
				asprintf (&option_name, "-rho_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->const_rho0[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
			}
				break;
				
			case 1:
			{
				asprintf(&option_name, "-eta0_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->const_eta0[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1 (PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",	option_name);
				}
				free(option_name);
				
				asprintf(&option_name, "-rho_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->const_rho0[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf(&option_name, "-Co_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->mises_tau_yield[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-Phi_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->dp_pressure_dependance[i]),&found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);CHKERRQ(ierr);
				}
				free (option_name);
				
				asprintf (&option_name, "-Tens_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->tens_cutoff[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-Hs_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->Hst_cutoff[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
			}
				break;
			case 2:
			{
				asprintf (&option_name, "-eta0_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->const_eta0[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-rho_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->const_rho0[i]), &found);CHKERRQ(ierr);
				asprintf (&option_name, "-Co_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->mises_tau_yield[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-Phi_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->dp_pressure_dependance[i]),&found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-Tens_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->tens_cutoff[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-Hs_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->Hst_cutoff[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-Co_inf_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->soft_Co_inf[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-Phi_inf_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->soft_phi_inf[i]), &found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-eps_min_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->soft_min_strain_cutoff[i]),&found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
				
				asprintf (&option_name, "-eps_max_%d", i);
				ierr = PetscOptionsGetReal(MODEL_NAME, option_name,&(rheology->soft_max_strain_cutoff[i]),&found);CHKERRQ(ierr);
				if (found == PETSC_FALSE) {
					SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,"Expected user to provide value to option %s \n",option_name);
				}
				free (option_name);
			}
				break;
			case 3:
			{
				SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "visco-elasto-plastic rheology not handled in 3D tatin \n");
			}
		}
	}
	
    /* bc type */
    data->boundary_conditon_type = GENEBC_FreeSlip;
	
	
    /* set initial values for model parameters */
	
	
    PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Gene3D"
PetscErrorCode ModelApplyBoundaryCondition_Gene3D(pTatinCtx user,void *ctx)
{
    ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
    PetscScalar zero = 0.0;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    switch (data->boundary_conditon_type)
	{
        case GENEBC_FreeSlip:
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            break;
			
        case GENEBC_NoSlip:
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            break;
			
        case GENEBC_FreeSlipFreeSurface:
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            break;
			
        case GENEBC_NoSlipFreeSurface:
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            break;
        default:
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
    PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Gene3D"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Gene3D(pTatinCtx c,void *ctx)
{
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
    PetscPrintf(PETSC_COMM_WORLD, "  NOT IMPLEMENTED \n", __FUNCT__);
	
    PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Gene3D"
PetscErrorCode ModelApplyInitialMeshGeometry_Gene3D(pTatinCtx c,void *ctx)
{
    ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    ierr =	DMDASetUniformCoordinates(c->stokes_ctx->dav, data->Ox, data->Lx, data->Oy, data->Ly, data->Oz, data->Lz); CHKERRQ(ierr);
	
    PetscFunctionReturn (0);
}

//=====================================================================================================================================

#undef __FUNCT__
#define __FUNCT__ "ModelSetMarkerIndexLayeredCake_Gene3D"
PetscErrorCode ModelSetMarkerIndexLayeredCake_Gene3D (pTatinCtx c,void *ctx)
/* define phase index on material points from a map file extruded in z direction */
{
    ModelGENE3DCtx *data = (ModelGENE3DCtx*)ctx;
    PetscInt i, nLayer;
    int p,n_mp_points;
    DataBucket db;
    DataField PField_std;
    int phase;
    PetscInt phaseLayer[LAYER_MAX];
    PetscScalar YLayer[LAYER_MAX + 1];
    char *option_name;
    PetscErrorCode ierr;
    PetscBool flg;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    /* define properties on material points */
    db = c->materialpoint_db;
    DataBucketGetDataFieldByName(db, MPntStd_classname, &PField_std);
    DataFieldGetAccess(PField_std);
    DataFieldVerifyAccess(PField_std, sizeof (MPntStd));
    DataBucketGetSizes(db, &n_mp_points, 0, 0);
	
    /* read layers from options */
    nLayer = 1;
    ierr = PetscOptionsGetInt(MODEL_NAME, "-nlayer", &nLayer, &flg);CHKERRQ(ierr);
    YLayer[0] = data->Oy;
    for (i=1; i<=nLayer; i++) {
        
		asprintf (&option_name, "-layer_y_%d", i);
		ierr = PetscOptionsGetReal(MODEL_NAME, option_name, &YLayer[i], &flg);CHKERRQ(ierr);
		if (flg == PETSC_FALSE) {
			/* NOTE - these error messages are useless if you don't include "model_Gene3D_" in the statement. I added &name[1] so that the "-" is skipped */
			SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"Expected user to provide value to option -model_Gene3D_%s \n",&option_name[1]);
		}
		free (option_name);
        
		if (YLayer[i] > YLayer[i-1]) {
			asprintf (&option_name, "-layer_phase_%d", i-1);
			ierr = PetscOptionsGetInt(MODEL_NAME, option_name, &phaseLayer[i-1],&flg);CHKERRQ(ierr);
			
			if (flg == PETSC_FALSE) {
				SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,"Expected user to provide value to option -model_Gene3D_%s \n",&option_name[1]);
			}
			free (option_name);
		} else {
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,"Layers must be entered so that layer_y[i] is larger than layer_y[i-1]\n");
		}
	}
	
	for (i=0; i<nLayer; i++) {
		PetscPrintf(PETSC_COMM_WORLD,"layer [%D] :  y coord range [%1.2e -- %1.2e] : phase index [%D] \n", i,YLayer[i],YLayer[i+1],phaseLayer[i]);
	}
	
    for (p = 0; p < n_mp_points; p++) {
		MPntStd *material_point;
		double *pos;
		
		DataFieldAccessPoint(PField_std, p, (void **) &material_point);
		MPntStdGetField_global_coord(material_point,&pos);
        
		for (i=0; i<nLayer; i++) {
			if ( (pos[1] >= YLayer[i]) && (pos[1] <= YLayer[i+1]) ) {
				phase = phaseLayer[i];
				break;
			}
		}
		/* check if the break was never executed */
		if (i==nLayer) {
			SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Unable to detect layer containing the marker with coordinates (%1.4e,%1.4e,%1.4e)",pos[0],pos[1],pos[2]);
		}
		
		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point, phase);
	}
	
    DataFieldRestoreAccess(PField_std);
    PetscFunctionReturn (0);
}

//===============================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "ModelSetMarkerIndexFromMap_Gene3D"
PetscErrorCode ModelSetMarkerIndexFromMap_Gene3D(pTatinCtx c,void *ctx)
/* define phase index on material points from a map file extruded in z direction */
{
    PetscErrorCode ierr;
    PhaseMap phasemap;
    PetscInt dir_0,dir_1,direction;
    DataBucket db;
    int p,n_mp_points;
    DataField PField_std;
    int phase_init, phase, phase_index, is_valid;
    char map_file[PETSC_MAX_PATH_LEN], *name;
    PetscBool flg;
    
	PetscFunctionBegin;
    PetscPrintf (PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    ierr = PetscOptionsGetString(MODEL_NAME,"-map_file",map_file,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (flg == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide a map file \n");
	}
    ierr = PetscOptionsGetInt(MODEL_NAME,"-extrude_dir",&direction,&flg);CHKERRQ(ierr);
    if (flg == PETSC_FALSE) {
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Expected user to provide an extrusion direction \n");
	}
    
    switch (direction){
        case 0:{
            dir_0 = 2;
            dir_1 = 1;
        }
            break;
            
        case 1:{
            dir_0 = 0;
            dir_1 = 2;
        }
            break;
            
        case 2:{
            dir_0 = 0;
            dir_1 = 1;
        }
            break;
    }
    
    asprintf(&name,"./inputdata/%s.pmap",map_file);
    PhaseMapLoadFromFile(name,&phasemap);
    free(name);
	
    asprintf(&name,"./inputdata/%s_phase_map.gp",map_file);
    PhaseMapViewGnuplot(name,phasemap);
    free(name);
	
	
    /* define properties on material points */
    db = c->materialpoint_db;
    DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
    DataFieldGetAccess(PField_std);
    DataFieldVerifyAccess(PField_std,sizeof (MPntStd));
	
    DataBucketGetSizes(db,&n_mp_points,0,0);
	
    for (p=0; p<n_mp_points; p++)
	{
		MPntStd *material_point;
		double position2D[2],*pos;
		
        
		DataFieldAccessPoint(PField_std, p, (void **) &material_point);
		MPntStdGetField_global_coord(material_point,&pos);
        
		position2D[0] = pos[dir_0];
		position2D[1] = pos[dir_1];
		
		MPntStdGetField_phase_index(material_point, &phase_init);
		
		PhaseMapGetPhaseIndex(phasemap, position2D, &phase_index);
		
		PhaseMapCheckValidity(phasemap, phase_index, &is_valid);
		//PetscPrintf(PETSC_COMM_WORLD,"Phase index : %d  is_valid %d \n", phase_index,is_valid);
		
		if (is_valid == 1) {			/* point located in the phase map */
			phase = phase_index;
		} else {
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,"marker outside the domain\n your phasemap is smaller than the domain \n please check your parameters and retry");
		}
		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point, phase);
		
	}
    PhaseMapDestroy(&phasemap);
    DataFieldRestoreAccess(PField_std);
	
    PetscFunctionReturn (0);
}


//======================================================================================================================================

#undef __FUNCT__
#define __FUNCT__ "ModelSetInitialStokesVariableOnMarker_Gene3D"
PetscErrorCode ModelSetInitialStokesVariableOnMarker_Gene3D(pTatinCtx c,void *ctx)
/* define properties on material points */
{
    int p, n_mp_points;
    DataBucket db;
    DataField PField_std, PField_stokes;
    int phase_index;
    RheologyConstants *rheology;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    rheology = &c->rheology_constants;
    db = c->materialpoint_db;
	
    DataBucketGetDataFieldByName(db, MPntStd_classname, &PField_std);
    DataFieldGetAccess(PField_std);
    DataFieldVerifyAccess(PField_std, sizeof (MPntStd));
    DataBucketGetDataFieldByName (db, MPntPStokes_classname, &PField_stokes);
    DataFieldGetAccess(PField_stokes);
    DataFieldVerifyAccess(PField_stokes, sizeof (MPntPStokes));
	
    DataBucketGetSizes(db, &n_mp_points, 0, 0);
	
    for (p = 0; p < n_mp_points; p++)
	{
		MPntStd *material_point;
		MPntPStokes *mpprop_stokes;
		
		DataFieldAccessPoint(PField_std, p, (void **) &material_point);
		DataFieldAccessPoint(PField_stokes, p, (void **) &mpprop_stokes);
		MPntStdGetField_phase_index (material_point, &phase_index);
		MPntPStokesSetField_eta_effective(mpprop_stokes,rheology->const_eta0[phase_index]);
		MPntPStokesSetField_density(mpprop_stokes,rheology->const_rho0[phase_index]);
	}
	
    DataFieldRestoreAccess(PField_std);
    DataFieldRestoreAccess(PField_stokes);
	
    PetscFunctionReturn (0);
}

//======================================================================================================================================

#undef __FUNCT__
#define __FUNCT__ "ModelGene3DInit"
PetscErrorCode ModelGene3DInit(DataBucket db)
{
    int                p,n_mp_points;
    DataField          PField_std;

    PetscFunctionBegin;
	
    DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
    DataFieldGetAccess(PField_std);
    DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
    DataBucketGetSizes(db,&n_mp_points,0,0);
	
    for (p=0; p<n_mp_points; p++) {
		int phase_index;
		MPntStd *material_point;
        
		DataFieldAccessPoint(PField_std,p,(void**)&material_point);
        
		MPntStdGetField_phase_index(material_point,&phase_index);
		phase_index = -1;
		MPntStdSetField_phase_index(material_point,phase_index);
	}
    DataFieldRestoreAccess(PField_std);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelGene3DCheckPhase"
PetscErrorCode ModelGene3DCheckPhase(DataBucket db,RheologyConstants *rheology)
{
    int                p,n_mp_points;
    DataField          PField_std;

    PetscFunctionBegin;
    DataBucketGetDataFieldByName (db,MPntStd_classname,&PField_std);
    DataFieldGetAccess(PField_std);
    DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
    DataBucketGetSizes(db,&n_mp_points,0,0);
	
    for (p=0; p<n_mp_points; p++) {
		int phase_index;
		double *pos;
		MPntStd *material_point;
        
		DataFieldAccessPoint(PField_std,p,(void**)&material_point);
		
		MPntStdGetField_phase_index(material_point,&phase_index);
		MPntStdGetField_global_coord(material_point,&pos);
		if (phase_index<0) {
			PetscPrintf(PETSC_COMM_SELF,"Phase of marker %D is uninitialized. Marker coor (%1.4e,%1.4e,%1.4e)\n",p,pos[0],pos[1],pos[2]);
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Marker phase is uninitialized");
		}
		if (phase_index>=rheology->nphases_active) {
			PetscPrintf(PETSC_COMM_SELF,"Phase of marker %D is larger than rheo->nphases_active = %D. Marker coor (%1.4e,%1.4e,%1.4e)\n",p,rheology->nphases_active, pos[0],pos[1],pos[2]);
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Marker phase is undefined");
		}
	}
    DataFieldRestoreAccess(PField_std);
	
    PetscFunctionReturn(0);
}

/*
 
 These models are required to
 1) set the value c->rheology_constants->nphases_active => although this seems to be done in ModelInitialize_Gene3D()
 2) all markers are assigned a phase index between [0 -- nphases_active-1]
 
 */
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Gene3D"
PetscErrorCode ModelApplyInitialMaterialGeometry_Gene3D(pTatinCtx c,void *ctx)
{
    PetscErrorCode ierr;
    ModelGENE3DCtx *data = (ModelGENE3DCtx *)ctx;

    PetscFunctionBegin;
	
	/* initalize all phase indices to -1 */
	ierr = ModelGene3DInit(c->materialpoint_db);CHKERRQ(ierr);
    switch (data->initial_geom)
	{
            /*Layered cake */
        case 0:
		{
			ierr = ModelSetMarkerIndexLayeredCake_Gene3D(c,ctx);CHKERRQ(ierr);
		}
            break;
            /*Extrude from Map along Z */
        case 1:
		{
			ierr = ModelSetMarkerIndexFromMap_Gene3D(c,ctx);CHKERRQ(ierr);
		}
            break;
            /*Read from CAD file */
        case 2:
		{
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Reading from CAD is not implemented yet \n");
		}
            break;
	}
	/* check all phase indices are between [0---rheo->max_phases-1] */
	ierr = ModelGene3DCheckPhase(c->materialpoint_db,&c->rheology_constants);CHKERRQ(ierr);
	
    ierr = ModelSetInitialStokesVariableOnMarker_Gene3D(c, ctx);CHKERRQ(ierr);
	
    PetscFunctionReturn (0);
}


//======================================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Gene3D"
PetscErrorCode ModelApplyUpdateMeshGeometry_Gene3D(pTatinCtx c,Vec X,void *ctx)
{
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
    PetscPrintf(PETSC_COMM_WORLD, "  NOT IMPLEMENTED \n", __FUNCT__);
	
    PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Gene3D"
PetscErrorCode ModelOutput_Gene3D(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
    ierr = pTatin3d_ModelOutput_MPntStd(c,prefix); CHKERRQ(ierr);
	
    PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Gene3D"
PetscErrorCode ModelDestroy_Gene3D(pTatinCtx c,void *ctx)
{
    ModelGENE3DCtx *data;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
    data = (ModelGENE3DCtx*)ctx;
    
    /* Free contents of structure */
	
    /* Free structure */
    ierr = PetscFree(data);CHKERRQ(ierr);
	
    PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Gene3D"
PetscErrorCode pTatinModelRegister_Gene3D(void)
{
    ModelGENE3DCtx *data;
    pTatinModel    m;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
	
    /* Allocate memory for the data structure for this model */
    ierr = PetscMalloc(sizeof(ModelGENE3DCtx),&data);CHKERRQ(ierr);
    ierr = PetscMemzero(data,sizeof(ModelGENE3DCtx));CHKERRQ(ierr);
	
    /* register user model */
    ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
    /* Set name, model select via -ptatin_model NAME */
    ierr = pTatinModelSetName(m,"Gene3D");CHKERRQ(ierr);
	
    /* Set model data */
    ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
    /* Set function pointers */
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_INIT,                  (void (*)(void)) ModelInitialize_Gene3D); CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_BC,              (void (*)(void)) ModelApplyBoundaryCondition_Gene3D); CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void)) ModelApplyMaterialBoundaryCondition_Gene3D);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void)) ModelApplyInitialMeshGeometry_Gene3D);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void)) ModelApplyInitialMaterialGeometry_Gene3D);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void)) ModelApplyUpdateMeshGeometry_Gene3D);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_OUTPUT,                (void (*)(void)) ModelOutput_Gene3D);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_DESTROY,               (void (*)(void)) ModelDestroy_Gene3D); CHKERRQ(ierr);
	
    /* Insert model into list */
    ierr = pTatinModelRegister(m); CHKERRQ(ierr);
	
    PetscFunctionReturn (0);
}

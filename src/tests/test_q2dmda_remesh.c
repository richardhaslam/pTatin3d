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
 **    filename:   test_q2dmda_remesh.c
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


#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "element_utils_q2.h"
#include "stokes_form_function.h"
#include "stokes_assembly.h"
#include "stokes_operators.h"
#include "mesh_quality_metrics.h"
#include "mesh_deformation.h"
#include "mesh_update.h"

#include "dmda_update_coords.h"
#include "dmda_remesh.h"



PetscErrorCode pTatin3d_remesh(void)
{
	PetscErrorCode  ierr;
	DM              multipys_pack,dav;
	pTatinCtx       user;
	PetscReal       x0,x1,y0,y1,z0,z1;
	PetscInt        M,N,P;
	Vec             X;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
	
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

	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	
	/////////////////
	x0 = y0 = -3.0;
	x1 = y1 = 3.0;
	z0 = 0.0;
	z1 = 1.0;
	ierr = DMDASetUniformCoordinates(dav, x0,x1, y0,y1, z0,z1);CHKERRQ(ierr);
	
	/* output */
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(user,X,"init");CHKERRQ(ierr);
	
	/* deform */
	ierr = MeshDeformation_Sinusodial_ZMAX(dav);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(user,X,"deformed");CHKERRQ(ierr);
	
	/* remesh */
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDARemeshSetUniformCoordinatesBetweenKLayers3d(dav,0,P);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(user,X,"remeshed");CHKERRQ(ierr);
	
	ierr = DMDABilinearizeQ2Elements(dav);CHKERRQ(ierr);
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(user,X,"bilinearized");CHKERRQ(ierr);
	
	/* test metrics */
	{
		PetscReal value;
		
		ierr = DMDAComputeMeshQualityMetric(dav,MESH_QUALITY_ASPECT_RATIO,&value);CHKERRQ(ierr);     PetscPrintf(PETSC_COMM_WORLD,"MESH_QUALITY_ASPECT_RATIO      : %1.4e \n", value);
		ierr = DMDAComputeMeshQualityMetric(dav,MESH_QUALITY_DISTORTION,&value);CHKERRQ(ierr);       PetscPrintf(PETSC_COMM_WORLD,"MESH_QUALITY_DISTORTION        : %1.4e \n", value);
		ierr = DMDAComputeMeshQualityMetric(dav,MESH_QUALITY_DIAGONAL_RATIO,&value);CHKERRQ(ierr);   PetscPrintf(PETSC_COMM_WORLD,"MESH_QUALITY_DIAGONAL_RATIO    : %1.4e \n", value);
        
		ierr = DMDAComputeMeshQualityMetric(dav,MESH_QUALITY_VERTEX_ANGLE,&value);CHKERRQ(ierr);     PetscPrintf(PETSC_COMM_WORLD,"MESH_QUALITY_VERTEX_ANGLE      : %1.4e \n", value);
		ierr = DMDAComputeMeshQualityMetric(dav,MESH_QUALITY_FACE_AREA_RATIO,&value);CHKERRQ(ierr);  PetscPrintf(PETSC_COMM_WORLD,"MESH_QUALITY_FACE_AREA_RATIO   : %1.4e \n", value);	
	}
	
	/////////////////
	
	
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = pTatin3dDestroyContext(&user);
	
	PetscFunctionReturn(0);
}

/*
 ./test_q2dmda_remesh.app -ptatin_model viscous_sinker -mx 64 -my 64 -mz 16
*/
int main( int argc,char **argv )
{
	PetscErrorCode ierr;
	
	ierr = pTatinInitialize(&argc,&argv,(char *)0,NULL);CHKERRQ(ierr);
	
	ierr = pTatin3d_remesh();CHKERRQ(ierr);
	
	ierr = pTatinFinalize();CHKERRQ(ierr);

	return 0;
}

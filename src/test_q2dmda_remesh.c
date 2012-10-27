
#define _GNU_SOURCE

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

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



#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_remesh"
PetscErrorCode pTatin3d_remesh(void)
{
	PetscErrorCode  ierr;
	DM              multipys_pack,dav,dap;
	pTatinCtx       user;
	PetscReal       x0,x1,y0,y1,z0,z1;
	PetscInt        M,N,P;
	Vec             X;
	
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
		
		ierr = DMDAComputeMeshQualityMetrics(dav,MESH_QUALITY_ASPECT_RATIO,&value);CHKERRQ(ierr);     PetscPrintf(PETSC_NULL,"MESH_QUALITY_ASPECT_RATIO    : %1.4e \n", value);
		ierr = DMDAComputeMeshQualityMetrics(dav,MESH_QUALITY_DISTORTION,&value);CHKERRQ(ierr);       PetscPrintf(PETSC_NULL,"MESH_QUALITY_DISTORTION      : %1.4e \n", value);
		ierr = DMDAComputeMeshQualityMetrics(dav,MESH_QUALITY_DIAGONAL_RATIO,&value);CHKERRQ(ierr);   PetscPrintf(PETSC_NULL,"MESH_QUALITY_DIAGONAL_RATIO  : %1.4e \n", value);
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
	
	PetscInitialize(&argc,&argv,(char *)0,0);
	
	ierr = pTatin3d_remesh();CHKERRQ(ierr);
	
	ierr = PetscFinalize();CHKERRQ(ierr);

	return 0;
}

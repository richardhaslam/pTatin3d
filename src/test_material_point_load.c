
static const char help[] = "Stokes solver using Q2-Pm1 mixed finite elements.\n"
"3D prototype of the (p)ragmatic version of Tatin. (pTatin3d_v0.0)\n\n";


#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "material_point_load.h"
#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "element_utils_q2.h"



#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_stokes"
PetscErrorCode pTatin3d_stokes(int argc,char **argv)
{
	PetscErrorCode ierr;
	DM              multipys_pack,dav,dap;
	pTatinCtx       user;

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
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	//ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	// LOAD FROM FILE
	{
		long int N1,N2,N3;
		double *coords_mp;
		double *eta_mp;
		int *phase_mp;
		
		ierr = MarkerCoordinatesLoadFromFile("filters/coords_markers.dat",&N1,&coords_mp);CHKERRQ(ierr);
		ierr = MarkerScalarFieldLoadFromFile("filters/phase_markers.dat",&N2,(void**)&phase_mp);CHKERRQ(ierr);
		ierr = MarkerScalarFieldLoadFromFile("filters/eta_markers.dat",&N3,(void**)&eta_mp);CHKERRQ(ierr);
		printf("N3 = %ld \n",N3);
	}
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);


	/* update markers = >> gauss points */
	{
		int               npoints;
		DataField         PField_std;
		DataField         PField_stokes;
		MPntStd           *mp_std;
		MPntPStokes       *mp_stokes;
		
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
		
		DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
		mp_std    = PField_std->data; /* should write a function to do this */
		mp_stokes = PField_stokes->data; /* should write a function to do this */
		
		ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(npoints,mp_std,mp_stokes,user->stokes_ctx->dav,user->stokes_ctx->volQ);CHKERRQ(ierr);
	}
	

	ierr = pTatin3dDestroyContext(&user);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	
	ierr = PetscInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	ierr = pTatinWriteOptionsFile(PETSC_NULL);CHKERRQ(ierr);
	
	ierr = pTatin3d_stokes(argc,argv);CHKERRQ(ierr);
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}

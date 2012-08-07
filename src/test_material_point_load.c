
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

#include "dmda_element_q2p1.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "material_point_point_location.h"


PetscErrorCode MaterialPointStdRemoval(DataBucket db,long int start,long int npoints,const int wil_key)
{
	DataField         PField_std;
	long int p,end;
	
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	printf("  loaded n_mp_points %ld \n", npoints );
	end = start+npoints;
	for (p=start; p<end; p++) {
		MPntStd     *material_point;

		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		if (material_point->wil==wil_key) {
			DataBucketRemovePointAtIndex(db,p);
			end--;
			p--;
		}
	}	
	DataFieldRestoreAccess(PField_std);
	
	printf("  kept   n_mp_points %ld \n", end-start );
	DataBucketSetSizes(db,end+1,-1); // -ve buffer val retains old value //
	
	PetscFunctionReturn(0);
}

PetscErrorCode MaterialPointStdInsertBasic(DataBucket db,DM da,long int start,long int npoints,double coords_mp[],int phase_mp[])
{
	DM cda;
	Vec gcoords;
	PetscScalar *LA_gcoords;
	double tolerance;
	int p,max_its;
	Truth use_nonzero_guess, monitor, log;
	PetscInt lmx,lmy,lmz;
	PetscInt nel,nen_u;
	MPntStd *local_list;
	const PetscInt *elnidx_u;
	DataField PField_std;
	PetscErrorCode ierr;
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	
	/* point location parameters */
	tolerance         = 1.0e-10;
	max_its           = 10;
	use_nonzero_guess = 0; //_FALSE;
	monitor           = 0; //_FALSE;
	log               = 0; //_FALSE;
	
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	printf("inserting %ld points from index %ld \n", npoints,start);
	for (p=start; p<start+npoints; p++) {
		MPntStd     *material_point;
		MPntStd     test_p;
		
		test_p.pid     = 0;
		test_p.coor[0] = coords_mp[3*p  ];
		test_p.coor[1] = coords_mp[3*p+1];
		test_p.coor[2] = coords_mp[3*p+2];
		test_p.phase   = phase_mp[p];
		test_p.wil     = -1;
		
		DataFieldInsertPoint(PField_std,p,(void*)&test_p);
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		
		InverseMappingDomain_3dQ2(tolerance, max_its,
															use_nonzero_guess, 
															monitor, log,
															(const double*)LA_gcoords, (const int)lmx,(const int)lmy,(const int)lmz, (const int*)elnidx_u,
															1, material_point );
	}
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	DataFieldRestoreAccess(PField_std);
		
	PetscFunctionReturn(0);
}


PetscErrorCode MaterialPointDataBasicLoadIntoListFromFile(DataBucket db,DM da,PetscBool append,const char coordfile[],const char phasefile[])
{
	long int N1,N2;
	double *coords_mp;
	int *phase_mp;
	long int p,start;
	int n_mp_points;
	PetscErrorCode ierr;
	
	/* read in from file */
	printf("reading files %s : %s \n", coordfile,phasefile);
	ierr = MarkerCoordinatesLoadFromFile(coordfile,&N1,&coords_mp);CHKERRQ(ierr);
	ierr = MarkerScalarFieldLoadFromFile(phasefile,&N2,(void**)&phase_mp);CHKERRQ(ierr);
	if (N1 != N2) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Coordinate file and marker file have different lengths");
	}

	if (append == PETSC_FALSE) {
		start = 0;
		DataBucketSetSizes(db,N1,-1); // -ve buffer val retains old value //
	} else {
		DataBucketGetSizes(db,&n_mp_points,0,0);
		start = n_mp_points - 1;
		n_mp_points = n_mp_points + N1;
		DataBucketSetSizes(db,n_mp_points,-1); // -ve buffer val retains old value //
	}
	
	ierr = MaterialPointStdInsertBasic(db,da,start,N1,coords_mp,phase_mp);CHKERRQ(ierr);
	ierr = MaterialPointStdRemoval(db,start,N1,-1);CHKERRQ(ierr);

	DataBucketGetSizes(db,&n_mp_points,0,0);
	ierr = SwarmMPntStd_AssignUniquePointIdentifiers(((PetscObject)da)->comm,db,0,n_mp_points);CHKERRQ(ierr);

	
	free(coords_mp);
	free(phase_mp);

	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_stokes"
PetscErrorCode pTatin3d_stokes(int argc,char **argv)
{
	PetscErrorCode ierr;
	DM              multipys_pack,dav,dap;
	pTatinCtx       user;
	Vec X;

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
	
	ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);

	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
#if 0
	// LOAD FROM FILE
	{
		long int N1,N2,N3;
		double *coords_mp;
		double *eta_mp;
		int *phase_mp;
		int p,n_mp_points;

		DataBucket db;
		DataField PField_std;

		/* read in from file */
		ierr = MarkerCoordinatesLoadFromFile("filters/coords_markers.dat",&N1,&coords_mp);CHKERRQ(ierr);
		ierr = MarkerScalarFieldLoadFromFile("filters/phase_markers.dat",&N2,(void**)&phase_mp);CHKERRQ(ierr);
		//ierr = MarkerScalarFieldLoadFromFile("filters/eta_markers.dat",&N3,(void**)&eta_mp);CHKERRQ(ierr);
		//printf("N3 = %ld \n",N3);

		db = user->materialpoint_db;
		DataBucketSetSizes(db,N1,-1); // -ve buffer val retains old value //

		
		/* update local coordinates and wil */
		{
			DM da,cda;
			Vec gcoords;
			PetscScalar *LA_gcoords;
			double tolerance;
			int max_its;
			Truth use_nonzero_guess, monitor, log;
			PetscInt lmx,lmy,lmz;
			PetscInt nel,nen_u;
			MPntStd *local_list;
			const PetscInt *elnidx_u;
			DataBucket db;
			DataField PField_std;
			
			/* setup for coords */
			da = dav;
			ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
			ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
			ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
			
			ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
			
			ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,&lmz);CHKERRQ(ierr);
			
			/* point location parameters */
			tolerance         = 1.0e-10;
			max_its           = 10;
			use_nonzero_guess = 0; //_FALSE;
			monitor           = 0; //_FALSE;
			log               = 0; //_FALSE;

			
			db = user->materialpoint_db;
			DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
			DataFieldGetAccess(PField_std);
			DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
			DataBucketGetSizes(db,&n_mp_points,0,0);

			for (p=0; p<n_mp_points; p++) {
				MPntStd     *material_point;
				MPntStd     test_p;
				
				test_p.pid = 0;
				test_p.coor[0] = coords_mp[3*p  ];
				test_p.coor[1] = coords_mp[3*p+1];
				test_p.coor[2] = coords_mp[3*p+2];
				test_p.phase = phase_mp[p];
				test_p.wil = -1;
				
				DataFieldInsertPoint(PField_std,p,(void*)&test_p);
				
				DataFieldAccessPoint(PField_std,p,   (void**)&material_point);

				InverseMappingDomain_3dQ2(tolerance, max_its,
																	use_nonzero_guess, 
																	monitor, log,
																	(const double*)LA_gcoords, (const int)lmx,(const int)lmy,(const int)lmz, (const int*)elnidx_u,
																	1, material_point );
			}
			
			ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
			DataFieldRestoreAccess(PField_std);
		}
		
		/* remove points not located */
		DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
		DataFieldGetAccess(PField_std);
		DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
		DataBucketGetSizes(db,&n_mp_points,0,0);
		
		printf("loaded n_mp_points %d \n", n_mp_points );
		for (p=0; p<n_mp_points; p++) {
			MPntStd     *material_point;
			DataFieldAccessPoint(PField_std,p,   (void**)&material_point);

			if (material_point->wil==-1) {
				DataBucketRemovePointAtIndex(db,p);
				n_mp_points--;
				p--;
			}
		}	
		DataFieldRestoreAccess(PField_std);

		printf("kept n_mp_points %d \n", n_mp_points );
		DataBucketSetSizes(db,n_mp_points,-1); // -ve buffer val retains old value //
		DataBucketGetSizes(db,&n_mp_points,0,0);
		//MPI_Allreduce(&npoints,&npoints_global_fin,1,MPI_INT,MPI_SUM,de->comm);
	}
#endif

//	ierr = MaterialPointDataBasicLoadIntoListFromFile(user->materialpoint_db,dav,PETSC_FALSE,"filters/coords_markers.dat","filters/phase_markers.dat");CHKERRQ(ierr);

//	ierr = MaterialPointDataBasicLoadIntoListFromFile(user->materialpoint_db,dav,PETSC_FALSE,"testdump/coords_test-0.dat","testdump/phase_test-0.dat");CHKERRQ(ierr);
//	ierr = MaterialPointDataBasicLoadIntoListFromFile(user->materialpoint_db,dav,PETSC_TRUE, "testdump/coords_test-2.dat","testdump/phase_test-2.dat");CHKERRQ(ierr);
//	ierr = MaterialPointDataBasicLoadIntoListFromFile(user->materialpoint_db,dav,PETSC_TRUE, "testdump/coords_test-4.dat","testdump/phase_test-4.dat");CHKERRQ(ierr);
	
	
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);


	/* update markers = >> gauss points */
#if 0	
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
#endif
	
	/* boundary conditions */
	ierr = pTatinModel_Output(user->model,user,X,"test");CHKERRQ(ierr);
	

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

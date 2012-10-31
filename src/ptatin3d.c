
#include "private/daimpl.h" 
#include "petsc.h"

#include "ptatin3d.h"
#include "ptatin3d_defs.h"

#include "element_type_Q2.h"
#include "dmda_element_q2p1.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "material_point_std_utils.h"
#include "ptatin_utils.h"
#include "ptatin3d_stokes.h"
#include "output_paraview.h"
#include "ptatin_models.h"
#include "dmda_checkpoint.h"

#include "private/ptatin_impl.h"


/* PHYSICS MESHES */

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_PhysCompStokesNew"
PetscErrorCode pTatin3d_PhysCompStokesNew(pTatinCtx user)
{
	PetscErrorCode  ierr;
	PhysCompStokes stokes;
	
	
	PetscFunctionBegin;
	ierr = PhysCompCreate_Stokes(&stokes);CHKERRQ(ierr);
	
	stokes->mx = user->mx;
	stokes->my = user->my;
	stokes->mz = user->mz;
	stokes->use_mf_stokes = user->use_mf_stokes;
	
	ierr = PhysCompCreateMesh_Stokes3d(stokes->mx,stokes->my,stokes->mz,stokes);CHKERRQ(ierr);
	ierr = PhysCompCreateBoundaryList_Stokes(stokes);CHKERRQ(ierr);
	ierr = PhysCompCreateVolumeQuadrature_Stokes(stokes);CHKERRQ(ierr);
	//	ierr = PhysCompCreateSurfaceQuadrature_Stokes(stokes);CHKERRQ(ierr);
	
	user->stokes_ctx = stokes;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_PhysCompStokesCreate"
PetscErrorCode pTatin3d_PhysCompStokesCreate(pTatinCtx user)
{
	PetscErrorCode  ierr;
	PhysCompStokes stokes;
	
	
	PetscFunctionBegin;

	if (user->restart_from_file) {
		/* load from file */
		char vname[PETSC_MAX_PATH_LEN];
		char pname[PETSC_MAX_PATH_LEN];		
		
		/* dav,dap */
		if (!StringEmpty(user->restart_prefix)) {
			sprintf(vname,"%s/ptat3dcpf.dmda-velocity_%s",user->restart_dir,user->restart_prefix);
			sprintf(pname,"%s/ptat3dcpf.dmda-pressure_%s",user->restart_dir,user->restart_prefix);
		} else {
			sprintf(vname,"%s/ptat3dcpf.dmda-velocity",user->restart_dir);
			sprintf(pname,"%s/ptat3dcpf.dmda-pressure",user->restart_dir);
		}
		PetscPrintf(PETSC_COMM_WORLD,"  reading %s \n", vname );
		PetscPrintf(PETSC_COMM_WORLD,"  reading %s \n", pname );
		
		ierr = pTatin3d_PhysCompStokesLoad(user,(const char*)vname,(const char*)pname);CHKERRQ(ierr);
	} else {
		/* create from data */
		ierr = pTatin3d_PhysCompStokesNew(user);CHKERRQ(ierr);
	}	
	
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_ModelOutput_VelocityPressure_Stokes"
PetscErrorCode pTatin3d_ModelOutput_VelocityPressure_Stokes(pTatinCtx ctx,Vec X,const char prefix[])
{
	PetscErrorCode ierr;
	char *name;
	DM stokes_pack;
	Vec UP;
	PetscLogDouble t0,t1;
	static int beenhere=0;
	static char *pvdfilename;
	PetscFunctionBegin;
	
	PetscGetTime(&t0);
	// PVD
	if (beenhere==0) {
		asprintf(&pvdfilename,"%s/timeseries_vp.pvd",ctx->outputpath);
		PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename );
		ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);

		beenhere = 1;
	}
	{
		char *vtkfilename;
		
		if (prefix) {
			asprintf(&vtkfilename, "%s_vp.pvts",prefix);
		} else {
			asprintf(&vtkfilename, "vp.pvts");
		}
		
		ierr = ParaviewPVDAppend(pvdfilename,ctx->time, vtkfilename, "");CHKERRQ(ierr);
		free(vtkfilename);
	}
	
	// PVTS + VTS
	if (prefix) {
		asprintf(&name,"%s_vp",prefix);
	} else {
		asprintf(&name,"vp");
	}
	
	PetscPrintf(PETSC_COMM_WORLD,"[[DESIGN FLAW]] %s: require better physics modularity to extract (u,p) <---| (X) \n", __FUNCT__ );
	
	stokes_pack = ctx->stokes_ctx->stokes_pack;
	UP = X;
	ierr = pTatinOutputParaViewMeshVelocityPressure(stokes_pack,UP,ctx->outputpath,name);CHKERRQ(ierr);
	free(name);
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_vp.(pvd,pvts,vts): CPU time %1.2e (sec) \n", __FUNCT__,prefix,t1-t0);
	
	PetscFunctionReturn(0);
}



/* MATERIAL POINTS */
#undef __FUNCT__
#define __FUNCT__ "SwarmDMDA3dDataExchangerCreate"
PetscErrorCode SwarmDMDA3dDataExchangerCreate(DM da,DataEx *_de)
{
	DataEx de;
	const PetscInt *neighborranks;
	PetscInt neighborranks2[27],neighborcount;
	PetscInt i;
	PetscLogDouble t0,t1;
	PetscMPIInt rank;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	MPI_Comm_rank(((PetscObject)da)->comm,&rank);
	ierr = DMDAGetNeighbors(da,&neighborranks);CHKERRQ(ierr);
	
	neighborcount = 0;
	for (i=0; i<27; i++) {
		neighborranks2[i] = -1;
		if ( (neighborranks[i]>=0) && (neighborranks[i]!=rank) ) {
			neighborranks2[neighborcount] = neighborranks[i];
			neighborcount++;
		}
	}
	
	PetscGetTime(&t0);
	de = DataExCreate(((PetscObject)da)->comm,0);
	//	de = DataExCreate(PETSC_COMM_WORLD,0);
	ierr = DataExTopologyInitialize(de);CHKERRQ(ierr);
	for (i=0; i<neighborcount; i++) {
		ierr = DataExTopologyAddNeighbour(de,neighborranks2[i]);CHKERRQ(ierr);
	}
	ierr = DataExTopologyFinalize(de);CHKERRQ(ierr);
	PetscGetTime(&t1);
	PetscPrintf(de->comm,"DataEx Communication setup time = %1.4e\n",t1-t0);
	
	*_de = de;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dCreateMaterialPoints"
PetscErrorCode pTatin3dCreateMaterialPoints(pTatinCtx ctx,DM dav)
{
	DataBucket     db;
	DataEx         ex;	
	PetscLogDouble t0,t1;
	PetscInt       lmx,lmy,lmz;
	PetscBool      flg;
	int            npoints;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* register marker structures here */
	PetscGetTime(&t0);
	DataBucketCreate(&db);
	DataBucketRegisterField(db,MPntStd_classname,    sizeof(MPntStd),PETSC_NULL);
	DataBucketRegisterField(db,MPntPStokes_classname,sizeof(MPntPStokes),PETSC_NULL);
	DataBucketFinalize(db);
	
	/* Choose type of projection (for eta and rho) */
	ctx->coefficient_projection_type = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-coefficient_projection_type",&ctx->coefficient_projection_type,&flg);CHKERRQ(ierr);
	switch (ctx->coefficient_projection_type) {
		case -1:
			PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using null projection\n");
			break;
		case 0:
			PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection\n");
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 not implemented");
			break;
		case 1:
			PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using Q1 projection\n");
			break;
		case 2:
			PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using Q2 projection\n");
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = Q2 not implemented");
			break;
		case 3:
			PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P1 projection\n");
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P1 not implemented");
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER," -coefficient_projection_type = {0,1,2} implying {P0,Q1,Q2}");
			break;
	}
	
	/* Choose coordinate layout for material points, -mp_layout 0,1 */
	/* set initial size */
  ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	DataBucketSetInitialSizes(db,lmx*lmy*lmz*2*2*2,1000);
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"[[Swarm initialization: %1.4lf sec]]\n", t1-t0);
	
	PetscGetTime(&t0);
	{
		/* defaults for lattice layout */
		PetscInt   Nxp[] = {2,2,2}; /* change with -lattice_layout_N{x,y,z} */
		PetscReal  perturb = 0.1;   /* change with -lattice_layout_perturb */
		/* defaults for random layout */
		PetscInt   nPerCell = 3*3*3; /* change with -random_layout_Np */
		
		PetscInt mplayout = 0;
		
		PetscOptionsGetInt(PETSC_NULL,"-mp_layout",&mplayout,PETSC_NULL);
		if (mplayout==0) {
			ierr = SwarmMPntStd_CoordAssignment_LatticeLayout3d(dav,Nxp,perturb,db);CHKERRQ(ierr);
		} else if (mplayout==1) {
			ierr = SwarmMPntStd_CoordAssignment_RandomLayout3d(dav,nPerCell,db);CHKERRQ(ierr);
		}
	}
	PetscGetTime(&t1);
	DataBucketGetSizes(db,&npoints,PETSC_NULL,PETSC_NULL);
	PetscPrintf(PETSC_COMM_WORLD,"[[Swarm->coordinate assignment: %d points : %1.4lf sec]]\n", npoints,t1-t0);
	
	
	/* create the data exchanger need for parallel particle movement */
	ierr = SwarmDMDA3dDataExchangerCreate(dav,&ex);CHKERRQ(ierr);
	
	ctx->materialpoint_db = db;
	ctx->materialpoint_ex = ex;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointCoordinateSetUp"
PetscErrorCode MaterialPointCoordinateSetUp(pTatinCtx ctx,DM da)
{
	PetscErrorCode ierr;
	DataBucket db;
	DM cda;
	Vec gcoords;
	PetscScalar *LA_gcoords;
	PetscInt nel,nen,e,i,p,n_mp_points;
	const PetscInt *elnidx;
	PetscScalar Ni_p[Q2_NODES_PER_EL_3D];
	PetscScalar elcoords[3*Q2_NODES_PER_EL_3D];
	DataField      PField_std;
	
	
  PetscFunctionBegin;
	db = ctx->materialpoint_db;
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		double      xp[3];
		double      *xi;
		
		DataFieldAccessPoint(PField_std,p,(void**)&material_point);
		
		MPntStdGetField_local_element_index(material_point,&e);
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		
		MPntStdGetField_local_coord(material_point,&xi);
		pTatin_ConstructNi_Q2_3D(xi,Ni_p);
		
		xp[0] = 0.0;
		xp[1] = 0.0;
		xp[2] = 0.0;
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			xp[0] = xp[0] + Ni_p[i] * elcoords[NSD*i  ];
			xp[1] = xp[1] + Ni_p[i] * elcoords[NSD*i+1];
			xp[2] = xp[2] + Ni_p[i] * elcoords[NSD*i+2];
		}
		
		MPntStdSetField_global_coord(material_point,xp);
		
	}
	
	DataFieldRestoreAccess(PField_std);
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_ModelOutput_MPntStd"
PetscErrorCode pTatin3d_ModelOutput_MPntStd(pTatinCtx ctx,const char prefix[])
{
	PetscErrorCode ierr;
	char *name;
	DM stokes_pack;
	Vec UP;
	PetscLogDouble t0,t1;
	static int beenhere=0;
	static char *pvdfilename;
	
	PetscFunctionBegin;
	
	PetscGetTime(&t0);
	// PVD
	if (beenhere==0) {
		asprintf(&pvdfilename,"%s/timeseries_mpoints_std.pvd",ctx->outputpath);
		PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename );
		ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
		
		beenhere = 1;
	}
	{
		char *vtkfilename;
		
		if (prefix) {
			asprintf(&vtkfilename, "%s_mpoints_std.pvtu",prefix);
		} else {
			asprintf(&vtkfilename, "mpoints_std.pvtu");
		}
		
		ierr = ParaviewPVDAppend(pvdfilename,ctx->time, vtkfilename, "");CHKERRQ(ierr);
		free(vtkfilename);
	}
	
	// PVTS + VTS
	if (prefix) {
		asprintf(&name,"%s_mpoints_std",prefix);
	} else {
		asprintf(&name,"mpoints_std");
	}
	
	ierr = SwarmOutputParaView_MPntStd(ctx->materialpoint_db,ctx->outputpath,name);CHKERRQ(ierr);
	free(name);
	
	PetscGetTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_mpoints_std.(pvd,pvtu,vtu): CPU time %1.2e (sec) \n", __FUNCT__,prefix,t1-t0);

	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dCreateContext"
PetscErrorCode pTatin3dCreateContext(pTatinCtx *ctx)
{
	PetscInt       e;
	pTatinCtx      user;
	PetscMPIInt    rank;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
  ierr = PetscNew(struct _p_pTatinCtx,&user);CHKERRQ(ierr);
	
	
	/* init */
	user->stokes_ctx = PETSC_NULL;
//	user->energy_ctx = PETSC_NULL;
//	user->coords_ctx = PETSC_NULL;
  
	user->pack     = PETSC_NULL; /* DM composite for velocity and pressure */
	
	/* set defaults */
	user->restart_from_file         = PETSC_FALSE;
	user->checkpoint_every          = 1000;
	user->checkpoint_every_nsteps   = 1000;
	user->checkpoint_every_ncpumins = 90.0;
		
	user->mx               = 4;
	user->my               = 4;
	user->mz               = 4;
	user->use_mf_stokes    = PETSC_FALSE;
	user->solverstatistics = PETSC_FALSE;
	
	user->continuation_m   = 1;
	user->continuation_M   = 1;
	user->coefficient_projection_type = 1; /* Q1 */
	
	/* time step control */
	user->nsteps           = 1;
	user->dt_max           = 1.0e30;
	user->dt_min           = 1.0e-30;
	user->dt               = user->dt_min;
	user->output_frequency = 1;
	user->time_max         = 1.0e32;
	user->time             = 0.0;
  user->step             = 0;
  user->dt_adv           = user->dt_min;
  
  ierr = RheologyConstantsInitialise(&user->rheology_constants);CHKERRQ(ierr);
	ierr = MaterialConstantsInitialize(&user->material_constants);CHKERRQ(ierr);
	ierr = PetscContainerCreate(PETSC_COMM_WORLD,&user->model_data);CHKERRQ(ierr);
  
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	if (rank==0) {
		pTatinGenerateFormattedTimestamp(user->formatted_timestamp);	
	}
	ierr = MPI_Bcast(user->formatted_timestamp,PETSC_MAX_PATH_LEN,MPI_INT,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	*ctx = user;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dDestroyContext"
PetscErrorCode pTatin3dDestroyContext(pTatinCtx *ctx)
{
	pTatinCtx user = *ctx;
	PetscInt e;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	if (user->materialpoint_ex) { /* ierr = DataExView(user->materialpoint_ex);CHKERRQ(ierr); */ ierr = DataExDestroy(user->materialpoint_ex);CHKERRQ(ierr); }
	if (user->materialpoint_db) { DataBucketDestroy(&user->materialpoint_db); }
	
	if (user->stokes_ctx) { ierr = PhysCompDestroy_Stokes(&user->stokes_ctx);CHKERRQ(ierr); }
	if (user->pack) {       ierr = DMDestroy(&user->pack);CHKERRQ(ierr); }
	
	/*
	 if (user->Q) { ierr = QuadratureStokesDestroy(&user->Q);CHKERRQ(ierr); }
	 for (e=0; e<QUAD_EDGES; e++) {
	 if (user->surfQ[e]) { ierr = SurfaceQuadratureStokesDestroy(&user->surfQ[e]);CHKERRQ(ierr); }
	 }
	 if (user->p_bclist) { ierr = BCListDestroy(&user->p_bclist);CHKERRQ(ierr); }
	 if (user->u_bclist) { ierr = BCListDestroy(&user->u_bclist);CHKERRQ(ierr); }
	 if (user->dap) { ierr = DMDestroy(&user->dap);CHKERRQ(ierr); }
	 if (user->dav) { ierr = DMDestroy(&user->dav);CHKERRQ(ierr); }
	 */

	if (user->material_constants) { DataBucketDestroy(&user->material_constants); }
	
	ierr = PetscContainerDestroy(&user->model_data);CHKERRQ(ierr);

	{
		char  logfile[PETSC_MAX_PATH_LEN];
		PetscViewer viewer;
		
		sprintf(logfile,"%s/ptatin.petsc.log_summary-%s",user->outputpath,user->formatted_timestamp);

		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,logfile,&viewer);CHKERRQ(ierr);
		ierr = PetscLogView(viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	}
		
	ierr = PetscFree(user);CHKERRQ(ierr);
	
	*ctx = PETSC_NULL;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinCtxGetModelData"
PetscErrorCode pTatinCtxGetModelData(pTatinCtx ctx,const char name[],void **data)
{
	PetscErrorCode ierr;
	PetscContainer container;
	
  PetscFunctionBegin;
	ierr = PetscObjectQuery((PetscObject)ctx->model_data,name,(PetscObject*)&container);CHKERRQ(ierr);
	if (!container) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"No data with name \"%s\" was composed with ctx->model_data",name);
	ierr = PetscContainerGetPointer(container,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinCtxAttachModelData"
PetscErrorCode pTatinCtxAttachModelData(pTatinCtx ctx,const char name[],void *data)
{
	PetscContainer container;
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
  ierr = PetscContainerCreate(PETSC_COMM_WORLD,&container);CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(container,(void*)data);CHKERRQ(ierr);
	
	ierr = PetscObjectCompose((PetscObject)ctx->model_data,name,(PetscObject)container);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dParseOptions"
PetscErrorCode pTatin3dParseOptions(pTatinCtx ctx)
{
  char           optionsfile[PETSC_MAX_PATH_LEN];
	PetscInt       mx,my,mz;
	PetscBool      flg;
	PetscErrorCode ierr;
	
	/* parse options */
	PetscOptionsGetInt(PETSC_NULL,"-mx",&ctx->mx,0);
	PetscOptionsGetInt(PETSC_NULL,"-my",&ctx->my,0);
	PetscOptionsGetInt(PETSC_NULL,"-mz",&ctx->mz,0);
	
	ierr = PetscOptionsGetBool(PETSC_NULL,"-use_mf_stokes",&ctx->use_mf_stokes,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(PETSC_NULL,"-with_statistics",&ctx->solverstatistics,&flg);CHKERRQ(ierr);
	
	flg = PETSC_FALSE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-output_path",ctx->outputpath,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg == PETSC_FALSE) { 
		sprintf(ctx->outputpath,"./output");
	}
	ierr = pTatinCreateDirectory(ctx->outputpath);CHKERRQ(ierr);
	
	/* checkpointing */
	PetscOptionsGetInt(PETSC_NULL,"-checkpoint_every",&ctx->checkpoint_every,&flg);
	PetscOptionsGetInt(PETSC_NULL,"-checkpoint_every_nsteps",&ctx->checkpoint_every_nsteps,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-checkpoint_every_ncpumins",&ctx->checkpoint_every_ncpumins,&flg);
	
	/* time stepping */
	PetscOptionsGetInt(PETSC_NULL,"-nsteps",&ctx->nsteps,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-dt_min",&ctx->dt_min,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-dt_max",&ctx->dt_max,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-time_max",&ctx->time_max,&flg);
	PetscOptionsGetInt(PETSC_NULL,"-output_frequency",&ctx->output_frequency,&flg);
	
	sprintf(optionsfile,"%s/ptatin.options-%s",ctx->outputpath,ctx->formatted_timestamp);
	ierr = pTatinWriteOptionsFile(optionsfile);CHKERRQ(ierr);

	sprintf(optionsfile,"%s/ptatin.options",ctx->outputpath);
	ierr = pTatinWriteOptionsFile(optionsfile);CHKERRQ(ierr);
	
//	ierr = pTatinModelLoad(ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelLoad"
PetscErrorCode pTatinModelLoad(pTatinCtx ctx)
{
	pTatinModel model;
	PetscBool flgname;
	char modelname[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;

	PetscFunctionBegin;
	
	flgname = PETSC_FALSE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-ptatin_model",modelname,PETSC_MAX_PATH_LEN-1,&flgname);CHKERRQ(ierr);
	if (flgname) {
		ierr = pTatinModelGetByName(modelname,&model);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  [pTatinModel]: -ptatin_model \"%s\" was detected\n",model->model_name);
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"  [pTatinModel]: -ptatin_model wasn't specified - running boring \"template\" model\n");
		ierr = pTatinModelGetByName("template",&model);CHKERRQ(ierr);
	}
	
	model->ptat_ctx = ctx;
	ctx->model = model;
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "pTatinGetMaterialPoints"
PetscErrorCode pTatinGetMaterialPoints(pTatinCtx ctx,DataBucket *db,DataEx *de)
{
	PetscErrorCode ierr;
	if (db) { *db = ctx->materialpoint_db; }
	if (de) { *de = ctx->materialpoint_ex; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinGetMaterialConstants"
PetscErrorCode pTatinGetMaterialConstants(pTatinCtx ctx,DataBucket *db)
{
	PetscErrorCode ierr;
	if (db) { *db = ctx->material_constants; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinGetModel"
PetscErrorCode pTatinGetModel(pTatinCtx ctx,pTatinModel *m)
{
	PetscErrorCode ierr;
	if (m) { *m = ctx->model; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinGetRheology"
PetscErrorCode pTatinGetRheology(pTatinCtx ctx,RheologyConstants **r)
{
	PetscErrorCode ierr;
	if (r) { *r = &ctx->rheology_constants; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinGetStokesContext"
PetscErrorCode pTatinGetStokesContext(pTatinCtx ctx,PhysCompStokes *s)
{
	PetscErrorCode ierr;
	if (s) { *s = ctx->stokes_ctx; }
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dContextLoad"
PetscErrorCode pTatin3dContextLoad(pTatinCtx *ctx,const char filename[])
{
	PetscViewer viewer;
	pTatinCtx cc;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	ierr = pTatin3dCreateContext(&cc);CHKERRQ(ierr);
	ierr = PetscMemzero(cc,sizeof(struct _p_pTatinCtx));CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
	
	ierr = PetscViewerBinaryRead(viewer,cc,sizeof(struct _p_pTatinCtx)/sizeof(char),PETSC_CHAR);
	
	ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	 
	/* zero out any pointers */
	cc->stokes_ctx = PETSC_NULL;
	cc->materialpoint_db = PETSC_NULL;
	cc->materialpoint_ex = PETSC_NULL;
	cc->material_constants = PETSC_NULL;
	cc->model = PETSC_NULL;
	cc->model_data = PETSC_NULL;
	
	*ctx = cc;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dContextSave"
PetscErrorCode pTatin3dContextSave(pTatinCtx ctx,const char filename[])
{
	PetscViewer viewer;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
	
	ierr = PetscViewerBinaryWrite(viewer,ctx,sizeof(struct _p_pTatinCtx)/sizeof(char),PETSC_CHAR,PETSC_FALSE);
	
	ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	
	 
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_PhysCompStokesLoad"
PetscErrorCode pTatin3d_PhysCompStokesLoad(pTatinCtx user,const char vname[],const char pname[])
{
	PhysCompStokes stokes;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PhysCompCreate_Stokes(&stokes);CHKERRQ(ierr);
	
	stokes->mx = user->mx;
	stokes->my = user->my;
	stokes->mz = user->mz;
	stokes->use_mf_stokes = user->use_mf_stokes;
	
	ierr = PhysCompLoadMesh_Stokes3d(stokes,vname,pname);CHKERRQ(ierr);
//	ierr = PhysCompCreateMesh_Stokes3d(stokes->mx,stokes->my,stokes->mz,stokes);CHKERRQ(ierr);
	
	ierr = PhysCompCreateBoundaryList_Stokes(stokes);CHKERRQ(ierr);
	ierr = PhysCompCreateVolumeQuadrature_Stokes(stokes);CHKERRQ(ierr);
	//	ierr = PhysCompCreateSurfaceQuadrature_Stokes(stokes);CHKERRQ(ierr);
	
	user->stokes_ctx = stokes;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dCheckpoint"
PetscErrorCode pTatin3dCheckpoint(pTatinCtx ctx,Vec X,const char prefix[])
{
	PetscViewer viewer;
	PetscErrorCode ierr;
	char start[PETSC_MAX_PATH_LEN];
	char f1[PETSC_MAX_PATH_LEN];
	char f2[PETSC_MAX_PATH_LEN];
	char f3[PETSC_MAX_PATH_LEN];
	
	PetscFunctionBegin;

	
	/* context */
	if (prefix) {
		sprintf(start,"%s/ptat3dcpf.ctx_%s",ctx->outputpath,prefix);
	} else {
		sprintf(start,"%s/ptat3dcpf.ctx",ctx->outputpath);
	}
	PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", start );
	ierr = pTatin3dContextSave(ctx,start);CHKERRQ(ierr);	
	
	/* dav,dap */
	if (prefix) {
		sprintf(f1,"%s/ptat3dcpf.dmda-velocity_%s",ctx->outputpath,prefix);
		sprintf(f2,"%s/ptat3dcpf.dmda-velocity-coords_%s",ctx->outputpath,prefix);
		sprintf(f3,"%s/ptat3dcpf.dmda-pressure_%s",ctx->outputpath,prefix);
	} else {
		sprintf(f1,"%s/ptat3dcpf.dmda-velocity",ctx->outputpath);
		sprintf(f2,"%s/ptat3dcpf.dmda-velocity-coords",ctx->outputpath);
		sprintf(f3,"%s/ptat3dcpf.dmda-pressure",ctx->outputpath);
	}
	PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", f1 );
	PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", f2 );
	PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", f3 );
	
	ierr = PhysCompSaveMesh_Stokes3d(ctx->stokes_ctx,f1,f3,f2);CHKERRQ(ierr);

	/* solution */
	if (prefix) {
		sprintf(f1,"%s/ptat3dcpf.dmda-X_%s",ctx->outputpath,prefix);
	} else {
		sprintf(f1,"%s/ptat3dcpf.dmda-X",ctx->outputpath);
	}
	PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", f1 );
	{
		PetscViewer viewer;
		
		ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD,f1,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
		ierr = VecView(X,viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	}
	
	/* material points */
	if (prefix) {
		sprintf(start,"%s/ptat3dcpf.markers_%s",ctx->outputpath,prefix);
	} else {
		sprintf(start,"%s/ptat3dcpf.markers",ctx->outputpath);
	}
	PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", start );
	DataBucketView(PETSC_COMM_WORLD,ctx->materialpoint_db,start,DATABUCKET_VIEW_BINARY);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dCheckpointManager"
PetscErrorCode pTatin3dCheckpointManager(pTatinCtx ctx,Vec X)
{
	PetscErrorCode ierr;
	PetscInt       checkpoint_every;
	PetscInt       checkpoint_every_nsteps,step;
	double         checkpoint_every_ncpumins, max_current_cpu_time, current_cpu_time;
	static double  last_cpu_time = 0.0;
	int            exists;
	char           prefix[PETSC_MAX_PATH_LEN],filetocheck[PETSC_MAX_PATH_LEN];
	
	PetscFunctionBegin;

	step                      = ctx->step;
	checkpoint_every          = ctx->checkpoint_every;
	checkpoint_every_nsteps   = ctx->checkpoint_every_nsteps;
	checkpoint_every_ncpumins = ctx->checkpoint_every_ncpumins;
	
	sprintf(prefix,"step%1.6d",step);

	
	/* -------------------------------------- */
	/* check one - this file has a fixed name */
	if (step%checkpoint_every==0) {
		char command[256];
		char file[256];
		

		sprintf(filetocheck,"%s/ptat3dcpf.ctx",ctx->outputpath);
		PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager: Writing\n");
		// call checkpoint routine //
		ierr = pTatin3dCheckpoint(ctx,X,PETSC_NULL);CHKERRQ(ierr);
	}
	
	/* -------------------------------------------------------------------- */
	/* check three - look at cpu time and decide if we need to write or not */
	PetscGetTime(&current_cpu_time);
	ierr = MPI_Allreduce(&current_cpu_time,&max_current_cpu_time,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	if (current_cpu_time > last_cpu_time + checkpoint_every_ncpumins) {
		char command[256];
		char file[256];
		
		sprintf(filetocheck,"%s/ptat3dcpf.ctx_step%1.6d",ctx->outputpath,step);
		FileExists(filetocheck,&exists);
		//PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager[checkpoint_every_ncpumins]: Checking for files %s\n",filetocheck);
		if (exists==0) {
			PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager[checkpoint_every_ncpumins]: Writing\n");
			// call checkpoint routine //
			ierr = pTatin3dCheckpoint(ctx,X,prefix);CHKERRQ(ierr);
		}
		
		last_cpu_time = max_current_cpu_time;
	}

	
	/* ----------------------------------------------------------------- */
	/* check two - these files have a file name related to the time step */
	if (step%checkpoint_every_nsteps==0) {
		char command[256];
		char file[256];
		
		sprintf(filetocheck,"%s/ptat3dcpf.ctx_step%1.6d",ctx->outputpath,step);
		FileExists(filetocheck,&exists);
		//PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager[checkpoint_every_nsteps]: Checking for files %s\n",filetocheck);
		if (exists==0) {
			PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager[checkpoint_every_nsteps]: Writing\n");
			// call checkpoint routine //
			ierr = pTatin3dCheckpoint(ctx,X,prefix);CHKERRQ(ierr);
		}
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinRestart_Initialize"
PetscErrorCode pTatinRestart_Initialize(pTatinCtx ctx,void *data)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinRestart_ApplyInitialMeshGeometry"
PetscErrorCode pTatinRestart_ApplyInitialMeshGeometry(pTatinCtx ctx,void *data)
{
	char name[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscFunctionBegin;

	if (!StringEmpty(ctx->restart_prefix)) {
		sprintf(name,"%s/ptat3dcpf.dmda-velocity-coords_%s",ctx->restart_dir,ctx->restart_prefix);
	} else {
		sprintf(name,"%s/ptat3dcpf.dmda-velocity-coords",ctx->restart_dir);
	}

	PetscPrintf(PETSC_COMM_WORLD,"Restart: [ApplyInitialMeshGeometry] from %s \n", name );
	ierr = DMDALoadCoordinatesFromFile(ctx->stokes_ctx->dav,name);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinRestart_ApplyInitialMaterialGeometry"
PetscErrorCode pTatinRestart_ApplyInitialMaterialGeometry(pTatinCtx ctx,void *data)
{
	char name[PETSC_MAX_PATH_LEN];
	DataBucket db;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	if (ctx->materialpoint_db) {
			DataBucketDestroy(&ctx->materialpoint_db);
	}
	
	if (!StringEmpty(ctx->restart_prefix)) {
		sprintf(name,"%s/ptat3dcpf.markers_%s",ctx->restart_dir,ctx->restart_prefix);
	} else {
		sprintf(name,"%s/ptat3dcpf.markers",ctx->restart_dir);
	}
	PetscPrintf(PETSC_COMM_WORLD,"Restart: [ApplyInitialMaterialGeometry] from %s \n", name );
	DataBucketLoadFromFile(PETSC_COMM_WORLD,name,DATABUCKET_VIEW_BINARY,&db);
	ctx->materialpoint_db = db;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinRestart_ApplyInitialSolution"
PetscErrorCode pTatinRestart_ApplyInitialSolution(pTatinCtx ctx,Vec X,void *data)
{
	char name[PETSC_MAX_PATH_LEN];
	Vec Xt;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if (!StringEmpty(ctx->restart_prefix)) {
		sprintf(name,"%s/ptat3dcpf.dmda-X_%s",ctx->restart_dir,ctx->restart_prefix);
	} else {
		sprintf(name,"%s/ptat3dcpf.dmda-X",ctx->restart_dir);
	}
	PetscPrintf(PETSC_COMM_WORLD,"Restart: [ApplyInitialSolution] from %s \n", name );
	ierr = DMDALoadGlobalVectorFromFile(ctx->pack,name,&Xt);CHKERRQ(ierr);
	ierr = VecCopy(Xt,X);CHKERRQ(ierr);
	ierr = VecDestroy(&Xt);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dRestart"
PetscErrorCode pTatin3dRestart(pTatinCtx ctx)
{
	PetscBool flg;
	pTatinCtx ctx2;
	PetscViewer viewer;
	PetscErrorCode ierr;
	char start[PETSC_MAX_PATH_LEN];
	char f1[PETSC_MAX_PATH_LEN];
	char f2[PETSC_MAX_PATH_LEN];
	char f3[PETSC_MAX_PATH_LEN];
	pTatinModel model;

	
	PetscFunctionBegin;
	
	ierr = PetscOptionsGetBool(PETSC_NULL,"-restart",&ctx->restart_from_file,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(PETSC_NULL,"-restart_prefix",ctx->restart_prefix,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg==PETSC_TRUE) {
		ctx->restart_from_file = PETSC_TRUE;
	}	
	if (!ctx->restart_from_file) {
		PetscPrintf(PETSC_COMM_WORLD,"pTatin3dRestart: Required to specify a suffix for your restart file via -restart_prefix\n");
		PetscPrintf(PETSC_COMM_WORLD,"pTatin3dRestart: Unable to restart job\n");
		PetscFunctionReturn(0);
	}

	sprintf(ctx->restart_dir,"%s",ctx->outputpath);
	ierr = PetscOptionsGetString(PETSC_NULL,"-restart_directory",ctx->restart_dir,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	
	/* context */
	if (!StringEmpty(ctx->restart_prefix)) {
		sprintf(start,"%s/ptat3dcpf.ctx_%s",ctx->restart_dir,ctx->restart_prefix);
	} else {
		sprintf(start,"%s/ptat3dcpf.ctx",ctx->restart_dir);
	}
	PetscPrintf(PETSC_COMM_WORLD,"  reading %s \n", start );
	ierr = pTatin3dContextLoad(&ctx2,start);CHKERRQ(ierr);	

	/* get a copy of the model pointer - or we load again */
//	model = ctx->model;
	ierr = PetscMemcpy(ctx,ctx2,sizeof(struct _p_pTatinCtx));CHKERRQ(ierr);
//	ctx->model = model;

	/* force these again */
	ctx->restart_from_file = PETSC_TRUE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-restart_prefix",ctx->restart_prefix,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	sprintf(ctx->restart_dir,"%s",ctx->outputpath);
	ierr = PetscOptionsGetString(PETSC_NULL,"-restart_directory",ctx->restart_dir,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	
	
	ierr = pTatin3dParseOptions(ctx);CHKERRQ(ierr);
	ierr = pTatinModelLoad(ctx);CHKERRQ(ierr);
	
	/* set new function pointers for loading model */
	ctx->model->FP_pTatinModel_ApplyInitialSolution = pTatinRestart_ApplyInitialSolution;
	ctx->model->FP_pTatinModel_ApplyInitialMeshGeometry = pTatinRestart_ApplyInitialMeshGeometry;
	ctx->model->FP_pTatinModel_ApplyInitialMaterialGeometry = pTatinRestart_ApplyInitialMaterialGeometry;
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DMCoarsenHierarchy2_DA"
PetscErrorCode  DMCoarsenHierarchy2_DA(DM da,PetscInt nlevels,DM dac[])
{
  PetscErrorCode ierr;
  PetscInt       i,n,*refx,*refy,*refz;
	
  PetscFunctionBegin;
  PetscValidHeaderSpecific(da,DM_CLASSID,1);
  if (nlevels < 0) SETERRQ(((PetscObject)da)->comm,PETSC_ERR_ARG_OUTOFRANGE,"nlevels cannot be negative");
  if (nlevels == 0) PetscFunctionReturn(0);
  PetscValidPointer(dac,3);
	
  /* Get refinement factors, defaults taken from the coarse DMDA */
  ierr = PetscMalloc3(nlevels,PetscInt,&refx,nlevels,PetscInt,&refy,nlevels,PetscInt,&refz);CHKERRQ(ierr);
  for (i=0; i<nlevels; i++) {
    ierr = DMDAGetRefinementFactor(da,&refx[i],&refy[i],&refz[i]);CHKERRQ(ierr);
  }
  n = nlevels;
  ierr = PetscOptionsGetIntArray(((PetscObject)da)->prefix,"-da_refine_hierarchy_x",refx,&n,PETSC_NULL);CHKERRQ(ierr);
  n = nlevels;
  ierr = PetscOptionsGetIntArray(((PetscObject)da)->prefix,"-da_refine_hierarchy_y",refy,&n,PETSC_NULL);CHKERRQ(ierr);
  n = nlevels;
  ierr = PetscOptionsGetIntArray(((PetscObject)da)->prefix,"-da_refine_hierarchy_z",refz,&n,PETSC_NULL);CHKERRQ(ierr);
	
	
	ierr = DMDASetRefinementFactor(da,refx[nlevels-1],refy[nlevels-1],refz[nlevels-1]);CHKERRQ(ierr);
  ierr = DMCoarsen(da,((PetscObject)da)->comm,&dac[0]);CHKERRQ(ierr);
  for (i=1; i<nlevels; i++) {
    ierr = DMDASetRefinementFactor(dac[i-1],refx[nlevels-1-i],refy[nlevels-1-i],refz[nlevels-1-i]);CHKERRQ(ierr);
    ierr = DMCoarsen(dac[i-1],((PetscObject)da)->comm,&dac[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree3(refx,refy,refz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


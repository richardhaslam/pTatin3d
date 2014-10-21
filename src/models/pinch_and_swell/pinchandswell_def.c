

#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "model_utils.h"
#include "dmda_bcs.h"
#include "data_bucket.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_iterator.h"
#include "mesh_deformation.h"
#include "mesh_update.h"
#include "dmda_remesh.h"
#include "output_material_points.h"
#include "material_constants.h"
#include "energy_output.h"
#include "element_utils_q2.h"
#include "dmda_element_q2p1.h"
#include "material_point_point_location.h"
#include "material_point_std_utils.h"


typedef enum { MATRIX_IDX=0, LAYER_IDX } MaterialDomain;

typedef struct {
    PetscReal x_bar,v_bar,t_bar,eta_bar,p_bar;
    PetscReal rhs_scale;
    PetscBool output_si;
    PetscViewer logviewer;
    PetscReal domain[3];
    PetscReal exx,ezz;
} ModelCtx;



#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_PAS"
PetscErrorCode ModelInitialize_PAS(pTatinCtx ptatinctx,void *modelctx)
{
	ModelCtx           *modeldata = (ModelCtx*)modelctx;
	RheologyConstants  *rheology;
	DataBucket         materialconstants;
    PetscInt           regionidx,midx;
    PetscReal          n_exp,F,preexp_A;
    PetscBool          flg;
	PetscErrorCode     ierr;
    
    
    /* --------------------------------- scaling params --------------------------------- */
    modeldata->v_bar   = 1.0;
    modeldata->x_bar   = 1.0;  /* m */
    modeldata->eta_bar = 1.0; /* eta */
    
    PetscOptionsGetReal(PETSC_NULL,"-model_charc_v",&modeldata->v_bar,0);
    PetscOptionsGetReal(PETSC_NULL,"-model_charc_x",&modeldata->x_bar,0);
    PetscOptionsGetReal(PETSC_NULL,"-model_charc_eta",&modeldata->eta_bar,0);
    
    modeldata->t_bar = modeldata->x_bar/modeldata->v_bar; /* time (sec) */
    modeldata->p_bar = modeldata->eta_bar * modeldata->v_bar / modeldata->x_bar; /* stress */
    
    modeldata->rhs_scale = modeldata->eta_bar * modeldata->v_bar / (modeldata->x_bar * modeldata->x_bar); /* [ eta.u/(L.L) ]^-1 */
    modeldata->rhs_scale = 1.0 / modeldata->rhs_scale;

    modeldata->output_si = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULL,"-model_output_si",&modeldata->output_si,0);
    
    /* --------------------------------- domain size --------------------------------- */
    modeldata->domain[0] = 1.0;
    modeldata->domain[1] = 1.0;
    modeldata->domain[2] = 1.0;
    {
        PetscInt nd = 3;
        
        flg = PETSC_FALSE;
        PetscOptionsGetRealArray(PETSC_NULL,"-model_domain_size",modeldata->domain,&nd,&flg);
        if (flg) {
            if (nd != 3) {
                SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"param:\"model domain size\" - status:\"missing entry\" => expected 3 entries, only found %D",nd);
            }
        }
    }
    
	/* --------------------------------- rheology prescription --------------------------------- */
	ierr = pTatinGetRheology(ptatinctx,&rheology);CHKERRQ(ierr);
	rheology->rheology_type = RHEOLOGY_VISCOUS;
    rheology->nphases_active = 2;
    
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global = 1.0e+25 / modeldata->eta_bar;
	rheology->eta_lower_cutoff_global = 1.0e-25 / modeldata->eta_bar;
    
    /* Material constant */
	ierr = pTatinGetMaterialConstants(ptatinctx,&materialconstants);CHKERRQ(ierr);
    ierr = MaterialConstantsSetDefaults(materialconstants);CHKERRQ(ierr);

    /* background */
    regionidx = MATRIX_IDX;
	MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS_2,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);
    
    F        = 1.0e21;
    preexp_A = 1.0;
    n_exp    = 1.0;
    
    MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,modeldata->rhs_scale * 3150.0 * 9.81 / 10.0);
    /* eta = F . pow(preexp_A,-1/n) . pow(e,1/n-1) . exp(E/nRT) */
    MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexp_A, F, 0.0, 0.0, n_exp, 1.0);
    MaterialConstantsScaleValues_ViscosityArrh(materialconstants,regionidx, modeldata->eta_bar, modeldata->p_bar );
    
    /* layer */
    for (midx=1; midx<(PetscInt)LAYER_IDX; midx++) {
        regionidx = midx;
        
        MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS_2,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);

        //F        = 0.5;
        //preexp_A = 2.0 * 4.75e11;
        F        = 4.75e11; /* mu_0 */
        preexp_A = 1.0;
        n_exp    = 4.0;
        
        MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,modeldata->rhs_scale * 3300.0 * 9.81 / 10.0);
        /* eta = F . pow(preexp_A,-1/n) . pow(e,1/n-1) . exp(E/nR(T+T0)) */
        MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexp_A, F, 0.0, 0.0, n_exp, 1.0); /* set T0 = 1 to ensure I can calc E/nR(T-T0) */
        MaterialConstantsScaleValues_ViscosityArrh(materialconstants,regionidx, modeldata->eta_bar, modeldata->p_bar );
    }
    
    /* Material constant */
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		ierr= MaterialConstantsSetFromOptions(materialconstants,"pas",regionidx,PETSC_FALSE);CHKERRQ(ierr);
	}
    
    /* output */
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		MaterialConstantsPrintAll(materialconstants,regionidx);
	}

    /* --------------------------------- logfile --------------------------------- */
    {
        char logfile[PETSC_MAX_PATH_LEN];
        
        sprintf(logfile,"%s/model.logfile",ptatinctx->outputpath);
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,logfile,&modeldata->logviewer);CHKERRQ(ierr);

        PetscViewerASCIIPrintf(modeldata->logviewer,"# Model logfile\n");
        PetscViewerASCIIPrintf(modeldata->logviewer,"#  pTatin3d scaling used:\n");
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    x*     = %1.4e (m) \n",modeldata->x_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    v*     = %1.4e (m/s) \n",modeldata->v_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    eta*   = %1.4e (Pa s) \n",modeldata->eta_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    t*     = %1.4e (sec) \n",modeldata->t_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    p*     = %1.4e (Pa) \n",modeldata->p_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    rho.g* = %1.4e (N/m/m) \n",modeldata->rhs_scale);
        PetscViewerASCIIPrintf(modeldata->logviewer,"\n");
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_PAS"
PetscErrorCode ModelApplyInitialMeshGeometry_PAS(pTatinCtx ptatinctx,void *modelctx)
{
	ModelCtx         *modeldata = (ModelCtx*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    PetscBool        flg;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* set initial velocity field */
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

	ierr = DMDASetUniformCoordinates(dav,-modeldata->domain[0]/2.0,modeldata->domain[0]/2.0, -modeldata->domain[1]/2.0,modeldata->domain[1]/2.0, -modeldata->domain[2]/2.0,modeldata->domain[2]/2.0);CHKERRQ(ierr);

    flg = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULL,"-model_domain_2d",&flg,0);
    if (flg) {
        ierr = pTatin3d_DefineVelocityMeshGeometryQuasi2D(ptatinctx);CHKERRQ(ierr);
    }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_PAS"
PetscErrorCode ModelApplyInitialMaterialGeometry_PAS(pTatinCtx c,void *ctx)
{
	ModelCtx         *data = (ModelCtx*)ctx;
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	DataBucket       materialpoint_db;
	DataBucket       materialconstants;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
	ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
	
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
    stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
    /* set region index for all materials */
    /* mantle - background */
    for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointSet_phase_index(mpX,p, MATRIX_IDX);CHKERRQ(ierr);
    }
    
    /* layer */
    srand(0);
	for (p=0; p<n_mpoints; p++) {
        double *position_p;
        double shift,dL;
        
		ierr = MaterialPointGet_global_coord(mpX,p,&position_p);CHKERRQ(ierr);
        
        dL = 0.33;
        shift = -1.0 + 2.0 * rand()/((double)(RAND_MAX));
        shift = shift * dL * (1.0/33.0);
        if (position_p[1] > -dL*0.5) {
            if (position_p[1] < dL*0.5 + shift) {
                ierr = MaterialPointSet_phase_index(mpX,p, LAYER_IDX);CHKERRQ(ierr);
            }
        }
	}
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);

    /* set initial eta,rho */
    for (p=0; p<n_mpoints; p++) {
        int ridx;
        
		ierr = MaterialPointGet_phase_index(mpX,p, &ridx);CHKERRQ(ierr);
        
        if (ridx == MATRIX_IDX) {
            /* background */
            ierr = MaterialPointSet_viscosity(mpX,  p, (1.0e21/20.0)/data->eta_bar);CHKERRQ(ierr);
            ierr = MaterialPointSet_density(mpX,    p, -data->rhs_scale * 3150.0 * 9.81);CHKERRQ(ierr);
        } else {
            /* layer */
            ierr = MaterialPointSet_viscosity(mpX,  p, (1.0e21/1.0)/data->eta_bar);CHKERRQ(ierr);
            ierr = MaterialPointSet_density(mpX,    p, -data->rhs_scale * 3300.0 * 9.81);CHKERRQ(ierr);
        }
    }
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

/* velocity bcs */
#undef __FUNCT__
#define __FUNCT__ "PAS_VelocityBC"
PetscErrorCode PAS_VelocityBC(BCList bclist,DM dav,pTatinCtx ptatinctx,ModelCtx *modeldata)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,EAST_FACE,1.0);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,WEST_FACE,-1.0);CHKERRQ(ierr);
 
   // ierr = DirichletBC_SetConstant(bclist,dav,WEST_FACE,0,0.0);CHKERRQ(ierr);

    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,0.0);CHKERRQ(ierr);
    
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,FRONT_FACE,0.0);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,BACK_FACE,0.0);CHKERRQ(ierr);
    
    //ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,NORTH_FACE,0.0);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_PAS"
PetscErrorCode ModelApplyBoundaryCondition_PAS(pTatinCtx ptatinctx,void *modelctx)
{
	ModelCtx         *modeldata = (ModelCtx*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;
    
	PetscFunctionBegin;
	/* Define velocity boundary conditions */
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = PAS_VelocityBC(stokes->u_bclist,dav,ptatinctx,modeldata);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_PAS"
PetscErrorCode ModelApplyBoundaryConditionMG_PAS(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx ptatinctx,void *modelctx)
{
	ModelCtx         *modeldata = (ModelCtx*)modelctx;
	PetscInt         n;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	/* Define velocity boundary conditions on each level within the MG hierarchy */
	for (n=0; n<nl; n++) {
		ierr = PAS_VelocityBC(bclist[n],dav[n],ptatinctx,modeldata);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_PAS"
PetscErrorCode ModelOutput_PAS(pTatinCtx ptatinctx,Vec X,const char prefix[],void *modelctx)
{
	ModelCtx         *modeldata = (ModelCtx*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    Vec              coords,velocity,pressure;
    DataBucket       materialpoint_db;
    
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    /* get the velocity mesh */
    ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
    stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

    if (modeldata->output_si) {
        
        /* get the coordinates of the velocity mesh and scale into SI units <note, local and ghosted coordinates should be scaled> */
        ierr = DMGetCoordinates(dav,&coords);CHKERRQ(ierr);
        ierr = VecScale(coords,modeldata->x_bar);CHKERRQ(ierr);
        ierr = DMGetCoordinatesLocal(dav,&coords);CHKERRQ(ierr);
        ierr = VecScale(coords,modeldata->x_bar);CHKERRQ(ierr);
        
        /* unscale vel, p */
        ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
        ierr = VecScale(velocity,modeldata->v_bar);CHKERRQ(ierr);
        ierr = VecScale(pressure,modeldata->p_bar);CHKERRQ(ierr);
        ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    }
    
	/* ---- Velocity-Pressure Mesh Output ---- */
    /* Light weight viewer: Only v is written out. v and coords are expressed as floats */
    ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);
    ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);
    
    if (modeldata->output_si) {
        /* undo the coordinate scaling of velocity mesh <note, local and ghosted coordinates should be scaled> */
        ierr = DMGetCoordinates(dav,&coords);CHKERRQ(ierr);
        ierr = VecScale(coords,1.0/modeldata->x_bar);CHKERRQ(ierr);
        ierr = DMGetCoordinatesLocal(dav,&coords);CHKERRQ(ierr);
        ierr = VecScale(coords,1.0/modeldata->x_bar);CHKERRQ(ierr);
        
        /* unscale vel, p */
        ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
        ierr = VecScale(pressure,1.0/modeldata->p_bar);CHKERRQ(ierr);
        ierr = VecScale(velocity,1.0/modeldata->v_bar);CHKERRQ(ierr);
        ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    }
 
	/* ---- Material Point Output ---- */
	/* [1] Basic viewer: Only reports coords, regionid and other internal data */
	/*
     ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
     */
	
	/* [2] Customized viewer: User defines specific fields they want to view - NOTE not .pvd file will be created */
    {
        DataBucket                materialpoint_db;
        const int                 nf = 2;
        const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes };
        char                      mp_file_prefix[256];
        
        ierr = pTatinGetMaterialPoints(ptatinctx,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
        sprintf(mp_file_prefix,"%s_mpoints",prefix);
        ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,ptatinctx->outputpath,mp_file_prefix);CHKERRQ(ierr);
    }
    
	/* [3] Customized marker->cell viewer: Marker data is projected onto the velocity mesh. User defines specific fields */
    {
        const int                    nf = 2;
        const MaterialPointVariable  mp_prop_list[] = { MPV_viscosity, MPV_density };
        
        ierr = pTatin3d_ModelOutput_MarkerCellFields(ptatinctx,nf,mp_prop_list,prefix);CHKERRQ(ierr);
    }
    
    /* SD3D output */
    ierr = pTatinGetMaterialPoints(ptatinctx,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_PAS"
PetscErrorCode pTatinModelRegister_PAS(void)
{
	ModelCtx        *data;
	pTatinModel     m;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
    
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"pas");CHKERRQ(ierr);
    
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_PAS);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_PAS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_PAS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_PAS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_PAS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_PAS);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_PAS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_PAS);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_PAS);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


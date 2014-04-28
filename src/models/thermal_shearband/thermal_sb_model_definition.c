

#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "swarm_fields.h"
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

typedef struct {
    PetscReal L_bar,V_bar,E_bar,T_bar,t_bar,eta_bar,P_bar;
    PetscReal srate_xx;
    PetscReal srate_yy;
    PetscReal srate_zz;
    PetscReal inclusion_radius;
} ThermalSBData;


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_ThermalSB"
PetscErrorCode ModelInitialize_ThermalSB(pTatinCtx ptatctx,void *modelctx)
{
	ThermalSBData  *modeldata = (ThermalSBData*)modelctx;
	PetscErrorCode ierr;

    modeldata->L_bar = 35.0e3;  /* m */
    modeldata->T_bar = 673.0;   /* K */
    modeldata->E_bar = 5.0e-14; /* strain rate */
    
    modeldata->V_bar = modeldata->E_bar * modeldata->L_bar; /* velocity */
    modeldata->t_bar = 1.0/modeldata->E_bar;
    
    modeldata->eta_bar = 1.0; /* todo */
    modeldata->P_bar = 1.0; /* todo */
    
    /* default strain rate for bcs */
    modeldata->srate_xx = -5.0e-14;
    modeldata->srate_yy = 0.0;
    modeldata->srate_zz = 0.0;

    
    /* scale strain rate */
    modeldata->srate_xx /= modeldata->E_bar;
    modeldata->srate_yy /= modeldata->E_bar;
    modeldata->srate_zz /= modeldata->E_bar;
    
    
    /* inclusion geometry */
    modeldata->inclusion_radius = 3.0e3;
    modeldata->inclusion_radius /= modeldata->L_bar;
    
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_ThermalSB"
PetscErrorCode ModelApplyInitialMeshGeometry_ThermalSB(pTatinCtx ptatinctx,void *modelctx)
{
	ThermalSBData    *modeldata = (ThermalSBData*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    PetscReal        Lx,Ly,Lz;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    
    /* dav -> velocity mesh */
    /* dap -> pressure mesh */
    /* sets cuboid mesh */
    Lx = 35.0 * 1.0e3; /* km */
    Ly = 30.0 * 1.0e3; /* km */
    Lz = 35.0 * 1.0e3; /* km */
	ierr = DMDASetUniformCoordinates(dav,0.0,Lx/modeldata->L_bar, 0.0,Ly/modeldata->L_bar, 0.0,Lz/modeldata->L_bar);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* velocity bcs */
#undef __FUNCT__
#define __FUNCT__ "ThermalSB_VelocityBC"
PetscErrorCode ThermalSB_VelocityBC(BCList bclist,DM dav,pTatinCtx ptatinctx,ThermalSBData *modelctx)
{
    ThermalSBData  *modeldata = (ThermalSBData*)modelctx;
    PetscReal      exx,eyy,ezz;
    PetscReal      Lx,Ly,vxE,vyN,gmin[3],gmax[3];
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
    exx = modeldata->srate_xx;
	//ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,exx,0);CHKERRQ(ierr);
    
    /*
     div(u) = exx + eyy + ezz
     
     2(vx)/Lx + 2(vy)/Ly + 0 = 0
     
     exx = 2.vx/Lx

     vE = exx.Lx/2
     vW = -vE
     vy = - 0.5 * exx * Ly
     
     vy = vx Ly/Lx on the top face
    */

    
    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    Lx = gmax[0] - gmin[0];
    Ly = gmax[1] - gmin[1];
    
    /* compute vx on the east face from exx */
    vxE = modeldata->srate_xx * Lx * 0.5;

	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&vxE);CHKERRQ(ierr);
    
    /* on the upper/north face */
    vyN = - 0.5 * exx * Ly;
    
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,NORTH_FACE,vyN);CHKERRQ(ierr);
    //ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,-vyN);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,0);CHKERRQ(ierr);

    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,FRONT_FACE,0);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,BACK_FACE,0);CHKERRQ(ierr);
    
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_ThermalSB"
PetscErrorCode ModelApplyBoundaryCondition_ThermalSB(pTatinCtx ptatinctx,void *modelctx)
{
	ThermalSBData         *modeldata = (ThermalSBData*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;
    
	
	PetscFunctionBegin;
    
	/* Define velocity boundary conditions */
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = ThermalSB_VelocityBC(stokes->u_bclist,dav,ptatinctx,modeldata);CHKERRQ(ierr);
    
	/* Define boundary conditions for any other physics */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_ThermalSB"
PetscErrorCode ModelApplyBoundaryConditionMG_ThermalSB(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx ptatinctx,void *modelctx)
{
	ThermalSBData    *modeldata = (ThermalSBData*)modelctx;
	PetscInt         n;
	PetscErrorCode   ierr;
	
	
	PetscFunctionBegin;
	
	/* Define velocity boundary conditions on each level within the MG hierarchy */
	for (n=0; n<nl; n++) {
		ierr = ThermalSB_VelocityBC(bclist[n],dav[n],ptatinctx,modeldata);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_ThermalSB"
PetscErrorCode ModelOutput_ThermalSB(pTatinCtx ptatinctx,Vec X,const char prefix[],void *modelctx)
{
	ThermalSBData    *modeldata = (ThermalSBData*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    Vec              coords,velocity,pressure;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    
    /* get the velocity mesh */
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    
    /* get the coordinates of the velocity mesh and scale into SI units <note, local and ghosted coordinates should be scaled> */
    ierr = DMDAGetCoordinates(dav,&coords);CHKERRQ(ierr);
    ierr = VecScale(coords,modeldata->L_bar);CHKERRQ(ierr);
    ierr = DMDAGetGhostedCoordinates(dav,&coords);CHKERRQ(ierr);
    ierr = VecScale(coords,modeldata->L_bar);CHKERRQ(ierr);
    
	/* unscale vel, p */
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    ierr = VecScale(velocity,modeldata->V_bar);CHKERRQ(ierr);
    ierr = VecScale(pressure,modeldata->P_bar);CHKERRQ(ierr);
    
    
	/* ---- Velocity-Pressure Mesh Output ---- */
	/* [1] Standard viewer: v,p written out as binary in double */
    ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);
	
    /* [2] Light weight viewer: Only v is written out. v and coords are expressed as floats */
    ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);

	/* [3] Write out v,p into PETSc Vec. These can be used to restart pTatin */
	/*
     ierr = pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
     */
    

    /* undo the coordinate scaling of velocity mesh <note, local and ghosted coordinates should be scaled> */
    ierr = DMDAGetCoordinates(dav,&coords);CHKERRQ(ierr);
    ierr = VecScale(coords,1.0/modeldata->L_bar);CHKERRQ(ierr);
    ierr = DMDAGetGhostedCoordinates(dav,&coords);CHKERRQ(ierr);
    ierr = VecScale(coords,1.0/modeldata->L_bar);CHKERRQ(ierr);

	/* unscale vel, p */
    ierr = VecScale(pressure,1.0/modeldata->P_bar);CHKERRQ(ierr);
    ierr = VecScale(velocity,1.0/modeldata->V_bar);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

    
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
        //const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
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
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_ThermalSB"
PetscErrorCode ModelApplyInitialSolution_ThermalSB(pTatinCtx ptatinctx,Vec X,void *modelctx)
{
	ThermalSBData    *modeldata = (ThermalSBData*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    Vec              velocity,pressure;
    PetscReal        Lx,Ly,vxR,vxL,vyT,vyB,gmin[3],gmax[3];
	DMDAVecTraverse3d_InterpCtx IntpCtx;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
    
    ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
    stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    
    /* interpolate vx and vy in the interior */
    /* velocity intial condition - background strain */
    ierr = VecZeroEntries(velocity);CHKERRQ(ierr);

    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    Lx = gmax[0] - gmin[0];
    Ly = gmax[1] - gmin[1];

    vxL = -modeldata->srate_xx * Lx * 0.5;
    vxR = -vxL;
    vxL = 0.0;
    
    /* on the upper/north face */
    vyT = -0.5 * modeldata->srate_xx * Ly;
    //vyB = -vyT;
    vyB = 0.0;
    
    ierr = DMDAVecTraverse3d_InterpCtxSetUp_X(&IntpCtx,(vxR-vxL)/(Lx),vxL,0.0);CHKERRQ(ierr);
    ierr = DMDAVecTraverse3d(dav,velocity,0,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
    
    ierr = DMDAVecTraverse3d_InterpCtxSetUp_Y(&IntpCtx,(vyT-vyB)/(Ly),vyB,0.0);CHKERRQ(ierr);
    ierr = DMDAVecTraverse3d(dav,velocity,1,DMDAVecTraverse3d_Interp,(void*)&IntpCtx);CHKERRQ(ierr);
    
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_ThermalSB"
PetscErrorCode ModelApplyInitialMaterialGeometry_ThermalSB(pTatinCtx c,void *ctx)
{
	ThermalSBData    *data = (ThermalSBData*)ctx;
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	DataBucket       materialpoint_db;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    PetscReal        Ox[3],gmin[3],gmax[3],inc_rad2;
	PetscErrorCode   ierr;
    
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
    stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointSet_phase_index(mpX,p, 0);CHKERRQ(ierr);
		ierr = MaterialPointSet_viscosity(mpX,  p, 1.0);CHKERRQ(ierr);
		ierr = MaterialPointSet_density(mpX,    p, 1.0);CHKERRQ(ierr);
	}

    inc_rad2 = data->inclusion_radius * data->inclusion_radius;

	for (p=0; p<n_mpoints; p++) {
        double *position_p,r2;
        
		ierr = MaterialPointGet_global_coord(mpX,p,&position_p);CHKERRQ(ierr);
        
        //Ox[0] = (double)(gmin[0] + 0.5 * (gmax[0] - gmin[0]));
        //Ox[1] = (double)(gmin[1] + 0.5 * (gmax[1] - gmin[1]));
        //Ox[2] = (double)(gmin[2] + 0.5 * (gmax[2] - gmin[2]));
        
        Ox[0] = Ox[1] = Ox[2] = 0.0;
        
        r2 = (position_p[0]-Ox[0])*(position_p[0]-Ox[0]);
        r2 += (position_p[1]-Ox[1])*(position_p[1]-Ox[1]);
        r2 += (position_p[2]-Ox[2])*(position_p[2]-Ox[2]);
        
        if (r2 < inc_rad2) {
            ierr = MaterialPointSet_phase_index(mpX,p, 1);CHKERRQ(ierr);
            ierr = MaterialPointSet_viscosity(mpX,  p, 10.0);CHKERRQ(ierr);
            ierr = MaterialPointSet_density(mpX,    p, 5.0);CHKERRQ(ierr);
        }
	}
	
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_ThermalSB"
PetscErrorCode ModelApplyUpdateMeshGeometry_ThermalSB(pTatinCtx c,Vec X,void *ctx)
{
	ThermalSBData    *data = (ThermalSBData*)ctx;
	PetscReal        step,gmin[3],gmax[3];
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	PetscErrorCode   ierr;
    
	/* fully lagrangian update */
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);

    ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);
    
    ierr = DMDAGetBoundingBox(dav,gmin,gmax);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates1D(dav,0,gmin[0],gmax[0]);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates1D(dav,1,gmin[1],gmax[1]);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates1D(dav,2,gmin[2],gmax[2]);CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
    
    PetscFunctionReturn(0);
}
    
#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_ThermalSB"
PetscErrorCode pTatinModelRegister_ThermalSB(void)
{
	ThermalSBData   *data;
	pTatinModel     m,model;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ThermalSBData),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ThermalSBData));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
    
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"thermal_sb");CHKERRQ(ierr);
    
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_ThermalSB);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_ThermalSB);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_ViscousSinker);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

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
 **    filename:   demo.c
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


typedef enum { M_IDX=0, C_IDX } MaterialDomain;

typedef struct {
    PetscReal x_bar,v_bar,t_bar,eta_bar,p_bar;
    PetscReal rhs_scale;
    PetscBool output_si;
    PetscViewer logviewer;
    PetscReal domain[3];
    PetscReal vx_bc,vz_bc;
    PetscReal exx_bc,ezz_bc;
    PetscInt  model_geom_type;
} ModelCtx;



PetscErrorCode ModelInitialize_PD(pTatinCtx ptatinctx,void *modelctx)
{
	ModelCtx           *modeldata = (ModelCtx*)modelctx;
	RheologyConstants  *rheology;
	DataBucket         materialconstants;
    PetscInt           regionidx,midx;
    PetscBool          flg;
    PetscReal          mmyr2ms;
	PetscErrorCode     ierr;
    
    
    /* --------------------------------- scaling params --------------------------------- */
    mmyr2ms = 1.0e-3/(3600.0*24.0*365.0);
    
    modeldata->v_bar   = 2.5*mmyr2ms;
    modeldata->x_bar   = 10.0*1.0e3;  /* m */
    modeldata->eta_bar = 1.0e22; /* eta */
    
    PetscOptionsGetReal(NULL,NULL,"-model_charc_v",&modeldata->v_bar,0);
    PetscOptionsGetReal(NULL,NULL,"-model_charc_x",&modeldata->x_bar,0);
    PetscOptionsGetReal(NULL,NULL,"-model_charc_eta",&modeldata->eta_bar,0);
    
    modeldata->t_bar = modeldata->x_bar/modeldata->v_bar; /* time (sec) */
    modeldata->p_bar = modeldata->eta_bar * modeldata->v_bar / modeldata->x_bar; /* stress */
    
    modeldata->rhs_scale = modeldata->eta_bar * modeldata->v_bar / (modeldata->x_bar * modeldata->x_bar); /* [ eta.u/(L.L) ]^-1 */
    modeldata->rhs_scale = 1.0 / modeldata->rhs_scale;
    modeldata->output_si = PETSC_FALSE;
    PetscOptionsGetBool(NULL,NULL,"-model_output_si",&modeldata->output_si,0);
    
    /* --------------------------------- model geometry type ------------------------- */
    modeldata->model_geom_type = 1;
    PetscOptionsGetInt(NULL,NULL,"-model_geom_type",&modeldata->model_geom_type,0);

    /* --------------------------------- domain size --------------------------------- */
    modeldata->domain[0] = 120.0e3;
    modeldata->domain[1] = 10.0e3;
    modeldata->domain[2] = 120.0e3;
    {
        PetscInt nd = 3;
        
        flg = PETSC_FALSE;
        PetscOptionsGetRealArray(NULL,NULL,"-model_domain_size",modeldata->domain,&nd,&flg);
        if (flg) {
            if (nd != 3) {
                SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"param:\"model domain size\" - status:\"missing entry\" => expected 3 entries, only found %D",nd);
            }
        }
    }
    
    modeldata->domain[0] = modeldata->domain[0]/modeldata->x_bar;
    modeldata->domain[1] = modeldata->domain[1]/modeldata->x_bar;
    modeldata->domain[2] = modeldata->domain[2]/modeldata->x_bar;
    
    
    /* --------------------------------- Dirichlet values for v --------------------------------- */
    modeldata->vx_bc = -2.5*mmyr2ms;
    modeldata->vz_bc = 0.0;
    
    modeldata->exx_bc = 0.0;
    modeldata->ezz_bc = 0.0;
    modeldata->vx_bc /= modeldata->v_bar; printf("vx_bc = %1.4e \n",modeldata->vx_bc);
    modeldata->vz_bc /= modeldata->v_bar;
    
    modeldata->exx_bc /= 1.0;
    modeldata->ezz_bc /= 1.0;
    
	/* --------------------------------- rheology prescription --------------------------------- */
	ierr = pTatinGetRheology(ptatinctx,&rheology);CHKERRQ(ierr);
	rheology->rheology_type = RHEOLOGY_VP_STD;
    rheology->nphases_active = 2;
    
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global = 1.0e+25 / modeldata->eta_bar;
	rheology->eta_lower_cutoff_global = 1.0e-25 / modeldata->eta_bar;
    
    /* Material constant */
	ierr = pTatinGetMaterialConstants(ptatinctx,&materialconstants);CHKERRQ(ierr);
    ierr = MaterialConstantsSetDefaults(materialconstants);CHKERRQ(ierr);
    
    /* background */
    regionidx = M_IDX;
    MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_CONSTANT,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);
    
    ierr = MaterialConstantsSetValues_ViscosityConst(materialconstants,regionidx,1.0e21);CHKERRQ(ierr);
    
    MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,2.700000e+03);
    
    /* layer */
    for (midx=1; midx<=(PetscInt)C_IDX; midx++) {
        regionidx = midx;
        
        MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_CONSTANT,PLASTIC_MISES,SOFTENING_NONE,DENSITY_CONSTANT);
        ierr = MaterialConstantsSetValues_ViscosityConst(materialconstants,regionidx,1.0e24);CHKERRQ(ierr);
        
        ierr = MaterialConstantsSetValues_PlasticMises(materialconstants,regionidx,1.000000e+08,1.000000e+8);CHKERRQ(ierr);
        
        ierr = MaterialConstantsSetValues_PlasticDP(materialconstants,regionidx,30.0*M_PI/180.0,30.0*M_PI/180.0,1.000000e+08,1.000000e+08,1.000000e+07,1.000000e+10);CHKERRQ(ierr);
        
        MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,2.700000e+03);
    }
    
    /* Material constant */
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		ierr = MaterialConstantsSetFromOptions(materialconstants,"mat_",regionidx,PETSC_FALSE);CHKERRQ(ierr);
	}
    
    /* output */
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		ierr = MaterialConstantsPrintAll(materialconstants,regionidx);CHKERRQ(ierr);
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
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
        ierr = MaterialConstantsScaleAll(materialconstants,regionidx,modeldata->x_bar,modeldata->v_bar,modeldata->t_bar,modeldata->eta_bar,1.0/modeldata->rhs_scale,modeldata->p_bar);CHKERRQ(ierr);
		ierr = MaterialConstantsPrintAll(materialconstants,regionidx);CHKERRQ(ierr);
    }
    
    
    PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMeshGeometry_PD(pTatinCtx ptatinctx,void *modelctx)
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

    switch (modeldata->model_geom_type) {
        case 1:
            ierr = DMDASetUniformCoordinates(dav,0.0,modeldata->domain[0]/2.0, -modeldata->domain[1]/2.0,modeldata->domain[1]/2.0, -modeldata->domain[2]/2.0,modeldata->domain[2]/2.0);CHKERRQ(ierr);
            break;
        case 2:
            ierr = DMDASetUniformCoordinates(dav,0.0,modeldata->domain[0]/2.0, -modeldata->domain[1]/2.0,modeldata->domain[1]/2.0, -modeldata->domain[2]/2.0,modeldata->domain[2]/2.0);CHKERRQ(ierr);
            break;
        case 3:
            ierr = DMDASetUniformCoordinates(dav,-modeldata->domain[0]/2.0,modeldata->domain[0]/2.0, -modeldata->domain[1]/2.0,modeldata->domain[1]/2.0, -modeldata->domain[2]/2.0,modeldata->domain[2]/2.0);CHKERRQ(ierr);
            break;
    }

    
    flg = PETSC_FALSE;
    PetscOptionsGetBool(NULL,NULL,"-model_domain_3d",&flg,0);
    if (!flg) {
        ierr = pTatin3d_DefineVelocityMeshGeometryQuasi2D(ptatinctx);CHKERRQ(ierr);
    }
	
	PetscFunctionReturn(0);
}

PetscErrorCode PD_MaterialGeometry_ex1(pTatinCtx ptatctx,ModelCtx *mctx)
{
	MPAccess         mpX;
	int              p,n_mpoints;
	DataBucket       materialpoint_db;
    PetscErrorCode   ierr;
    
	ierr = pTatinGetMaterialPoints(ptatctx,&materialpoint_db,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
    for (p=0; p<n_mpoints; p++) {
        PetscBool in_mantle;
        double *position_p;
        double xp,yp,H,dnx,dny;
        
		ierr = MaterialPointGet_global_coord(mpX,p,&position_p);CHKERRQ(ierr);
        /* convert into km */
        xp = position_p[0] * mctx->domain[0]*0.5 * mctx->x_bar;
        yp = (position_p[1] + mctx->domain[1]*0.5) * mctx->x_bar;
        H = mctx->domain[1] * mctx->x_bar;
        
        in_mantle = PETSC_FALSE;
        
        if ((H-yp) > 22.5*1.0e3) {
            in_mantle = PETSC_TRUE;
        }
        
        
        dny = (0.25/3.0) * 30.0e3;
        dnx = 2.0 * dny;
        if (fabs(xp) < 1.0*dnx) {
            if ((H-yp) > (22.5*1.0e3-dny)) {
                in_mantle = PETSC_TRUE;
            }
        }
        
        if (!in_mantle) {
            ierr = MaterialPointSet_phase_index(mpX,p, C_IDX);CHKERRQ(ierr);
        }
    }
    
    {
        PetscBool qp_proj = PETSC_FALSE;
        
        PetscOptionsGetBool(NULL,NULL,"-qp_proj",&qp_proj,0);
        if (qp_proj) {
            for (p=0; p<n_mpoints/27; p++) {
                int pc = p*27 + 13;
                int pk,region;
                
                ierr = MaterialPointGet_phase_index(mpX,pc,&region);CHKERRQ(ierr);
                for (pk=p*27; pk<p*27+27; pk++) {
                    ierr = MaterialPointSet_phase_index(mpX,pk,region);CHKERRQ(ierr);
                }
            }
        }
    }
    
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

PetscErrorCode PD_MaterialGeometry_ex2(pTatinCtx ptatctx,ModelCtx *mctx)
{
	MPAccess         mpX;
	int              p,n_mpoints;
	DataBucket       materialpoint_db;
    PetscErrorCode   ierr;
    PetscReal        Ly = 10.0e3;
    
	ierr = pTatinGetMaterialPoints(ptatctx,&materialpoint_db,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
    for (p=0; p<n_mpoints; p++) {
        PetscBool in_mantle;
        double *position_p;
        double xp,yp,H,dnx,dny;
        
		ierr = MaterialPointGet_global_coord(mpX,p,&position_p);CHKERRQ(ierr);
        
        /* convert into km */
        xp = position_p[0] * mctx->domain[0]*0.5 * mctx->x_bar;
        yp = (position_p[1] + mctx->domain[1]*0.5) * mctx->x_bar;
        H = mctx->domain[1] * mctx->x_bar;
        
        in_mantle = PETSC_FALSE;
        
        dny = 2.0*(0.25/3.0) * Ly;
        dnx = 8.0 * dny;
        
        if (fabs(xp) < 1.0*dnx) {
            if ((H-yp) > (Ly-dny)) {
                in_mantle = PETSC_TRUE;
            }
        }
        
        if (!in_mantle) {
            ierr = MaterialPointSet_phase_index(mpX,p, C_IDX);CHKERRQ(ierr);
        }
    }
    
    {
        PetscBool qp_proj = PETSC_FALSE;
        
        PetscOptionsGetBool(NULL,NULL,"-qp_proj",&qp_proj,0);
        if (qp_proj) {
            for (p=0; p<n_mpoints/27; p++) {
                int pc = p*27 + 13;
                int pk,region;
                
                
                ierr = MaterialPointGet_phase_index(mpX,pc,&region);CHKERRQ(ierr);
                for (pk=p*27; pk<p*27+27; pk++) {
                    ierr = MaterialPointSet_phase_index(mpX,pk,region);CHKERRQ(ierr);
                }
            }
        }
    }
    
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
    
	PetscFunctionReturn(0);
}

PetscErrorCode PD_MaterialGeometry_ex3(pTatinCtx ptatctx,ModelCtx *mctx)
{
	MPAccess         mpX;
	int              p,n_mpoints;
	DataBucket       materialpoint_db;
    PetscErrorCode   ierr;
    
	ierr = pTatinGetMaterialPoints(ptatctx,&materialpoint_db,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
    for (p=0; p<n_mpoints; p++) {
        PetscBool in_mantle;
        double *position_p;
        double xp,yp,zp,H,dnx,dny;
        
		ierr = MaterialPointGet_global_coord(mpX,p,&position_p);CHKERRQ(ierr);
        /* convert into km */
        xp = position_p[0] * mctx->x_bar;
        yp = (position_p[1] + mctx->domain[1]*0.5) * mctx->x_bar;
        zp = position_p[2] * mctx->x_bar;
        H = mctx->domain[1] * mctx->x_bar;
        
        in_mantle = PETSC_FALSE;
        
        if ((H-yp) > 22.5*1.0e3) {
            in_mantle = PETSC_TRUE;
        }
        
        
        dny = (0.25/3.0) * 30.0e3;
        dnx = 2.0 * dny;
        if (fabs(xp) < 1.0*dnx) {
            if ((H-yp) > (22.5*1.0e3-dny)) {
                if (zp > 0.0) {
                    in_mantle = PETSC_TRUE;
                }
            }
        }
        
        if (!in_mantle) {
            ierr = MaterialPointSet_phase_index(mpX,p, C_IDX);CHKERRQ(ierr);
        }
    }
    
    {
        PetscBool qp_proj = PETSC_FALSE;
        
        PetscOptionsGetBool(NULL,NULL,"-qp_proj",&qp_proj,0);
        if (qp_proj) {
            for (p=0; p<n_mpoints/27; p++) {
                int pc = p*27 + 13;
                int pk,region;
                
                ierr = MaterialPointGet_phase_index(mpX,pc,&region);CHKERRQ(ierr);
                for (pk=p*27; pk<p*27+27; pk++) {
                    ierr = MaterialPointSet_phase_index(mpX,pk,region);CHKERRQ(ierr);
                }
            }
        }
    }
    
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_PD(pTatinCtx c,void *ctx)
{
	ModelCtx         *data = (ModelCtx*)ctx;
	MPAccess         mpX;
	int              p,n_mpoints;
	DataBucket       materialpoint_db;
	DataBucket       materialconstants;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;
    PetscReal        gravity[3];
	DataField        PField_visc_const,PField_dens_const;
    
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
	ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
	
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
    
    gravity[0] = 0.0;
    gravity[1] = -9.8;
    gravity[2] = 0.0;
    ierr = PhysCompStokesSetGravityVector(stokes,gravity);CHKERRQ(ierr);
    
    stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
    /* set region index for all materials */
    /* mantle - background */
    for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointSet_phase_index(mpX,p, M_IDX);CHKERRQ(ierr);
    }
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
    switch (data->model_geom_type) {
        case 1:
            /* Mises layer over viscous layer */
            ierr = PD_MaterialGeometry_ex1(c,data);CHKERRQ(ierr);
            break;
        case 2:
            /* Mises layer with viscous notch */
            ierr = PD_MaterialGeometry_ex2(c,data);CHKERRQ(ierr);
            break;
        case 3:
            /* Mises layer with viscous notch */
            ierr = PD_MaterialGeometry_ex3(c,data);CHKERRQ(ierr);
            break;
    }
    
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
	DataBucketGetDataFieldByName(materialconstants,MaterialConst_DensityConst_classname,&PField_dens_const);
    DataFieldGetAccess(PField_dens_const);
	DataBucketGetDataFieldByName(materialconstants,MaterialConst_ViscosityConst_classname,&PField_visc_const);
    DataFieldGetAccess(PField_visc_const);
    
    /* set initial eta,rho */
    for (p=0; p<n_mpoints; p++) {
        int ridx;
        MaterialConst_DensityConst   *mc_density;
        MaterialConst_ViscosityConst *mc_visc;
        
        
		ierr = MaterialPointGet_phase_index(mpX,p, &ridx);CHKERRQ(ierr);
        
        DataFieldAccessPoint(PField_dens_const,ridx,(void**)&mc_density);
        DataFieldAccessPoint(PField_visc_const,ridx,(void**)&mc_visc);
        
        ierr = MaterialPointSet_density(mpX,    p, mc_density->density);CHKERRQ(ierr);
        ierr = MaterialPointSet_viscosity(mpX,  p, mc_visc->eta0);CHKERRQ(ierr);
    }
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
	DataFieldRestoreAccess(PField_dens_const);
	DataFieldRestoreAccess(PField_visc_const);
    
	PetscFunctionReturn(0);
}

/* velocity bcs */
PetscErrorCode PD_VelocityBC(BCList bclist,DM dav,pTatinCtx ptatinctx,ModelCtx *modeldata)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,EAST_FACE, modeldata->vx_bc);CHKERRQ(ierr);
    //ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,WEST_FACE,-modeldata->vx_bc);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,WEST_FACE,0.0);CHKERRQ(ierr);
    
    // ierr = DirichletBC_SetConstant(bclist,dav,WEST_FACE,0,0.0);CHKERRQ(ierr);
    
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,0.0);CHKERRQ(ierr);
    
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,FRONT_FACE, modeldata->vz_bc);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,BACK_FACE, -modeldata->vz_bc);CHKERRQ(ierr);
    
    //ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,NORTH_FACE,0.0);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyBoundaryCondition_PD(pTatinCtx ptatinctx,void *modelctx)
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
	
	ierr = PD_VelocityBC(stokes->u_bclist,dav,ptatinctx,modeldata);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyBoundaryConditionMG_PD(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx ptatinctx,void *modelctx)
{
	ModelCtx         *modeldata = (ModelCtx*)modelctx;
	PetscInt         n;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	/* Define velocity boundary conditions on each level within the MG hierarchy */
	for (n=0; n<nl; n++) {
		ierr = PD_VelocityBC(bclist[n],dav[n],ptatinctx,modeldata);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode ModelOutput_PD(pTatinCtx ptatinctx,Vec X,const char prefix[],void *modelctx)
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
    //ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(ptatinctx,X,prefix);CHKERRQ(ierr);
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
        const int                 nf = 3;
        const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
        char                      mp_file_prefix[256];
        
        ierr = pTatinGetMaterialPoints(ptatinctx,&materialpoint_db,NULL);CHKERRQ(ierr);
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
    ierr = pTatinGetMaterialPoints(ptatinctx,&materialpoint_db,NULL);CHKERRQ(ierr);
    
    
	PetscFunctionReturn(0);
}

PetscErrorCode ModelDestroy_PD(pTatinCtx c,void *ctx)
{
	ModelCtx *data = (ModelCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* Free contents of structure */
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode pTatinModelRegister_PD(void)
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
	ierr = pTatinModelSetName(m,"pd");CHKERRQ(ierr);
    
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_PD);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_PD);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_PD);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_PD);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_PD);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_PD);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_PD);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_PD);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_PD);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


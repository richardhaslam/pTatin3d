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
 **    filename:   sd3d_definition.c
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


typedef enum { CASE_1A=0, CASE_1B, CASE_2A, CASE_2B, CASE_2C } SD3DModelType;

typedef enum { MANTLE_IDX=0, SLAB_IDX, SLAB_BASE_IDX, SLAB_EDGE_IDX, SLAB_EDGE_CENTER_IDX, SLAB_FRONT_FACE_IDX, SLAB_BACK_FACE_IDX, PLATE_IDX } SD3MaterialDomain;

typedef struct {
    PetscReal x_bar,v_bar,t_bar,eta_bar,p_bar;
    PetscReal rhs_scale;
    PetscReal slab_length;
    SD3DModelType model_type;
    PetscBool output_si;
    PetscViewer logviewer;
    PetscReal Lc,tc,vc;
} SD3DCtx;



#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_SD3D"
PetscErrorCode ModelInitialize_SD3D(pTatinCtx ptatinctx,void *modelctx)
{
	SD3DCtx            *modeldata = (SD3DCtx*)modelctx;
	RheologyConstants  *rheology;
	DataBucket         materialconstants;
    PetscInt           regionidx,midx;
    PetscReal          n_exp,F,preexp_A;
    PetscInt           mtype;
	PetscErrorCode     ierr;
    
    /* scaling params */
    //modeldata->v_bar   = 0.0000031709792; //1.0e-10; /* m/s */ /* 1.0e-11 */
    //modeldata->v_bar   = 1.0e-11; //1.0e-10; /* m/s */ /* 1.0e-11 */
    modeldata->v_bar   = 1.0e-11; // v +2.8533e-02 ; p +4.1781e-01
    modeldata->v_bar   = 1.0e-10;
    modeldata->x_bar   = 660.0 * 1.0e3;  /* m */
    modeldata->eta_bar = 10.0e24; /* eta ~ 1e22 Pa s */
    
    // preferred
    modeldata->v_bar   = 1.0e-10;
    modeldata->x_bar   = 660.0 * 1.0e3;
    modeldata->eta_bar = 1.0e24;

    // 2d
    //modeldata->v_bar   = 1.0e-9;
    //modeldata->x_bar   = 1.0e6;
    //modeldata->eta_bar = 1.0e21;
    
    modeldata->t_bar = modeldata->x_bar/modeldata->v_bar; /* time (sec) */
    modeldata->p_bar = modeldata->eta_bar * modeldata->v_bar / modeldata->x_bar; /* stress */
    
    modeldata->rhs_scale = modeldata->eta_bar * modeldata->v_bar / (modeldata->x_bar * modeldata->x_bar); /* [ eta.u/(L.L) ]^-1 */
    modeldata->rhs_scale = 1.0 / modeldata->rhs_scale;

    modeldata->output_si = PETSC_FALSE;
    PetscOptionsGetBool(NULL,NULL,"-model_output_si",&modeldata->output_si,0);
    
    modeldata->model_type = CASE_1A;
    modeldata->slab_length = 250.0;
    
    mtype = 0;
    PetscOptionsGetInt(NULL,NULL,"-model_sd3d_mtype",&mtype,NULL);
    switch (mtype) {
        case 0:
            modeldata->model_type = CASE_1A;
            PetscPrintf(PETSC_COMM_WORLD,"  [model sd3d]: CASE 1A\n");
            break;
        case 1:
            modeldata->model_type = CASE_1B;
            PetscPrintf(PETSC_COMM_WORLD,"  [model sd3d]: CASE 1B\n");
            break;
        case 2:
            modeldata->model_type = CASE_2A;
            PetscPrintf(PETSC_COMM_WORLD,"  [model sd3d]: CASE 2A\n");
            break;
        case 3:
            modeldata->model_type = CASE_2B;
            PetscPrintf(PETSC_COMM_WORLD,"  [model sd3d]: CASE 2B\n");
            break;
        case 4:
            modeldata->model_type = CASE_2C;
            PetscPrintf(PETSC_COMM_WORLD,"  [model sd3d]: CASE 2C\n");
            PetscOptionsGetReal(NULL,NULL,"-model_sd3d_case2c_slab_ly",&modeldata->slab_length,0);
            PetscPrintf(PETSC_COMM_WORLD,"  [model sd3d]: slab length = %1.4e (km)\n",modeldata->slab_length);
            break;
    }
    
	/* Rheology prescription */
	ierr = pTatinGetRheology(ptatinctx,&rheology);CHKERRQ(ierr);
	rheology->rheology_type = RHEOLOGY_VP_STD;
    rheology->nphases_active = 5;
    
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global = 1.0e+25 / modeldata->eta_bar;
	rheology->eta_lower_cutoff_global = 1.0e+18 / modeldata->eta_bar;
    
    /* Material constant */
	ierr = pTatinGetMaterialConstants(ptatinctx,&materialconstants);CHKERRQ(ierr);
    ierr = MaterialConstantsSetDefaults(materialconstants);CHKERRQ(ierr);

    /* background */
    regionidx = MANTLE_IDX;
	MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS_2,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);
    
    if ((modeldata->model_type == CASE_1A) || (modeldata->model_type == CASE_2A)) {
        F        = 1.0e21;
        preexp_A = 1.0;
        n_exp    = 1.0;
    }
    if ((modeldata->model_type == CASE_1B) || (modeldata->model_type == CASE_2B) || (modeldata->model_type == CASE_2C)) {
        //F        = 0.5;
        //preexp_A = 2.0 * 4.54e10;
        F        = 4.54e10; /* mu_0 */
        preexp_A = 1.0;
        n_exp    = 3.0;
    }
    
    MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,modeldata->rhs_scale * 3150.0);
    /* eta = F . pow(preexp_A,-1/n) . pow(e,1/n-1) . exp(E/nRT) */
    MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexp_A, F, 0.0, 0.0, n_exp, 1.0);
    MaterialConstantsScaleValues_ViscosityArrh(materialconstants,regionidx, modeldata->eta_bar, modeldata->p_bar );
    
    /* slab */
    for (midx=1; midx<(PetscInt)PLATE_IDX; midx++) {
        regionidx = midx;
        
        MaterialConstantsSetValues_MaterialType(materialconstants,regionidx,VISCOUS_ARRHENIUS_2,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);

        //F        = 0.5;
        //preexp_A = 2.0 * 4.75e11;
        F        = 4.75e11; /* mu_0 */
        preexp_A = 1.0;
        n_exp    = 4.0;
        
        MaterialConstantsSetValues_DensityConst(materialconstants,regionidx,modeldata->rhs_scale * 3300.0);
        /* eta = F . pow(preexp_A,-1/n) . pow(e,1/n-1) . exp(E/nR(T+T0)) */
        MaterialConstantsSetValues_ViscosityArrh(materialconstants,regionidx,preexp_A, F, 0.0, 0.0, n_exp, 1.0); /* set T0 = 1 to ensure I can calc E/nR(T-T0) */
        MaterialConstantsScaleValues_ViscosityArrh(materialconstants,regionidx, modeldata->eta_bar, modeldata->p_bar );
    }
    
    /* Material constant */
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		ierr= MaterialConstantsSetFromOptions(materialconstants,"sd3d",regionidx,PETSC_FALSE);CHKERRQ(ierr);
	}
    
    /* output */
	for (regionidx=0; regionidx<rheology->nphases_active; regionidx++) {
		MaterialConstantsPrintAll(materialconstants,regionidx);
	}

    /* benchmark characteristic scales */
    modeldata->Lc = 80.0 * 1.0e3; /* (m) */
    modeldata->tc = 22.56 * 1.0e6 * (60.0 * 60.0 * 24.0 * 365.0); /* (s) */
    modeldata->vc = modeldata->Lc / modeldata->tc; /* (m/s) */
    
    {
        char logfile[PETSC_MAX_PATH_LEN];
        
        sprintf(logfile,"%s/sd3d.logfile",ptatinctx->outputpath);
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,logfile,&modeldata->logviewer);CHKERRQ(ierr);

        PetscViewerASCIIPrintf(modeldata->logviewer,"# Slab detachment 3D logfile\n");
        PetscViewerASCIIPrintf(modeldata->logviewer,"#  pTatin3d scaling used:\n");
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    x*     = %1.4e (m) \n",modeldata->x_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    v*     = %1.4e (m/s) \n",modeldata->v_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    eta*   = %1.4e (Pa s) \n",modeldata->eta_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    t*     = %1.4e (sec) \n",modeldata->t_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    p*     = %1.4e (Pa) \n",modeldata->p_bar);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    rho.g* = %1.4e (N/m/m) \n",modeldata->rhs_scale);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#  Benchmark scales:\n");
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    Lc = %1.4e (m)\n",modeldata->Lc);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    tc = %1.4e (s)\n",modeldata->tc);
        PetscViewerASCIIPrintf(modeldata->logviewer,"#    vc = %1.4e (m/s)\n",modeldata->vc);
        PetscViewerASCIIPrintf(modeldata->logviewer,"# ----------------------------------------------------------------------------------------------------------------- \n");
//        PetscViewerASCIIPrintf(modeldata->logviewer,"# step | time (ptat3d)    | time (sec)       | Ys (m)            | Yn (m)            | Xm (m)            | Ymin (m)          | Ymax (m)          | Vrms (m/s)        | Phi (Pa/s)        | volume (m^3)\n");
  
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","step");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","time (ptat3d)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","time (sec)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Ys (m)");

        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Xn_edge (m)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Yn_edge (m)");
        
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Xn_back (m)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Yn_back (m)");

        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Zn_front (m)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Yn_front (m)");
        
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Ymin (m)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Ymax (m)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Vrms (m/s)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Phi (W)");
        PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","volume (m^3)");
        
        PetscViewerASCIIPrintf(modeldata->logviewer,"\n");
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_SD3D"
PetscErrorCode ModelApplyInitialMeshGeometry_SD3D(pTatinCtx ptatinctx,void *modelctx)
{
	SD3DCtx          *modeldata = (SD3DCtx*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    PetscReal        Lx,Ly,Lz;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* set initial velocity field */
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    
    Lx = 500.0 * 1.0e3; /* km */
    Ly = 660.0 * 1.0e3; /* km */
    Lz = 500.0 * 1.0e3; /* km */
	ierr = DMDASetUniformCoordinates(dav,0.0,Lx/modeldata->x_bar, 0.0,Ly/modeldata->x_bar, 0.0,Lz/modeldata->x_bar);CHKERRQ(ierr);
  {
    PetscReal gvec[] = { 0.0, -9.81, 0.0 };
    ierr = PhysCompStokesSetGravityVector(stokes,gvec);CHKERRQ(ierr);
  }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SD3D_InsertSlabEdge"
PetscErrorCode SD3D_InsertSlabEdge(DataBucket materialconstants_db,DM dav,DataBucket materialpoint_db,SD3DCtx *modeldata)
{
	MPAccess         mpX;
	MaterialConst_DensityConst      *DensityConst_data;
	DataField                       PField_DensityConst;
    PetscReal  slab_length,Ly_slab,Lz_slab,dys,dzs;
    PetscInt   nyp,nzp,j,k;
	double         tolerance;
	int            max_its;
	PetscBool      use_nonzero_guess,monitor;
	DM             cda;
	Vec            gcoords;
	PetscScalar    *LA_gcoords;
	const PetscInt *elnidx_u;
	PetscInt       nel,nen_u;
	PetscInt       lmx,lmy,lmz;
    
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
    DataBucketGetDataFieldByName(materialconstants_db,MaterialConst_DensityConst_classname,&PField_DensityConst);
	DensityConst_data      = (MaterialConst_DensityConst*)PField_DensityConst->data;
    

    slab_length = modeldata->slab_length;

    /* 
     Slab edge plane is located at
     
     x = 500.0-40.0
     660.0 - slab_length - 80.0 <= y <= 660.0 - 80.0
     0 <= z <= 250.0

    */
    
    //dye = ;
    //dze = ;
    
    Ly_slab = slab_length;
    Lz_slab = 250.0;
    nyp = 40;
    nzp = 40;
    
    dys = Ly_slab/((PetscReal)nyp);
    dzs = Lz_slab/((PetscReal)nzp);
    
	tolerance         = 1.0e-10;
	max_its           = 10;
	use_nonzero_guess = PETSC_FALSE;
	monitor           = PETSC_FALSE;
    
	/* setup for coords */
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dav,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
    
    for (k=0; k<nzp; k++) {
        for (j=0; j<nyp; j++) {
            PetscBool point_on_edge;
            double xslab_pos[3];
            int n_mpoints_orig;
            MPntStd mp_std;
            
            /* km */
            xslab_pos[0] = 500.0 - 40.0;
            xslab_pos[1] = 660.0 - slab_length - 80.0 + j*dys;
            xslab_pos[2] = k*dzs;
            
            /* non-dim units */
            xslab_pos[0] = xslab_pos[0] * 1.0e3 / modeldata->x_bar;
            xslab_pos[1] = xslab_pos[1] * 1.0e3 / modeldata->x_bar;
            xslab_pos[2] = xslab_pos[2] * 1.0e3 / modeldata->x_bar;
     
            /* check if position is INSIDE the local FE mesh */
            mp_std.coor[0] = xslab_pos[0];
            mp_std.coor[1] = xslab_pos[1];
            mp_std.coor[2] = xslab_pos[2];
            
            InverseMappingDomain_3dQ2(tolerance,max_its,
                                      use_nonzero_guess,
                                      monitor,
                                      (const PetscReal*)LA_gcoords, (const PetscInt)lmx,(const PetscInt)lmy,(const PetscInt)lmz, (const PetscInt*)elnidx_u,
                                      1, &mp_std );
            
            point_on_edge = PETSC_FALSE;
            if (mp_std.wil != -1) {
                point_on_edge = PETSC_TRUE;
            }
            
            if (point_on_edge) {
                int pidx;

                DataBucketGetSizes(materialpoint_db,&n_mpoints_orig,0,0);
                DataBucketSetSizes(materialpoint_db,n_mpoints_orig+1,-1);
                pidx = n_mpoints_orig;
                
                ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
                
                ierr = MaterialPointSet_global_coord(mpX,pidx,mp_std.coor);CHKERRQ(ierr);
                ierr = MaterialPointSet_local_coord(mpX,pidx,mp_std.xi);CHKERRQ(ierr);
                ierr = MaterialPointSet_local_element_index(mpX,pidx,mp_std.wil);CHKERRQ(ierr);
                
                ierr = MaterialPointSet_phase_index(mpX,pidx, SLAB_EDGE_IDX);CHKERRQ(ierr);

                ierr = MaterialPointSet_viscosity(mpX,  pidx, 1.0e22/modeldata->eta_bar);CHKERRQ(ierr);
                ierr = MaterialPointSet_density(mpX,    pidx, DensityConst_data[SLAB_EDGE_IDX].density);CHKERRQ(ierr);
                
                ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);

                //printf("inserting point at %1.4f %1.4f %1.4f [%d] %1.4f %1.4f %1.4f \n",mp_std.coor[0],mp_std.coor[1],mp_std.coor[2],mp_std.wil,mp_std.xi[0],mp_std.xi[1],mp_std.xi[2]);
            }
        }
    }
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_SD3D"
PetscErrorCode ModelApplyInitialMaterialGeometry_SD3D(pTatinCtx c,void *ctx)
{
	SD3DCtx          *data = (SD3DCtx*)ctx;
	MPAccess         mpX;
	int              p,n_mpoints;
	DataBucket       materialpoint_db;
	DataBucket       materialconstants;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;
    PetscReal        slab_length;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
    
	ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
	
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
    stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
    slab_length = data->slab_length;

    /* set region index for all materials */
    /* mantle - background */
    for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointSet_phase_index(mpX,p, MANTLE_IDX);CHKERRQ(ierr);
    }
    
    /* slab */
	for (p=0; p<n_mpoints; p++) {
        double *position_p,pos_si_p[3];
        
		ierr = MaterialPointGet_global_coord(mpX,p,&position_p);CHKERRQ(ierr);
        pos_si_p[0] = position_p[0] * data->x_bar / 1.0e3;
        pos_si_p[1] = position_p[1] * data->x_bar / 1.0e3;
        pos_si_p[2] = position_p[2] * data->x_bar / 1.0e3;
        
        if (pos_si_p[1] >= 660.0-80.0) {
            ierr = MaterialPointSet_phase_index(mpX,p, SLAB_IDX);CHKERRQ(ierr);
        }
        if (pos_si_p[2] <= 250.0) {
            if (pos_si_p[1] >= 660.0 - slab_length - 80.0) {
                if (pos_si_p[0] >= 500.0-40.0) {
                    ierr = MaterialPointSet_phase_index(mpX,p, SLAB_IDX);CHKERRQ(ierr);
                }
            }
        }
        
	}
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
    //ierr = SD3D_InsertSlabEdge(materialconstants,dav,materialpoint_db,data);CHKERRQ(ierr);

    { /* slab surface */
        PetscInt Nxp2[] = { 50, 50 };
        PetscReal vert_coord[4*3];
        
        vert_coord[3*0+0] = 0.0;                   vert_coord[3*0+2] = 0.0;
        vert_coord[3*1+0] = 500.0e3/data->x_bar;   vert_coord[3*1+2] = 0.0;
        vert_coord[3*2+0] = 0.0;                   vert_coord[3*2+2] = 500.0e3/data->x_bar;
        vert_coord[3*3+0] = 500.0e3/data->x_bar;   vert_coord[3*3+2] = 500.0e3/data->x_bar;
        
        vert_coord[3*0+1] = vert_coord[3*1+1] = vert_coord[3*2+1] = vert_coord[3*3+1] = (660.0-80.0)*1.0e3/data->x_bar;
        
        ierr = SwarmMPntStd_CoordAssignment_InsertWithinPlane(materialpoint_db,dav,Nxp2,PLATE_IDX,vert_coord);CHKERRQ(ierr);
    }

    { /* slab bottom */
        PetscInt Nxp2[] = { 5, 25 };
        PetscReal vert_coord[4*3];
        
        vert_coord[3*0+0] = 500.0e3/data->x_bar;             vert_coord[3*0+2] = 0.0;
        vert_coord[3*1+0] = (500.0-40.0)*1.0e3/data->x_bar;  vert_coord[3*1+2] = 0.0;
        vert_coord[3*2+0] = 500.0e3/data->x_bar;             vert_coord[3*2+2] = 250.0e3/data->x_bar;
        vert_coord[3*3+0] = (500.0-40.0)*1.0e3/data->x_bar;  vert_coord[3*3+2] = 250.0e3/data->x_bar;
        
        vert_coord[3*0+1] = vert_coord[3*1+1] = vert_coord[3*2+1] = vert_coord[3*3+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;
        
        ierr = SwarmMPntStd_CoordAssignment_InsertWithinPlane(materialpoint_db,dav,Nxp2,SLAB_BASE_IDX,vert_coord);CHKERRQ(ierr);
    }

    { /* slab front face */
        PetscInt Nxp2[] = { 11, 25 };
        PetscReal vert_coord[4*3];
        
        Nxp2[1] = (PetscInt)(slab_length/10.0);
        
        vert_coord[3*0+0] = (500.0-40.0)*1.0e3/data->x_bar;  vert_coord[3*0+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;
        vert_coord[3*1+0] = 500.0e3/data->x_bar;             vert_coord[3*1+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;
        vert_coord[3*2+0] = (500.0-40.0)*1.0e3/data->x_bar;  vert_coord[3*2+1] = (660.0-80.0)*1.0e3/data->x_bar;
        vert_coord[3*3+0] = 500.0e3/data->x_bar;             vert_coord[3*3+1] = (660.0-80.0)*1.0e3/data->x_bar;
        
        vert_coord[3*0+2] = vert_coord[3*1+2] = vert_coord[3*2+2] = vert_coord[3*3+2] = 250.0e3/data->x_bar;
        
        ierr = SwarmMPntStd_CoordAssignment_InsertWithinPlane(materialpoint_db,dav,Nxp2,SLAB_FRONT_FACE_IDX,vert_coord);CHKERRQ(ierr);
    }

    { /* slab back face */
        PetscInt Nxp2[] = { 2, 25 };
        PetscReal vert_coord[4*3];
        
        Nxp2[1] = (PetscInt)(slab_length/10.0);
        
        vert_coord[3*0+0] = (500.0-40.0)*1.0e3/data->x_bar;  vert_coord[3*0+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;
        vert_coord[3*1+0] = (500.0-39.999)*1.0e3/data->x_bar; vert_coord[3*1+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;
        vert_coord[3*2+0] = (500.0-40.0)*1.0e3/data->x_bar;  vert_coord[3*2+1] = (660.0-80.0)*1.0e3/data->x_bar;
        vert_coord[3*3+0] = (500.0-39.999)*1.0e3/data->x_bar; vert_coord[3*3+1] = (660.0-80.0)*1.0e3/data->x_bar;
        
        vert_coord[3*0+2] = vert_coord[3*1+2] = vert_coord[3*2+2] = vert_coord[3*3+2] = 0.0e3/data->x_bar;
        
        ierr = SwarmMPntStd_CoordAssignment_InsertWithinPlane(materialpoint_db,dav,Nxp2,SLAB_BACK_FACE_IDX,vert_coord);CHKERRQ(ierr);
    }

    { /* slab face */
        PetscInt Nxp2[] = { 2*25, 2*25 };
        PetscReal vert_coord[4*3];

        Nxp2[0] = (PetscInt)(slab_length/10.0);

        vert_coord[3*0+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;     vert_coord[3*0+2] = 0.0;
        vert_coord[3*1+1] = (660.0-80.0)*1.0e3/data->x_bar;                 vert_coord[3*1+2] = 0.0;
        vert_coord[3*2+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;     vert_coord[3*2+2] = 250.0e3/data->x_bar;
        vert_coord[3*3+1] = (660.0-80.0)*1.0e3/data->x_bar;                 vert_coord[3*3+2] = 250.0e3/data->x_bar;
        
        vert_coord[3*0+0] = vert_coord[3*1+0] = vert_coord[3*2+0] = vert_coord[3*3+0] = (500.0-40.0)*1.0e3/data->x_bar;
        
        ierr = SwarmMPntStd_CoordAssignment_InsertWithinPlane(materialpoint_db,dav,Nxp2,SLAB_EDGE_IDX,vert_coord);CHKERRQ(ierr);
    }
    
    { /* slab center plane */
        PetscInt Nxp2[] = { 2*25, 2*25 };
        PetscReal vert_coord[4*3];
        
        Nxp2[0] = (PetscInt)(slab_length/10.0);
        
        vert_coord[3*0+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;     vert_coord[3*0+2] = 0.0;
        vert_coord[3*1+1] = (660.0-80.0)*1.0e3/data->x_bar;                 vert_coord[3*1+2] = 0.0;
        vert_coord[3*2+1] = (660.0-80.0-slab_length)*1.0e3/data->x_bar;     vert_coord[3*2+2] = 250.0e3/data->x_bar;
        vert_coord[3*3+1] = (660.0-80.0)*1.0e3/data->x_bar;                 vert_coord[3*3+2] = 250.0e3/data->x_bar;
        
        vert_coord[3*0+0] = vert_coord[3*1+0] = vert_coord[3*2+0] = vert_coord[3*3+0] = 500.0*1.0e3/data->x_bar;
        
        ierr = SwarmMPntStd_CoordAssignment_InsertWithinPlane(materialpoint_db,dav,Nxp2,SLAB_EDGE_CENTER_IDX,vert_coord);CHKERRQ(ierr);
    }
    
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    /* set initial eta,rho */
    for (p=0; p<n_mpoints; p++) {
        int ridx;
        
		ierr = MaterialPointGet_phase_index(mpX,p, &ridx);CHKERRQ(ierr);
        
        if (ridx == MANTLE_IDX) {
            /* mantle - background */
            ierr = MaterialPointSet_viscosity(mpX,  p, 1.0e21/data->eta_bar);CHKERRQ(ierr);
            ierr = MaterialPointSet_density(mpX,    p, data->rhs_scale * 3150.0);CHKERRQ(ierr);
        } else {
            /* slab */
            ierr = MaterialPointSet_viscosity(mpX,  p, 1.0e22/data->eta_bar);CHKERRQ(ierr);
            ierr = MaterialPointSet_density(mpX,    p, data->rhs_scale * 3300.0);CHKERRQ(ierr);
        }
    }
    ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

/* velocity bcs */
#undef __FUNCT__
#define __FUNCT__ "SD3D_VelocityBC"
PetscErrorCode SD3D_VelocityBC(BCList bclist,DM dav,pTatinCtx ptatinctx,SD3DCtx *modeldata)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,EAST_FACE,0.0);CHKERRQ(ierr);
 
    ierr = DirichletBC_SetConstant(bclist,dav,WEST_FACE,0,0.0);CHKERRQ(ierr);
    ierr = DirichletBC_SetConstant(bclist,dav,WEST_FACE,1,0.0);CHKERRQ(ierr);
    ierr = DirichletBC_SetConstant(bclist,dav,WEST_FACE,2,0.0);CHKERRQ(ierr);

    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,SOUTH_FACE,0.0);CHKERRQ(ierr);
    
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,FRONT_FACE,0.0);CHKERRQ(ierr);
    ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,BACK_FACE,0.0);CHKERRQ(ierr);
    
    switch (modeldata->model_type) {
        case CASE_1A:
            ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,NORTH_FACE,0.0);CHKERRQ(ierr);
            break;
        case CASE_1B:
            ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,NORTH_FACE,0.0);CHKERRQ(ierr);
            break;
    }
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_SD3D"
PetscErrorCode ModelApplyBoundaryCondition_SD3D(pTatinCtx ptatinctx,void *modelctx)
{
	SD3DCtx          *modeldata = (SD3DCtx*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;
    
	PetscFunctionBegin;
	/* Define velocity boundary conditions */
	ierr = pTatinGetStokesContext(ptatinctx,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = SD3D_VelocityBC(stokes->u_bclist,dav,ptatinctx,modeldata);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_SD3D"
PetscErrorCode ModelApplyBoundaryConditionMG_SD3D(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx ptatinctx,void *modelctx)
{
	SD3DCtx          *modeldata = (SD3DCtx*)modelctx;
	PetscInt         n;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	/* Define velocity boundary conditions on each level within the MG hierarchy */
	for (n=0; n<nl; n++) {
		ierr = SD3D_VelocityBC(bclist[n],dav[n],ptatinctx,modeldata);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SD3DOutput_ComputeDs"
PetscErrorCode SD3DOutput_ComputeDs(DataBucket db,PetscReal *Ds)
{
    PetscReal gmin[3],gmax[3];
    PetscErrorCode ierr;
    
    ierr = MPntStdComputeBoundingBoxInRangeInRegion(db,NULL,NULL,SLAB_BASE_IDX,gmin,gmax);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"slab base range: gmin %1.4e %1.4e %1.4e\n",gmin[0],gmin[1],gmin[2]);
    PetscPrintf(PETSC_COMM_WORLD,"slab base range: gmax %1.4e %1.4e %1.4e\n",gmax[0],gmax[1],gmax[2]);
    *Ds = gmin[1]; /* min y coordinate of the markers on the slab base */

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SD3DOutput_ComputeWDn_minX"
PetscErrorCode SD3DOutput_ComputeWDn_minX(DataBucket db,PetscReal *W,PetscReal *Dn)
{
    PetscBool      coord_mask[] = { PETSC_TRUE, PETSC_FALSE, PETSC_FALSE };
    PetscReal      w,coord[3];
    PetscReal      gmin[3],gmax[3];
    int            pidx_w;
    PetscMPIInt    rank,rank_w;
    double         shared_coord[3];
	MPAccess       mpX;
    PetscErrorCode ierr;
    
    ierr = MPntStdComputeBoundingBoxInRangeInRegion(db,NULL,NULL,SLAB_EDGE_IDX,gmin,gmax);CHKERRQ(ierr);
    w = gmax[0];
    PetscPrintf(PETSC_COMM_WORLD,"slab edge range: gmax(x) %1.4e\n",w);

    coord[0] = w;
    coord[1] = 0.0;
    coord[2] = 0.0;
    ierr = MPntStdIdentifyFromPosition(db,coord,coord_mask,SLAB_EDGE_IDX,1.0e-16,&pidx_w,&rank_w);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"slab edge identified: pidx_w,rank_w %d %d\n",pidx_w,rank_w);

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
    shared_coord[0] = shared_coord[1] = shared_coord[2] = 1.0e32;
    if (rank == rank_w) {
        double  *pos_p;

		ierr = MaterialPointGet_global_coord(mpX,pidx_w,&pos_p);CHKERRQ(ierr);
        shared_coord[0] = pos_p[0];
        shared_coord[1] = pos_p[1];
        shared_coord[2] = pos_p[2];
    }
	ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);

    ierr =  MPI_Bcast(shared_coord,3,MPI_DOUBLE,rank_w,PETSC_COMM_WORLD);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"slab edge range: W,Dn %1.4e %1.4e\n",w,shared_coord[1]);
    
    *W = w;
    *Dn = shared_coord[1];
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SD3DOutput_ComputeWDn_backface_minX"
PetscErrorCode SD3DOutput_ComputeWDn_backface_minX(DataBucket db,PetscReal *W,PetscReal *Dn)
{
    PetscBool      coord_mask[] = { PETSC_TRUE, PETSC_FALSE, PETSC_FALSE };
    PetscReal      w,coord[3];
    PetscReal      gmin[3],gmax[3];
    int            pidx_w;
    PetscMPIInt    rank,rank_w;
    double         shared_coord[3];
	MPAccess       mpX;
    PetscErrorCode ierr;
    
    ierr = MPntStdComputeBoundingBoxInRangeInRegion(db,NULL,NULL,SLAB_BACK_FACE_IDX,gmin,gmax);CHKERRQ(ierr);
    w = gmax[0];
    PetscPrintf(PETSC_COMM_WORLD,"slab edge (back face) range: gmax(x) %1.4e\n",w);
    
    coord[0] = w;
    coord[1] = 0.0;
    coord[2] = 0.0;
    ierr = MPntStdIdentifyFromPosition(db,coord,coord_mask,SLAB_BACK_FACE_IDX,1.0e-16,&pidx_w,&rank_w);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"slab edge (back face) identified: pidx_w,rank_w %d %d\n",pidx_w,rank_w);
    
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
    shared_coord[0] = shared_coord[1] = shared_coord[2] = 1.0e32;
    if (rank == rank_w) {
        double  *pos_p;
        
		ierr = MaterialPointGet_global_coord(mpX,pidx_w,&pos_p);CHKERRQ(ierr);
        shared_coord[0] = pos_p[0];
        shared_coord[1] = pos_p[1];
        shared_coord[2] = pos_p[2];
    }
	ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);
    
    ierr =  MPI_Bcast(shared_coord,3,MPI_DOUBLE,rank_w,PETSC_COMM_WORLD);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"slab edge (back face) range: W,Dn %1.4e %1.4e\n",w,shared_coord[1]);
    
    *W = w;
    *Dn = shared_coord[1];
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SD3DOutput_ComputeWDn_frontface_minZ"
PetscErrorCode SD3DOutput_ComputeWDn_frontface_minZ(DataBucket db,PetscReal *W,PetscReal *Dn)
{
    PetscBool      coord_mask[] = { PETSC_FALSE, PETSC_FALSE, PETSC_TRUE };
    PetscReal      w,coord[3];
    PetscReal      gmin[3],gmax[3];
    int            pidx_w;
    PetscMPIInt    rank,rank_w;
    double         shared_coord[3];
	MPAccess       mpX;
    PetscErrorCode ierr;
    
    ierr = MPntStdComputeBoundingBoxInRangeInRegion(db,NULL,NULL,SLAB_FRONT_FACE_IDX,gmin,gmax);CHKERRQ(ierr);
    w = gmin[2];
    PetscPrintf(PETSC_COMM_WORLD,"slab edge (front face) range: gmin(z) %1.4e\n",w);
    
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = w;
    ierr = MPntStdIdentifyFromPosition(db,coord,coord_mask,SLAB_FRONT_FACE_IDX,1.0e-16,&pidx_w,&rank_w);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"slab edge (front face) identified: pidx_w,rank_w %d %d\n",pidx_w,rank_w);
    
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
    shared_coord[0] = shared_coord[1] = shared_coord[2] = 1.0e32;
    if (rank == rank_w) {
        double  *pos_p;
        
		ierr = MaterialPointGet_global_coord(mpX,pidx_w,&pos_p);CHKERRQ(ierr);
        shared_coord[0] = pos_p[0];
        shared_coord[1] = pos_p[1];
        shared_coord[2] = pos_p[2];
    }
	ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);
    
    ierr =  MPI_Bcast(shared_coord,3,MPI_DOUBLE,rank_w,PETSC_COMM_WORLD);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"slab edge (front face) range: W,Dn %1.4e %1.4e\n",w,shared_coord[1]);
    
    *W = w;
    *Dn = shared_coord[1];
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SD3DOutput_ComputeZrange"
PetscErrorCode SD3DOutput_ComputeZrange(DM dav,PetscReal range[])
{
    PetscReal gmin[3],gmax[3];
    PetscErrorCode ierr;
    
    ierr =  DMDAComputeBoundingBoxBoundaryFace(dav,NORTH_FACE,gmin,gmax);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"topo range: %1.4e %1.4e\n",gmin[1],gmax[1]);
    range[0] = gmin[1];
    range[1] = gmax[1];
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SD3DOutput_ComputeVrms"
PetscErrorCode SD3DOutput_ComputeVrms(DM dav,Vec v,PetscReal *volume,PetscReal *vrms)
{
    PetscReal v2,vol,VRMS;
    PetscErrorCode ierr;
    
    ierr =  StokesComputeVRMS(dav,v,&v2,&vol);CHKERRQ(ierr);
	VRMS = PetscSqrtReal(v2/vol);
    PetscPrintf(PETSC_COMM_WORLD,"vol,int v.v,vrms: %1.4e %1.4e %1.4e\n",vol,v2,VRMS);
    *volume = vol;
    *vrms = VRMS;
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SD3DOutput_ComputeViscousDissipation"
PetscErrorCode SD3DOutput_ComputeViscousDissipation(pTatinCtx ptatinctx,PhysCompStokes stokes,Vec X,PetscReal *_value)
{
    DM              stokes_pack,dav,dap;
    Vec             velocity,pressure;
	Vec             Uloc,Ploc;
	PetscScalar     *LA_Uloc,*LA_Ploc;
    Quadrature      volQ;
    PetscReal       value[3];
    PetscErrorCode  ierr;
    
	stokes_pack = stokes->stokes_pack;
    volQ        = stokes->volQ;
    
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
    ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    
	ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);

	ierr = pTatin_EvaluateRheologyNonlinearities(ptatinctx,dav,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
	
    ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);

    value[0] = value[1] = value[2] = 0.0;
    //ierr =  StokesComputeViscousDissipation(dav,dap,velocity,pressure,volQ,0,&value[0]);CHKERRQ(ierr);
    ierr =  StokesComputeViscousDissipation(dav,dap,velocity,pressure,volQ,1,&value[1]);CHKERRQ(ierr);
    //ierr =  StokesComputeViscousDissipation(dav,dap,velocity,pressure,volQ,2,&value[2]);CHKERRQ(ierr);

    //printf("sigma:e %+1.10e \n",value[0]);
    //printf("tau:e   %+1.10e \n",value[1]);
    //printf("-p:e    %+1.10e \n",value[2]);
    
    *_value = value[1]; /* tau_ij : e_ij */
    
    ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    
    PetscPrintf(PETSC_COMM_WORLD,"Phi: %1.4e\n",*value);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_SD3D"
PetscErrorCode ModelOutput_SD3D(pTatinCtx ptatinctx,Vec X,const char prefix[],void *modelctx)
{
	SD3DCtx          *modeldata = (SD3DCtx*)modelctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
    Vec              coords,velocity,pressure;
    DataBucket       materialpoint_db;
    PetscReal        Ds,Dn,W,Dnx,Wx,Dnz,Wz,range[2],Vrms,Phi,e_bar,volume,x_bar3;
    
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

    Ds = Dn = W = 0.0;
    Dnx = Wx = Dnz = Wz = 0.0;
    range[0] = range[1] = 0.0;
    Vrms = Phi = 0.0;
    
    ierr = SD3DOutput_ComputeDs(materialpoint_db,&Ds);CHKERRQ(ierr);

    ierr = SD3DOutput_ComputeWDn_minX(materialpoint_db,&W,&Dn);CHKERRQ(ierr);
    //Dn = 660.0e3/modeldata->x_bar - Dn;
    //W = 500.0e3/modeldata->x_bar - W;

    ierr = SD3DOutput_ComputeWDn_backface_minX(materialpoint_db,&Wx,&Dnx);CHKERRQ(ierr);
    ierr = SD3DOutput_ComputeWDn_frontface_minZ(materialpoint_db,&Wz,&Dnz);CHKERRQ(ierr);
    
    
    ierr = SD3DOutput_ComputeZrange(dav,range);CHKERRQ(ierr);
    
    ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    ierr = SD3DOutput_ComputeVrms(dav,velocity,&volume,&Vrms);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    
    ierr = SD3DOutput_ComputeViscousDissipation(ptatinctx,stokes,X,&Phi);CHKERRQ(ierr);
    e_bar = modeldata->v_bar / modeldata->x_bar; /* 1/s */
    x_bar3 = modeldata->x_bar * modeldata->x_bar * modeldata->x_bar;
    
    /* 
     Notes:
     - dissipation must be scaled the following scales: characteristic pressure, characteristic strain rate, characteristic volume
     - dissipation isn't normalized by the volume, THUS, to define dissipation in the full domain we need to multiply the result here (half domain) by 2
     
    */
    /*
    PetscViewerASCIIPrintf(modeldata->logviewer,"%.6D %1.12e %1.12e %+1.12e %+1.12e %+1.12e %+1.12e %+1.12e %+1.12e %+1.12e %+1.12e\n",
                           ptatinctx->step, ptatinctx->time, ptatinctx->time*modeldata->t_bar,
                           Ds*modeldata->x_bar,       Dn*modeldata->x_bar,               W*modeldata->x_bar,
                           range[0]*modeldata->x_bar, range[1]*modeldata->x_bar,
                           Vrms*modeldata->v_bar,     2.0 * Phi*modeldata->p_bar*e_bar*x_bar3, volume*x_bar3);
    */

    /*
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","step");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","time (ptat3d)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","time (sec)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Ys (m)");
     
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Xn_edge (m)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Yn_edge (m)");
     
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Xn_back (m)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Yn_back (m)");
     
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Zn_front (m)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Yn_front (m)");
     
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Ymin (m)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Ymax (m)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Vrms (m/s)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","Phi (W)");
     PetscViewerASCIIPrintf(modeldata->logviewer,"| %16.20s ","volume (m^3)");
    */
    
    PetscViewerASCIIPrintf(modeldata->logviewer,"%19.6D ",ptatinctx->step);
    PetscViewerASCIIPrintf(modeldata->logviewer,"%8.12e ",ptatinctx->time);
    PetscViewerASCIIPrintf(modeldata->logviewer,"%8.12e ",ptatinctx->time*modeldata->t_bar);

    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",Ds*modeldata->x_bar);

    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",W*modeldata->x_bar);
    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",Dn*modeldata->x_bar);

    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",Wx*modeldata->x_bar);
    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",Dnx*modeldata->x_bar);

    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",Wz*modeldata->x_bar);
    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",Dnz*modeldata->x_bar);

    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",range[0]*modeldata->x_bar);
    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",range[1]*modeldata->x_bar);

    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",Vrms*modeldata->v_bar);
    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",2.0 * Phi*modeldata->p_bar*e_bar*x_bar3);
    PetscViewerASCIIPrintf(modeldata->logviewer,"%+1.11e ",volume*x_bar3);
    
    PetscViewerASCIIPrintf(modeldata->logviewer,"\n");
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_SD3D"
PetscErrorCode ModelDestroy_SD3D(pTatinCtx c,void *ctx)
{
	SD3DCtx *data = (SD3DCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* Free contents of structure */
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_SD3D"
PetscErrorCode pTatinModelRegister_SD3D(void)
{
	SD3DCtx         *data;
	pTatinModel     m;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(SD3DCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(SD3DCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
    
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"sd3d");CHKERRQ(ierr);
    
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_SD3D);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_SD3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_SD3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_SD3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_SD3D);CHKERRQ(ierr);
	//ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_ThermalSB);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_SD3D);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_SD3D);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


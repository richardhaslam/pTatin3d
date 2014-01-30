
/*
 
 Developed by 
 Xiaolin Mao [xlmao@gps.caltech.edu]
 
 Defines the iso-viscous models described here
 
 @article {GJI:GJI23,
   author = {Blankenbach, B. and Busse, F. and Christensen, U. and Cserepes, L. and Gunkel, D. and Hansen, U. and Harder, H. and Jarvis, G. and Koch, M. and Marquart, G. and Moore, D. and Olson, P. and Schmeling, H. and Schnaubelt, T.},
   title = {A benchmark comparison for mantle convection codes},
   journal = {Geophysical Journal International},
   volume = {98},
   number = {1},
   pages = {23--38},
   year = {1989},
 }
 
*/


#define _GNU_SOURCE
#include "petsc.h"
#include "petscmath.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_utils.h"
#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_iterator.h"
#include "mesh_update.h"
#include "output_material_points.h"
#include "material_point_std_utils.h"
#include "energy_output.h"
#include "geometry_object.h"
#include "ptatin3d_energy.h"
#include "model_utils.h"

#include "model_convection2d_ctx.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Thermal_Convection2d"
PetscErrorCode ModelInitialize_Thermal_Convection2d(pTatinCtx c,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	PetscReal         cm_yr2m_s;
	RheologyConstants *rheology;
	DataBucket        materialconstants;
	PetscInt          regionidx;
	PetscReal         DeltaT;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatin3d_DefineVelocityMeshQuasi2D(c);CHKERRQ(ierr);
	
	cm_yr2m_s = 1.0e-2 / ( 365.25 * 24.0 * 60.0 * 60.0 ) ;
	
	/*box geometry*/
	data->Lx = 10.0e5;
	data->Ly = 0.0e5;
	data->Lz = 0.01e5;
	data->Ox = 0.0;
	data->Oy = -10.0e5;
	data->Oz = 0.0e5;
	
	/*material constants*/	
	data->diffusivity[0] = 1.0e-6;
	data->alpha[0]       = 2.5e-2; /* 2.5e-5 - table 1 */
	data->H[0]           = 0.0; /* heat generation */
	
	/* density */
	data->rho[0] = 4000.0;
		
	/* temperature */
	data->Ttop = 0.0;
	data->Tbot = 1.0; //	
	
	/* viscosity */
	data->eta[0] = 1.0e21; /* 2.5e19 Pa.s */
	
	/* scaling bar */
	data->length_bar    = 1000.0 * 1.0e3;
	data->viscosity_bar = 1.0e21;
	data->velocity_bar  = 1.0e-6;
	DeltaT              = 1.0e3;
	
	ierr = pTatinGetRheology(c,&rheology);CHKERRQ(ierr);
	rheology->rheology_type = RHEOLOGY_VP_STD;
	rheology->nphases_active = 1;
	rheology->apply_viscosity_cutoff_global = PETSC_TRUE;
	rheology->eta_upper_cutoff_global = 1.0e+26;
	rheology->eta_lower_cutoff_global = 1.0e+20;
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_Lx",&data->Lx,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_Ox",&data->Ox,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_Ly",&data->Ly,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_Oy",&data->Oy,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_Lz",&data->Lz,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_Oz",&data->Oz,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_velocity",&data->velocity,PETSC_NULL);CHKERRQ(ierr); 	
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_length_bar",&data->length_bar,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_viscosity_bar",&data->viscosity_bar,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_velocity_bar",&data->velocity_bar,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_Ttop",&data->Ttop,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_Tbot",&data->Tbot,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_eta0",&data->eta[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_rho0",&data->rho[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_diffusivity0",&data->diffusivity[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_alpha0",&data->alpha[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_conv2d_H0",&data->H[0],PETSC_NULL);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"[thermal_convection2d]: Ra = %1.4e \n",data->alpha[0]*10*DeltaT*data->Lx*data->Lx*data->Lx/(data->diffusivity[0]*data->eta[0]));
	
	data->velocity = data->velocity * cm_yr2m_s;
	
	ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
	MaterialConstantsSetDefaults(materialconstants);
	
	MaterialConstantsSetValues_MaterialType(materialconstants,0,VISCOUS_CONSTANT,PLASTIC_NONE,SOFTENING_NONE,DENSITY_BOUSSINESQ);
	MaterialConstantsSetValues_ViscosityConst(materialconstants,0,data->eta[0]);
	MaterialConstantsSetValues_DensityBoussinesq(materialconstants,0,data->rho[0],data->alpha[0],0.0);
	MaterialConstantsSetValues_DensityConst(materialconstants,0,data->rho[0]);
	
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) { 
		MaterialConstantsPrintAll(materialconstants,regionidx);
	} 
	
	/* Compute additional scaling parameters */
	data->time_bar      = data->length_bar / data->velocity_bar;
	data->pressure_bar  = data->viscosity_bar/data->time_bar;
	data->density_bar   = data->pressure_bar / data->length_bar;
	
	PetscPrintf(PETSC_COMM_WORLD,"[thermal_convection2d]:  during the solve scaling will be done using \n");
	PetscPrintf(PETSC_COMM_WORLD,"  L*    : %1.4e [m]\n", data->length_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  U*    : %1.4e [m.s^-1]\n", data->velocity_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  t*    : %1.4e [s]\n", data->time_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  eta*  : %1.4e [Pa.s]\n", data->viscosity_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  rho*  : %1.4e [kg.m^-3]\n", data->density_bar );
	PetscPrintf(PETSC_COMM_WORLD,"  P*    : %1.4e [Pa]\n", data->pressure_bar );
	
	PetscPrintf(PETSC_COMM_WORLD,"[thermal_convection2d]:  parameters \n");
	PetscPrintf(PETSC_COMM_WORLD,"  alpha    : %1.4e [K^-1]\n", data->alpha[0] );
	PetscPrintf(PETSC_COMM_WORLD,"  eta_upper_cutoff   : %1.4e [Pas]\n", rheology->eta_upper_cutoff_global );
	PetscPrintf(PETSC_COMM_WORLD,"  eta_lower_cutoff   : %1.4e [Pas]\n", rheology->eta_lower_cutoff_global );
	
	/* scale viscosity cutoff */
	rheology->eta_lower_cutoff_global = rheology->eta_lower_cutoff_global / data->viscosity_bar;
	rheology->eta_upper_cutoff_global = rheology->eta_upper_cutoff_global / data->viscosity_bar;
	/* scale length */
	data->Lx = data->Lx / data->length_bar;
	data->Ly = data->Ly / data->length_bar;
	data->Lz = data->Lz / data->length_bar;
	data->Ox = data->Ox / data->length_bar;
	data->Oy = data->Oy / data->length_bar;
	data->Oz = data->Oz / data->length_bar;
	
	/* scale velocity */
	data->velocity = data->velocity/data->velocity_bar;
	
	/* scale rho0 */
	data->rho[0] = data->rho[0]/data->density_bar;
	
	/* scale diffusivity and alpha */
	data->diffusivity[0] = data->diffusivity[0]/data->length_bar/data->length_bar*data->time_bar;
	/* Since you scaled temperature by \Delta T = 1000, alphs should be increased by a factor 1000 */
	//data->alpha[0]       = data->alpha[0] * DeltaT;
#if 1
	// scale material properties
	for (regionidx=0; regionidx<rheology->nphases_active;regionidx++) {
		MaterialConstantsScaleAll(materialconstants,regionidx,data->length_bar,data->velocity_bar,data->time_bar,data->viscosity_bar,data->density_bar,data->pressure_bar);
	}

	ierr = PetscOptionsInsertString("-activate_energy");CHKERRQ(ierr);
#endif	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Thermal_Convection2d_VelocityBC"
PetscErrorCode Thermal_Convection2d_VelocityBC(BCList bclist,DM dav,pTatinCtx c,ModelThermal_Convection2dCtx *data)
{
	PetscScalar    zero = 0.0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,EAST_FACE);CHKERRQ(ierr);
	
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,WEST_FACE);CHKERRQ(ierr);
	
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,NORTH_FACE);CHKERRQ(ierr);
	
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,SOUTH_FACE);CHKERRQ(ierr);
	
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	//ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,BACK_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,FRONT_FACE);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Thermal_Convection2d"
PetscErrorCode ModelApplyBoundaryCondition_Thermal_Convection2d(pTatinCtx c,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	PetscScalar    zero = 0.0,velocity;
	PetscBool      active_energy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = Thermal_Convection2d_VelocityBC(c->stokes_ctx->u_bclist,c->stokes_ctx->dav,c,data);CHKERRQ(ierr);
	
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		
		PetscReal      val_T;
		PhysCompEnergy energy;
		BCList         bclist;
		DM             daT;
		
		ierr   = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		daT    = energy->daT;
		bclist = energy->T_bclist;
		
		val_T = data->Tbot;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);
		
		val_T = data->Ttop;
		ierr = DMDABCListTraverse3d(bclist,daT,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&val_T);CHKERRQ(ierr);		
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Thermal_Convection2d"
PetscErrorCode ModelApplyBoundaryConditionMG_Thermal_Convection2d(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	PetscInt       n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = Thermal_Convection2d_VelocityBC(bclist[n],dav[n],c,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Thermal_Convection2d"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Thermal_Convection2d(pTatinCtx c,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Thermal_Convection2d"
PetscErrorCode ModelApplyInitialMeshGeometry_Thermal_Convection2d(pTatinCtx c,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, data->Ox,data->Lx, data->Oy,data->Ly, data->Oz, data->Lz);CHKERRQ(ierr);
	ierr = pTatin3d_DefineVelocityMeshGeometryQuasi2D(c);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Thermal_Convection2d"
PetscErrorCode ModelApplyInitialMaterialGeometry_Thermal_Convection2d(pTatinCtx c,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_pls,PField_stokes, PField_energy;
	int                    phase;
	PetscBool              active_energy;
	MPAccess               mpX;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetMaterialPoints(c,&db,PETSC_NULL);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd       *material_point;
		MPntPStokes   *mpprop_stokes;
		double        *position,ycoord,xcoord,zcoord;
		
		DataFieldAccessPoint(PField_std,p,(void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		MPntStdGetField_global_coord(material_point,&position);
		
		xcoord = position[0] * data->length_bar;
		ycoord = position[1] * data->length_bar;
		zcoord = position[2] * data->length_bar;
		
		phase = 0;
		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,phase);
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	
	ierr = MaterialPointGetAccess(db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mp_points; p++) {
		double kappa,H;
		
		ierr = MaterialPointGet_phase_index(mpX,p,&phase);CHKERRQ(ierr);
		
		kappa = data->diffusivity[0];
		H     = data->H[0];
		ierr = MaterialPointSet_diffusivity(mpX,p,kappa);CHKERRQ(ierr);
		ierr = MaterialPointSet_heat_source(mpX,p,H);CHKERRQ(ierr);
	}
	ierr = MaterialPointRestoreAccess(db,&mpX);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse3d_PerturbationThermal_Convection2d"
PetscBool DMDAVecTraverse3d_PerturbationThermal_Convection2d(PetscScalar position[], PetscScalar *val, void *ctx) 
{
	PetscScalar x,y,z;
	PetscReal  *coeffs;
	PetscBool  impose;
	PetscReal  L_bar,Tbot,age,Ttop,Ly,Lx,Oy,Ox,length_bar,diffusivity;
	PetscReal  cm_yr2m_s;
	
	cm_yr2m_s = 1.0e-2 / ( 365.25 * 24.0 * 60.0 * 60.0 ) ;
	
	/* get coordinates */
	x = position[0];
	y = position[1];
	z = position[2];
	/* fetch user data */
	coeffs = (PetscReal*)ctx;
	Ttop       = coeffs[0]; 
	Tbot       = coeffs[1]; 
	Oy         = coeffs[2];
	Ly         = coeffs[3];
	Ox         = coeffs[4];
	Lx         = coeffs[5];
	
	*val = Ttop+(Tbot-Ttop)*(Ly-y)/(Ly-Oy)+0.01*sin(2.0*3.1415926*(x-Ox)/(Lx-Ox));
	
	impose = PETSC_TRUE;
	
	return impose;
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_Thermal_Convection2d"
PetscErrorCode ModelApplyInitialSolution_Thermal_Convection2d(pTatinCtx c,Vec X, void *ctx)
{
	
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	DM           stokes_pack,dau,dap,daT;
	Vec          velocity,pressure,temperature;
	PetscReal    MeshMin[3],MeshMax[3],domain_height;
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	DMDAVecTraverse3d_InterpCtx                  IntpCtx;
	PetscReal      vals[6];
	PetscScalar    zero = 0.0;
	PhysCompEnergy energy;
	PetscBool      active_energy;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	stokes_pack = c->stokes_ctx->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
	ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
	
	ierr = DMDAGetBoundingBox(dau,MeshMin,MeshMax);CHKERRQ(ierr);
	domain_height = MeshMax[1] - MeshMin[1];
	
	HPctx.surface_pressure = 0.0;
	HPctx.ref_height = domain_height;
	HPctx.ref_N      = c->stokes_ctx->my-1;
	HPctx.grav       = 10.0;
	HPctx.rho        = data->rho[0];
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%e %d %e]]\n", HPctx.ref_height,HPctx.ref_N,HPctx.rho);
	
	ierr = DMDAVecTraverseIJK(dap,pressure,0,DMDAVecTraverseIJK_HydroStaticPressure_v2,     (void*)&HPctx);CHKERRQ(ierr); /* P = P0 + a.x + b.y + c.z, modify P0 (idx=0) */
	ierr = DMDAVecTraverseIJK(dap,pressure,2,DMDAVecTraverseIJK_HydroStaticPressure_dpdy_v2,(void*)&HPctx);CHKERRQ(ierr); /* P = P0 + a.x + b.y + c.z, modify b  (idx=2) */
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	
	if (active_energy) {
		ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(c,&temperature,PETSC_NULL);CHKERRQ(ierr);
		daT  = energy->daT;
		
		vals[0] = data->Ttop; //top temperature K
		vals[1] = data->Tbot; //bottom temperature K
		vals[2] = data->Oy;
		vals[3] = data->Ly;
		vals[4] = data->Ox;
		vals[5] = data->Lx;
		
		ierr = DMDAVecTraverse3d(daT,temperature,0, DMDAVecTraverse3d_PerturbationThermal_Convection2d, (void*)vals);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialStokesVariableMarkers_Thermal_Convection2d"
PetscErrorCode ModelApplyInitialStokesVariableMarkers_Thermal_Convection2d(pTatinCtx c,Vec X,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	DM                stokes_pack,dau,dap;
	PhysCompStokes    stokes;
	Vec               Uloc,Ploc;
	PetscScalar       *LA_Uloc,*LA_Ploc;
	DataField         PField;
	MaterialConst_MaterialType *truc;
	PetscInt        regionidx;
	PetscErrorCode  ierr;

	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	ierr = pTatin_EvaluateRheologyNonlinearities(c,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Thermal_Convection2d"
PetscErrorCode ModelApplyUpdateMeshGeometry_Thermal_Convection2d(pTatinCtx c,Vec X,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	PetscReal      step;
	PhysCompStokes stokes;
	Vec            velocity,pressure;
	DM             stokes_pack,dav,dap;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
#if 0	
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	//	ierr = UpdateMeshGeometry_FullLagrangian(dav,velocity,step);CHKERRQ(ierr);
	ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
#endif	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Thermal_Convection2d"
PetscErrorCode ModelOutput_Thermal_Convection2d(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	PetscBool      active_energy;
	DataBucket     materialpoint_db;
	DM             stokes_pack,dau,dap;
	Vec            velocity,pressure;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

#if 0	
	{
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
		//  Write out just the stokes variable?
		//  const int nf = 1;
		//  const MaterialPointField mp_prop_list[] = { MPField_Stokes };
		//
		//  Write out just std, stokes and plastic variables
		const int nf = 4;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
		char mp_file_prefix[256];
		
		sprintf(mp_file_prefix,"%s_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
#endif	
	/*		
	 {
	 const int                   nf = 4;
	 const MaterialPointVariable mp_prop_list[] = { MPV_viscosity, MPV_density,MPV_inv2stress, MPV_plastic_strain }; 
	 
	 ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
	 }	
	 */	

	/* standard viewer */
	ierr = pTatinContextValid_Energy(c,&active_energy);CHKERRQ(ierr);
	if (active_energy) {
		PhysCompEnergy energy;
		Vec            temperature;
		
		ierr = pTatinGetContext_Energy(c,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(c,&temperature,PETSC_NULL);CHKERRQ(ierr);
		
		ierr = pTatin3d_ModelOutput_Temperature_Energy(c,temperature,prefix);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Thermal_Convection2d"
PetscErrorCode ModelDestroy_Thermal_Convection2d(pTatinCtx c,void *ctx)
{
	ModelThermal_Convection2dCtx *data = (ModelThermal_Convection2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Thermal_Convection2d"
PetscErrorCode pTatinModelRegister_Thermal_Convection2d(void)
{
	ModelThermal_Convection2dCtx *data;
	pTatinModel    m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelThermal_Convection2dCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelThermal_Convection2dCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"convection2d");CHKERRQ(ierr);
	
	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_STOKES_VARIABLE_MARKERS,   (void (*)(void))ModelApplyInitialStokesVariableMarkers_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Thermal_Convection2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Thermal_Convection2d);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

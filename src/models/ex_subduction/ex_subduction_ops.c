/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      ex_subduction_ops.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id:$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

static char model_help[] = "\
|-------------------------------------------------------------------------------------------------------------------------- \n\
| \n\
| ModelName: ExSubduction [option -ptatin_model ex_subduction] \n\
| \n\
| ModelDescription: Prototype model for studying subduction \n\
| \n\
| ModelDetails: \n\
|   The default setup is defined by the model labelled as \"Case 1\" from the following paper: \n\
|    \"A benchmark comparison of spontaneous subduction models - Towards a free surface\"                                   \n\
|     Schmeling, Babeyko, Enns, Faccenna, Funiciello, Gerya, Golabek, Grigull, Kaus, Morra, Schmalholz, van Hunen           \n\
|		 Physics of the Earth and Planetary Interiors 171 (2008) 198–223                                                        \n\
|   See Figure 1 in this paper for the geometric and material paramters used.                                               \n\
|   Note that our \"Case 1\" setup employs a free surface, and not a `sticky-air' layer.                                    \n\
| \n\
| ModelOptions: \n\
|   -exsubduction_mode2d : Defines a two geometry for the domain (one-element thick) \n\
|   -exsubduction_finite_plate  : Indicated you want a plate which doesn't extend across the entire domain in z \n\
|   -exsubduction_finite_plate_width_factor VALUE : Width of slab is defined as VALUE * Lz  \n\
| \n\
|--------------------------------------------------------------------------------------------------------------------------\n";

/*
 [1] Model description:
 
 This is intended to be a prototype model for studying subduction processes.
 The default setup is defined by the model labelled as "Case 1" from the following paper: 
   "A benchmark comparison of spontaneous subduction models - Towards a free surface"
   Schmeling, Babeyko, Enns, Faccenna, Funiciello, Gerya, Golabek, Grigull, Kaus, Morra, Schmalholz, van Hunen
   Physics of the Earth and Planetary Interiors 171 (2008) 198–223
 See Figure 1 in this paper for the geometric and material paramters used.
 Note that our "Case 1" setup employs a free surface, and not a `sticky-air' layer.
 
 [2] Information about non-dimensionalization:
 
 We scaled the input (dimensional) using the following scales;
 
 option: --elv_eta =  1.0e23
 option: --elv_length =  1000.0e3
 option: --elv_velocity =  0.000000000317098
 <StokesScales: ND>
	 eta          1e+23
	 length       1000000.0
	 velocity     3.17098e-10
	 time         3.15359920277e+15
	 strain-rate  3.17098e-16
	 pressure     31709800.0
	 rhs-scale    0.0315359920277
 </StokesScales>
 <StokesScales: Convert user units into scaled units>
	 eta         -> eta         x  1e-23
	 length      -> length      x  1e-06
	 velocity    -> velocity    x  3153599202.77
	 time        -> time        x  3.17098e-16
	 strain-rate -> strain-rate x  3.15359920277e+15
	 pressure    -> pressure    x  3.15359920277e-08
 </StokesScales>
 
 
 [3] Example usage:
 
 <2d model>
 ./ptatin_driver_linear_ts.app -options_file test_option_files/test_basic_mfgmg.opts  -ptatin_model ex_subduction -output_path tsub -mx 92 -my 48 -dau_nlevels 3 -exsubduction_mode2d -fieldsplit_u_mg_levels_est_ksp_max_it 4 -fieldsplit_u_mg_levels_ksp_max_it 4

 <3d, finite plate model>
 ./ptatin_driver_linear_ts.app -options_file test_option_files/test_basic_mfgmg.opts  -ptatin_model ex_subduction -output_path tsub -mx 16 -my 12 -mz 16  -dau_nlevels 3 -fieldsplit_u_mg_levels_est_ksp_max_it 4 -fieldsplit_u_mg_levels_ksp_max_it 4 -lattice_layout_Nx 6 -lattice_layout_Ny 6 -lattice_layout_Nz 6 -exsubduction_finite_plate -exsubduction_finite_plate_width_factor 0.5 -mp_popctrl_np_upper 1000 -constant_dt 1.0e-5 -fieldsplit_u_ksp_rtol 1.0e-1 -ksp_atol 1.0e-4 -nsteps 3000
 
*/


const double char_eta=         1.0e+23;
const double char_length=      1000000.0;
const double char_velocity=    3.17098e-10;
const double char_time=        3.15359920277e+15;
const double char_strainrate=  3.17098e-16;
const double char_pressure=    31709800.0;
const double char_rhsscale=    0.0315359920277;



#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_utils.h"
#include "ptatin_models.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "energy_output.h"

#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "dmda_iterator.h"
#include "geometry_object.h"
#include "geometry_object_evaluator.h"
#include "rheology.h"
#include "material_constants.h"
#include "mesh_update.h"

#include "ex_subduction_ctx.h"

#define geom_eps 1.0e-8

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_ExSubduction"
PetscErrorCode ModelInitialize_ExSubduction(pTatinCtx c,void *ctx)
{
	ExSubductionCtx *data = (ExSubductionCtx*)ctx;
	PetscBool       flg;
	PetscReal       x0[3],Lx[3],cox[3];
	GeometryObject  G,A,B;
	DataBucket      materialconstants;
	PetscErrorCode  ierr;
	PetscBool       mode_2d = PETSC_FALSE;
	PetscBool       finite_plate = PETSC_FALSE;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	PetscPrintf(PETSC_COMM_WORLD,"%s\n",model_help);
	
	data->domain[0] = 3000.0e3 / char_length;
	data->domain[1] = 700.0e3 / char_length;
	data->domain[2] = 3000.0e3  / char_length;
	
	PetscOptionsGetBool(PETSC_NULL,"-exsubduction_mode2d",&mode_2d,PETSC_NULL);
	if (mode_2d) {
		c->mz           = 1;
		data->domain[2] = 50.0e3 / char_length;
		PetscOptionsInsertString("-da_refine_z 1");
	}

	data->dip = 90.0;
	PetscOptionsGetReal(PETSC_NULL,"-exsubduction_dip",&data->dip,PETSC_NULL);
	
	
	data->eta[ RegionId_Mantle ]      = 1.0e21/char_eta;     data->rho[ RegionId_Mantle ]      = 3200.0*char_rhsscale;
	data->eta[ RegionId_Slab ]        = 1.0e23/char_eta;     data->rho[ RegionId_Slab ]        = 3300.0*char_rhsscale;
	data->eta[ RegionId_LowerMantle ] = 1.0e21/char_eta;     data->rho[ RegionId_LowerMantle ] = 3200.0*char_rhsscale;


	/* define geometry here so everyone can use it */
	ierr = GeometryObjectCreate("mantle",&G);CHKERRQ(ierr);
	x0[0] = -geom_eps;                  x0[1] = -geom_eps;                  x0[2] = -geom_eps;
	Lx[0] = data->domain[0] + geom_eps; Lx[1] = data->domain[1] + geom_eps; Lx[2] = data->domain[2] + geom_eps; 
	ierr = GeometryObjectSetType_BoxCornerReference(G,x0,Lx);CHKERRQ(ierr);
	data->go[RegionId_Mantle] = G;

	ierr = GeometryObjectCreate("lowermantle",&G);CHKERRQ(ierr);
	x0[0] = -geom_eps;                  x0[1] = -geom_eps;                  x0[2] = -geom_eps;
	Lx[0] = data->domain[0] + geom_eps; Lx[1] = 100.0e3 / char_length;      Lx[2] = data->domain[2] + geom_eps; 
	ierr = GeometryObjectSetType_BoxCornerReference(G,x0,Lx);CHKERRQ(ierr);
	data->go[RegionId_LowerMantle] = G;
	
	ierr = GeometryObjectCreate("plate",&A);CHKERRQ(ierr);
	x0[0] = 1050.0e3 / char_length;   x0[1] = 600.0e3 / char_length;   x0[2] = 0.0;
	Lx[0] = 2000.0e3 / char_length;   Lx[1] = 100.0e3 / char_length;   Lx[2] = data->domain[2];
	ierr = GeometryObjectSetType_BoxCornerReference(A,x0,Lx);CHKERRQ(ierr);

	ierr = GeometryObjectCreate("slabtip",&B);CHKERRQ(ierr);
	x0[0] = 1050.0e3 / char_length;   x0[1] = 650.0e3 / char_length;   x0[2] = 0.5 * data->domain[2];
	Lx[0] = 100.0e3 / char_length;    Lx[1] = 400.0e3 / char_length;   Lx[2] = data->domain[2];
	ierr = GeometryObjectSetType_Box(B,x0,Lx);CHKERRQ(ierr);
	ierr = GeometryObjectRotate(B,ROTATE_AXIS_Z,(data->dip - 90.0));CHKERRQ(ierr);
	
	ierr = GeometryObjectCreate("slab",&G);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_SetOperationDefault(G,GeomSet_Union,A,B);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&A);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&B);CHKERRQ(ierr);
	data->go[RegionId_Slab] = G;
	data->ngo = 3;
	
	PetscOptionsGetBool(PETSC_NULL,"-exsubduction_finite_plate",&finite_plate,PETSC_NULL);
	if (finite_plate) {
		GeometryObject mask;
		PetscReal      factor;
		
		factor = 0.25;
		PetscOptionsGetReal(PETSC_NULL,"-exsubduction_finite_plate_width_factor",&factor,PETSC_NULL);
		
		
		ierr = GeometryObjectCreate("mask",&mask);CHKERRQ(ierr);
		x0[0] = 0.0e3 / char_length;      x0[1] =   0.0e3 / char_length;   x0[2] = 0.0;
		Lx[0] = 3000.0e3 / char_length;   Lx[1] = 700.0e3 / char_length;   Lx[2] = factor * data->domain[2];
		ierr = GeometryObjectSetType_BoxCornerReference(mask,x0,Lx);CHKERRQ(ierr);
		
		
		ierr = GeometryObjectCreate("slab",&G);CHKERRQ(ierr);
		ierr = GeometryObjectSetType_SetOperationDefault(G,GeomSet_Intersection,data->go[RegionId_Slab],mask);CHKERRQ(ierr);
		ierr = GeometryObjectDestroy(&data->go[RegionId_Slab]);CHKERRQ(ierr);
		ierr = GeometryObjectDestroy(&mask);CHKERRQ(ierr);
		data->go[RegionId_Slab] = G;
	}
	
	
	

	
	/* material properties */
	ierr = pTatinGetMaterialConstants(c,&materialconstants);CHKERRQ(ierr);
	MaterialConstantsSetDefaults(materialconstants);
	/* water */
	MaterialConstantsSetValues_MaterialType(materialconstants,0,VISCOUS_CONSTANT,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);		
	MaterialConstantsSetValues_ViscosityConst(materialconstants,0,1.0);

	/* lava types */
	MaterialConstantsSetValues_MaterialType(materialconstants,1,VISCOUS_CONSTANT,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);		
	MaterialConstantsSetValues_ViscosityConst(materialconstants,1,1.0);
	MaterialConstantsSetValues_MaterialType(materialconstants,2,VISCOUS_CONSTANT,PLASTIC_NONE,SOFTENING_NONE,DENSITY_CONSTANT);		
	MaterialConstantsSetValues_ViscosityConst(materialconstants,2,1.0);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ExSubduction_VelocityBC"
PetscErrorCode ExSubduction_VelocityBC(BCList bclist,DM dav,pTatinCtx c,ExSubductionCtx *data)
{
	PetscReal val_V;
	PetscReal basal_val_V[2];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	val_V = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&val_V);CHKERRQ(ierr);

	val_V = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&val_V);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&val_V);CHKERRQ(ierr);

	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&val_V);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&val_V);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_ExSubduction"
PetscErrorCode ModelApplyBoundaryConditionMG_ExSubduction(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
	PetscInt        n;
	ExSubductionCtx *data = (ExSubductionCtx*)ctx;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = ExSubduction_VelocityBC(bclist[n],dav[n],c,data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_ExSubduction"
PetscErrorCode ModelApplyBoundaryCondition_ExSubduction(pTatinCtx c,void *ctx)
{
	ExSubductionCtx *data = (ExSubductionCtx*)ctx;
	PhysCompStokes            stokes;
	DM                        stokes_pack,dav,dap;
	PhysCompEnergy            energy;
	PetscBool                 energy_active;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* velocity */
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = ExSubduction_VelocityBC(stokes->u_bclist,dav,c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_ExSubduction"
PetscErrorCode ModelApplyInitialMeshGeometry_ExSubduction(pTatinCtx c,void *ctx)
{
	ExSubductionCtx *data = (ExSubductionCtx*)ctx;
	PhysCompStokes       stokes;
	DM                   stokes_pack,dav,dap;
	PetscErrorCode       ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMDASetUniformCoordinates(dav,0.0,data->domain[0],0.0,data->domain[1],0.0,data->domain[2]);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_ExSubduction"
PetscErrorCode ModelApplyInitialMaterialGeometry_ExSubduction(pTatinCtx c,void *ctx)
{
	ExSubductionCtx *data = (ExSubductionCtx*)ctx;
	PetscInt       k;
	int            p,n_mp_points;
	DataBucket     material_point_db;
	MPAccess       mpX;
	int            inside;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	ierr = pTatinGetMaterialPoints(c,&material_point_db,PETSC_NULL);CHKERRQ(ierr);
	DataBucketGetSizes(material_point_db,&n_mp_points,0,0);
	ierr = MaterialPointGetAccess(material_point_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mp_points; p++) {
		double  *position;
		int     phase_index;
		
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		ierr = MaterialPointGet_global_coord(mpX,p,&position);CHKERRQ(ierr);

		inside = 0;
		ierr = GeometryObjectPointInside(data->go[0],position,&inside);CHKERRQ(ierr);
		if (inside == 1) {
			ierr = MaterialPointSet_phase_index(mpX,p,(int)RegionId_Mantle);CHKERRQ(ierr);
			ierr = MaterialPointSet_viscosity(mpX,p, data->eta[0]     );CHKERRQ(ierr);
			ierr = MaterialPointSet_density(mpX,p,  -data->rho[0]*9.8 );CHKERRQ(ierr);
		} else {
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Point appears to be outside domain");
		}
			
		for (k=1; k<data->ngo; k++) {
			inside = 0;
			ierr = GeometryObjectPointInside(data->go[k],position,&inside);CHKERRQ(ierr);
			if (inside == 1) {
				ierr = MaterialPointSet_phase_index(mpX,p,(int)k);CHKERRQ(ierr);
				ierr = MaterialPointSet_viscosity(mpX,p, data->eta[k]     );CHKERRQ(ierr);
				ierr = MaterialPointSet_density(mpX,p,  -data->rho[k]*9.8 );CHKERRQ(ierr);
			}
		}

	}
	ierr = MaterialPointRestoreAccess(material_point_db,&mpX);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_ExSubduction"
PetscErrorCode ModelApplyUpdateMeshGeometry_ExSubduction(pTatinCtx c,Vec X,void *ctx)
{
	ExSubductionCtx  *data = (ExSubductionCtx*)ctx;
	PetscErrorCode   ierr;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
  ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,c->dt);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_ExSubduction"
PetscErrorCode ModelOutput_ExSubduction(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ExSubductionCtx  *data = (ExSubductionCtx*)ctx;
	PetscBool        output_markers;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(c,X,prefix);CHKERRQ(ierr);

	output_markers = PETSC_TRUE;
/*
	if (output_markers) {
		DataBucket  materialpoint_db;
		const int   nf = 1;
		const MaterialPointField mp_prop_list[] = { MPField_Stokes };
		//  Write out just std, stokes and plastic variables
		//const int nf = 4;
		//const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
		char mp_file_prefix[256];
		
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
		sprintf(mp_file_prefix,"%s_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
*/
	if (output_markers) {
		ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_ExSubduction"
PetscErrorCode ModelDestroy_ExSubduction(pTatinCtx c,void *ctx)
{
	ExSubductionCtx *data = (ExSubductionCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_ExSubduction"
PetscErrorCode pTatinModelRegister_ExSubduction(void)
{
	ExSubductionCtx  *data;
	pTatinModel      m,model;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ExSubductionCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ExSubductionCtx));CHKERRQ(ierr);
	
	/* set initial values for model parameters */
	PetscMemzero(data->go,sizeof(GeometryObject)*10);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"ex_subduction");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_ExSubduction);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_ExSubduction);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_ExSubduction);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_ExSubduction);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_ExSubduction);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_ExSubduction);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_ExSubduction);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_ExSubduction);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

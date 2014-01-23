
#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "dmda_element_q2p1.h"
#include "mesh_update.h"
#include "ptatin_models.h"
#include "model_utils.h"
#include "ptatin_utils.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "geometry_object.h"

#include "iplus_ctx.h"


/* ------------- [iplus] C. Meriaux slab geometry ------------- */
/*
plate
 x0[0] = 0.0 / char_length;   x0[1] = 0.385 / char_length;   x0[2] = 0.15 / char_length;
 Lx[0] = 0.5 / char_length;   Lx[1] = 0.015 / char_length;   Lx[2] = 0.3 / char_length;

 tip
 x0[0] = 0.5 / char_length;   x0[1] = 0.385 / char_length;   x0[2] = 0.15 / char_length;
 Lx[0] = 0.06 / char_length;  Lx[1] = 0.015 / char_length;   Lx[2] = 0.3 / char_length;
 
*/
#undef __FUNCT__
#define __FUNCT__ "iPLUS_CreateSlabGeometry_Standard"
PetscErrorCode iPLUS_CreateSlabGeometry_Standard(iPLUSCtx *data)
{
	GeometryObject plate,tip;
	PetscReal      x0[3],Lx[3];
	PetscErrorCode ierr;
	
	ierr = GeometryObjectCreate("plate",&plate);CHKERRQ(ierr);
	x0[0] = 0.00;   x0[1] = 0.385;   x0[2] = 0.15;
	Lx[0] = 0.50;   Lx[1] = 0.015;   Lx[2] = 0.30;
	ierr = GeometryObjectSetType_BoxCornerReference(plate,x0,Lx);CHKERRQ(ierr);
	
	ierr = GeometryObjectCreate("slab_tip",&tip);CHKERRQ(ierr);
	x0[0] = 0.50;    x0[1] = 0.385+0.015*0.5;  x0[2] = 0.30;
	Lx[0] = 0.012;   Lx[1] = 0.12;             Lx[2] = 0.30;
	//ierr = GeometryObjectSetType_BoxCornerReference(tip,x0,Lx);CHKERRQ(ierr);
	//ierr = GeometryObjectRotate(tip,ROTATE_AXIS_Z,-34.0);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_Box(tip,x0,Lx);CHKERRQ(ierr);
	ierr = GeometryObjectRotate(tip,ROTATE_AXIS_Z,90-34.0);CHKERRQ(ierr);
	
	ierr = GeometryObjectCreate("slab",&data->slab_geometry);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_SetOperationDefault(data->slab_geometry,GeomSet_Union,plate,tip);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&plate);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&tip);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* ------------- Schellart, G3, 2008 slab geometry ------------- */
#undef __FUNCT__
#define __FUNCT__ "iPLUS_CreateSlabGeometry_Schellart_G3_2008"
PetscErrorCode iPLUS_CreateSlabGeometry_Schellart_G3_2008(iPLUSCtx *data)
{
	GeometryObject plate,tip;
	PetscReal      x0[3],Lx[3],tip_thickness,tip_length,tip_angle;
	PetscErrorCode ierr;
	
	/* Note - 2.5 cm of the slab is bent downwards into the mantle */
	ierr = GeometryObjectCreate("plate",&plate);CHKERRQ(ierr);
	x0[0] = 0.31+0.025;   x0[1] = 0.367;   x0[2] = 0.225;
	Lx[0] = 0.55-0.025;   Lx[1] = 0.013;   Lx[2] = 0.15;
	ierr = GeometryObjectSetType_BoxCornerReference(plate,x0,Lx);CHKERRQ(ierr);
	
	ierr = GeometryObjectCreate("slab_tip",&tip);CHKERRQ(ierr);
	tip_thickness  = 0.013;
	tip_length = 0.025;
	x0[0] = 0.31+0.025;        x0[1] = 0.38 - 0.5*tip_thickness;  x0[2] = 0.30;
	Lx[0] = 2.0*tip_length;    Lx[1] = tip_thickness;             Lx[2] = 0.15;
	ierr = GeometryObjectSetType_Box(tip,x0,Lx);CHKERRQ(ierr);
	/* Wouter indicates the slab tip has an angle between 15 - 30 degress [paragraph 16] */
	tip_angle = 30.0;
	PetscOptionsGetReal(PETSC_NULL,"-iplus_schellart_g3_2008_tip_angle",&tip_angle,PETSC_NULL);
	ierr = GeometryObjectRotate(tip,ROTATE_AXIS_Z, tip_angle);CHKERRQ(ierr);
	
	ierr = GeometryObjectCreate("slab",&data->slab_geometry);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_SetOperationDefault(data->slab_geometry,GeomSet_Union,plate,tip);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&plate);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&tip);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUS_DefineSlabMaterial"
PetscErrorCode iPLUS_DefineSlabMaterial(DM dav,DataBucket materialpoint_db,iPLUSCtx *data)
{
	MPAccess       mpX;
	int            p,n_mpoints;
	PetscErrorCode ierr;
	
	
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		int     inside_slab;
		double  *position;
		
		
		ierr = MaterialPointGet_global_coord(mpX,p,&position);CHKERRQ(ierr);
		
		inside_slab = 0;
		ierr = GeometryObjectPointInside(data->slab_geometry,position,&inside_slab);CHKERRQ(ierr);
		if (inside_slab == 1) {
			ierr = MaterialPointSet_phase_index(mpX,p,iPLUSMatSlab);CHKERRQ(ierr);
			ierr = MaterialPointSet_viscosity(mpX,p,data->slab_eta);CHKERRQ(ierr);
			ierr = MaterialPointSet_density(mpX,p,-GRAVITY*data->slab_rho);CHKERRQ(ierr);
		}
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUS_CreateSlabGeometry"
PetscErrorCode iPLUS_CreateSlabGeometry(iPLUSCtx *data)
{
	PetscErrorCode ierr;
	PetscBool      g3_slab = PETSC_FALSE;
	
	PetscOptionsGetBool(PETSC_NULL,"-iplus_slab_type_schellart_g3_2008",&g3_slab,PETSC_NULL);
	if (!g3_slab) {
		ierr = iPLUS_CreateSlabGeometry_Standard(data);CHKERRQ(ierr);
	} else {
		ierr = iPLUS_CreateSlabGeometry_Schellart_G3_2008(data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

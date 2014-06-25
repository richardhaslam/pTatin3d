
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
        /* set default tip length */
        tip_length = 0.025;
        /* check for input parameter -iplus_schellart_g3_2008_tip_length from file or command line */
        PetscOptionsGetReal(PETSC_NULL,"-iplus_schellart_g3_2008_tip_length",&tip_length,PETSC_NULL);
	
	/* Note - 2.5 cm of the slab is bent downwards into the mantle */
	ierr = GeometryObjectCreate("plate",&plate);CHKERRQ(ierr);
	x0[0] = 0.31+tip_length;   x0[1] = 0.367;   x0[2] = 0.225;
	Lx[0] = 0.55-tip_length;   Lx[1] = 0.013;   Lx[2] = 0.15;
	ierr = GeometryObjectSetType_BoxCornerReference(plate,x0,Lx);CHKERRQ(ierr);
	
	ierr = GeometryObjectCreate("slab_tip",&tip);CHKERRQ(ierr);
	tip_thickness  = 0.013;
	x0[0] = 0.31+tip_length;        x0[1] = 0.38 - 0.5*tip_thickness;  x0[2] = 0.30;
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

/* ------------- Li & Ribe, JGR, 2012 slab geometry ------------- */
#undef __FUNCT__
#define __FUNCT__ "iPLUS_CreateSlabGeometry_LiRibe_JGR_2012"
PetscErrorCode iPLUS_CreateSlabGeometry_LiRibe_JGR_2012(iPLUSCtx *data)
{
	GeometryObject plate,tail,nose;
	GeometryObject l_segment,otip,itip,annulus,clip_vert,clip_angle,clipper;

	PetscReal      x0[3],Lx[3],tip_thickness,tip_length,tip_angle;
	PetscReal      l,L,theta_0,h,W,x_0,y_0,z_0;
	PetscReal      R,r,xT[3];
	PetscBool      arcuate_slab_ends = PETSC_FALSE;
	PetscErrorCode ierr;

	x_0 = 0.31+0.025;   /* plate offset in x direction - subducting end */
	y_0 = 0.367; /* base of plate in y direction */
	z_0 = 0.225; /* plate edge offset in z direction */ 
	h  = 0.013;  /* plate thickness */
	L  = 0.55-0.025;   /* plate length */
	W  = 0.15;   /* plate width */
	
	l = 4.0 * h; /* input from paper, use default of l/h = 4 */
	PetscOptionsGetReal(PETSC_NULL,"-iplus_liribe_jgr_2012_l",&l,PETSC_NULL);
        

	theta_0 = 60.0; /* input from paper */
	PetscOptionsGetReal(PETSC_NULL,"-iplus_liribe_jgr_2012_theta0",&theta_0,PETSC_NULL);
	
	arcuate_slab_ends = PETSC_FALSE;
	PetscOptionsGetBool(PETSC_NULL,"-iplus_liribe_jgr_2012_arcuate_slab_ends",&arcuate_slab_ends,PETSC_NULL);
	
	
	/* horizontal plate ------------------------- */	
	ierr = GeometryObjectCreate("plate",&plate);CHKERRQ(ierr);
	x0[0] = x_0;  x0[1] = y_0;   x0[2] = z_0;
	Lx[0] = L;    Lx[1] = h;     Lx[2] = W;
	ierr = GeometryObjectSetType_BoxCornerReference(plate,x0,Lx);CHKERRQ(ierr);

	
	/* arcuate slab tail ------------------------- */
	ierr = GeometryObjectCreate("slab_tail",&tail);CHKERRQ(ierr);
	x0[0] = x_0 + L;
	x0[1] = y_0 + 0.5*h;
	x0[2] = z_0 + 0.5*W;
	ierr = GeometryObjectSetType_Cylinder(tail,x0,0.5*h,W,ROTATE_AXIS_Z);CHKERRQ(ierr);
	
	
	/* arcuate trench - l_segment ------------------------- */
	r = l / (theta_0 * M_PI / 180.0); /* theta.r = segment_length */
	R = r + h; /* outer radius = inner + slab_thickness */
	
	/* a) create annulus */
	ierr = GeometryObjectCreate("otip",&otip);CHKERRQ(ierr);
	x0[0] = x_0;
	x0[1] = y_0 + h - R;
	x0[2] = z_0 + 0.5*W;
	ierr = GeometryObjectSetType_Cylinder(otip,x0,R,W,ROTATE_AXIS_Z);CHKERRQ(ierr);

	ierr = GeometryObjectCreate("itip",&itip);CHKERRQ(ierr);
	x0[0] = x_0;
	x0[1] = y_0 + h - R;
	x0[2] = z_0 + 0.5*W;
	ierr = GeometryObjectSetType_Cylinder(itip,x0,r,W,ROTATE_AXIS_Z);CHKERRQ(ierr);

	ierr = GeometryObjectCreate("annulus",&annulus);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_SetOperationDefault(annulus,GeomSet_Complement,otip,itip);CHKERRQ(ierr);
	
	/* b) create clippers */
	ierr = GeometryObjectCreate("clip_vert",&clip_vert);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_HalfSpace(clip_vert,x0,SIGN_NEGATIVE,ROTATE_AXIS_X);CHKERRQ(ierr);

	ierr = GeometryObjectCreate("clip_angle",&clip_angle);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_HalfSpace(clip_angle,x0,SIGN_POSITIVE,ROTATE_AXIS_X);CHKERRQ(ierr);
	ierr = GeometryObjectRotate(clip_angle,ROTATE_AXIS_Z,theta_0);CHKERRQ(ierr);

	ierr = GeometryObjectCreate("clipper",&clipper);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_SetOperationDefault(clipper,GeomSet_Intersection,clip_vert,clip_angle);CHKERRQ(ierr);
	
	/* merge objects */
	ierr = GeometryObjectCreate("l_segment",&l_segment);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_SetOperationDefault(l_segment,GeomSet_Intersection,annulus,clipper);CHKERRQ(ierr);

	/* arcuate slab nose ------------------------- */
	/* radius = R - 0.5h */
	xT[0] = x0[0] - (R-0.5*h)*sin(theta_0 * M_PI/180.0);
	xT[1] = x0[1] + (R-0.5*h)*cos(theta_0 * M_PI/180.0);
	xT[2] = x0[2];
	
	ierr = GeometryObjectCreate("slab_nose",&nose);CHKERRQ(ierr);
	ierr = GeometryObjectSetType_Cylinder(nose,xT,0.5*h,W,ROTATE_AXIS_Z);CHKERRQ(ierr);

	
	/* merge plate + slab (l_segment) */
	if (!arcuate_slab_ends) {
		ierr = GeometryObjectCreate("slab",&data->slab_geometry);CHKERRQ(ierr);
		ierr = GeometryObjectSetType_SetOperationDefault(data->slab_geometry,GeomSet_Union,plate,l_segment);CHKERRQ(ierr);
		
	} else {
		GeometryObject tmp_slab,tmp_slab1;
		
		/* merge plate + slab */
		ierr = GeometryObjectCreate("slab_1",&tmp_slab);CHKERRQ(ierr);
		ierr = GeometryObjectSetType_SetOperationDefault(tmp_slab,GeomSet_Union,plate,l_segment);CHKERRQ(ierr);

		/* merge nose */
		ierr = GeometryObjectCreate("slab_1a",&tmp_slab1);CHKERRQ(ierr);
		ierr = GeometryObjectSetType_SetOperationDefault(tmp_slab1,GeomSet_Union,tmp_slab,nose);CHKERRQ(ierr);
		
		/* merge tail */
		ierr = GeometryObjectCreate("slab",&data->slab_geometry);CHKERRQ(ierr);
		ierr = GeometryObjectSetType_SetOperationDefault(data->slab_geometry,GeomSet_Union,tmp_slab1,tail);CHKERRQ(ierr);

		ierr = GeometryObjectDestroy(&tmp_slab);CHKERRQ(ierr);
		ierr = GeometryObjectDestroy(&tmp_slab1);CHKERRQ(ierr);
	}
	
	
	ierr = GeometryObjectDestroy(&tail);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&nose);CHKERRQ(ierr);

	ierr = GeometryObjectDestroy(&l_segment);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&annulus);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&clip_vert);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&clip_angle);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&clipper);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&otip);CHKERRQ(ierr);
	ierr = GeometryObjectDestroy(&itip);CHKERRQ(ierr);
		
	ierr = GeometryObjectDestroy(&plate);CHKERRQ(ierr);
	
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
	PetscBool      wouter_slab  = PETSC_FALSE;
	PetscBool      ribe_slab = PETSC_FALSE;
	
	PetscOptionsGetBool(PETSC_NULL,"-iplus_slab_type_schellart_g3_2008",&wouter_slab,PETSC_NULL);
	PetscOptionsGetBool(PETSC_NULL,"-iplus_slab_type_liribe_jgr_2012",&ribe_slab,PETSC_NULL);
	if (wouter_slab) {
		ierr = iPLUS_CreateSlabGeometry_Schellart_G3_2008(data);CHKERRQ(ierr);
	}
	if (ribe_slab) {
		ierr = iPLUS_CreateSlabGeometry_LiRibe_JGR_2012(data);CHKERRQ(ierr);
	}
	
	if ( (!wouter_slab) && (!ribe_slab) ) {
		ierr = iPLUS_CreateSlabGeometry_Standard(data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}


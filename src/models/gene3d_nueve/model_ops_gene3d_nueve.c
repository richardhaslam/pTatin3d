

#define _GNU_SOURCE
#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "phase_map.h"

#include "model_gene3d_ctx.h"

#define MODEL_NAME "Gene3DNueve_"

PetscErrorCode ReportOptionMissing(const char optprefix[],const char description[],const char opt_name[],const char defaults[])
{
	char fulloptionname[512];
	
	if (optprefix) {
		sprintf(fulloptionname,"-%s%s",optprefix,&opt_name[1]);
	} else {
		sprintf(fulloptionname,"%s",opt_name);
	}
	if (defaults) {
		PetscPrintf(PETSC_COMM_WORLD,"Missing Option: Expected user to specify \"%s\" via option %s : defaults {%s}\n",description,fulloptionname,defaults);
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Missing Option: Expected user to specify \"%s\" via option %s\n",description,fulloptionname);
	}
	SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"<< Essential User Command Line Option Missing >>\n");	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelGene3DNueve_ParseMaterialConstants_RHEOLOGY_VISCOUS"
PetscErrorCode ModelGene3DNueve_ParseMaterialConstants_RHEOLOGY_VISCOUS(DataBucket material_constants,const char model_name[],const PetscInt region_id)
{
	PetscErrorCode ierr;

	PetscFunctionBegin;
	ierr = MaterialConstantsSetFromOptions_ViscosityConst(material_constants,model_name,region_id,PETSC_TRUE);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelGene3DNueve_ParseMaterialConstants_RHEOLOGY_VP_STD"
PetscErrorCode ModelGene3DNueve_ParseMaterialConstants_RHEOLOGY_VP_STD(DataBucket material_constants,const char model_name[],const PetscInt region_id)
{
	PetscInt visc_type,plastic_type,softening_type;
	PetscBool found;
	char option_name[256];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* viscous type */
	sprintf(option_name,"-viscous_%d",region_id);
  ierr = PetscOptionsGetInt(model_name,option_name,&visc_type,&found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER,"Expected user to provide value for viscous type for region %D; via -viscous_%D\n",region_id,region_id);
	}
	
	/* plastic type */
	sprintf(option_name,"-plastic_%d",region_id);
  ierr = PetscOptionsGetInt(model_name,option_name,&plastic_type,&found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER,"Expected user to provide value for plastic type for region %D; via -plastic_%D\n",region_id,region_id);
	}
	
  /* softening  type */
	sprintf(option_name,"-softening_%d",region_id);
  ierr = PetscOptionsGetInt(model_name,option_name,&softening_type,&found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER,"Expected user to provide value for softening type for region %D; via -softening_%D\n",region_id,region_id);
	}
	
  /* density type */

	
	
	switch (visc_type) {
			
		case VISCOUS_CONSTANT:
			ierr = MaterialConstantsSetFromOptions_ViscosityConst(material_constants,model_name,region_id,PETSC_TRUE);CHKERRQ(ierr);
			break;
			
		case VISCOUS_FRANKK:
			printf("TODO\n");
			break;
			
		case VISCOUS_ARRHENIUS:
			printf("TODO\n");
			break;
	}
		
	switch (plastic_type) {
			
		case PLASTIC_NONE:
			printf("TODO\n");
			break;
			
		case PLASTIC_MISES:
			printf("TODO\n");
			break;
			
		case PLASTIC_DP:
			printf("TODO\n");
			break;
	}
	
	switch (softening_type) {
			
		case SOFTENING_NONE:
			printf("TODO\n");
			break;
			
		case SOFTENING_LINEAR:
			printf("TODO\n");
			break;
			
		case SOFTENING_EXPONENTIAL:
			printf("TODO\n");
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Gene3DNueve"
PetscErrorCode ModelInitialize_Gene3DNueve(pTatinCtx c,void *ctx)
{
  ModelGene3DNueveCtx    *data = (ModelGene3DNueveCtx*)ctx;
  PetscBool         flg, found;
  char              *option_name;
  PetscInt          r,nregions,opts_bcs,rheol_type;
  PetscErrorCode    ierr;
	
  PetscFunctionBegin;
	
	
  PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
  /* model geometry */
  PetscPrintf(PETSC_COMM_WORLD,"reading model initial geometry from options\n");
  ierr = PetscOptionsGetInt(PETSC_NULL,"-initial_region_geom_type",(PetscInt*)&data->initial_geom,&found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Type of region index initialisation","-initial_region_geom_type","1=pmap, 2=CAD");CHKERRQ(ierr);
	}
	
  /* box geometry */
  PetscPrintf(PETSC_COMM_WORLD, "reading box geometry from options\n");
	
  ierr = PetscOptionsGetReal(PETSC_NULL, "-Lx", &data->Lx, &found);CHKERRQ(ierr);
  if (found == PETSC_FALSE) {
		ierr = ReportOptionMissing(0,"Model length in x","-Lx",0);CHKERRQ(ierr);
	}
	
  ierr = PetscOptionsGetReal(PETSC_NULL, "-Ly", &data->Ly, &found);CHKERRQ(ierr);
  if (found == PETSC_FALSE) {
		ierr = ReportOptionMissing(0,"Model length in y","-Ly",0);CHKERRQ(ierr);
	}
	
  ierr = PetscOptionsGetReal(PETSC_NULL, "-Lz", &data->Lz, &found);CHKERRQ(ierr);
  if (found == PETSC_FALSE) {
		ierr = ReportOptionMissing(0,"Model length in z","-Lz",0);CHKERRQ(ierr);
	}
	
	
  ierr = PetscOptionsGetReal(PETSC_NULL, "-Ox", &data->Ox, &found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Model origin in x","-Ox",0);CHKERRQ(ierr);
	}
	
  ierr = PetscOptionsGetReal(PETSC_NULL, "-Oy", &data->Oy, &found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Model origin in y","-Oy",0);CHKERRQ(ierr);
	}
	
  ierr = PetscOptionsGetReal(PETSC_NULL, "-Oz", &data->Oz, &found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Model origin in z","-Oz",0);CHKERRQ(ierr);
	}
	
	
  /* bc type */
	opts_bcs = GENEBC_FreeSlip;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-bc_type",&opts_bcs,&found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Boundary condition type","-bc_type",0);CHKERRQ(ierr);
	}
  data->boundary_conditon_type = opts_bcs;
	
		
  /* material properties */
  PetscPrintf(PETSC_COMM_WORLD,"reading material properties from options\n");
  ierr = PetscOptionsGetInt(MODEL_NAME, "-nphase", &nregions, &found);CHKERRQ(ierr);
  if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(MODEL_NAME,"Maximum number of material regions","-nphase",0);CHKERRQ(ierr);
	}
	
	/* set the active phases on the rheo structure */
	DataBucketSetSizes(c->material_constants,nregions,-1);
	
	ierr = PetscOptionsGetInt(MODEL_NAME,"-rheol",&rheol_type,&found);CHKERRQ(ierr);
	if (found == PETSC_FALSE){
		ierr = ReportOptionMissing(MODEL_NAME,"Rheology type","-rheol","0=viscous , 4=visco_plastic");CHKERRQ(ierr);
	}
	
	
	switch  (rheol_type) {
		
		case RHEOLOGY_VISCOUS:
			for (r=0; r<nregions; r++) {
				/* this should live in rheology.c */
				ierr = ModelGene3DNueve_ParseMaterialConstants_RHEOLOGY_VISCOUS(c->material_constants,MODEL_NAME,r);CHKERRQ(ierr);
			}
			break;

		case RHEOLOGY_VP_STD:
			for (r=0; r<nregions; r++) {
				/* this should live in rheology.c */
				ierr = ModelGene3DNueve_ParseMaterialConstants_RHEOLOGY_VP_STD(c->material_constants,MODEL_NAME,r);CHKERRQ(ierr);
			}
			break;
			
		case RHEOLOGY_VISCO_PLASTIC:
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "-rheol RHEOLOGY_VISCO_PLASTIC");
			break;
			
		case RHEOLOGY_VISCO_PLASTIC_STRAIN_WEAKENING:
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "-rheol RHEOLOGY_VISCO_PLASTIC_STRAIN_WEAKENING");
			break;
			
		case RHEOLOGY_VISCO_ELASTIC_PLASTIC:
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "-rheol RHEOLOGY_VISCO_ELASTIC_PLASTIC");
			break;
			
	}
	
  PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Gene3DNueve"
PetscErrorCode ModelApplyBoundaryCondition_Gene3DNueve(pTatinCtx user,void *ctx)
{
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscScalar zero = 0.0;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
  switch (data->boundary_conditon_type) {
			
    case GENEBC_FreeSlip:
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      break;
			
    case GENEBC_NoSlip:
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_KMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      break;
			
    case GENEBC_FreeSlipFreeSurface:
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      break;
			
    case GENEBC_NoSlipFreeSurface:
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_IMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      ierr = DMDABCListTraverse3d(user->stokes_ctx->u_bclist,user->stokes_ctx->dav, DMDABCList_JMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
      break;
			
    default:
      break;
	}
	
  PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Gene3DNueve"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Gene3DNueve(pTatinCtx c,void *ctx)
{
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
  PetscPrintf(PETSC_COMM_WORLD, "  NOT IMPLEMENTED \n", __FUNCT__);
	
  PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Gene3DNueve"
PetscErrorCode ModelApplyInitialMeshGeometry_Gene3DNueve(pTatinCtx c,void *ctx)
{
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
  ierr =	DMDASetUniformCoordinates(c->stokes_ctx->dav, data->Ox, data->Lx, data->Oy, data->Ly, data->Oz, data->Lz); CHKERRQ(ierr);
	
  PetscFunctionReturn (0);
}

/* Define eta,rho on material points */
#undef __FUNCT__
#define __FUNCT__ "Gene3DNueve_MaterialPointSetInitialStokesVariables"
PetscErrorCode Gene3DNueve_MaterialPointSetInitialStokesVariables(pTatinCtx c,void *ctx)
{
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscInt e, p, n_mp_points;
  DataBucket material_point_db;
  DataBucket material_constants_db;
  DataField PField_std, PField_stokes, PField_ViscConst;
	MaterialConst_ViscosityConst *material_visc_const;
  int region_id, i;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
	
  material_point_db = c->materialpoint_db;
  material_constants_db = c->material_constants;

	/* Get material point data */
  DataBucketGetDataFieldByName(material_point_db, MPntStd_classname, &PField_std);
  DataFieldGetAccess(PField_std);
  DataFieldVerifyAccess(PField_std, sizeof (MPntStd));

  DataBucketGetDataFieldByName (material_point_db, MPntPStokes_classname, &PField_stokes);
  DataFieldGetAccess(PField_stokes);
  DataFieldVerifyAccess(PField_stokes, sizeof (MPntPStokes));
	
  DataBucketGetSizes(material_point_db,&n_mp_points,0,0);

	/* Get material constants data */
  DataBucketGetDataFieldByName(material_constants_db, MaterialConst_ViscosityConst_classname, &PField_ViscConst);
	material_visc_const = (MaterialConst_ViscosityConst*)PField_ViscConst->data;
	
  for (p=0; p<n_mp_points; p++) {
		MPntStd     *mpprop_std;
		MPntPStokes *mpprop_stokes;
		
		DataFieldAccessPoint(PField_std, p, (void**)&mpprop_std);
		DataFieldAccessPoint(PField_stokes, p, (void**)&mpprop_stokes);

		MPntStdGetField_phase_index (mpprop_std, &region_id);
		
		MPntPStokesSetField_eta_effective(mpprop_stokes, material_visc_const[region_id].eta0 );
		MPntPStokesSetField_density(mpprop_stokes,       2.0 );
	}
	
  DataFieldRestoreAccess(PField_std);
  DataFieldRestoreAccess(PField_stokes);
	
  PetscFunctionReturn (0);
}

/* define phase index on material points from a map file extruded in z direction */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSetRegionIndexFromMap"
PetscErrorCode MaterialPointSetRegionIndexFromMap(pTatinCtx c,void *ctx)
{
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscErrorCode ierr;
  PhaseMap phasemap;
  PetscInt p, n_mp_points,dir_0,dir_1,direction;
  DataBucket db;
  DataField PField_std;
  int phase_init, phase, phase_index, is_valid;
  char map_file[PETSC_MAX_PATH_LEN], *name;
  PetscBool flg;
  
	PetscFunctionBegin;
  PetscPrintf (PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
  ierr = PetscOptionsGetString(MODEL_NAME,"-map_file",map_file,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (flg == PETSC_FALSE) {
		ierr = ReportOptionMissing(MODEL_NAME,"Name of pmap file","-map_file","HINT: Don't include the extension .pmap in the file name");CHKERRQ(ierr);
	}
  ierr = PetscOptionsGetInt(MODEL_NAME,"-extrude_dir",&direction,&flg);CHKERRQ(ierr);
 if (flg == PETSC_FALSE) {
	 ierr = ReportOptionMissing(MODEL_NAME,"Direction to extrude the pmap","-extrude_dir","0=x, 1=y, 2=z");CHKERRQ(ierr);
	}	
	
  switch (direction) {
    case 0:
		{
			dir_0 = 2;
			dir_1 = 1;
    } 
			break;  
			
    case 1:
		{
			dir_0 = 0;
			dir_1 = 2;
    } 
			break;   
			
    case 2:
		{
			dir_0 = 0;
			dir_1 = 1;
    } 
			break;  
  } 
	
  asprintf(&name,"./inputdata/%s.pmap",map_file);
  PhaseMapLoadFromFile(name,&phasemap);
  free(name);
	
  /* define properties on material points */
  db = c->materialpoint_db;
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  DataFieldVerifyAccess(PField_std,sizeof (MPntStd));
	
  DataBucketGetSizes(db,&n_mp_points,0,0);
	
  for (p=0; p<n_mp_points; p++) {
		MPntStd *material_point;
		double position2D[2],*pos;
		

		DataFieldAccessPoint(PField_std, p, (void **) &material_point);
		MPntStdGetField_global_coord(material_point,&pos);

		position2D[0] = pos[dir_0];
		position2D[1] = pos[dir_1];
		
		MPntStdGetField_phase_index(material_point, &phase_init);
		
		PhaseMapGetPhaseIndex(phasemap, position2D, &phase_index);
		
		PhaseMapCheckValidity(phasemap, phase_index, &is_valid);
		//PetscPrintf(PETSC_COMM_WORLD,"Phase index : %d  is_valid %d \n", phase_index,is_valid);
		
		if (is_valid == 1) {			/* point located in the phase map */
			phase = phase_index;
		} else {
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,"marker outside the domain\n your phasemap is smaller than the domain \n please check your parameters and retry");
		}
		MPntStdSetField_phase_index(material_point, phase);
		
	}
  PhaseMapDestroy(&phasemap);	
  DataFieldRestoreAccess(PField_std);
	
  PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointInitializeRegionIndex"
PetscErrorCode MaterialPointInitializeRegionIndex(DataBucket db,const int init_region_id)
{
  PetscInt           p,n_mp_points;
  DataField          PField_std;
  PetscErrorCode     ierr;
	
  PetscFunctionBegin;
	
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
  DataBucketGetSizes(db,&n_mp_points,0,0);
	
  for (p=0; p<n_mp_points; p++) {
		int     region_index;
		MPntStd *material_point;

		DataFieldAccessPoint(PField_std,p,(void**)&material_point);

		region_index = init_region_id;
		MPntStdSetField_phase_index(material_point,region_index);
	}
  DataFieldRestoreAccess(PField_std);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointCheckRegionIndexBounds"
PetscErrorCode MaterialPointCheckRegionIndexBounds(DataBucket db,PetscInt lower_regionid,PetscInt upper_regionid)
{
  PetscInt           p,n_mp_points;
  DataField          PField_std;
  PetscErrorCode     ierr;
  PetscFunctionBegin;
	
  DataBucketGetDataFieldByName (db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
  DataBucketGetSizes(db,&n_mp_points,0,0);
	
  for (p=0; p<n_mp_points; p++) {
		int    region_index;
		double *pos;
		MPntStd *material_point;

		DataFieldAccessPoint(PField_std,p,(void**)&material_point);
		
		MPntStdGetField_phase_index(material_point,&region_index);
		MPntStdGetField_global_coord(material_point,&pos);
		if (region_index < lower_regionid) {
			PetscPrintf(PETSC_COMM_SELF,"RegionId of material point %D below lower bound (%D): marker coor (%1.4e,%1.4e,%1.4e)\n",p,lower_regionid,pos[0],pos[1],pos[2]);
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Material point region_id is uninitialized");
		}
		if (region_index >= upper_regionid) {
			PetscPrintf(PETSC_COMM_SELF,"RegionId of material point %D above upper bound (%D): marker coor (%1.4e,%1.4e,%1.4e)\n",p,upper_regionid,pos[0],pos[1],pos[2]);
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Material point region_id is uninitialized");
		}
	}
  DataFieldRestoreAccess(PField_std);
	
  PetscFunctionReturn(0);
}

/*
These models are required to
 1) set the value c->rheology_constants->nphases_active => although this seems to be done in ModelInitialize_Gene3DNueve()
 2) all markers are assigned a phase index between [0 -- nphases_active-1]
*/
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Gene3DNueve"
PetscErrorCode ModelApplyInitialMaterialGeometry_Gene3DNueve(pTatinCtx c,void *ctx)
{
	int                 nregions;
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscErrorCode      ierr;
  PetscFunctionBegin;
	
	/* Onitalize all region indices on markers to be -1 */
	ierr = MaterialPointInitializeRegionIndex(c->materialpoint_db,-1);CHKERRQ(ierr);
  switch (data->initial_geom) {
      
    case 0: /* Layer cake <obsolete> */
		{
		}
      break;
      
    case 1: /* Extrude map along x,y,z directions */
		{
			ierr = MaterialPointSetRegionIndexFromMap(c,ctx);CHKERRQ(ierr);
		}
      break;
      
    case 2: /* Read from CAD file */
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Reading region index from CAD is not implemented yet\n");
		}
      break;
	}
	/* Check all phase indices are between [0---rheo->max_phases-1] */
	DataBucketGetSizes(c->material_constants,&nregions,0,0);
	ierr = MaterialPointCheckRegionIndexBounds(c->materialpoint_db,0,nregions);CHKERRQ(ierr);

	/* Set eta0, rho0 */
  ierr = Gene3DNueve_MaterialPointSetInitialStokesVariables(c,ctx);CHKERRQ(ierr);
	
  PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Gene3DNueve"
PetscErrorCode ModelApplyUpdateMeshGeometry_Gene3DNueve(pTatinCtx c,Vec X,void *ctx)
{
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
  PetscPrintf(PETSC_COMM_WORLD, "  NOT IMPLEMENTED \n", __FUNCT__);
	
  PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Gene3DNueve"
PetscErrorCode ModelOutput_Gene3DNueve(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
  ierr = pTatin3d_ModelOutput_MPntStd(c,prefix); CHKERRQ(ierr);
	
  PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Gene3DNueve"
PetscErrorCode ModelDestroy_Gene3DNueve(pTatinCtx c,void *ctx)
{
  ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
  /* Free contents of structure */
	
  /* Free structure */
  ierr = PetscFree(data);CHKERRQ(ierr);
	
  PetscFunctionReturn (0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Gene3DNueve"
PetscErrorCode pTatinModelRegister_Gene3DNueve(void)
{
  ModelGene3DNueveCtx *data;
  pTatinModel m, model;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
	
  /* Allocate memory for the data structure for this model */
  ierr = PetscMalloc(sizeof(ModelGene3DNueveCtx),&data);CHKERRQ(ierr);
  ierr = PetscMemzero(data,sizeof(ModelGene3DNueveCtx));CHKERRQ(ierr);
	
  /* register user model */
  ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	
  /* Set name, model select via -ptatin_model NAME */
  ierr = pTatinModelSetName(m,"Gene3DNueve");CHKERRQ(ierr);
	
  /* Set model data */
  ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
  /* Set function pointers */
  ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_INIT,                  (void (*)(void)) ModelInitialize_Gene3DNueve); CHKERRQ(ierr);
  ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_BC,              (void (*)(void)) ModelApplyBoundaryCondition_Gene3DNueve); CHKERRQ(ierr);
  ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void)) ModelApplyMaterialBoundaryCondition_Gene3DNueve);CHKERRQ(ierr);
  ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void)) ModelApplyInitialMeshGeometry_Gene3DNueve);CHKERRQ(ierr);
  ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void)) ModelApplyInitialMaterialGeometry_Gene3DNueve);CHKERRQ(ierr);
  ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void)) ModelApplyUpdateMeshGeometry_Gene3DNueve);CHKERRQ(ierr);
  ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_OUTPUT,                (void (*)(void)) ModelOutput_Gene3DNueve);CHKERRQ(ierr);
  ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_DESTROY,               (void (*)(void)) ModelDestroy_Gene3DNueve); CHKERRQ(ierr);
	
  /* Insert model into list */
  ierr = pTatinModelRegister(m); CHKERRQ(ierr);
	
  PetscFunctionReturn (0);
}

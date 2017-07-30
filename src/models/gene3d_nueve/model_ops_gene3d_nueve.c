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
 **    filename:   model_ops_gene3d_nueve.c
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

/*
 Developed by Laetitia Le Pourhiet [laetitia.le_pourhiet@upmc.fr]
 */

#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_bcs.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "data_bucket.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "cartgrid.h"

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
	SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"<< Essential User Command Line Option Missing >>");
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
	PetscInt                   visc_type,plastic_type,softening_type;
	PetscBool                  found;
	char                       option_name[256];
	MaterialConst_MaterialType *data;
	DataField                  PField;
	PetscErrorCode             ierr;
	
	PetscFunctionBegin;
    
    
	/* Load in the material types */
	DataBucketGetDataFieldByName(material_constants,MaterialConst_MaterialType_classname,&PField);
	DataFieldGetAccess(PField);
	DataFieldAccessPoint(PField,region_id,(void**)&data);
	
	/* viscous type */
	sprintf(option_name,"-viscous_%d",region_id);
    ierr = PetscOptionsGetInt(NULL,model_name,option_name,&visc_type,&found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		//SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER,"Expected user to provide value for viscous type for region %D; via -viscous_%D\n",region_id,region_id);
		ierr = ReportOptionMissing(model_name,"Viscous type for region",option_name,"0 = CONSTANT, 1 = FRANKK, 2 = ARRHENIUS");CHKERRQ(ierr);
	}
    
	MaterialConst_MaterialTypeSetField_visc_type(data,visc_type); /* use setter */
	
	/* plastic type */
	sprintf(option_name,"-plastic_%d",region_id);
    ierr = PetscOptionsGetInt(NULL,model_name,option_name,&plastic_type,&found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		//SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER,"Expected user to provide value for plastic type for region %D; via -plastic_%D\n",region_id,region_id);
		ierr = ReportOptionMissing(model_name,"Plastic type for region",option_name,"0 = NONE, 1 = MISES, 2 = DP");CHKERRQ(ierr);
	}
	
	MaterialConst_MaterialTypeSetField_plastic_type(data,plastic_type);
	
    /* softening  type */
	sprintf(option_name,"-softening_%d",region_id);
    ierr = PetscOptionsGetInt(NULL,model_name,option_name,&softening_type,&found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		//SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER,"Expected user to provide value for softening type for region %D; via -softening_%D\n",region_id,region_id);
		ierr = ReportOptionMissing(model_name,"Softening type for region",option_name,"0 = NONE, 1 = LINEAR, 2 = EXPONENTIAL");CHKERRQ(ierr);
	}
	
	MaterialConst_MaterialTypeSetField_softening_type(data,softening_type);
	
    /* density type */
    
	
	DataFieldRestoreAccess(PField);
	
	
	switch (visc_type) {
			
		case VISCOUS_CONSTANT:
			ierr = MaterialConstantsSetFromOptions_ViscosityConst(material_constants,model_name,region_id,PETSC_TRUE);CHKERRQ(ierr);
			break;
			
		case VISCOUS_FRANKK:
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"-viscous FRANKK currently unsupported");
			break;
			
		case VISCOUS_ARRHENIUS:
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"-viscous ARRHENIUS currently unsupported");
			break;
	}
    
	switch (plastic_type) {
			
		case PLASTIC_NONE:
			break;
			
		case PLASTIC_MISES:
			ierr = MaterialConstantsSetFromOptions_PlasticMises(material_constants,model_name,region_id,PETSC_TRUE);CHKERRQ(ierr);
			break;
			
		case PLASTIC_DP:
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"-plastic DP currently unsupported");
			break;
	}
	
	switch (softening_type) {
			
		case SOFTENING_NONE:
			break;
			
		case SOFTENING_LINEAR:
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"-softening LINEAR currently unsupported");
			break;
			
		case SOFTENING_EXPONENTIAL:
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"-softening EXPONENTIAL currently unsupported");
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Gene3DNueve"
PetscErrorCode ModelInitialize_Gene3DNueve(pTatinCtx c,void *ctx)
{
    ModelGene3DNueveCtx    *data = (ModelGene3DNueveCtx*)ctx;
    PetscBool         found;
    PetscInt          r,nregions,opts_bcs,rheol_type;
    PetscErrorCode    ierr;
	
    PetscFunctionBegin;
	
	
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    /* model geometry */
    PetscPrintf(PETSC_COMM_WORLD,"reading model initial geometry from options\n");
    ierr = PetscOptionsGetInt(NULL,NULL,"-initial_region_geom_type",(PetscInt*)&data->initial_geom,&found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Type of region index initialisation","-initial_region_geom_type","1=pmap, 2=CAD");CHKERRQ(ierr);
	}
	
    /* box geometry */
    PetscPrintf(PETSC_COMM_WORLD, "reading box geometry from options\n");
	
    ierr = PetscOptionsGetReal(NULL,NULL, "-Lx", &data->Lx, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE) {
		ierr = ReportOptionMissing(0,"Model length in x","-Lx",0);CHKERRQ(ierr);
	}
	
    ierr = PetscOptionsGetReal(NULL,NULL, "-Ly", &data->Ly, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE) {
		ierr = ReportOptionMissing(0,"Model length in y","-Ly",0);CHKERRQ(ierr);
	}
	
    ierr = PetscOptionsGetReal(NULL,NULL, "-Lz", &data->Lz, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE) {
		ierr = ReportOptionMissing(0,"Model length in z","-Lz",0);CHKERRQ(ierr);
	}
	
	
    ierr = PetscOptionsGetReal(NULL,NULL, "-Ox", &data->Ox, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Model origin in x","-Ox",0);CHKERRQ(ierr);
	}
	
    ierr = PetscOptionsGetReal(NULL,NULL, "-Oy", &data->Oy, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Model origin in y","-Oy",0);CHKERRQ(ierr);
	}
	
    ierr = PetscOptionsGetReal(NULL,NULL, "-Oz", &data->Oz, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Model origin in z","-Oz",0);CHKERRQ(ierr);
	}
	
	
    /* bc type */
	opts_bcs = GENEBC_FreeSlip;
	ierr = PetscOptionsGetInt(NULL,NULL,"-bc_type",&opts_bcs,&found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(0,"Boundary condition type","-bc_type",0);CHKERRQ(ierr);
	}
    data->boundary_conditon_type = opts_bcs;
	
    
    /* material properties */
    PetscPrintf(PETSC_COMM_WORLD,"reading material properties from options\n");
    ierr = PetscOptionsGetInt(NULL,MODEL_NAME, "-nphase", &nregions, &found);CHKERRQ(ierr);
    if (found == PETSC_FALSE)	{
		ierr = ReportOptionMissing(MODEL_NAME,"Maximum number of material regions","-nphase",0);CHKERRQ(ierr);
	}
	
	/* set the active phases on the rheo structure */
	DataBucketSetSizes(c->material_constants,nregions,-1);
	
	ierr = PetscOptionsGetInt(NULL,MODEL_NAME,"-rheol",&rheol_type,&found);CHKERRQ(ierr);
	if (found == PETSC_FALSE){
		ierr = ReportOptionMissing(MODEL_NAME,"Rheology type","-rheol","0=viscous , 4=visco_plastic");CHKERRQ(ierr);
	}
	/* This needs to be stored ELSEWHERE */
	c->rheology_constants.rheology_type = rheol_type;
	
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
			
            //		case RHEOLOGY_VISCO_PLASTIC:
            //			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "-rheol RHEOLOGY_VISCO_PLASTIC");
            //			break;
			
            //		case RHEOLOGY_VISCO_PLASTIC_STRAIN_WEAKENING:
            //			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "-rheol RHEOLOGY_VISCO_PLASTIC_STRAIN_WEAKENING");
            //			break;
			
            //		case RHEOLOGY_VISCO_ELASTIC_PLASTIC:
            //			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "-rheol RHEOLOGY_VISCO_ELASTIC_PLASTIC");
            //			break;
			
	}
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelGene3DNueve_ApplyBoundaryCondition"
PetscErrorCode ModelGene3DNueve_ApplyBoundaryCondition(DM dav,BCList u_bclist,ModelGene3DNueveCtx *data)
{
    PetscScalar zero = 0.0;
	PetscReal ext_value;
	PetscBool flg;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    switch (data->boundary_conditon_type) {
			
        case GENEBC_FreeSlip:
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            break;
			
        case GENEBC_NoSlip:
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            break;
			
        case GENEBC_FreeSlipFreeSurface:
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            
			/* basement */
			ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            
			
			ierr = PetscOptionsGetReal(NULL,NULL,"-extension_x",&ext_value,&flg);CHKERRQ(ierr);
			if (flg) {
				ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &ext_value);CHKERRQ(ierr);
				ext_value = -1.0 * ext_value;
				ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &ext_value);CHKERRQ(ierr);
			}
			
			ierr = PetscOptionsGetReal(NULL,NULL,"-extension_z",&ext_value,&flg);CHKERRQ(ierr);
			if (flg) {
				ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMAX_LOC, 2,BCListEvaluator_constant, (void *) &ext_value);CHKERRQ(ierr);
				ext_value = -1.0 * ext_value;
				ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_KMIN_LOC, 2,BCListEvaluator_constant, (void *) &ext_value);CHKERRQ(ierr);
			}
			
			break;
			
        case GENEBC_NoSlipFreeSurface:
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_IMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMAX_LOC, 0,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMAX_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMAX_LOC, 2,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
            break;
			
		case GENEBC_ShearFreeSlipFreeSurface:
		{
			PetscReal exz_value,zero;
			
			exz_value = 0.25;
			ierr = PetscOptionsGetReal(NULL,NULL,"-shear_exz",&exz_value,&flg);CHKERRQ(ierr);
            /*
             // all sides //
             ierr = DirichletBC_ApplyStrainRateExz_b(u_bclist,dav,exz_value);CHKERRQ(ierr);
             // basement //
             zero = 0.0;
             ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
             */
			
			// east/west sides - y is unconstrained //
			ierr = DirichletBC_ApplyStrainRateExz_c(u_bclist,dav,exz_value);CHKERRQ(ierr);
            
			/* front/back sides */
			ierr = DirichletBC_FreeSlip(u_bclist,dav,FRONT_FACE);CHKERRQ(ierr);
			ierr = DirichletBC_FreeSlip(u_bclist,dav,BACK_FACE);CHKERRQ(ierr);
			
			// basement //
			zero = 0.0;
			ierr = DMDABCListTraverse3d(u_bclist,dav, DMDABCList_JMIN_LOC, 1,BCListEvaluator_constant, (void *) &zero);CHKERRQ(ierr);
			
		}
			break;
			
        default:
            break;
	}
	
    PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Gene3DNueve"
PetscErrorCode ModelApplyBoundaryCondition_Gene3DNueve(pTatinCtx user,void *ctx)
{
    ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
	BCList              u_bclist;
	DM                  dav;
    PetscErrorCode      ierr;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
    
	u_bclist = user->stokes_ctx->u_bclist;
	dav      = user->stokes_ctx->dav;
	
	ierr = ModelGene3DNueve_ApplyBoundaryCondition(dav,u_bclist,data);CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Gene3DNueve"
PetscErrorCode ModelApplyBoundaryConditionMG_Gene3DNueve(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
    ModelGene3DNueveCtx *data = (ModelGene3DNueveCtx*)ctx;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		ierr = ModelGene3DNueve_ApplyBoundaryCondition(dav[n],bclist[n],data);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Gene3DNueve"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Gene3DNueve(pTatinCtx c,void *ctx)
{
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
    PetscPrintf(PETSC_COMM_WORLD, "  NOT IMPLEMENTED \n", __FUNCT__);
	
    PetscFunctionReturn(0);
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
	
    PetscFunctionReturn(0);
}

/* Define eta,rho on material points */
#undef __FUNCT__
#define __FUNCT__ "Gene3DNueve_MaterialPointSetInitialStokesVariables"
PetscErrorCode Gene3DNueve_MaterialPointSetInitialStokesVariables(pTatinCtx c,void *ctx)
{
    int p, n_mp_points;
    DataBucket material_point_db;
    DataBucket material_constants_db;
    DataField PField_std, PField_stokes, PField_ViscConst;
	MaterialConst_ViscosityConst *material_visc_const;
    int region_id;
	
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
		MPntPStokesSetField_density(mpprop_stokes,       0.0 );
	}
	
    DataFieldRestoreAccess(PField_std);
    DataFieldRestoreAccess(PField_stokes);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointSetRegion_MyFunction"
PetscErrorCode MaterialPointSetRegion_MyFunction(pTatinCtx c,void *ctx)
{
    int p, n_mp_points;
    DataBucket db;
    DataField PField_std;
    int phase_init, phase;
    PetscReal arg,alpha, beta,gamma,delta;
	
	PetscFunctionBegin;
    PetscPrintf (PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    /* define properties on material points */
    db = c->materialpoint_db;
    DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
    DataFieldGetAccess(PField_std);
    DataFieldVerifyAccess(PField_std,sizeof (MPntStd));
	
    DataBucketGetSizes(db,&n_mp_points,0,0);
    
	alpha = 1.0;
	beta = 1.2;
	gamma = 1.0/6.0;
	delta = 0.0;
    
	PetscOptionsGetReal(NULL,NULL,"-int_alpha",&alpha,0);
	PetscOptionsGetReal(NULL,NULL,"-int_beta",&beta,0);
	PetscOptionsGetReal(NULL,NULL,"-int_gamma",&gamma,0);
	PetscOptionsGetReal(NULL,NULL,"-int_delta",&delta,0);
	
    for (p=0; p<n_mp_points; p++) {
		MPntStd *material_point;
		double *pos,xp,yp,zp,height;
		
		
		DataFieldAccessPoint(PField_std, p, (void **) &material_point);
		MPntStdGetField_global_coord(material_point,&pos);
		
		MPntStdGetField_phase_index(material_point, &phase_init);
        
		xp = pos[0] - 6.0;
		yp = pos[1];
		zp = pos[2] - 6.0;
		
		arg = ( beta * cos(gamma*M_PI*(xp + delta * zp)) - zp );
		height = 0.5 * exp( -alpha * arg * arg ) + 1.5;
        
		phase = 0;
		if (yp >= height ) {
			phase = 1;
		}
		
		MPntStdSetField_phase_index(material_point, phase);
	}
    DataFieldRestoreAccess(PField_std);
	
    PetscFunctionReturn(0);
}


/* define phase index on material points from a map file extruded in z direction */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSetRegionIndexFromMap"
PetscErrorCode MaterialPointSetRegionIndexFromMap(pTatinCtx c,void *ctx)
{
    PetscErrorCode ierr;
    CartGrid phasemap;
    PetscInt dir_0,dir_1,direction;
    int p, n_mp_points;
    DataBucket db;
    DataField PField_std;
    int phase_init, phase, phase_index;
    char map_file[PETSC_MAX_PATH_LEN], *name;
    PetscBool flg,phasefound;
    
	PetscFunctionBegin;
    PetscPrintf (PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    ierr = PetscOptionsGetString(NULL,MODEL_NAME,"-map_file",map_file,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (flg == PETSC_FALSE) {
		ierr = ReportOptionMissing(MODEL_NAME,"Name of pmap file","-map_file","HINT: Don't include the extension .pmap in the file name");CHKERRQ(ierr);
	}
    ierr = PetscOptionsGetInt(NULL,MODEL_NAME,"-extrude_dir",&direction,&flg);CHKERRQ(ierr);
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
      default :
        SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"-extrude_dir %d not valid",direction);
    }

    asprintf(&name,"./inputdata/%s.pmap",map_file);
    ierr = CartGridCreate(&phasemap);CHKERRQ(ierr);
    ierr = CartGridSetFilename(phasemap,map_file);CHKERRQ(ierr);
    ierr = CartGridSetUp(phasemap);CHKERRQ(ierr);
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
		
		ierr = CartGridGetValue(phasemap,position2D,&phase_index,&phasefound);CHKERRQ(ierr);
		
		if (phasefound == 1) {			/* point located in the phase map */
			phase = phase_index;
		} else {
			SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,"marker outside the domain\n your phasemap is smaller than the domain \n please check your parameters and retry");
		}
		MPntStdSetField_phase_index(material_point, phase);
		
	}
    ierr = CartGridDestroy(&phasemap);CHKERRQ(ierr);
    DataFieldRestoreAccess(PField_std);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointSetRegionIndexFromExtrudedMap"
PetscErrorCode MaterialPointSetRegionIndexFromExtrudedMap(DataBucket db,const char map_filename[],PetscInt direction,double c0,double c1)
{
    CartGrid  phasemap;
    PetscInt  dir_0,dir_1,dir_normal;
    int p, n_mp_points;
    DataField PField_std;
    int       phase_init,phase,phase_index;
    PetscBool phasefound;
    PetscErrorCode ierr;
    
	PetscFunctionBegin;
	
  switch (direction) {
    case 0:
      {
        dir_0 = 2;
        dir_1 = 1;
        dir_normal = 0;
      }
      break;

    case 1:
      {
        dir_0 = 0;
        dir_1 = 2;
        dir_normal = 1;
      }
      break;

    case 2:
      {
        dir_0 = 0;
        dir_1 = 1;
        dir_normal = 2;
      }
      break;
    default :
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"direction %d not valid",direction);
  }
	
    ierr = CartGridCreate(&phasemap);CHKERRQ(ierr);
    ierr = CartGridSetFilename(phasemap,map_filename);CHKERRQ(ierr);
    ierr = CartGridSetUp(phasemap);CHKERRQ(ierr);
	
    /* define properties on material points */
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
		
		if (pos[dir_normal] < c0) continue;
		if (pos[dir_normal] > c1) continue;
		
		
		MPntStdGetField_phase_index(material_point, &phase_init);
		
		ierr = CartGridGetValue(phasemap,position2D,&phase_index,&phasefound);CHKERRQ(ierr);
		
		if (phasefound) {			/* point located in the phase map */
			phase = phase_index;
		} else { /* Don't error if the point is outside the phase map - catch this error outside this function */
			phase = -1;
		}
		MPntStdSetField_phase_index(material_point, phase);
		
	}
    DataFieldRestoreAccess(PField_std);
    ierr = CartGridDestroy(&phasemap);CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelGene3DNueve_MaterialPointSetRegionIndexFromMultipleExtrudedMap"
PetscErrorCode ModelGene3DNueve_MaterialPointSetRegionIndexFromMultipleExtrudedMap(pTatinCtx c,void *ctx)
{
    PetscErrorCode ierr;
    DataBucket     db;
	PetscInt       n,nmaps,direction;
    char           opt_name[PETSC_MAX_PATH_LEN],map_file[PETSC_MAX_PATH_LEN],name[PETSC_MAX_PATH_LEN];
	PetscReal      range[2];
	PetscInt       nfetched;
    PetscBool      flg;
    
	PetscFunctionBegin;
    PetscPrintf (PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
    
    db = c->materialpoint_db;
	
	
    ierr = PetscOptionsGetInt(NULL,MODEL_NAME,"-num_map_files",&nmaps,&flg);CHKERRQ(ierr);
    if (flg == PETSC_FALSE) {
		ierr = ReportOptionMissing(MODEL_NAME,"Number of pmap files to read","-num_map_files",0);CHKERRQ(ierr);
	}
    
	direction = 0;
    ierr = PetscOptionsGetInt(NULL,MODEL_NAME,"-extrude_dir",&direction,&flg);CHKERRQ(ierr);
	if ( (flg == PETSC_FALSE) || (direction>2) || (direction<0)) {
		ierr = ReportOptionMissing(MODEL_NAME,"Direction to extrude the pmaps","-extrude_dir","0=x, 1=y, 2=z");CHKERRQ(ierr);
	}
	
	for (n=0; n<nmaps; n++) {
		/* Get file name */
		sprintf(opt_name,"-map_file_%d",n);
		ierr = PetscOptionsGetString(NULL,MODEL_NAME,opt_name,map_file,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
		if (flg == PETSC_FALSE) {
			ierr = ReportOptionMissing(MODEL_NAME,"Name of pmap file",opt_name,"HINT: Don't include the extension .pmap in the file name");CHKERRQ(ierr);
		}
        
		sprintf(name,"./inputdata/%s.pmap",map_file);
		
		/* Get extent of map */
		sprintf(opt_name,"-map_range_%d",n);
		nfetched = 2;
		ierr = PetscOptionsGetRealArray(NULL,MODEL_NAME,opt_name,range,&nfetched,&flg);CHKERRQ(ierr);
		if (!flg) {
			ierr = ReportOptionMissing(MODEL_NAME,"Range of pmap file in the normal direction",opt_name,"e.g. 0.0,1.0");CHKERRQ(ierr);
		}
		if (nfetched != 2) {
			ierr = ReportOptionMissing(MODEL_NAME,"Two values for pmap range in the normal direction",opt_name,"e.g. 0.0,1.0");CHKERRQ(ierr);
		}
		
		ierr = MaterialPointSetRegionIndexFromExtrudedMap(db,name,direction,range[0],range[1]);CHKERRQ(ierr);
	}
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointInitializeRegionIndex"
PetscErrorCode MaterialPointInitializeRegionIndex(DataBucket db,const int init_region_id)
{
    int                p,n_mp_points;
    DataField          PField_std;
	
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
    int                p,n_mp_points;
    DataField          PField_std;

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
            
        case GENE_LayeredCake: /* Layer cake <obsolete> */
		{
			ierr = MaterialPointSetRegion_MyFunction(c,ctx);CHKERRQ(ierr);
		}
            break;
            
        case GENE_ExtrudeFromMap: /* Extrude map along x,y,z directions */
		{
			ierr = MaterialPointSetRegionIndexFromMap(c,ctx);CHKERRQ(ierr);
		}
            break;
            
        case GENE_ReadFromCAD: /* Read from CAD file */
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Reading region index from CAD is not implemented yet\n");
		}
            break;
			
		case GENE_ExtrudeByPartsFromMap:
		{
			ierr = ModelGene3DNueve_MaterialPointSetRegionIndexFromMultipleExtrudedMap(c,ctx);CHKERRQ(ierr);
		}
	}
	/* Check all phase indices are between [0---rheo->max_phases-1] */
	DataBucketGetSizes(c->material_constants,&nregions,0,0);
	ierr = MaterialPointCheckRegionIndexBounds(c->materialpoint_db,0,nregions);CHKERRQ(ierr);
    
	/* Set eta0, rho0 */
    ierr = Gene3DNueve_MaterialPointSetInitialStokesVariables(c,ctx);CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Gene3DNueve"
PetscErrorCode ModelApplyUpdateMeshGeometry_Gene3DNueve(pTatinCtx c,Vec X,void *ctx)
{
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
    PetscPrintf(PETSC_COMM_WORLD, "  NOT IMPLEMENTED \n", __FUNCT__);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Gene3DNueve"
PetscErrorCode ModelOutput_Gene3DNueve(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
    PetscPrintf(PETSC_COMM_WORLD, "[[%s]]\n", __FUNCT__);
	
    ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
    ierr = pTatin3d_ModelOutput_MPntStd(c,prefix); CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
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
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Gene3DNueve"
PetscErrorCode pTatinModelRegister_Gene3DNueve(void)
{
    ModelGene3DNueveCtx *data;
    pTatinModel m;
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
	ierr =  pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_BCMG,            (void (*)(void)) ModelApplyBoundaryConditionMG_Gene3DNueve);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void)) ModelApplyMaterialBoundaryCondition_Gene3DNueve);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void)) ModelApplyInitialMeshGeometry_Gene3DNueve);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void)) ModelApplyInitialMaterialGeometry_Gene3DNueve);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void)) ModelApplyUpdateMeshGeometry_Gene3DNueve);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_OUTPUT,                (void (*)(void)) ModelOutput_Gene3DNueve);CHKERRQ(ierr);
    ierr =	pTatinModelSetFunctionPointer(m, PTATIN_MODEL_DESTROY,               (void (*)(void)) ModelDestroy_Gene3DNueve); CHKERRQ(ierr);
	
    /* Insert model into list */
    ierr = pTatinModelRegister(m); CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}

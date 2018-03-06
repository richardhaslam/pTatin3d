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
 **    filename:   phys_comp_energy.c
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
 
 Physics componenet: Energy (advection-diffusion equation)
 
 Require functions:
 
 PhysCompCreateMesh_Energy()
 PhysCompCreateBoundaryList_Energy()
 PhysCompCreateVolumeQuadrature_Energy()
 
 PhysCompCreate_Energy() : 
 PhysCompDestroy_Energy() : 
 
 PhysCompsNew_Energy() : Allocates new structure
 PhysCompLoad_Energy() : Loads from checkpoint file 
 PhysCompSave_Energy() : Saves to checkpoint file

 pTatin3d_PhysCompCreate_Energy() : Reads checkpoint data or calls PhysCompsNew_Energy()
 
*/


#include "petsc.h"
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "data_bucket.h"
#include "dmda_bcs.h"
#include "element_utils_q1.h"
#include "dmdae.h"
#include "dmda_element_q1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"
#include "QPntVolCoefEnergy_def.h"
#include "phys_comp_energy.h"
#include "dmda_checkpoint.h"
#include "cJSON.h"
#include "cjson_utils.h"


#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreate_Energy"
PetscErrorCode PhysCompCreate_Energy(PhysCompEnergy *E)
{
	PetscErrorCode ierr;
	PhysCompEnergy energy;
	
	PetscFunctionBegin;
  *E = NULL;
	ierr = PetscMalloc(sizeof(struct _p_PhysCompEnergy),&energy);CHKERRQ(ierr);
	ierr = PetscMemzero(energy,sizeof(struct _p_PhysCompEnergy));CHKERRQ(ierr);
	*E = energy;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompDestroy_Energy"
PetscErrorCode PhysCompDestroy_Energy(PhysCompEnergy *E)
{
	PetscErrorCode ierr;
	PhysCompEnergy ctx;
	
	PetscFunctionBegin;
	
	if (!E) {PetscFunctionReturn(0);}
	ctx = *E;
	
	//	for (e=0; e<HEX_FACES; e++) {
	//		if (ctx->surfQ[e]) { ierr = SurfaceQuadratureDestroy(&ctx->surfQ[e]);CHKERRQ(ierr); }
	//	}
	if (ctx->volQ) { ierr = QuadratureDestroy(&ctx->volQ);CHKERRQ(ierr); }
	if (ctx->T_bclist) { ierr = BCListDestroy(&ctx->T_bclist);CHKERRQ(ierr); }
  if (ctx->daT) {
    ierr = DMDestroyDMDAE(ctx->daT);CHKERRQ(ierr);
    ierr = DMDestroy(&ctx->daT);CHKERRQ(ierr);
  }
  if (ctx->u_minus_V) { ierr = VecDestroy(&ctx->u_minus_V);CHKERRQ(ierr); }
  if (ctx->Told) {      ierr = VecDestroy(&ctx->Told);CHKERRQ(ierr); }
  if (ctx->Xold) {      ierr = VecDestroy(&ctx->Xold);CHKERRQ(ierr); }
	if (ctx) { ierr = PetscFree(ctx);CHKERRQ(ierr); }
	
	*E = NULL;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateMesh_Energy"
PetscErrorCode PhysCompCreateMesh_Energy(PhysCompEnergy E,DM dav,PetscInt mx,PetscInt my, PetscInt mz,PetscInt mesh_generator_type)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	E->energy_mesh_type = mesh_generator_type;

	switch (mesh_generator_type) {
		DMDAE dae;
			
		case 0:
			PetscPrintf(PETSC_COMM_WORLD,"PhysCompCreateMesh_Energy: Generating standard Q1 DMDA\n");
			E->mx = mx;
			E->my = my;
			E->mz = mz;
			
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only overlapping and nested supported {1,2} ");
			break;
		
		case 1:
			PetscPrintf(PETSC_COMM_WORLD,"PhysCompCreateMesh_Energy: Generating overlapping Q1 DMDA\n");
			if (!dav) {
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Require valid DM dav");
			}
			ierr = DMDACreateOverlappingQ1FromQ2(dav,1,&E->daT);CHKERRQ(ierr);
			ierr = DMGetDMDAE(E->daT,&dae);CHKERRQ(ierr);
			E->mx = dae->mx;
			E->my = dae->my;
			E->mz = dae->mz;
			break;
		
		case 2:
			PetscPrintf(PETSC_COMM_WORLD,"PhysCompCreateMesh_Energy: Generating nested Q1 DMDA\n");
			if (!dav) {
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Require valid DM dav");
			}
			ierr = DMDACreateNestedQ1FromQ2(dav,1,&E->daT);CHKERRQ(ierr);
			ierr = DMGetDMDAE(E->daT,&dae);CHKERRQ(ierr);
			E->mx = dae->mx;
			E->my = dae->my;
			E->mz = dae->mz;
			break;

		default:
			E->energy_mesh_type = 1;

			PetscPrintf(PETSC_COMM_WORLD,"PhysCompCreateMesh_Energy: Generating overlapping Q1 DMDA\n");
			if (!dav) {
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Require valid DM dav");
			}
			ierr = DMDACreateOverlappingQ1FromQ2(dav,1,&E->daT);CHKERRQ(ierr);
			ierr = DMGetDMDAE(E->daT,&dae);CHKERRQ(ierr);
			E->mx = dae->mx;
			E->my = dae->my;
			E->mz = dae->mz;
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateBoundaryList_Energy"
PetscErrorCode PhysCompCreateBoundaryList_Energy(PhysCompEnergy E)
{
	DM daT;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	daT = E->daT;
	if (!daT) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"daT must be set"); }
	ierr = DMDABCListCreate(daT,&E->T_bclist);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateVolumeQuadrature_Energy"
PetscErrorCode PhysCompCreateVolumeQuadrature_Energy(PhysCompEnergy E)
{
	PetscInt dim, np_per_dim, ncells;
	DMDAE dae;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	np_per_dim = 2;
	dim = 3;

	ierr = DMGetDMDAE(E->daT,&dae);CHKERRQ(ierr);
	ncells = dae->lmx * dae->lmy * dae->lmz;
	
	ierr = VolumeQuadratureCreate_GaussLegendreEnergy(dim,np_per_dim,ncells,&E->volQ);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompNew_Energy"
PetscErrorCode PhysCompNew_Energy(DM dav,PetscInt mx,PetscInt my, PetscInt mz,PetscInt mesh_generator_type,PhysCompEnergy *E)
{
	PetscErrorCode  ierr;
	PhysCompEnergy  energy;
	DM              cda;
	
	PetscFunctionBegin;

	ierr = PhysCompCreate_Energy(&energy);CHKERRQ(ierr);

	ierr = PhysCompCreateMesh_Energy(energy,dav,mx,my,mz,mesh_generator_type);CHKERRQ(ierr);
	ierr = PhysCompCreateBoundaryList_Energy(energy);CHKERRQ(ierr);
	ierr = PhysCompCreateVolumeQuadrature_Energy(energy);CHKERRQ(ierr);
	
	/* create aux vectors */
	ierr = DMCreateGlobalVector(energy->daT,&energy->Told);CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(energy->daT,&cda);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(cda,&energy->u_minus_V);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(cda,&energy->Xold);CHKERRQ(ierr);
	
	*E = energy;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompLoad_Energy"
PetscErrorCode PhysCompLoad_Energy(void)
{
	PetscFunctionBegin;
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Currently unavailable");
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompSave_Energy"
PetscErrorCode PhysCompSave_Energy(void)
{
	PetscFunctionBegin;
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Currently unavailable");
	PetscFunctionReturn(0);
}

/* quadrature */
#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureCreate_GaussLegendreEnergy"
PetscErrorCode VolumeQuadratureCreate_GaussLegendreEnergy(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature)
{
	Quadrature Q;
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
	
	ierr = QuadratureCreate(&Q);CHKERRQ(ierr);
	Q->dim  = nsd;
	Q->type = VOLUME_QUAD;
	
	/*PetscPrintf(PETSC_COMM_WORLD,"VolumeQuadratureCreate_GaussLegendreEnergy:\n");*/
	switch (np_per_dim) {
		case 1:
			/*PetscPrintf(PETSC_COMM_WORLD,"\tUsing 1 pnt Gauss Legendre quadrature\n");*/
			//QuadratureCreateGauss_1pnt_3D(&ngp,gp_xi,gp_weight);
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"This will result in a rank-deficient operator");
			break;
			
		case 2:
			/*PetscPrintf(PETSC_COMM_WORLD,"\tUsing 2x2 pnt Gauss Legendre quadrature\n");*/
			QuadratureCreateGauss_2pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
			
		case 3:
			/*PetscPrintf(PETSC_COMM_WORLD,"\tUsing 3x3 pnt Gauss Legendre quadrature\n");*/
			QuadratureCreateGauss_3pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
			
		default:
			/*PetscPrintf(PETSC_COMM_WORLD,"\tUsing 3x3 pnt Gauss Legendre quadrature\n");*/
			QuadratureCreateGauss_3pnt_3D(&Q->npoints,&Q->q_xi_coor,&Q->q_weight);
			break;
	}
	
	Q->n_elements = ncells;
	if (ncells != 0) {
		
		DataBucketCreate(&Q->properties_db);
		DataBucketRegisterField(Q->properties_db,QPntVolCoefEnergy_classname, sizeof(QPntVolCoefEnergy),NULL);
		DataBucketFinalize(Q->properties_db);
		
		DataBucketSetInitialSizes(Q->properties_db,Q->npoints*ncells,1);
		
    /*
    // Note: This call will hang if any rank contained zero elements
		DataBucketView(PETSC_COMM_WORLD, Q->properties_db,"GaussLegendre EnergyCoefficients",DATABUCKET_VIEW_STDOUT);
    */
	}
	
	*quadrature = Q;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureGetAllCellData_Energy"
PetscErrorCode VolumeQuadratureGetAllCellData_Energy(Quadrature Q,QPntVolCoefEnergy *coeffs[])
{
	QPntVolCoefEnergy *quadraturepoint_data;
  DataField          PField;
	PetscFunctionBegin;
	
	DataBucketGetDataFieldByName(Q->properties_db, QPntVolCoefEnergy_classname ,&PField);
	quadraturepoint_data = PField->data;
	*coeffs = quadraturepoint_data;
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VolumeQuadratureGetCellData_Energy"
PetscErrorCode VolumeQuadratureGetCellData_Energy(Quadrature Q,QPntVolCoefEnergy coeffs[],PetscInt cidx,QPntVolCoefEnergy *cell[])
{
  PetscFunctionBegin;
	if (cidx>=Q->n_elements) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_SIZ,"cidx > max cells");
	}
	
	*cell = &coeffs[cidx*Q->npoints];
  PetscFunctionReturn(0);
}

/* material points */
#undef __FUNCT__  
#define __FUNCT__ "PhysCompAddMaterialPointCoefficients_Energy"
PetscErrorCode PhysCompAddMaterialPointCoefficients_Energy(DataBucket db)
{
	PetscFunctionBegin;
	/* register marker structures here */
	DataBucketRegisterField(db,MPntPEnergy_classname,sizeof(MPntPEnergy),NULL);

	PetscFunctionReturn(0);
}

/* unneeded */
#undef __FUNCT__
#define __FUNCT__ "PhysCompCheckpointLoad_Energy"
PetscErrorCode PhysCompCheckpointLoad_Energy(const char prefix[],PhysCompEnergy e)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PhysCompCheckpointWrite_Energy"
PetscErrorCode PhysCompCheckpointWrite_Energy(PhysCompEnergy e,PetscBool write_dmda,const char path[],const char prefix[])
{
  char jfilename[PETSC_MAX_PATH_LEN],vfilename[3][PETSC_MAX_PATH_LEN],daprefix[PETSC_MAX_PATH_LEN];
  MPI_Comm       comm;
  PetscMPIInt    commsize,commrank;
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)e->daT,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  
  if (path && prefix) SETERRQ(comm,PETSC_ERR_SUP,"Only support for path/file or prefix_file");
  
  daprefix[0] = '\0';
  if (write_dmda) {
    if (path) {
      PetscSNPrintf(daprefix,PETSC_MAX_PATH_LEN-1,"%s/physcomp_energy",path);
    } else {
      PetscSNPrintf(daprefix,PETSC_MAX_PATH_LEN-1,"%s_physcomp_energy",prefix);
    }
    ierr = DMDACheckpointWrite(e->daT,daprefix);CHKERRQ(ierr);
  }
  
  if (path) {
    PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s/physcomp_energy.json",path);
    PetscSNPrintf(vfilename[0],PETSC_MAX_PATH_LEN-1,"%s/physcomp_energy_u_minus_V.pbvec",path);
    PetscSNPrintf(vfilename[1],PETSC_MAX_PATH_LEN-1,"%s/physcomp_energy_Told.pbvec",path);
    PetscSNPrintf(vfilename[2],PETSC_MAX_PATH_LEN-1,"%s/physcomp_energy_Xold.pbvec",path);
  } else {
    PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s_physcomp_energy.json",prefix);
    PetscSNPrintf(vfilename[0],PETSC_MAX_PATH_LEN-1,"%s_physcomp_energy_u_minus_V.pbvec",prefix);
    PetscSNPrintf(vfilename[1],PETSC_MAX_PATH_LEN-1,"%s_physcomp_energy_Told.pbvec",prefix);
    PetscSNPrintf(vfilename[2],PETSC_MAX_PATH_LEN-1,"%s_physcomp_energy_Xold.pbvec",prefix);
  }
  
  if (e->u_minus_V) {
    ierr = DMDAWriteVectorToFile(e->u_minus_V,vfilename[0],PETSC_FALSE);CHKERRQ(ierr);
  }
  if (e->Told) {
    ierr = DMDAWriteVectorToFile(e->Told,vfilename[1],PETSC_FALSE);CHKERRQ(ierr);
  }
  if (e->Xold) {
    ierr = DMDAWriteVectorToFile(e->Xold,vfilename[2],PETSC_FALSE);CHKERRQ(ierr);
  }
  
  if (commrank == 0) {
    cJSON *jso_file = NULL,*jso_energy = NULL,*content;
    
    /* create json meta data file */
    jso_file = cJSON_CreateObject();
    
    jso_energy = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_file,"PhysCompEnergy",jso_energy);
    
    content = cJSON_CreateInt((int)e->mx);         cJSON_AddItemToObject(jso_energy,"mx",content);
    content = cJSON_CreateInt((int)e->my);         cJSON_AddItemToObject(jso_energy,"my",content);
    content = cJSON_CreateInt((int)e->mz);         cJSON_AddItemToObject(jso_energy,"mz",content);

    content = cJSON_CreateInt((int)e->energy_mesh_type);   cJSON_AddItemToObject(jso_energy,"meshType",content);
    content = cJSON_CreateNumber((double)e->time);         cJSON_AddItemToObject(jso_energy,"time",content);
    content = cJSON_CreateNumber((double)e->dt);           cJSON_AddItemToObject(jso_energy,"timeStepSize",content);

    if (write_dmda) {
      cJSON *jso_dmda;
      char subdmfilename[PETSC_MAX_PATH_LEN];
      
      jso_dmda = cJSON_CreateObject();
      cJSON_AddItemToObject(jso_energy,"sub_dmda",jso_dmda);
      
      PetscSNPrintf(subdmfilename,PETSC_MAX_PATH_LEN-1,"%s_dmda.json",daprefix);
      content = cJSON_CreateString(subdmfilename);       cJSON_AddItemToObject(jso_dmda,"fileName",content);
      content = cJSON_CreateString("json-meta");         cJSON_AddItemToObject(jso_dmda,"dataFormat",content);
    }

    {
      cJSON *jso_petscvec;
      
      jso_petscvec = cJSON_CreateObject();
      cJSON_AddItemToObject(jso_energy,"u_minus_V",jso_petscvec);
      
      content = cJSON_CreateString(vfilename[0]);   cJSON_AddItemToObject(jso_petscvec,"fileName",content);
      content = cJSON_CreateString("petsc-binary"); cJSON_AddItemToObject(jso_petscvec,"dataFormat",content);
    }
    {
      cJSON *jso_petscvec;
      
      jso_petscvec = cJSON_CreateObject();
      cJSON_AddItemToObject(jso_energy,"Told",jso_petscvec);
      
      content = cJSON_CreateString(vfilename[1]);   cJSON_AddItemToObject(jso_petscvec,"fileName",content);
      content = cJSON_CreateString("petsc-binary"); cJSON_AddItemToObject(jso_petscvec,"dataFormat",content);
    }
    {
      cJSON *jso_petscvec;
      
      jso_petscvec = cJSON_CreateObject();
      cJSON_AddItemToObject(jso_energy,"Xold",jso_petscvec);
      
      content = cJSON_CreateString(vfilename[2]);   cJSON_AddItemToObject(jso_petscvec,"fileName",content);
      content = cJSON_CreateString("petsc-binary"); cJSON_AddItemToObject(jso_petscvec,"dataFormat",content);
    }
    
    /* write json meta data file */
    {
      FILE *fp;
      char *jbuff = cJSON_Print(jso_file);
      
      fp = fopen(jfilename,"w");
      if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to open file %s",jfilename);
      fprintf(fp,"%s\n",jbuff);
      fclose(fp);
      free(jbuff);
    }
    
    cJSON_Delete(jso_file);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PhysCompLoad2_Energy"
PetscErrorCode PhysCompLoad2_Energy(DM dav,const char jfilename[],PhysCompEnergy *E)
{
  PetscErrorCode ierr;
  PetscMPIInt rank;
  char pathtovec[PETSC_MAX_PATH_LEN];
  cJSON *jfile = NULL,*jphys = NULL;
  PetscInt mx,my,mz,mesh_generator_type;
  PetscReal time,dt;
  MPI_Comm comm;
  PetscBool found;
  
  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)dav,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  if (rank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open JSON file \"%s\"",jfilename);
    jphys = cJSON_GetObjectItem(jfile,"PhysCompEnergy");
  }

  /* query JSON file for input parameters */
  ierr = cJSONGetPetscInt(comm,jphys,"mx",&mx,&found);CHKERRQ(ierr);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"mx\"");
  
  ierr = cJSONGetPetscInt(comm,jphys,"my",&my,&found);CHKERRQ(ierr);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"my\"");
  
  ierr = cJSONGetPetscInt(comm,jphys,"mz",&mz,&found);CHKERRQ(ierr);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"mz\"");
  
  ierr = cJSONGetPetscInt(comm,jphys,"meshType",&mesh_generator_type,&found);CHKERRQ(ierr);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"meshType\"");

  ierr = cJSONGetPetscReal(comm,jphys,"time",&time,&found);CHKERRQ(ierr);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"time\"");
  
  ierr = cJSONGetPetscReal(comm,jphys,"timeStepSize",&dt,&found);CHKERRQ(ierr);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"timeStepSize\"");

  /* 
   This function creates the DMDA for temperature.
   I elect to perform a self creation, rather than creating via DMDACheckpointLoad() as
   this particular DMDA possess extra content associated with finite elements in the
   form of an attached struct.
  */
  ierr = PhysCompNew_Energy(dav,mx,my,mz,mesh_generator_type,E);CHKERRQ(ierr);
  
  /* reset parameters in energy struct */
  (*E)->time = time;
  (*E)->dt = dt;
  
  /* query file for state vector filenames - load vectors into energy struct */
  {
    cJSON *jso_petscvec = NULL;
    
    /* u - V */
    if (jphys) { jso_petscvec = cJSON_GetObjectItem(jphys,"u_minus_V"); }
    if (!jso_petscvec) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to locate key \"u_minus_V\"");
    
    ierr = cJSONGetPetscString(comm,jso_petscvec,"fileName",pathtovec,&found);CHKERRQ(ierr);
    if (found) { ierr = VecLoadFromFile((*E)->u_minus_V,pathtovec);CHKERRQ(ierr);
    } else       SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"fileName\"");

    /* T^{k} */
    if (jphys) { jso_petscvec = cJSON_GetObjectItem(jphys,"Told"); }
    if (!jso_petscvec) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to locate key \"Told\"");
    
    ierr = cJSONGetPetscString(comm,jso_petscvec,"fileName",pathtovec,&found);CHKERRQ(ierr);
    if (found) { ierr = VecLoadFromFile((*E)->Told,pathtovec);CHKERRQ(ierr);
    } else       SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"fileName\"");

    /* X^{k} */
    if (jphys) { jso_petscvec = cJSON_GetObjectItem(jphys,"Xold"); }
    if (!jso_petscvec) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to locate key \"Xold\"");
    
    ierr = cJSONGetPetscString(comm,jso_petscvec,"fileName",pathtovec,&found);CHKERRQ(ierr);
    if (found) { ierr = VecLoadFromFile((*E)->Xold,pathtovec);CHKERRQ(ierr);
    } else       SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"fileName\"");
  }
  
  if (rank == 0) {
    cJSON_Delete(jfile);
  }
  
  PetscFunctionReturn(0);
}

/* unneeded */
#undef __FUNCT__
#define __FUNCT__ "PhysCompSave2_Energy"
PetscErrorCode PhysCompSave2_Energy(PhysCompEnergy e)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

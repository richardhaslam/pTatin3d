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
 **    filename:   dmda_checkpoint.c
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


#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>
#include "dmda_update_coords.h"
#include "dmda_checkpoint.h"
#include "cjson_utils.h"

PetscErrorCode DMDALoadGlobalVectorFromFile(DM da,const char name[],Vec *da_x)
{
	PetscErrorCode ierr;
	PetscViewer    v;
	MPI_Comm       comm;
	Vec            xn;
	

	PetscFunctionBegin;
	PetscObjectGetComm( (PetscObject)da, &comm );
	if (da == NULL) SETERRQ(comm,PETSC_ERR_USER, "da is NULL");
	
	ierr = DMCreateGlobalVector(da, &xn);CHKERRQ(ierr);

	ierr = PetscViewerCreate(comm,&v);CHKERRQ(ierr);
	ierr = PetscViewerSetType(v,PETSCVIEWERBINARY);CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(v,FILE_MODE_READ);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO
	ierr = PetscViewerBinarySetMPIIO(v);CHKERRQ(ierr);
#endif
	ierr = PetscViewerFileSetName(v,name);CHKERRQ(ierr);
	ierr = VecLoad(xn,v); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
	
	/* 
	 putain - VecLoadIntoVector inserts the option below into the command line.
	 This will screw shit up if you load in vectors with different block sizes.
	 */
	ierr = PetscOptionsClearValue(NULL,"-vecload_block_size");CHKERRQ(ierr);
	
	*da_x = xn;
	
	PetscFunctionReturn(0);	
}

PetscErrorCode DMDALoadCoordinatesFromFile(DM da,const char name[])
{
	PetscErrorCode ierr;
	DM             cda;
	Vec            coords,da_coords;
	
	
	PetscFunctionBegin;
	if (da == NULL) SETERRQ( PetscObjectComm((PetscObject)da),PETSC_ERR_USER, "da is NULL" );
	
	/* make sure the vector is present */
	ierr = DMDASetUniformCoordinates(da, 0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMDALoadGlobalVectorFromFile(cda,name,&coords);CHKERRQ(ierr);
	
	/* set the global coordinates */
	ierr = DMGetCoordinates(da,&da_coords);CHKERRQ(ierr);
	ierr = VecCopy(coords,da_coords);CHKERRQ(ierr);
	
	/* make sure the local coordinates are upto date */
	ierr = DMDAUpdateGhostedCoordinates(da);CHKERRQ(ierr);
	
	ierr = VecDestroy(&coords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);		
}

PetscErrorCode DMDAWriteVectorToFile(Vec x,const char name[],PetscBool zip_file)
{
	char fieldname[PETSC_MAX_PATH_LEN];
	PetscViewer viewer;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	if (zip_file) {
		sprintf(fieldname,"%s.gz",name);
	} else {
		sprintf(fieldname,"%s",name);
	}
	
  ierr = PetscViewerCreate(PetscObjectComm((PetscObject)x),&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO	
	ierr = PetscViewerBinarySetMPIIO(viewer);CHKERRQ(ierr);
#endif
	ierr = PetscViewerFileSetName(viewer,fieldname);CHKERRQ(ierr);
	
	ierr = VecView(x,viewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode VecLoadFromFile(Vec x,const char name[])
{
	PetscErrorCode ierr;
	PetscViewer    v;
	
	
	PetscFunctionBegin;
	
	ierr = PetscViewerCreate(PetscObjectComm((PetscObject)x),&v);CHKERRQ(ierr);
	ierr = PetscViewerSetType(v,PETSCVIEWERBINARY);CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(v,FILE_MODE_READ);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO
	ierr = PetscViewerBinarySetMPIIO(v);CHKERRQ(ierr);
#endif
	ierr = PetscViewerFileSetName(v,name);CHKERRQ(ierr);
	ierr = VecLoad(x,v); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&v); CHKERRQ(ierr);
	
	/* 
	 putain - VecLoadIntoVector inserts the option below into the command line.
	 This will screw shit up if you load in vectors with different block sizes.
	 */
	ierr = PetscOptionsClearValue(NULL,"-vecload_block_size");CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

PetscErrorCode DMDACheckpointWrite(DM da,const char jprefix[])
{
  char jfilename[PETSC_MAX_PATH_LEN],cfilename[PETSC_MAX_PATH_LEN];
  MPI_Comm       comm;
  PetscMPIInt    commsize,commrank;
  PetscInt       i,dim,M,N,P,m,n,p,ndof,sw;
  PetscInt       refine[3];
  DMBoundaryType bc_type[3];
  DMDAStencilType  st;
  const char *prefix;
  const PetscInt *_lri,*_lrj,*_lrk;
  int *lri,*lrj,*lrk;
  Vec coords = NULL;
  PetscErrorCode ierr;
  
  
  ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  
  ierr = DMGetOptionsPrefix(da,&prefix);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,&dim,&M,&N,&P,&m,&n,&p,&ndof,&sw,&bc_type[0],&bc_type[1],&bc_type[2],&st);CHKERRQ(ierr);
  ierr = DMDAGetRefinementFactor(da,&refine[0],&refine[1],&refine[2]);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRanges(da,&_lri,&_lrj,&_lrk);CHKERRQ(ierr);
  
  ierr = PetscMalloc1(m,&lri);CHKERRQ(ierr);
  ierr = PetscMalloc1(n,&lrj);CHKERRQ(ierr);
  ierr = PetscMalloc1(p,&lrk);CHKERRQ(ierr);
  for (i=0; i<m; i++) { lri[i] = (int)_lri[i]; }
  for (i=0; i<n; i++) { lrj[i] = (int)_lrj[i]; }
  for (i=0; i<p; i++) { lrk[i] = (int)_lrk[i]; }

  if (jprefix) {
    PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s_dmda.json",jprefix);
    PetscSNPrintf(cfilename,PETSC_MAX_PATH_LEN-1,"%s_dmda_coords.pbvec",jprefix);
  } else {
    PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"dmda.json");
    PetscSNPrintf(cfilename,PETSC_MAX_PATH_LEN-1,"dmda_coords.pbvec");
  }
  
  ierr = DMGetCoordinates(da,&coords);CHKERRQ(ierr);
  if (coords) {
    PetscViewer v;
    
    ierr = PetscViewerCreate(PetscObjectComm((PetscObject)da),&v);CHKERRQ(ierr);
    ierr = PetscViewerSetType(v,PETSCVIEWERBINARY);CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(v,FILE_MODE_WRITE);CHKERRQ(ierr);
#ifdef PTATIN_USE_MPIIO
    ierr = PetscViewerBinarySetMPIIO(v);CHKERRQ(ierr);
#endif
    ierr = PetscViewerFileSetName(v,cfilename);CHKERRQ(ierr);
    ierr = VecView(coords,v);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&v);CHKERRQ(ierr);
  }
  
  if (commrank == 0) {
    cJSON *jso_file = NULL,*jso_dm = NULL,*ja,*jso_part,*content,*obj;
    
    /* create json meta data file */
    jso_file = cJSON_CreateObject();
    
    jso_dm = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_file,"DMDA",jso_dm);
    
    if (prefix) {
      content = cJSON_CreateString((char*)prefix);  cJSON_AddItemToObject(jso_dm,"optionsPrefix",content);
    }
    content = cJSON_CreateInt((int)dim);          cJSON_AddItemToObject(jso_dm,"dim",content);
    content = cJSON_CreateInt((int)ndof);         cJSON_AddItemToObject(jso_dm,"ndof",content);
    content = cJSON_CreateInt((int)sw);           cJSON_AddItemToObject(jso_dm,"stencilWidth",content);
    switch (st) {
      case DMDA_STENCIL_STAR:
        content = cJSON_CreateString("star");     cJSON_AddItemToObject(jso_dm,"stencilType",content);
        break;
      case DMDA_STENCIL_BOX:
        content = cJSON_CreateString("box");      cJSON_AddItemToObject(jso_dm,"stencilType",content);
        break;
      default:
        SETERRQ(comm,PETSC_ERR_USER,"Unknown stenctil type detected");
        break;
    }			
    
    ja = cJSON_CreateArray();
    cJSON_AddItemToObject(jso_dm,"directions",ja);

    if (dim >= 1) {
      obj = cJSON_CreateObject();
      content = cJSON_CreateInt((int)M);            cJSON_AddItemToObject(obj,"points",content);
      content = cJSON_CreateInt((int)refine[0]);     cJSON_AddItemToObject(obj,"refinementFactor",content);
      switch (bc_type[0]) {
        case DM_BOUNDARY_NONE:
          content = cJSON_CreateString("none");     cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_GHOSTED:
          content = cJSON_CreateString("ghosted");  cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_MIRROR:
          content = cJSON_CreateString("mirror");   cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_PERIODIC:
          content = cJSON_CreateString("periodic"); cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        default:
          SETERRQ(comm,PETSC_ERR_USER,"Unknown bounaryType type detected");
          break;
      }
      cJSON_AddItemToArray(ja,obj);
    }
    if (dim >= 2) {
      obj = cJSON_CreateObject();
      content = cJSON_CreateInt((int)N);            cJSON_AddItemToObject(obj,"points",content);
      content = cJSON_CreateInt((int)refine[1]);     cJSON_AddItemToObject(obj,"refinementFactor",content);
      switch (bc_type[1]) {
        case DM_BOUNDARY_NONE:
          content = cJSON_CreateString("none");     cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_GHOSTED:
          content = cJSON_CreateString("ghosted");  cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_MIRROR:
          content = cJSON_CreateString("mirror");   cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_PERIODIC:
          content = cJSON_CreateString("periodic"); cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        default:
          SETERRQ(comm,PETSC_ERR_USER,"Unknown bounaryType type detected");
          break;
      }
      cJSON_AddItemToArray(ja,obj);
    }
    if (dim == 3) {
      obj = cJSON_CreateObject();
      content = cJSON_CreateInt((int)P);            cJSON_AddItemToObject(obj,"points",content);
      content = cJSON_CreateInt((int)refine[2]);     cJSON_AddItemToObject(obj,"refinementFactor",content);
      switch (bc_type[2]) {
        case DM_BOUNDARY_NONE:
          content = cJSON_CreateString("none");     cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_GHOSTED:
          content = cJSON_CreateString("ghosted");  cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_MIRROR:
          content = cJSON_CreateString("mirror");   cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        case DM_BOUNDARY_PERIODIC:
          content = cJSON_CreateString("periodic"); cJSON_AddItemToObject(obj,"boundaryType",content);
          break;
        default:
          SETERRQ(comm,PETSC_ERR_USER,"Unknown bounaryType type detected");
          break;
      }
      cJSON_AddItemToArray(ja,obj);
    }

    {
      cJSON *jso_c;
      
      jso_c = cJSON_CreateObject();
      cJSON_AddItemToObject(jso_dm,"coordinate",jso_c);

      if (coords) {
        content = cJSON_CreateString(cfilename); cJSON_AddItemToObject(jso_c,"dmCoordinateFile",content);
        content = cJSON_CreateString("petsc-binary"); cJSON_AddItemToObject(jso_c,"dataFormat",content);
      } else {
        content = cJSON_CreateString("null"); cJSON_AddItemToObject(jso_c,"dmCoordinateFile",content);
      }
    }
    
    jso_part = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_dm,"partition",jso_part);
    content = cJSON_CreateInt((int)commsize);  cJSON_AddItemToObject(jso_part,"commSize",content);
    
    ja = cJSON_CreateArray();
    cJSON_AddItemToObject(jso_part,"directions",ja);

    if (dim >= 1) {
      obj = cJSON_CreateObject();
      content = cJSON_CreateInt((int)m);                      cJSON_AddItemToObject(obj,"ranks",content);
      content = cJSON_CreateIntArray((const int*)lri,(int)m); cJSON_AddItemToObject(obj,"pointsPerRank",content);
      cJSON_AddItemToArray(ja,obj);
    }
    if (dim >= 2) {
      obj = cJSON_CreateObject();
      content = cJSON_CreateInt((int)n);                      cJSON_AddItemToObject(obj,"ranks",content);
      content = cJSON_CreateIntArray((const int*)lrj,(int)n); cJSON_AddItemToObject(obj,"pointsPerRank",content);
      cJSON_AddItemToArray(ja,obj);
    }
    if (dim == 3) {
      obj = cJSON_CreateObject();
      content = cJSON_CreateInt((int)p);                      cJSON_AddItemToObject(obj,"ranks",content);
      content = cJSON_CreateIntArray((const int*)lrk,(int)p); cJSON_AddItemToObject(obj,"pointsPerRank",content);
      cJSON_AddItemToArray(ja,obj);
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
  
  ierr = PetscFree(lri);CHKERRQ(ierr);
  ierr = PetscFree(lrj);CHKERRQ(ierr);
  ierr = PetscFree(lrk);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode _DMDAParseLayout_JSON(cJSON *jdmda,MPI_Comm comm,
                      char prefix[],
                      PetscInt *dim,
                      DMBoundaryType bc_type[],
                      DMDAStencilType *stencil_type,
                      PetscInt M[],
                      PetscInt *ndof,
                      PetscInt *sw)
{
  PetscErrorCode ierr;
  cJSON *dirlist = NULL,*d_k = NULL;
  PetscBool found,match;
  int k,nf;
  char stencilTypeName[PETSC_MAX_PATH_LEN];
  char boundaryTypeName[PETSC_MAX_PATH_LEN];
  
  ierr = cJSONGetPetscString(comm,jdmda,"optionsPrefix",prefix,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscInt(comm,jdmda,"dim",dim,&found);CHKERRQ(ierr); if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"dim\"");
  ierr = cJSONGetPetscInt(comm,jdmda,"ndof",ndof,&found);CHKERRQ(ierr); if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"ndof\"");
  ierr = cJSONGetPetscInt(comm,jdmda,"stencilWidth",sw,&found);CHKERRQ(ierr); if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"stencilWidth\"");
  ierr = cJSONGetPetscString(comm,jdmda,"stencilType",stencilTypeName,&found);CHKERRQ(ierr); if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"stencilType\"");
  PetscStrcmp(stencilTypeName,"star",&match);   if (match) { *stencil_type = DMDA_STENCIL_STAR; goto foundValidStencil; }
  PetscStrcmp(stencilTypeName,"box",&match);    if (match) { *stencil_type = DMDA_STENCIL_BOX;  goto foundValidStencil; }
  if (!match) SETERRQ(comm,PETSC_ERR_USER,"Failed to match any valid DMDAStencilType");
  foundValidStencil:
  
  nf = (int)(*dim);
  if (jdmda) {
    dirlist = cJSON_GetObjectItem(jdmda,"directions");
    if (!dirlist) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to locate key \"directions\"\n");
    nf = cJSON_GetArraySize(dirlist);
    d_k = cJSON_GetArrayItemRoot(dirlist);
  }
  
  for (k=0; k<nf; k++) {
    ierr = cJSONGetPetscInt(comm,d_k,"points",&M[k],&found);CHKERRQ(ierr); if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"points\"");
    
    ierr = cJSONGetPetscString(comm,d_k,"boundaryType",boundaryTypeName,&found);CHKERRQ(ierr); if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"boundaryType\"");
    PetscStrcmp(boundaryTypeName,"none",&match);        if (match) { bc_type[k] = DM_BOUNDARY_NONE;     goto foundValidBoundary; }
    PetscStrcmp(boundaryTypeName,"ghosted",&match);     if (match) { bc_type[k] = DM_BOUNDARY_GHOSTED;  goto foundValidBoundary; }
    PetscStrcmp(boundaryTypeName,"mirror",&match);      if (match) { bc_type[k] = DM_BOUNDARY_MIRROR;   goto foundValidBoundary; }
    PetscStrcmp(boundaryTypeName,"periodic",&match);    if (match) { bc_type[k] = DM_BOUNDARY_PERIODIC; goto foundValidBoundary; }
    if (!match) SETERRQ(comm,PETSC_ERR_USER,"Failed to match any valid DMBoundaryType");
    foundValidBoundary:
    
    if (d_k) { d_k = cJSON_GetArrayItemNext(d_k); }
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode _DMDAParsePartition_JSON(cJSON *jdmda,MPI_Comm comm,PetscInt dim,
                                     PetscInt mr[],
                                     PetscInt **lri,PetscInt **lrj,PetscInt **lrk,
                                     PetscMPIInt *csize)
{
  PetscErrorCode ierr;
  cJSON *part = NULL,*dirlist = NULL,*d_k = NULL;
  PetscBool found;
  int k,nf;
  PetscInt commsize,nvalues;
  
  if (jdmda) {
    part = cJSON_GetObjectItem(jdmda,"partition");
    if (!part) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to locate key \"partition\"\n");
  }
  
  ierr = cJSONGetPetscInt(comm,part,"commSize",&commsize,&found);CHKERRQ(ierr); if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"commSize\"");
  *csize = (PetscMPIInt)commsize;
  
  nf = (int)dim;
  if (part) {
    dirlist = cJSON_GetObjectItem(part,"directions");
    if (!dirlist) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to locate key \"directions\"\n");
    nf = cJSON_GetArraySize(dirlist);
    d_k = cJSON_GetArrayItemRoot(dirlist);
  }
  
  for (k=0; k<nf; k++) {
    ierr = cJSONGetPetscInt(comm,d_k,"ranks",&mr[k],&found);CHKERRQ(ierr);
    if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"ranks\"");

    if (k == 0) {
      ierr = cJSONGetPetscIntArray(comm,d_k,"pointsPerRank",&nvalues,lri,&found);CHKERRQ(ierr);
      if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"pointsPerRank\"");
    }
    else if (k == 1) {
      ierr = cJSONGetPetscIntArray(comm,d_k,"pointsPerRank",&nvalues,lrj,&found);CHKERRQ(ierr);
      if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"pointsPerRank\"");
    }
    else if (k == 2) {
      ierr = cJSONGetPetscIntArray(comm,d_k,"pointsPerRank",&nvalues,lrk,&found);CHKERRQ(ierr);
      if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"pointsPerRank\"");
    } else {
      SETERRQ(comm,PETSC_ERR_USER,"\"direction\" can only contain at most three entries \n");
    }
    
    if (d_k) { d_k = cJSON_GetArrayItemNext(d_k); }
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode _DMDAParseRefine_JSON(cJSON *jdmda,MPI_Comm comm,PetscInt dim,
                                     PetscInt refine[])
{
  PetscErrorCode ierr;
  cJSON *dirlist = NULL,*d_k = NULL;
  PetscBool found;
  int k,nf;
  
  nf = (int)dim;
  if (jdmda) {
    dirlist = cJSON_GetObjectItem(jdmda,"directions");
    if (!dirlist) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to locate key \"directions\"\n");
    nf = cJSON_GetArraySize(dirlist);
    d_k = cJSON_GetArrayItemRoot(dirlist);
  }
  
  for (k=0; k<nf; k++) {
    ierr = cJSONGetPetscInt(comm,d_k,"refinementFactor",&refine[k],&found);CHKERRQ(ierr); if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"refinementFactor\"");
    
    if (d_k) { d_k = cJSON_GetArrayItemNext(d_k); }
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode DMDACheckpointLoad(MPI_Comm comm,const char jfilename[],DM *_da)
{
  PetscErrorCode ierr;
  PetscMPIInt nproc,rank,csize_file = -1;
  cJSON *jfile = NULL,*jdmda = NULL,*jdmdacoor = NULL;
  char optionsprefix[PETSC_MAX_PATH_LEN],cfilename[PETSC_MAX_PATH_LEN];
  DM da = NULL;
  PetscInt dim,ndof,sw,M[3],mr[3],*lri,*lrj,*lrk,refine[3];
  DMBoundaryType bc_type[3];
  DMDAStencilType stencil_type;
  PetscBool found,has_coords;
  
  
  ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  if (rank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to open JSON file \"%s\"",jfilename);
    jdmda = cJSON_GetObjectItem(jfile,"DMDA");
  }

  dim = 0;
  bc_type[0] = bc_type[1] = bc_type[2] = DM_BOUNDARY_NONE;
  stencil_type = DMDA_STENCIL_STAR;
  M[0] = M[1] = M[2] = 0;
  mr[0] = mr[1] = mr[2] = PETSC_DECIDE;
  ndof = 1;
  sw = 0;
  lri = lrj = lrk = NULL;
  
  refine[0] = refine[1] = refine[2] = 2;
  optionsprefix[0] = '\0';
  
  ierr = _DMDAParseLayout_JSON(jdmda,comm,optionsprefix,&dim,bc_type,&stencil_type,M,&ndof,&sw);CHKERRQ(ierr);
  ierr = _DMDAParsePartition_JSON(jdmda,comm,dim,mr,&lri,&lrj,&lrk,&csize_file);CHKERRQ(ierr);
  ierr = _DMDAParseRefine_JSON(jdmda,comm,dim,refine);CHKERRQ(ierr);

  if (csize_file != nproc) {
    /* reset all information regarding the partition */
    PetscPrintf(comm,"[DMDACheckpointLoad]: Ignoring partition information from file. Current comm.size = %d, whilst data was written using comm.size = %d\n",(int)nproc,(int)csize_file);
    mr[0] = mr[1] = mr[2] = PETSC_DECIDE;
    
    ierr = PetscFree(lri);CHKERRQ(ierr);
    ierr = PetscFree(lrj);CHKERRQ(ierr);
    ierr = PetscFree(lrk);CHKERRQ(ierr);
    lri = lrj = lrk = NULL;
  }
  
  switch (dim) {
    case 1:
      SETERRQ(comm,PETSC_ERR_SUP,"1D requires code adjustment: Only tested for 3D");
      break;
    case 2:
      SETERRQ(comm,PETSC_ERR_SUP,"2D requires testing: Only tested for 3D");
      ierr = DMDACreate2d(comm,bc_type[0],bc_type[1],stencil_type,M[0],M[1],mr[0],mr[1],ndof,sw,lri,lrj,&da);CHKERRQ(ierr);
      ierr = DMSetUp(da);CHKERRQ(ierr);
      break;
    case 3:
      ierr = DMDACreate3d(comm,bc_type[0],bc_type[1],bc_type[2],stencil_type,M[0],M[1],M[2],mr[0],mr[1],mr[2],ndof,sw,lri,lrj,lrk,&da);CHKERRQ(ierr);
      ierr = DMSetUp(da);CHKERRQ(ierr);
      break;
      
    default:
      SETERRQ(comm,PETSC_ERR_SUP,"Spatial dimension must be 1, 2 or 3");
      break;
  }
  
  /* load coordinate file if found */
  cfilename[0] = '\0';
  has_coords = PETSC_FALSE;
  if (rank == 0) {
    jdmdacoor = cJSON_GetObjectItem(jdmda,"coordinate");
    has_coords = PETSC_TRUE;
  }
  ierr = MPI_Bcast(&has_coords,1,MPIU_BOOL,0,comm);CHKERRQ(ierr);
  
  if (has_coords) {
    ierr = cJSONGetPetscString(comm,jdmdacoor,"dmCoordinateFile",cfilename,&found);CHKERRQ(ierr);
    if (!found) SETERRQ(comm,PETSC_ERR_USER,"Failed to locate key \"dmCoordinateFile\"");
    
    ierr = DMDALoadCoordinatesFromFile(da,cfilename);CHKERRQ(ierr);
  }
  
  /* set prefix and refinement factor */
  if (da) {
    if (optionsprefix[0]) {
      ierr = DMSetOptionsPrefix(da,optionsprefix);CHKERRQ(ierr);
    }
    ierr = DMDASetRefinementFactor(da,refine[0],refine[1],refine[2]);CHKERRQ(ierr);
  }
  
  ierr = PetscFree(lri);CHKERRQ(ierr);
  ierr = PetscFree(lrj);CHKERRQ(ierr);
  ierr = PetscFree(lrk);CHKERRQ(ierr);
  
  if (jfile) { cJSON_Delete(jfile); }
  
  *_da = da;
  
  PetscFunctionReturn(0);
}

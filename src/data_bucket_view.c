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
 **    filename:   data_bucket_view.c
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
#include <cJSON.h>
#include <mpiio_blocking.h>
#include <data_bucket.h>
#include <geometry_object_parse.h>
#include <data_bucket_view.h>


#undef __FUNCT__
#define __FUNCT__ "DataBucketWriteMeta"
PetscErrorCode DataBucketWriteMeta(FILE *fp,long int len,DataBucket db,PetscBool open)
{
  if (fp) {
    if (open) {
      fprintf(fp,"{\n\"DataBucket\": {\n");
      //fprintf(fp,"  \"Object\":        \"DataBucket\",\n");
      fprintf(fp,"  \"DBVersion\":     \"1.0.0\",\n");
      fprintf(fp,"  \"DBStorageType\": \"native\",\n");
      //fprintf(fp,"  \"DBLength\":      %ld,\n",len);
      fprintf(fp,"  \"DBLength\":      %d,\n",db->L);
      fprintf(fp,"  \"DBBuffer\":      %d,\n",db->buffer);
      fprintf(fp,"  \"DBAllocated\":   %d,\n",db->allocated);
      fprintf(fp,"  \"DBNFields\":     %d,\n",db->nfields);
    } else {
      fprintf(fp,"}\n}\n");
    }
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DataBucketWriteMPIMeta"
PetscErrorCode DataBucketWriteMPIMeta(FILE *fp,DataBucket db,MPI_Comm comm)
{
  PetscMPIInt p,rank,commsize;
  int L,allocated,buffer,*list;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);

	DataBucketGetSizes(db,&L,&buffer,&allocated);
  
  list = NULL;
  if (rank == 0) {
    PetscMalloc(sizeof(int)*commsize,&list);
  }
  
  ierr = MPI_Gather(&L,1,MPI_INT,list,1,MPI_INT,0,comm);CHKERRQ(ierr);
  if (rank == 0) {
    fprintf(fp,"  \"DBRanks\":       %d,\n",commsize);
    fprintf(fp,"  \"DBRankLength\":  [");
    for (p=0; p<commsize-1; p++) {
      fprintf(fp,"%d,",list[p]);
    }
    fprintf(fp,"%d],\n",list[commsize-1]);
  }

  ierr = MPI_Gather(&buffer,1,MPI_INT,list,1,MPI_INT,0,comm);CHKERRQ(ierr);
  if (rank == 0) {
    fprintf(fp,"  \"DBRankBuffer\":  [");
    for (p=0; p<commsize-1; p++) {
      fprintf(fp,"%d,",list[p]);
    }
    fprintf(fp,"%d],\n",list[commsize-1]);
  }

  ierr = MPI_Gather(&allocated,1,MPI_INT,list,1,MPI_INT,0,comm);CHKERRQ(ierr);
  if (rank == 0) {
    fprintf(fp,"  \"DBRankAllocated\":  [");
    for (p=0; p<commsize-1; p++) {
      fprintf(fp,"%d,",list[p]);
    }
    fprintf(fp,"%d],\n",list[commsize-1]);
  }

  
  
  if (list) { PetscFree(list); }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DataBucketWriteAllFieldMeta"
PetscErrorCode DataBucketWriteAllFieldMeta(FILE *fp,long int len,DataBucket db,const char fieldprefix[],MPI_Comm comm)
{
  char           filename[PETSC_MAX_PATH_LEN];
  int            f;
  PetscMPIInt    rank;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (rank == 0) {
    for (f=0; f<db->nfields; f++) {
      fprintf(fp,"\n");
      
      fprintf(fp,"    \"DataField\":  {\n");
      fprintf(fp,"      \"Name\":        \"%s\",\n",db->field[f]->name);
      fprintf(fp,"      \"Index\":       %d,\n",f);
      fprintf(fp,"      \"AtomicSize\":  %zu,\n",db->field[f]->atomic_size);
      fprintf(fp,"      \"Length\":      %ld,\n",len);
      fprintf(fp,"      \"RegisterationFunction\":  \"%s\",\n",db->field[f]->registeration_function);

      if (fieldprefix) {
        PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s%s_db.bin",fieldprefix,db->field[f]->name);
      } else {
        PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s_db.bin",db->field[f]->name);
      }
      fprintf(fp,"      \"DataFile\":    \"%s\"\n",filename);

      if (f == (db->nfields-1)) { fprintf(fp,"    }\n\n");}
      else {                      fprintf(fp,"    },\n");}
    }
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DataBucketWriteFieldDataFiles"
PetscErrorCode DataBucketWriteFieldDataFiles(DataBucket db,const char fieldprefix[],MPI_Comm comm)
{
  char           filename[PETSC_MAX_PATH_LEN];
  int            f;
  PetscMPIInt    rank;
  PetscErrorCode ierr;
  FILE           *fp;
  
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  for (f=0; f<db->nfields; f++) {
    DataField df;
    long int  len;
    
    if (fieldprefix) {
      PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s%s_db.bin",fieldprefix,db->field[f]->name);
    } else {
      PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s_db.bin",db->field[f]->name);
    }

    df = db->field[f];
    /* 
      Note:
        We ONLY write the valid range of the data field and
        excluding writing anything in the buffer.
        Upon loading, the same buffer size will be allocated,
        however it's contents will be zeroed
    */
    //len = (long int)df->L;
    len = (long int)db->L;
    
    fp = NULL;
    if (rank == 0) {
      fp = fopen(filename,"w");
      if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to open file (write) %s",filename);
    }
    ierr = MPIWrite_Blocking(fp,df->data,len,df->atomic_size,0,PETSC_FALSE,comm);CHKERRQ(ierr);
    if (rank == 0) { fclose(fp); }
  }
  PetscFunctionReturn(0);
}

void _cJSON_GetObjectValue_intarray(cJSON *cj,const char name[],int *found,int *nv,int *vals[])
{
  cJSON *list = NULL;
  int   k,len,*_v;
  
  *found = cJSON_True;
  list = cJSON_GetObjectItem(cj,name);
  if (!list) {
    *nv = 0;
    *found = cJSON_False;
    return;
  }
  
  len = cJSON_GetArraySize(list);
  _v = malloc(sizeof(int)*len);
  for (k=0; k<len; k++) {
    cJSON *list_entry;
    
    list_entry = cJSON_GetArrayItem(list,k);
    _v[k] = list_entry->valueint;
  }
  
  *nv = len;
  *vals = _v;
}

#undef __FUNCT__
#define __FUNCT__ "DataBucketReadDataFieldMeta"
PetscErrorCode DataBucketReadDataFieldMeta(DataBucket db,cJSON *root,MPI_Comm comm)
{
  PetscErrorCode ierr;
  char *j_name = NULL,*j_regfunc = NULL,*j_datafname;
  int j_asize = 0;
  size_t atomicsize,slen = 0;
  int found;
  char name[PETSC_MAX_PATH_LEN],func[PETSC_MAX_PATH_LEN];
  DataField gfield;
  long int len;
  FILE *fp;
  PetscMPIInt rank;
  
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (rank == 0) {
    cJSON_GetObjectValue_char(root,"Name",&found,&j_name); //printf("Name %s\n",j_name);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"DataField: Missing field \"Name\"");
    
    cJSON_GetObjectValue_int(root,"AtomicSize",&found,&j_asize); //printf("AtomicSize %d\n",j_asize);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"DataField: Missing field \"AtomicSize\"");

    //cJSON_GetObjectValue_int(root,"Length",&found,&j_length); printf("Length %d\n",j_length);

    cJSON_GetObjectValue_char(root,"RegisterationFunction",&found,&j_regfunc); //printf("Function %s\n",j_regfunc);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"DataField: Missing field \"RegisterationFunction\"");

    cJSON_GetObjectValue_char(root,"DataFile",&found,&j_datafname); //printf("DataFile %s\n",j_datafname);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"DataField: Missing field \"DataFile\"");
  }

  PetscMemzero(name,sizeof(char)*PETSC_MAX_PATH_LEN);
  if (rank == 0) {
    ierr = PetscStrcpy(name,j_name);CHKERRQ(ierr);
    ierr = PetscStrlen(name,&slen);CHKERRQ(ierr);
    //printf("PetscStrlen %d \n",slen);
  }
  ierr = MPI_Bcast(name,PETSC_MAX_PATH_LEN,MPI_CHAR,0,comm);CHKERRQ(ierr);
  //printf("[bc] name %s \n",name);
  
  ierr = MPI_Bcast(&j_asize,1,MPI_INT,0,comm);CHKERRQ(ierr);
  atomicsize = (size_t)j_asize;
  //printf("[bc] j_asize %d \n",j_asize);
  
  PetscMemzero(func,sizeof(char)*PETSC_MAX_PATH_LEN);
  if (rank == 0) {
    ierr = PetscStrcpy(func,j_regfunc);CHKERRQ(ierr);
    ierr = PetscStrlen(func,&slen);CHKERRQ(ierr);
  }
  ierr = MPI_Bcast(func,PETSC_MAX_PATH_LEN,MPI_CHAR,0,comm);CHKERRQ(ierr);
  //printf("[bc] func %s \n",func);
  
	DataBucketRegisterField(db,name,atomicsize,&gfield);
  /* overide registered function name */
  free(gfield->registeration_function);
	asprintf(&gfield->registeration_function, "%s", func );
  
  /* 
    Note:
      We ONLY read the valid range of the data field and
      excluding writing anything in the buffer.
      After the valid range was been loaded, the buffer
      space will contain zero entries.
  */
  //len = (long int)gfield->L;
  len = (long int)db->L;
  fp = NULL;
  if (rank == 0) {
    fp = fopen(j_datafname,"r");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to open file (read) %s",j_datafname);
  }
  ierr = MPIRead_Blocking(fp,&gfield->data,len,gfield->atomic_size,0,PETSC_FALSE,comm);CHKERRQ(ierr);

  if (rank == 0) { fclose(fp); }
  
  PetscFunctionReturn(0);
}

/* Load / register field specification first */
/* Allocate space */
/* Load contents of data fields */
#undef __FUNCT__
#define __FUNCT__ "DataBucketReadMeta"
PetscErrorCode DataBucketReadMeta(DataBucket db,cJSON *jfile,MPI_Comm comm)
{
  cJSON *root,*dfitem;
  int   found,j_length = 9,j_buffer = 0,j_alloc = 0,j_ranks = 0,j_nfields = 0;
  int   count,f,nv,*list1 = NULL,*list2 = NULL,*list3 = NULL;
  cJSON **df_list = NULL;
  int length,buffer,allocated;
  PetscMPIInt r,commsize,rank,avg,rem;
  PetscErrorCode ierr;

  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  /* check DB is fresh and no fields have been registered */
  if (db->nfields != 0) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"The DataBucket object should not have any fields registered prior to reading data from a file");
  }
  
  if (jfile) {
    root = cJSON_GetArrayItemRoot(jfile);
  }
  
  if (jfile) {
    cJSON_GetObjectValue_int(root,"DBLength",&found,&j_length); //printf("len %d\n",j_length);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"DataBucket: Missing field \"DBLength\"");

    cJSON_GetObjectValue_int(root,"DBBuffer",&found,&j_buffer); //printf("buff %d\n",j_buffer);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"DataBucket: Missing field \"DBBuffer\"");

    cJSON_GetObjectValue_int(root,"DBAllocated",&found,&j_alloc); //printf("all %d\n",j_alloc);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"DataBucket: Missing field \"DBAllocated\"");

    cJSON_GetObjectValue_int(root,"DBRanks",&found,&j_ranks); //printf("nr %d\n",j_ranks);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"DataBucket: Missing field \"DBRanks\"");
    
    //if (j_ranks > 1) { // read values from all ranks //
      _cJSON_GetObjectValue_intarray(root,"DBRankLength",&found,&nv,&list1);
      _cJSON_GetObjectValue_intarray(root,"DBRankBuffer",&found,&nv,&list2);
      _cJSON_GetObjectValue_intarray(root,"DBRankAllocated",&found,&nv,&list3);
    //}
  }
 
  ierr = MPI_Bcast(&j_ranks,1,MPI_INT,0,comm);CHKERRQ(ierr);

  length = buffer = allocated = 0;
  if (commsize == 1) {
    if (j_ranks > 1) { PetscPrintf(comm,"[DataBucketReadMeta] Serializing parallel DataBucket data\n"); }
    for (r=0; r<j_ranks; r++) {
      length += list1[r];
      buffer += list2[r];
      allocated += list3[r];
    }
    PetscPrintf(PETSC_COMM_SELF,"  [rank 0] L %d ; b %d ; a %d \n",length,buffer,allocated);
  } else if (commsize != j_ranks) {
    int *_list1,*_list2,*_list3;
    
    /* define a new partition */
    PetscPrintf(comm,"[DataBucketReadMeta] comm.size differs from comm.size of the parallel DataBucket data\n");
    PetscPrintf(comm,"[DataBucketReadMeta] Defining a new partition\n");

    PetscMalloc1(j_ranks,&_list1);
    PetscMalloc1(j_ranks,&_list2);
    PetscMalloc1(j_ranks,&_list3);
    for (r=0; r<j_ranks; r++) {
      _list1[r] = _list2[r] = _list3[r] = 0;
    }
    
    if (rank == 0) {
      PetscPrintf(PETSC_COMM_SELF,"[DataBucketReadMeta] Stored partition\n");
      for (r=0; r<j_ranks; r++) {
        _list1[r] = list1[r];
        _list2[r] = list2[r];
        _list3[r] = list3[r];
        PetscPrintf(PETSC_COMM_SELF,"  [rank %d] L %d ; b %d ; a %d \n",r,list1[r],list2[r],list3[r]);
      }
    }
    
    PetscPrintf(comm,"[DataBucketReadMeta] New partition\n");
    
    ierr = MPI_Bcast(_list1,j_ranks,MPI_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(_list2,j_ranks,MPI_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(_list3,j_ranks,MPI_INT,0,comm);CHKERRQ(ierr);
    
    if (j_ranks > commsize) {
      avg = j_ranks / commsize;
      rem = j_ranks - avg * commsize;
    } else {
      avg = commsize / j_ranks;
      rem = commsize - avg * j_ranks;
      PetscPrintf(comm,"[DataBucketReadMeta] comm.size matches comm.size of the parallel DataBucket data\n");
      SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"[DataBucketReadMeta] Cannot load a DataBucket object using a communicator with more ranks than that used to write the data [current %d : stored %d]",(int)commsize,j_ranks);
    }
    
    for (r=0; r<commsize; r++) {
      int k,start,end,sum[3];
      
      start = r * avg;
      end = start + avg;
      if (r == (commsize-1)) end += rem;
      sum[0] = sum[1] = sum[2] = 0;
      for (k=start; k<end; k++) {
        sum[0] += _list1[k];
        sum[1] += _list2[k];
        sum[2] += _list3[k];
      }
      if (rank == r) {
        length = sum[0];
        buffer = sum[1];
        allocated = sum[2];
      }
      PetscPrintf(comm,"  [rank %d] L %d ; b %d ; a %d \n",r,sum[0],sum[1],sum[2]);
    }
    
    PetscFree(_list3);
    PetscFree(_list2);
    PetscFree(_list1);
  } else {
    PetscPrintf(comm,"[DataBucketReadMeta] comm.size matches comm.size of the parallel DataBucket data\n");
    if (rank == 0) {
      for (r=0; r<j_ranks; r++) {
        PetscPrintf(PETSC_COMM_SELF,"  [rank %d] L %d ; b %d ; a %d \n",r,list1[r],list2[r],list3[r]);
      }
    }
    ierr = MPI_Scatter(list1,1,MPI_INT,&length,1,MPI_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Scatter(list2,1,MPI_INT,&buffer,1,MPI_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Scatter(list3,1,MPI_INT,&allocated,1,MPI_INT,0,comm);CHKERRQ(ierr);
  }

  db->L = length;
  db->buffer = buffer;
  db->allocated = allocated;
  
  //printf("[rank %d] L %d ; b %d ; a %d \n",rank,length,buffer,allocated);
  
  dfitem = NULL;
  if (jfile) {
    cJSON_GetObjectValue_int(root,"DBNFields",&found,&j_nfields); //printf("nf %d\n",j_nfields);
  }

  ierr = MPI_Bcast(&j_nfields,1,MPI_INT,0,comm);CHKERRQ(ierr);
  PetscMalloc1(j_nfields,&df_list);
  for (f=0; f<j_nfields; f++) {
    df_list[f] = NULL;
  }
  
  if (jfile) {
    dfitem = cJSON_GetObjectItem(root,"DataField");
    if (!dfitem) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"No DataItem's were found in DataBucket meta file");
    df_list[0] = dfitem;
  }
  
  count = 1;
  for (f=1; f<j_nfields; f++) {
    if (dfitem && dfitem->next) {
      dfitem = dfitem->next;
      df_list[f] = dfitem;
      count++;
    }
  }
  ierr = MPI_Bcast(&count,1,MPI_INT,0,comm);CHKERRQ(ierr);
  if (count != j_nfields) {
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"Number of DataField's found (%D) does not match number expected (%D)",count,j_nfields);
  }
  
  ierr = MPI_Barrier(comm);CHKERRQ(ierr);
  for (f=0; f<j_nfields; f++) {
//    for (f=1; f<2; f++) {
    ierr = DataBucketReadDataFieldMeta(db,df_list[f],comm);CHKERRQ(ierr);
  }

  
  PetscFree(df_list);
  if (list3) { free(list3); }
  if (list2) { free(list2); }
  if (list1) { free(list1); }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DataBucketLoad_Native"
PetscErrorCode DataBucketLoad_Native(DataBucket db,const char filename[],MPI_Comm comm)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  cJSON          *jfile = NULL;

  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  // open json file
  if (rank == 0) {
    cJSON_FileView(filename,&jfile);
    if (!jfile) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open JSON file: %s",filename);
  }
  
  // load meta
  ierr = DataBucketReadMeta(db,jfile,comm);CHKERRQ(ierr);
  
  // load heavy data
  
  
  if (jfile) {
    cJSON_Delete(jfile);
  }
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DataBucketView_Native"
PetscErrorCode DataBucketView_Native(DataBucket db,const char prefix[],MPI_Comm comm)
{
  int         _len;
  long int    len,len_g;
  FILE        *fp;
  PetscMPIInt rank;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  DataBucketGetSizes(db,&_len,NULL,NULL);
  len = (long int)_len;
  len_g = 0;
  ierr = MPI_Reduce(&len,&len_g,1,MPI_LONG,MPI_SUM,0,comm);CHKERRQ(ierr);
  
  fp = NULL;
  if (rank == 0) {
    fp = fopen("db.meta","w");
    if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open binary file: XXX");
  }
  ierr = DataBucketWriteMeta(fp,len_g,db,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DataBucketWriteMPIMeta(fp,db,comm);CHKERRQ(ierr);
  
  ierr = DataBucketWriteAllFieldMeta(fp,len_g,db,prefix,comm);CHKERRQ(ierr);
  
  ierr = DataBucketWriteFieldDataFiles(db,prefix,comm);CHKERRQ(ierr);

  ierr = DataBucketWriteMeta(fp,len_g,db,PETSC_FALSE);CHKERRQ(ierr);
  if (fp) { fclose(fp); }
  
  PetscFunctionReturn(0);
}

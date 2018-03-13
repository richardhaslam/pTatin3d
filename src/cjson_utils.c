
#include <stdio.h>
#include <stdlib.h>
#include <petsc.h>
#include <cjson_utils.h>

/*
 If we only want to go the next item in the list (in sequence),
 it is pointless to use the function cJSON_GetArrayItem().
 This function traverses the entire tree from the root upon each call.
 
 Usage:
 gobject = cJSON_GetObjectItem(jfile,"GeometryObjectList");
 ngobj = cJSON_GetArraySize(gobject);
 
 gobj_k = cJSON_GetArrayItemRoot(gobject);
 for (k=0; k<ngobj; k++) {
   gobj_k = cJSON_GetArrayItemNext(gobj_k);
 }
*/
cJSON* cJSON_GetArrayItemRoot(cJSON *gobject)
{
  if (gobject) {
    return cJSON_GetArrayItem(gobject,0);
  } else {
    return NULL;
  }
}

cJSON* cJSON_GetArrayItemNext(cJSON *gobj_k)
{
  if (gobj_k) {
    return gobj_k->next;
  } else {
    return NULL;
  }
}

void cJSON_FileView(const char filename[],cJSON **jf)
{
  FILE  *fp = NULL;
  char  *data;
  long  len;
  cJSON *jfile;
  
  fp = fopen(filename,"rb");
  if (!fp) {
    *jf = NULL;
    return;
  }
  
  fseek(fp,0,SEEK_END);
  
  len = ftell(fp);
  fseek(fp,0,SEEK_SET);
  
  data = malloc(sizeof(char)*(len+1));
  if (fread(data,1,len,fp) != len) {
    printf("fread error (%s:%d). Exiting ungracefully.\n",__FILE__,__LINE__);
    exit(1);
  }
  fclose(fp);
  
  jfile = NULL;
  jfile = cJSON_Parse(data);
  free(data);
  
  if (!jfile) {
    printf("***************************** cJSON PARSING ERROR *****************************\n%s\n**********************************************************\n",cJSON_GetErrorPtr());
  }
  
  *jf = jfile;
}

void cJSON_GetObjectValue_bool(cJSON *cj,const char name[],int *found,int *val)
{
  cJSON *obj = NULL;
  
  *found = cJSON_True;
  obj = cJSON_GetObjectItem(cj,name);
  if (!obj) {
    *val = -1;
    *found = cJSON_False;
    return;
  }
  
  *val = (int)obj->type;
}

void cJSON_GetObjectValue_int(cJSON *cj,const char name[],int *found,int *val)
{
  cJSON *obj = NULL;
  
  *found = cJSON_True;
  obj = cJSON_GetObjectItem(cj,name);
  if (!obj) {
    *val = 0;
    *found = cJSON_False;
    return;
  }
  
  *val = obj->valueint;
}

void cJSON_GetObjectValue_intarray(cJSON *cj,const char name[],int *found,int *nv,int vals[])
{
  cJSON *list = NULL;
  int   k,len;
  
  *found = cJSON_True;
  list = cJSON_GetObjectItem(cj,name);
  if (!list) {
    *nv = 0;
    *found = cJSON_False;
    return;
  }
  
  len = cJSON_GetArraySize(list);
  for (k=0; k<len; k++) {
    cJSON *list_entry;
    
    list_entry = cJSON_GetArrayItem(list,k);
    vals[k] = list_entry->valueint;
  }
  
  *nv = len;
}

void cJSON_GetObjectValue_double(cJSON *cj,const char name[],int *found,double *val)
{
  cJSON *obj = NULL;
  
  *found = cJSON_True;
  obj = cJSON_GetObjectItem(cj,name);
  if (!obj) {
    *val = 0.0;
    *found = cJSON_False;
    return;
  }
  
  *val = obj->valuedouble;
}

void cJSON_GetObjectValue_doublearray(cJSON *cj,const char name[],int *found,int *nv,double vals[])
{
  cJSON *list = NULL;
  int   k,len;
  
  *found = cJSON_True;
  list = cJSON_GetObjectItem(cj,name);
  if (!list) {
    *nv = 0;
    *found = cJSON_False;
    return;
  }
  
  len = cJSON_GetArraySize(list);
  for (k=0; k<len; k++) {
    cJSON *list_entry;
    
    list_entry = cJSON_GetArrayItem(list,k);
    vals[k] = list_entry->valuedouble;
  }
  
  *nv = len;
}

void cJSON_GetObjectValue_char(cJSON *cj,const char name[],int *found,char **val)
{
  cJSON *obj = NULL;
  
  *found = cJSON_True;
  obj = cJSON_GetObjectItem(cj,name);
  if (!obj) {
    *val = NULL;
    *found = cJSON_False;
    return;
  }
  *val = obj->valuestring;
}

cJSON *cJSON_CreateInt(int num)
{
  cJSON *blob;
  
  blob = cJSON_CreateNumber((double)num);
  return(blob);
}

cJSON *cJSON_CreateDouble(double num)
{
  cJSON *blob;
  
  blob = cJSON_CreateNumber((double)num);
  return(blob);
}

void cJSON_GetObjectValue_arraylength(cJSON *cj,const char name[],int *found,int *nv)
{
  cJSON *list = NULL;
  int   len;
  
  *found = cJSON_True;
  list = cJSON_GetObjectItem(cj,name);
  if (!list) {
    *nv = 0;
    *found = cJSON_False;
    return;
  }
  
  len = cJSON_GetArraySize(list);
  *nv = len;
}

#undef __FUNCT__
#define __FUNCT__ "cJSONGetPetscInt"
PetscErrorCode cJSONGetPetscInt(MPI_Comm comm,cJSON *cj,const char name[],PetscInt *val,PetscBool *found)
{
  PetscMPIInt commsize,commrank;
  int cfound = cJSON_False;
  int ival = 0;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  *found = PETSC_FALSE;
  if (commrank == 0) {
    cJSON_GetObjectValue_int(cj,name,&cfound,&ival);
    *val = (PetscInt)ival;
  }
  if (commsize > 1) {
    /* broadcast value and found */
    ierr = MPI_Bcast(val,1,MPIU_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(&cfound,1,MPI_INT,0,comm);CHKERRQ(ierr);
  }
  if (cfound == cJSON_True) { *found  = PETSC_TRUE; }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "cJSONGetPetscReal"
PetscErrorCode cJSONGetPetscReal(MPI_Comm comm,cJSON *cj,const char name[],PetscReal *val,PetscBool *found)
{
  PetscMPIInt commsize,commrank;
  int cfound = cJSON_False;
  double dval = 0.0;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  *found = PETSC_FALSE;
  if (commrank == 0) {
    cJSON_GetObjectValue_double(cj,name,&cfound,&dval);
    *val = (PetscReal)dval;
  }
  if (commsize > 1) {
    /* broadcast value and found */
    ierr = MPI_Bcast(val,1,MPIU_REAL,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(&cfound,1,MPI_INT,0,comm);CHKERRQ(ierr);
  }
  if (cfound == cJSON_True) { *found  = PETSC_TRUE; }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "cJSONGetPetscRealArray"
PetscErrorCode cJSONGetPetscRealArray(MPI_Comm comm,cJSON *cj,const char name[],PetscInt *nvalues,PetscReal **values,PetscBool *found)
{
  PetscMPIInt commsize,commrank;
  int cfound = cJSON_False;
  int k,nvals = 0;
  double *dvals = NULL;
  PetscReal *_values = NULL;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  *found = PETSC_FALSE;
  if (commrank == 0) {
    
    cJSON_GetObjectValue_arraylength(cj,name,&cfound,&nvals);
    dvals = (double*)malloc(sizeof(double)*(nvals + 1));
    if (cfound == cJSON_True) {
      cJSON_GetObjectValue_doublearray(cj,name,&cfound,&nvals,dvals);
    }
  }

  /* broadcast nvals,ivals, and found */
  ierr = MPI_Bcast(&nvals,1,MPI_INT,0,comm);CHKERRQ(ierr);
  *nvalues = (PetscInt)nvals;
  
  ierr = PetscMalloc1(*nvalues,&_values);CHKERRQ(ierr);
  ierr = PetscMemzero(_values,sizeof(PetscReal)*(*nvalues));CHKERRQ(ierr);
  
  if (commrank == 0) {
    for (k=0; k<nvals; k++) {
      _values[k] = (PetscReal)dvals[k];
    }
  }
  ierr = MPI_Bcast(_values,*nvalues,MPIU_REAL,0,comm);CHKERRQ(ierr);
  ierr = MPI_Bcast(&cfound,1,MPI_INT,0,comm);CHKERRQ(ierr);
  
  *values = _values;
  if (cfound == cJSON_True) { *found  = PETSC_TRUE; }
  if (dvals) { free(dvals); }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "cJSONGetPetscIntArray"
PetscErrorCode cJSONGetPetscIntArray(MPI_Comm comm,cJSON *cj,const char name[],PetscInt *nvalues,PetscInt **values,PetscBool *found)
{
  PetscMPIInt commsize,commrank;
  int cfound = cJSON_False;
  int k,nvals = 0;
  int *ivals = NULL;
  PetscInt *_values = NULL;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  *found = PETSC_FALSE;
  if (commrank == 0) {

    cJSON_GetObjectValue_arraylength(cj,name,&cfound,&nvals);
    ivals = (int*)malloc(sizeof(int)*(nvals + 1));
    
    if (cfound == cJSON_True) {
      cJSON_GetObjectValue_intarray(cj,name,&cfound,&nvals,ivals);
    }
  }

  /* broadcast nvals,ivals, and found */
  ierr = MPI_Bcast(&nvals,1,MPI_INT,0,comm);CHKERRQ(ierr);
  *nvalues = (PetscInt)nvals;

  ierr = PetscMalloc1(*nvalues,&_values);CHKERRQ(ierr);
  ierr = PetscMemzero(_values,sizeof(PetscInt)*(*nvalues));CHKERRQ(ierr);

  if (commrank == 0) {
    for (k=0; k<nvals; k++) {
      _values[k] = (PetscInt)ivals[k];
    }
  }
  ierr = MPI_Bcast(_values,*nvalues,MPIU_INT,0,comm);CHKERRQ(ierr);
  ierr = MPI_Bcast(&cfound,1,MPI_INT,0,comm);CHKERRQ(ierr);
  
  *values = _values;
  if (cfound == cJSON_True) { *found  = PETSC_TRUE; }
  if (ivals) { free(ivals); }

  PetscFunctionReturn(0);
}

/*
 A "petsc string" is simply char val[PETSC_MAX_PATH_LEN];
*/
#undef __FUNCT__
#define __FUNCT__ "cJSONGetPetscString"
PetscErrorCode cJSONGetPetscString(MPI_Comm comm,cJSON *cj,const char name[],char val[],PetscBool *found)
{
  PetscMPIInt commsize,commrank;
  int cfound = cJSON_False;
  PetscInt k;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  *found = PETSC_FALSE;
  for (k=0; k<PETSC_MAX_PATH_LEN; k++) {
    val[k] = '\0';
  }
  if (commrank == 0) {
    char *cval;

    cJSON_GetObjectValue_char(cj,name,&cfound,&cval);
    if (cfound == cJSON_True) {
      ierr = PetscSNPrintf(val,PETSC_MAX_PATH_LEN-1,"%s",cval);CHKERRQ(ierr);
    }
  }
  if (commsize > 1) {
    /* broadcast value and found */
    ierr = MPI_Bcast(val,PETSC_MAX_PATH_LEN,MPI_CHAR,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(&cfound,1,MPI_INT,0,comm);CHKERRQ(ierr);
  }
  if (cfound == cJSON_True) { *found  = PETSC_TRUE; }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "cJSONGetPetscBool"
PetscErrorCode cJSONGetPetscBool(MPI_Comm comm,cJSON *cj,const char name[],PetscBool *val,PetscBool *found)
{
  PetscMPIInt commsize,commrank;
  int cfound = cJSON_False;
  int ival = 0;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  *found = PETSC_FALSE;
  if (commrank == 0) {
    cJSON_GetObjectValue_bool(cj,name,&cfound,&ival);
  }
  if (commsize > 1) {
    /* broadcast value and found */
    ierr = MPI_Bcast(&ival,1,MPIU_INT,0,comm);CHKERRQ(ierr);
    ierr = MPI_Bcast(&cfound,1,MPI_INT,0,comm);CHKERRQ(ierr);
  }
  if (ival == cJSON_True) { *val  = PETSC_TRUE; }
  else                    { *val  = PETSC_FALSE; }
  if (cfound == cJSON_True) { *found  = PETSC_TRUE; }
  PetscFunctionReturn(0);
}

/*
#undef __FUNCT__
#define __FUNCT__ "SETERRQ_JSONKEY"
PetscErrorCode SETERRQ_JSONKEY(MPI_Comm comm,const char keyname[])
{
  SETERRQ1(comm,PETSC_ERR_FILE_UNEXPECTED,"Failed to locate essential JSON key \"%s\"",keyname);
  PetscFunctionReturn(0);
}
*/

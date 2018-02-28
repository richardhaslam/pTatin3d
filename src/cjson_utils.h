
#ifndef __cjson_utils_h__
#define __cjson_utils_h__

#include <cJSON.h>

void cJSON_FileView(const char filename[],cJSON **jf);

void cJSON_GetObjectValue_bool(cJSON *cj,const char name[],int *found,int *val);
void cJSON_GetObjectValue_int(cJSON *cj,const char name[],int *found,int *val);
void cJSON_GetObjectValue_intarray(cJSON *cj,const char name[],int *found,int *nv,int vals[]);
void cJSON_GetObjectValue_double(cJSON *cj,const char name[],int *found,double *val);
void cJSON_GetObjectValue_doublearray(cJSON *cj,const char name[],int *found,int *nv,double vals[]);
void cJSON_GetObjectValue_char(cJSON *cj,const char name[],int *found,char **val);

cJSON* cJSON_GetArrayItemRoot(cJSON *gobject);
cJSON* cJSON_GetArrayItemNext(cJSON *gobj_k);

cJSON *cJSON_CreateInt(int num);
cJSON *cJSON_CreateDouble(double num);

void cJSON_GetObjectValue_arraylength(cJSON *cj,const char name[],int *found,int *nv);

PetscErrorCode cJSONGetPetscInt(MPI_Comm comm,cJSON *cj,const char name[],PetscInt *val,PetscBool *found);
PetscErrorCode cJSONGetPetscReal(MPI_Comm comm,cJSON *cj,const char name[],PetscReal *val,PetscBool *found);
PetscErrorCode cJSONGetPetscIntArray(MPI_Comm comm,cJSON *cj,const char name[],PetscInt *nvalues,PetscInt **values,PetscBool *found);
PetscErrorCode cJSONGetPetscString(MPI_Comm comm,cJSON *cj,const char name[],char val[],PetscBool *found);
PetscErrorCode cJSONGetPetscBool(MPI_Comm comm,cJSON *cj,const char name[],PetscBool *val,PetscBool *found);
/*
PetscErrorCode SETERRQ_JSONKEY(MPI_Comm comm,const char keyname[]);
*/
#define SETERRQ_JSONKEY(comm,keyname) SETERRQ1((comm),PETSC_ERR_FILE_UNEXPECTED,"Failed to locate essential JSON key \"%s\"",(keyname));


#endif

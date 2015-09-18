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
 **    filename:   cJSON_helpers.c
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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "cJSON.h"
#include "cJSON_helpers.h"

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
  fread(data,1,len,fp);
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
  
  *val = obj->type;
}

void cJSON_GetObjectValue_logical(cJSON *cj,const char name[],int *found,int *val)
{
  char *logical = NULL,*upper;
  const char *positive[] = { "YES", "TRUE",  "PETSC_TRUE",  "Y", "T" };
  const char *negative[] = { "NO",  "FALSE", "PETSC_FALSE", "N", "F" };
  int len,cmp,k;
  
  cJSON_GetObjectValue_char(cj,name,found,&logical);
  if (!(*found)) {
    *val = -1;
    return;
  }
  
  asprintf(&upper,"%s",logical);
  len = strlen(upper);
  for (k=0; k<len; k++) {
    upper[k] = toupper(logical[k]);
  }
  
  /* compare with yes/no, t/f, True/False, PETSC_TRUE/PETSC_FALSE, y/n */
  for (k=0; k<5; k++) {
    cmp = strcmp(upper,positive[k]);
    if (cmp == 0) { *val = 1; break; }
  }
  for (k=0; k<5; k++) {
    cmp = strcmp(upper,negative[k]);
    if (cmp == 0) { *val = 0; break; }
  }
  free(upper);
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

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
 **    filename:   geometry_object_parse.h
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

#ifndef __geometry_object_parse_h__
#define __geometry_object_parse_h__

#include "petsc.h"
#include "geometry_object.h"
#include "cJSON.h"

/* cJSON extensions */
void cJSON_FileView(const char filename[],cJSON **jf);
void cJSON_GetObjectValue_bool(cJSON *cj,const char name[],int *found,int *val);
void cJSON_GetObjectValue_int(cJSON *cj,const char name[],int *found,int *val);
void cJSON_GetObjectValue_intarray(cJSON *cj,const char name[],int *found,int *nv,int vals[]);
void cJSON_GetObjectValue_double(cJSON *cj,const char name[],int *found,double *val);
void cJSON_GetObjectValue_doublearray(cJSON *cj,const char name[],int *found,int *nv,double vals[]);
void cJSON_GetObjectValue_char(cJSON *cj,const char name[],int *found,char **val);
cJSON* cJSON_GetArrayItemRoot(cJSON *gobject);
cJSON* cJSON_GetArrayItemNext(cJSON *gobj_k);

PetscErrorCode GeometryObjectLoadJSON(const char filename[],PetscInt *n,GeometryObject **golist);
PetscErrorCode GeometryObjectListParseJSON(cJSON *jfile,int *_nlist,GeometryObject **_list);
PetscErrorCode GeometryObjectQueryFromJSONIsSetOperation(cJSON *obj,PetscBool *isset);
PetscErrorCode GeometryObjectPrimitiveLoadFromJSON(cJSON *obj,GeometryObject *g);

PetscErrorCode GeometryObjectParseDetermineAxis(const char axisname[],GeomRotateAxis *a);
PetscErrorCode GeometryObjectParseDetermineSign(const char name[],GeomSign *a);
PetscErrorCode GeometryObjectParseDetermineSetOperatorType(const char name[],GeomSetOperatorType *a);

#endif

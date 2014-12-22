/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012,
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics,
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      geometry_object_parse.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

#include "stdio.h"
#include "stdlib.h"
#include "petsc.h"
#include "geometry_object.h"
#include "cJSON.h"
#include "geometry_object_parse.h"


#undef __FUNCT__
#define __FUNCT__ "GeometryObjectLoadJSON"
PetscErrorCode GeometryObjectLoadJSON(const char filename[],PetscInt *n,GeometryObject **golist)
{
    PetscErrorCode ierr;
    cJSON          *jfile;
    int            nlist;
    GeometryObject *list;
    
    /* perform file validation */
    
    /* open file in a cJSON object */
    cJSON_FileView(filename,&jfile);
    if (!jfile) {
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Failed to open JSON file: %s",filename);
    }
    
    /* parse file */
    ierr = GeometryObjectListParseJSON(jfile,&nlist,&list);CHKERRQ(ierr);
 
    *n      = (PetscInt)nlist;
    *golist = list;

    /* filter list if required */
    {
        cJSON          *gofilterobject;
        int            k,cnt,nobj;
        GeometryObject *filter_list;
        
        gofilterobject = cJSON_GetObjectItem(jfile,"GeometryObjectListFilter");
        if (gofilterobject) {
            PetscBool isall;
            
            nobj = cJSON_GetArraySize(gofilterobject);
            
            filter_list = malloc(sizeof(GeometryObject)*nobj);
            
            PetscPrintf(PETSC_COMM_WORLD,"<GeometryObjectListFilter>\n");
            PetscPrintf(PETSC_COMM_WORLD,"...Filtering geometry objects with name\n");
            cnt = 0;
            for (k=0; k<nobj; k++) {
                char           *name1;
                GeometryObject go;
                
                name1 = cJSON_GetArrayItem(gofilterobject,k)->valuestring;
                ierr = GeometryObjectFindByName(list,name1,&go);CHKERRQ(ierr);
                
                PetscStrcmp(name1,"all",&isall);
                if (isall) {
                    PetscPrintf(PETSC_COMM_WORLD,"    \"all\" --> every object found in JSON file will be returned\n");
                    free(filter_list);
                    break;
                }
                
                if (!go) {
                    SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObjectListFilter: Failed to locate object with name \"%s\"\n",name1);
                }
                PetscPrintf(PETSC_COMM_WORLD,"    \"%s\"\n",name1);
                
                filter_list[cnt] = go;
                go->ref_cnt++;
                
                cnt++;
            }
            
            /* zip through and Destroy all entries in original list */
            if (!isall) {
            
                for (k=0; k<nlist; k++) {
                    ierr = GeometryObjectDestroy(&list[k]);CHKERRQ(ierr);
                }
                free(list);
                
                *n      = (PetscInt)nobj;
                *golist = filter_list;
            }
            
        }
    }
    
    cJSON_Delete(jfile);
    
    
    PetscFunctionReturn(0);
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


/*
 Geometry file must contain the following array/list
 
 GeometryObjectList: [
   {
     "name": "crust"
     "type": "GeomType_Sphere",
   }
 ]

 Optionally you may also specify

 GeometryObjectListFilter: [
   "name1", "name2", ...
 ]
 
*/
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectListParseJSON"
PetscErrorCode GeometryObjectListParseJSON(cJSON *jfile,int *_nlist,GeometryObject **_list)
{
    PetscErrorCode ierr;
    cJSON          *gobject = NULL;
    cJSON          *gobj_k;
    int            k,ngobj,found;
    int            np;
    GeometryObject *primitive_list,*tmp;
    PetscBool      *constructed_list;
    int            missing_objects,recursion;
    
    gobject = cJSON_GetObjectItem(jfile,"GeometryObjectList");
    if (!gobject) {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"<GeometryObjectParser> Failed to locate list named \"GeometryObjectList\" in JSON file");
    }
    
    np = 1;
    primitive_list = malloc(sizeof(GeometryObject)*np);
    primitive_list[0] = NULL;
    
    /* parse non-Set objects first */
    PetscPrintf(PETSC_COMM_WORLD,"<GeometryObjectList>\n");
    PetscPrintf(PETSC_COMM_WORLD,"...Processing primitive geometry objects\n");

    ngobj = cJSON_GetArraySize(gobject);
    PetscMalloc(sizeof(PetscBool)*ngobj,&constructed_list);
    
    for (k=0; k<ngobj; k++) {
        PetscBool      isset;
        GeometryObject G;

        constructed_list[k] = PETSC_FALSE;
        gobj_k = cJSON_GetArrayItem(gobject,k);

        ierr = GeometryObjectQueryFromJSONIsSetOperation(gobj_k,&isset);CHKERRQ(ierr);
        if (isset) { continue; }
        
        G = NULL;
        ierr = GeometryObjectPrimitiveLoadFromJSON(gobj_k,&G);CHKERRQ(ierr);
        
        /* re-size list and insert G */
        if (G != NULL) {
            PetscPrintf(PETSC_COMM_WORLD,"   Parsed: %s (%s)\n",G->name,GeomTypeNames[(int)(G->type)]);
            
            tmp = realloc(primitive_list,sizeof(GeometryObject)*(np+1));
            primitive_list = tmp;
            primitive_list[np-1] = G;
            primitive_list[np] = NULL;
            np++;

            constructed_list[k] = PETSC_TRUE;
        }
    }
    
    recursion = 0;
GeometryObjectSetOperators_Scan:
    PetscPrintf(PETSC_COMM_WORLD,"...Processing primitive set operator geometry objects [recursive level %d]\n",recursion);
    
    missing_objects = 0;
    for (k=0; k<ngobj; k++) {
        cJSON          *cjl_AB;
        char           *name,*name1,*name2;
        GeometryObject g1,g2;
        
        if (constructed_list[k]) { continue; }
        
        
        /* object is a set operation */
        gobj_k = cJSON_GetArrayItem(gobject,k);

        /* fetch name */
        cJSON_GetObjectValue_char(gobj_k,"name",&found,&name);
        if (!name) {
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject: Requires field \"name\"");
        }

        /* fetch fields */
        cjl_AB = cJSON_GetObjectItem(gobj_k,"fields");
        if (!cjl_AB) {
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.SetOperator: Requires field \"fields\"");
        }
        if (cJSON_GetArraySize(cjl_AB) != 2) {
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.SetOperator.fields: Requires array \"fields\" be of length 2");
        }

        name1 = cJSON_GetArrayItem(cjl_AB,0)->valuestring;
        name2 = cJSON_GetArrayItem(cjl_AB,1)->valuestring;

        PetscPrintf(PETSC_COMM_WORLD,"   Set operator \"%s\" depends on [ \"%s\" , \"%s\" ] \n",name,name1,name2);
        
        ierr = GeometryObjectFindByName(primitive_list,name1,&g1);CHKERRQ(ierr);
        ierr = GeometryObjectFindByName(primitive_list,name2,&g2);CHKERRQ(ierr);
        
        
        /* generate object */
        if ( (g1 != NULL) && (g2 != NULL) ) {
            char                *setopname;
            GeometryObject      G;
            double              dummy[] = { 0.0, 0.0, 0.0 };
            GeomSetOperatorType setoptype;
            
            /* fetch operator */
            cJSON_GetObjectValue_char(gobj_k,"operator",&found,&setopname);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.SetOperator: Requires field \"operator\"");
            }
            ierr = GeometryObjectParseDetermineSetOperatorType(setopname,&setoptype);CHKERRQ(ierr);
            
            PetscPrintf(PETSC_COMM_WORLD,"   ... generatoring geom object\n");
            ierr = GeometryObjectCreate(name,&G);CHKERRQ(ierr);
            ierr = GeometryObjectParseDetermineSetOperatorType(setopname,&setoptype);CHKERRQ(ierr);
            ierr = GeometryObjectSetType_SetOperation(G,setoptype,dummy,g1,g2);CHKERRQ(ierr);
            
            tmp = realloc(primitive_list,sizeof(GeometryObject)*(np+1));
            primitive_list = tmp;
            primitive_list[np-1] = G;
            primitive_list[np] = NULL;
            np++;
            
            constructed_list[k] = PETSC_TRUE;
        } else {
            /* count number of missing object */
            if (g1) {PetscPrintf(PETSC_COMM_WORLD,"      + \"%s\" Found\n",name1);}
            else    {PetscPrintf(PETSC_COMM_WORLD,"      - \"%s\" Not found\n",name1);}
            if (g2) {PetscPrintf(PETSC_COMM_WORLD,"      + \"%s\" Found\n",name2);}
            else    {PetscPrintf(PETSC_COMM_WORLD,"      - \"%s\" Not found\n",name1);}
            missing_objects++;
        }
    }
    if (missing_objects != 0) {
        if (recursion < 20) {
            PetscPrintf(PETSC_COMM_WORLD,"...Missing #%d dependent objects\n",missing_objects);
            recursion++;
            goto GeometryObjectSetOperators_Scan;
        } else {
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"JSON.GeometryObject parsing required more than 20 levels of recursion - aborting\n");
        }
    }
    
    PetscFree(constructed_list);

    /* return result */
    *_nlist = np-1;
    *_list = primitive_list;
    
    PetscFunctionReturn(0);
}



/* ------------------------------------------------------------------------------------------- */
/* Implementations */

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectQueryFromJSONIsSetOperation"
PetscErrorCode GeometryObjectQueryFromJSONIsSetOperation(cJSON *obj,PetscBool *isset)
{
    int       found;
    char      *type;
    PetscBool same;
    
    cJSON_GetObjectValue_char(obj,"type",&found,&type);
    if (found == cJSON_False) {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"JSON.GeometryObject requires field \"type\"");
    }
    
    PetscStrcmp(type,GeomTypeNames[(int)GeomType_SetOperation],&same);

    *isset = PETSC_FALSE;
    if (same) {
        *isset = PETSC_TRUE;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectParseDetermineAxis"
PetscErrorCode GeometryObjectParseDetermineAxis(const char axisname[],GeomRotateAxis *a)
{
    PetscBool same;
    
    *a = ROTATE_AXIS_UNDEFINED;
    
    same = PETSC_FALSE; PetscStrcmp(axisname,"x",&same); if (same) { *a = ROTATE_AXIS_X; }
    same = PETSC_FALSE; PetscStrcmp(axisname,"X",&same); if (same) { *a = ROTATE_AXIS_X; }
    same = PETSC_FALSE; PetscStrcmp(axisname,GeomRotateAxisNames[0],&same); if (same) { *a = ROTATE_AXIS_X; }

    same = PETSC_FALSE; PetscStrcmp(axisname,"y",&same); if (same) { *a = ROTATE_AXIS_Y; }
    same = PETSC_FALSE; PetscStrcmp(axisname,"Y",&same); if (same) { *a = ROTATE_AXIS_Y; }
    same = PETSC_FALSE; PetscStrcmp(axisname,GeomRotateAxisNames[1],&same); if (same) { *a = ROTATE_AXIS_Y; }

    same = PETSC_FALSE; PetscStrcmp(axisname,"z",&same); if (same) { *a = ROTATE_AXIS_Z; }
    same = PETSC_FALSE; PetscStrcmp(axisname,"Z",&same); if (same) { *a = ROTATE_AXIS_Z; }
    same = PETSC_FALSE; PetscStrcmp(axisname,GeomRotateAxisNames[2],&same); if (same) { *a = ROTATE_AXIS_Z; }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectParseDetermineSign"
PetscErrorCode GeometryObjectParseDetermineSign(const char name[],GeomSign *a)
{
    PetscBool same;
    
    *a = SIGN_UNDEFINED;
    
    same = PETSC_FALSE; PetscStrcmp(name,"+",&same); if (same) { *a = SIGN_POSITIVE; }
    same = PETSC_FALSE; PetscStrcmp(name,"positive",&same); if (same) { *a = SIGN_POSITIVE; }
    same = PETSC_FALSE; PetscStrcmp(name,GeomSignNames[0],&same); if (same) { *a = SIGN_POSITIVE; }
    
    same = PETSC_FALSE; PetscStrcmp(name,"-",&same); if (same) { *a = SIGN_NEGATIVE; }
    same = PETSC_FALSE; PetscStrcmp(name,"negative",&same); if (same) { *a = SIGN_NEGATIVE; }
    same = PETSC_FALSE; PetscStrcmp(name,GeomSignNames[1],&same); if (same) { *a = SIGN_NEGATIVE; }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectParseDetermineSetOperatorType"
PetscErrorCode GeometryObjectParseDetermineSetOperatorType(const char name[],GeomSetOperatorType *a)
{
    PetscBool same;
    
    *a = SIGN_UNDEFINED;
    
    same = PETSC_FALSE; PetscStrcmp(name,"+",&same); if (same) { *a = GeomSet_Union; }
    same = PETSC_FALSE; PetscStrcmp(name,"positive",&same); if (same) { *a = GeomSet_Union; }
    same = PETSC_FALSE; PetscStrcmp(name,GeomSignNames[0],&same); if (same) { *a = GeomSet_Union; }
    
    same = PETSC_FALSE; PetscStrcmp(name,"-",&same); if (same) { *a = SIGN_NEGATIVE; }
    same = PETSC_FALSE; PetscStrcmp(name,"negative",&same); if (same) { *a = SIGN_NEGATIVE; }
    same = PETSC_FALSE; PetscStrcmp(name,GeomSignNames[1],&same); if (same) { *a = SIGN_NEGATIVE; }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPrimitiveLoadFromJSON"
PetscErrorCode GeometryObjectPrimitiveLoadFromJSON(cJSON *obj,GeometryObject *g)
{
    char *type;
    char *name;
    GeomType gtype;
    int i,found;
    GeometryObject go;
    PetscErrorCode ierr;
    
    /* fetch type */
    cJSON_GetObjectValue_char(obj,"type",&found,&type);
    if (!type) {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject: Requires field \"type\"");
    }

    gtype = GeomType_NULL;
    for (i=0; i<((int)GeomType_NULL); i++) {
        PetscBool same = PETSC_FALSE;
        
        PetscStrcmp(type,GeomTypeNames[i],&same);
        if (same) {
            gtype = (GeomType)i;
            break;
        }
    }
    if (gtype == GeomType_NULL) {
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject: Failed to map \"type\" field (%s) to a valid GeomType",type);
    }
    
    if (gtype == GeomType_SetOperation) {
        printf("NNNNN\n");
        *g = NULL;
        PetscFunctionReturn(0);
    }
    
    /* fetch name */
    cJSON_GetObjectValue_char(obj,"name",&found,&name);
    if (!name) {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject: Requires field \"name\"");
    }

    ierr = GeometryObjectCreate(name,&go);CHKERRQ(ierr);
    
    switch (gtype) {
            
        case GeomType_Box:
        {
            cJSON  *cjl_c,*cjl_l;
            double cx[3],lx[3];
            int    nv,cref;
            
            cJSON_GetObjectValue_bool(obj,"corner_ref",&found,&cref);
            
            cjl_c = cJSON_GetObjectItem(obj,"centroid");
            if (!cjl_c) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Box: Requires field \"centroid\"");
            }
            if (cJSON_GetArraySize(cjl_c) != 3) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Box.centroid: Requires array \"centroid\" be of length 3");
            }
            cJSON_GetObjectValue_doublearray(obj,"centroid",&found,&nv,cx);
            
            cjl_l = cJSON_GetObjectItem(obj,"length");
            if (!cjl_l) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject:Box: Requires field \"length\"");
            }
            if (cJSON_GetArraySize(cjl_l) != 3) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Box.length: Requires array \"length\" be of length 3");
            }
            cJSON_GetObjectValue_doublearray(obj,"length",&found,&nv,lx);
            
            /* construct */
            if (cref == cJSON_True) {
                ierr = GeometryObjectSetType_BoxCornerReference(go,cx,lx);CHKERRQ(ierr);
            } else {
                ierr = GeometryObjectSetType_Box(go,cx,lx);CHKERRQ(ierr);
            }
        }
            break;

        case GeomType_Cylinder:
        {
            cJSON  *cjl_c;
            double cx[3],length,rad;
            int    nv;
            char   *axisname;
            GeomRotateAxis axis;
            
            cjl_c = cJSON_GetObjectItem(obj,"centroid");
            if (!cjl_c) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Cylinder: Requires field \"centroid\"");
            }
            if (cJSON_GetArraySize(cjl_c) != 3) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Cylinder.centroid: Requires array \"centroid\" be of length 3");
            }
            cJSON_GetObjectValue_doublearray(obj,"centroid",&found,&nv,cx);
        
            cJSON_GetObjectValue_double(obj,"radius",&found,&rad);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Cylinder: Requires field \"radius\"");
            }

            cJSON_GetObjectValue_double(obj,"length",&found,&length);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Cylinder: Requires field \"length\"");
            }

            cJSON_GetObjectValue_char(obj,"axis",&found,&axisname);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Cylinder: Requires field \"axis\"");
            }
            ierr = GeometryObjectParseDetermineAxis(axisname,&axis);CHKERRQ(ierr);
            
            /* construct */
            ierr = GeometryObjectSetType_Cylinder(go,cx,rad,length,axis);CHKERRQ(ierr);
        }
            break;

        case GeomType_Sphere:
        {
            cJSON  *cjl_c;
            double cx[3],rad;
            int    nv;

            cJSON_GetObjectValue_double(obj,"radius",&found,&rad);

            cjl_c = cJSON_GetObjectItem(obj,"origin");
            if (!cjl_c) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Sphere: Requires field \"origin\"");
            }
            if (cJSON_GetArraySize(cjl_c) != 3) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Sphere.origin: Requires array \"origin\" be of length 3");
            }
            cJSON_GetObjectValue_doublearray(obj,"origin",&found,&nv,cx);
            
            cJSON_GetObjectValue_double(obj,"radius",&found,&rad);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Sphere: Requires field \"radius\"");
            }
            
            /* construct */
            ierr = GeometryObjectSetType_Sphere(go,cx,rad);CHKERRQ(ierr);
        }
            break;

        case GeomType_EllipticCylinder:
        {
            cJSON  *cjl_c;
            double cx[3],length,radA,radB;
            int    nv;
            char   *axisname;
            GeomRotateAxis axis;
            
            cjl_c = cJSON_GetObjectItem(obj,"centroid");
            if (!cjl_c) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.EllipticCylinder: Requires field \"centroid\"");
            }
            if (cJSON_GetArraySize(cjl_c) != 3) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.EllipticCylinder.centroid: Requires array \"centroid\" be of length 3");
            }
            cJSON_GetObjectValue_doublearray(obj,"centroid",&found,&nv,cx);
            
            cJSON_GetObjectValue_double(obj,"radiusA",&found,&radA);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.EllipticCylinder: Requires field \"radiusA\"");
            }

            cJSON_GetObjectValue_double(obj,"radiusB",&found,&radB);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.EllipticCylinder: Requires field \"radiusB\"");
            }
            
            cJSON_GetObjectValue_double(obj,"length",&found,&length);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.EllipticCylinder: Requires field \"length\"");
            }
            
            cJSON_GetObjectValue_char(obj,"axis",&found,&axisname);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.EllipticCylinder: Requires field \"axis\"");
            }
            ierr = GeometryObjectParseDetermineAxis(axisname,&axis);CHKERRQ(ierr);
            
            /* construct */
            ierr = GeometryObjectSetType_EllipticCylinder(go,cx,radA,radB,length,axis);CHKERRQ(ierr);
        }
            break;

        case GeomType_Ellipsoid:
        {
            cJSON  *cjl_c;
            double cx[3],radA,radB,radC;
            int    nv;
            
            cjl_c = cJSON_GetObjectItem(obj,"centroid");
            if (!cjl_c) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Ellipsoid: Requires field \"centroid\"");
            }
            if (cJSON_GetArraySize(cjl_c) != 3) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Ellipsoid.centroid: Requires array \"centroid\" be of length 3");
            }
            cJSON_GetObjectValue_doublearray(obj,"centroid",&found,&nv,cx);
            
            cJSON_GetObjectValue_double(obj,"radiusA",&found,&radA);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Ellipsoid: Requires field \"radiusA\"");
            }
            
            cJSON_GetObjectValue_double(obj,"radiusB",&found,&radB);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Ellipsoid: Requires field \"radiusB\"");
            }

            cJSON_GetObjectValue_double(obj,"radiusC",&found,&radC);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.Ellipsoid: Requires field \"radiusC\"");
            }

            /* construct */
            ierr = GeometryObjectSetType_Ellipsoid(go,cx,radA,radB,radC);CHKERRQ(ierr);
        }
            break;

        case GeomType_InfLayer:
        {
            cJSON  *cjl_c;
            double cx[3],thickness;
            int nv;
            GeomRotateAxis axis;
            char           *axisname;

            cjl_c = cJSON_GetObjectItem(obj,"centroid");
            if (!cjl_c) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.InfLayer: Requires field \"centroid\"");
            }
            if (cJSON_GetArraySize(cjl_c) != 3) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.InfLayer.centroid: Requires array \"centroid\" be of length 3");
            }
            cJSON_GetObjectValue_doublearray(obj,"centroid",&found,&nv,cx);
            
            cJSON_GetObjectValue_double(obj,"thickness",&found,&thickness);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.InfLayer: Requires field \"thickness\"");
            }
            
            cJSON_GetObjectValue_char(obj,"axis",&found,&axisname);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.InfLayer: Requires field \"axis\"");
            }
            ierr = GeometryObjectParseDetermineAxis(axisname,&axis);CHKERRQ(ierr);

            /* construct */
            ierr = GeometryObjectSetType_InfLayer(go,cx,thickness,axis);CHKERRQ(ierr);
        }
            break;
            
        case GeomType_SetOperation:
        {
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"JSON.GeometryObject.Set: Should be generated by other functionality");
        }
            break;

        case GeomType_HalfSpace:
        {
            cJSON  *cjl_c;
            double cx[3];
            int nv;
            GeomSign sign;
            GeomRotateAxis axis;
            char   *signname,*axisname;
            
            cjl_c = cJSON_GetObjectItem(obj,"centroid");
            if (!cjl_c) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.HalfSpace: Requires field \"centroid\"");
            }
            if (cJSON_GetArraySize(cjl_c) != 3) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.HalfSpace.centroid: Requires array \"centroid\" be of length 3");
            }
            cJSON_GetObjectValue_doublearray(obj,"centroid",&found,&nv,cx);
            
            cJSON_GetObjectValue_char(obj,"axis",&found,&axisname);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.HalfSpace: Requires field \"axis\"");
            }
            ierr = GeometryObjectParseDetermineAxis(axisname,&axis);CHKERRQ(ierr);

            cJSON_GetObjectValue_char(obj,"sign",&found,&signname);
            if (found == cJSON_False) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"JSON.GeometryObject.HalfSpace: Requires field \"sign\"");
            }
            ierr = GeometryObjectParseDetermineSign(signname,&sign);CHKERRQ(ierr);
            
            /* construct */
            ierr = GeometryObjectSetType_HalfSpace(go,cx,sign,axis);CHKERRQ(ierr);
        }
            break;
    }

    *g = go;
    
    PetscFunctionReturn(0);
}

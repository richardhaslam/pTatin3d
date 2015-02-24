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
 **    filename:   cartgrid.c
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

#include <ctype.h>
#include <petsc.h>
#include "ptatin3d.h"
#include "cartgrid.h"


const char *PMapTypeNames[] = { "udef", "inmem", "outofcore", 0 };
const char *PMapDataTypeNames[] = { "udef", "int", "long", "float", "double", "short", "char", 0 };


PetscErrorCode PMapGetIndex_InMem(PMap map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index);
PetscErrorCode PMapGetValue_InMem(PMap map,PetscReal xp[],void *value,PMapValueFound *found);

PetscErrorCode PMapGetValue_OutOfCore(PMap map,PetscReal xp[],void *value,PMapValueFound *found);



#undef __FUNCT__
#define __FUNCT__ "PMapCreate"
PetscErrorCode PMapCreate(PMap *map)
{
	PMap p;
    
    PetscFunctionBegin;
	PetscNew(&p);
    
    p->type = PMAP_INMEM;
    p->getindex = PMapGetIndex_InMem;
    p->getvalue = PMapGetValue_InMem;
    p->destroy = NULL;
    
	*map = p;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapDestroy"
PetscErrorCode PMapDestroy(PMap *map)
{
	PMap p;
    
    PetscFunctionBegin;
	if (!map) { PetscFunctionReturn(0); }
	p = *map;
	
	if (p->data) {
		PetscFree(p->data);
		p->data = NULL;
	}

    if (p->type == PMAP_OUTOFCORE) {
        if (p->data_fp) {
            fclose(p->data_fp);
            p->data_fp = NULL;
        }
    }
    
    PetscFree(p);
    
	*map = NULL;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetType"
PetscErrorCode PMapSetType(PMap map,PMapType t)
{
    PetscFunctionBegin;
    map->type = t;
    
    switch (t) {
            
        case PMAP_INMEM:
            map->getindex = PMapGetIndex_InMem;
            map->getvalue = PMapGetValue_InMem;
            map->destroy = NULL;
            break;

        case PMAP_OUTOFCORE:
            map->getindex = PMapGetIndex_InMem;
            map->getvalue = PMapGetValue_OutOfCore;
            map->destroy = NULL;
            break;
            
        default:
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetDim"
PetscErrorCode PMapSetDim(PMap map,PetscInt dim)
{
    PetscFunctionBegin;
    map->dim = dim;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetSizes"
PetscErrorCode PMapSetSizes(PMap map,PetscInt m,PetscInt n,PetscInt p)
{
    PetscFunctionBegin;
    map->mx = m;
    map->my = n;
    map->mz = p;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetDataType"
PetscErrorCode PMapSetDataType(PMap map,PMapDataType type)
{
    PetscFunctionBegin;
    map->data_type = type;
    switch (type) {
        case PMAP_DTYPE_UDEF:
            map->bytes = -1;
            break;

        case PMAP_INT:
            map->bytes = sizeof(int);
            break;
            
        case PMAP_LONG:
            map->bytes = sizeof(long int);
            break;

        case PMAP_FLOAT:
            map->bytes = sizeof(float);
            break;

        case PMAP_DOUBLE:
            map->bytes = sizeof(double);
            break;

        case PMAP_SHORT:
            map->bytes = sizeof(short);
            break;
        
        case PMAP_CHAR:
            map->bytes = sizeof(char);
            break;

        default:
            map->bytes = -1;
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetDomain"
PetscErrorCode PMapSetDomain(PMap map,PetscReal xr[],PetscReal yr[],PetscReal zr[])
{
    PetscFunctionBegin;
    if (xr) { map->range_x[0] = xr[0]; map->range_x[1] = xr[1]; }
    if (yr) { map->range_y[0] = yr[0]; map->range_y[1] = yr[1]; }
    if (zr) { map->range_z[0] = zr[0]; map->range_z[1] = zr[1]; }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetFilename"
PetscErrorCode PMapSetFilename(PMap map,const char fname[])
{
    PetscFunctionBegin;
    PetscSNPrintf(map->metadatafile_name,PETSC_MAX_PATH_LEN-1,"%s",fname);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetDataFilename"
PetscErrorCode PMapSetDataFilename(PMap map,const char fname[])
{
    PetscFunctionBegin;
    PetscSNPrintf(map->datafile_name,PETSC_MAX_PATH_LEN-1,"%s",fname);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetUp_InMem"
PetscErrorCode PMapSetUp_InMem(PMap map)
{
    FILE *fp_data = NULL;
    
    PetscFunctionBegin;
    map->start[0] = 0;
    map->start[1] = 0;
    map->start[2] = 0;

    map->end[0] = map->mx;
    map->end[1] = map->my;
    map->end[2] = map->mz;

	/* allocate data */
	PetscMalloc(map->bytes* map->mx * map->my * map->mz,&map->data);
    PetscMemzero(map->data,map->bytes* map->mx * map->my * map->mz);
    
    fp_data = fopen(map->datafile_name,"r");
    if (!fp_data) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open %s",map->datafile_name);
    fread(map->data,map->bytes,map->mx * map->my * map->mz,fp_data);
    fclose(fp_data);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetUp_OutOfCore"
PetscErrorCode PMapSetUp_OutOfCore(PMap map)
{
    FILE *fp_data = NULL;

    PetscFunctionBegin;
    map->start[0] = 0;
    map->start[1] = 0;
    map->start[2] = 0;
    
    map->end[0] = map->mx;
    map->end[1] = map->my;
    map->end[2] = map->mz;

    fp_data = fopen(map->datafile_name,"r");
    if (!fp_data) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open %s",map->datafile_name);
    
    map->data_fp = fp_data;
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapGetIndex_InMem2d"
PetscErrorCode PMapGetIndex_InMem2d(PMap map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index)
{
    PetscFunctionBegin;
	if (i < map->start[0]) { *index = -1; PetscFunctionReturn(0); }
	if (j < map->start[1]) { *index = -1; PetscFunctionReturn(0); }
	if (i >= map->end[0])  { *index = -1; PetscFunctionReturn(0); }
	if (j >= map->end[1])  { *index = -1; PetscFunctionReturn(0); }
    
	*index = i + j * map->mx;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapGetIndex_InMem3d"
PetscErrorCode PMapGetIndex_InMem3d(PMap map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index)
{
    PetscFunctionBegin;
	if (i < map->start[0]) { *index = -1; PetscFunctionReturn(0); }
	if (j < map->start[1]) { *index = -1; PetscFunctionReturn(0); }
	if (k < map->start[2]) { *index = -1; PetscFunctionReturn(0); }
	if (i >= map->end[0])  { *index = -1; PetscFunctionReturn(0); }
	if (j >= map->end[1])  { *index = -1; PetscFunctionReturn(0); }
	if (k >= map->end[2])  { *index = -1; PetscFunctionReturn(0); }
    
	*index = i + j * map->mx + k * map->mx * map->my;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapGetIndex_InMem"
PetscErrorCode PMapGetIndex_InMem(PMap map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    switch (map->dim) {
        case 2:
            ierr = PMapGetIndex_InMem2d(map,i,j,-1,index);CHKERRQ(ierr);
            break;
        case 3:
            ierr = PMapGetIndex_InMem3d(map,i,j,k,index);CHKERRQ(ierr);
            break;
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapGetValue_InMem2d"
PetscErrorCode PMapGetValue_InMem2d(PMap map,PetscReal xp[],void *value,PMapValueFound *found)
{
	PetscInt i,j,index;
    void     *value_i;
	
    PetscFunctionBegin;
	if (xp[0] < map->range_x[0]) { PetscFunctionReturn(0); }
	if (xp[0] > map->range_x[1]) { PetscFunctionReturn(0); }
    
	if (xp[1] < map->range_y[0]) { PetscFunctionReturn(0); }
	if (xp[1] > map->range_y[1]) { PetscFunctionReturn(0); }
	
	i = (xp[0] - map->range_x[0])/map->dx;
	j = (xp[1] - map->range_y[0])/map->dy;
	if (i == map->mx) { i--; }
	if (j == map->my) { j--; }
	
	PMapGetIndex_InMem2d(map,i,j,-1,&index);
	if (index < 0) { PetscFunctionReturn(0); }
    
    *found = PMAP_INSIDE;
    
    value_i = (void*)( (char*)map->data + index * map->bytes );
    PetscMemcpy(value,value_i,map->bytes);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapGetValue_InMem3d"
PetscErrorCode PMapGetValue_InMem3d(PMap map,PetscReal xp[],void *value,PMapValueFound *found)
{
	PetscInt i,j,k,index;
    void     *value_i;
	
    PetscFunctionBegin;
	if (xp[0] < map->range_x[0]) { PetscFunctionReturn(0); }
	if (xp[0] > map->range_x[1]) { PetscFunctionReturn(0); }

	if (xp[1] < map->range_y[0]) { PetscFunctionReturn(0); }
	if (xp[1] > map->range_y[1]) { PetscFunctionReturn(0); }

	if (xp[2] < map->range_z[0]) { PetscFunctionReturn(0); }
	if (xp[2] > map->range_z[1]) { PetscFunctionReturn(0); }
	
	i = (xp[0] - map->range_x[0])/map->dx;
	j = (xp[1] - map->range_y[0])/map->dy;
	k = (xp[2] - map->range_z[0])/map->dz;
	if (i == map->mx) { i--; }
	if (j == map->my) { j--; }
	if (k == map->mz) { k--; }
	
	PMapGetIndex_InMem3d(map,i,j,k,&index);
	if (index < 0) { PetscFunctionReturn(0); }
    
    *found = PMAP_INSIDE;
    
    value_i = (void*)( (char*)map->data + index * map->bytes );
    PetscMemcpy(value,value_i,map->bytes);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapGetValue_InMem"
PetscErrorCode PMapGetValue_InMem(PMap map,PetscReal xp[],void *value,PMapValueFound *found)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    switch (map->dim) {
        case 2:
            ierr = PMapGetValue_InMem2d(map,xp,value,found);CHKERRQ(ierr);
            break;
        case 3:
            ierr = PMapGetValue_InMem3d(map,xp,value,found);CHKERRQ(ierr);
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapGetValue_OutOfCore"
PetscErrorCode PMapGetValue_OutOfCore(PMap map,PetscReal xp[],void *value,PMapValueFound *found)
{
	PetscInt i,j,k,index;
    long     offset;
	
    PetscFunctionBegin;
	if (xp[0] < map->range_x[0]) { PetscFunctionReturn(0); }
	if (xp[0] > map->range_x[1]) { PetscFunctionReturn(0); }
    
	if (xp[1] < map->range_y[0]) { PetscFunctionReturn(0); }
	if (xp[1] > map->range_y[1]) { PetscFunctionReturn(0); }
    
	i = (xp[0] - map->range_x[0])/map->dx;
	j = (xp[1] - map->range_y[0])/map->dy;
	if (i == map->mx) { i--; }
	if (j == map->my) { j--; }
	
    if (map->dim == 3) {
        if (xp[2] < map->range_z[0]) { PetscFunctionReturn(0); }
        if (xp[2] > map->range_z[1]) { PetscFunctionReturn(0); }
        k = (xp[2] - map->range_z[0])/map->dz;
        if (k == map->mz) { k--; }
        
        PMapGetIndex_InMem3d(map,i,j,k,&index);
    } else {
        PMapGetIndex_InMem2d(map,i,j,-1,&index);
    }
	if (index < 0) { PetscFunctionReturn(0); }
    
    *found = PMAP_INSIDE;

    /* reset position */
    fseek(map->data_fp,0,SEEK_SET);
    /* displace position */
    offset = (long)(index * map->bytes);
    fseek(map->data_fp,offset,SEEK_CUR);
    
    /* fread */
    fread(value,map->bytes,1,map->data_fp);
 
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapGetValue"
PetscErrorCode PMapGetValue(PMap map,PetscReal xp[],void *value,PMapValueFound *found)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
	(*found) = PMAP_OUTSIDE;
    ierr = map->getvalue(map,xp,value,found);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "trimright"
PetscErrorCode trimright(char str[])
{
    size_t len;
    int    c,flg;
    
    PetscFunctionBegin;
    len = strlen(str);
    for (c=len-1; c>=0; c--) {
        flg = isspace((int)str[c]);
        if (flg != 0) {
            str[c] = '\0';
        } else {
            break;
        }
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "trimnewline"
PetscErrorCode trimnewline(char str[])
{
    size_t len;
    PetscFunctionBegin;
    len = strlen(str);
    if (str[len-1] == '\n') {
        str[len-1] = '\0';
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapSetUp"
PetscErrorCode PMapSetUp(PMap map)
{
    FILE           *fp;
    int            t,dim,mx,my,mz;
    double         xr[2],yr[2],zr[2],dx,dy,dz;
    char           str[60];
    int            matched;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
	/* open file to parse */
    fp = fopen(map->metadatafile_name,"r");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open %s",map->metadatafile_name);
    
    /* pmap type */
    if( fgets (str,60,fp) != NULL) {
        PMapType type;
        
        trimnewline(str);
        trimright(str);
        type = PMAP_TYPE_UDEF;
        for (t=1; t<4; t++) {
            matched = strcmp(str,PMapTypeNames[t]);
            if (matched == 0) {
                type = (PMapType)t;
                break;
            }
        }
        ierr = PMapSetType(map,type);CHKERRQ(ierr);
        printf("pmaptype: %s\n",PMapTypeNames[(int)map->type]);
    }
    
    /* data_type */
    if( fgets (str,60,fp) != NULL) {
        PMapDataType data_type;
        
        trimnewline(str);
        trimright(str);
        data_type = PMAP_DTYPE_UDEF;
        for (t=1; t<5; t++) {
            matched = strcmp(str,PMapDataTypeNames[t]);
            if (matched == 0) {
                data_type = (PMapDataType)t;
                break;
            }
        }
        ierr = PMapSetDataType(map,data_type);CHKERRQ(ierr);
        printf("datatype: %s\n",PMapDataTypeNames[(int)map->data_type]);
    }
    
    /* datafile_name */
    if( fgets (str,60,fp) != NULL) {
        trimnewline(str);
        trimright(str);
        printf("datafile: %s\n",str);
        ierr = PMapSetDataFilename(map,str);CHKERRQ(ierr);
    }
    
    /* read header: dim */
    fscanf(fp,"%d\n",&dim);
    map->dim = (PetscInt)dim;
    
    if (dim == 2) {
        /* read header: mx,my */
        fscanf(fp,"%d %d\n",&mx,&my);
        mz = 1;
        /* read header: x0,x1,y0,y1 */
        fscanf(fp,"%lf %lf %lf %lf\n",&xr[0],&xr[1],&yr[0],&yr[1]);
        zr[0] = -1.0e32;
        zr[1] =  1.032;
        printf("[%d x %d] --> [%1.4e,%1.4e] x [%1.4e,%1.4e]\n",mx,my,xr[0],xr[1],yr[0],yr[1]);
        
        dx = (xr[1] - xr[0])/(double)(mx);
        dy = (yr[1] - yr[0])/(double)(my);
        dz = 0.0;
    } else if (dim == 3) {
        /* read header: mx,my,mz */
        fscanf(fp,"%d %d %d\n",&mx,&my,&mz);
        /* read header: x0,x1,y0,y1,z0,z1 */
        fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",&xr[0],&xr[1],&yr[0],&yr[1],&zr[0],&zr[1]);
        printf("[%d x %d x %d] --> [%1.2e,%1.2e] x [%1.2e,%1.2e] x [%1.2e,%1.2e]\n",mx,my,mz,xr[0],xr[1],yr[0],yr[1],zr[0],zr[1]);

        dx = (xr[1] - xr[0])/(double)(mx);
        dy = (yr[1] - yr[0])/(double)(my);
        dz = (zr[1] - zr[0])/(double)(mz);
    } else {
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"dim must be specified");
    }
    
    map->mx = (PetscInt)mx;
    map->my = (PetscInt)my;
    map->mz = (PetscInt)mz;
    
    map->dx = (PetscReal)dx;
    map->dy = (PetscReal)dy;
    map->dz = (PetscReal)dz;
    
    map->range_x[0] = (PetscReal)xr[0];  map->range_x[1] = (PetscReal)xr[1];
    map->range_y[0] = (PetscReal)yr[0];  map->range_y[1] = (PetscReal)yr[1];
    map->range_z[0] = (PetscReal)zr[0];  map->range_z[1] = (PetscReal)zr[1];
    
    fclose(fp);
    
    switch (map->type) {
        case PMAP_INMEM:
            ierr = PMapSetUp_InMem(map);CHKERRQ(ierr);
            break;
        case PMAP_OUTOFCORE:
            ierr = PMapSetUp_OutOfCore(map);CHKERRQ(ierr);
            break;
        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Must specify a valid CartGrid type");
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapViewMetaData"
PetscErrorCode PMapViewMetaData(PMap map)
{
    FILE *fp;
    
    PetscFunctionBegin;
    fp = fopen(map->metadatafile_name,"w");
    if (!fp)  SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open %s",map->metadatafile_name);
    
    fprintf(fp,"%s\n",PMapTypeNames[(int)map->type]);
    fprintf(fp,"%s\n",PMapDataTypeNames[(int)map->data_type]);
    fprintf(fp,"%s\n",map->datafile_name);

    fprintf(fp,"%d\n",map->dim);
    if (map->dim == 2) {
        fprintf(fp,"%d %d\n",map->mx,map->my);
        fprintf(fp,"%1.4e %1.4e %1.4e %1.4e\n",
                map->range_x[0],map->range_x[1],
                map->range_y[0],map->range_y[1]);
    } else {
        fprintf(fp,"%d %d %d\n",map->mx,map->my,map->mz);
        fprintf(fp,"%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n",
                map->range_x[0],map->range_x[1],
                map->range_y[0],map->range_y[1],
                map->range_z[0],map->range_z[1]);
    }
    fclose(fp);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PMapViewPV"
PetscErrorCode PMapViewPV(PMap map,const char filename[])
{
	FILE     *fp = NULL;
	PetscInt i,j,k;
	
    PetscFunctionBegin;
	/* open file to parse */
    fp = fopen(filename,"w");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open %s",filename);

    fprintf(fp,"<?xml version=\"1.0\"?>\n");
    fprintf(fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    if (map->dim == 2) {
        fprintf(fp,"  <ImageData WholeExtent=\"%d %d %d %d %d %d\"\n",0,map->mx,0,map->my,0,1);
        fprintf(fp,"   Origin=\"%1.4e %1.4e 0.0\"\n",map->range_x[0],map->range_y[0]);
        fprintf(fp,"   Spacing=\"%1.2e %1.2e %1.2e\">\n",map->dx,map->dy,1.0e-3);
        fprintf(fp,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",0,map->mx,0,map->my,0,1);
    } else {
        fprintf(fp,"  <ImageData WholeExtent=\"%d %d %d %d %d %d\"\n",0,map->mx,0,map->my,0,map->mz);
        fprintf(fp,"   Origin=\"%1.4e %1.4e %1.4e\"\n",map->range_x[0],map->range_y[0],map->range_z[0]);
        fprintf(fp,"   Spacing=\"%1.2e %1.2e %1.2e\">\n",map->dx,map->dy,map->dz);
        fprintf(fp,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",0,map->mx,0,map->my,0,map->mz);
    }

    fprintf(fp,"      <CellData>\n");
    switch (map->data_type) {
        case PMAP_INT:
            fprintf(fp,"      <DataArray type=\"Int32\" Name=\"values\" format=\"ascii\">\n");
            break;
        case PMAP_LONG:
            fprintf(fp,"      <DataArray type=\"Int64\" Name=\"values\" format=\"ascii\">\n");
            break;
        case PMAP_FLOAT:
            fprintf(fp,"      <DataArray type=\"Float32\" Name=\"values\" format=\"ascii\">\n");
            break;
        case PMAP_DOUBLE:
            fprintf(fp,"      <DataArray type=\"Float64\" Name=\"values\" format=\"ascii\">\n");
            break;
        case PMAP_SHORT:
            fprintf(fp,"      <DataArray type=\"Int16\" Name=\"values\" format=\"ascii\">\n");
            break;
        case PMAP_CHAR:
            fprintf(fp,"      <DataArray type=\"Int8\" Name=\"values\" format=\"ascii\">\n");
            break;

        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"A valid VTK data type must be specified");
            break;
    }

    switch (map->type) {
        case PMAP_INMEM:
            
            for (k=0; k<map->mz; k++) {
                for (j=0; j<map->my; j++) {
                    for (i=0; i<map->mx; i++) {
                        PetscInt index;
                        void     *data_i;
                        
                        index = i + j * map->mx + k * map->mx * map->my;
                        
                        data_i = (void*)( (char*)map->data + index * map->bytes );
                        
                        switch (map->data_type) {
                            case PMAP_INT:
                                fprintf(fp,"%d ",*((int*)data_i));
                                break;
                            case PMAP_LONG:
                                fprintf(fp,"%ld ",*((long int*)data_i));
                                break;
                            case PMAP_FLOAT:
                                fprintf(fp,"%1.4e ",*((float*)data_i));
                                break;
                            case PMAP_DOUBLE:
                                fprintf(fp,"%1.4e ",*((double*)data_i));
                                break;
                            case PMAP_SHORT:
                                fprintf(fp,"%hd ",*((short*)data_i));
                                break;
                            case PMAP_CHAR:
                                fprintf(fp,"%hd ",(short)*((char*)data_i));
                                break;
                            case PMAP_DTYPE_UDEF:
                                break;
                        }
                    }
                }
            }
            break;
            
        case PMAP_OUTOFCORE:
        {
            void *data_i;

            data_i = malloc(sizeof(map->bytes));
            /* reset position */
            fseek(map->data_fp,0,SEEK_SET);

            for (k=0; k<map->mz; k++) {
                for (j=0; j<map->my; j++) {
                    for (i=0; i<map->mx; i++) {

                        /* fread */
                        fread(data_i,map->bytes,1,map->data_fp);

                        switch (map->data_type) {
                            case PMAP_INT:
                                fprintf(fp,"%d ",*((int*)data_i));
                                break;
                            case PMAP_LONG:
                                fprintf(fp,"%ld ",*((long int*)data_i));
                                break;
                            case PMAP_FLOAT:
                                fprintf(fp,"%1.4e ",*((float*)data_i));
                                break;
                            case PMAP_DOUBLE:
                                fprintf(fp,"%1.4e ",*((double*)data_i));
                                break;
                            case PMAP_SHORT:
                                fprintf(fp,"%hd ",*((short*)data_i));
                                break;
                            case PMAP_CHAR:
                                fprintf(fp,"%hd ",(short)*((char*)data_i));
                                break;
                            case PMAP_DTYPE_UDEF:
                                break;
                        }
                    }
                }
            }
            free(data_i);
        }
            break;
            
        case PMAP_TYPE_UDEF:
            break;
    }
    fprintf(fp,"\n");
    fprintf(fp,"      </DataArray>\n");
    fprintf(fp,"      </CellData>\n");
    
    fprintf(fp,"    </Piece>\n");
    fprintf(fp,"  </ImageData>\n");
    fprintf(fp,"</VTKFile>\n");
	
    fclose(fp);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinCtxAttachPMap"
PetscErrorCode pTatinCtxAttachPMap(pTatinCtx ctx,PMap map)
{
	PetscErrorCode ierr;
    PetscFunctionBegin;
	ierr = pTatinCtxAttachModelData(ctx,"pmap",(void*)map);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinCtxGetPMap"
PetscErrorCode pTatinCtxGetPMap(pTatinCtx ctx,PMap *map)
{
	void           *mymap;
	PetscErrorCode ierr;
	
    PetscFunctionBegin;
	ierr = pTatinCtxGetModelData(ctx,"pmap",&mymap);CHKERRQ(ierr);
	*map = (PMap)mymap;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ex1"
PetscErrorCode ex1(void)
{
    PMap           cg;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    {
        double t[] = {1.0,2.0,3.3,4.4,5.5,6.6, 7.7,8.8,9.9,10.1,11.2,12.3 };
        FILE *fp;
        
        fp = fopen("mydata.bin","w");
        if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open mydata.bin");
        
        fwrite(t,sizeof(double),12,fp);
        fclose(fp);        
    }
    
    ierr = PMapCreate(&cg);CHKERRQ(ierr);
    ierr = PMapSetFilename(cg,"test2.dat");CHKERRQ(ierr);
    ierr = PMapSetUp(cg);CHKERRQ(ierr);
    
    {
        int            k;
        double         pos[] = {1.0,60.0};
        double         value;
        PMapValueFound found;
        
        value = 0.0;
        for (k=0; k<10; k++) {
            ierr = PMapGetValue(cg,pos,&value,&found);CHKERRQ(ierr);
        }
        printf("val %1.6e [%1.4e,%1.4e]\n",value,pos[0],pos[1]);
    }
    ierr = PMapViewPV(cg,"test.vti");CHKERRQ(ierr);
    ierr = PMapDestroy(&cg);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ex2"
PetscErrorCode ex2(void)
{
    PetscReal      xr[2],yr[2],zr[2];
    PMap           cgin,cg;
    int            i,j,k;
    long int       *data;
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    ierr = PMapCreate(&cgin);CHKERRQ(ierr);
    ierr = PMapSetFilename(cgin,"test3.pmap");CHKERRQ(ierr);
    ierr = PMapSetDataFilename(cgin,"test.bin");CHKERRQ(ierr);
    ierr = PMapSetType(cgin,PMAP_INMEM);CHKERRQ(ierr);
    ierr = PMapSetDim(cgin,3);CHKERRQ(ierr);
    ierr = PMapSetDataType(cgin,PMAP_LONG);CHKERRQ(ierr);

    xr[0] = 0.0;  xr[1] = 10.0;
    yr[0] = 0.0;  yr[1] = 13.0;
    zr[0] = 0.0;  zr[1] = 20.0;
    ierr = PMapSetDomain(cgin,xr,yr,zr);CHKERRQ(ierr);
    ierr = PMapSetSizes(cgin,30,40,50);CHKERRQ(ierr);

    ierr = PMapViewMetaData(cgin);CHKERRQ(ierr);

    data = malloc(sizeof(long int)*cgin->mx*cgin->my*cgin->mz);
    for (k=0; k<cgin->mz; k++) {
        for (j=0; j<cgin->my; j++) {
            for (i=0; i<cgin->mx; i++) {
                int index;
                
                index = i+j*cgin->mx+k*cgin->mx*cgin->my;

                data[index] = 0;
                if (i <cgin->mx/2) { data[index] = 1; }
                if (j <cgin->my/2) { data[index] = 2; }
                if (k <cgin->mz/2) { data[index] = 3; }
            }
        }
    }
    
    {
        FILE *fp;
        
        fp = fopen("test.bin","w");
        if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open test.bin");
        fwrite(data,sizeof(long int),cgin->mx*cgin->my*cgin->mz,fp);
        fclose(fp);
    }

    free(data);
    
    ierr = PMapDestroy(&cgin);CHKERRQ(ierr);
    
    ierr = PMapCreate(&cg);CHKERRQ(ierr);
    ierr = PMapSetFilename(cg,"test3.pmap");CHKERRQ(ierr);
    ierr = PMapSetUp(cg);CHKERRQ(ierr);
    ierr = PMapViewPV(cg,"test.vti");CHKERRQ(ierr);
    ierr = PMapDestroy(&cg);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#if 0
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
    PetscErrorCode ierr;
    
    ierr = PetscInitialize(&argc,&args,(char*)0,NULL);CHKERRQ(ierr);

    //ierr = ex1();CHKERRQ(ierr);
    ierr = ex2();CHKERRQ(ierr);
    
    
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
}
#endif

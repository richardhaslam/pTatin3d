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


const char *CartGridTypeNames[] = { "udef", "inmem", "outofcore", 0 };
const char *CartGridDataTypeNames[] = { "udef", "int", "long", "float", "double", "short", "char", 0 };


PetscErrorCode CartGridGetIndex_InMem(CartGrid map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index);
PetscErrorCode CartGridGetValue_InMem(CartGrid map,PetscReal xp[],void *value,PetscBool *found);
PetscErrorCode CartGridGetValue_OutOfCore(CartGrid map,PetscReal xp[],void *value,PetscBool *found);


#undef __FUNCT__
#define __FUNCT__ "CartGridCreate"
PetscErrorCode CartGridCreate(CartGrid *map)
{
  PetscErrorCode ierr;
  CartGrid       p;

  PetscFunctionBegin;
  *map = NULL;
  ierr = PetscNew(&p);CHKERRQ(ierr);

  p->type = CARTGRID_INMEM;
  p->getindex = CartGridGetIndex_InMem;
  p->getvalue = CartGridGetValue_InMem;
  p->destroy = NULL;

  *map = p;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridDestroy"
PetscErrorCode CartGridDestroy(CartGrid *map)
{
  CartGrid p;

  PetscFunctionBegin;
  if (!map) { PetscFunctionReturn(0); }
  p = *map;

  if (p->data) {
    PetscFree(p->data);
    p->data = NULL;
  }

  if (p->type == CARTGRID_OUTOFCORE) {
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
#define __FUNCT__ "CartGridSetType"
PetscErrorCode CartGridSetType(CartGrid map,CartGridType t)
{
    PetscFunctionBegin;
    map->type = t;
    
    switch (t) {
            
        case CARTGRID_INMEM:
            map->getindex = CartGridGetIndex_InMem;
            map->getvalue = CartGridGetValue_InMem;
            map->destroy = NULL;
            break;

        case CARTGRID_OUTOFCORE:
            map->getindex = CartGridGetIndex_InMem;
            map->getvalue = CartGridGetValue_OutOfCore;
            map->destroy = NULL;
            break;
            
        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"A valid CartGrid type {inmem,outofcore} must be specified");
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridSetDim"
PetscErrorCode CartGridSetDim(CartGrid map,PetscInt dim)
{
    PetscFunctionBegin;
    map->dim = dim;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridSetSizes"
PetscErrorCode CartGridSetSizes(CartGrid map,PetscInt m,PetscInt n,PetscInt p)
{
    PetscFunctionBegin;
    map->mx = m;
    map->my = n;
    map->mz = p;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridSetDataType"
PetscErrorCode CartGridSetDataType(CartGrid map,CartGridDataType type)
{
    PetscFunctionBegin;
    map->data_type = type;
    switch (type) {
        case CARTGRID_DTYPE_UDEF:
            map->bytes = -1;
            break;

        case CARTGRID_INT:
            map->bytes = sizeof(int);
            break;
            
        case CARTGRID_LONG:
            map->bytes = sizeof(long int);
            break;

        case CARTGRID_FLOAT:
            map->bytes = sizeof(float);
            break;

        case CARTGRID_DOUBLE:
            map->bytes = sizeof(double);
            break;

        case CARTGRID_SHORT:
            map->bytes = sizeof(short);
            break;
        
        case CARTGRID_CHAR:
            map->bytes = sizeof(char);
            break;

        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"A valid data type {char,short,int,long,float,double} must be specified");
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridSetDomain"
PetscErrorCode CartGridSetDomain(CartGrid map,PetscReal xr[],PetscReal yr[],PetscReal zr[])
{
    PetscFunctionBegin;
    if (xr) { map->range_x[0] = xr[0]; map->range_x[1] = xr[1]; }
    if (yr) { map->range_y[0] = yr[0]; map->range_y[1] = yr[1]; }
    if (zr) { map->range_z[0] = zr[0]; map->range_z[1] = zr[1]; }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridSetFilename"
PetscErrorCode CartGridSetFilename(CartGrid map,const char fname[])
{
    PetscFunctionBegin;
    PetscSNPrintf(map->metadatafile_name,PETSC_MAX_PATH_LEN-1,"%s",fname);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridSetDataFilename"
PetscErrorCode CartGridSetDataFilename(CartGrid map,const char fname[])
{
    PetscFunctionBegin;
    PetscSNPrintf(map->datafile_name,PETSC_MAX_PATH_LEN-1,"%s",fname);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridSetUp_InMem"
PetscErrorCode CartGridSetUp_InMem(CartGrid map)
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
#define __FUNCT__ "CartGridSetUp_OutOfCore"
PetscErrorCode CartGridSetUp_OutOfCore(CartGrid map)
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
#define __FUNCT__ "CartGridGetIndex_InMem2d"
PetscErrorCode CartGridGetIndex_InMem2d(CartGrid map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index)
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
#define __FUNCT__ "CartGridGetIndex_InMem3d"
PetscErrorCode CartGridGetIndex_InMem3d(CartGrid map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index)
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
#define __FUNCT__ "CartGridGetIndex_InMem"
PetscErrorCode CartGridGetIndex_InMem(CartGrid map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    switch (map->dim) {
        case 2:
            ierr = CartGridGetIndex_InMem2d(map,i,j,-1,index);CHKERRQ(ierr);
            break;
        case 3:
            ierr = CartGridGetIndex_InMem3d(map,i,j,k,index);CHKERRQ(ierr);
            break;
        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Dimension must be 2 or 3");
            break;
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridGetValue_InMem2d"
PetscErrorCode CartGridGetValue_InMem2d(CartGrid map,PetscReal xp[],void *value,PetscBool *found)
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
	
	CartGridGetIndex_InMem2d(map,i,j,-1,&index);
	if (index < 0) { PetscFunctionReturn(0); }
    
    *found = PETSC_TRUE;
    
    value_i = (void*)( (char*)map->data + index * map->bytes );
    PetscMemcpy(value,value_i,map->bytes);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridGetValue_InMem3d"
PetscErrorCode CartGridGetValue_InMem3d(CartGrid map,PetscReal xp[],void *value,PetscBool *found)
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
	
	CartGridGetIndex_InMem3d(map,i,j,k,&index);
	if (index < 0) { PetscFunctionReturn(0); }
    
    *found = PETSC_TRUE;
    
    value_i = (void*)( (char*)map->data + index * map->bytes );
    PetscMemcpy(value,value_i,map->bytes);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridGetValue_InMem"
PetscErrorCode CartGridGetValue_InMem(CartGrid map,PetscReal xp[],void *value,PetscBool *found)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    switch (map->dim) {
        case 2:
            ierr = CartGridGetValue_InMem2d(map,xp,value,found);CHKERRQ(ierr);
            break;
        case 3:
            ierr = CartGridGetValue_InMem3d(map,xp,value,found);CHKERRQ(ierr);
            break;
        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Dimension must be 2 or 3");
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridGetValue_OutOfCore"
PetscErrorCode CartGridGetValue_OutOfCore(CartGrid map,PetscReal xp[],void *value,PetscBool *found)
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
        
        CartGridGetIndex_InMem3d(map,i,j,k,&index);
    } else {
        CartGridGetIndex_InMem2d(map,i,j,-1,&index);
    }
	if (index < 0) { PetscFunctionReturn(0); }
    
    *found = PETSC_TRUE;

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
#define __FUNCT__ "CartGridGetValue"
PetscErrorCode CartGridGetValue(CartGrid map,PetscReal xp[],void *value,PetscBool *found)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  (*found) = PETSC_FALSE;
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
#define __FUNCT__ "CartGridSetUp"
PetscErrorCode CartGridSetUp(CartGrid map)
{
    FILE           *fp;
    int            t,dim,mx,my,mz;
    double         xr[2],yr[2],zr[2],dx,dy,dz;
    char           str[256];
    int            matched;
    const int      STRMAX = 256;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
	/* open file to parse */
    fp = fopen(map->metadatafile_name,"r");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open %s",map->metadatafile_name);

    /* header */
ParseHeader:
    if (!fgets(str,STRMAX,fp)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_READ,"fgets() failed");
    if (str[0] == '#') { goto ParseHeader; }
    if (str[0] == '!') { goto ParseHeader; }
    if (str[0] == '%') { goto ParseHeader; }
    if (str[0] == '#') { goto ParseHeader; }
    
    
    /* CartGrid type */
    {
        CartGridType type;
        
        trimnewline(str);
        trimright(str);
        type = CARTGRID_TYPE_UDEF;
        for (t=1; t<3; t++) {
            matched = strcmp(str,CartGridTypeNames[t]);
            if (matched == 0) {
                type = (CartGridType)t;
                break;
            }
        }
        ierr = CartGridSetType(map,type);CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD,"CartGrid:type: %s\n",CartGridTypeNames[(int)map->type]);
    }
    
    /* data_type */
    if( fgets(str,STRMAX,fp) != NULL) {
        CartGridDataType data_type;
        
        trimnewline(str);
        trimright(str);
        data_type = CARTGRID_DTYPE_UDEF;
        for (t=1; t<7; t++) {
            matched = strcmp(str,CartGridDataTypeNames[t]);
            if (matched == 0) {
                data_type = (CartGridDataType)t;
                break;
            }
        }
        ierr = CartGridSetDataType(map,data_type);CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD,"CartGrid:datatype: %s\n",CartGridDataTypeNames[(int)map->data_type]);
    }
    
    /* datafile_name */
    if( fgets(str,STRMAX,fp) != NULL) {
        trimnewline(str);
        trimright(str);
        PetscPrintf(PETSC_COMM_WORLD,"CartGrid:datafile: %s\n",str);
        ierr = CartGridSetDataFilename(map,str);CHKERRQ(ierr);
    }
    
    /* read header: dim */
    if (fscanf(fp,"%d\n",&dim) < 1) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
    map->dim = (PetscInt)dim;
    
    if (dim == 2) {
        /* read header: mx,my */
        if (fscanf(fp,"%d %d",&mx,&my) < 2) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
        mz = 1;
        /* read header: x0,x1,y0,y1 */
        if (fscanf(fp,"%lf %lf %lf %lf",&xr[0],&xr[1],&yr[0],&yr[1]) < 4) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
        zr[0] = -1.0e32;
        zr[1] =  1.032;
        printf("[%d x %d] --> [%1.4e,%1.4e] x [%1.4e,%1.4e]\n",mx,my,xr[0],xr[1],yr[0],yr[1]);
        
        dx = (xr[1] - xr[0])/(double)(mx);
        dy = (yr[1] - yr[0])/(double)(my);
        dz = 0.0;
    } else if (dim == 3) {
        /* read header: mx,my,mz */
        if (fscanf(fp,"%d %d %d",&mx,&my,&mz) < 3) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
        /* read header: x0,x1,y0,y1,z0,z1 */
        if (fscanf(fp,"%lf %lf %lf %lf %lf %lf",&xr[0],&xr[1],&yr[0],&yr[1],&zr[0],&zr[1]) < 6) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
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
        case CARTGRID_INMEM:
            ierr = CartGridSetUp_InMem(map);CHKERRQ(ierr);
            break;
        case CARTGRID_OUTOFCORE:
            ierr = CartGridSetUp_OutOfCore(map);CHKERRQ(ierr);
            break;
        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"A valid type {inmem,outofcore} must be specified");
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridViewMetaData"
PetscErrorCode CartGridViewMetaData(CartGrid map)
{
    FILE *fp;
    
    PetscFunctionBegin;
    fp = fopen(map->metadatafile_name,"w");
    if (!fp)  SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open %s",map->metadatafile_name);
    
    fprintf(fp,"%s\n",CartGridTypeNames[(int)map->type]);
    fprintf(fp,"%s\n",CartGridDataTypeNames[(int)map->data_type]);
    fprintf(fp,"%s\n",map->datafile_name);

    fprintf(fp,"%d\n",map->dim);
    if (map->dim == 2) {
        fprintf(fp,"%d %d\n",(int)map->mx,(int)map->my);
        fprintf(fp,"%1.4e %1.4e %1.4e %1.4e\n",
                map->range_x[0],map->range_x[1],
                map->range_y[0],map->range_y[1]);
    } else {
        fprintf(fp,"%d %d %d\n",(int)map->mx,(int)map->my,(int)map->mz);
        fprintf(fp,"%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n",
                map->range_x[0],map->range_x[1],
                map->range_y[0],map->range_y[1],
                map->range_z[0],map->range_z[1]);
    }
    fclose(fp);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CartGridViewPV"
PetscErrorCode CartGridViewPV(CartGrid map,const char filename[])
{
	FILE     *fp = NULL;
	PetscInt i,j,k;
	
    PetscFunctionBegin;
	/* open file to parse */
    fp = fopen(filename,"w");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open %s",filename);

    fprintf(fp,"<?xml version=\"1.0\"?>\n");
#ifdef WORDSIZE_BIGENDIAN
    fprintf(fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
    fprintf(fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
    if (map->dim == 2) {
        fprintf(fp,"  <ImageData WholeExtent=\"%d %d %d %d %d %d\"\n",0,(int)map->mx,0,(int)map->my,0,1);
        fprintf(fp,"   Origin=\"%1.4e %1.4e 0.0\"\n",map->range_x[0],map->range_y[0]);
        fprintf(fp,"   Spacing=\"%1.2e %1.2e %1.2e\">\n",map->dx,map->dy,1.0e-3);
        fprintf(fp,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",0,(int)map->mx,0,(int)map->my,0,1);
    } else {
        fprintf(fp,"  <ImageData WholeExtent=\"%d %d %d %d %d %d\"\n",0,(int)map->mx,0,(int)map->my,0,(int)map->mz);
        fprintf(fp,"   Origin=\"%1.4e %1.4e %1.4e\"\n",map->range_x[0],map->range_y[0],map->range_z[0]);
        fprintf(fp,"   Spacing=\"%1.2e %1.2e %1.2e\">\n",map->dx,map->dy,map->dz);
        fprintf(fp,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",0,(int)map->mx,0,(int)map->my,0,(int)map->mz);
    }

    fprintf(fp,"      <CellData>\n");
    switch (map->data_type) {
        case CARTGRID_INT:
            fprintf(fp,"      <DataArray type=\"Int32\" Name=\"values\" format=\"ascii\">\n");
            break;
        case CARTGRID_LONG:
            fprintf(fp,"      <DataArray type=\"Int64\" Name=\"values\" format=\"ascii\">\n");
            break;
        case CARTGRID_FLOAT:
            fprintf(fp,"      <DataArray type=\"Float32\" Name=\"values\" format=\"ascii\">\n");
            break;
        case CARTGRID_DOUBLE:
            fprintf(fp,"      <DataArray type=\"Float64\" Name=\"values\" format=\"ascii\">\n");
            break;
        case CARTGRID_SHORT:
            fprintf(fp,"      <DataArray type=\"Int16\" Name=\"values\" format=\"ascii\">\n");
            break;
        case CARTGRID_CHAR:
            fprintf(fp,"      <DataArray type=\"Int8\" Name=\"values\" format=\"ascii\">\n");
            break;
        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"A valid VTK data type must be specified");
            break;
    }

    switch (map->type) {
        case CARTGRID_INMEM:
            
            for (k=0; k<map->mz; k++) {
                for (j=0; j<map->my; j++) {
                    for (i=0; i<map->mx; i++) {
                        PetscInt index;
                        void     *data_i;
                        
                        index = i + j * map->mx + k * map->mx * map->my;
                        
                        data_i = (void*)( (char*)map->data + index * map->bytes );
                        
                        switch (map->data_type) {
                            case CARTGRID_INT:
                                fprintf(fp,"%d ",*((int*)data_i));
                                break;
                            case CARTGRID_LONG:
                                fprintf(fp,"%ld ",*((long int*)data_i));
                                break;
                            case CARTGRID_FLOAT:
                                fprintf(fp,"%1.4e ",*((float*)data_i));
                                break;
                            case CARTGRID_DOUBLE:
                                fprintf(fp,"%1.4e ",*((double*)data_i));
                                break;
                            case CARTGRID_SHORT:
                                fprintf(fp,"%hd ",*((short*)data_i));
                                break;
                            case CARTGRID_CHAR:
                                fprintf(fp,"%hd ",(short)*((char*)data_i));
                                break;
                            case CARTGRID_DTYPE_UDEF:
                                break;
                            default:
                                SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"A valid type {char,short,int,long,float,double} must be specified");
                                break;
                        }
                    }
                }
            }
            break;
            
        case CARTGRID_OUTOFCORE:
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
                            case CARTGRID_INT:
                                fprintf(fp,"%d ",*((int*)data_i));
                                break;
                            case CARTGRID_LONG:
                                fprintf(fp,"%ld ",*((long int*)data_i));
                                break;
                            case CARTGRID_FLOAT:
                                fprintf(fp,"%1.4e ",*((float*)data_i));
                                break;
                            case CARTGRID_DOUBLE:
                                fprintf(fp,"%1.4e ",*((double*)data_i));
                                break;
                            case CARTGRID_SHORT:
                                fprintf(fp,"%hd ",*((short*)data_i));
                                break;
                            case CARTGRID_CHAR:
                                fprintf(fp,"%hd ",(short)*((char*)data_i));
                                break;
                            case CARTGRID_DTYPE_UDEF:
                                break;
                            default:
                                SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"A valid type {char,short,int,long,float,double} must be specified");
                                break;
                        }
                    }
                }
            }
            free(data_i);
        }
            break;
            
        case CARTGRID_TYPE_UDEF:
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
#define __FUNCT__ "pTatinCtxAttachCartGrid"
PetscErrorCode pTatinCtxAttachCartGrid(pTatinCtx ctx,CartGrid map)
{
	PetscErrorCode ierr;
    PetscFunctionBegin;
	ierr = pTatinCtxAttachModelData(ctx,"CartGrid",(void*)map);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinCtxGetCartGrid"
PetscErrorCode pTatinCtxGetCartGrid(pTatinCtx ctx,CartGrid *map)
{
	void           *mymap;
	PetscErrorCode ierr;
	
    PetscFunctionBegin;
	ierr = pTatinCtxGetModelData(ctx,"CartGrid",&mymap);CHKERRQ(ierr);
	*map = (CartGrid)mymap;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ex1"
PetscErrorCode ex1(void)
{
    CartGrid       cg;
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
    
    ierr = CartGridCreate(&cg);CHKERRQ(ierr);
    ierr = CartGridSetFilename(cg,"test2.dat");CHKERRQ(ierr);
    ierr = CartGridSetUp(cg);CHKERRQ(ierr);
    
    {
        int       k;
        double    pos[] = {1.0,60.0};
        double    value;
        PetscBool found;
        
        value = 0.0;
        for (k=0; k<10; k++) {
            ierr = CartGridGetValue(cg,pos,&value,&found);CHKERRQ(ierr);
        }
        printf("val %1.6e [%1.4e,%1.4e]\n",value,pos[0],pos[1]);
    }
    ierr = CartGridViewPV(cg,"test.vti");CHKERRQ(ierr);
    ierr = CartGridDestroy(&cg);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ex2"
PetscErrorCode ex2(void)
{
    PetscReal      xr[2],yr[2],zr[2];
    CartGrid       cgin,cg;
    int            i,j,k;
    long int       *data;
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    ierr = CartGridCreate(&cgin);CHKERRQ(ierr);
    ierr = CartGridSetFilename(cgin,"test3.CartGrid");CHKERRQ(ierr);
    ierr = CartGridSetDataFilename(cgin,"test.bin");CHKERRQ(ierr);
    ierr = CartGridSetType(cgin,CARTGRID_INMEM);CHKERRQ(ierr);
    ierr = CartGridSetDim(cgin,3);CHKERRQ(ierr);
    ierr = CartGridSetDataType(cgin,CARTGRID_LONG);CHKERRQ(ierr);

    xr[0] = 0.0;  xr[1] = 10.0;
    yr[0] = 0.0;  yr[1] = 13.0;
    zr[0] = 0.0;  zr[1] = 20.0;
    ierr = CartGridSetDomain(cgin,xr,yr,zr);CHKERRQ(ierr);
    ierr = CartGridSetSizes(cgin,30,40,50);CHKERRQ(ierr);

    ierr = CartGridViewMetaData(cgin);CHKERRQ(ierr);

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
    
    ierr = CartGridDestroy(&cgin);CHKERRQ(ierr);
    
    ierr = CartGridCreate(&cg);CHKERRQ(ierr);
    ierr = CartGridSetFilename(cg,"test3.CartGrid");CHKERRQ(ierr);
    ierr = CartGridSetUp(cg);CHKERRQ(ierr);
    ierr = CartGridViewPV(cg,"test.vti");CHKERRQ(ierr);
    ierr = CartGridDestroy(&cg);CHKERRQ(ierr);
    
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

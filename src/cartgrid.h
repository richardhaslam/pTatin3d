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
 **    filename:   cartgrid.h
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

#ifndef __PMAP_H__
#define __PMAP_H__

typedef enum { PMAP_OUTSIDE=0, PMAP_INSIDE=1 } PMapValueFound;

typedef enum { PMAP_ENDIAN_LITTLE , PMPA_ENDIAN_BIG } PMapEndian;

typedef enum { PMAP_TYPE_UDEF=0, PMAP_INMEM, PMAP_OUTOFCORE } PMapType; /* def, blk, ooc */
extern const char *PMapTypeNames[];

typedef enum { PMAP_DTYPE_UDEF=0, PMAP_INT, PMAP_LONG, PMAP_FLOAT, PMAP_DOUBLE, PMAP_SHORT, PMAP_CHAR } PMapDataType; /* int, lng, flt, dbl */
extern const char *PMapDataTypeNames[];

typedef struct _p_PMap *PMap;

struct _p_PMap {
    PetscInt     dim;
    PetscReal    range_x[2],range_y[2],range_z[2];
    PMapDataType data_type;
    size_t       bytes;
    PetscInt     start[3],end[3];
    PetscInt     mx,my,mz;
    PMapType     type;
    PMapEndian   endian;
	PetscReal    dx,dy,dz;
    char         metadatafile_name[PETSC_MAX_PATH_LEN];
    char         datafile_name[PETSC_MAX_PATH_LEN];
    FILE         *data_fp;
	void         *data;
    PetscErrorCode (*getindex)(PMap,PetscInt,PetscInt,PetscInt,PetscInt*);
    PetscErrorCode (*getvalue)(PMap,PetscReal*,void*,PMapValueFound*);
    PetscErrorCode (*destroy)(PMap);
};

PetscErrorCode PMapCreate(PMap *map);

PetscErrorCode PMapSetType(PMap map,PMapType t);
PetscErrorCode PMapSetDim(PMap map,PetscInt dim);
PetscErrorCode PMapSetSizes(PMap map,PetscInt m,PetscInt n,PetscInt p);
PetscErrorCode PMapSetDataType(PMap map,PMapDataType type);
PetscErrorCode PMapSetDomain(PMap map,PetscReal xr[],PetscReal yr[],PetscReal zr[]);
PetscErrorCode PMapSetFilename(PMap map,const char fname[]);
PetscErrorCode PMapSetDataFilename(PMap map,const char fname[]);
PetscErrorCode PMapViewMetaData(PMap map);

PetscErrorCode PMapSetUp(PMap map);

PetscErrorCode PMapDestroy(PMap *map);


PetscErrorCode PMapGetIndex(PMap map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index);

PetscErrorCode PMapGetValue(PMap map,PetscReal xp[],void *value,PMapValueFound*);

PetscErrorCode PMapGetDataRange(PMap map,void *min,void *max);

PetscErrorCode PMapView(PMap map);
PetscErrorCode PMapViewPV(PMap map,const char filename[]);


#endif


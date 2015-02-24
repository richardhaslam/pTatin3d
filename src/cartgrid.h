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

#ifndef __CARTGRID_H__
#define __CARTGRID_H__

typedef enum { CARTGRID_ENDIAN_LITTLE , PMPA_ENDIAN_BIG } CartGridEndian;

typedef enum { CARTGRID_TYPE_UDEF=0, CARTGRID_INMEM, CARTGRID_OUTOFCORE } CartGridType; /* def, blk, ooc */
extern const char *CartGridTypeNames[];

typedef enum { CARTGRID_DTYPE_UDEF=0, CARTGRID_INT, CARTGRID_LONG, CARTGRID_FLOAT, CARTGRID_DOUBLE, CARTGRID_SHORT, CARTGRID_CHAR } CartGridDataType; /* int, lng, flt, dbl */
extern const char *CartGridDataTypeNames[];

typedef struct _p_CartGrid *CartGrid;

struct _p_CartGrid {
    PetscInt     dim;
    PetscReal    range_x[2],range_y[2],range_z[2];
    CartGridDataType data_type;
    size_t       bytes;
    PetscInt     start[3],end[3];
    PetscInt     mx,my,mz;
    CartGridType     type;
    CartGridEndian   endian;
	PetscReal    dx,dy,dz;
    char         metadatafile_name[PETSC_MAX_PATH_LEN];
    char         datafile_name[PETSC_MAX_PATH_LEN];
    FILE         *data_fp;
	void         *data;
    PetscErrorCode (*getindex)(CartGrid,PetscInt,PetscInt,PetscInt,PetscInt*);
    PetscErrorCode (*getvalue)(CartGrid,PetscReal*,void*,PetscBool*);
    PetscErrorCode (*destroy)(CartGrid);
};

PetscErrorCode CartGridCreate(CartGrid *map);

PetscErrorCode CartGridSetType(CartGrid map,CartGridType t);
PetscErrorCode CartGridSetDim(CartGrid map,PetscInt dim);
PetscErrorCode CartGridSetSizes(CartGrid map,PetscInt m,PetscInt n,PetscInt p);
PetscErrorCode CartGridSetDataType(CartGrid map,CartGridDataType type);
PetscErrorCode CartGridSetDomain(CartGrid map,PetscReal xr[],PetscReal yr[],PetscReal zr[]);
PetscErrorCode CartGridSetFilename(CartGrid map,const char fname[]);
PetscErrorCode CartGridSetDataFilename(CartGrid map,const char fname[]);
PetscErrorCode CartGridViewMetaData(CartGrid map);

PetscErrorCode CartGridSetUp(CartGrid map);

PetscErrorCode CartGridDestroy(CartGrid *map);

PetscErrorCode CartGridGetIndex(CartGrid map,PetscInt i,PetscInt j,PetscInt k,PetscInt *index);

PetscErrorCode CartGridGetValue(CartGrid map,PetscReal xp[],void *value,PetscBool*);

PetscErrorCode CartGridView(CartGrid map);
PetscErrorCode CartGridViewPV(CartGrid map,const char filename[]);


#endif


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
 **    filename:   quantity.h
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

#ifndef _ptatin_quantity_h__
#define _ptatin_quantity_h__

#include <petsc.h>

typedef struct _p_Quantity *Quantity;

struct _p_Quantity {
  char      name[PETSC_MAX_PATH_LEN],unit[PETSC_MAX_PATH_LEN];
  PetscReal magnitude,recip_magnitude;
};

extern Quantity *PTATIN_SI_QLIST;
extern Quantity *PTATIN_GEO_QLIST;

#define QuantityGetRecipMagnitude(Q) (Q)->magnitude
#define QuantityGetMagnitude(Q) (Q)->magnitude
#define QuantityGetUnit(Q) (Q)->unit

#define QuantityApplyND(Q,X) (X) * (Q)->recip_magnitude
#define QuantityApplyD(Q,X) (X) * (Q)->magnitude

#define pQuantitySIGetRecipMagnitude(type) PTATIN_SI_QLIST[(type)]->magnitude
#define pQuantitySIGetMagnitude(type)      PTATIN_SI_QLIST[(type)]->magnitude
#define pQuantitySIGetUnit(type)           PTATIN_SI_QLIST[(type)]->unit

#define pQuantityGeoGetRecipMagnitude(type) PTATIN_GEO_QLIST[(type)]->magnitude
#define pQuantityGeoGetMagnitude(type)      PTATIN_GEO_QLIST[(type)]->magnitude
#define pQuantityGeoGetUnit(type)           PTATIN_GEO_QLIST[(type)]->unit

#define pQConvert_GeoUnits2ModelUnits(type,val) QuantityApplyND(PTATIN_GEO_QLIST[(type)],(val))
#define pQConvert_ModelUnits2GeoUnits(type,val) QuantityApplyD(PTATIN_GEO_QLIST[(type)],(val))
#define pQConvert_SIUnits2ModelUnits(type,val)  QuantityApplyND(PTATIN_SI_QLIST[(type)],(val))
#define pQConvert_ModelUnits2SIUnits(type,val)  QuantityApplyD(PTATIN_SI_QLIST[(type)],(val))

typedef enum { UNITS_NATIVE = 0, UNITS_SI, UNITS_GEO } ptatinInputUnits;

typedef enum {
  QLength = 0,
  QTime,
  QVelocity,
  QStress,
  QViscosity,
  QStrainRate,
  QDensity,
  QAcceleration,
  QNull
} QuantityType;
extern const char *QuantityTypeNames[];
extern const char *QuantityTypeUnits[];
extern const char *ptatinInputUnitsNames[];


/*
void UnitsApplyScaling(Units u,double y,double *ystar);
void UnitsApplyScalingToArray(Units u,int n,double y[],double ystar[]);
void UnitsApplyInverseScaling(Units u,double y,double *ystar);
void UnitsApplyInverseScalingToArray(Units u,int n,double y[],double ystar[]);

void UnitsConvertSI2Kilo(double si,double *ksi);
void UnitsConvertSI2Mega(double si,double *ksi);
void UnitsConvertSI2Giga(double si,double *ksi);
void UnitsConvertSI2Milli(double si,double *ksi);
void UnitsConvertSI2Micro(double si,double *ksi);
void UnitsConvertSI2Pico(double si,double *ksi);

void UnitsCreate(Units *_u,const char quantity_name[],const char unit_name[],double scale);
void UnitsSetScale(Units u,double scale);
void UnitsView(Units u);

void UnitsCreateSIVelocity(Units *u);
void UnitsCreateSILength(Units *u);
void UnitsCreateSITime(Units *u);
void UnitsCreateSIStress(Units *u);
void UnitsCreateSIViscosity(Units *u);
void UnitsCreateSITemperature(Units *u);
*/

PetscErrorCode QuantityCreate(Quantity *q);
PetscErrorCode QuantitySetValues(Quantity q,const char name[],const char unit[],PetscReal mag);
PetscErrorCode QuantitySetMagnitude(Quantity q,PetscReal mag);
PetscErrorCode QuantitySetUnit(Quantity q,const char unit[]);
PetscErrorCode QuantityView(Quantity q);

extern PetscErrorCode pQConvert_GeoUnits2ModelUnits_Vec(QuantityType type,Vec x);
extern PetscErrorCode pQConvert_GeoUnits2ModelUnits_ArrayInPlace(QuantityType type,PetscInt len,PetscReal vals[]);
extern PetscErrorCode pQConvert_GeoUnits2ModelUnits_Array(QuantityType type,PetscInt len,PetscReal ain[],PetscReal aout[]);

extern PetscErrorCode pQConvert_ModelUnits2GeoUnits_Vec(QuantityType type,Vec x);
extern PetscErrorCode pQConvert_ModelUnits2GeoUnits_ArrayInPlace(QuantityType type,PetscInt len,PetscReal vals[]);
extern PetscErrorCode pQConvert_ModelUnits2GeoUnits_Array(QuantityType type,PetscInt len,PetscReal ain[],PetscReal aout[]);

extern PetscErrorCode pQConvert_SIUnits2ModelUnits_Vec(QuantityType type,Vec x);
extern PetscErrorCode pQConvert_SIUnits2ModelUnits_ArrayInPlace(QuantityType type,PetscInt len,PetscReal vals[]);
extern PetscErrorCode pQConvert_SIUnits2ModelUnits_Array(QuantityType type,PetscInt len,PetscReal ain[],PetscReal aout[]);

extern PetscErrorCode pQConvert_ModelUnits2SIUnits_Vec(QuantityType type,Vec x);
extern PetscErrorCode pQConvert_ModelUnits2SIUnits_ArrayInPlace(QuantityType type,PetscInt len,PetscReal vals[]);
extern PetscErrorCode pQConvert_ModelUnits2SIUnits_Array(QuantityType type,PetscInt len,PetscReal ain[],PetscReal aout[]);

PetscErrorCode ptatinQuantityCreate(void);
PetscErrorCode ptatinQuantityDestroy(void);

#endif


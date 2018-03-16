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
 **    filename:   quantity.c
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

#include <petsc.h>
#include <quantity.h>

const char *QuantityTypeNames[] = { "length",  "time", "velocity", "stress", "viscosity", "strain-rate", "density",  "acceleration", 0 };
const char *QuantityTypeUnits[] = { "[m]",     "[s]",  "[m/s]",    "[Pa]",   "[Pa s]",    "[1/s]",       "[kg/m^3]", "[m/s^2]",      0 };
const char *ptatinInputUnitsNames[] = { "native", "si", "geo", "ptatinInputUnits", "PTATIN_INPUT_UNITS_", 0 };

/* global variables */
Quantity *PTATIN_SI_QLIST = NULL;
Quantity *PTATIN_GEO_QLIST = NULL;


#undef __FUNCT__
#define __FUNCT__ "QuantityView"
PetscErrorCode QuantityView(Quantity q)
{
  PetscPrintf(PETSC_COMM_WORLD,"  \"%20s\" %1.4e %s \n",q->name,q->magnitude,q->unit);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QuantityCreate"
PetscErrorCode QuantityCreate(Quantity *q)
{
  Quantity       qq;
  PetscErrorCode ierr;
  
  ierr = PetscMalloc1(1,&qq);CHKERRQ(ierr);
  ierr = PetscMemzero(qq,sizeof(struct _p_Quantity));
  qq->magnitude = 1.0;
  qq->recip_magnitude = 1.0;
  *q = qq;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QuantitySetValues"
PetscErrorCode QuantitySetValues(Quantity q,const char name[],const char unit[],PetscReal mag)
{
  PetscSNPrintf(q->name,PETSC_MAX_PATH_LEN,"%s",name);
  PetscSNPrintf(q->unit,PETSC_MAX_PATH_LEN,"%s",unit);
  q->magnitude = mag;
  q->recip_magnitude = 1.0/mag;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QuantitySetMagnitude"
PetscErrorCode QuantitySetMagnitude(Quantity q,PetscReal mag)
{
  q->magnitude = mag;
  q->recip_magnitude = 1.0/mag;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QuantitySetUnit"
PetscErrorCode QuantitySetUnit(Quantity q,const char unit[])
{
  PetscSNPrintf(q->unit,PETSC_MAX_PATH_LEN,"%s",unit);
  PetscFunctionReturn(0);
}

/* Conversion routines */
// [1] geo -> model
// [2] model -> geo
// [3] si -> model
// [4] model -> si

// [1] geo -> model
/*  #define QuantityXYZ_ConvertFromGeoUnit2ModelUnit(g) QuantityApplyND(geo_list[QXYZ],(g)) */
#undef __FUNCT__
#define __FUNCT__ "pQConvert_GeoUnits2ModelUnits_Vec"
inline PetscErrorCode pQConvert_GeoUnits2ModelUnits_Vec(QuantityType type,Vec x)
{
  PetscErrorCode ierr;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  ierr = VecScale(x,QuantityGetRecipMagnitude(PTATIN_GEO_QLIST[type]));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pQConvert_GeoUnits2ModelUnits_ArrayInPlace"
inline PetscErrorCode pQConvert_GeoUnits2ModelUnits_ArrayInPlace(QuantityType type,PetscInt len,PetscReal vals[])
{
  PetscInt i;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  for (i=0; i<len; i++) {
    vals[i] = QuantityApplyND(PTATIN_GEO_QLIST[type],vals[i]);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pQConvert_GeoUnits2ModelUnits_Array"
inline PetscErrorCode pQConvert_GeoUnits2ModelUnits_Array(QuantityType type,PetscInt len,PetscReal ain[],PetscReal aout[])
{
  PetscInt i;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  for (i=0; i<len; i++) {
    aout[i] = QuantityApplyND(PTATIN_GEO_QLIST[type],ain[i]);
  }
  PetscFunctionReturn(0);
}

// [2] model -> geo
/*  #define QuantityXYZ_ConvertFromModelUnit2GeoUnit(m) QuantityApplyD(geo_list[QXYZ],(m)) */
#undef __FUNCT__
#define __FUNCT__ "pQConvert_ModelUnits2GeoUnits_Vec"
inline PetscErrorCode pQConvert_ModelUnits2GeoUnits_Vec(QuantityType type,Vec x)
{
  PetscErrorCode ierr;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  ierr = VecScale(x,QuantityGetMagnitude(PTATIN_GEO_QLIST[type]));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pQConvert_ModelUnits2GeoUnits_ArrayInPlace"
inline PetscErrorCode pQConvert_ModelUnits2GeoUnits_ArrayInPlace(QuantityType type,PetscInt len,PetscReal vals[])
{
  PetscInt i;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  for (i=0; i<len; i++) {
    vals[i] = QuantityApplyD(PTATIN_GEO_QLIST[type],vals[i]);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pQConvert_ModelUnits2GeoUnits_Array"
inline PetscErrorCode pQConvert_ModelUnits2GeoUnits_Array(QuantityType type,PetscInt len,PetscReal ain[],PetscReal aout[])
{
  PetscInt i;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  for (i=0; i<len; i++) {
    aout[i] = QuantityApplyD(PTATIN_GEO_QLIST[type],ain[i]);
  }
  PetscFunctionReturn(0);
}

// [3] si -> model
/* #define QuantityXYZ_ConvertFromSI2ModelUnit(g)      QuantityApplyND(list[QXYZ],(g)) */
/*
#undef __FUNCT__
#define __FUNCT__ "pQConvert_SIUnits2ModelUnits"
inline PetscErrorCode pQConvert_SIUnits2ModelUnits(QuantityType type,PetscReal *a)
{
  if (type > QNull) printf("unknown type detected\n");
  if (type == QNull) printf("NULL type detected\n");
  *a = QuantityApplyND(PTATIN_SI_QLIST[type],*a);
  PetscFunctionReturn(0);
}
*/

#undef __FUNCT__
#define __FUNCT__ "pQConvert_SIUnits2ModelUnits_Vec"
inline PetscErrorCode pQConvert_SIUnits2ModelUnits_Vec(QuantityType type,Vec x)
{
  PetscErrorCode ierr;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  ierr = VecScale(x,QuantityGetRecipMagnitude(PTATIN_SI_QLIST[type]));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pQConvert_SIUnits2ModelUnits_ArrayInPlace"
inline PetscErrorCode pQConvert_SIUnits2ModelUnits_ArrayInPlace(QuantityType type,PetscInt len,PetscReal vals[])
{
  PetscInt i;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  for (i=0; i<len; i++) {
    vals[i] = QuantityApplyND(PTATIN_SI_QLIST[type],vals[i]);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pQConvert_SIUnits2ModelUnits_Array"
inline PetscErrorCode pQConvert_SIUnits2ModelUnits_Array(QuantityType type,PetscInt len,PetscReal ain[],PetscReal aout[])
{
  PetscInt i;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  for (i=0; i<len; i++) {
    aout[i] = QuantityApplyND(PTATIN_SI_QLIST[type],ain[i]);
  }
  PetscFunctionReturn(0);
}

// [4] model -> si
/* #define QuantityXYZ_ConvertFromModelUnit2SI(m)      QuantityApplyD(list[QXYZ],(m)) */
#undef __FUNCT__
#define __FUNCT__ "pQConvert_ModelUnits2SIUnits_Vec"
inline PetscErrorCode pQConvert_ModelUnits2SIUnits_Vec(QuantityType type,Vec x)
{
  PetscErrorCode ierr;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  ierr = VecScale(x,QuantityGetMagnitude(PTATIN_SI_QLIST[type]));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pQConvert_ModelUnits2SIUnits_ArrayInPlace"
inline PetscErrorCode pQConvert_ModelUnits2SIUnits_ArrayInPlace(QuantityType type,PetscInt len,PetscReal vals[])
{
  PetscInt i;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  for (i=0; i<len; i++) {
    vals[i] = QuantityApplyD(PTATIN_SI_QLIST[type],vals[i]);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pQConvert_ModelUnits2SIUnits_Array"
inline PetscErrorCode pQConvert_ModelUnits2SIUnits_Array(QuantityType type,PetscInt len,PetscReal ain[],PetscReal aout[])
{
  PetscInt i;
  if (type > QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"unknown quantity detected");
  if (type == QNull) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"NULL type detected");
  for (i=0; i<len; i++) {
    aout[i] = QuantityApplyD(PTATIN_SI_QLIST[type],ain[i]);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ptatinQuantitySetupSIUnitList_UsingLVEScales"
PetscErrorCode ptatinQuantitySetupSIUnitList_UsingLVEScales(PetscReal LengScale,PetscReal VeloScale,PetscReal ViscScale)
{
  PetscReal Tstar;
  
  // geo scales //
  QuantitySetMagnitude(PTATIN_SI_QLIST[QLength],LengScale);
  QuantitySetMagnitude(PTATIN_SI_QLIST[QVelocity],VeloScale);
  QuantitySetMagnitude(PTATIN_SI_QLIST[QViscosity],ViscScale);
  
  QuantitySetMagnitude(PTATIN_SI_QLIST[QTime],      LengScale/VeloScale);
  QuantitySetMagnitude(PTATIN_SI_QLIST[QStrainRate],VeloScale/LengScale);
  QuantitySetMagnitude(PTATIN_SI_QLIST[QStress],    ViscScale/QuantityGetMagnitude(PTATIN_SI_QLIST[QTime]));
  
  // stress = F/m^2 = (kg.m/s^2) / m^2 = kg.m^3/s^2 //
  Tstar = QuantityGetMagnitude(PTATIN_SI_QLIST[QTime]);
  //QuantitySetMagnitude(list[QDensity],     QuantityGetMagnitude(list[QStress])*Tstar*Tstar/(LengScale*LengScale*LengScale * LengScale*LengScale*LengScale) );
  //QuantitySetMagnitude(list[QAcceleration],VeloScale/QuantityGetMagnitude(list[QTime]));
  
  QuantitySetMagnitude(PTATIN_SI_QLIST[QAcceleration],VeloScale / QuantityGetMagnitude(PTATIN_SI_QLIST[QTime]));
  QuantitySetMagnitude(PTATIN_SI_QLIST[QDensity],     ViscScale * Tstar / (LengScale * LengScale) );
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ptatinQuantitySetupGeoUnitList"
PetscErrorCode ptatinQuantitySetupGeoUnitList(void)
{
  PetscReal secondsPerYear,secondsPerkyr;
  
  /* geo units */
  /* reset scale factors */
  QuantitySetMagnitude(PTATIN_GEO_QLIST[QLength],QuantityGetMagnitude(PTATIN_SI_QLIST[QLength])/1.0e3);
  QuantitySetUnit(PTATIN_GEO_QLIST[QLength],"[km]");
  
  secondsPerYear = 60.0 * 60.0 * 24.0 * 365.0;
  secondsPerkyr = secondsPerYear * 1.0e3;
  QuantitySetMagnitude(PTATIN_GEO_QLIST[QTime],QuantityGetMagnitude(PTATIN_SI_QLIST[QTime])/secondsPerkyr);
  QuantitySetUnit(PTATIN_GEO_QLIST[QTime],"[kyr]");
  
  QuantitySetMagnitude(PTATIN_GEO_QLIST[QVelocity],QuantityGetMagnitude(PTATIN_SI_QLIST[QVelocity])/(1.0e-2/secondsPerYear));
  QuantitySetUnit(PTATIN_GEO_QLIST[QVelocity],"[cm/yr]");
  
  QuantitySetMagnitude(PTATIN_GEO_QLIST[QStress],QuantityGetMagnitude(PTATIN_SI_QLIST[QStress])/1.0e6);
  QuantitySetUnit(PTATIN_GEO_QLIST[QStress],"[MPa]");
  
  QuantitySetMagnitude(PTATIN_GEO_QLIST[QStrainRate],QuantityGetMagnitude(PTATIN_SI_QLIST[QStrainRate]) * secondsPerkyr);
  QuantitySetUnit(PTATIN_GEO_QLIST[QStrainRate],"[1/kyr]");
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ptatinQuantitySetup_UsingLVEScales"
PetscErrorCode ptatinQuantitySetup_UsingLVEScales(ptatinInputUnits units)
{
  PetscInt       i,n;
  PetscReal      Ls,Vs,Es,geo_time,geo_length;
  PetscReal      secondsPerYear,secondsPerkyr;
  PetscErrorCode ierr;
  
  if (units == UNITS_GEO) PetscPrintf(PETSC_COMM_WORLD,"pTatin: \"geo\" units are <length, velocity, time, stress, strain-rate> are <km, cm/yr, kyr, MPa, 1/kyr>\n");
  Ls = 1.0;
  Vs = 1.0;
  Es = 1.0;
  PetscOptionsGetReal(NULL,NULL,"-ptatin_q_length_mag",&Ls,NULL);
  PetscOptionsGetReal(NULL,NULL,"-ptatin_q_velocity_mag",&Vs,NULL);
  PetscOptionsGetReal(NULL,NULL,"-ptatin_q_viscosity_mag",&Es,NULL);

  secondsPerYear = 60.0 * 60.0 * 24.0 * 365.0;
	secondsPerkyr = secondsPerYear * 1.0e3;
  geo_time = secondsPerkyr; /* kyr */

  geo_length = 1.0e3; /* km */
  
  switch (units) {
    case UNITS_NATIVE:
      Ls = 1.0;
      Vs = 1.0;
      Es = 1.0;
      break;
    case UNITS_SI:
      break;

    case UNITS_GEO: /* convert length, vel and visc into SI */
      /* convert input length to SI (km -> m) */
      Ls = Ls * geo_length;

      /* convert input velocty to SI (cm/yr -> m/s) */
      Vs = Vs * 1.0e-2 / geo_time;
      break;
      
    default:
      break;
  }
  
  ierr = ptatinQuantitySetupSIUnitList_UsingLVEScales(Ls,Vs,Es);CHKERRQ(ierr);
  if (units != UNITS_NATIVE) {
    ierr = ptatinQuantitySetupGeoUnitList();CHKERRQ(ierr);
  }

  n = (PetscInt)QNull;
  PetscPrintf(PETSC_COMM_WORLD,"[ptatin Quantities (SI units)]\n");
  for (i=0; i<n; i++) {
    ierr = QuantityView(PTATIN_SI_QLIST[i]);CHKERRQ(ierr);
  }
  
  PetscPrintf(PETSC_COMM_WORLD,"[ptatin Quantities (geo units)]\n");
  for (i=0; i<n; i++) {
    ierr = QuantityView(PTATIN_GEO_QLIST[i]);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ptatinQuantityCreate"
PetscErrorCode ptatinQuantityCreate(void)
{
  PetscInt         i,n;
  Quantity         *si_list,*geo_list;
  ptatinInputUnits value;
  PetscErrorCode   ierr;
  
  n = (PetscInt)QNull;
  ierr = PetscMalloc1(n,&si_list);CHKERRQ(ierr);
  ierr = PetscMalloc1(n,&geo_list);CHKERRQ(ierr);
  for (i=0; i<n; i++) {
    si_list[i]  = NULL;
    geo_list[i] = NULL;
  }

  for (i=0; i<n; i++) {
    ierr = QuantityCreate(&si_list[i]);CHKERRQ(ierr);
    ierr = QuantitySetValues(si_list[i],QuantityTypeNames[i],QuantityTypeUnits[i],1.0);CHKERRQ(ierr);

    ierr = QuantityCreate(&geo_list[i]);CHKERRQ(ierr);
    ierr = QuantitySetValues(geo_list[i],QuantityTypeNames[i],QuantityTypeUnits[i],1.0);CHKERRQ(ierr);
  }
  
  PTATIN_SI_QLIST = si_list;
  PTATIN_GEO_QLIST = geo_list;
  
  value = UNITS_NATIVE;
  ierr = PetscOptionsGetEnum(NULL,NULL,"-ptatin_input_units",ptatinInputUnitsNames,(PetscEnum*)&value,NULL);CHKERRQ(ierr);
  if (value == UNITS_NATIVE) PetscPrintf(PETSC_COMM_WORLD,"pTatin: Using \"native\" units for input parameters\n");
  if (value == UNITS_SI) PetscPrintf(PETSC_COMM_WORLD,"pTatin: Using \"si\" units for input parameters\n");
  if (value == UNITS_GEO) PetscPrintf(PETSC_COMM_WORLD,"pTatin: Using \"geo\" units for input parameters\n");

  ierr = ptatinQuantitySetup_UsingLVEScales(value);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ptatinQuantityDestroy"
PetscErrorCode ptatinQuantityDestroy(void)
{
  PetscInt i,n;

  n = (PetscInt)QNull;
  if (PTATIN_SI_QLIST) {
    for (i=0; i<n; i++) {
      if (PTATIN_SI_QLIST[i]) { PetscFree(PTATIN_SI_QLIST[i]); }
    }
    PetscFree(PTATIN_SI_QLIST);
    PTATIN_SI_QLIST = NULL;
  }
  if (PTATIN_GEO_QLIST) {
    for (i=0; i<n; i++) {
      if (PTATIN_GEO_QLIST[i]) { PetscFree(PTATIN_GEO_QLIST[i]); }
    }
    PetscFree(PTATIN_GEO_QLIST);
    PTATIN_GEO_QLIST = NULL;
  }
  
  PetscFunctionReturn(0);
}


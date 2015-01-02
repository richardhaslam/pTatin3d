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
 **    filename:   units.h
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

#ifndef _ptatin_units_h__
#define _ptatin_units_h__

typedef struct _p_Units *Units;

struct _p_Units {
	char   *quantity_name;
	char   *unit_name;
	double characteristic_scale;
};


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

#endif


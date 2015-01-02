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
 **    filename:   phase_map.h
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

#ifndef __PHASE_MAP_H__
#define __PHASE_MAP_H__

typedef enum { PHASE_MAP_POINT_OUTSIDE=-1, PHASE_MAP_POINT_INSIDE=1 } PhaseMapLocationIndicator;

typedef struct _p_PhaseMap *PhaseMap;
struct _p_PhaseMap {
	double x0,y0,x1,y1;
	double dx,dy;
	int mx,my;
	int nphases;
	int *data;
};

void PhaseMapCreate(PhaseMap *map);
void PhaseMapDestroy(PhaseMap *map);
void PhaseMapGetIndex(PhaseMap pm,const int i,const int j, int *index);
void PhaseMapLoadFromFile(const char filename[],PhaseMap *map);
void PhaseMapGetPhaseIndex(PhaseMap phasemap,double xp[],int *phase);
void PhaseMapCheckValidity(PhaseMap phasemap,int phase,int *is_valid);
void PhaseMapGetMaxPhases(PhaseMap phasemap,int *maxphase);
void PhaseMapViewGnuplot(const char filename[],PhaseMap phasemap);

PetscErrorCode pTatinCtxAttachPhaseMap(pTatinCtx ctx,PhaseMap map);
PetscErrorCode pTatinCtxGetPhaseMap(pTatinCtx ctx,PhaseMap *map);


#endif


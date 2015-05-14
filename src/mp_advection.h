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
 **    filename:   mp_advection.h
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

#ifndef __ptatin_mp_advection_h__
#define __ptatin_mp_advection_h__

typedef enum { RK_ORDER_1=1, RK_ORDER_2, RK_ORDER_3, RK_ORDER_4 } RKOrder;

PetscErrorCode SwarmUpdatePosition_MPntStd_Euler(DM da,Vec velocity,PetscReal step,int npoints,MPntStd marker[]);
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep(DM da,Vec velocity,PetscReal *step);
PetscErrorCode SwarmUpdateProperties_MPntStd(DataBucket db,pTatinCtx ctx,Vec X);
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates(DataBucket materialpoints,DM dav,Vec velocity,PetscReal dt);
PetscErrorCode MaterialPointStd_UpdateLocalCoordinates(DataBucket materialpoints,DM dav);
PetscErrorCode MaterialPointStd_Removal(DataBucket materialpoints);
PetscErrorCode SwarmUpdatePosition_Communication_Generic(DataBucket db,DM da,DataEx de);
PetscErrorCode MaterialPointStd_UpdateCoordinates(DataBucket materialpoints,DM dav,DataEx de);

PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates_RK2(DataBucket materialpoints,DataEx de,DM dav,Vec velocity,PetscReal dt);
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates_RK3(DataBucket materialpoints,DataEx de,DM dav,Vec velocity,PetscReal dt);
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates_RK4(DataBucket materialpoints,DataEx de,DM dav,Vec velocity,PetscReal dt);
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinatesRK(DataBucket materialpoints,RKOrder order,DataEx de,DM dav,Vec velocity,PetscReal dt);

#endif


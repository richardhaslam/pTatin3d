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
 **    Filename:      mp_advection.h
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

#ifndef __ptatin_mp_advection_h__
#define __ptatin_mp_advection_h__

PetscErrorCode SwarmUpdatePosition_MPntStd_Euler(DM da,Vec velocity,PetscReal step,int npoints,MPntStd marker[]);
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep(DM da,Vec velocity,PetscReal *step);
PetscErrorCode SwarmUpdateProperties_MPntStd(DataBucket db,pTatinCtx ctx,Vec X);
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates(DataBucket materialpoints,DM dav,Vec velocity,PetscReal dt);
PetscErrorCode MaterialPointStd_UpdateLocalCoordinates(DataBucket materialpoints,DM dav);
PetscErrorCode MaterialPointStd_Removal(DataBucket materialpoints);
PetscErrorCode SwarmUpdatePosition_Communication_Generic(DataBucket db,DM da,DataEx de);
PetscErrorCode MaterialPointStd_UpdateCoordinates(DataBucket materialpoints,DM dav,DataEx de);


#endif

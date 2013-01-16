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
 **    Filename:      material_point_std_utils.h
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

#ifndef __MATERIAL_POINT_STD_UTILS_H__
#define __MATERIAL_POINT_STD_UTILS_H__

#include "swarm_fields.h"
#include "data_exchanger.h"
#include "MPntStd_def.h"


PetscErrorCode SwarmMPntStd_AssignUniquePointIdentifiers(MPI_Comm comm,DataBucket db,int start_pid,int end_pid);
PetscErrorCode SwarmMPntStd_CoordAssignment_LatticeLayout3d(DM da,PetscInt Nxp[],PetscReal perturb,DataBucket db);
PetscErrorCode SwarmMPntStd_CoordAssignment_RandomLayout3d(DM da,PetscInt nPerCell,DataBucket db);
PetscErrorCode SwarmOutputParaView_MPntStd(DataBucket db,const char path[],const char prefix[]);

#endif


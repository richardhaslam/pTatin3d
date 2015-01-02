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
 **    filename:   material_point_load.h
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



#ifndef __ptatin_material_point_load_h__
#define __ptatin_material_point_load_h__

PetscErrorCode MarkerCoordinatesLoadFromFile(const char name[],long int *length,double **coords);
PetscErrorCode MarkerScalarFieldLoadFromFile(const char name[],long int *length,void **field);

PetscErrorCode MaterialPointStdRemoval(DataBucket db,long int start,long int npoints,const int wil_key);
PetscErrorCode MaterialPointStdInsertBasic(DataBucket db,DM da,long int start,long int npoints,double coords_mp[],int phase_mp[]);
PetscErrorCode MaterialPointDataBasicLoadIntoListFromFile(DataBucket db,DM da,PetscBool append,const char coordfile[],const char phasefile[]);

PetscErrorCode SwarmDataWriteToPetscVec(DataBucket db,const char suffix[]);
PetscErrorCode SwarmDataLoadFromPetscVec(DataBucket db,const char suffix[]);

PetscErrorCode SwarmDataWriteToPetscVec_MPntStd(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz);
PetscErrorCode SwarmDataWriteToPetscVec_MPntPStokes(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz);
PetscErrorCode SwarmDataWriteToPetscVec_MPntPStokesPl(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz);
PetscErrorCode SwarmDataWriteToPetscVec_MPntPEnergy(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz);

PetscErrorCode SwarmDataLoadFromPetscVec_MPntStd(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz);
PetscErrorCode SwarmDataLoadFromPetscVec_MPntPStokes(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz);
PetscErrorCode SwarmDataLoadFromPetscVec_MPntPStokesPl(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz);
PetscErrorCode SwarmDataLoadFromPetscVec_MPntPEnergy(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz);

#endif


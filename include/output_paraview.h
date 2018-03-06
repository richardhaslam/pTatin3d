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
 **    filename:   output_paraview.h
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

#ifndef __ptatin_output_paraview_h__
#define __ptatin_output_paraview_h__

#include "petsc.h"
#include "petscdm.h"

PetscErrorCode pTatinGenerateParallelVTKName(const char prefix[],const char suffix[],char **name);
PetscErrorCode pTatinGenerateVTKName(const char prefix[],const char suffix[],char **name);

PetscErrorCode ParaviewPVDOpen(const char pvdfilename[]);
PetscErrorCode ParaviewPVDAppend(const char pvdfilename[],double time,const char datafile[], const char DirectoryName[]);

PetscErrorCode pTatinOutputParaViewMeshVelocityPressure(DM pack,Vec X,const char path[],const char prefix[]);
PetscErrorCode pTatinOutputMeshVelocityPressureVTS_v0(DM pack,Vec X,const char name[]);
PetscErrorCode pTatinOutputMeshVelocityPressureVTS_v0_binary(DM pack,Vec X,const char name[]);
PetscErrorCode pTatinOutputMeshVelocityPressureVTS_v0_binary_gz(DM pack,Vec X,const char name[]);
PetscErrorCode pTatinOutputMeshVelocityPressurePVTS(DM pack,const char prefix[],const char name[]);

PetscErrorCode pTatinOutputLiteParaViewMeshVelocity(DM pack,Vec X,const char path[],const char prefix[]);
PetscErrorCode pTatinOutputLiteMeshVelocityVTS_v0_binary(DM pack,Vec X,const char name[]);
PetscErrorCode pTatinOutputLiteMeshVelocityPVTS(DM pack,const char prefix[],const char name[]);

PetscErrorCode DAQ2PieceExtendForGhostLevelZero( FILE *vtk_fp, int indent_level, DM dau, const char local_file_prefix[] );
PetscErrorCode pTatinOutputLiteMeshVelocitySlicedPVTS(DM pack,const char path[],const char prefix[]);

#endif


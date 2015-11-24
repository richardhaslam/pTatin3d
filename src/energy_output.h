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
 **    filename:   energy_output.h
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

#ifndef __energy_output_h__
#define __energy_output_h__

PetscErrorCode pTatin3d_ModelOutput_Temperature_Energy(pTatinCtx ctx,Vec X,const char prefix[]);
PetscErrorCode pTatinOutputParaViewMeshEnergy(DM daT,Quadrature Q,Vec X,const char path[],const char prefix[]);
PetscErrorCode pTatinOutputMeshEnergyPVTS(DM daT,Quadrature Q,const char prefix[],const char name[]);
PetscErrorCode DAQ1PieceExtendForGhostLevelZero(FILE *vtk_fp,int indent_level,DM da,const char local_file_prefix[]);

#endif

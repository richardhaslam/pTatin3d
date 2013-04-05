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
 **    Filename:      ptatin_log.h
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

#ifndef __ptatin_log_h__
#define __ptatin_log_h__

PetscErrorCode pTatinLogOpenFile(pTatinCtx ctx);
PetscErrorCode pTatinLogCloseFile(pTatinCtx ctx);
PetscErrorCode pTatinLogHeader(pTatinCtx ctx);

PetscErrorCode pTatinLogBasic(pTatinCtx ctx);
PetscErrorCode pTatinLogBasicKSP(pTatinCtx ctx,const char kspname[],KSP ksp);
PetscErrorCode pTatinLogBasicSNES(pTatinCtx ctx,const char snesname[],SNES snes);
PetscErrorCode pTatinLogBasicStokesSolution(pTatinCtx ctx,DM pack,Vec X);
PetscErrorCode pTatinLogBasicStokesSolutionResiduals(pTatinCtx ctx,SNES snes,DM pack,Vec X);
PetscErrorCode pTatinLogBasicDMDA(pTatinCtx ctx,const char dmname[],DM dm);
PetscErrorCode pTatinLogBasicMaterialPoints(pTatinCtx ctx,const char mpname[],DataBucket db);
PetscErrorCode pTatinLogBasicCPUtime(pTatinCtx ctx,const char component_description[],double time);
PetscErrorCode pTatinLogNote(pTatinCtx ctx,const char comment[]);
PetscErrorCode pTatinLogNote2(pTatinCtx ctx,const char comment1[],const char comment2[]);

#endif


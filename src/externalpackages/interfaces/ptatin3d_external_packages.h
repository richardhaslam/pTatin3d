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
 **    filename:   ptatin3d_external_packages.h
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

#ifndef __ptatin_external_packages_h__
#define __ptatin_external_packages_h__

PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_SPMA(pTatinCtx pctx,Vec X);


PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_FastScape_V3(pTatinCtx pctx,Vec X,
                                                                   PetscInt refinement_factor,
                                                                   PetscReal Lstar,PetscReal Vstar,
                                                                   PetscReal dt_mechanical,PetscReal dt_spm,
                                                                   PetscInt _law,PetscReal _m,PetscReal _kf,PetscReal _kd,PetscInt _bc);


#endif


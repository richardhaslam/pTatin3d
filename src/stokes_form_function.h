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
 **    Filename:      stokes_form_function.h
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

#ifndef __ptatin_stokes_form_function_h__
#define __ptatin_stokes_form_function_h__

PetscErrorCode FormFunction_Stokes(SNES snes,Vec X,Vec F,void *ctx);
PetscErrorCode FormFunction_Stokes_QuasiNewtonX(SNES snes,Vec X,Vec F,void *ctx);
PetscErrorCode MF_Stokes(Vec X,Vec Y,void *ctx);

PetscErrorCode MatStokesJijGetContext(Mat J,void **data);
PetscErrorCode MatStokesJijUpdateGlobalFields(Mat J,Vec u,Vec p,Vec x);
PetscErrorCode MatStokesJijUpdateLocalFields(Mat J,Vec u,Vec p,Vec x);

PetscErrorCode MatCreateStokesJux(pTatinCtx ctx,void *Jctx,Mat *_Jux);
PetscErrorCode MatCreateStokesJpx(pTatinCtx ctx,void *Jctx,Mat *_Jpx);
PetscErrorCode MatCreateStokesJuu(pTatinCtx ctx,void *Jctx,Mat *_Juu);

#endif


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
 **    filename:   pswarm_impl.h
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

#ifndef __PSWARM_IMPL_H__
#define __PSWARM_IMPL_H__

#include <petsc.h>
#include <petscviewer.h>
#include <petscdm.h>
#include <petsc-private/petscimpl.h>

#include <pswarm.h>
#include <data_bucket.h>
#include <data_exchanger.h>

typedef enum { PSW_TS_UNINIT=0, PSW_TS_STALE, PSW_TS_INSYNC } PSwarmStateType;

typedef struct _PSwarmOps *PSwarmOps;
struct _PSwarmOps {
    PetscErrorCode (*view)(PSwarm,PetscViewer);
    PetscErrorCode (*view_pv)(PSwarm,const char*);
    PetscErrorCode (*advect)(PSwarm,DM,Vec);
    PetscErrorCode (*field_update_finitestrain)(PSwarm,DM,Vec);
    PetscErrorCode (*field_update_ptt)(PSwarm,DM,DM,Vec,Vec,PetscReal);
};

struct _p_PSwarm {
    PETSCHEADER(struct _PSwarmOps);
    DataBucket db;
    DataEx     de;
    pTatinCtx  pctx;
    PSwarmStateType state;
    PetscBool setup;
    PetscBool db_set_by_user,de_set_by_user;
};

#endif

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
 **    filename:   pswarm.h
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

#ifndef __PSWARM_H__
#define __PSWARM_H__

PETSC_EXTERN PetscClassId PSWARM_CLASSID;

typedef struct _p_PSwarm* PSwarm;

typedef enum { PSWARM_FU_NULL=0,
    PSWARM_FU_ADVECT,
    PSWARM_FU_STRESS,
    PSWARM_FU_STRAINRATE,
    PSWARM_FU_FINITESTRAIN,
    PSWARM_FU_PTT,
    PSWARM_FU_Pressure
} PSwarmFieldUpdateType;

typedef enum { PSWARM_TM_EULERIAN=0, PSWARM_TM_LAGRANGIAN } PSwarmTransportModeType;
typedef enum { PSWARM_ADV_RK1=0, PSWARM_ADV_RK2, PSWARM_ADV_RK4 } PSwarmAdvectopmRKType;


PetscErrorCode PSwarmInitializePackage(void);

PetscErrorCode PSwarmSetOptionsPrefix(PSwarm ps,const char prefix[]);
PetscErrorCode PSwarmCreate(MPI_Comm comm,PSwarm *ps);
PetscErrorCode PSwarmView(PSwarm ps);
PetscErrorCode PSwarmDestroy(PSwarm *ps);
PetscErrorCode PSwarmSetUp(PSwarm ps);
PetscErrorCode PSwarmSetFromOptions(PSwarm ps);
PetscErrorCode PSwarmFieldUpdateAll(PSwarm ps);
PetscErrorCode PSwarmCreateFromPtatinCtx(pTatinCtx pctx,PSwarm *ps);
PetscErrorCode PSwarmDefineCommTopologyFromDMDA(PSwarm ps,DM dm);
PetscErrorCode PSwarmCreateMultipleInstances(MPI_Comm comm,PSwarm **pslist);
PetscErrorCode PSwarmViewInfo(PSwarm ps);
PetscErrorCode PSwarmSetRegionIndex(PSwarm ps,PetscInt ridx);

PetscErrorCode PSwarmSetPtatinCtx(PSwarm ps,pTatinCtx pctx);
PetscErrorCode PSwarmSetDataBucket(PSwarm ps,DataBucket db);
PetscErrorCode PSwarmSetDataExchanger(PSwarm ps,DataEx de);
PetscErrorCode PSwarmSetTransportModeType(PSwarm ps,PSwarmTransportModeType type);
PetscErrorCode PSwarmSetFieldUpdateType(PSwarm ps,PSwarmFieldUpdateType type);

PetscErrorCode PSwarmAttachStateVecVelocityPressure(PSwarm ps,Vec x);
PetscErrorCode PSwarmAttachStateVecTemperature(PSwarm ps,Vec x);

PetscErrorCode PSwarmGetDataBucket(PSwarm ps,DataBucket *db);

PetscErrorCode _PSwarmFieldUpdate_Advect(PSwarm ps,DM dmv,Vec v);
PetscErrorCode _PSwarmFieldUpdate_FiniteStrain(PSwarm ps,DM dmv,Vec v);
PetscErrorCode _PSwarmFieldUpdate_PressTempTime(PSwarm ps,DM dmp,DM dmt,Vec p,Vec t,PetscReal time);

PetscErrorCode PSwarmFieldUpdate_Advect(PSwarm ps);
PetscErrorCode PSwarmFieldUpdate_FiniteStrain(PSwarm ps);
PetscErrorCode PSwarmFieldUpdate_PressTempTime(PSwarm ps);

#endif

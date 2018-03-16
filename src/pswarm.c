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
 **    filename:   pswarm.c
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
#include <ptatin3d.h>
#include <private/ptatin_impl.h>
#include <MPntStd_def.h>
#include <material_point_std_utils.h>
#include <dmda_element_q2p1.h>
#include <mp_advection.h>
#include <ptatin3d.h>
#include <ptatin_utils.h>
#include <ptatin3d_stokes.h>
#include <ptatin3d_energy.h>
#include <pswarm.h>
#include <private/pswarm_impl.h>
#include <mpiio_blocking.h>
#include <model_utils.h>
#include <element_utils_q2.h>


PetscClassId PSWARM_CLASSID;

const char PSWARM_COMPOSED_STATE_VELPRES[] = "PSWarmStateVector_VP";
const char PSWARM_COMPOSED_STATE_TEMP[]    = "PSWarmStateVector_T";

PetscErrorCode SwarmDMDA3dDataExchangerCreate(DM da,DataEx *_de);

PetscErrorCode PSwarmView(PSwarm ps,PSwarmViewType type);
PetscErrorCode PSwarmDestroy(PSwarm *ps);
PetscErrorCode _PSwarmFieldUpdate_AdvectEulerian(PSwarm ps,DM dmv,Vec v);
PetscErrorCode _PSwarmFieldUpdate_AdvectLagrangian(PSwarm ps,DM dmv,Vec v);

/*
 
 PSwarmCreate()
 PSwarmSetOptionsPrefix()
 PSwarmSetPtatinCtx()
 PSwarmAttachStateVecVelocityPressure()
 PSwarmAttachStateVecTemperature()
 ..
 PSwarmDestroy()

 PSwarmCreateFromPtatinCtx()
 PSwarmSetOptionsPrefix()
 PSwarmAttachStateVecVelocityPressure()
 PSwarmAttachStateVecTemperature()
 
 PSwarmDestroy()

 
*/

#undef __FUNCT__
#define __FUNCT__ "PSwarmInitializePackage"
PetscErrorCode PSwarmInitializePackage(void)
{
    PetscErrorCode   ierr;
    static PetscBool pswarm_registered = PETSC_FALSE;

    PetscFunctionBegin;
    if (!pswarm_registered) {
        ierr = PetscClassIdRegister("Particle Swarm Mangement",&PSWARM_CLASSID);CHKERRQ(ierr);
        pswarm_registered = PETSC_TRUE;
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetOptionsPrefix"
PetscErrorCode PSwarmSetOptionsPrefix(PSwarm ps,const char prefix[])
{
    PetscObjectSetOptionsPrefix((PetscObject)ps,prefix);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmCreate"
PetscErrorCode PSwarmCreate(MPI_Comm comm,PSwarm *ps)
{
    PetscErrorCode ierr;
    PSwarm         p;

    PetscFunctionBegin;
    PetscValidPointer(ps,2);
    *ps = NULL;

    ierr = PSwarmInitializePackage();CHKERRQ(ierr);
    
    ierr = PetscHeaderCreate(p, PSWARM_CLASSID, "PSwarm", "Particle Swarm Manager", "PSwarm", comm, PSwarmDestroy, PSwarmView);CHKERRQ(ierr);
    ierr = PetscMemzero(p->ops, sizeof(struct _PSwarmOps));CHKERRQ(ierr);

    p->state = PSW_TS_UNINIT;

    DataBucketCreate(&p->db);
    DataBucketRegisterField(p->db,MPntStd_classname,sizeof(MPntStd),NULL);
    DataBucketFinalize(p->db);
    
    ierr = PSwarmSetTransportModeType(p,PSWARM_TM_LAGRANGIAN);CHKERRQ(ierr);
    
    p->setup = PETSC_FALSE;
    p->db_set_by_user = PETSC_FALSE;
    p->de_set_by_user = PETSC_FALSE;
    p->transport_mode = PSWARM_TM_LAGRANGIAN;
    p->advection_type = PSWARM_ADV_RK1;
    p->pvdopen = PETSC_FALSE;
  
    *ps = p;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmCreateFromPtatinCtx"
PetscErrorCode PSwarmCreateFromPtatinCtx(pTatinCtx pctx,PSwarm *ps)
{
    PetscErrorCode ierr;
    MPI_Comm       comm;
    
    PetscFunctionBegin;
    ierr = PetscObjectGetComm((PetscObject)pctx->pack,&comm);CHKERRQ(ierr);
    ierr = PSwarmCreate(comm,ps);CHKERRQ(ierr);
    ierr = PSwarmSetPtatinCtx(*ps,pctx);CHKERRQ(ierr);
    ierr = PSwarmSetDataBucket(*ps,pctx->materialpoint_db);CHKERRQ(ierr);
    ierr = PSwarmSetDataExchanger(*ps,pctx->materialpoint_ex);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetPtatinCtx"
PetscErrorCode PSwarmSetPtatinCtx(PSwarm ps,pTatinCtx pctx)
{
    ps->pctx = pctx;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetDataBucket"
PetscErrorCode PSwarmSetDataBucket(PSwarm ps,DataBucket db)
{
    if (ps->db && !ps->db_set_by_user) {
        DataBucketDestroy(&ps->db);
    }
    ps->db = db;
    ps->db_set_by_user = PETSC_TRUE;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmGetDataBucket"
PetscErrorCode PSwarmGetDataBucket(PSwarm ps,DataBucket *db)
{
    if (db) { *db = ps->db; }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetDataExchanger"
PetscErrorCode PSwarmSetDataExchanger(PSwarm ps,DataEx de)
{
    PetscErrorCode ierr;
    if (ps->de && !ps->de_set_by_user) {
        ierr = DataExDestroy(ps->de);CHKERRQ(ierr);
    }
    ps->de = de;
    ps->de_set_by_user = PETSC_TRUE;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmDefineCommTopologyFromDMDA"
PetscErrorCode PSwarmDefineCommTopologyFromDMDA(PSwarm ps,DM dm)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = SwarmDMDA3dDataExchangerCreate(dm,&ps->de);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmAttachStateVecVelocityPressure"
PetscErrorCode PSwarmAttachStateVecVelocityPressure(PSwarm ps,Vec x)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    if (x) {
        Vec Xtmp;
      
        ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_VELPRES,(PetscObject*)&Xtmp);CHKERRQ(ierr);
        if (Xtmp) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector X=(u,p) already attached to PSwarm");

        ierr = PetscObjectCompose((PetscObject)ps,PSWARM_COMPOSED_STATE_VELPRES,(PetscObject)x);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmAttachStateVecTemperature"
PetscErrorCode PSwarmAttachStateVecTemperature(PSwarm ps,Vec x)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    if (x) {
        Vec Xtmp;
      
        ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_TEMP,(PetscObject*)&Xtmp);CHKERRQ(ierr);
        if (Xtmp) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector T already attached to PSwarm");

        ierr = PetscObjectCompose((PetscObject)ps,PSWARM_COMPOSED_STATE_TEMP,(PetscObject)x);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmViewInfo"
PetscErrorCode PSwarmViewInfo(PSwarm ps)
{
  const char *prefix;
  int n_points;
  MPI_Comm comm;
  PetscErrorCode ierr;
  
  comm = PetscObjectComm((PetscObject)ps);
  DataBucketGetSizes(ps->db,&n_points,NULL,NULL);

  ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);
  if (prefix) {
    PetscPrintf(comm,"PSwarm [%s] --------------------\n",prefix);
  } else {
    PetscPrintf(comm,"PSwarm [%s] --------------------\n",prefix);
  }
  PetscPrintf(comm,"  npoints %D\n",n_points);
  PetscPrintf(comm,"  Transport mode: \n");
  if (ps->transport_mode == PSWARM_TM_EULERIAN) {
    PetscPrintf(comm,"    Eulerian\n");
  }
  if (ps->transport_mode == PSWARM_TM_LAGRANGIAN) {
    PetscPrintf(comm,"    Lagrangian\n");
  }
  PetscPrintf(comm,"  Update methods: \n");
  if (ps->ops->field_update_finitestrain) {
    PetscPrintf(comm,"    Finite strain\n");
  }
  if (ps->ops->field_update_ptt) {
    PetscPrintf(comm,"    Pressure-Temperature-Time\n");
  }
  if (ps->ops->field_update_pressure) {
    PetscPrintf(comm,"    Pressure\n");
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmDestroy"
PetscErrorCode PSwarmDestroy(PSwarm *ps)
{
    PetscErrorCode ierr;
    PSwarm         p;
    
    PetscFunctionBegin;
    if (!ps) PetscFunctionReturn(0);
    p = *ps;
    if (!p) PetscFunctionReturn(0);
  
    if (!p->db_set_by_user) {
        if (p->db) { DataBucketDestroy(&p->db); }
    }
    if (!p->de_set_by_user) {
        if (p->de) { ierr = DataExDestroy(p->de);CHKERRQ(ierr); }
    }

    PetscHeaderDestroy(ps);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetTransportModeType"
PetscErrorCode PSwarmSetTransportModeType(PSwarm ps,PSwarmTransportModeType type)
{
    PetscFunctionBegin;
    ps->transport_mode = type;
    switch (type) {
        case PSWARM_TM_EULERIAN:
            ps->ops->advect = NULL;//_PSwarmFieldUpdate_AdvectEulerian;
            break;
        case PSWARM_TM_LAGRANGIAN:
            ps->ops->advect = _PSwarmFieldUpdate_AdvectLagrangian;
            break;
            
        default:
            ps->ops->advect = NULL;
            break;
    }
    
    PetscFunctionReturn(0);
}

/* Pressure update functionality */
#undef __FUNCT__
#define __FUNCT__ "PSwarmUpdate_Pressure"
PetscErrorCode PSwarmUpdate_Pressure(PSwarm ps,DM dmv,DM dmp,Vec pressure)
{
  PetscErrorCode ierr;
  DataField datafield_tracers,datafield;
  MPntStd *tracer;
  double *tracer_pressure;
  int n_tracers,p;
  PetscReal elcoords[3*Q2_NODES_PER_EL_3D],elp[P_BASIS_FUNCTIONS];
  Vec gcoords,pressure_l;
  const PetscScalar *LA_gcoords;
  const PetscScalar *LA_pfield;
  PetscInt nel,nen_u,nen_p,k;
  const PetscInt *elnidx_u,*elnidx_p;
  
  if (ps->state == PSW_TS_STALE) {
    /* update local coordinates and perform communication */
		ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
    ps->state = PSW_TS_INSYNC;
  }
  
  DataBucketGetDataFieldByName(ps->db,MPntStd_classname,&datafield_tracers);
	DataFieldGetEntries(datafield_tracers,(void**)&tracer);

  DataBucketGetDataFieldByName(ps->db,"pressure",&datafield);
	DataFieldGetEntries(datafield,(void**)&tracer_pressure);
  
  ierr = DMDAGetElements_pTatinQ2P1(dmv,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dmv,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);

  ierr = DMDAGetElements_pTatinQ2P1(dmp,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dmp,&pressure_l);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dmp,pressure,INSERT_VALUES,pressure_l);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dmp,pressure,INSERT_VALUES,pressure_l);CHKERRQ(ierr);
  ierr = VecGetArrayRead(pressure_l,&LA_pfield);CHKERRQ(ierr);
  
	DataBucketGetSizes(ps->db,&n_tracers,NULL,NULL);
  for (p=0; p<n_tracers; p++) {
    double        *xi_p;
    PetscReal     NIp[P_BASIS_FUNCTIONS],pressure_p;
    PetscInt      eidx;
    
    xi_p = tracer[p].xi;
    eidx = (PetscInt)tracer[p].wil;

    ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*eidx],(PetscScalar*)LA_gcoords);CHKERRQ(ierr);
    ierr = DMDAGetScalarElementField(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*eidx],(PetscScalar*)LA_pfield);CHKERRQ(ierr);
    
    ConstructNi_pressure(xi_p,elcoords,NIp);

    pressure_p = 0.0;
    for (k=0; k<P_BASIS_FUNCTIONS; k++) {
      pressure_p += NIp[k] * elp[k];
    }
    tracer_pressure[p] = pressure_p;
  }
  ierr = VecRestoreArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(pressure_l,&LA_pfield);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dmp,&pressure_l);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetFieldType_Pressure"
PetscErrorCode PSwarmSetFieldType_Pressure(PSwarm ps)
{
  BTruth found;
  
  DataBucketQueryDataFieldByName(ps->db,"pressure",&found);
  if (!found) {
    DataBucketRegisterField(ps->db,"pressure",sizeof(double),NULL);
    DataBucketFinalize(ps->db);
  }
  ps->ops->field_update_pressure = PSwarmUpdate_Pressure;
  
  PetscFunctionReturn(0);
}


/*
 [1] Register data
 [2] Define any updates for history variables
*/
#undef __FUNCT__
#define __FUNCT__ "PSwarmSetFieldUpdateType"
PetscErrorCode PSwarmSetFieldUpdateType(PSwarm ps,PSwarmFieldUpdateType type)
{
    PetscErrorCode ierr;
  
    PetscFunctionBegin;
    switch (type) {
        case PSWARM_FU_NULL:
            break;
        case PSWARM_FU_ADVECT:
            SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"Set transport mode via PSwarmSetTransportModeType()");
            break;
        case PSWARM_FU_FINITESTRAIN:
            ps->ops->field_update_finitestrain = NULL;
            break;
        case PSWARM_FU_PTT:
            ps->ops->field_update_ptt = NULL;
            break;
        case PSWARM_FU_Pressure:
            ierr = PSwarmSetFieldType_Pressure(ps);CHKERRQ(ierr);
            break;
        
        default:
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_PSwarmFieldUpdate_AdvectEulerian"
PetscErrorCode _PSwarmFieldUpdate_AdvectEulerian(PSwarm ps,DM dmv,Vec v)
{
    ps->state = PSW_TS_INSYNC;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_PSwarmFieldUpdate_AdvectLagrangian"
PetscErrorCode _PSwarmFieldUpdate_AdvectLagrangian(PSwarm ps,DM dmv,Vec v)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = MaterialPointStd_UpdateGlobalCoordinates(ps->db,dmv,v,ps->pctx->dt);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmFieldUpdate_Advect"
PetscErrorCode PSwarmFieldUpdate_Advect(PSwarm ps)
{
    PetscErrorCode ierr;
    PhysCompStokes stokes;
    DM             dmv,dmstokes;
    
    PetscFunctionBegin;
    ierr = PSwarmSetUp(ps);CHKERRQ(ierr);

    ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);

    if (ps->state == PSW_TS_STALE) {
        /* update local coordinates and perform communication */
		ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
        ps->state = PSW_TS_INSYNC;
    }
    
    if (ps->ops->advect) {
        Vec X,velocity,pressure;
        
        ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_VELPRES,(PetscObject*)&X);CHKERRQ(ierr);
        if (!X) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector X=(u,p) was not provided. User must call PSwarmAttachStateVecVelocityPressure()");
        
        ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
        
        ierr = ps->ops->advect(ps,dmv,velocity);CHKERRQ(ierr);

        ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
        
        /* flag as being stale => local coordinates need updating */
        ps->state = PSW_TS_STALE;
    } else {
        SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"FieldUpdate(Advect) was not activated");
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmFieldUpdate_FiniteStrain"
PetscErrorCode PSwarmFieldUpdate_FiniteStrain(PSwarm ps)
{
    PetscErrorCode ierr;
    PhysCompStokes stokes;
    DM             dmv,dmstokes;
    
    PetscFunctionBegin;
    ierr = PSwarmSetUp(ps);CHKERRQ(ierr);

    ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);

    if (ps->state == PSW_TS_STALE) {
        /* update local coordinates and perform communication */
		ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
        ps->state = PSW_TS_INSYNC;
    }
    
    if (ps->ops->field_update_finitestrain) {
        Vec X,velocity,pressure;
        
        ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_VELPRES,(PetscObject*)&X);CHKERRQ(ierr);
        if (!X) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector X=(u,p) was not provided. User must call PSwarmAttachStateVecVelocityPressure()");

		ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);

        ierr = ps->ops->field_update_finitestrain(ps,dmv,velocity);CHKERRQ(ierr);

		ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
    } else {
        SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"FieldUpdate(FiniteStrain) was not activated");
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmFieldUpdate_PressTempTime"
PetscErrorCode PSwarmFieldUpdate_PressTempTime(PSwarm ps)
{
    PetscErrorCode ierr;
    PhysCompStokes stokes;
    PhysCompEnergy energy;
    DM             dmv,dmp,dmstokes,dmT;

    PetscFunctionBegin;
    ierr = PSwarmSetUp(ps);CHKERRQ(ierr);

    ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMs(stokes,&dmv,&dmp);CHKERRQ(ierr);

    if (ps->state == PSW_TS_STALE) {
        /* update local coordinates and perform communication */
		ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
        ps->state = PSW_TS_INSYNC;
    }

    if (ps->ops->field_update_ptt) {
        Vec X,velocity,pressure,temperature;
        
		ierr = pTatinGetContext_Energy(ps->pctx,&energy);CHKERRQ(ierr);
		dmT = energy->daT;
        ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_VELPRES,(PetscObject*)&X);CHKERRQ(ierr);
        ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_TEMP,(PetscObject*)&temperature);CHKERRQ(ierr);
        if (!X) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector X=(u,p) was not provided. User must call PSwarmAttachStateVecVelocityPressure()");
        if (!temperature) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector T was not provided. User must call PSwarmAttachStateVecTemperature()");

		ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
        
        ierr = ps->ops->field_update_ptt(ps,dmp,dmT,pressure,temperature,ps->pctx->time);CHKERRQ(ierr);
    
		ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
    } else {
        SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"FieldUpdate(PressTempTime) was not activated");
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmFieldUpdateAll"
PetscErrorCode PSwarmFieldUpdateAll(PSwarm ps)
{
    PetscErrorCode ierr;
    PhysCompStokes stokes;
    PhysCompEnergy energy;
    DM             dmstokes,dmv,dmp,dmT;
    Vec            X,velocity,pressure,temperature;
    
    PetscFunctionBegin;
    
    ierr = PSwarmSetUp(ps);CHKERRQ(ierr);
    
    ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMs(stokes,&dmv,&dmp);CHKERRQ(ierr);
    ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_VELPRES,(PetscObject*)&X);CHKERRQ(ierr);
    
    ierr = pTatinGetContext_Energy(ps->pctx,&energy);CHKERRQ(ierr);
    temperature = NULL;
    dmT = NULL;
    if (energy) {
        dmT = energy->daT;
        ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_TEMP,(PetscObject*)&temperature);CHKERRQ(ierr);
    }
  
    if (ps->state == PSW_TS_STALE) {
        ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
        ps->state = PSW_TS_INSYNC;
    }
    
    if (ps->ops->field_update_finitestrain) {
        if (!X) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector X=(u,p) was not provided. User must call PSwarmAttachStateVecVelocityPressure()");
        ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
        ierr = ps->ops->field_update_finitestrain(ps,dmv,velocity);CHKERRQ(ierr);
        ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
    }
    
    if (ps->ops->field_update_ptt) {
        if (!X) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector X=(u,p) was not provided. User must call PSwarmAttachStateVecVelocityPressure()");
        ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
        ierr = ps->ops->field_update_ptt(ps,dmp,dmT,pressure,temperature,ps->pctx->time);CHKERRQ(ierr);
        ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
    }

    if (ps->ops->field_update_pressure) {
      if (!X) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector X=(u,p) was not provided. User must call PSwarmAttachStateVecVelocityPressure()");
      ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
      ierr = ps->ops->field_update_pressure(ps,dmv,dmp,pressure);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
    }

    /* position must always be the last state variable to be updated */
    if (ps->ops->advect) {
        if (!X) SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"State vector X=(u,p) was not provided. User must call PSwarmAttachStateVecVelocityPressure()");
        ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
        ierr = ps->ops->advect(ps,dmv,velocity);CHKERRQ(ierr);
        /* flag as being stale => local coordinates need updating */
        ps->state = PSW_TS_STALE;
        ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
    }
  
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmMPntStd_CoordAssignment_RestrictedLatticeLayout"
PetscErrorCode SwarmMPntStd_CoordAssignment_RestrictedLatticeLayout(DataBucket db,DM da,PetscReal xmin[],PetscReal xmax[],PetscInt Nxp[],PetscReal perturb)
{
  DataField    PField;
  PetscInt     e;
  Vec          gcoords;
  PetscScalar  *LA_coords;
  PetscScalar  el_coords[Q2_NODES_PER_EL_3D*NSD];
  int          np_cur;
  PetscInt     nel,nen;
  const PetscInt     *elnidx;
  PetscInt     p,k,pi,pj,pk;
  PetscReal    dxi,deta,dzeta;
  long int     np_local;
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  
  ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
  
  //DataBucketSetSizes(db,100,-1);
  DataBucketSetInitialSizes(db,100,-1);
  DataBucketGetSizes(db,&np_cur,NULL,NULL);
  
  if (perturb < 0.0) {
    SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_USER,"Cannot use a negative perturbation");
  }
  if (perturb > 1.0) {
    SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_USER,"Cannot use a perturbation greater than 1.0");
  }
  
  /* setup for coords */
  ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
  
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField);
  DataFieldGetAccess(PField);
  
  dxi    = 2.0/(PetscReal)Nxp[0];
  deta   = 2.0/(PetscReal)Nxp[1];
  dzeta  = 2.0/(PetscReal)Nxp[2];
  
  p = 0;
  for (e=0; e<nel; e++) {
    /* get coords for the element */
    ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_coords);CHKERRQ(ierr);
    
    for (pk=0; pk<Nxp[2]; pk++) {
      for (pj=0; pj<Nxp[1]; pj++) {
        for (pi=0; pi<Nxp[0]; pi++) {
          MPntStd *marker;
          double xip[NSD],xip_shift[NSD],xip_rand[NSD],xp_rand[NSD],Ni[Q2_NODES_PER_EL_3D];
          
          xip[0] = -1.0 + dxi    * (pi + 0.5);
          xip[1] = -1.0 + deta   * (pj + 0.5);
          xip[2] = -1.0 + dzeta  * (pk + 0.5);
          
          /* random between -0.5 <= shift <= 0.5 */
          xip_shift[0] = 1.0*(rand()/((double)RAND_MAX)) - 0.5;
          xip_shift[1] = 1.0*(rand()/((double)RAND_MAX)) - 0.5;
          xip_shift[2] = 1.0*(rand()/((double)RAND_MAX)) - 0.5;
          
          xip_rand[0] = xip[0] + perturb * dxi    * xip_shift[0];
          xip_rand[1] = xip[1] + perturb * deta   * xip_shift[1];
          xip_rand[2] = xip[2] + perturb * dzeta  * xip_shift[2];
          
          pTatin_ConstructNi_Q2_3D(xip_rand,Ni);
          
          xp_rand[0] = xp_rand[1] = xp_rand[2] = 0.0;
          for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
            xp_rand[0] += Ni[k] * el_coords[NSD*k+0];
            xp_rand[1] += Ni[k] * el_coords[NSD*k+1];
            xp_rand[2] += Ni[k] * el_coords[NSD*k+2];
          }
          
          DataFieldAccessPoint(PField,p,(void**)&marker);
          
          if ( (xp_rand[0] >= xmin[0]) && (xp_rand[0] <= xmax[0]) ) {
            if ( (xp_rand[1] >= xmin[1]) && (xp_rand[1] <= xmax[1]) ) {
              if ( (xp_rand[2] >= xmin[2]) && (xp_rand[2] <= xmax[2]) ) {
                marker->coor[0] = xp_rand[0];
                marker->coor[1] = xp_rand[1];
                marker->coor[2] = xp_rand[2];
                
                marker->xi[0] = xip_rand[0];
                marker->xi[1] = xip_rand[1];
                marker->xi[2] = xip_rand[2];
                
                marker->wil    = e;
                marker->pid    = 0;
                
                p++;
              }
            }
          }
          
          /* Re-size bucket */
          if (p == np_cur) {
            DataBucketSetSizes(db,np_cur+100,-1);
            DataBucketGetSizes(db,&np_cur,NULL,NULL);
          }
          
        }
      }
    }
  }
  DataFieldRestoreAccess(PField);
  ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
  
  np_local = p;
  DataBucketSetSizes(db,np_local,-1);
  ierr = SwarmMPntStd_AssignUniquePointIdentifiers(PetscObjectComm((PetscObject)da),db,0,np_local);CHKERRQ(ierr);
  
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetUpCoords_FillDM"
PetscErrorCode PSwarmSetUpCoords_FillDM(PSwarm ps)
{
    PetscErrorCode ierr;
    PhysCompStokes stokes;
    DM             dmv;
    PetscInt       Nxp[] = {1,1,1}; /* change with -lattice_layout_N{x,y,z} */
    
    PetscFunctionBegin;
    
    ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
    
    ierr = SwarmMPntStd_CoordAssignment_LatticeLayout3d(dmv,Nxp,0.0,ps->db);CHKERRQ(ierr);
  
    ps->state = PSW_TS_INSYNC;
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetUpCoords_FillDMWithinBoundingBox"
PetscErrorCode PSwarmSetUpCoords_FillDMWithinBoundingBox(PSwarm ps)
{
  PetscErrorCode ierr;
  PhysCompStokes stokes;
  DM             dmv;
  PetscInt       nn,Nxp[] = {1,1,1}; /* change with -lattice_layout_N{x,y,z} */
  PetscReal      xmin[3],xmax[3];
  const char     *prefix;
  
  PetscFunctionBegin;
  
  ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
  
  ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);

  xmin[0] = -1.0e32;
  xmin[1] = -1.0e32;
  xmin[2] = -1.0e32;

  xmax[0] = 1.0e32;
  xmax[1] = 1.0e32;
  xmax[2] = 1.0e32;

  nn = 3;
	PetscOptionsGetIntArray(NULL,prefix,"-pswarm_lattice_nx",Nxp,&nn,NULL);
  nn = 3;
  PetscOptionsGetRealArray(NULL,prefix,"-pswarm_lattice_min",xmin,&nn,NULL);
  nn = 3;
  PetscOptionsGetRealArray(NULL,prefix,"-pswarm_lattice_max",xmax,&nn,NULL);

  ierr = SwarmMPntStd_CoordAssignment_RestrictedLatticeLayout(ps->db,dmv,xmin,xmax,Nxp,0.0);CHKERRQ(ierr);
  
  ps->state = PSW_TS_INSYNC;
  
  PetscFunctionReturn(0);
}

/*
 There is a possibility that this routine may create duplicate points, e.g.
 two ranks may define a particle with identical coordinates. 
 In general, for usage with passive swarms, this is likely to not be a problem.
 When using this method to define a deformation mesh, duplicate points may be 
 problematic. As a work around, I explicitly assign ALL points a pid value given
 by i + j*nx + k*nx*ny. In this way, associating a particle coordinate with a 
 point in a structured mesh using pid will be safe.
*/
#undef __FUNCT__
#define __FUNCT__ "PSwarmSetUpCoords_FillBox"
PetscErrorCode PSwarmSetUpCoords_FillBox(PSwarm ps)
{
  PetscErrorCode ierr;
  PhysCompStokes stokes;
  DM             dmv;
  PetscInt       nn,nlist,Nxp[] = {2,2,2}; /* change with -lattice_layout_N{x,y,z} */
  PetscReal      xmin[3],xmax[3],*coorlist;
  const char     *prefix;
  PetscBool      found;
  PetscMPIInt    rank;
  MPI_Comm       comm;
  PetscInt       ii,jj,kk,c;
  PetscReal      dx[3],damin[3],damax[3],elmin[3],elmax[3],coor[3];
  
  PetscFunctionBegin;
  
  ierr = PetscObjectGetComm((PetscObject)ps,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
  
  ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);
  
  xmin[0] = -1.0e32;
  xmin[1] = -1.0e32;
  xmin[2] = -1.0e32;
  
  xmax[0] = 1.0e32;
  xmax[1] = 1.0e32;
  xmax[2] = 1.0e32;
  
  nn = 3;
	PetscOptionsGetIntArray(NULL,prefix,"-pswarm_box_nx",Nxp,&nn,NULL);
  
  if (Nxp[0] <= 1) SETERRQ(comm,PETSC_ERR_USER,"Nxp[0] must be greater than 1, use -pswarm_box_nx");
  if (Nxp[1] <= 1) SETERRQ(comm,PETSC_ERR_USER,"Nxp[1] must be greater than 1, use -pswarm_box_nx");
  if (Nxp[2] <= 1) SETERRQ(comm,PETSC_ERR_USER,"Nxp[2] must be greater than 1, use -pswarm_box_nx");
  
  nn = 3;
  PetscOptionsGetRealArray(NULL,prefix,"-pswarm_box_min",xmin,&nn,&found);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Must specify box min extent via -pswarm_box_min");
  
  nn = 3;
  PetscOptionsGetRealArray(NULL,prefix,"-pswarm_box_max",xmax,&nn,&found);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Must specify box max extent via -pswarm_box_max");

  /* Create array of coordinates */
  /* Two pass: first count, then allocate and fill */
  dx[0] = (xmax[0]-xmin[0])/((PetscReal)Nxp[0]-1);
  dx[1] = (xmax[1]-xmin[1])/((PetscReal)Nxp[1]-1);
  dx[2] = (xmax[2]-xmin[2])/((PetscReal)Nxp[2]-1);

  /*ierr = DMDAGetLocalBoundingBox(dmv,damin,damax);CHKERRQ(ierr);*/
  ierr = DMDAComputeQ2ElementBoundingBox(dmv,elmin,elmax);CHKERRQ(ierr);
  ierr = DMDAComputeQ2LocalBoundingBox(dmv,damin,damax);CHKERRQ(ierr);
  damin[0] -= elmin[0] * 1.0e-6;
  damin[1] -= elmin[1] * 1.0e-6;
  damin[2] -= elmin[2] * 1.0e-6;

  damax[0] += elmax[0] * 1.0e-6;
  damax[1] += elmax[1] * 1.0e-6;
  damax[2] += elmax[2] * 1.0e-6;

  c = 0;
  for (kk=0; kk<Nxp[2]; kk++) {
    for (jj=0; jj<Nxp[1]; jj++) {
      for (ii=0; ii<Nxp[0]; ii++) {
        coor[0] = xmin[0] + ii * dx[0];
        coor[1] = xmin[1] + jj * dx[1];
        coor[2] = xmin[2] + kk * dx[2];
        
        if ( (coor[2] < damin[2]) || (coor[2] > damax[2]) ) continue;
        if ( (coor[1] < damin[1]) || (coor[1] > damax[1]) ) continue;
        if ( (coor[0] < damin[0]) || (coor[0] > damax[0]) ) continue;
        
        c++;
      }
    }
  }

  nlist = c;
  PetscMalloc(sizeof(PetscReal)*3*nlist,&coorlist);

  c = 0;
  for (kk=0; kk<Nxp[2]; kk++) {
    for (jj=0; jj<Nxp[1]; jj++) {
      for (ii=0; ii<Nxp[0]; ii++) {

        coor[0] = xmin[0] + ii * dx[0];
        coor[1] = xmin[1] + jj * dx[1];
        coor[2] = xmin[2] + kk * dx[2];
        
        if ( (coor[2] < damin[2]) || (coor[2] > damax[2]) ) continue;
        if ( (coor[1] < damin[1]) || (coor[1] > damax[1]) ) continue;
        if ( (coor[0] < damin[0]) || (coor[0] > damax[0]) ) continue;

        coorlist[c*3+0] = coor[0];
        coorlist[c*3+1] = coor[1];
        coorlist[c*3+2] = coor[2];
        c++;
      }
    }
  }

  ierr = SwarmMPntStd_CoordAssignment_InsertFromList(ps->db,dmv,nlist,coorlist,0,PETSC_TRUE);CHKERRQ(ierr);
  
  PetscFree(coorlist);
  
  /* Traverse points, examine coordinates, set pid based on xp,yp,zp */
  {
    DataField df;
    int p,npoints;
    MPntStd *points;
    
    DataBucketGetSizes(ps->db,&npoints,NULL,NULL);
    DataBucketGetDataFieldByName(ps->db,MPntStd_classname,&df);
    DataFieldGetEntries(df,(void**)&points);
    for (p=0; p<npoints; p++) {
      long int ii,jj,kk;
      
      ii = (long int)( (1.0e-12*elmin[0] + points[p].coor[0] - xmin[0])/ dx[0]);
      jj = (long int)( (1.0e-12*elmin[1] + points[p].coor[1] - xmin[1])/ dx[1]);
      kk = (long int)( (1.0e-12*elmin[2] + points[p].coor[2] - xmin[2])/ dx[2]);
      
      points[p].pid = ii + jj * Nxp[0] + kk * Nxp[0] * Nxp[1];
    }
  }
  
  if (rank == 0) {
    char filename[PETSC_MAX_PATH_LEN];
    
    PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/deformation_grid_ref.vts",ps->pctx->outputpath);
    ierr = pSwarmParaViewMeshDeformationBaseVTS(xmin,dx,Nxp,filename);CHKERRQ(ierr);
  }
  
  ps->state = PSW_TS_INSYNC;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetUpCoords_FromUserList"
PetscErrorCode PSwarmSetUpCoords_FromUserList(PSwarm ps)
{
  PetscErrorCode ierr;
  PhysCompStokes stokes;
  DM             dmv;
  const char     *prefix;
  PetscBool      found;
  MPI_Comm       comm;
  PetscReal      *coorlist,*coorx,*coory,*coorz;
  PetscInt       i,nlist,nlistsize;
  
  PetscFunctionBegin;
  
  ierr = PetscObjectGetComm((PetscObject)ps,&comm);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);

  PetscOptionsGetInt(NULL,prefix,"-pswarm_coor_n",&nlistsize,&found);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Must specify number of coordinates via -pswarm_coor_n");

  PetscMalloc(sizeof(PetscReal)*3*nlistsize,&coorlist);
  PetscMalloc(sizeof(PetscReal)*nlistsize,&coorx);
  PetscMalloc(sizeof(PetscReal)*nlistsize,&coory);
  PetscMalloc(sizeof(PetscReal)*nlistsize,&coorz);
  
  nlist = nlistsize;
  PetscOptionsGetRealArray(NULL,prefix,"-pswarm_coor_x",coorx,&nlist,&found);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Must specify x coordinates via -pswarm_coor_x");
  if (nlist != nlistsize) SETERRQ1(comm,PETSC_ERR_USER,"Must specify %D x coordinates",nlistsize);
  
  nlist = nlistsize;
  PetscOptionsGetRealArray(NULL,prefix,"-pswarm_coor_y",coory,&nlist,&found);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Must specify y coordinates via -pswarm_coor_y");
  if (nlist != nlistsize) SETERRQ1(comm,PETSC_ERR_USER,"Must specify %D y coordinates",nlistsize);

  nlist = nlistsize;
  PetscOptionsGetRealArray(NULL,prefix,"-pswarm_coor_z",coorz,&nlist,&found);
  if (!found) SETERRQ(comm,PETSC_ERR_USER,"Must specify z coordinates via -pswarm_coor_z");
  if (nlist != nlistsize) SETERRQ1(comm,PETSC_ERR_USER,"Must specify %D z coordinates",nlistsize);

  for (i=0; i<nlistsize; i++) {
    coorlist[3*i+0] = coorx[i];
    coorlist[3*i+1] = coory[i];
    coorlist[3*i+2] = coorz[i];
  }
  
  ierr = SwarmMPntStd_CoordAssignment_InsertFromList(ps->db,dmv,nlistsize,coorlist,0,PETSC_TRUE);CHKERRQ(ierr);

  PetscFree(coorlist);
  PetscFree(coorx);
  PetscFree(coory);
  PetscFree(coorz);
  
  ps->state = PSW_TS_INSYNC;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetUpCoords"
PetscErrorCode PSwarmSetUpCoords(PSwarm ps)
{
  PetscErrorCode ierr;
  const char *prefix;
  PetscInt type;
  MPI_Comm comm;
  PetscFunctionBegin;
    
  ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);
  
  type = 1;
  PetscOptionsGetInt(NULL,prefix,"-pswarm_coord_layout",&type,NULL);

  comm = PetscObjectComm((PetscObject)ps);
  switch (type) {
    case 0:
      PetscPrintf(comm,"[PSwarmSetUpCoords_FillDM]\n");
      ierr = PSwarmSetUpCoords_FillDM(ps);CHKERRQ(ierr);
      break;

    case 1:
      PetscPrintf(comm,"[PSwarmSetUpCoords_FillDMWithinBoundingBox]\n");
      ierr = PSwarmSetUpCoords_FillDMWithinBoundingBox(ps);CHKERRQ(ierr);
      break;

    case 2:
      PetscPrintf(comm,"[PSwarmSetUpCoords_FillBox]\n");
      ierr = PSwarmSetUpCoords_FillBox(ps);CHKERRQ(ierr);
      break;

    case 3:
      PetscPrintf(comm,"[PSwarmSetUpCoords_FromUserList]\n");
      ierr = PSwarmSetUpCoords_FromUserList(ps);CHKERRQ(ierr);
      break;
      
    case 4:
      //ierr = PSwarmSetUpCoords_FillDMWithGeometryObject(ps);CHKERRQ(ierr);
      break;

    default:
      break;
  }

  ps->state = PSW_TS_INSYNC;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetUp"
PetscErrorCode PSwarmSetUp(PSwarm ps)
{
    PetscErrorCode ierr;
    PhysCompStokes stokes;
    DM             dmv;

    PetscFunctionBegin;
    
    if (ps->setup) PetscFunctionReturn(0);
    
    if (!ps->pctx) {
        SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_USER,"Must provide a valid pTatinCtx");
    }
    
    if (!ps->de) {
        ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
        ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
        ierr = PSwarmDefineCommTopologyFromDMDA(ps,dmv);CHKERRQ(ierr);
    }

    ierr = PSwarmSetUpCoords(ps);CHKERRQ(ierr);

    {
        const char *prefix;
        PetscInt   ridx;
        PetscBool  isactive;
      
      
        ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);
        ridx = 0;
        isactive = PETSC_FALSE;
        ierr = PetscOptionsGetInt(NULL,prefix,"-pswarm_region_index",&ridx,&isactive);CHKERRQ(ierr);
        if (isactive) { ierr = PSwarmSetRegionIndex(ps,ridx);CHKERRQ(ierr); }
    }
  
    ps->setup = PETSC_TRUE;
    
    PetscFunctionReturn(0);
}

/* set from options */
#undef __FUNCT__
#define __FUNCT__ "PSwarmSetFromOptions"
PetscErrorCode PSwarmSetFromOptions(PSwarm ps)
{
    PetscErrorCode ierr;
    PetscBool      isactive;
  
    PetscFunctionBegin;
    ierr = PetscObjectOptionsBegin((PetscObject)ps);CHKERRQ(ierr);
    ierr = PetscOptionsHead(PetscOptionsObject,"PSwarm options");CHKERRQ(ierr);
    
    isactive = PETSC_FALSE;
    ierr = PetscOptionsBool("-pswarm_transport_mode_eulerian","Transport mode set to Eulerian","PSwarmSetTransportModeType",isactive,&isactive,0);CHKERRQ(ierr);
    if (isactive) { ierr = PSwarmSetTransportModeType(ps,PSWARM_TM_EULERIAN);CHKERRQ(ierr); }

    isactive = PETSC_FALSE;
    ierr = PetscOptionsBool("-pswarm_hvar_finite_strain","Activate the tracking of finite strain","PSwarmSetFieldUpdateType",isactive,&isactive,0);CHKERRQ(ierr);
    if (isactive) { ierr = PSwarmSetFieldUpdateType(ps,PSWARM_FU_FINITESTRAIN);CHKERRQ(ierr); }

    isactive = PETSC_FALSE;
    ierr = PetscOptionsBool("-pswarm_pressure","Activate the tracking of pressure","PSwarmSetFieldUpdateType",isactive,&isactive,0);CHKERRQ(ierr);
    if (isactive) { ierr = PSwarmSetFieldUpdateType(ps,PSWARM_FU_Pressure);CHKERRQ(ierr); }

    isactive = PETSC_FALSE;
    ierr = PetscOptionsBool("-pswarm_view","View PSwarm info","PSwarmView",isactive,&isactive,0);CHKERRQ(ierr);
    if (isactive) { ierr = PSwarmViewInfo(ps);CHKERRQ(ierr); }

    ierr = PetscOptionsTail();CHKERRQ(ierr);
    ierr = PetscOptionsEnd();CHKERRQ(ierr);
  
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmCreateMultipleInstances"
PetscErrorCode PSwarmCreateMultipleInstances(MPI_Comm comm,PSwarm **pslist)
{
    PetscErrorCode ierr;
    PSwarm         *plist;
    PetscInt       k,nswarms;
    PetscInt       max = 20;
    PetscBool      found;
    char           *namelist[20];
    
    PetscFunctionBegin;
    ierr = PetscOptionsGetStringArray(NULL,NULL,"-pswarm_list",namelist,&max,&found);CHKERRQ(ierr);
    nswarms = max;
    
    PetscMalloc(sizeof(PSwarm)*(nswarms+1),&plist);
    PetscMemzero(plist,sizeof(PSwarm)*(nswarms+1));
    if (!found) {
      *pslist = plist;
      PetscFunctionReturn(0);
    }
  
    for (k=0; k<nswarms; k++) {
        char prefix[PETSC_MAX_PATH_LEN];
        
        ierr = PSwarmCreate(comm,&plist[k]);CHKERRQ(ierr);
        
        PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"%s_",namelist[k]);
        ierr = PSwarmSetOptionsPrefix(plist[k],prefix);CHKERRQ(ierr);
        
        //ierr = PSwarmSetPtatinCtx(plist[k],ctx);CHKERRQ(ierr);
        //if (X) ierr = PSwarmAttachStateVecVelocityPressure(plist[k],X);CHKERRQ(ierr);
        //if (T) ierr = PSwarmAttachStateVecTemperature(plist[k],T);CHKERRQ(ierr);
        //ierr = PSwarmSetFromOptions(plist[k]);CHKERRQ(ierr);
    }
    
    *pslist = plist;
    
    for (k=0; k<nswarms; k++) {
        PetscFree(namelist[k]);
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmCoordinatesSetSynchronization"
PetscErrorCode PSwarmCoordinatesSetSynchronization(PSwarm ps,PetscBool val)
{
  if (val) { ps->state = PSW_TS_INSYNC; }
  else { ps->state = PSW_TS_STALE; }
  PetscFunctionReturn(0);
}


/* ----------------- */
/* VTU headers writers */

/* defaults for int, float, double */
void PSWarmArray_VTUWriteBinaryAppendedHeader_int(FILE *vtk_fp,const char name[],int *offset,const int N)
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"%s\" format=\"appended\"  offset=\"%d\" />\n",name,*offset);
  *offset = *offset + sizeof(int) + N * sizeof(int);
}

void PSWarmArray_VTUWriteBinaryAppendedHeader_float(FILE *vtk_fp,const char name[],int *offset,const int N)
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\"  offset=\"%d\" />\n",name,*offset);
  *offset = *offset + sizeof(int) + N * sizeof(float);
}

void PSWarmArray_VTUWriteBinaryAppendedHeader_double(FILE *vtk_fp,const char name[],int *offset,const int N)
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"appended\"  offset=\"%d\" />\n",name,*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);
}

/* MPStd specific for phase, pid */
void MPntStd_VTUWriteBinaryAppendedHeader_phase(FILE *vtk_fp,int *offset,const int N)
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"phase\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(int);
}

void MPntStd_VTUWriteBinaryAppendedHeader_pid(FILE *vtk_fp,int *offset,const int N)
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int64\" Name=\"index\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(long int);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarm_VTUWriteBinaryAppendedHeaderAllFields"
PetscErrorCode PSwarm_VTUWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,DataBucket db,int npoints,int *byte_offset)
{
  BTruth found;
	
	PetscFunctionBegin;
  
  { /* MPStd */
    MPntStd_VTUWriteBinaryAppendedHeader_pid(vtk_fp,byte_offset,(const int)npoints);
    MPntStd_VTUWriteBinaryAppendedHeader_phase(vtk_fp,byte_offset,(const int)npoints);
  }
  
  { /* Pressure */
    DataBucketQueryDataFieldByName(db,"pressure",&found);
    if (found) {
      PSWarmArray_VTUWriteBinaryAppendedHeader_double(vtk_fp,"pressure",byte_offset,(const int)npoints);
    }
  }
  
	PetscFunctionReturn(0);
}

/* data writers */
/* defaults for int, float, double */
void PSwarmArray_VTUWriteBinaryAppendedData_int(FILE *vtk_fp,const int N,int data[])
{
  int length;
  size_t atomic_size;
  
  atomic_size = sizeof(int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite(&length,sizeof(int),1,vtk_fp);
  fwrite(data,atomic_size,N,vtk_fp);
}

void PSwarmArray_VTUWriteBinaryAppendedData_float(FILE *vtk_fp,const int N,double data[])
{
  int length;
  size_t atomic_size;
  
  atomic_size = sizeof(float);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite(&length,sizeof(int),1,vtk_fp);
  fwrite(data,atomic_size,N,vtk_fp);
}

void PSwarmArray_VTUWriteBinaryAppendedData_double(FILE *vtk_fp,const int N,double data[])
{
  int length;
  size_t atomic_size;
  
  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite(&length,sizeof(int),1,vtk_fp);
  fwrite(data,atomic_size,N,vtk_fp);
}

/* MPStd specific for phase, pid */
void MPntStd_VTUWriteBinaryAppendedData_phase(FILE *vtk_fp,const int N,const MPntStd points[])
{
  int p,length;
  size_t atomic_size;
  
  atomic_size = sizeof(int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite(&points[p].phase,atomic_size,1,vtk_fp);
  }
}

void MPntStd_VTUWriteBinaryAppendedData_pid(FILE *vtk_fp,const int N,const MPntStd points[])
{
  int p,length;
  size_t atomic_size;
  
  atomic_size = sizeof(long int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite(&points[p].pid,atomic_size,1,vtk_fp);
  }
}

#undef __FUNCT__
#define __FUNCT__ "PSwarm_VTUWriteBinaryAppendedDataAllFields"
PetscErrorCode PSwarm_VTUWriteBinaryAppendedDataAllFields(FILE *vtk_fp,DataBucket db)
{
	int npoints;
  BTruth found;
	DataField datafield;
  
	PetscFunctionBegin;
  
	DataBucketGetSizes(db,&npoints,NULL,NULL);
  
  { /* MPStd */
    MPntStd *tracer;

    DataBucketGetDataFieldByName(db,MPntStd_classname,&datafield);
    DataFieldGetEntries(datafield,(void**)&tracer);

    MPntStd_VTUWriteBinaryAppendedData_pid(vtk_fp,(const int)npoints,tracer);
    MPntStd_VTUWriteBinaryAppendedData_phase(vtk_fp,(const int)npoints,tracer);
  }
  
  { /* Pressure */
    DataBucketQueryDataFieldByName(db,"pressure",&found);
    if (found) {
      double *data;
      
      DataBucketGetDataFieldByName(db,"pressure",&datafield);
      DataFieldGetEntries(datafield,(void**)&data);

      PSwarmArray_VTUWriteBinaryAppendedData_double(vtk_fp,(const int)npoints,data);
    }
  }
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarm_PVTUWriteBinaryAppendedHeaderAllFields"
PetscErrorCode PSwarm_PVTUWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,DataBucket db)
{
  BTruth found;
  
  /* MPStd */
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Int64\" Name=\"index\" NumberOfComponents=\"1\"/>\n");
  
  fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Int32\" Name=\"phase\" NumberOfComponents=\"1\"/>\n");
  
  /* pressure */
  DataBucketQueryDataFieldByName(db,"pressure",&found);
  if (found) {
    fprintf(vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\"/>\n");
  }

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmView_VTUXML_binary_appended"
PetscErrorCode PSwarmView_VTUXML_binary_appended(PSwarm ps,const char name[])
{
	FILE *vtk_fp;
	PetscInt k;
	int npoints;
	PetscLogDouble t0,t1;
	DataField PField;
	int byte_offset,length;
  DataBucket db;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscTime(&t0);CHKERRQ(ierr);
	
	if ((vtk_fp = fopen (name,"w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
	}

	db = ps->db;
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	
	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
	
#ifdef WORDSIZE_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	fprintf( vtk_fp, "\t<UnstructuredGrid>\n" );
	
	DataBucketGetSizes(db,&npoints,NULL,NULL);
	fprintf( vtk_fp, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",npoints,npoints );
	
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<Cells>\n");
	
	byte_offset = 0;
	
	// connectivity //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
  byte_offset = byte_offset + sizeof(int) + npoints * sizeof(int);
	
	// offsets //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
  byte_offset = byte_offset + sizeof(int) + npoints * sizeof(int);
	
	// types //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
  byte_offset = byte_offset + sizeof(int) + npoints * sizeof(unsigned char);
	
	fprintf( vtk_fp, "\t\t\t</Cells>\n");
	
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<CellData>\n");
	fprintf( vtk_fp, "\t\t\t</CellData>\n");
	fprintf( vtk_fp, "\n");
	
	fprintf( vtk_fp, "\t\t\t<Points>\n");
	
	/* coordinates */
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
  byte_offset = byte_offset + sizeof(int) + npoints * 3 * sizeof(double);
	
	fprintf( vtk_fp, "\t\t\t</Points>\n");
	fprintf( vtk_fp, "\n");
	
	/* point data BEGIN */
	fprintf( vtk_fp, "\t\t\t<PointData>\n");
	/* auto generated shit for the header goes here */
	{
		ierr = PSwarm_VTUWriteBinaryAppendedHeaderAllFields(vtk_fp,db,npoints,&byte_offset);CHKERRQ(ierr);
	}
	fprintf( vtk_fp, "\t\t\t</PointData>\n");
	fprintf( vtk_fp, "\n");
	/* point data END */
	
	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</UnstructuredGrid>\n");
	
	/* WRITE APPENDED DATA HERE */
	fprintf( vtk_fp,"\t<AppendedData encoding=\"raw\">\n");
	fprintf( vtk_fp,"_");
	
	/* connectivity, offsets, types, coords */
	////////////////////////////////////////////////////////
	/* write connectivity */
	length = sizeof(int)*npoints;
	fwrite( &length,sizeof(int),1,vtk_fp);
	for (k=0; k<npoints; k++) {
		int idx = k;
		fwrite( &idx, sizeof(int),1, vtk_fp );
	}
	////////////////////////////////////////////////////////
	/* write offset */
	length = sizeof(int)*npoints;
	fwrite( &length,sizeof(int),1,vtk_fp);
	for (k=0; k<npoints; k++) {
		int idx = k+1;
		fwrite( &idx, sizeof(int),1, vtk_fp );
	}
	////////////////////////////////////////////////////////
	/* write types */
	length = sizeof(unsigned char)*npoints;
	fwrite( &length,sizeof(int),1,vtk_fp);
	for (k=0; k<npoints; k++) {
		unsigned char idx = 1; /* VTK_VERTEX */
		fwrite( &idx, sizeof(unsigned char),1, vtk_fp );
	}
	////////////////////////////////////////////////////////
	/* write coordinates */
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	
	length = sizeof(double)*npoints*3;
	fwrite( &length,sizeof(int),1,vtk_fp);
	for (k=0; k<npoints; k++) {
		MPntStd *marker;
		double  *coor;
		double  coords_k[] = {0.0, 0.0, 0.0};
		
		DataFieldAccessPoint(PField,k,(void**)&marker);
		MPntStdGetField_global_coord(marker,&coor);
		coords_k[0] = coor[0];
		coords_k[1] = coor[1];
		coords_k[2] = coor[2];
		
		fwrite( coords_k, sizeof(double), 3, vtk_fp );
	}
	DataFieldRestoreAccess(PField);
	
	/* auto generated shit for the marker data goes here */
	{
		ierr = PSwarm_VTUWriteBinaryAppendedDataAllFields(vtk_fp,db);CHKERRQ(ierr);
	}
	
	fprintf( vtk_fp,"\n\t</AppendedData>\n");
	
	fprintf( vtk_fp, "</VTKFile>\n");
	
	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}
	
	ierr = PetscTime(&t1);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/* deprecated */
#undef __FUNCT__
#define __FUNCT__ "_PSwarmViewMPntStd"
PetscErrorCode _PSwarmViewMPntStd(PSwarm ps)
{
	PetscErrorCode ierr;
  pTatinCtx ctx;
	static PetscBool beenhere=PETSC_FALSE;
	char pvdfilename[PETSC_MAX_PATH_LEN];
  const char *prefix;
  char vtkfilename[PETSC_MAX_PATH_LEN];
	char name[PETSC_MAX_PATH_LEN];
	
	PetscFunctionBegin;
  
	ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);
  ctx = ps->pctx;
	// PVD
  if (prefix) {
    PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_%spswarm.pvd",ctx->outputpath,prefix);
    PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1, "step%D_%spswarm.pvtu",ctx->step,prefix);
  } else {
    PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_pswarm.pvd",ctx->outputpath);
    PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1, "step%D_pswarm.pvtu",ctx->step);
  }
  if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename ); }
  ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time,vtkfilename,"");CHKERRQ(ierr);
  beenhere = PETSC_TRUE;
	
	// PVTS + VTS
  if (prefix) {
    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"step%D_%spswarm",ctx->step,prefix);
  } else {
    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"step%D_pswarm",ctx->step);
  }
	
	ierr = SwarmOutputParaView_MPntStd(ps->db,ctx->outputpath,name);CHKERRQ(ierr);
  
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_PSwarmViewStd"
PetscErrorCode _PSwarmViewStd(PSwarm ps,const char fileprefix[])
{
  PhysCompStokes stokes;
  DM             dmv,dmstokes;
  PetscErrorCode ierr;
  
  ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
  
  /*
   We could update the state here to ensure that particles which have left the
   sub-domain are not plotted. I'm not sure if this is really crucial or not.
   */
  /*
   if (ps->state == PSW_TS_STALE) {
   ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
   ps->state = PSW_TS_INSYNC;
   }
   */
  ierr = _PSwarmViewMPntStd(ps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmViewParaview_VTU"
PetscErrorCode PSwarmViewParaview_VTU(PSwarm ps,const char path[],const char stepprefix[],const char petscprefix[])
{
	char *vtkfilename,filename[PETSC_MAX_PATH_LEN],basename[PETSC_MAX_PATH_LEN];
  int n_points;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

  if (petscprefix) { PetscSNPrintf(basename,PETSC_MAX_PATH_LEN-1, "%s_%spswarm",stepprefix,petscprefix); }
  else {             PetscSNPrintf(basename,PETSC_MAX_PATH_LEN-1, "%s_pswarm",stepprefix); }
	
	ierr = pTatinGenerateParallelVTKName(basename,"vtu",&vtkfilename);CHKERRQ(ierr);
	if (path) { PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/%s",path,vtkfilename); }
	else {      PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"./%s",vtkfilename); }
  
  DataBucketGetSizes(ps->db,&n_points,NULL,NULL);
  if (n_points > 0) {
    ierr = PSwarmView_VTUXML_binary_appended(ps,filename);CHKERRQ(ierr);
  }
	free(vtkfilename);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "__PSwarmViewParaview_PVTU"
PetscErrorCode __PSwarmViewParaview_PVTU(DataBucket db,const char filename[],const char fileprefix[],PetscMPIInt nplist[])
{
	PetscMPIInt nproc;
	FILE *vtk_fp;
	PetscInt i;
	char *sourcename;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	if ((vtk_fp = fopen (filename,"w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",filename);
	}
	
	/* (VTK) generate pvts header */
	fprintf(vtk_fp,"<?xml version=\"1.0\"?>\n");
	
#ifdef WORDSIZE_BIGENDIAN
	fprintf(vtk_fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(vtk_fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	/* define size of the nodal mesh based on the cell DM */
	fprintf(vtk_fp,"  <PUnstructuredGrid GhostLevel=\"0\">\n" ); /* note overlap = 0 */
	
	/* DUMP THE CELL REFERENCES */
	fprintf(vtk_fp,"    <PCellData>\n");
	fprintf(vtk_fp,"    </PCellData>\n");
	
	///////////////
	fprintf(vtk_fp,"    <PPoints>\n");
	fprintf(vtk_fp,"      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	fprintf(vtk_fp,"    </PPoints>\n");
	///////////////
	
	///////////////
  fprintf(vtk_fp,"    <PPointData>\n");

  ierr = PSwarm_PVTUWriteBinaryAppendedHeaderAllFields(vtk_fp,db);CHKERRQ(ierr);
  
  fprintf(vtk_fp,"    </PPointData>\n");
	///////////////
	
	/* write out the parallel information */
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&nproc);CHKERRQ(ierr);
	for (i=0; i<nproc; i++) {
    int i32;
    
    if (nplist[i] != 0) {
      
      PetscMPIIntCast(i,&i32);
      if (asprintf(&sourcename,"%s-subdomain%1.5d.vtu",fileprefix,i32) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
      fprintf(vtk_fp,"    <Piece Source=\"%s\"/>\n",sourcename);
      free(sourcename);
    }
	}
	
	/* close the file */
	fprintf(vtk_fp,"  </PUnstructuredGrid>\n");
	fprintf(vtk_fp,"</VTKFile>\n");
	
	if (vtk_fp != NULL){
		fclose(vtk_fp);
		vtk_fp = NULL;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmViewParaview_PVTU"
PetscErrorCode PSwarmViewParaview_PVTU(DataBucket db,const char path[],const char stepprefix[],const char petscprefix[])
{
	char *vtufilename,fileprefix[PETSC_MAX_PATH_LEN],fileprefix2[PETSC_MAX_PATH_LEN];
	PetscMPIInt commsize,rank;
  int n_points,*nplist;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
  
  if (petscprefix) {
    PetscSNPrintf(fileprefix,PETSC_MAX_PATH_LEN-1,"%s_%spswarm",stepprefix,petscprefix);
  } else {
    PetscSNPrintf(fileprefix,PETSC_MAX_PATH_LEN-1,"%s_pswarm",stepprefix);
  }
  
  if (path) {
    PetscSNPrintf(fileprefix2,PETSC_MAX_PATH_LEN-1,"%s/%s",path,fileprefix);
  } else {
    PetscSNPrintf(fileprefix2,PETSC_MAX_PATH_LEN-1,"./%s",fileprefix);
  }
  
	ierr = pTatinGenerateVTKName(fileprefix2,"pvtu",&vtufilename);CHKERRQ(ierr);
	
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  DataBucketGetSizes(db,&n_points,NULL,NULL);
  if (rank == 0) {
    PetscMalloc(sizeof(int)*commsize,&nplist);
  }
  ierr = MPI_Gather(&n_points,1,MPI_INT,nplist,1,MPI_INT,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  if (rank == 0) {
    ierr = __PSwarmViewParaview_PVTU(db,vtufilename,fileprefix,nplist);CHKERRQ(ierr);
    PetscFree(nplist);
	}
  free(vtufilename);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmViewParaview_PVD"
PetscErrorCode PSwarmViewParaview_PVD(PSwarm ps,const char path[],const char stepprefix[],const char petscprefix[])
{
	PetscErrorCode ierr;
  pTatinCtx ctx;
	char pvdfilename[PETSC_MAX_PATH_LEN];
  const char *prefix;
  char vtkfilename[PETSC_MAX_PATH_LEN];
	
	PetscFunctionBegin;
  
	ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);
  ctx = ps->pctx;
  
  if (prefix) { PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_%spswarm.pvd",path,petscprefix); }
  else { PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_pswarm.pvd",path); }
	if (!ps->pvdopen) {
		PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename );
		ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
		
		ps->pvdopen = PETSC_TRUE;
	}
  if (prefix) { PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1, "%s_%spswarm.pvtu",stepprefix,petscprefix); }
  else {        PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1, "%s_pswarm.pvtu",stepprefix); }
  
  ierr = ParaviewPVDAppend(pvdfilename,ctx->time,vtkfilename,stepprefix);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
 STEPPREFIX = "step00000"
 
 filename for .pvd
 OUTPUTPATH/timeseries_PETSCPREFIX_pswarm.pvtu
 
 filenames for pvd file references
 STEPPREFIX_PETSCPREFIX_pswarm.pvtu
 
 filename for .pvtu
 OUTPUTPATH/STEPPREFIX_PETSCPREFIX_pswarm.pvtu

 filenames for .pvtu file references
 STEPPREFIX_PETSCPREFIX_pswarm-subdomainRANK.vtu

 
*/
#undef __FUNCT__
#define __FUNCT__ "PSwarmView_PerRank"
PetscErrorCode PSwarmView_PerRank(PSwarm ps)
{
  PhysCompStokes stokes;
  DM             dmv,dmstokes;
  char           stepprefix[PETSC_MAX_PATH_LEN],pvoutputdir[PETSC_MAX_PATH_LEN];
  const char     *petscprefix;
  PetscBool      found;
  PetscErrorCode ierr;
  
  ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
  
  /*
   We could update the state here to ensure that particles which have left the
   sub-domain are not plotted. I'm not sure if this is really crucial or not.
   */
  /*
   if (ps->state == PSW_TS_STALE) {
   ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
   ps->state = PSW_TS_INSYNC;
   }
   */
  
	ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&petscprefix);CHKERRQ(ierr);
  PetscSNPrintf(stepprefix,PETSC_MAX_PATH_LEN-1,"step%D",ps->pctx->step);
  
  ierr = PetscSNPrintf(pvoutputdir,PETSC_MAX_PATH_LEN-1,"%s/step%D",ps->pctx->outputpath,ps->pctx->step);CHKERRQ(ierr);
  PetscTestDirectory(pvoutputdir,'w',&found);
  if (!found) { ierr = pTatinCreateDirectory(pvoutputdir);CHKERRQ(ierr); }

  ierr = PSwarmViewParaview_PVD(ps,ps->pctx->outputpath,stepprefix,petscprefix);CHKERRQ(ierr);
  ierr = PSwarmViewParaview_VTU(ps,pvoutputdir,stepprefix,petscprefix);CHKERRQ(ierr);
  ierr = PSwarmViewParaview_PVTU(ps->db,pvoutputdir,stepprefix,petscprefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetRegionIndex"
PetscErrorCode PSwarmSetRegionIndex(PSwarm ps,PetscInt ridx)
{
  MPntStd *tracer;
  DataField datafield;
  int p,np;
  
  DataBucketGetDataFieldByName(ps->db,MPntStd_classname,&datafield);
  DataFieldGetEntries(datafield,(void**)&tracer);
  DataBucketGetSizes(ps->db,&np,NULL,NULL);
  for (p=0; p<np; p++) {
    tracer[p].phase = ridx;
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pSwarmParaViewMeshDeformationBaseVTS"
PetscErrorCode pSwarmParaViewMeshDeformationBaseVTS(PetscReal xmin[],PetscReal dx[],PetscInt nx[],const char name[])
{
	PetscInt i,j,k;
	FILE*	vtk_fp = NULL;
  float xp,yp,zp;
	
	PetscFunctionBegin;
	if ((vtk_fp = fopen (name,"w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
	}
	
	/* VTS HEADER - OPEN */
#ifdef WORDSIZE_BIGENDIAN
	fprintf(vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	fprintf(vtk_fp,"  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0,nx[0]-1, 0,nx[1]-1,0,nx[2]-1);
	fprintf(vtk_fp,"    <Piece Extent=\"%d %d %d %d %d %d\">\n", 0,nx[0]-1, 0,nx[1]-1,0,nx[2]-1);
	
	/* VTS COORD DATA */
	fprintf(vtk_fp,"    <Points>\n");
	fprintf(vtk_fp,"      <DataArray Name=\"coords\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for (k=0; k<nx[2]; k++) {
		for (j=0; j<nx[1]; j++) {
			for (i=0; i<nx[0]; i++) {
        xp = xmin[0] + dx[0] * i;
        yp = xmin[1] + dx[1] * j;
        zp = xmin[2] + dx[2] * k;
				fprintf(vtk_fp,"      %1.6e %1.6e %1.6e\n",xp,yp,zp);
			}
		}
	}
	fprintf(vtk_fp,"      </DataArray>\n");
	fprintf(vtk_fp,"    </Points>\n");
	
	/* VTS CELL DATA */
	fprintf(vtk_fp,"    <CellData>\n");
	
	fprintf(vtk_fp,"      <DataArray Name=\"index\" type=\"Int32\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	fprintf(vtk_fp,"      ");
	for (k=0; k<nx[2]; k++) {
		for (j=0; j<nx[1]; j++) {
			for (i=0; i<nx[0]; i++) {
				fprintf(vtk_fp,"%d ", 0 );
			}
		}
	}
	fprintf(vtk_fp,"\n");
	fprintf(vtk_fp,"      </DataArray>\n");
	
	fprintf(vtk_fp,"    </CellData>\n");
	
	/* VTS NODAL DATA */
	fprintf(vtk_fp,"    <PointData>\n");
	fprintf(vtk_fp,"    </PointData>\n");
	
	/* VTS HEADER - CLOSE */
	fprintf(vtk_fp,"    </Piece>\n");
	fprintf(vtk_fp,"  </StructuredGrid>\n");
	fprintf(vtk_fp,"</VTKFile>\n");
	fclose(vtk_fp);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmViewSingletonParaview_PVD"
PetscErrorCode PSwarmViewSingletonParaview_PVD(PSwarm ps,const char path[],const char stepprefix[],const char petscprefix[])
{
	PetscErrorCode ierr;
  pTatinCtx ctx;
	char pvdfilename[PETSC_MAX_PATH_LEN];
  const char *prefix;
  char vtkfilename[PETSC_MAX_PATH_LEN];
	
	PetscFunctionBegin;
  
	ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);
  ctx = ps->pctx;
  
  if (prefix) { PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_%spswarm.pvd",path,petscprefix); }
  else { PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_pswarm.pvd",path); }
	if (!ps->pvdopen) {
		PetscPrintf(PetscObjectComm((PetscObject)ps),"  writing pvdfilename %s \n", pvdfilename );
		ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
		
		ps->pvdopen = PETSC_TRUE;
	}
  if (prefix) { PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1, "%s_%spswarm.vtu",stepprefix,petscprefix); }
  else {        PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1, "%s_pswarm.vtu",stepprefix); }
  
  ierr = ParaviewPVDAppend(pvdfilename,ctx->time,vtkfilename,stepprefix);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSingleton_VTUWriteBinaryAppendedDataAllFields"
PetscErrorCode PSwarmSingleton_VTUWriteBinaryAppendedDataAllFields(FILE *vtk_fp,DataBucket db,MPI_Comm comm)
{
	int i,npoints,npoints_g = 0;
  BTruth found;
	DataField datafield;
  PetscMPIInt rank;
  size_t atomic_size;
  PetscErrorCode ierr;
  
	PetscFunctionBegin;
  
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	DataBucketGetSizes(db,&npoints,NULL,NULL);
  ierr = MPI_Reduce(&npoints,&npoints_g,1,MPI_INT,MPI_SUM,0,comm);CHKERRQ(ierr);
  
  { /* MPStd - pid */
    MPntStd *tracer;
    long int *buffer;
    
    DataBucketGetDataFieldByName(db,MPntStd_classname,&datafield);
    DataFieldGetEntries(datafield,(void**)&tracer);
    
    PetscMalloc(sizeof(long int)*npoints,&buffer);
    for (i=0; i<npoints; i++) {
      buffer[i] = tracer[i].pid;
    }
    
    atomic_size = sizeof(long int);
    if (rank == 0) {
      int length;

      length = (int)( atomic_size * ((size_t)npoints_g) );
      fwrite(&length,sizeof(int),1,vtk_fp);
    }

    ierr = MPIWrite_Blocking(vtk_fp,(void*)buffer,npoints,atomic_size,0,PETSC_TRUE,comm);CHKERRQ(ierr);

    PetscFree(buffer);
  }

  { /* MPStd - region */
    MPntStd *tracer;
    int *buffer;
    
    DataBucketGetDataFieldByName(db,MPntStd_classname,&datafield);
    DataFieldGetEntries(datafield,(void**)&tracer);
    
    PetscMalloc(sizeof(int)*npoints,&buffer);
    for (i=0; i<npoints; i++) {
      buffer[i] = tracer[i].phase;
    }
    
    atomic_size = sizeof(int);
    if (rank == 0) {
      int length;
      
      length = (int)( atomic_size * ((size_t)npoints_g) );
      fwrite(&length,sizeof(int),1,vtk_fp);
    }
    
    ierr = MPIWrite_Blocking(vtk_fp,(void*)buffer,npoints,atomic_size,0,PETSC_TRUE,comm);CHKERRQ(ierr);
    
    PetscFree(buffer);
  }
  
  { /* Pressure */
    DataBucketQueryDataFieldByName(db,"pressure",&found);
    if (found) {
      double *data;
      
      DataBucketGetDataFieldByName(db,"pressure",&datafield);
      DataFieldGetEntries(datafield,(void**)&data);
      
      atomic_size = sizeof(double);
      if (rank == 0) {
        int length;
        
        length = (int)( atomic_size * ((size_t)npoints_g) );
        fwrite(&length,sizeof(int),1,vtk_fp);
      }

      ierr = MPIWrite_Blocking(vtk_fp,(void*)data,npoints,atomic_size,0,PETSC_TRUE,comm);CHKERRQ(ierr);

    }
  }
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmViewSingleton_VTUXML_binary_appended"
PetscErrorCode PSwarmViewSingleton_VTUXML_binary_appended(PSwarm ps,const char name[])
{
	FILE *vtk_fp;
	PetscInt k;
	int npoints,npoints_g;
	PetscLogDouble t0,t1;
	DataField PField;
	int byte_offset,length;
  DataBucket db;
  PetscMPIInt rank;
  MPI_Comm comm;
  double *buffer;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscTime(&t0);CHKERRQ(ierr);
  
  comm = PetscObjectComm((PetscObject)ps);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  vtk_fp = NULL;
  if (rank == 0) {
    if ((vtk_fp = fopen (name,"w")) == NULL)  {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
    }
  }
  
	db = ps->db;
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	
  if (rank == 0) {
    fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
    
  #ifdef WORDSIZE_BIGENDIAN
    fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  #else
    fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  #endif
    
    fprintf( vtk_fp, "\t<UnstructuredGrid>\n" );
  }
	
	DataBucketGetSizes(db,&npoints,NULL,NULL);
  ierr = MPI_Reduce(&npoints,&npoints_g,1,MPI_INT,MPI_SUM,0,comm);CHKERRQ(ierr);
  
  if (rank == 0) {
    fprintf( vtk_fp, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",npoints_g,npoints_g );
    
    fprintf( vtk_fp, "\n");
    fprintf( vtk_fp, "\t\t\t<Cells>\n");
  }
	
	byte_offset = 0;
	
  if (rank == 0) {
    // connectivity //
    fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
    byte_offset = byte_offset + sizeof(int) + npoints_g * sizeof(int);
    
    // offsets //
    fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
    byte_offset = byte_offset + sizeof(int) + npoints_g * sizeof(int);
    
    // types //
    fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
    byte_offset = byte_offset + sizeof(int) + npoints_g * sizeof(unsigned char);
    
    fprintf( vtk_fp, "\t\t\t</Cells>\n");
    
    fprintf( vtk_fp, "\n");
    fprintf( vtk_fp, "\t\t\t<CellData>\n");
    fprintf( vtk_fp, "\t\t\t</CellData>\n");
    fprintf( vtk_fp, "\n");
    
    fprintf( vtk_fp, "\t\t\t<Points>\n");
    
    /* coordinates */
    fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
    byte_offset = byte_offset + sizeof(int) + npoints_g * 3 * sizeof(double);
    
    fprintf( vtk_fp, "\t\t\t</Points>\n");
    fprintf( vtk_fp, "\n");
  }

	if (rank == 0) {
    /* point data BEGIN */
    fprintf( vtk_fp, "\t\t\t<PointData>\n");
    /* auto generated shit for the header goes here */
    {
      ierr = PSwarm_VTUWriteBinaryAppendedHeaderAllFields(vtk_fp,db,npoints_g,&byte_offset);CHKERRQ(ierr);
    }
    fprintf( vtk_fp, "\t\t\t</PointData>\n");
    fprintf( vtk_fp, "\n");
    /* point data END */
  }
	
	if (rank == 0) {
    fprintf( vtk_fp, "\t\t</Piece>\n");
    fprintf( vtk_fp, "\t</UnstructuredGrid>\n");
    
    /* WRITE APPENDED DATA HERE */
    fprintf( vtk_fp,"\t<AppendedData encoding=\"raw\">\n");
    fprintf( vtk_fp,"_");
  }

	if (rank == 0) {
    /* connectivity, offsets, types, coords */
    ////////////////////////////////////////////////////////
    /* write connectivity */
    length = sizeof(int)*npoints_g;
    fwrite( &length,sizeof(int),1,vtk_fp);
    for (k=0; k<npoints_g; k++) {
      int idx = k;
      fwrite( &idx, sizeof(int),1, vtk_fp );
    }
    
    ////////////////////////////////////////////////////////
    /* write offset */
    length = sizeof(int)*npoints_g;
    fwrite( &length,sizeof(int),1,vtk_fp);
    for (k=0; k<npoints_g; k++) {
      int idx = k+1;
      fwrite( &idx, sizeof(int),1, vtk_fp );
    }

    ////////////////////////////////////////////////////////
    /* write types */
    length = sizeof(unsigned char)*npoints_g;
    fwrite( &length,sizeof(int),1,vtk_fp);
    for (k=0; k<npoints_g; k++) {
      unsigned char idx = 1; // VTK_VERTEX //
      fwrite( &idx, sizeof(unsigned char),1, vtk_fp );
    }
  }
  
	////////////////////////////////////////////////////////
	/* write coordinates */
	PetscMalloc(sizeof(double)*npoints*3,&buffer);
  
  DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
  for (k=0; k<npoints; k++) {
		MPntStd *marker;
		double  *coor;
		
		DataFieldAccessPoint(PField,k,(void**)&marker);
		MPntStdGetField_global_coord(marker,&coor);
		buffer[3*k+0] = coor[0];
		buffer[3*k+1] = coor[1];
		buffer[3*k+2] = coor[2];
	}
	DataFieldRestoreAccess(PField);

  if (rank == 0) {
    length = sizeof(double)*npoints_g*3;
    fwrite( &length,sizeof(int),1,vtk_fp);
  }
  ierr = MPIWrite_Blocking(vtk_fp,(void*)buffer,npoints,3*sizeof(double),0,PETSC_TRUE,comm);CHKERRQ(ierr);
  PetscFree(buffer);
  
	/* auto generated shit for the marker data goes here */
	{
    ierr = PSwarmSingleton_VTUWriteBinaryAppendedDataAllFields(vtk_fp,db,comm);CHKERRQ(ierr);
  }
	
  if (rank == 0) {
    fprintf( vtk_fp,"\n\t</AppendedData>\n");
    fprintf( vtk_fp, "</VTKFile>\n");
    
    if( vtk_fp!= NULL ) {
      fclose( vtk_fp );
      vtk_fp = NULL;
    }
  }
	
	ierr = PetscTime(&t1);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmViewSingletonParaview_VTU"
PetscErrorCode PSwarmViewSingletonParaview_VTU(PSwarm ps,const char path[],const char stepprefix[],const char petscprefix[])
{
	char vtkfilename[PETSC_MAX_PATH_LEN],filename[PETSC_MAX_PATH_LEN],basename[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
  
  if (petscprefix) { PetscSNPrintf(basename,PETSC_MAX_PATH_LEN-1, "%s_%spswarm",stepprefix,petscprefix); }
  else {             PetscSNPrintf(basename,PETSC_MAX_PATH_LEN-1, "%s_pswarm",stepprefix); }
	
  PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1, "%s.vtu",basename);
	if (path) { PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/%s",path,vtkfilename); }
	else {      PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"./%s",vtkfilename); }
  
  ierr = PSwarmViewSingleton_VTUXML_binary_appended(ps,filename);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmView_Singleton"
PetscErrorCode PSwarmView_Singleton(PSwarm ps)
{
  PhysCompStokes stokes;
  DM             dmv,dmstokes;
  char           stepprefix[PETSC_MAX_PATH_LEN],pvoutputdir[PETSC_MAX_PATH_LEN];
  const char     *petscprefix;
  PetscBool      found;
  PetscErrorCode ierr;
  
  ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
  
	ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&petscprefix);CHKERRQ(ierr);
  PetscSNPrintf(stepprefix,PETSC_MAX_PATH_LEN-1,"step%D",ps->pctx->step);
  
  ierr = PetscSNPrintf(pvoutputdir,PETSC_MAX_PATH_LEN-1,"%s/step%D",ps->pctx->outputpath,ps->pctx->step);CHKERRQ(ierr);
  PetscTestDirectory(pvoutputdir,'w',&found);
  if (!found) { ierr = pTatinCreateDirectory(pvoutputdir);CHKERRQ(ierr); }

  ierr = PSwarmViewSingletonParaview_PVD(ps,ps->pctx->outputpath,stepprefix,petscprefix);CHKERRQ(ierr);
  
  ierr = PSwarmViewSingletonParaview_VTU(ps,pvoutputdir,stepprefix,petscprefix);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmView"
PetscErrorCode PSwarmView(PSwarm ps,PSwarmViewType type)
{
  PetscErrorCode ierr;
  
  if (type == PSW_VT_PERRANK) {
    ierr = PSwarmView_PerRank(ps);CHKERRQ(ierr);
  } else if (type == PSW_VT_SINGLETON) {
    ierr = PSwarmView_Singleton(ps);CHKERRQ(ierr);
  } else {
    SETERRQ(PetscObjectComm((PetscObject)ps),PETSC_ERR_SUP,"Unknown PSwarmViewType detected");
  }
  
  PetscFunctionReturn(0);
}




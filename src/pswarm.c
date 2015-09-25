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
#include <ptatin3d.c>
#include <ptatin3d_stokes.h>
#include <ptatin3d_energy.h>
#include <pswarm.h>
#include <pswarm_impl.h>


PetscClassId PSWARM_CLASSID;

const char PSWARM_COMPOSED_STATE_VELPRES[] = "PSWarmStateVector_VP";
const char PSWARM_COMPOSED_STATE_TEMP[]    = "PSWarmStateVector_T";

PetscErrorCode SwarmDMDA3dDataExchangerCreate(DM da,DataEx *_de);

PetscErrorCode PSwarmView(PSwarm ps,PetscViewer viewer);
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
    
    PetscHeaderCreate(p, _p_PSwarm, struct _PSwarmOps, PSWARM_CLASSID, "PSwarm", "Particle Swarm Manager", "PSwarm", comm, PSwarmDestroy, PSwarmView);
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
        ierr = PetscObjectCompose((PetscObject)ps,PSWARM_COMPOSED_STATE_TEMP,(PetscObject)x);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_PSwarmViewMPntStd"
PetscErrorCode _PSwarmViewMPntStd(PSwarm ps)
{
	PetscErrorCode ierr;
    pTatinCtx ctx;
	char *name;
	static int beenhere=0;
	static char *pvdfilename;
    const char *prefix;
	
	PetscFunctionBegin;

	ierr = PetscObjectGetOptionsPrefix((PetscObject)ps,&prefix);CHKERRQ(ierr);
    ctx = ps->pctx;
	// PVD
	if (beenhere == 0) {
        if (prefix) {
            asprintf(&pvdfilename,"%s/timeseries_%spswarm_std.pvd",ctx->outputpath,prefix);
        } else {
            asprintf(&pvdfilename,"%s/timeseries_pswarm_std.pvd",ctx->outputpath);
        }
		PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename );
		ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
		
		beenhere = 1;
	}
	{
		char *vtkfilename;

		if (prefix) {
            asprintf(&vtkfilename, "%d_%spswarm_std.pvtu",ctx->step,prefix);
		} else {
            asprintf(&vtkfilename, "%d_pswarm_std.pvtu",ctx->step);
        }
        
		ierr = ParaviewPVDAppend(pvdfilename,ctx->time,vtkfilename,"");CHKERRQ(ierr);
		free(vtkfilename);
	}
	
	// PVTS + VTS
    if (prefix) {
        asprintf(&name,"%d_%spswarm_std",ctx->step,prefix);
    } else {
        asprintf(&name,"%d_pswarm_std",ctx->step);
    }
	
	ierr = SwarmOutputParaView_MPntStd(ps->db,ctx->outputpath,name);CHKERRQ(ierr);
  
	free(name);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmView"
PetscErrorCode PSwarmView(PSwarm ps,PetscViewer viewer)
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
#define __FUNCT__ "PSwarmViewInfo"
PetscErrorCode PSwarmViewInfo(PSwarm ps)
{
  PetscErrorCode ierr;
  int n_points;
  
  DataBucketGetSizes(ps->db,&n_points,0,0);

  PetscPrintf(PETSC_COMM_WORLD,"PSwarm --------------------\n");
  PetscPrintf(PETSC_COMM_WORLD,"  npoints %D\n",n_points);
  PetscPrintf(PETSC_COMM_WORLD,"  Transport mode: \n");
  if (ps->transport_mode == PSWARM_TM_EULERIAN) {
    PetscPrintf(PETSC_COMM_WORLD,"    Eulerian\n");
  }
  if (ps->transport_mode == PSWARM_TM_LAGRANGIAN) {
    PetscPrintf(PETSC_COMM_WORLD,"    Lagrangian\n");
  }
  PetscPrintf(PETSC_COMM_WORLD,"  Update methods: \n");
  if (ps->ops->field_update_finitestrain) {
    PetscPrintf(PETSC_COMM_WORLD,"    Finite strain\n");
  }
  if (ps->ops->field_update_ptt) {
    PetscPrintf(PETSC_COMM_WORLD,"    Pressure-Temperature-Time\n");
  }
  if (ps->ops->field_update_pressure) {
    PetscPrintf(PETSC_COMM_WORLD,"    Pressure\n");
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
            ps->ops->advect = _PSwarmFieldUpdate_AdvectEulerian;
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
  
  if (ps->state == PSW_TS_STALE) {
    /* update local coordinates and perform communication */
		ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
    ps->state = PSW_TS_INSYNC;
  }
  
  DataBucketGetDataFieldByName(ps->db,MPntStd_classname,&datafield_tracers);
	DataFieldGetEntries(datafield_tracers,(void**)&tracer);

  DataBucketGetDataFieldByName(ps->db,"pressure",&datafield);
	DataFieldGetEntries(datafield,(void**)&tracer_pressure);
  
	DataBucketGetSizes(ps->db,&n_tracers,0,0);
  for (p=0; p<n_tracers; p++) {
    tracer_pressure[p] = tracer[p].coor[1];
  }
  
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
    if (energy) {
        dmT = energy->daT;
        ierr = PetscObjectQuery((PetscObject)ps,PSWARM_COMPOSED_STATE_TEMP,(PetscObject*)&temperature);CHKERRQ(ierr);
    }

    ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
    
    if (ps->state == PSW_TS_STALE) {
		ierr = MaterialPointStd_UpdateCoordinates(ps->db,dmv,ps->de);CHKERRQ(ierr);
        ps->state = PSW_TS_INSYNC;
    }
    
    if (ps->ops->field_update_finitestrain) {
        ierr = ps->ops->field_update_finitestrain(ps,dmv,velocity);CHKERRQ(ierr);
    }
    
    if (ps->ops->field_update_ptt) {
        ierr = ps->ops->field_update_ptt(ps,dmp,dmT,pressure,temperature,ps->pctx->time);CHKERRQ(ierr);
    }

    if (ps->ops->field_update_pressure) {
      ierr = ps->ops->field_update_pressure(ps,dmv,dmp,pressure);CHKERRQ(ierr);
    }

    /* position must always be the last state variable to be updated */
    if (ps->ops->advect) {
        ierr = ps->ops->advect(ps,dmv,velocity);CHKERRQ(ierr);
        /* flag as being stale => local coordinates need updating */
        ps->state = PSW_TS_STALE;
    }
    ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetUpCoords_FillDM"
PetscErrorCode PSwarmSetUpCoords_FillDM(PSwarm ps)
{
    PetscErrorCode ierr;
    PhysCompStokes stokes;
    DM             dmv;
    PetscInt       lmx,lmy,lmz;
    PetscInt       Nxp[] = {1,1,1}; /* change with -lattice_layout_N{x,y,z} */
    
    PetscFunctionBegin;
    
    ierr = pTatinGetStokesContext(ps->pctx,&stokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMs(stokes,&dmv,NULL);CHKERRQ(ierr);
    
    ierr = DMDAGetLocalSizeElementQ2(dmv,&lmx,&lmy,&lmz);CHKERRQ(ierr);
    
    ierr = SwarmMPntStd_CoordAssignment_LatticeLayout3d(dmv,Nxp,0.0,ps->db);CHKERRQ(ierr);
    
    ps->state = PSW_TS_INSYNC;
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmSetUpCoords"
PetscErrorCode PSwarmSetUpCoords(PSwarm ps)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    
    ierr = PSwarmSetUpCoords_FillDM(ps);CHKERRQ(ierr);
    //ierr = PSwarmSetUpCoords_FillDMWithGeometryObject(ps);CHKERRQ(ierr);
    //ierr = PSwarmSetUpCoords_FromList(ps);CHKERRQ(ierr);
    
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
    ierr = PetscOptionsHead("PSwarm options");CHKERRQ(ierr);
    
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
    ierr = PetscOptionsBool("-pswarm_view","View PSwaem info","PSwarmView",isactive,&isactive,0);CHKERRQ(ierr);
    if (isactive) { ierr = PSwarmViewInfo(ps);CHKERRQ(ierr); }

    ierr = PetscOptionsTail();CHKERRQ(ierr);
    ierr = PetscOptionsEnd();CHKERRQ(ierr);
  
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarmCreateMultipleInstances"
PetscErrorCode PSwarmCreateMultipleInstances(pTatinCtx ctx,Vec X,Vec T,PSwarm **pslist)
{
    PetscErrorCode ierr;
    PSwarm         *plist;
    PetscInt       k,nswarms;
    MPI_Comm       comm;
    PetscInt       max = 20;
    PetscBool      found;
    char           *namelist[20];
    
    PetscFunctionBegin;
    ierr = PetscOptionsGetStringArray(NULL,"-pswarm_list",namelist,&max,&found);CHKERRQ(ierr);
    if (!found) {
        *pslist = NULL;
        PetscFunctionReturn(0);
    }
    
    nswarms = max;
    
    PetscMalloc(sizeof(PSwarm)*(nswarms+1),&plist);
    PetscMemzero(plist,sizeof(PSwarm)*(nswarms+1));

    ierr = PetscObjectGetComm((PetscObject)ctx->pack,&comm);CHKERRQ(ierr);

    for (k=0; k<nswarms; k++) {
        char prefix[PETSC_MAX_PATH_LEN];
        
        ierr = PSwarmCreate(comm,&plist[k]);CHKERRQ(ierr);
        
        PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN,"%s_",namelist[k]);
        ierr = PSwarmSetOptionsPrefix(plist[k],prefix);CHKERRQ(ierr);
        
        ierr = PSwarmSetPtatinCtx(plist[k],ctx);CHKERRQ(ierr);
        ierr = PSwarmAttachStateVecVelocityPressure(plist[k],X);CHKERRQ(ierr);
        ierr = PSwarmAttachStateVecTemperature(plist[k],T);CHKERRQ(ierr);
        ierr = PSwarmSetFromOptions(plist[k]);CHKERRQ(ierr);
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
/* headers writers */
void PSWarmArray_VTKWriteBinaryAppendedHeader_int(FILE *vtk_fp,const char name[],int *offset,const int N)
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"%s\" format=\"appended\"  offset=\"%d\" />\n",name,*offset);
  *offset = *offset + sizeof(int) + N * sizeof(int);
}

void PSWarmArray_VTKWriteBinaryAppendedHeader_double(FILE *vtk_fp,const char name[],int *offset,const int N)
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"appended\"  offset=\"%d\" />\n",name,*offset);
  *offset = *offset + sizeof(int) + N * sizeof(double);
}
void MPntStd_VTKWriteBinaryAppendedHeader_Phase(FILE *vtk_fp,int *offset,const int N)
{
  fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"phase\" format=\"appended\"  offset=\"%d\" />\n",*offset);
  *offset = *offset + sizeof(int) + N * sizeof(int);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarm_VTKWriteBinaryAppendedHeaderAllFields"
PetscErrorCode PSwarm_VTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,DataBucket db,int *byte_offset)
{
	int npoints;
  BTruth found;
	
	PetscFunctionBegin;
  
	DataBucketGetSizes(db,&npoints,NULL,NULL);
  
  { /* MPStd */
    MPntStd_VTKWriteBinaryAppendedHeader_Phase(vtk_fp,byte_offset,(const int)npoints);
  }
  
  { /* Pressure */
    DataBucketQueryDataFieldByName(db,"pressure",&found);
    if (found) {
      PSWarmArray_VTKWriteBinaryAppendedHeader_double(vtk_fp,"pressure",byte_offset,(const int)npoints);
    }
  }
  
	PetscFunctionReturn(0);
}

/* data writers */
void MPntStd_VTKWriteBinaryAppendedData_Phase(FILE *vtk_fp,const int N,const MPntStd points[])
{
  int p,length;
  size_t atomic_size;
  
  atomic_size = sizeof(int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  for(p=0;p<N;p++) {
    fwrite( &points[p].phase,atomic_size,1,vtk_fp);
  }
}

void PSwarmArray_VTKWriteBinaryAppendedData_int(FILE *vtk_fp,const int N,int data[])
{
  int length;
  size_t atomic_size;
  
  atomic_size = sizeof(int);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  fwrite( data,atomic_size,N,vtk_fp);
}

void PSwarmArray_VTKWriteBinaryAppendedData_double(FILE *vtk_fp,const int N,double data[])
{
  int length;
  size_t atomic_size;
  
  atomic_size = sizeof(double);
  length = (int)( atomic_size * ((size_t)N) );
  fwrite( &length,sizeof(int),1,vtk_fp);
  fwrite( data,atomic_size,N,vtk_fp);
}

#undef __FUNCT__
#define __FUNCT__ "PSwarm_VTKWriteBinaryAppendedDataAllFields"
PetscErrorCode PSwarm_VTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,DataBucket db)
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

    MPntStd_VTKWriteBinaryAppendedData_Phase(vtk_fp,(const int)npoints,tracer);
  }
  
  { /* Pressure */
    DataBucketQueryDataFieldByName(db,"pressure",&found);
    if (found) {
      double *data;
      
      DataBucketGetDataFieldByName(db,"pressure",&datafield);
      DataFieldGetEntries(datafield,(void**)&data);

      PSwarmArray_VTKWriteBinaryAppendedData_double(vtk_fp,(const int)npoints,data);
    }
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
		ierr = PSwarm_VTKWriteBinaryAppendedHeaderAllFields(vtk_fp,db,&byte_offset);CHKERRQ(ierr);
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
		ierr = PSwarm_VTKWriteBinaryAppendedDataAllFields(vtk_fp,db);CHKERRQ(ierr);
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

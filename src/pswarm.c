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
#include <ptatin3d_stokes.h>
#include <ptatin3d_energy.h>
#include <pswarm.h>
#include <pswarm_impl.h>


const char PSWARM_COMPOSED_STATE_VELPRES[] = "PSWarmStateVector_VP";
const char PSWARM_COMPOSED_STATE_TEMP[]    = "PSWarmStateVector_T";

PetscErrorCode SwarmDMDA3dDataExchangerCreate(DM da,DataEx *_de);

PetscErrorCode PSwarmView(PSwarm ps,PetscViewer viewer);
PetscErrorCode PSwarmDestroy(PSwarm *ps);
PetscErrorCode _PSwarmFieldUpdate_AdvectEulerian(PSwarm ps,DM dmv,Vec v);


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
#define __FUNCT__ "PSwarmView"
PetscErrorCode PSwarmView(PSwarm ps,PetscViewer viewer)
{
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
    switch (type) {
        case PSWARM_TM_EULERIAN:
            ps->ops->advect = _PSwarmFieldUpdate_AdvectEulerian;
            break;
        case PSWARM_TM_LAGRANGIAN:
            ps->ops->advect = NULL;
            break;
            
        default:
            ps->ops->advect = NULL;
            break;
    }
    
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
    ierr = PetscOptionsHead("PSwarm options");CHKERRQ(ierr);
    
    isactive = PETSC_FALSE;
    ierr = PetscOptionsBool("-pswarm_transport_mode_eulerian","Transport mode set to Eulerian","PSwarmSetTransportModeType",isactive,&isactive,0);CHKERRQ(ierr);
    if (isactive) { ierr = PSwarmSetTransportModeType(ps,PSWARM_TM_EULERIAN);CHKERRQ(ierr); }

    isactive = PETSC_FALSE;
    ierr = PetscOptionsBool("-pswarm_hvar_finite_strain","Activate the tracking of finite strain","PSwarmSetFieldUpdateType",isactive,&isactive,0);CHKERRQ(ierr);
    if (isactive) { ierr = PSwarmSetFieldUpdateType(ps,PSWARM_FU_FINITESTRAIN);CHKERRQ(ierr); }

    
    ierr = PetscOptionsTail();CHKERRQ(ierr);
    
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

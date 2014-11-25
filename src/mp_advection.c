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
 **    Filename:      mp_advection.c
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


#define PTAT3D_PROFILE_SwarmUpdatePosition

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "petscdm.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "data_bucket.h"
#include "data_exchanger.h"
#include "MPntStd_def.h"
#include "ptatin3d_stokes.h"
#include "element_utils_q2.h"
#include "dmda_element_q2p1.h"
#include "material_point_point_location.h"
#include "mp_advection.h"

PetscLogEvent PTATIN_MaterialPointAdvGlobalCoordUpdate;
PetscLogEvent PTATIN_MaterialPointAdvLocalCoordUpdate;
PetscLogEvent PTATIN_MaterialPointAdvCommunication;
PetscLogEvent PTATIN_MaterialPointAdvRemoval;


#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_AdvectEuler"
PetscErrorCode MaterialPointStd_AdvectEuler(DM da,Vec velocity,PetscReal step,int npoints,MPntStd marker[])
{
	Vec             Lvelocity;
	PetscScalar     *LA_velocity;
	PetscScalar     el_velocity[Q2_NODES_PER_EL_3D*NSD];
	PetscInt        e,i;
	PetscScalar     Ni_p[Q2_NODES_PER_EL_3D],vel_p[NSD];
	PetscInt        nel,nen_u;
	const PetscInt  *elnidx_u;
	PetscInt        vel_el_lidx[U_BASIS_FUNCTIONS*3];
	int             p,wil;
	PetscErrorCode  ierr;

    
	PetscFunctionBegin;
	/* scatter velocity to local vector */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p = &marker[p];
		
		wil   = marker_p->wil;
		e     = wil;
		if (wil < 0) { SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Point[%d] has wil_e < 0", wil ); }
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = DMDAGetVectorElementFieldQ2_3D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		P3D_ConstructNi_Q2_3D(marker_p->xi,Ni_p);
		
		vel_p[0] = vel_p[1] = vel_p[2] = 0.0;
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			vel_p[0] += Ni_p[i] * el_velocity[NSD*i+0];
			vel_p[1] += Ni_p[i] * el_velocity[NSD*i+1];
			vel_p[2] += Ni_p[i] * el_velocity[NSD*i+2];
		}
		
		marker_p->coor[0] = marker_p->coor[0] + step * vel_p[0];
		marker_p->coor[1] = marker_p->coor[1] + step * vel_p[1];
		marker_p->coor[2] = marker_p->coor[2] + step * vel_p[2];
	}
	
    ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_ComputeCourantStep"
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep(DM da,Vec velocity,PetscReal *step)
{
	Vec             Lvelocity, gcoords;
	PetscScalar     *LA_velocity, *LA_coords;
	PetscScalar     el_coords[Q2_NODES_PER_EL_3D*NSD];
	PetscScalar     el_velocity[Q2_NODES_PER_EL_3D*NSD];
	PetscInt        e;
	PetscInt        nel,nen_u,ii;
	const PetscInt  *elnidx_u;
	DM              cda;
	PetscReal       dt_min_local, dt_min;
	MPI_Comm        comm;
	PetscReal       hx,hy,hz,vavg[NSD],dtx,dty,dtz;
	PetscReal       xmin,xmax,ymin,ymax,zmin,zmax,coor;
	PetscErrorCode  ierr;
	
    
	PetscFunctionBegin;
    /* setup for coords */
    ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
    ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	/* setup velocity */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	dt_min_local = 1.0e32;
	for (e=0; e<nel; e++) {
		
		/* get coords for the element */
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx_u[nen_u*e],LA_coords);CHKERRQ(ierr);
		/* get velocity */
		ierr = DMDAGetVectorElementFieldQ2_3D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		/* find min,max x,y,z */
		xmin = ymin = zmin =  1.0e32;
		xmax = ymax = zmax = -1.0e32;
		vavg[0] = vavg[1] = vavg[2] = 0.0;
		for (ii=0; ii<Q2_NODES_PER_EL_3D; ii++) {
			coor = el_coords[3*ii+0];
			if (coor < xmin) { xmin = coor; }
			if (coor > xmax) { xmax = coor; }
			
			coor = el_coords[3*ii+1];
			if (coor < ymin) { ymin = coor; }
			if (coor > ymax) { ymax = coor; }
			
			coor = el_coords[3*ii+2];
			if (coor < zmin) { zmin = coor; }
			if (coor > zmax) { zmax = coor; }
			
			vavg[0] = vavg[0] + el_velocity[3*ii+0];
			vavg[1] = vavg[1] + el_velocity[3*ii+1];
			vavg[2] = vavg[2] + el_velocity[3*ii+2];
		}
		vavg[0] = vavg[0]/(PetscReal)(Q2_NODES_PER_EL_3D);
		vavg[1] = vavg[1]/(PetscReal)(Q2_NODES_PER_EL_3D);
		vavg[2] = vavg[2]/(PetscReal)(Q2_NODES_PER_EL_3D);
		
		hx = xmax - xmin;
		hy = ymax - ymin;
		hz = zmax - zmin;
		
		dtx = fabs( hx / vavg[0] );
		dty = fabs( hy / vavg[1] );
		dtz = fabs( hz / vavg[2] );
		
		if (dtx < dt_min_local) { dt_min_local = dtx; }
		if (dty < dt_min_local) { dt_min_local = dty; }
		if (dtz < dt_min_local) { dt_min_local = dtz; }
	}
	
    ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
    ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	
	ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&dt_min_local,&dt_min,1,MPIU_REAL,MPI_MIN,comm);CHKERRQ(ierr);
	
	*step = dt_min;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdateProperties_MPntStd"
PetscErrorCode SwarmUpdateProperties_MPntStd(DataBucket db,pTatinCtx ctx,Vec X)
{
	BTruth found;
	
	PetscFunctionBegin;
	DataBucketQueryDataFieldByName(db,MPntStd_classname,&found);
	if (found == BFALSE) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cannot find DataField with name %s \n", MPntStd_classname );
	}
	PetscFunctionReturn(0);
}

/* ADVECT MARKERS */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateGlobalCoordinates"
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates(DataBucket materialpoints,DM dav,Vec velocity,PetscReal dt)
{
	PetscErrorCode ierr;
	int            npoints;
	MPntStd        *mp_std;
	DataField      PField;
	
	
	PetscFunctionBegin;
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvGlobalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
	ierr = MaterialPointStd_AdvectEuler(dav,velocity,dt,npoints,mp_std);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvGlobalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* UPDATE local coordinates */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateLocalCoordinates"
PetscErrorCode MaterialPointStd_UpdateLocalCoordinates(DataBucket materialpoints,DM dav)
{
	PetscErrorCode ierr;
	int            npoints;
	MPntStd        *mp_std;
	DataField      PField;
	PetscReal      tolerance;
	int            max_its;
	PetscBool      use_nonzero_guess,monitor;
	DM             cda;
	Vec            gcoords;
	PetscScalar    *LA_gcoords;
	const PetscInt *elnidx_u;
	PetscInt       nel,nen_u;
	PetscInt       lmx,lmy,lmz;
	

	PetscFunctionBegin;
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvLocalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
	/* get marker fields */
	DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dav,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	
	/* point location parameters */
	tolerance         = 1.0e-10;
	max_its           = 10;
	use_nonzero_guess = PETSC_TRUE;
	monitor           = PETSC_FALSE;
	
	InverseMappingDomain_3dQ2(tolerance, max_its,
                              use_nonzero_guess,
                              monitor,
                              (const PetscReal*)LA_gcoords, (const PetscInt)lmx,(const PetscInt)lmy,(const PetscInt)lmz, elnidx_u,
                              npoints, mp_std );
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvLocalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* remove all points which didn't find a home */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_Removal"
PetscErrorCode MaterialPointStd_Removal(DataBucket materialpoints)
{
	int            p,npoints,escaped;
	MPntStd        *mp_std;
	DataField      PField;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvRemoval,0,0,0,0);CHKERRQ(ierr);
	/* get marker fields */
	DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
	escaped = 0;
	for (p=0; p<npoints; p++) {
		if (mp_std[p].wil == -1) {
			escaped++;
		}
	}
	
	/* remove points which left processor */
/*
	if (escaped != 0) {
		PetscPrintf(PETSC_COMM_SELF,"  *** MPntStd removal: Identified %d points which are not contained on subdomain (after communication) \n", escaped );
	}
*/
	if (escaped != 0) {
		for (p=0; p<npoints; p++) {
			if (mp_std[p].wil == -1) {
				//printf("############################ Point %d not located in domain (%1.6e , %1.6e) ############################# \n",p,mp_std[p].coor[0],mp_std[p].coor[1]);
				
				/* kill point */
				DataBucketRemovePointAtIndex(materialpoints,p);
				DataBucketGetSizes(materialpoints,&npoints,0,0); /* you need to update npoints as the list size decreases! */
				p--; /* check replacement point */
			}
		}
	}
	ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvRemoval,0,0,0,0);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_Communication_Generic"
PetscErrorCode SwarmUpdatePosition_Communication_Generic(DataBucket db,DM da,DataEx de)
{
	DataField      PField_std;
	DataField      *PField_all;
	int            f,nfields;
	size_t         sizeof_marker_contents;
	int            p,npoints;
	void           *recv_data;
	void           *data_p;
	PetscMPIInt    n,neighborcount, *neighborranks2;
	PetscInt       recv_length;
	PetscInt       npoints_accepted;
	PetscMPIInt    rank,size;
	MPntStd        *marker_std;
	PetscErrorCode ierr;
	
    
	PetscFunctionBegin;
	/* communucate */
	ierr = MPI_Comm_size(PetscObjectComm((PetscObject)da),&size);CHKERRQ(ierr);
	if (size == 1) {
		PetscFunctionReturn(0);
	}
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvCommunication,0,0,0,0);CHKERRQ(ierr);
	
	ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)da),&rank);CHKERRQ(ierr);
	
	neighborcount  = de->n_neighbour_procs;
	neighborranks2 = de->neighbour_procs;
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataBucketGetDataFields(db,&nfields,&PField_all);
	for (f=1; f<nfields; f++) {
		DataFieldGetAccess(PField_all[f]);
	}
	
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	DataBucketGetSizes(db,&npoints,0,0);
	
	/* figure out how many points left processor */
	ierr = DataExInitializeSendCount(de);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		PetscBool onproc;
		MPntStd   *marker;
		
		DataFieldAccessPoint(PField_std,p,(void**)&marker);
		onproc = PETSC_TRUE;
		if (marker->wil == -1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc == PETSC_FALSE) {
			for (n=0; n<neighborcount; n++) {
				ierr = DataExAddToSendCount( de, neighborranks2[n], 1 );CHKERRQ(ierr);
			}
		}
	}
	ierr = DataExFinalizeSendCount(de);CHKERRQ(ierr);
	
	DataFieldRestoreAccess(PField_std);
	
	/* pack points which left processor */
	marker_std    = PField_std->data; /* should write a function to do this */
	
	sizeof_marker_contents = 0;
	sizeof_marker_contents = sizeof_marker_contents + sizeof(MPntStd);
	for (f=1; f<nfields; f++) {
		sizeof_marker_contents += PField_all[f]->atomic_size;
	}
	
	ierr = DataExPackInitialize(de,sizeof_marker_contents);CHKERRQ(ierr);
	
	/* allocate a temporary buffer whihc is large enough to store all marker fields from an individual point,p */
	ierr = PetscMalloc(sizeof_marker_contents,&data_p);CHKERRQ(ierr);
	
	for (p=0; p<npoints; p++) {
		PetscBool   onproc;
		MPntStd     *marker_p;
		void        *marker_p_prop;
		size_t      atomic_size,offset;
		
		/* access fields from the bucket */
		marker_p     = &marker_std[p];
		
		onproc = PETSC_TRUE;
		if (marker_p->wil == -1) {
			onproc = PETSC_FALSE;
			/* pack together */
			//(void*)((char*)field->data + index*field->atomic_size)
			
			/* copy all fields from the bucket into a temporary buffer */
			PetscMemcpy((void*)data_p,                         marker_p,    sizeof(MPntStd));
			
			offset = sizeof(MPntStd);
			for (f=1; f<nfields; f++) {
				atomic_size = PField_all[f]->atomic_size;
				DataFieldAccessPoint(PField_all[f],p,&marker_p_prop);
				
				PetscMemcpy((void*)((char*)data_p+offset),marker_p_prop,atomic_size);
				offset = offset + atomic_size;
			}
		}
		
		if (onproc == PETSC_FALSE) {
			for (n=0; n<neighborcount; n++) {
				ierr = DataExPackData( de, neighborranks2[n], 1,(void*)data_p );CHKERRQ(ierr);
			}
		}
	}		
	for (f=1; f<nfields; f++) {
		DataFieldRestoreAccess(PField_all[f]);
	}
	ierr = DataExPackFinalize(de);CHKERRQ(ierr);
	
	/* remove points which left processor */
	DataBucketGetSizes(db,&npoints,0,0);
	DataFieldGetAccess(PField_std);
	for (p=0; p<npoints; p++) {
		PetscBool onproc;
		MPntStd   *marker_p;
		
		DataFieldAccessPoint(PField_std,p,(void**)&marker_p);
		onproc = PETSC_TRUE;
		if (marker_p->wil == -1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc == PETSC_FALSE) { 
			/* kill point */
			DataBucketRemovePointAtIndex(db,p);
			DataBucketGetSizes(db,&npoints,0,0); /* you need to update npoints as the list size decreases! */
			p--; /* check replacement point */
		}
	}		
	DataFieldRestoreAccess(PField_std);
	
	// START communicate //
	ierr = DataExBegin(de);CHKERRQ(ierr);
	ierr = DataExEnd(de);CHKERRQ(ierr);
	// END communicate //
	
	// receive, if i own them, add new points to list //
	ierr = DataExGetRecvData( de, &recv_length, (void**)&recv_data );CHKERRQ(ierr);
    /*
	{
		PetscInt totalsent;
		ierr = MPI_Allreduce(&recv_length,&totalsent,1,MPIU_INT,MPI_SUM,de->comm);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  DataEx: total points sent = %D \n", totalsent);
	}
	*/
	/* update the local coordinates and cell owner for all recieved points */
	{
		DM cda;
		Vec gcoords;
		PetscScalar *LA_gcoords;
		PetscReal tolerance;
		PetscInt max_its;
		PetscBool use_nonzero_guess,monitor;
		PetscInt lmx,lmy,lmz;
		PetscInt nel,nen_u;
		MPntStd *marker_p;
		const PetscInt *elnidx_u;
		
		/* setup for coords */
		ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
		ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
		
		ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		
		/* point location parameters */
		tolerance         = 1.0e-10;
		max_its           = 10;
		use_nonzero_guess = PETSC_FALSE; /* for markers sent across processors, it is necessary to NOT use the last known values! */
		monitor           = PETSC_FALSE;
		
		/* 
		 NOTE: My algorithm will break down if you try to use an non-zero initial guess
		 if the marker has been sent to another cpu. This could be fixed, however if would
		 require additional logic to be added into the function InverseMappingDomain_2dQ2()
		 
		 It seems more reasonable to simply assume that if the marker was sent to another cpu,
		 the value of wil currently stored on the marker is completely meaningless and thus we
		 should ALWAYS use a zero initial guess (i.e. assume we know nothing)
		 */
		
		for (p=0; p<recv_length; p++) {
			marker_p = (MPntStd*)( (char*)recv_data + p*(sizeof_marker_contents) );
			
			InverseMappingDomain_3dQ2(tolerance, max_its,
                                      use_nonzero_guess,
                                      monitor,
                                      (const PetscReal*)LA_gcoords, (const PetscInt)lmx,(const PetscInt)lmy,(const PetscInt)lmz, elnidx_u,
                                      1, marker_p );
		}
		
		ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	}
	
	/* accepte all points living locally */
	npoints_accepted = 0;
	for (p=0; p<recv_length; p++) {
		PetscBool   onproc;
		MPntStd     *marker_p;
		size_t      offset;
		
		offset   = 0;
		marker_p = (MPntStd*)(     (char*)recv_data + p*(sizeof_marker_contents) + offset);
		
		onproc = PETSC_TRUE;
		if (marker_p->wil == -1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc == PETSC_TRUE) {
			int end;
			
			DataBucketAddPoint(db);
			DataBucketGetSizes(db,&end,0,0);
			end = end - 1;
			
			offset = 0;
			for (f=0; f<nfields; f++) {
				void *data_p = (void*)( (char*)recv_data + p*(sizeof_marker_contents) + offset );
				
				DataFieldInsertPoint(PField_all[f], end, (void*)data_p );
				
				offset = offset + PField_all[f]->atomic_size;
			}
			npoints_accepted++;
		}
	}	
	
	ierr = PetscFree(data_p);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvCommunication,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateCoordinates"
PetscErrorCode MaterialPointStd_UpdateCoordinates(DataBucket materialpoints,DM dav,DataEx de)
{
	PetscErrorCode ierr;
	PetscLogDouble t0,t1;
#ifdef PTAT3D_PROFILE_SwarmUpdatePosition
	PetscLogDouble tlocal[3],tgmax[3],tgmin[3];
	long int       npoints_global_init,npoints_global_fin;
#endif
    
	PetscFunctionBegin;
	PetscTime(&t0);
	ierr = MaterialPointStd_UpdateLocalCoordinates(materialpoints,dav);CHKERRQ(ierr);
	PetscTime(&t1);
	tlocal[0] = t1 - t0;
	
	PetscTime(&t0);
#ifdef PTAT3D_PROFILE_SwarmUpdatePosition
	DataBucketGetGlobalSizes(de->comm,materialpoints,&npoints_global_init,NULL,NULL);
#endif
	ierr = SwarmUpdatePosition_Communication_Generic(materialpoints,dav,de);CHKERRQ(ierr);
#ifdef PTAT3D_PROFILE_SwarmUpdatePosition
	DataBucketGetGlobalSizes(de->comm,materialpoints,&npoints_global_fin,NULL,NULL);
#endif
	PetscTime(&t1);
	tlocal[1] = t1 - t0;

	PetscTime(&t0);
	ierr = MaterialPointStd_Removal(materialpoints);CHKERRQ(ierr);
	PetscTime(&t1);
	tlocal[2] = t1 - t0;

#ifdef PTAT3D_PROFILE_SwarmUpdatePosition
	ierr = MPI_Allreduce(tlocal,tgmax,3,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(tlocal,tgmin,3,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
    
	PetscPrintf(PETSC_COMM_WORLD,"=========== MaterialPointStd_UpdateCoordinates ============= \n");
	PetscPrintf(PETSC_COMM_WORLD,"      Number of material points                     %ld (init) / %ld (final)\n",npoints_global_init,npoints_global_fin);
	PetscPrintf(PETSC_COMM_WORLD,"    	MaterialPointStd_UpdateLocalCoordinates   %10.2e (sec) efficiency %1.1lf%%\n",tgmax[0],100.0 - 100.0*(tgmax[0]-tgmin[0])/tgmax[0]);
	PetscPrintf(PETSC_COMM_WORLD,"    	SwarmUpdatePosition_Communication_Generic %10.2e (sec) efficiency %1.1lf%%\n",tgmax[1],100.0 - 100.0*(tgmax[1]-tgmin[1])/tgmax[1]);
	PetscPrintf(PETSC_COMM_WORLD,"    	MaterialPointStd_Removal                  %10.2e (sec) efficiency %1.1lf%%\n",tgmax[2],100.0 - 100.0*(tgmax[2]-tgmin[2])/tgmax[2]);
	PetscPrintf(PETSC_COMM_WORLD,"=========== ================================== ============= \n");
#endif
    
	PetscFunctionReturn(0);
}

/* ADVECT MARKERS USING RK2,RK3,RK4 */
typedef struct {
    double coor[3];
    double k1[3];
    double k2[3];
} MaterialPointRK2;
const char MaterialPointRK2ClassName[] = "MaterialPointRK2";

typedef struct {
    double coor[3];
    double k1[3];
    double k2[3];
    double k3[3];
} MaterialPointRK3;
const char MaterialPointRK3ClassName[] = "MaterialPointRK3";

typedef struct {
    double coor[3];
    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
} MaterialPointRK4;
const char MaterialPointRK4ClassName[] = "MaterialPointRK4";

#undef __FUNCT__
#define __FUNCT__ "MaterialPointRK4_Interp"
PetscErrorCode MaterialPointRK4_Interp(DM da,Vec velocity,RKOrder order,PetscInt stage,int npoints,MPntStd marker[],
                                       MaterialPointRK2 markerrk2[],
                                       MaterialPointRK3 markerrk3[],
                                       MaterialPointRK4 markerrk4[])
{
	Vec             Lvelocity;
	PetscScalar     *LA_velocity;
	PetscScalar     el_velocity[Q2_NODES_PER_EL_3D*NSD];
	PetscInt        e,i;
	PetscScalar     Ni_p[Q2_NODES_PER_EL_3D],vel_p[NSD];
	PetscInt        nel,nen_u;
	const PetscInt  *elnidx_u;
	PetscInt        vel_el_lidx[U_BASIS_FUNCTIONS*3];
	int             p,wil;
	PetscErrorCode  ierr;
    
    
	PetscFunctionBegin;
	/* scatter velocity to local vector */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p = &marker[p];
		
		wil   = marker_p->wil;
		e     = wil;
		if (wil < 0) { SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Point[%d] has wil_e < 0", wil ); }
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = DMDAGetVectorElementFieldQ2_3D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		P3D_ConstructNi_Q2_3D(marker_p->xi,Ni_p);
		
		vel_p[0] = vel_p[1] = vel_p[2] = 0.0;
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			vel_p[0] += Ni_p[i] * el_velocity[NSD*i+0];
			vel_p[1] += Ni_p[i] * el_velocity[NSD*i+1];
			vel_p[2] += Ni_p[i] * el_velocity[NSD*i+2];
		}
        
        switch (order) {
                
            case RK_ORDER_2:
                switch (stage) {
                    case 1:
                        markerrk2[p].k1[0] = vel_p[0];
                        markerrk2[p].k1[1] = vel_p[1];
                        markerrk2[p].k1[2] = vel_p[2];
                        break;
                    case 2:
                        markerrk2[p].k2[0] = vel_p[0];
                        markerrk2[p].k2[1] = vel_p[1];
                        markerrk2[p].k2[2] = vel_p[2];
                        break;
                }
                break;
                
                
            case RK_ORDER_3:
                switch (stage) {
                    case 1:
                        markerrk3[p].k1[0] = vel_p[0];
                        markerrk3[p].k1[1] = vel_p[1];
                        markerrk3[p].k1[2] = vel_p[2];
                        break;
                    case 2:
                        markerrk3[p].k2[0] = vel_p[0];
                        markerrk3[p].k2[1] = vel_p[1];
                        markerrk3[p].k2[2] = vel_p[2];
                        break;
                    case 3:
                        markerrk3[p].k3[0] = vel_p[0];
                        markerrk3[p].k3[1] = vel_p[1];
                        markerrk3[p].k3[2] = vel_p[2];
                        break;
                }
                break;
                
            case RK_ORDER_4:
                switch (stage) {
                    case 1:
                        markerrk4[p].k1[0] = vel_p[0];
                        markerrk4[p].k1[1] = vel_p[1];
                        markerrk4[p].k1[2] = vel_p[2];
                        break;
                    case 2:
                        markerrk4[p].k2[0] = vel_p[0];
                        markerrk4[p].k2[1] = vel_p[1];
                        markerrk4[p].k2[2] = vel_p[2];
                        break;
                    case 3:
                        markerrk4[p].k3[0] = vel_p[0];
                        markerrk4[p].k3[1] = vel_p[1];
                        markerrk4[p].k3[2] = vel_p[2];
                        break;
                    case 4:
                        markerrk4[p].k4[0] = vel_p[0];
                        markerrk4[p].k4[1] = vel_p[1];
                        markerrk4[p].k4[2] = vel_p[2];
                        break;
                }
                break;
        }
	}
	
    ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_Communication_Generic_2"
PetscErrorCode SwarmUpdatePosition_Communication_Generic_2(DataBucket db1,DataBucket db2,DM da,DataEx de)
{
	DataField      PField_std;
	size_t         sizeof_marker_contents,sizeof_db1_contents,sizeof_db2_contents;
	int            p,npoints;
	void           *recv_data;
	PetscMPIInt    n,neighborcount,*neighborranks2;
	PetscInt       recv_length;
	PetscInt       npoints_accepted;
	PetscMPIInt    rank,size;
	MPntStd        *marker_std;
    void           *dbuf,*dbuf1,*dbuf2;
	PetscErrorCode ierr;
	
    
	PetscFunctionBegin;
    
	ierr = MPI_Comm_size(PetscObjectComm((PetscObject)da),&size);CHKERRQ(ierr);
	if (size == 1) {
		PetscFunctionReturn(0);
	}
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvCommunication,0,0,0,0);CHKERRQ(ierr);
	
	ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)da),&rank);CHKERRQ(ierr);
	
	neighborcount  = de->n_neighbour_procs;
	neighborranks2 = de->neighbour_procs;
	
    DataBucketCreatePackedArray(db1,&sizeof_db1_contents,&dbuf1);
    DataBucketCreatePackedArray(db2,&sizeof_db2_contents,&dbuf2);
    sizeof_marker_contents = sizeof_db1_contents + sizeof_db2_contents;
    PetscMalloc(sizeof_marker_contents,&dbuf);
    
	DataBucketGetDataFieldByName(db1,MPntStd_classname,&PField_std);
    DataFieldGetEntries(PField_std,(void**)&marker_std);
    
	DataBucketGetSizes(db1,&npoints,0,0);
	
	/* figure out how many points left processor */
	ierr = DataExInitializeSendCount(de);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		if (marker_std[p].wil == -1) {
			for (n=0; n<neighborcount; n++) {
				ierr = DataExAddToSendCount( de, neighborranks2[n], 1 );CHKERRQ(ierr);
			}
		}
	}
	ierr = DataExFinalizeSendCount(de);CHKERRQ(ierr);
	
    DataFieldRestoreEntries(PField_std,(void**)&marker_std);
	
	ierr = DataExPackInitialize(de,sizeof_marker_contents);CHKERRQ(ierr);
	
    DataFieldGetEntries(PField_std,(void**)&marker_std);
	for (p=0; p<npoints; p++) {
		MPntStd     *marker_p;
		
		/* access fields from the bucket */
		marker_p     = &marker_std[p];
		
		if (marker_p->wil == -1) {
            
            DataBucketFillPackedArray(db1,p,dbuf1);
            DataBucketFillPackedArray(db2,p,dbuf2);
            
            PetscMemcpy(dbuf,dbuf1,sizeof_db1_contents);
            PetscMemcpy((void*)((char*)dbuf+sizeof_db1_contents),dbuf2,sizeof_db2_contents);
            
			for (n=0; n<neighborcount; n++) {
				ierr = DataExPackData( de, neighborranks2[n], 1,(void*)dbuf );CHKERRQ(ierr);
            }
        }
	}
    DataFieldRestoreEntries(PField_std,(void**)&marker_std);
    
	ierr = DataExPackFinalize(de);CHKERRQ(ierr);
	
	/* remove points which left processor */
	DataBucketGetSizes(db1,&npoints,0,0);
	DataFieldGetAccess(PField_std);
	for (p=0; p<npoints; p++) {
		MPntStd   *marker_p;
		
		DataFieldAccessPoint(PField_std,p,(void**)&marker_p);
		if (marker_p->wil == -1) {
			/* kill point */
			DataBucketRemovePointAtIndex(db1,p);
			DataBucketRemovePointAtIndex(db2,p);
			DataBucketGetSizes(db1,&npoints,0,0); /* you need to update npoints as the list size decreases! */
			p--; /* check replacement point */
		}
	}
	DataFieldRestoreAccess(PField_std);
	
	// START communicate //
	ierr = DataExBegin(de);CHKERRQ(ierr);
	ierr = DataExEnd(de);CHKERRQ(ierr);
	// END communicate //
    
	// receive, if i own them, add new points to list //
	ierr = DataExGetRecvData( de, &recv_length, (void**)&recv_data );CHKERRQ(ierr);
    
	/* update the local coordinates and cell owner for all recieved points */
	{
		DM             cda;
		Vec            gcoords;
		PetscScalar    *LA_gcoords;
		PetscReal      tolerance;
		PetscInt       max_its;
		PetscBool      use_nonzero_guess,monitor;
		PetscInt       lmx,lmy,lmz,nel,nen_u;
		MPntStd        *marker_p;
		const PetscInt *elnidx_u;
		
		/* setup for coords */
		ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
		ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
		
		ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		
		/* point location parameters */
		tolerance         = 1.0e-10;
		max_its           = 10;
		use_nonzero_guess = PETSC_FALSE; /* for markers sent across processors, it is necessary to NOT use the last known values! */
		monitor           = PETSC_FALSE;
		
		for (p=0; p<recv_length; p++) {
			marker_p = (MPntStd*)( (char*)recv_data + p*(sizeof_marker_contents) );
			
			InverseMappingDomain_3dQ2(tolerance, max_its,
                                      use_nonzero_guess,
                                      monitor,
                                      (const PetscReal*)LA_gcoords, (const PetscInt)lmx,(const PetscInt)lmy,(const PetscInt)lmz, elnidx_u,
                                      1, marker_p );
		}
		
		ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	}
	
	/* accept all points living locally */
	npoints_accepted = 0;
	for (p=0; p<recv_length; p++) {
        MPntStd *marker_p;
		void  *data1;
		void  *data2;
		
		marker_p = (MPntStd*)( (char*)recv_data + p*(sizeof_marker_contents));
		data1    = (void*)(    (char*)recv_data + p*(sizeof_marker_contents));
		data2    = (void*)(    (char*)recv_data + p*(sizeof_marker_contents) + sizeof_db1_contents);
		
		if (marker_p->wil != -1) {
			int end;
			
			DataBucketAddPoint(db1);
			DataBucketAddPoint(db2);
			DataBucketGetSizes(db1,&end,0,0);
			end = end - 1;
			
            DataBucketInsertPackedArray(db1,end,data1);
            DataBucketInsertPackedArray(db2,end,data2);
            
			npoints_accepted++;
		}
	}
	
	ierr = PetscFree(dbuf);CHKERRQ(ierr);
    DataBucketDestroyPackedArray(db1,&dbuf1);
    DataBucketDestroyPackedArray(db2,&dbuf2);
	ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvCommunication,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_Removal_2"
PetscErrorCode MaterialPointStd_Removal_2(DataBucket materialpoints,DataBucket db2)
{
	int            p,npoints,escaped;
	MPntStd        *mp_std;
	DataField      PField;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvRemoval,0,0,0,0);CHKERRQ(ierr);
	/* get marker fields */
	DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
	escaped = 0;
	for (p=0; p<npoints; p++) {
		if (mp_std[p].wil == -1) {
			escaped++;
		}
	}
	
	/* remove points which left processor */
	if (escaped != 0) {
		for (p=0; p<npoints; p++) {
			if (mp_std[p].wil == -1) {
				/* kill point */
				DataBucketRemovePointAtIndex(materialpoints,p);
				DataBucketRemovePointAtIndex(db2,p);
				DataBucketGetSizes(materialpoints,&npoints,0,0); /* you need to update npoints as the list size decreases! */
				p--; /* check replacement point */
			}
		}
	}
	ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvRemoval,0,0,0,0);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateCoordinates_2"
PetscErrorCode MaterialPointStd_UpdateCoordinates_2(DataBucket materialpoints,DataBucket db2,DM dav,DataEx de)
{
	PetscErrorCode ierr;
    
	PetscFunctionBegin;

	ierr = MaterialPointStd_UpdateLocalCoordinates(materialpoints,dav);CHKERRQ(ierr);
	ierr = SwarmUpdatePosition_Communication_Generic_2(materialpoints,db2,dav,de);CHKERRQ(ierr);
	ierr = MaterialPointStd_Removal_2(materialpoints,db2);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateGlobalCoordinates_RK4"
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates_RK4(DataBucket materialpoints,DataEx de,DM dav,Vec velocity,PetscReal dt)
{
	PetscErrorCode   ierr;
	int              p,npoints;
	MPntStd          *mp_std;
	MaterialPointRK4 *mp_rk4;
	DataField        PField,PFieldRK4;
	DataBucket       dbrk4;
    
	
	PetscFunctionBegin;
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvGlobalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
    
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
    
    DataBucketCreate(&dbrk4);
	DataBucketRegisterField(dbrk4,MaterialPointRK4ClassName,sizeof(MaterialPointRK4),NULL);
	DataBucketFinalize(dbrk4);
	DataBucketSetSizes(dbrk4,npoints,100);
    
    DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
    DataBucketGetDataFieldByName(dbrk4, MaterialPointRK4ClassName ,&PFieldRK4);
	mp_rk4 = PFieldRK4->data;
    for (p=0; p<npoints; p++) {
        mp_rk4[p].coor[0] = mp_std[p].coor[0];
        mp_rk4[p].coor[1] = mp_std[p].coor[1];
        mp_rk4[p].coor[2] = mp_std[p].coor[2];
    }
    
    /* stage 1 : k1 = v(x_n) */
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_4,1,npoints,mp_std,NULL,NULL,mp_rk4);CHKERRQ(ierr);
    
    /* stage 2 : k2 = v(x_n + 0.5 * dt * k1) */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk4[p].coor[0] + 0.5 * dt * mp_rk4[p].k1[0];
        mp_std[p].coor[1] = mp_rk4[p].coor[1] + 0.5 * dt * mp_rk4[p].k1[1];
        mp_std[p].coor[2] = mp_rk4[p].coor[2] + 0.5 * dt * mp_rk4[p].k1[2];
    }
    ierr = MaterialPointStd_UpdateCoordinates_2(materialpoints,dbrk4,dav,de);CHKERRQ(ierr);
    
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	mp_std = PField->data;
	mp_rk4 = PFieldRK4->data;
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_4,2,npoints,mp_std,NULL,NULL,mp_rk4);CHKERRQ(ierr);
    
    /* stage 3 : k3 = v(x_n + 0.5 * dt * k2) */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk4[p].coor[0] + 0.5 * dt * mp_rk4[p].k2[0];
        mp_std[p].coor[1] = mp_rk4[p].coor[1] + 0.5 * dt * mp_rk4[p].k2[1];
        mp_std[p].coor[2] = mp_rk4[p].coor[2] + 0.5 * dt * mp_rk4[p].k2[2];
    }
    ierr = MaterialPointStd_UpdateCoordinates_2(materialpoints,dbrk4,dav,de);CHKERRQ(ierr);
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	mp_std = PField->data;
	mp_rk4 = PFieldRK4->data;
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_4,3,npoints,mp_std,NULL,NULL,mp_rk4);CHKERRQ(ierr);
    
    /* stage 4 : k4 = v(x_n + dt * k3) */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk4[p].coor[0] + dt * mp_rk4[p].k3[0];
        mp_std[p].coor[1] = mp_rk4[p].coor[1] + dt * mp_rk4[p].k3[1];
        mp_std[p].coor[2] = mp_rk4[p].coor[2] + dt * mp_rk4[p].k3[2];
    }
    ierr = MaterialPointStd_UpdateCoordinates_2(materialpoints,dbrk4,dav,de);CHKERRQ(ierr);
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	mp_std = PField->data;
	mp_rk4 = PFieldRK4->data;
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_4,4,npoints,mp_std,NULL,NULL,mp_rk4);CHKERRQ(ierr);
    
    /* final step */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk4[p].coor[0] + (1.0/6.0) * dt * ( mp_rk4[p].k1[0] + 2.0*mp_rk4[p].k2[0] + 2.0*mp_rk4[p].k3[0] + mp_rk4[p].k4[0]);
        mp_std[p].coor[1] = mp_rk4[p].coor[1] + (1.0/6.0) * dt * ( mp_rk4[p].k1[1] + 2.0*mp_rk4[p].k2[1] + 2.0*mp_rk4[p].k3[1] + mp_rk4[p].k4[1]);
        mp_std[p].coor[2] = mp_rk4[p].coor[2] + (1.0/6.0) * dt * ( mp_rk4[p].k1[2] + 2.0*mp_rk4[p].k2[2] + 2.0*mp_rk4[p].k3[2] + mp_rk4[p].k4[2]);
    }
    
	DataBucketDestroy(&dbrk4);
    
    ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvGlobalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateGlobalCoordinates_RK2"
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates_RK2(DataBucket materialpoints,DataEx de,DM dav,Vec velocity,PetscReal dt)
{
	PetscErrorCode    ierr;
	int               p,npoints;
	MPntStd           *mp_std;
	MaterialPointRK2  *mp_rk;
	DataField         PField,PFieldRK;
	DataBucket        dbrk;
	
    
	PetscFunctionBegin;
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvGlobalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
    
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
    
    DataBucketCreate(&dbrk);
	DataBucketRegisterField(dbrk,MaterialPointRK2ClassName,sizeof(MaterialPointRK2),NULL);
	DataBucketFinalize(dbrk);
	DataBucketSetSizes(dbrk,npoints,100);
    
    DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
    DataBucketGetDataFieldByName(dbrk, MaterialPointRK2ClassName ,&PFieldRK);
	mp_rk = PFieldRK->data;
    for (p=0; p<npoints; p++) {
        mp_rk[p].coor[0] = mp_std[p].coor[0];
        mp_rk[p].coor[1] = mp_std[p].coor[1];
        mp_rk[p].coor[2] = mp_std[p].coor[2];
    }
    
    /* stage 1 : k1 = v(x_n) */
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_2,1,npoints,mp_std,mp_rk,NULL,NULL);CHKERRQ(ierr);
    
    /* stage 2 : k2 = v(x_n + 0.5 * dt * k1) */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk[p].coor[0] + 0.5 * dt * mp_rk[p].k1[0];
        mp_std[p].coor[1] = mp_rk[p].coor[1] + 0.5 * dt * mp_rk[p].k1[1];
        mp_std[p].coor[2] = mp_rk[p].coor[2] + 0.5 * dt * mp_rk[p].k1[2];
    }
    ierr = MaterialPointStd_UpdateCoordinates_2(materialpoints,dbrk,dav,de);CHKERRQ(ierr);
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	mp_std = PField->data;
	mp_rk  = PFieldRK->data;
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_2,2,npoints,mp_std,mp_rk,NULL,NULL);CHKERRQ(ierr);
    
    /* final step */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk[p].coor[0] + dt * mp_rk[p].k2[0];
        mp_std[p].coor[1] = mp_rk[p].coor[1] + dt * mp_rk[p].k2[1];
        mp_std[p].coor[2] = mp_rk[p].coor[2] + dt * mp_rk[p].k2[2];
    }
    
	DataBucketDestroy(&dbrk);
    
    ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvGlobalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateGlobalCoordinates_RK3"
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinates_RK3(DataBucket materialpoints,DataEx de,DM dav,Vec velocity,PetscReal dt)
{
	PetscErrorCode    ierr;
	int               p,npoints;
	MPntStd           *mp_std;
	MaterialPointRK3  *mp_rk;
	DataField         PField,PFieldRK;
	DataBucket        dbrk;
	
    
	PetscFunctionBegin;
	ierr = PetscLogEventBegin(PTATIN_MaterialPointAdvGlobalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
    
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
    
    DataBucketCreate(&dbrk);
	DataBucketRegisterField(dbrk,MaterialPointRK3ClassName,sizeof(MaterialPointRK3),NULL);
	DataBucketFinalize(dbrk);
	DataBucketSetSizes(dbrk,npoints,100);
    
    DataBucketGetDataFieldByName(materialpoints, MPntStd_classname ,&PField);
	mp_std = PField->data;
	
    DataBucketGetDataFieldByName(dbrk, MaterialPointRK3ClassName ,&PFieldRK);
	mp_rk = PFieldRK->data;
    for (p=0; p<npoints; p++) {
        mp_rk[p].coor[0] = mp_std[p].coor[0];
        mp_rk[p].coor[1] = mp_std[p].coor[1];
        mp_rk[p].coor[2] = mp_std[p].coor[2];
    }
    
    /* stage 1 : k1 = v(x_n) */
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_3,1,npoints,mp_std,NULL,mp_rk,NULL);CHKERRQ(ierr);
    
    /* stage 2 : k2 = v(x_n + 0.5 * dt * k1) */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk[p].coor[0] + 0.25 * dt * mp_rk[p].k1[0];
        mp_std[p].coor[1] = mp_rk[p].coor[1] + 0.25 * dt * mp_rk[p].k1[1];
        mp_std[p].coor[2] = mp_rk[p].coor[2] + 0.25 * dt * mp_rk[p].k1[2];
    }
    ierr = MaterialPointStd_UpdateCoordinates_2(materialpoints,dbrk,dav,de);CHKERRQ(ierr);
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	mp_std = PField->data;
	mp_rk  = PFieldRK->data;
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_3,2,npoints,mp_std,NULL,mp_rk,NULL);CHKERRQ(ierr);
    
    /* stage 3 : k2 = v(x_n -(2/3) * dt * k1 + k2) */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk[p].coor[0] - (2.0/3.0) * dt * mp_rk[p].k1[0] + dt * mp_rk[p].k2[0];
        mp_std[p].coor[1] = mp_rk[p].coor[1] - (2.0/3.0) * dt * mp_rk[p].k1[1] + dt * mp_rk[p].k2[1];
        mp_std[p].coor[2] = mp_rk[p].coor[2] - (2.0/3.0) * dt * mp_rk[p].k1[2] + dt * mp_rk[p].k2[2];
    }
    ierr = MaterialPointStd_UpdateCoordinates_2(materialpoints,dbrk,dav,de);CHKERRQ(ierr);
    DataBucketGetSizes(materialpoints,&npoints,NULL,NULL);
	mp_std = PField->data;
	mp_rk  = PFieldRK->data;
    ierr = MaterialPointRK4_Interp(dav,velocity,RK_ORDER_3,3,npoints,mp_std,NULL,mp_rk,NULL);CHKERRQ(ierr);
    
    /* final step */
    for (p=0; p<npoints; p++) {
        mp_std[p].coor[0] = mp_rk[p].coor[0] + (1.0/3.0) * dt * (mp_rk[p].k1[0] + mp_rk[p].k2[0] + mp_rk[p].k3[0]);
        mp_std[p].coor[1] = mp_rk[p].coor[1] + (1.0/3.0) * dt * (mp_rk[p].k1[1] + mp_rk[p].k2[1] + mp_rk[p].k3[1]);
        mp_std[p].coor[2] = mp_rk[p].coor[2] + (1.0/3.0) * dt * (mp_rk[p].k1[2] + mp_rk[p].k2[2] + mp_rk[p].k3[2]);
    }
    
	DataBucketDestroy(&dbrk);
    
    ierr = PetscLogEventEnd(PTATIN_MaterialPointAdvGlobalCoordUpdate,0,0,0,0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStd_UpdateGlobalCoordinatesRK"
PetscErrorCode MaterialPointStd_UpdateGlobalCoordinatesRK(DataBucket materialpoints,RKOrder order,DataEx de,DM dav,Vec velocity,PetscReal dt)
{
    PetscErrorCode ierr;
    
    switch (order) {
        case RK_ORDER_1:
            ierr = MaterialPointStd_UpdateGlobalCoordinates(materialpoints,dav,velocity,dt);CHKERRQ(ierr);
            break;
        case RK_ORDER_2:
            ierr = MaterialPointStd_UpdateGlobalCoordinates_RK2(materialpoints,de,dav,velocity,dt);CHKERRQ(ierr);
            break;
        case RK_ORDER_3:
            ierr = MaterialPointStd_UpdateGlobalCoordinates_RK3(materialpoints,de,dav,velocity,dt);CHKERRQ(ierr);
            break;
        case RK_ORDER_4:
            ierr = MaterialPointStd_UpdateGlobalCoordinates_RK4(materialpoints,de,dav,velocity,dt);CHKERRQ(ierr);
            break;
        default:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Runge-Kutta material point advection only supports RK1, RK2, RK3, RK4");
            break;
    }
    
	PetscFunctionReturn(0);
}



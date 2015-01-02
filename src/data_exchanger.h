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
 **    filename:   data_exchanger.h
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



#ifndef __TATIN_DATA_EXCHANGER_H__
#define __TATIN_DATA_EXCHANGER_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

typedef enum { DEOBJECT_INITIALIZED=0, DEOBJECT_FINALIZED, DEOBJECT_STATE_UNKNOWN } DEObjectState;

typedef struct _p_DataEx* DataEx;
struct  _p_DataEx {
	PetscInt    instance;
	MPI_Comm    comm;
	PetscMPIInt rank;
	
	PetscMPIInt  n_neighbour_procs;
	PetscMPIInt *neighbour_procs; /* [n_neighbour_procs] */
	PetscInt    *messages_to_be_sent; /* [n_neighbour_procs] */
	PetscInt    *message_offsets; /* [n_neighbour_procs] */
	PetscInt    *messages_to_be_recvieved; /* [n_neighbour_procs] */
	size_t       unit_message_size;
	void        *send_message;
	PetscInt     send_message_length;
	void        *recv_message;
	PetscInt     recv_message_length;
	int         *send_tags, *recv_tags;
	PetscInt     total_pack_cnt;
	PetscInt    *pack_cnt; /* [n_neighbour_procs] */
	DEObjectState     topology_status;
	DEObjectState     message_lengths_status;
	DEObjectState     packer_status;
	DEObjectState     communication_status;
	
	MPI_Status  *_stats;
	MPI_Request *_requests;
};


/* OBJECT_STATUS */
//#define OBJECT_INITIALIZED    0
//#define OBJECT_FINALIZED      1
//#define OBJECT_STATE_UNKNOWN  2

extern const char *status_names[];

DataEx DataExCreate(MPI_Comm comm,const PetscInt count);
PetscErrorCode DataExView(DataEx d);
PetscErrorCode DataExDestroy(DataEx d);
PetscErrorCode DataExTopologyInitialize(DataEx d);
PetscErrorCode DataExTopologyAddNeighbour(DataEx d,const PetscMPIInt proc_id);
PetscErrorCode DataExTopologyFinalize(DataEx d);
PetscErrorCode DataExInitializeSendCount(DataEx de);
PetscErrorCode DataExAddToSendCount(DataEx de,const PetscMPIInt proc_id,const PetscInt count);
PetscErrorCode DataExFinalizeSendCount(DataEx de);
PetscErrorCode DataExPackInitialize(DataEx de,size_t unit_message_size);
PetscErrorCode DataExPackData(DataEx de,PetscMPIInt proc_id,PetscInt n,void *data);
PetscErrorCode DataExPackFinalize(DataEx de);
PetscErrorCode DataExBegin(DataEx de);
PetscErrorCode DataExEnd(DataEx de);
PetscErrorCode DataExGetSendData(DataEx de,PetscInt *length,void **send);
PetscErrorCode DataExGetRecvData(DataEx de,PetscInt *length,void **recv);
PetscErrorCode DataExTopologyGetNeighbours(DataEx de,PetscMPIInt *n,PetscInt *neigh[]);


#endif


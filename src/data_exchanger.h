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
 **    Filename:      data_exchanger.h
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
	int         instance;
	MPI_Comm    comm;
	int         rank;
	
	int     n_neighbour_procs;
	int    *neighbour_procs; /* [n_neighbour_procs] */
	int    *messages_to_be_sent; /* [n_neighbour_procs] */
	int    *message_offsets; /* [n_neighbour_procs] */
	int    *messages_to_be_recvieved; /* [n_neighbour_procs] */
	size_t  unit_message_size;
	void   *send_message;
	int     send_message_length;
	void   *recv_message;
	int     recv_message_length;
	int     *send_tags, *recv_tags;
	int     total_pack_cnt;
	int    *pack_cnt; /* [n_neighbour_procs] */
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

DataEx DataExCreate( MPI_Comm comm, const int count );
PetscErrorCode DataExView( DataEx d );
PetscErrorCode DataExDestroy( DataEx d );
PetscErrorCode DataExTopologyInitialize( DataEx d );
PetscErrorCode DataExTopologyAddNeighbour( DataEx d, const int proc_id );
PetscErrorCode DataExTopologyFinalize( DataEx d );
PetscErrorCode DataExInitializeSendCount( DataEx de );
PetscErrorCode DataExAddToSendCount( DataEx de, const int proc_id, const int count );
PetscErrorCode DataExFinalizeSendCount( DataEx de );
PetscErrorCode DataExPackInitialize( DataEx de, size_t unit_message_size );
PetscErrorCode DataExPackData( DataEx de, int proc_id, int n, void *data );
PetscErrorCode DataExPackFinalize( DataEx de );
PetscErrorCode DataExBegin( DataEx de );
PetscErrorCode DataExEnd( DataEx de );
PetscErrorCode DataExGetSendData( DataEx de, int *length, void **send );
PetscErrorCode DataExGetRecvData( DataEx de, int *length, void **recv );



#endif


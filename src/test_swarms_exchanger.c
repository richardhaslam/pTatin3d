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
 **    Filename:      test_swarms_exchanger.c
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

/*
 Information and usage:
 
 To compile as a stand alone example, you need to
 0) Install MPI and PETSc
 1) Create the object file for data_bucket.c and data_exchanger.c
 2) Compile test_swarms.c and link against data_bucket.o data_exchanger.o
 
 [step 1]
 ${MPI_DIR}/bin/mpicc -O0 -g -c data_bucket.c -I. -I${MPI_DIR}/include
 ${MPI_DIR}/bin/mpicc -O0 -g -c data_exchanger.c -I. -I${MPI_DIR}/include
 [step 2]
 ${MPI_DIR}/bin/mpicc -O0 -g -o test_swarms test_swarms.c data_bucket.o data_exchanger.o -I. -I${PETSC_DIR}/include -L${PETSC_DIR}/lib -lpetsc
 
 /Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/bin/mpicc -O0 -g -c data_bucket.c -I. -I/Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/include -I/Users/dmay/software/petsc-3.5.1/include
 /Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/bin/mpicc -O0 -g -c data_exchanger.c -I. -I/Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/include -I/Users/dmay/software/petsc-3.5.1/include
 /Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/bin/mpicc -O0 -g -o test_swarms_exchanger test_swarms_exchanger.c data_bucket.o data_exchanger.o -I. -I/Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/include -I/Users/dmay/software/petsc-3.5.1/include -L/Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/lib -lpetsc
*/

#include <data_bucket.h>
#include <data_exchanger.h>

/* Define structures for storing material point data */
typedef struct {
	double    coor[3];
	long int  pid;
	float     viscosity;
    short int region_id;
} MaterialPoint;

typedef struct {
	float plastic_strain;
    char  failure_type;
} MaterialPointVP;

/* Define names used to identify data structure for material point data */
const char MaterialPointClassName[]   = "DBF_MaterialPoint";
const char MaterialPointVPClassName[] = "DBF_MaterialPointVP";



#undef __FUNCT__
#define __FUNCT__ "SwarmTest_Communication1"
PetscErrorCode SwarmTest_Communication1(MPI_Comm comm)
{
    PetscErrorCode ierr;
    DataEx          data_exchanger;
	DataBucket      db;
	DataField       dbField;
    MaterialPoint   *mp_list,*mp_list_recv;
    int             nproc,rank,i,k,init_size,L_local;
    long int        L_global;
    PetscInt        num_items_recv;
    
    
    ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    
    if (nproc != 2) {
        SETERRQ(comm,PETSC_ERR_SUP,"Example only valid if nproc = 2");
    }
    
    if (rank) {
        printf("[[%s]]\n",__FUNCTION__);
    }
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    
    DataBucketCreate(&db);
	DataBucketRegisterField(db,MaterialPointClassName,sizeof(MaterialPoint),NULL);
	DataBucketFinalize(db);
    
    /* Each processor will define a different number of entries */
    if (rank == 0) { init_size = 7; }
    if (rank == 1) { init_size = 3; }
	DataBucketSetSizes(db,init_size,-1);
    
	DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
	DataFieldGetAccess(dbField);
    /* Get number of entries */
	DataBucketGetSizes(db,&L_local,NULL,NULL);
    
    for (k=0; k<L_local; k++) {
		MaterialPoint *mp_k;
        
		DataFieldAccessPoint(dbField,k,(void**)&mp_k);
        
        /* Set values in the k'th entry */
        mp_k->coor[0] = (double)3*k + 1.1 + (double)(rank*10);
        mp_k->coor[1] = (double)3*k + 2.2 + (double)(rank*10);
        mp_k->coor[2] = (double)3*k + 3.3 + (double)(rank*10);
        
        mp_k->pid       = (long int)k + (double)(rank*100);
        mp_k->viscosity = (float)(k+10) + (float)(rank*10);
        mp_k->region_id = (short int)(rank+1);
    }
    
	DataFieldRestoreAccess(dbField);
    
    /* Examine result */
    DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
    DataFieldGetEntries(dbField,(void**)&mp_list);
	DataBucketGetSizes(db,&L_local,NULL,NULL);
    
    ierr = MPI_Barrier(comm);CHKERRQ(ierr); fflush(stdout);
    if (rank) { printf("[Initial data bucket (pre communication)]\n"); fflush(stdout);} ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    for (i=0; i<nproc; i++) {
        ierr = MPI_Barrier(comm);CHKERRQ(ierr);
        if (rank == i) {
            for (k=0; k<L_local; k++) {
                printf("  [rank %.2d] [%.2d] pid = %.4ld ; viscosity = %1.4e\n",rank,k,mp_list[k].pid,mp_list[k].viscosity); fflush(stdout);
            }
        }
    } fflush(stdout);

    DataFieldRestoreEntries(dbField,(void**)&mp_list);
    
    /* Get total number of entries (summed over all ranks in comm */
	DataBucketGetGlobalSizes(comm,db,&L_global,NULL,NULL);
    
    /* Report parallel information about data bucket */
    DataBucketView(comm,db,"Material Point Coefficients",DATABUCKET_VIEW_STDOUT);
    
    
    
    /*
     Create the data exchanger
     
     DataEx DataExCreate(MPI_Comm comm,const PetscInt count)
       const PetscInt count : Index assigned to data exchange, this is used to generate unique message tags
                              if multiple data exchangers are to be used simulataneously.
    */
    data_exchanger = DataExCreate(comm,0);
    
    
    /* 
     Define communication topology using 3 step process
       [1] Call DataExTopologyInitialize()
       [2] Call DataExTopologyAddNeighbour() multiple times
       [3] DataExTopologyFinalize()
     
     PetscErrorCode DataExTopologyAddNeighbour(DataEx d,const PetscMPIInt proc_id)
       const PetscMPIInt proc_id : The index of a process in communicator used when DataEx was created which is
                                   to labelled as a "neighour" of the current processor.
    */
    ierr = DataExTopologyInitialize(data_exchanger);CHKERRQ(ierr);
    
    /* Configure exchange such that rank0 will send data to rank1 (and vice-versa) */
    if (rank == 0) {  ierr = DataExTopologyAddNeighbour(data_exchanger,1);CHKERRQ(ierr); }
    if (rank == 1) {  ierr = DataExTopologyAddNeighbour(data_exchanger,0);CHKERRQ(ierr); }
    
    ierr = DataExTopologyFinalize(data_exchanger);CHKERRQ(ierr);
    
    /* 
     Send material points with indices {1,3,6} FROM rank 0 TO rank 1
     Send material points with indices {0}     FROM rank 1 TO rank 0
     
     Communication occurs in four stages;
       (a) Establish number of items to be sent and where they are to be sent to
       (b) Pack all data to be sent (order is not important)
       (c) Send data
       (d) Get the buffer containing the data recieved
    */
    
    /*
     Phase a) Establish send counts
     
     PetscErrorCode DataExInitializeSendCount(DataEx de)
     PetscErrorCode DataExAddToSendCount(DataEx de,const PetscMPIInt proc_id,const PetscInt count)
     PetscErrorCode DataExFinalizeSendCount(DataEx de)
    */
    /* Initialze all counters, flush all memory buffers */
    ierr = DataExInitializeSendCount(data_exchanger);CHKERRQ(ierr);
    if (rank == 0) { /* send three items form rank 0 to rank 1 */
        ierr = DataExAddToSendCount(data_exchanger,1,3);CHKERRQ(ierr);
    }
    if (rank == 1) { /* send one item from rank 1 to rank 0 */
        ierr = DataExAddToSendCount(data_exchanger,0,1);CHKERRQ(ierr);
    }
    /* Finalize all counters - no additional items can be sent unless DataExInitializeSendCount() is called again */
    ierr = DataExFinalizeSendCount(data_exchanger);CHKERRQ(ierr);
    
    /*
     Phase b) Pack data into buffers to be sent

     PetscErrorCode DataExPackInitialize(DataEx de,size_t unit_message_size)
     PetscErrorCode DataExPackData(DataEx de,PetscMPIInt proc_id,PetscInt n,void *data)
     PetscErrorCode DataExPackFinalize(DataEx de)
    */
    /* Flush buffers and indicate what the atomic size is of each item to be sent */
    ierr = DataExPackInitialize(data_exchanger,sizeof(MaterialPoint));CHKERRQ(ierr);
    
    /* Pack data into buffers - the number of items packed must match the total number determined from calls to DataExAddToSendCount() */
    DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
    DataFieldGetEntries(dbField,(void**)&mp_list);
    if (rank == 0) {
        /* Items to be sent can be packed one at a time, or all at once. Here we pack markers one at a time */
        /* Insert markers 1, 3, 6 (to be sent to rank 1) */
        ierr = DataExPackData(data_exchanger,1,1,&mp_list[1]);CHKERRQ(ierr);
        ierr = DataExPackData(data_exchanger,1,1,&mp_list[3]);CHKERRQ(ierr);
        ierr = DataExPackData(data_exchanger,1,1,&mp_list[6]);CHKERRQ(ierr);
    }
    if (rank == 1) {
        /* Insert markers 0 (to be sent to rank 0) */
        ierr = DataExPackData(data_exchanger,0,1,&mp_list[0]);CHKERRQ(ierr);
    }
    DataFieldRestoreEntries(dbField,(void**)&mp_list);
    
    /* 
     [DATASTRUCTURE - APPLICATION SPECIFIC] NOTE
     
     Users should consider if they want to remove the items being sent from their data structure.
     This is an application specific choice.
     For a typical material point implementation, sending a marker to another processor indicates
     that the material point has physically moved from one part of the domain to another, and thus
     "ownership" of the material point has changed - hence we should remove the material points which
     are being sent.
    */
    if (rank == 0) {
        /* I remove them in reverse order for convienence */
        DataBucketRemovePointAtIndex(db,6);
        DataBucketRemovePointAtIndex(db,3);
        DataBucketRemovePointAtIndex(db,1);
    }
    if (rank == 1) {
        DataBucketRemovePointAtIndex(db,0);
    }
    
    
    /* Finalize packing - no further data can be added into the buffer after this call */
    ierr = DataExPackFinalize(data_exchanger);CHKERRQ(ierr);
    
    /*
     Phase c) Send the data

     PetscErrorCode DataExBegin(DataEx de)
     PetscErrorCode DataExEnd(DataEx de)
    */
    ierr = DataExBegin(data_exchanger);CHKERRQ(ierr);
    ierr = DataExEnd(data_exchanger);CHKERRQ(ierr);
    
    /*
     Phase d) Fetch the recv buffer and do something with it

     PetscErrorCode DataExGetRecvData(DataEx de,PetscInt *length,void **recv)
    */
    ierr = DataExGetRecvData(data_exchanger,&num_items_recv,(void**)&mp_list_recv);CHKERRQ(ierr);
    
    
    /*
     [DATASTRUCTURE - APPLICATION SPECIFIC] NOTE
     
     The newly received data needs to be inserted into the users data structure.
     We will insert new material points at the end of the array within our data bucket
    */
	DataBucketGetSizes(db,&L_local,NULL,NULL);
	DataBucketSetSizes(db,L_local + num_items_recv,-1);
    DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
    for (k=0; k<num_items_recv; k++) {
        DataFieldInsertPoint(dbField,L_local+k,&mp_list_recv[k]);
    }

    
    /* View result */
    DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
    DataFieldGetEntries(dbField,(void**)&mp_list);
	DataBucketGetSizes(db,&L_local,NULL,NULL);
    
    ierr = MPI_Barrier(comm);CHKERRQ(ierr); fflush(stdout);
    if (rank) { printf("[Final data bucket (post communication)]\n"); fflush(stdout); } ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    for (i=0; i<nproc; i++) {
        ierr = MPI_Barrier(comm);CHKERRQ(ierr);
        if (rank == i) {
            for (k=0; k<L_local; k++) {
                printf("  [*rank %.2d] [%.2d] pid = %.4ld ; viscosity = %1.4e\n",rank,k,mp_list[k].pid,mp_list[k].viscosity); fflush(stdout);
            }
        }
    } fflush(stdout);
    
    DataFieldRestoreEntries(dbField,(void**)&mp_list);
    
    
    
    /* View information about data_exchanger */
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    ierr = DataExView(data_exchanger);CHKERRQ(ierr);
    
    ierr = DataExDestroy(data_exchanger);CHKERRQ(ierr);
	DataBucketDestroy(&db);
    
    PetscFunctionReturn(0);
}

/* 
 Communication example when data bucket contains two different data types.
 Strategy adopted is simple, but slightly awkard.
 The strategy is as follows:
   (i)   Setup the send counts
   (ii)  Pack,send,unpack the individual fields
   (iii) This requires us to tag points which are to be removed and this
         requires we send/add received entries BEFORE removing points from the data bucket
*/
#undef __FUNCT__
#define __FUNCT__ "SwarmTest_Communication2"
PetscErrorCode SwarmTest_Communication2(MPI_Comm comm)
{
    PetscErrorCode ierr;
    DataEx          data_exchanger;
	DataBucket      db;
	DataField       dbField,dbFieldVP;
    int             nproc,rank,i,k,init_size,L_local;
    long int        L_global;
    PetscInt        num_items_recv;
    MaterialPoint   *mp,*mp_recv;
    MaterialPointVP *mpv,*mpv_recv;
    
    
    ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    
    if (nproc != 2) {
        SETERRQ(comm,PETSC_ERR_SUP,"Example only valid if nproc = 2");
    }
    
    if (rank) {
        printf("[[%s]]\n",__FUNCTION__);
    }
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    
    DataBucketCreate(&db);
	DataBucketRegisterField(db,MaterialPointClassName,  sizeof(MaterialPoint),  NULL);
	DataBucketRegisterField(db,MaterialPointVPClassName,sizeof(MaterialPointVP),NULL);
	DataBucketFinalize(db);
    
    /* Each processor will define a different number of entries */
    if (rank == 0) { init_size = 7; }
    if (rank == 1) { init_size = 3; }
	DataBucketSetSizes(db,init_size,-1);
	DataBucketGetSizes(db,&L_local,NULL,NULL);
    
	DataBucketGetDataFieldByName(db,MaterialPointClassName,  &dbField);
    DataBucketGetDataFieldByName(db,MaterialPointVPClassName,&dbFieldVP);
    DataFieldGetEntries(dbField,  (void**)&mp);
    DataFieldGetEntries(dbFieldVP,(void**)&mpv);

    for (k=0; k<L_local; k++) {
        /* Set values of standard material point data */
        mp[k].coor[0] = (double)3*k + 1.1 + (double)(rank*10);
        mp[k].coor[1] = (double)3*k + 2.2 + (double)(rank*10);
        mp[k].coor[2] = (double)3*k + 3.3 + (double)(rank*10);
        
        mp[k].pid       = (long int)k + (double)(rank*100);
        mp[k].viscosity = (float)(k+10) + (float)(rank*10);
        mp[k].region_id = (short int)(rank+1);

        /* Set values of visco-plastic material point data */
        mpv[k].plastic_strain = 0.1 * ( (float)(k+10) + (double)(rank*100) );
        mpv[k].failure_type   = (char)(rank+1);
    }

    DataFieldRestoreEntries(dbFieldVP,(void**)&mpv);
    DataFieldRestoreEntries(dbField,(void**)&mp);


    /* Examine result */
	DataBucketGetSizes(db,&L_local,NULL,NULL);

	DataBucketGetDataFieldByName(db,MaterialPointClassName,  &dbField);
    DataBucketGetDataFieldByName(db,MaterialPointVPClassName,&dbFieldVP);
    DataFieldGetEntries(dbField,  (void**)&mp);
    DataFieldGetEntries(dbFieldVP,(void**)&mpv);
    
    ierr = MPI_Barrier(comm);CHKERRQ(ierr); fflush(stdout);
    if (rank) { printf("[Initial data bucket (pre communication)]\n"); fflush(stdout); } ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    for (i=0; i<nproc; i++) {
        ierr = MPI_Barrier(comm);CHKERRQ(ierr);
        if (rank == i) {
            for (k=0; k<L_local; k++) {
                printf("  [rank %.2d] [%.2d] { pid = %.4ld ; viscosity = %1.4e ; e_plastic = %1.4e }\n",
                       rank,k,mp[k].pid,mp[k].viscosity,mpv[k].plastic_strain); fflush(stdout);
            }
        }
    } fflush(stdout);
    DataFieldRestoreEntries(dbField,(void**)&mpv);
    DataFieldRestoreEntries(dbField,(void**)&mp);
    
    /* Get total number of entries (summed over all ranks in comm */
	DataBucketGetGlobalSizes(comm,db,&L_global,NULL,NULL);
    
    /* Report parallel information about data bucket */
    DataBucketView(comm,db,"Material Point Coefficients",DATABUCKET_VIEW_STDOUT);
    
    
    
    /* Create the data exchanger */
    data_exchanger = DataExCreate(comm,0);
    
    
    /* Define communication topology */
    ierr = DataExTopologyInitialize(data_exchanger);CHKERRQ(ierr);
    
    if (rank == 0) {  ierr = DataExTopologyAddNeighbour(data_exchanger,1);CHKERRQ(ierr); }
    if (rank == 1) {  ierr = DataExTopologyAddNeighbour(data_exchanger,0);CHKERRQ(ierr); }
    
    ierr = DataExTopologyFinalize(data_exchanger);CHKERRQ(ierr);
    
    /* Phase a) Establish send counts */
    ierr = DataExInitializeSendCount(data_exchanger);CHKERRQ(ierr);
    if (rank == 0) { /* send three items form rank 0 to rank 1 */
        ierr = DataExAddToSendCount(data_exchanger,1,3);CHKERRQ(ierr);
    }
    if (rank == 1) { /* send one item from rank 1 to rank 0 */
        ierr = DataExAddToSendCount(data_exchanger,0,1);CHKERRQ(ierr);
    }
    ierr = DataExFinalizeSendCount(data_exchanger);CHKERRQ(ierr);

    /* ---------------------------------------------------------- */
    /* Phase b) Pack data into buffers to be sent [MaterialPoint] */
    ierr = DataExPackInitialize(data_exchanger,sizeof(MaterialPoint));CHKERRQ(ierr);
    
    DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
    DataFieldGetEntries(dbField,(void**)&mp);
    if (rank == 0) {
        /* Insert markers 1, 3, 6 (to be sent to rank 1) */
        ierr = DataExPackData(data_exchanger,1,1,&mp[1]);CHKERRQ(ierr);
        ierr = DataExPackData(data_exchanger,1,1,&mp[3]);CHKERRQ(ierr);
        ierr = DataExPackData(data_exchanger,1,1,&mp[6]);CHKERRQ(ierr);
        
        /* tag points to be removed */
        mp[1].region_id = -10;
        mp[3].region_id = -10;
        mp[6].region_id = -10;
    }
    if (rank == 1) {
        /* Insert markers 0 (to be sent to rank 0) */
        ierr = DataExPackData(data_exchanger,0,1,&mp[0]);CHKERRQ(ierr);

        /* tag points to be removed */
        mp[0].region_id = -10;
    }
    DataFieldRestoreEntries(dbField,(void**)&mp);
    
    /* Finalize packing - no further data can be added into the buffer after this call */
    ierr = DataExPackFinalize(data_exchanger);CHKERRQ(ierr);
    
    /* Phase c) Send the data */
    ierr = DataExBegin(data_exchanger);CHKERRQ(ierr);
    ierr = DataExEnd(data_exchanger);CHKERRQ(ierr);
    
    /* Phase d) Fetch the recv buffer and do something with it */
    ierr = DataExGetRecvData(data_exchanger,&num_items_recv,(void**)&mp_recv);CHKERRQ(ierr);

    /* Increase size of data bucket */
	DataBucketGetSizes(db,&L_local,NULL,NULL);
	DataBucketSetSizes(db,L_local + num_items_recv,-1);
    
    /* Insert new values for MaterialPoint */
    DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
    for (k=0; k<num_items_recv; k++) {
        DataFieldInsertPoint(dbField,L_local+k,&mp_recv[k]);
    }

    
    
    /* ------------------------------------------------------------ */
    /* Phase b) Pack data into buffers to be sent [MaterialPointVP] */
    ierr = DataExPackInitialize(data_exchanger,sizeof(MaterialPointVP));CHKERRQ(ierr);
    
    DataBucketGetDataFieldByName(db,MaterialPointVPClassName,&dbFieldVP);
    DataFieldGetEntries(dbFieldVP,(void**)&mpv);
    if (rank == 0) {
        /* Insert markers 1, 3, 6 (to be sent to rank 1) */
        ierr = DataExPackData(data_exchanger,1,1,&mpv[1]);CHKERRQ(ierr);
        ierr = DataExPackData(data_exchanger,1,1,&mpv[3]);CHKERRQ(ierr);
        ierr = DataExPackData(data_exchanger,1,1,&mpv[6]);CHKERRQ(ierr);
    }
    if (rank == 1) {
        /* Insert markers 0 (to be sent to rank 0) */
        ierr = DataExPackData(data_exchanger,0,1,&mpv[0]);CHKERRQ(ierr);
    }
    DataFieldRestoreEntries(dbFieldVP,(void**)&mpv);
    
    /* Finalize packing - no further data can be added into the buffer after this call */
    ierr = DataExPackFinalize(data_exchanger);CHKERRQ(ierr);
    
    /* Phase c) Send the data */
    ierr = DataExBegin(data_exchanger);CHKERRQ(ierr);
    ierr = DataExEnd(data_exchanger);CHKERRQ(ierr);
    
    /* Phase d) Fetch the recv buffer and do something with it */
    ierr = DataExGetRecvData(data_exchanger,&num_items_recv,(void**)&mpv_recv);CHKERRQ(ierr);

    /* Insert new values for MaterialPointVP */
    DataBucketGetDataFieldByName(db,MaterialPointVPClassName,&dbFieldVP);
    for (k=0; k<num_items_recv; k++) {
        DataFieldInsertPoint(dbFieldVP,L_local+k,&mpv_recv[k]);
    }
    
    
    /*
     [DATASTRUCTURE - APPLICATION SPECIFIC] NOTE
     
     Users should consider if they want to remove the items being sent from their data structure.
     This is an application specific choice.
     For a typical material point implementation, sending a marker to another processor indicates
     that the material point has physically moved from one part of the domain to another, and thus
     "ownership" of the material point has changed - hence we should remove the material points which
     are being sent.
     */
	DataBucketGetSizes(db,&L_local,NULL,NULL);
	DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
    DataFieldGetEntries(dbField,(void**)&mp);
    for (k=0; k<L_local; k++) {
        if (mp[k].region_id == -10) {
            DataBucketRemovePointAtIndex(db,k);
        }
    }
    DataFieldRestoreEntries(dbField,(void**)&mp);
    
    
    /* View result */
	DataBucketGetSizes(db,&L_local,NULL,NULL);

	DataBucketGetDataFieldByName(db,MaterialPointClassName,  &dbField);
    DataBucketGetDataFieldByName(db,MaterialPointVPClassName,&dbFieldVP);
    DataFieldGetEntries(dbField,  (void**)&mp);
    DataFieldGetEntries(dbFieldVP,(void**)&mpv);
    ierr = MPI_Barrier(comm);CHKERRQ(ierr); fflush(stdout);
    if (rank) { printf("[Final data bucket (post communication)]\n"); fflush(stdout); } ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    for (i=0; i<nproc; i++) {
        ierr = MPI_Barrier(comm);CHKERRQ(ierr);
        if (rank == i) {
            for (k=0; k<L_local; k++) {
                printf("  [*rank %.2d] [%.2d] { pid = %.4ld ; viscosity = %1.4e ; e_plastic = %1.4e }\n",
                       rank,k,mp[k].pid,mp[k].viscosity,mpv[k].plastic_strain); fflush(stdout);
            }
        }
    } fflush(stdout);
    DataFieldRestoreEntries(dbFieldVP,(void**)&mpv);
    DataFieldRestoreEntries(dbField,(void**)&mp);
    
    /* View information about data_exchanger */
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    ierr = DataExView(data_exchanger);CHKERRQ(ierr);
    
    ierr = DataExDestroy(data_exchanger);CHKERRQ(ierr);
	DataBucketDestroy(&db);
    
    PetscFunctionReturn(0);
}

/*
 Communication example when data bucket contains two different data types.
 Strategy adopted is packs/sends/unpacks all entries in the data bucket at once.
 */
#undef __FUNCT__
#define __FUNCT__ "SwarmTest_Communication3"
PetscErrorCode SwarmTest_Communication3(MPI_Comm comm)
{
    PetscErrorCode ierr;
    DataEx          data_exchanger;
	DataBucket      db;
	DataField       dbField,dbFieldVP;/*,*dbAllFields;*/
    int             nproc,rank,i,k,init_size,L_local;
    long int        L_global;
    PetscInt        num_items_recv;
    MaterialPoint   *mp;
    MaterialPointVP *mpv;
/*    int             nfields; */
	size_t          sizeof_marker_contents;
	void            /**data_buffer,*/*recv_buffer,*dbuf;
    
    
    ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    
    if (nproc != 2) {
        SETERRQ(comm,PETSC_ERR_SUP,"Example only valid if nproc = 2");
    }
    
    if (rank) {
        printf("[[%s]]\n",__FUNCTION__);
    }
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    
    DataBucketCreate(&db);
	DataBucketRegisterField(db,MaterialPointClassName,  sizeof(MaterialPoint),  NULL);
	DataBucketRegisterField(db,MaterialPointVPClassName,sizeof(MaterialPointVP),NULL);
	DataBucketFinalize(db);
    
    /* Each processor will define a different number of entries */
    if (rank == 0) { init_size = 7; }
    if (rank == 1) { init_size = 3; }
	DataBucketSetSizes(db,init_size,-1);
	DataBucketGetSizes(db,&L_local,NULL,NULL);
    
	DataBucketGetDataFieldByName(db,MaterialPointClassName,  &dbField);
    DataBucketGetDataFieldByName(db,MaterialPointVPClassName,&dbFieldVP);
    DataFieldGetEntries(dbField,  (void**)&mp);
    DataFieldGetEntries(dbFieldVP,(void**)&mpv);
    
    for (k=0; k<L_local; k++) {
        /* Set values of standard material point data */
        mp[k].coor[0] = (double)3*k + 1.1 + (double)(rank*10);
        mp[k].coor[1] = (double)3*k + 2.2 + (double)(rank*10);
        mp[k].coor[2] = (double)3*k + 3.3 + (double)(rank*10);
        
        mp[k].pid       = (long int)k + (double)(rank*100);
        mp[k].viscosity = (float)(k+10) + (float)(rank*10);
        mp[k].region_id = (short int)(rank+1);
        
        /* Set values of visco-plastic material point data */
        mpv[k].plastic_strain = 0.1 * ( (float)(k+10) + (double)(rank*100) );
        mpv[k].failure_type   = (char)(rank+1);
    }
    
    DataFieldRestoreEntries(dbFieldVP,(void**)&mpv);
    DataFieldRestoreEntries(dbField,(void**)&mp);
    
    
    /* Examine result */
	DataBucketGetSizes(db,&L_local,NULL,NULL);
    
	DataBucketGetDataFieldByName(db,MaterialPointClassName,  &dbField);
    DataBucketGetDataFieldByName(db,MaterialPointVPClassName,&dbFieldVP);
    DataFieldGetEntries(dbField,  (void**)&mp);
    DataFieldGetEntries(dbFieldVP,(void**)&mpv);
    
    ierr = MPI_Barrier(comm);CHKERRQ(ierr); fflush(stdout);
    if (rank) { printf("[Initial data bucket (pre communication)]\n"); fflush(stdout); } ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    for (i=0; i<nproc; i++) {
        ierr = MPI_Barrier(comm);CHKERRQ(ierr);
        if (rank == i) {
            for (k=0; k<L_local; k++) {
                printf("  [rank %.2d] [%.2d] { pid = %.4ld ; viscosity = %1.4e ; e_plastic = %1.4e }\n",
                       rank,k,mp[k].pid,mp[k].viscosity,mpv[k].plastic_strain); fflush(stdout);
            }
        }
    } fflush(stdout);
    DataFieldRestoreEntries(dbField,(void**)&mpv);
    DataFieldRestoreEntries(dbField,(void**)&mp);
    
    /* Get total number of entries (summed over all ranks in comm */
	DataBucketGetGlobalSizes(comm,db,&L_global,NULL,NULL);
    
    /* Report parallel information about data bucket */
    DataBucketView(comm,db,"Material Point Coefficients",DATABUCKET_VIEW_STDOUT);
    
    /* Create the data exchanger */
    data_exchanger = DataExCreate(comm,0);
    
    
    /* Define communication topology */
    ierr = DataExTopologyInitialize(data_exchanger);CHKERRQ(ierr);
    
    if (rank == 0) {  ierr = DataExTopologyAddNeighbour(data_exchanger,1);CHKERRQ(ierr); }
    if (rank == 1) {  ierr = DataExTopologyAddNeighbour(data_exchanger,0);CHKERRQ(ierr); }
    
    ierr = DataExTopologyFinalize(data_exchanger);CHKERRQ(ierr);
    
    /* Phase a) Establish send counts */
    ierr = DataExInitializeSendCount(data_exchanger);CHKERRQ(ierr);
    if (rank == 0) { /* send three items form rank 0 to rank 1 */
        ierr = DataExAddToSendCount(data_exchanger,1,3);CHKERRQ(ierr);
    }
    if (rank == 1) { /* send one item from rank 1 to rank 0 */
        ierr = DataExAddToSendCount(data_exchanger,0,1);CHKERRQ(ierr);
    }
    ierr = DataExFinalizeSendCount(data_exchanger);CHKERRQ(ierr);
    
    /* ---------------------------------------------------------- */
    /* compute size of merged data bucket entries */
    /*
    DataBucketGetDataFields(db,&nfields,&dbAllFields);
    sizeof_marker_contents = 0;
    for (f=0; f<nfields; f++) {
        size_t asize;
        
        DataFieldGetAtomicSize(dbAllFields[f],&asize);
        sizeof_marker_contents += asize;
    }
    */

	/* allocate a temporary buffer whihc is large enough to store all marker fields from an individual point,p */
	//ierr = PetscMalloc(sizeof_marker_contents,&data_buffer);CHKERRQ(ierr);
    DataBucketCreatePackedArray(db,&sizeof_marker_contents,&dbuf);
    
    /* Phase b) Pack data into buffers to be sent [MaterialPoint + MaterialPointVP] */
    ierr = DataExPackInitialize(data_exchanger,sizeof_marker_contents);CHKERRQ(ierr);
    
    {
        int nps,pidx[3];
        
        
        if (rank == 0) {
            nps = 3;
            pidx[0] = 1;
            pidx[1] = 3;
            pidx[2] = 6;
            
        }
        if (rank == 1) {
            nps = 1;
            pidx[0] = 0;
        }
     
/*
        for (k=0; k<nps; k++) {
            void *data0,*data;
            size_t asize,offset;
            
            DataFieldGetAtomicSize(dbAllFields[0],&asize);
            DataFieldGetEntries(dbAllFields[0],&data0);
			PetscMemcpy((void*)data_buffer,    data0 + pidx[k]*asize,    asize);

            offset = asize;
            for (f=1; f<nfields; f++) {
                DataFieldGetAtomicSize(dbAllFields[f],&asize);
                DataFieldGetEntries(dbAllFields[f],&data);
                
                PetscMemcpy((void*)((char*)data_buffer + offset),    data + pidx[k]*asize,    asize);
				offset = offset + asize;
            }
            
            if (rank == 0) {
                ierr = DataExPackData(data_exchanger,1,1,(void*)data_buffer);CHKERRQ(ierr);
            }
            if (rank == 1) {
                ierr = DataExPackData(data_exchanger,0,1,(void*)data_buffer);CHKERRQ(ierr);
            }
        }
*/
  
        for (k=0; k<nps; k++) {
            DataBucketFillPackedArray(db,pidx[k],dbuf);

            if (rank == 0) {
                ierr = DataExPackData(data_exchanger,1,1,(void*)dbuf);CHKERRQ(ierr);
            }
            if (rank == 1) {
                ierr = DataExPackData(data_exchanger,0,1,(void*)dbuf);CHKERRQ(ierr);
            }
        }
    }
    
    /* Finalize packing - no further data can be added into the buffer after this call */
    ierr = DataExPackFinalize(data_exchanger);CHKERRQ(ierr);
    
    
    /*
     [DATASTRUCTURE - APPLICATION SPECIFIC] NOTE
     Rempove points who were sent
    */
    if (rank == 0) {
        /* I remove them in reverse order for convienence */
        DataBucketRemovePointAtIndex(db,6);
        DataBucketRemovePointAtIndex(db,3);
        DataBucketRemovePointAtIndex(db,1);
    }
    if (rank == 1) {
        DataBucketRemovePointAtIndex(db,0);
    }

    
    /* Phase c) Send the data */
    ierr = DataExBegin(data_exchanger);CHKERRQ(ierr);
    ierr = DataExEnd(data_exchanger);CHKERRQ(ierr);
    
    /* Phase d) Fetch the recv buffer and do something with it */
    ierr = DataExGetRecvData(data_exchanger,&num_items_recv,(void**)&recv_buffer);CHKERRQ(ierr);
    
    /* Increase size of data bucket */
	DataBucketGetSizes(db,&L_local,NULL,NULL);
	DataBucketSetSizes(db,L_local + num_items_recv,-1);
    
    /* Insert new values for MaterialPoint */
    /*
    for (k=0; k<num_items_recv; k++) {
        void *data_p;
        size_t asize,offset;
        
        offset = 0;
        for (f=0; f<nfields; f++) {
            data_p = (void*)( (char*)recv_buffer + k*(sizeof_marker_contents) + offset );
            
            DataFieldInsertPoint(dbAllFields[f], L_local + k, (void*)data_p );
            
            DataFieldGetAtomicSize(dbAllFields[f],&asize);
            offset = offset + asize;
        }
    }
    */
    for (k=0; k<num_items_recv; k++) {
        void *data_p = (void*)( (char*)recv_buffer + k*(sizeof_marker_contents) );
        
        DataBucketInsertPackedArray(db,L_local+k,data_p);
    }
    
    /* View result */
	DataBucketGetSizes(db,&L_local,NULL,NULL);
    
	DataBucketGetDataFieldByName(db,MaterialPointClassName,  &dbField);
    DataBucketGetDataFieldByName(db,MaterialPointVPClassName,&dbFieldVP);
    DataFieldGetEntries(dbField,  (void**)&mp);
    DataFieldGetEntries(dbFieldVP,(void**)&mpv);
    ierr = MPI_Barrier(comm);CHKERRQ(ierr); fflush(stdout);
    if (rank) { printf("[Final data bucket (post communication)]\n"); fflush(stdout); } ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    for (i=0; i<nproc; i++) {
        ierr = MPI_Barrier(comm);CHKERRQ(ierr);
        if (rank == i) {
            for (k=0; k<L_local; k++) {
                printf("  [*rank %.2d] [%.2d] { pid = %.4ld ; viscosity = %1.4e ; e_plastic = %1.4e }\n",
                       rank,k,mp[k].pid,mp[k].viscosity,mpv[k].plastic_strain); fflush(stdout);
            }
        }
    } fflush(stdout);
    DataFieldRestoreEntries(dbFieldVP,(void**)&mpv);
    DataFieldRestoreEntries(dbField,(void**)&mp);
    
    /* View information about data_exchanger */
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    ierr = DataExView(data_exchanger);CHKERRQ(ierr);
    
    DataBucketDestroyPackedArray(db,&dbuf);
	//ierr = PetscFree(data_buffer);CHKERRQ(ierr);
    ierr = DataExDestroy(data_exchanger);CHKERRQ(ierr);
	DataBucketDestroy(&db);
    
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int nargs,char **args)
{
    int test_id;
    
    PetscErrorCode ierr;
    ierr = PetscInitialize(&nargs,&args,(char*)0,NULL);CHKERRQ(ierr);
    
    if (nargs == 2) {
        test_id = atoi(args[1]);
    } else {
        test_id = 1;
    }
    
    /* serial tests */
    switch (test_id) {
            
        case 1:
            ierr = SwarmTest_Communication1(PETSC_COMM_WORLD);CHKERRQ(ierr);
            break;
        case 2:
            ierr = SwarmTest_Communication2(PETSC_COMM_WORLD);CHKERRQ(ierr);
            break;
        case 3:
            ierr = SwarmTest_Communication3(PETSC_COMM_WORLD);CHKERRQ(ierr);
            break;
            
    }

    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}

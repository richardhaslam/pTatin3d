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
 **    filename:   test_swarms.c
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

/*
   Information and usage:

   To compile as a stand alone example, you need to
   0) Install MPI
   1) Create the object file for data_bucket.c
   2) Compile test_swarms.c and link against data_bucket.o

   [step 1]
   ${MPI_DIR}/bin/mpicc -O0 -g -c data_bucket.c -I. -I${MPI_DIR}/include
   [step 2]
   ${MPI_DIR}/bin/mpicc -O0 -g -o test_swarms test_swarms.c data_bucket.o -I. -I${MPI_DIR}/include -L${MPI_DIR}/lib -lmpi

   /Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/bin/mpicc -O0 -g -c data_bucket.c -I. -I/Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/include
   /Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/bin/mpicc -O0 -g -o test_swarms test_swarms.c data_bucket.o -I. -I/Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/include -L/Users/dmay/software/petsc-3.5.1/arch-darwin-c-mpi-debug/lib -lmpi
   */


#include <data_bucket.h>

/* Define structure for storing material point data */
typedef struct {
  double    coor[3];
  long int  pid;
  float     viscosity;
  short int region_id;
} MaterialPoint;

/* Define name used to identify data structure for material point data */
const char MaterialPointClassName[] = "DBF_MaterialPoint";



int SwarmTest_Initialization1(void)
{
  DataBucket db;
  int        L;


  printf("[[%s]]\n",__FUNCTION__);

  /* Create object */
  DataBucketCreate(&db);

  /* Load a specific data type into data bucket, specifying textual name and its size in bytes */
  DataBucketRegisterField(db,MaterialPointClassName,sizeof(MaterialPoint),NULL);

  DataBucketFinalize(db);

  /* Set number of instances of all data within the data bucket (active size, total size) */
  DataBucketSetSizes(db,10,-1);

  /* Get number of entries */
  DataBucketGetSizes(db,&L,NULL,NULL);

  /* Report information about the bucket and its contents */
  DataBucketView(MPI_COMM_SELF,db,"Material Point Coefficients",DATABUCKET_VIEW_STDOUT);

  /* Destroy object */
  DataBucketDestroy(&db);

  return(0);
}

int SwarmTest_Initialization2(void)
{
  DataBucket db;
  int        L;


  printf("[[%s]]\n",__FUNCTION__);

  /* Create object */
  DataBucketCreate(&db);

  /* Load a specific data type into data bucket, specifying textual name and its size in bytes */
  DataBucketRegisterField(db,"MP_vx",sizeof(double),NULL);
  DataBucketRegisterField(db,"MP_vy",sizeof(double),NULL);
  DataBucketRegisterField(db,"MP_pid",sizeof(short int),NULL);

  DataBucketFinalize(db);

  /* Set number of instances of all data within the data bucket (active size, total size) */
  DataBucketSetSizes(db,10,-1);

  /* Get number of entries */
  DataBucketGetSizes(db,&L,NULL,NULL);

  /* Report information about the bucket and its contents */
  DataBucketView(MPI_COMM_SELF,db,"Material Point Data Arrays",DATABUCKET_VIEW_STDOUT);

  /* Destroy object */
  DataBucketDestroy(&db);

  return(0);
}

int SwarmTest_AccessPatterns1(void)
{
  DataBucket db;
  DataField  dbField;
  int        k,L;


  printf("[[%s]]\n",__FUNCTION__);

  DataBucketCreate(&db);
  DataBucketRegisterField(db,MaterialPointClassName,sizeof(MaterialPoint),NULL);
  DataBucketFinalize(db);
  DataBucketSetSizes(db,10,-1);

  /* Get number of entries */
  DataBucketGetSizes(db,&L,NULL,NULL);

  /* Set data into the data bucket */
  /*   a) Fetch data type from bucket using textual name */
  DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
  /*   b) Get access to data type. Calls to DataFieldGetAccess() must be paired with a call to DataFieldRestoreAccess() */
  DataFieldGetAccess(dbField);
  /*   c) Optional: Perform check to ensure data field requested matches size of what you expect */
  DataFieldVerifyAccess(dbField,sizeof(MaterialPoint));

  for (k=0; k<L; k++) {
    MaterialPoint *mp_k;

    /* From the data field object, get pointer to the k'th entry of type MaterialPoint */
    DataFieldAccessPoint(dbField,k,(void**)&mp_k);

    /* Set values in the k'th entry */
    mp_k->coor[0] = (double)3*k + 1.1;
    mp_k->coor[1] = (double)3*k + 2.2;
    mp_k->coor[2] = (double)3*k + 3.3;

    mp_k->pid       = (long int)k;
    mp_k->viscosity = (float)(k*10);
    mp_k->region_id = 2;
  }

  /*   d) Restore access to data type */
  DataFieldRestoreAccess(dbField);

  /* Examine result */
  DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
  DataFieldGetAccess(dbField);
  for (k=0; k<L; k++) {
    MaterialPoint *mp_k;

    DataFieldAccessPoint(dbField,k,(void**)&mp_k);
    printf("  [%.2d] pid = %ld ; viscosity = %1.4e\n",k,mp_k->pid,mp_k->viscosity);
  }
  DataFieldRestoreAccess(dbField);

  DataBucketDestroy(&db);

  return(0);
}

int SwarmTest_AccessPatterns2(void)
{
  DataBucket    db;
  DataField     dbField;
  MaterialPoint *mp_list;
  int           k,L;


  printf("[[%s]]\n",__FUNCTION__);

  DataBucketCreate(&db);
  DataBucketRegisterField(db,MaterialPointClassName,sizeof(MaterialPoint),NULL);
  DataBucketFinalize(db);
  DataBucketSetSizes(db,10,-1);

  /* Get number of entries */
  DataBucketGetSizes(db,&L,NULL,NULL);

  /* Set data into the data bucket */
  /*   a) Fetch data type from bucket using textual name */
  DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
  /*   b) Get access to data type. Calls to DataFieldGetAccess() must be paired with a call to DataFieldRestoreAccess() */
  DataFieldGetAccess(dbField);
  /*   c) Optional: Perform check to ensure data field requested matches size of what you expect */
  DataFieldVerifyAccess(dbField,sizeof(MaterialPoint));

  for (k=0; k<L; k++) {
    MaterialPoint *mp_k;

    /* From the data field object, get pointer to the k'th entry of type MaterialPoint */
    DataFieldAccessPoint(dbField,k,(void**)&mp_k);

    /* Set values in the k'th entry */
    mp_k->coor[0] = (double)3*k + 1.1;
    mp_k->coor[1] = (double)3*k + 2.2;
    mp_k->coor[2] = (double)3*k + 3.3;

    mp_k->pid       = (long int)k;
    mp_k->viscosity = (float)(k*10);
    mp_k->region_id = 2;
  }

  /*   d) Restore access to data type */
  DataFieldRestoreAccess(dbField);

  /* Examine result */
  DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);

  /*
     Get access to complete array of all entries (rather than accessing them one by one).
     Less safe, recommended if you are sure you know what you are doing.
     Calls to DataFieldGetEntries() must be paired with calls to DataFieldRestoreEntries().
     */
  DataFieldGetEntries(dbField,(void**)&mp_list);
  for (k=0; k<L; k++) {
    printf("  [%.2d] pid = %ld ; viscosity = %1.4e\n",k,mp_list[k].pid,mp_list[k].viscosity);
  }
  DataFieldRestoreEntries(dbField,(void**)&mp_list);

  DataBucketDestroy(&db);

  return(0);
}

int SwarmTest_AccessPatterns3(void)
{
  DataBucket db;
  DataField  dbField;
  double     *data_d;
  short int  *data_si;
  int        k,L;


  printf("[[%s]]\n",__FUNCTION__);

  /* Create object */
  DataBucketCreate(&db);

  /* Load a specific data type into data bucket, specifying textual name and its size in bytes */
  DataBucketRegisterField(db,"MP_vx",sizeof(double),NULL);
  DataBucketRegisterField(db,"MP_vy",sizeof(double),NULL);
  DataBucketRegisterField(db,"MP_pid",sizeof(short int),NULL);

  DataBucketFinalize(db);

  DataBucketSetSizes(db,5,-1);

  DataBucketGetSizes(db,&L,NULL,NULL);

  /* Set array entries : vx */
  DataBucketGetDataFieldByName(db,"MP_vx",&dbField);
  DataFieldGetEntries(dbField,(void**)&data_d);
  for (k=0; k<L; k++) {
    data_d[k] = k + 1.1;
  }
  DataFieldRestoreEntries(dbField,(void**)&data_d);

  /* Set array entries : vy */
  DataBucketGetDataFieldByName(db,"MP_vy",&dbField);
  DataFieldGetEntries(dbField,(void**)&data_d);
  for (k=0; k<L; k++) {
    data_d[k] = k + 12.1;
  }
  DataFieldRestoreEntries(dbField,(void**)&data_d);

  /* Set array entries : pid */
  DataBucketGetDataFieldByName(db,"MP_pid",&dbField);
  DataFieldGetEntries(dbField,(void**)&data_si);
  for (k=0; k<L; k++) {
    data_si[k] = k + 1;
  }
  DataFieldRestoreEntries(dbField,(void**)&data_si);

  DataBucketGetDataFieldByName(db,"MP_vx",&dbField);
  DataFieldGetEntries(dbField,(void**)&data_d);
  printf("vx = {\n");
  for (k=0; k<L; k++) {
    printf("  %1.4e\n",data_d[k]);
  } printf("}\n");
  DataFieldRestoreEntries(dbField,(void**)&data_d);

  DataBucketGetDataFieldByName(db,"MP_vy",&dbField);
  DataFieldGetEntries(dbField,(void**)&data_d);
  printf("vy = {\n");
  for (k=0; k<L; k++) {
    printf("  %1.4e\n",data_d[k]);
  } printf("}\n");
  DataFieldRestoreEntries(dbField,(void**)&data_d);

  DataBucketGetDataFieldByName(db,"MP_pid",&dbField);
  DataFieldGetEntries(dbField,(void**)&data_si);
  printf("pid = {\n");
  for (k=0; k<L; k++) {
    printf("  %d\n",data_si[k]);
  } printf("}\n");
  DataFieldRestoreEntries(dbField,(void**)&data_si);

  DataBucketDestroy(&db);

  return(0);
}

int SwarmTest_LengthManipulations1(void)
{
  DataBucket    db;
  DataField     dbField;
  MaterialPoint *mp_list;
  int           k,L,Lnew;


  printf("[[%s]]\n",__FUNCTION__);

  DataBucketCreate(&db);
  DataBucketRegisterField(db,MaterialPointClassName,sizeof(MaterialPoint),NULL);
  DataBucketFinalize(db);
  DataBucketSetSizes(db,6,-1);

  DataBucketGetSizes(db,&L,NULL,NULL);

  DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);

  /* Set data into the data bucket */
  DataFieldGetAccess(dbField);
  DataBucketGetSizes(db,&L,NULL,NULL);
  for (k=0; k<L; k++) {
    MaterialPoint *mp_k;

    /* From the data field object, get pointer to the k'th entry of type MaterialPoint */
    DataFieldAccessPoint(dbField,k,(void**)&mp_k);

    /* Set values in the k'th entry */
    mp_k->coor[0] = (double)3*k + 1.1;
    mp_k->coor[1] = (double)3*k + 2.2;
    mp_k->coor[2] = (double)3*k + 3.3;

    mp_k->pid       = (long int)k;
    mp_k->viscosity = (float)(k*10);
    mp_k->region_id = 2;
  }
  DataFieldRestoreAccess(dbField);

  /* Examine result */
  DataFieldGetEntries(dbField,(void**)&mp_list);
  DataBucketGetSizes(db,&Lnew,NULL,NULL);
  printf("[a] Initial entries\n");
  for (k=0; k<L; k++) {
    printf("  [%.2d] pid = %ld ; viscosity = %1.4e\n",k,mp_list[k].pid,mp_list[k].viscosity);
  }
  DataFieldRestoreEntries(dbField,(void**)&mp_list);

  /*
     Resize list - extend

     void DataBucketSetSizes(DataBucket db,int new_size,int buffer)
     int new_size : Set a new size of the list
     int buffer : Indicate a new buffer size - or keep old buffer size by passing in -1

Note: New entries are placed at the end of the list
*/
  DataBucketSetSizes(db,10,20);
  DataBucketGetSizes(db,&Lnew,NULL,NULL);

  /* Set values on new points */
  DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
  DataFieldGetEntries(dbField,(void**)&mp_list);
  for (k=L; k<Lnew; k++) {
    mp_list[k].coor[0] = (double)3*k + 10.1;
    mp_list[k].coor[1] = (double)3*k + 20.2;
    mp_list[k].coor[2] = (double)3*k + 30.3;

    mp_list[k].pid       = (long int)k;
    mp_list[k].viscosity = (float)(k*100);
    mp_list[k].region_id = 3;
  }
  DataFieldRestoreEntries(dbField,(void**)&mp_list);

  /* Examine result */
  DataFieldGetEntries(dbField,(void**)&mp_list);
  DataBucketGetSizes(db,&Lnew,NULL,NULL);
  printf("[b] Current entries (extended)\n");
  for (k=0; k<Lnew; k++) {
    printf("  [%.2d] pid = %ld ; viscosity = %1.4e\n",k,mp_list[k].pid,mp_list[k].viscosity);
  }
  DataFieldRestoreEntries(dbField,(void**)&mp_list);

  /*
     Resize list - remove the k'th point

     void DataBucketRemovePointAtIndex(DataBucket db,int idx)
     int idx : Index of entry which will be removed

Note: When an entry is removed, it is replaced by the entry at the end of the array
and the total size of the list is decreased by one
*/
  DataBucketRemovePointAtIndex(db,2);
  DataBucketRemovePointAtIndex(db,4);
  DataBucketRemovePointAtIndex(db,7);

  DataBucketGetDataFieldByName(db,MaterialPointClassName,&dbField);
  DataFieldGetEntries(dbField,(void**)&mp_list);
  DataBucketGetSizes(db,&Lnew,NULL,NULL);
  printf("[c] Current entries (2,4,8 removed)\n");
  for (k=0; k<Lnew; k++) {
    printf("  [%.2d] pid = %ld ; viscosity = %1.4e\n",k,mp_list[k].pid,mp_list[k].viscosity);
  }
  DataFieldRestoreEntries(dbField,(void**)&mp_list);

  /* Report information about the bucket and its contents */
  DataBucketView(MPI_COMM_SELF,db,"Material Point Data Arrays",DATABUCKET_VIEW_STDOUT);

  DataBucketDestroy(&db);

  return(0);
}

int SwarmTest_Parallel1(MPI_Comm comm)
{
  DataBucket    db;
  DataField     dbField;
  MaterialPoint *mp_list;
  int           nproc,rank,i,k,init_size,L_local;
  long int      L_global;


  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);
  if (rank) {
    printf("[[%s]]\n",__FUNCTION__);
  }
  MPI_Barrier(comm);

  DataBucketCreate(&db);
  DataBucketRegisterField(db,MaterialPointClassName,sizeof(MaterialPoint),NULL);
  DataBucketFinalize(db);

  /* Each processor will define a different number of entries */
  init_size = (rank + 1) * 3;
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
  for (i=0; i<nproc; i++) {
    MPI_Barrier(comm);
    if (rank == i) {
      for (k=0; k<L_local; k++) {
        printf("  [rank %.2d] [%.2d] pid = %ld ; viscosity = %1.4e\n",rank,k,mp_list[k].pid,mp_list[k].viscosity);
      }
    }
  }
  DataFieldRestoreEntries(dbField,(void**)&mp_list);

  /* Get total number of entries (summed over all ranks in comm */
  DataBucketGetGlobalSizes(comm,db,&L_global,NULL,NULL);

  /* Report parallel information about data bucket */
  DataBucketView(comm,db,"Material Point Coefficients",DATABUCKET_VIEW_STDOUT);

  DataBucketDestroy(&db);

  return(0);
}

int main(int nargs,char **args)
{
  int test_id;


  MPI_Init(&nargs,&args);

  if (nargs == 2) {
    test_id = atoi(args[1]);
  } else {
    test_id = 1;
  }

  switch (test_id) {

    /* serial tests */

    /* Demonstrate how to define a data bucket */
    case 1:
      SwarmTest_Initialization1();
      break;
    case 2:
      SwarmTest_Initialization2();
      break;

      /* Demonstrate how to access entries contained within a data bucket */
    case 10:
      SwarmTest_AccessPatterns1();
      break;
    case 11:
      SwarmTest_AccessPatterns2();
      break;
    case 12:
      SwarmTest_AccessPatterns3();
      break;

      /* Demonstrate how to manipulate the number of entries within a data bucket */
    case 20:
      SwarmTest_LengthManipulations1();
      break;

      /* parallel tests */
    case 30:
      SwarmTest_Parallel1(MPI_COMM_WORLD);
      break;
  }

  MPI_Finalize();

  return 0;
}

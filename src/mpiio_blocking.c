
#include <stdio.h>
#include <mpi.h>
#include <petsc.h>
#include "mpiio_blocking.h"

/*

 MPIWrite_Blocking()
 
 Description:
   Blocking MPI output to disk implementation in which data is streamed to root, 
   and only root writes to disk.
   Each MPI-rank sends its data to root, where it is written to disk.
   The data written to disk will be in order of the rank which defined it.

 Collective over comm

 Input:
   fp   - file pointer
   data - the data (defined on each rank) to be written into the file referenced by fp
   len  - the number of entries in data
   size - the size in bytes of each entry in data[]
   root - the rank responsible for collecting data and writing data
   comm - the MPI communicator
 
 Notes:
   * The file must be opened prior to calling this function
   * The value of fp is only required to be a valid file pointer on root
   * The values data[] is not altered
   * The implementation does not require O(p) storage
 
*/
PetscErrorCode MPIWrite_Blocking(FILE *fp,void *data,long int len,size_t size,int root,PetscBool skip_header,MPI_Comm comm)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,commsize,r;
  long int       len_max,len_r,len_total,len_total_bytes,bytes;
  void           *rbuffer;
  int            tagI,tagD;
  MPI_Request    request,requestS;
  MPI_Status     status;
  

  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  if (rank == root) {
    if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"[MPIWrite_Blocking] File error: File pointer is NULL on root");
  }

  ierr = MPI_Reduce(&len,&len_total,1,MPI_LONG,MPI_SUM,root,comm);CHKERRQ(ierr);
  len_total_bytes = len_total * size;
  ierr = MPI_Reduce(&len,&len_max,1,MPI_LONG,MPI_MAX,root,comm);CHKERRQ(ierr);

  rbuffer = NULL;
  if (rank == root) {
    rbuffer = malloc(size*len_max);
    if (!rbuffer) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"[MPIWrite_Blocking] Memory allocation error: could not allocate auxiliary buffer");
  }
  

  if (rank == root) {
    //if (rank == root) printf("<win 1> filepos %ld \n",ftell(fp));

    /* write total length into file */
    if (!skip_header) {
      fwrite(&len_total_bytes,sizeof(long int),1,fp);
    }
    /*printf("[write] total bytes %ld \n",len_total_bytes);*/

    bytes = 0;
    for (r=0; r<commsize; r++) {
      tagI = 2*r + 0;
      tagD = 2*r + 1;

      //printf("r[%d] byte write %ld \n",r,bytes);
      if (r == root) {
        /* write data */
        len_r = len;
        //if (rank == root) printf("<win> filepos %ld \n",ftell(fp));
        fwrite(data,size,len_r,fp);
      } else {
        /* flush old data */
        memset(rbuffer,0,size*len_max);
        
        /* recv length */
        ierr = MPI_Irecv(&len_r,1,MPI_LONG,r,tagI,comm,&request);CHKERRQ(ierr);
        ierr = MPI_Wait(&request,&status);CHKERRQ(ierr);

        /* recv data */
        //printf("size_t = %1.4e\n",1.0e-9*((double)(len_r*size)));
        ierr = MPI_Irecv(rbuffer,len_r*size,MPI_CHAR,r,tagD,comm,&request);CHKERRQ(ierr);
        ierr = MPI_Wait(&request,&status);CHKERRQ(ierr);

        //for (i=0; i<len_r; i++) {
        //  printf("%d %1.4e\n",i,((double*)rbuffer)[i]);
        //}

        /* write data */
        //if (rank == root) printf("<win> filepos %ld \n",ftell(fp));
        fwrite(rbuffer,size,len_r,fp);
      }
      bytes += len_r * size;
    }
  } else {
    tagI = 2*rank + 0;
    tagD = 2*rank + 1;
    /* send length then data */
    ierr = MPI_Isend(&len,1,MPI_LONG,root,tagI,comm,&requestS);CHKERRQ(ierr);
    ierr = MPI_Isend(data,len*size,MPI_CHAR,root,tagD,comm,&requestS);CHKERRQ(ierr);
  }

  ierr = MPI_Barrier(comm);CHKERRQ(ierr);
  if (rank == root) { free(rbuffer); }
  
  PetscFunctionReturn(0);
}

/*
 
 MPIRead_Blocking()
 
 Description:
   Blocking MPI input implementation in which data is streamed from root to all ranks in comm.
   Only root reads data from disk.
   Each MPI-rank will receive a chunk of data from root.
 
 Collective over comm
 
 Input:
   fp   - file pointer
   len  - the number of entries to be read from the file
   size - the size in bytes of each entry
   root - the rank responsible for reading data from disk and broadcasting the data
   comm - the MPI communicator
   data - where the bytes read from disk will be stored

 Notes:
   * The file must be opened prior to calling this function
   * The value of fp is only required to be a valid file pointer on root
   * If a valid pointer for data is provided, MPIRead_Blocking() will use the existing memory
     space. Otherwise, a new allocation will be performed
   * If sum_{all ranks} len does not equal the number of bytes written to disk, no data will be 
     read from the file (or scattered)
   * The implementation does not require O(p) storage

*/
PetscErrorCode MPIRead_Blocking(FILE *fp,void **data,long int len,size_t size,int root,PetscBool skip_header,MPI_Comm comm)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,commsize,r;
  long int       len_total,len_total_bytes,len_total_bytes_est;
  void           *sbuffer,*buffer = NULL;
  int            tagI,tagD;
  MPI_Request    request,requestS;
  MPI_Status     status,statusS;
  long int       ipack[2],ipackr[2];
  long           byte_offset;

  
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  if (!data) {
    if (!rank) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"[MPIRead_Blocking] Pointer error: A valid pointer to store data must be provided");
  }
  
  /* check file pointer is valid */
  if (rank == root) {
    if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"[MPIRead_Blocking] File error: File pointer is NULL on root");
  }

  //if (rank == root) printf("<in> filepos %ld \n",ftell(fp));
  
  ipack[0] = 0;
  ipack[1] = len;

  /* check input sums to num bytes in file */
  ierr = MPI_Allreduce(&len,&len_total,1,MPI_LONG,MPI_SUM,comm);CHKERRQ(ierr);
  len_total_bytes_est = len_total * size;

  len_total_bytes = 0;
  if (rank == root) {
    if (!skip_header) {
      /* read header */
      {
        size_t nread;
        nread =  fread(&len_total_bytes,sizeof(long int),1,fp);
        if (nread != 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_READ,"fread() failure.");
      }
      /*printf("[read] total bytes = %ld\n",len_total_bytes);*/
    }
    //printf("[read] total len = %ld\n",len_total_bytes/size);
  }

  if (!skip_header) {
    ierr = MPI_Bcast(&len_total_bytes,1,MPI_LONG,root,comm);CHKERRQ(ierr);
    if (len_total_bytes != len_total_bytes_est) {
      if (!rank) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_USER,"[MPIRead_Blocking] File error: Sizes don't match. File contains %ld bytes, your lengths sum to %ld bytes",len_total_bytes,len_total_bytes_est);
    }
  }

  /*
   User pointer is not null - assume it is long enough and re-use it, 
   otherwise allocate additional memory
  */
  if (*data) {
    buffer = *data;
  } else {
    buffer = malloc(size*len);
  }
  
  /*
  buffer = malloc(size*len);
  memset(buffer,0,size*len);
  */
  if (rank == root) {

    byte_offset = 0;
    for (r=0; r<commsize; r++) {
      sbuffer = NULL;
      
      //printf("r[%d]: %ld byte offset \n",r,byte_offset);

      if (r == root) {
        ipackr[0] = ipack[0];
        ipackr[1] = ipack[1];
        //printf("r[%d] requested [%ld,%ld] \n",r,ipackr[0],ipackr[0]+ipackr[1]);
        {
          size_t nread;
          nread = fread(buffer,size,ipackr[1],fp);
          if (nread != (size_t) ipackr[1]) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_READ,"fread() failure.");
        }
      } else {
        /* get size */
        tagI = 2*r;
        ierr = MPI_Irecv(ipackr,2,MPI_LONG,r,tagI,comm,&request);CHKERRQ(ierr);
        ierr = MPI_Wait(&request,&status);CHKERRQ(ierr);
        
        //printf("r[%d] requested [%ld,%ld] \n",r,ipackr[0],ipackr[0]+ipackr[1]);
        
        /* malloc */
        sbuffer = malloc(size*ipackr[1]);
        memset(sbuffer,0,size*ipackr[1]);

        /* read */
        {
          size_t nread;
          nread =  fread(sbuffer,size,ipackr[1],fp);
          if (nread != (size_t) ipackr[1]) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_READ,"fread() failure.");
        }
        
        /* send - must wait until MPI_Isend is finished before it is safe to free sbuffer */
        tagD = 2*r + 1;
        ierr = MPI_Isend(sbuffer,ipackr[1]*size,MPI_CHAR,r,tagD,comm,&requestS);CHKERRQ(ierr);
        ierr = MPI_Wait(&requestS,&statusS);CHKERRQ(ierr);
        //printf("r0 --> %d with tag %d\n",r,tagD);
        
        free(sbuffer);
      }

      byte_offset += ipackr[1]*size;
    }
    
  }
  if (rank != root) {
    tagI = 2*rank;
    ierr = MPI_Isend(ipack,2,MPI_LONG,root,tagI,comm,&requestS);CHKERRQ(ierr);
  }
  
  if (rank != root) {
    tagD = 2*rank + 1;
    //printf("r%d <-- r0 with tag %d\n",rank,tagD);
    ierr = MPI_Irecv(buffer,len*size,MPI_CHAR,root,tagD,comm,&request);CHKERRQ(ierr);
    ierr = MPI_Wait(&request,&status);CHKERRQ(ierr);
  }
  
  ierr = MPI_Barrier(comm);CHKERRQ(ierr);

  //if (rank == root) printf("<out> filepos %ld \n",ftell(fp));
  
  if (data) { *data = buffer; }
  
  PetscFunctionReturn(0);
}

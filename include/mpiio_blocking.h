
#ifndef __mpiio_blocking_h__
#define __mpiio_blocking_h__

#include <petsc.h>

PetscErrorCode MPIWrite_Blocking(FILE *fp,void *data,long int len,size_t size,int root,PetscBool skip_header,MPI_Comm comm);
PetscErrorCode MPIRead_Blocking(FILE *fp,void **data,long int len,size_t size,int root,PetscBool skip_header,MPI_Comm comm);

#endif


#ifndef __ptable_h__
#define __ptable_h__

typedef enum { PTABLE_UNINIT = -1, PTABLE_DENSE, PTABLE_SPARSE } PTableType;

typedef struct _p_PTable *PTable;

struct _p_PTable {
  PTableType  type;
  MPI_Comm    comm;
  long int    size;
  long int    start,end;
  PTable      cache;
  void        *data;
  PetscBool   issetup;
  PetscMPIInt commsize;
  long int    *start_all,*end_all;
  /* operations */
  PetscErrorCode (*create)(PTable);
  PetscErrorCode (*destroy)(PTable);
  PetscErrorCode (*setvalue)(PTable,long int,long int);
  PetscErrorCode (*getvalues)(PTable,long int,long int*,const long int**);
  PetscErrorCode (*hasvalues)(PTable,long int,PetscBool*);
  PetscErrorCode (*sync)(PTable);
  PetscBool issynchronized;
};


PetscErrorCode PTableCreate(MPI_Comm comm,PTable *p);
PetscErrorCode PTableDestroy(PTable *p);
PetscErrorCode PTableSetType(PTable p,PTableType type);
PetscErrorCode PTableSetRange(PTable p,long int start,long int end);
PetscErrorCode PTableGetRange(PTable p,long int *start,long int *end);
PetscErrorCode PTableSetValue(PTable p,long int row,long int val_j);
PetscErrorCode PTableGetValues(PTable p,long int row,long int *len,const long int *vals[]);
PetscErrorCode PTableHasValues(PTable p,long int row,PetscBool *found);
PetscErrorCode PTableGetNumberOfEntriesLocal(PTable p,long int *nentries);
PetscErrorCode PTableSynchronize(PTable p);
PetscErrorCode PTableSetup(PTable p);
PetscErrorCode PTableFlattenIntoIS(PTable p,PetscBool ignore_empty_slots,IS *fis,IS *tis);

#endif

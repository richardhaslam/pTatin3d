
/* ----------------------------------------------------------------------------------------------------------------
 
 PTable: A parallel table implementation for long int datatype's.
         PTable stores the following type of data
         row -> (value_0, value_1, ...)
         where row, value_i are all of type long int.
         value_i != value_j
 
 - PTable supports two storage mechanisms; dense (PTABLE_DENSE) and sparse (PTABLE_SPARSE)
 Support for sparse storage is provided by the khash library.
 
 - Two insertition modes were explored, 
   a) PTableValues_Insert_default(), and
     + PTableValues_Insert_default() checks for the existence of "value" in row, and only inserts a new
       entry if "value" does not exact. It uses exact memory allocation (no over allocation).
       This implementation is found to become very slow for large numbers of values (say > 1e4) due to the
       brute force search for existence employed.
   b) PTableValues_Insert_greedy()
     + PTableValues_Insert_greedy() skips the check for existence of "value" and simply inserts all values
       into row. Periodically, the row values are (i) sorted, (ii) duplicates removed and the
       (iii) memory allocated is adjusted to be exact. This is also enforced during PTableSynchronize().
       See _PTableSortPurgeValues().
 
 - In parallel, a sparse PTable is created to cache values which live on another MPI rank. 
   The contents of this table are scattered to the correct MPI rank, and the table subsequently flushed
   during PTableSynchronize().
 
   ---------------------------------------------------------------------------------------------------------------- */

#include <petsc.h>
#include <khash.h>
#include <data_exchanger.h>
#include <ptable.h>

typedef struct {
  long int len,allocated_len;
  long int *values; /* this can be a contiguous array defined by [len,allocated_len,value0,value1,...,] */
} PTableValues;

typedef struct {
  long int     nentries;
  PTableValues *row;
} PTable_DENSE;

/* Name associated with this khash table is "ptable64" */
KHASH_MAP_INIT_INT64(ptable64, PTableValues*)

typedef struct {
  long int          nentries;
  khash_t(ptable64) *hash;
  long int          *key_list,key_list_len;
  PetscBool         key_list_sorted;
} PTable_SPARSE;


PetscErrorCode PTableSparse_PrepareKeyList(PTable_SPARSE *sparse);

#define PT_SWAP(a,b,t) {t=a;a=b;b=t;}

#define PT_MEDIAN3(v,a,b,c)                        \
(v[a]<v[b]                                    \
? (v[b]<v[c]                                 \
? b                                       \
: (v[a]<v[c] ? c : a))                    \
: (v[c]<v[b]                                 \
? b                                       \
: (v[a]<v[c] ? a : c)))

#define PT_MEDIAN(v,right) PT_MEDIAN3(v,right/4,right/2,right/4*3)


/*
 A simple version of quicksort; taken from Kernighan and Ritchie, page 87.
 Assumes 0 origin for v, number of elements = right+1 (right is index of
 right-most member).
 */
static void PetscSortLongInt_Private(long int *v,long int right)
{
  long int i,j;
  long int pivot,tmp;

  if (right <= 1) {
    if (right == 1) {
      if (v[0] > v[1]) PT_SWAP(v[0],v[1],tmp);
    }
    return;
  }
  i = PT_MEDIAN(v,right);          /* Choose a pivot */
  PT_SWAP(v[0],v[i],tmp);          /* Move it out of the way */
  pivot = v[0];
  for (i=0,j=right+1;; ) {
    while (++i < j && v[i] <= pivot) ; /* Scan from the left */
    while (v[--j] > pivot) ;           /* Scan from the right */
    if (i >= j) break;
    PT_SWAP(v[i],v[j],tmp);
  }
  PT_SWAP(v[0],v[j],tmp);          /* Put pivot back in place. */
  PetscSortLongInt_Private(v,j-1);
  PetscSortLongInt_Private(v+j+1,right-(j+1));
}

PetscErrorCode PetscSortedRemoveDupsLongInt(long int *n,long int ii[])
{
  long int i,s = 0,N = *n, b = 0;

  PetscFunctionBegin;
  for (i=0; i<N-1; i++) {
    if (PetscUnlikely(ii[b+s+1] < ii[b])) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Input array is not sorted");
    if (ii[b+s+1] != ii[b]) {
      ii[b+1] = ii[b+s+1]; b++;
    } else s++;
  }
  *n = N - s;
  PetscFunctionReturn(0);
}


static inline PetscErrorCode PTableValues_Insert_default(PTableValues *ptv,long int v)
{
  long int  i;
  PetscBool found = PETSC_FALSE;
  long int  *tmp;

  for (i=0; i<ptv->len; i++) {
    if (ptv->values[i] == v) {
      found = PETSC_TRUE;
      break;
    }
  }
  if (found == PETSC_TRUE) PetscFunctionReturn(0);

  tmp = (long int*)realloc(ptv->values,sizeof(long int)*(ptv->len+1));
  ptv->allocated_len = ptv->len+1;
  ptv->values = tmp;
  ptv->values[ptv->len] = v;
  ptv->len++;
  
  PetscFunctionReturn(0);
}

PetscErrorCode PTableValues_SortPurge(PTableValues *ptv)
{
  long int  *tmp;

  //printf("sort\n");
  PetscSortLongInt_Private(ptv->values,ptv->len-1);
  PetscSortedRemoveDupsLongInt(&ptv->len,ptv->values);
  tmp = (long int*)realloc(ptv->values,sizeof(long int)*(ptv->len));
  ptv->values = tmp;
  ptv->allocated_len = ptv->len;
  PetscFunctionReturn(0);
}

static inline PetscErrorCode PTableValues_Insert_greedy(PTableValues *ptv,long int v)
{
  long int  *tmp;
  static long int counter = 0;

  tmp = (long int*)realloc(ptv->values,sizeof(long int)*(ptv->len+1));
  ptv->allocated_len = ptv->len+1;
  ptv->values = tmp;
  ptv->values[ptv->len] = v;
  ptv->len++;

  counter++;

  if (counter == 100000) {
    PTableValues_SortPurge(ptv);
    counter = 0;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PTableValues_Insert(PTableValues *ptv,long int v)
{
  PTableValues_Insert_default(ptv,v);
  //PTableValues_Insert_greedy(ptv,v);
  PetscFunctionReturn(0);
}

PetscErrorCode PTableValues_New(PTableValues *ptv)
{
  ptv->len = 0;
  ptv->allocated_len = 1;
  //ierr = PetscMalloc1(1,dense->row[i].values);CHKERRQ(ierr);
  ptv->values = (long int*)malloc(sizeof(long int));
  ptv->values[0] = 0;
  PetscFunctionReturn(0);
}

PetscErrorCode PTableValues_Destroy(PTableValues *ptv)
{
  if (!ptv) PetscFunctionReturn(0);
  free(ptv->values);
  ptv->values = NULL;
  ptv->len = 0;
  ptv->allocated_len = 0;
  PetscFunctionReturn(0);
}


PetscErrorCode _PTableSortPurgeValues(PTable table)
{
  PetscErrorCode ierr;
  long int       k;

  if (table->type == PTABLE_DENSE) {
    PTable_DENSE *dense = (PTable_DENSE*)table->data;

    for (k=0; k<table->size; k++) {
      PTableValues *row = NULL;

      row = &dense->row[k];
      ierr = PTableValues_SortPurge(row);CHKERRQ(ierr);
    }
  } else {
    PTable_SPARSE *sparse = (PTable_SPARSE*)table->data;

    ierr = PTableSparse_PrepareKeyList(sparse);CHKERRQ(ierr);

    /* traverse and pluck */
    for (k=0; k<sparse->key_list_len; k++) {
      long int     key;
      khiter_t     iterator;
      PTableValues *row = NULL;

      key = sparse->key_list[k];
      iterator = kh_get(ptable64, sparse->hash, key);
      row = kh_value(sparse->hash, iterator);
      ierr = PTableValues_SortPurge(row);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/* ========================================== */
/* implementation: dense */
/* ========================================== */
PetscErrorCode PTableSetValue_DENSE(PTable p,long int r,long int v)
{
  PTable_DENSE   *dense = (PTable_DENSE*)p->data;
  PetscErrorCode ierr;
  long int       rid;
  PTableValues   *row;

  rid = r - p->start;
  row = &dense->row[rid];
  ierr = PTableValues_Insert(row,v);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PTableHasValues_DENSE(PTable p,long int r,PetscBool *found)
{
  *found = PETSC_TRUE;
  PetscFunctionReturn(0);
}

PetscErrorCode PTableGetValues_DENSE(PTable p,long int r,long int *len,const long int *v[])
{
  PTable_DENSE   *dense = (PTable_DENSE*)p->data;
  long int       rid;

  rid = r - p->start;
  if (len) {
    *len = dense->row[ rid ].len;
  }
  if (v) {
    *v = (const long int*)dense->row[ rid ].values;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PTableCreate_DENSE(PTable p)
{
  PTable_DENSE   *dense;
  long int       i;
  PetscErrorCode ierr;

  ierr = PetscMalloc1(1,&dense);CHKERRQ(ierr);
  p->data = (void*)dense;
  dense->nentries = p->size;
  ierr = PetscMalloc1(dense->nentries,&dense->row);CHKERRQ(ierr);
  for (i=0; i<dense->nentries; i++) {
    dense->row[i].len = 0;
    dense->row[i].allocated_len = 1;
    //ierr = PetscMalloc1(1,dense->row[i].values);CHKERRQ(ierr);
    dense->row[i].values = (long int*)malloc(sizeof(long int));
    dense->row[i].values[0] = 0;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PTableDestroy_DENSE(PTable p)
{
  PTable_DENSE   *dense = (PTable_DENSE*)p->data;
  long int       i;
  PetscErrorCode ierr;

  for (i=0; i<dense->nentries; i++) {
    ierr = PTableValues_Destroy(&dense->row[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(dense->row);CHKERRQ(ierr);
  ierr = PetscFree(dense);CHKERRQ(ierr);
  p->data = NULL;
  PetscFunctionReturn(0);
}

/* ========================================== */
/* implementation: sparse */
/* ========================================== */

PetscErrorCode PTableSparse_PrepareKeyList(PTable_SPARSE *sparse)
{
  PetscErrorCode ierr;
  khiter_t       k0,k1,k;
  khint_t        len;
  long int       key,cnt;

  if (sparse->key_list) {
    ierr = PetscFree(sparse->key_list);CHKERRQ(ierr);
  }

  k0 = kh_begin(sparse->hash);
  k1 = kh_end(sparse->hash);
  len = kh_size(sparse->hash);

  ierr = PetscMalloc1(len,&sparse->key_list);CHKERRQ(ierr);
  ierr = PetscMemzero(sparse->key_list,sizeof(long int)*len);CHKERRQ(ierr);
  cnt = 0;
  for (k=k0; k<k1; k++) {
    if (kh_exist(sparse->hash, k)) {
      key = kh_key(sparse->hash, k);
      if (key >= 0) {
        sparse->key_list[cnt] = key;
        cnt++;
      }
    }
  }
  sparse->key_list_len = cnt;

  /* sort keys */
  PetscSortLongInt_Private(sparse->key_list,cnt-1);
  sparse->key_list_sorted = PETSC_TRUE;
  PetscFunctionReturn(0);
}

PetscErrorCode PTableSetValue_SPARSE(PTable p,long int r,long int v)
{
  PTable_SPARSE  *sparse = (PTable_SPARSE*)p->data;
  PetscErrorCode ierr;
  PTableValues   *row = NULL;
  khiter_t       iterator;
  long int       key;

  key = r - p->start; /* not sure it is a clever idea to offset this */

  iterator = kh_get(ptable64, sparse->hash, key);
  if (!kh_exist(sparse->hash, iterator)) {
    int ret;

    ierr = PetscMalloc1(1,&row);CHKERRQ(ierr);
    ierr = PTableValues_New(row);CHKERRQ(ierr);

    iterator = kh_put(ptable64, sparse->hash, key, &ret);
    (kh_value(sparse->hash, iterator)) = row;
    sparse->nentries++;
  } else {
    row = kh_value(sparse->hash, iterator);
  }
  if (!row) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Key %ld maps to a NULL row entry in the table. An initialization error has occurred",key);
  }
  ierr = PTableValues_Insert(row,v);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PTableHasValues_SPARSE(PTable p,long int r,PetscBool *found)
{
  PTable_SPARSE  *sparse = (PTable_SPARSE*)p->data;
  khiter_t       iterator;
  long int       key;

  key = r - p->start;
  iterator = kh_get(ptable64, sparse->hash, key);
  if (!kh_exist(sparse->hash, iterator)) {
    *found  = PETSC_FALSE;
  } else {
    *found = PETSC_TRUE;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PTableGetValues_SPARSE(PTable p,long int r,long int *len,const long int *v[])
{
  PTable_SPARSE  *sparse = (PTable_SPARSE*)p->data;
  PTableValues   *row = NULL;
  khiter_t       iterator;
  long int       key;

  key = r - p->start;

  iterator = kh_get(ptable64, sparse->hash, key);
  //if (!kh_exist(sparse->hash, iterator)) {
  //  SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Key %ld does not exist in table. Use PTableQueryValues() if you want to double-check prior to calling PTableGetValues()",key);
  //}
  if (!kh_exist(sparse->hash, iterator)) {
    if (len) { *len = 0; }
    if (v) {   *v = NULL; }
  }

  row = kh_value(sparse->hash, iterator);
  if (!row) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Key %ld maps to a NULL row entry in the table. An initialization error has occurred",key);
  }

  if (len) {
    *len = row->len;
  }
  if (v) {
    *v = (const long int*)row->values;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PTableCreate_SPARSE(PTable p)
{
  PTable_SPARSE  *sparse;
  PetscErrorCode ierr;
  int            ret;
  khiter_t       iterator;
  long int       key = -1;

  ierr = PetscMalloc1(1,&sparse);CHKERRQ(ierr);
  p->data = (void*)sparse;
  sparse->nentries = 0;

  sparse->hash = kh_init(ptable64);
  sparse->nentries = 0;

  /* insert key into slot 1 */
  iterator = kh_put(ptable64, sparse->hash, key, &ret);
  {
    PTableValues *row_dummy;

    ierr = PetscMalloc1(1,&row_dummy);CHKERRQ(ierr);
    ierr = PTableValues_New(row_dummy);CHKERRQ(ierr);

    (kh_value(sparse->hash, iterator)) = row_dummy;
  }

  sparse->key_list = NULL;
  sparse->key_list_sorted = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode PTableDestroy_SPARSE(PTable p)
{
  PTable_SPARSE  *sparse = (PTable_SPARSE*)p->data;
  PTableValues   *row = NULL;
  khiter_t       iterator;
  long int       key;
  PetscErrorCode ierr;
  khiter_t k,k0,k1;

  if (sparse->key_list) { ierr = PetscFree(sparse->key_list);CHKERRQ(ierr); }

  k0 = kh_begin(sparse->hash);
  k1 = kh_end(sparse->hash);

  for (k=k0; k<k1; k++) {
    key = kh_key(sparse->hash, k);
    if (kh_exist(sparse->hash, k)) {
      iterator = kh_get(ptable64, sparse->hash, key);
      row = kh_value(sparse->hash, iterator);
      ierr = PTableValues_Destroy(row);CHKERRQ(ierr);
      ierr = PetscFree(row);CHKERRQ(ierr);
    }
  }
  kh_destroy(ptable64, sparse->hash);
  sparse->hash = NULL;

  ierr = PetscFree(sparse);CHKERRQ(ierr);
  p->data = NULL;
  PetscFunctionReturn(0);
}

/* ========================================== */
/* =============== PTable API =============== */
/* ========================================== */

PetscErrorCode PTableView_Self(PTable p)
{
  long int k,i,nentries;
  PetscBool has_row;

  if (!p) PetscFunctionReturn(0);

  printf("PTableView:\n");
  PetscPrintf(p->comm,"  commsize: %D\n",(PetscInt)p->commsize);
  if (p->type == PTABLE_DENSE ) {  printf("  Type: \"dense\"\n"); }
  if (p->type == PTABLE_SPARSE ) { printf("  Type: \"sparse\"\n"); }
  printf("  start: %ld\n",p->start);
  printf("  end:   %ld\n",p->end);
  printf("  size:  %ld\n",p->size);
  PTableGetNumberOfEntriesLocal(p,&nentries);
  printf("  nentries:  %ld\n",nentries);
  //printf("nentries %ld ",p->end);
  for (k=p->start; k<p->end; k++) {
    PTableHasValues(p,k,&has_row);
    if (!has_row) {
      //printf("[r%.6ld] : \"empty\"\n",k);
    } else {
      long int       rl;
      const long int *rv;

      printf("row %.6ld ",k);
      PTableGetValues(p,k,&rl,&rv);
      printf(": nentries %.6ld : ",rl);

      // viewer for short arrary of values
      if (rl < 20) {
        for (i=0; i<rl-1; i++) {
          printf("%ld,",rv[i]);
        }
        printf("%ld",rv[i]);
        printf("\n");
      }
      // viewer for longer arrary of values
      if (rl >= 20) {
        long int nnr = rl/20;
        
        for (long int kk=0; kk<nnr; kk++) {
          printf("\n    ");
          for (long int ii=0; ii<20; ii++) {
            printf("%ld,",rv[20*kk+ii]);
          }
        }
        if (rl - nnr*20 > 0) { /* non zero remainder */
          printf("\n    ");
          for (long int ii=nnr*20; ii<rl-1; ii++) {
            printf("%ld,",rv[ii]);
          }
          printf("%ld\n",rv[rl-1]);
        }
      }

    }
  }

  PetscFunctionReturn(0);
}

/* collective */
PetscErrorCode PTableView(PTable p)
{
  PetscMPIInt commrank;
  PetscErrorCode ierr;
  long int k,i;
  PetscBool has_row;

  if (!p) PetscFunctionReturn(0);

  ierr = MPI_Comm_rank(p->comm,&commrank);CHKERRQ(ierr);

  PetscPrintf(p->comm,"PTableView:\n");
  if (p->type == PTABLE_DENSE ) {  PetscPrintf(p->comm,"  Type: \"dense\"\n"); }
  if (p->type == PTABLE_SPARSE ) { PetscPrintf(p->comm,"  Type: \"sparse\"\n"); }
  //PetscSynchronizedPrintf(p->comm,"  [%d] start: %ld\n",commrank,p->start);
  //PetscSynchronizedPrintf(p->comm,"  [%d] end:   %ld\n",commrank,p->end);
  PetscPrintf(p->comm,"  commsize: %D\n",(PetscInt)p->commsize);
  if (p->commsize > 1) {
    PetscPrintf(p->comm,"  ranges:\n    [");
    for (k=0; k<p->commsize; k++) {
      PetscPrintf(p->comm,"%ld,",p->start_all[k]);
    }
    PetscPrintf(p->comm,"%ld]\n",p->end_all[p->commsize-1]);
  } else {
    PetscPrintf(p->comm,"  ranges:\n    [%ld,%ld]\n",p->start,p->end);
  }
  //PetscSynchronizedPrintf(p->comm,"  [%d] size:  %ld\n",commrank,p->size);
  //PTableGetNumberOfEntriesLocal(p,&nentries);
  //PetscSynchronizedPrintf(p->comm,"  [%d] nentries:  %ld\n",commrank,nentries);
  //PetscSynchronizedFlush(p->comm,PETSC_STDOUT);
  for (k=p->start; k<p->end; k++) {
    PTableHasValues(p,k,&has_row);
    if (!has_row) {
      //printf("[r%.6ld] : \"empty\"\n",k);
    } else {
      long int       rl;
      const long int *rv;

      PetscSynchronizedPrintf(p->comm,"row %.6ld ",k);
      PTableGetValues(p,k,&rl,&rv);
      PetscSynchronizedPrintf(p->comm,": nentries %.6ld : ",rl);

      // viewer for short arrary of values
      if (rl < 20) {
        for (i=0; i<rl-1; i++) {
          PetscSynchronizedPrintf(p->comm,"%ld,",rv[i]);
        }
        PetscSynchronizedPrintf(p->comm,"%ld",rv[i]);
        PetscSynchronizedPrintf(p->comm,"\n");
      }
      // viewer for longer arrary of values
      if (rl >= 20) {
        long int nnr = rl/20;

        for (long int kk=0; kk<nnr; kk++) {
          PetscSynchronizedPrintf(p->comm,"\n    ");
          for (long int ii=0; ii<20; ii++) {
            PetscSynchronizedPrintf(p->comm,"%ld,",rv[20*kk+ii]);
          }
        }
        if (rl - nnr*20 > 0) { /* non zero remainder */
          PetscSynchronizedPrintf(p->comm,"\n    ");
          for (long int ii=nnr*20; ii<rl-1; ii++) {
            PetscSynchronizedPrintf(p->comm,"%ld,",rv[ii]);
          }
          PetscSynchronizedPrintf(p->comm,"%ld\n",rv[rl-1]);
        }
      }

    }
  }
  PetscSynchronizedFlush(p->comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}

PetscErrorCode PTableViewLite(PTable p)
{
  PetscMPIInt commrank;
  PetscErrorCode ierr;
  long int k,nentries,nentries_g;
  PetscBool has_row;

  if (!p) PetscFunctionReturn(0);

  ierr = MPI_Comm_rank(p->comm,&commrank);CHKERRQ(ierr);

  PetscPrintf(p->comm,"PTableViewLite:\n");
  if (p->type == PTABLE_DENSE ) {  PetscPrintf(p->comm,"  Type: \"dense\"\n"); }
  if (p->type == PTABLE_SPARSE ) { PetscPrintf(p->comm,"  Type: \"sparse\"\n"); }
  PetscPrintf(p->comm,"  commsize: %D\n",(PetscInt)p->commsize);
  PTableGetNumberOfEntriesLocal(p,&nentries);
  ierr = MPI_Reduce(&nentries,&nentries_g,1,MPI_LONG,MPI_SUM,0,p->comm);CHKERRQ(ierr);
  PetscPrintf(p->comm,"  nentries:  %ld\n",nentries_g);
  if (p->commsize > 1) {
    PetscPrintf(p->comm,"  ranges:\n    [");
    for (k=0; k<p->commsize; k++) {
      PetscPrintf(p->comm,"%ld,",p->start_all[k]);
    }
    PetscPrintf(p->comm,"%ld]\n",p->end_all[p->commsize-1]);
  } else {
    PetscPrintf(p->comm,"  ranges:\n    [%ld,%ld]\n",p->start,p->end);
  }
  for (k=p->start; k<p->end; k++) {
    PTableHasValues(p,k,&has_row);
    if (!has_row) {
      //printf("[r%.6ld] : \"empty\"\n",k);
    } else {
      long int       rl;
      const long int *rv;

      PTableGetValues(p,k,&rl,&rv);
      PetscSynchronizedPrintf(p->comm,"row %.6ld : nentries %.6ld\n",k,rl);
    }
  }
  PetscSynchronizedFlush(p->comm,PETSC_STDOUT);

  PetscFunctionReturn(0);
}

PetscErrorCode PTableDestroy(PTable *pt)
{
  PetscErrorCode ierr;
  PTable p;

  if (!pt) PetscFunctionReturn(0);
  p = *pt;

  if (p->start_all) { ierr = PetscFree(p->start_all);CHKERRQ(ierr); }
  if (p->end_all) { ierr = PetscFree(p->end_all);CHKERRQ(ierr); }
  if (p->cache) {
    ierr = PTableDestroy(&p->cache);CHKERRQ(ierr);
  }
  ierr = p->destroy(p);CHKERRQ(ierr);
  ierr = PetscFree(p);CHKERRQ(ierr);
  *pt = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode PTableCreate(MPI_Comm comm,PTable *p)
{
  PTable _p;
  PetscErrorCode ierr;


  ierr = PetscMalloc1(1,&_p);CHKERRQ(ierr);
  ierr = PetscMemzero(_p,sizeof(struct _p_PTable));CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&_p->commsize);CHKERRQ(ierr);

  _p->comm = comm;
  _p->start = -2;
  _p->end   = -1;
  _p->size  = 0;
  _p->type = PTABLE_UNINIT;
  _p->issetup = PETSC_FALSE;
  _p->cache = NULL;
  _p->start_all = NULL;
  _p->end_all   = NULL;
  _p->issynchronized = PETSC_FALSE;

  if (_p->commsize > 1) {
    ierr = PTableCreate(PETSC_COMM_SELF,&_p->cache);CHKERRQ(ierr);
    ierr = PTableSetType(_p->cache,PTABLE_SPARSE);CHKERRQ(ierr);
    ierr = PetscMalloc1(_p->commsize,&_p->start_all);CHKERRQ(ierr);
    ierr = PetscMalloc1(_p->commsize,&_p->end_all);CHKERRQ(ierr);
  }

  *p = _p;
  PetscFunctionReturn(0);
}

/* 
 Take all entries which are stored in p->cache
 Pack and send them to the correct rank
 Remove entries from p->cache upon recv
*/
PetscErrorCode PTableSync_MPI(PTable p)
{
  PetscErrorCode ierr;
  PTable         cache;
  PTable_SPARSE  *sparse;
  PetscMPIInt    commrank,j;
  long int       k,key;
  DataEx         de;
  long int       *pack_buffer,*recv_buffer;
  PetscInt       recv_length,cnt;


  ierr = MPI_Comm_rank(p->comm,&commrank);CHKERRQ(ierr);
  cache = p->cache;

  /* should do a sort/purge on the cached table */
  ierr = _PTableSortPurgeValues(cache);CHKERRQ(ierr);
  p->cache->issynchronized = PETSC_TRUE; /* cache is sequential so it always up to date - furthermore, one should never access it via GetValues */

  sparse = (PTable_SPARSE*)cache->data;
  ierr = PTableSparse_PrepareKeyList(sparse);CHKERRQ(ierr);

  //for (k=0; k<sparse->key_list_len; k++) {
  //  printf("r %d : k %ld : key %ld \n",commrank,k,sparse->key_list[k]);
  //}

  /* =================== DSDE =================== */
  // [dsde:create]
  de = DataExCreate(p->comm,0);

  // [dsde:topology]
  ierr = DataExTopologyInitialize(de);CHKERRQ(ierr);

  /* This is a hack and simple done because I was lazy */
  for (j=0; j<p->commsize; j++) {
    PetscMPIInt rank_j = (PetscMPIInt)j;
    if (rank_j != commrank) {
      ierr = DataExTopologyAddNeighbour(de,rank_j);CHKERRQ(ierr);
    }
  }
  ierr = DataExTopologyFinalize(de);CHKERRQ(ierr);

  // [dsde:send counts]
  ierr = DataExInitializeSendCount(de);CHKERRQ(ierr);

  for (k=0; k<sparse->key_list_len; k++) {
    long int len;
    PetscMPIInt target_rank;

    key = sparse->key_list[k];

    //ierr = PTableGetValues(cache,key,&len,NULL);CHKERRQ(ierr);
    ierr = cache->getvalues(cache,key,&len,NULL);CHKERRQ(ierr);

    target_rank = -1;
    for (j=0; j<p->commsize; j++) {
      if (key >= p->start_all[j]) {
        if (key < p->end_all[j]) {
          target_rank = j;
          break;
        }
      }
    }

    ierr = DataExAddToSendCount(de, target_rank, len+2 );CHKERRQ(ierr);

  }

  ierr = DataExFinalizeSendCount(de);CHKERRQ(ierr);

  // [dsde:pack]
  ierr = DataExPackInitialize(de,sizeof(long int));CHKERRQ(ierr);

  for (k=0; k<sparse->key_list_len; k++) {
    long int len;
    const long int *vals;
    PetscMPIInt target_rank;

    key = sparse->key_list[k];

    //ierr = PTableGetValues(cache,key,&len,&vals);CHKERRQ(ierr);
    ierr = cache->getvalues(cache,key,&len,&vals);CHKERRQ(ierr);

    target_rank = -1;
    for (j=0; j<p->commsize; j++) {
      if (key >= p->start_all[j]) {
        if (key < p->end_all[j]) {
          target_rank = j;
          break;
        }
      }
    }


    ierr = PetscMalloc1(len+2,&pack_buffer);CHKERRQ(ierr);

    pack_buffer[0] = key;
    pack_buffer[1] = len;
    //if (commrank == 0){ printf("[p] key %ld : len %ld : v \n",key,len); }

    ierr = PetscMemcpy(&pack_buffer[2],vals,len*sizeof(long int));CHKERRQ(ierr);
    //for (j=2; j<len+2; j++) {
    //  if (commrank == 0){ printf("[p] key %ld : len %ld : v %ld \n",key,len,pack_buffer[j]); }
    //}

    for (j=0; j<len+2; j++) {
      ierr = DataExPackData( de, target_rank, 1,(void*)&pack_buffer[j] );CHKERRQ(ierr);
    }

    ierr = PetscFree(pack_buffer);CHKERRQ(ierr);
  }

  ierr = DataExPackFinalize(de);CHKERRQ(ierr);

  // [send-recv]
  ierr = DataExBegin(de);CHKERRQ(ierr);
  ierr = DataExEnd(de);CHKERRQ(ierr);


  // [un-pack]
  ierr = DataExGetRecvData( de, &recv_length, (void**)&recv_buffer );CHKERRQ(ierr);

  /*
  if (commrank == 1) {
    for (k=0; k<recv_length; k++) {
      printf("i %ld val %ld\n",k,recv_buffer[k]);
    }
  }
  */

  cnt = 0;
  while (cnt < recv_length) {
    long int len;

    key = recv_buffer[cnt]; cnt++;
    len = recv_buffer[cnt]; cnt++;
    for (k=0; k<len; k++) {
      long int v = recv_buffer[cnt];

      //if (commrank == 1){ printf("[up] key %ld : len %ld : v %ld\n",key,len,v); }
      ierr = PTableSetValue(p,key,v);CHKERRQ(ierr);
      cnt++;
    }
  }



  ierr = DataExDestroy(de);CHKERRQ(ierr);

  /* =================== Flush cache =================== */


  for (k=0; k<sparse->key_list_len; k++) {
    khiter_t iterator;

    key = sparse->key_list[k];
    iterator = kh_get(ptable64, sparse->hash, key);
    if (!kh_exist(sparse->hash, iterator)) {
      SETERRQ2(p->comm,PETSC_ERR_USER,"k %ld : key %ld - a valid value should have been found...\n",k,key);
    } else {
      PTableValues *row = NULL;

      row = kh_value(sparse->hash, iterator);
      //printf("   --> row %p\n",row);
      if (row) {
        free(row->values);
        ierr = PetscFree(row);CHKERRQ(ierr);
      }
      kh_del(ptable64, sparse->hash, iterator);
    }
  }
  sparse->nentries = (long int)kh_size(sparse->hash);
  sparse->nentries = sparse->nentries - 1; /* subtract one here to account for the dummy value with key = -1 */


  PetscFunctionReturn(0);
}

/* collective */
PetscErrorCode PTableSynchronize(PTable p)
{
  PetscErrorCode ierr;

  /* should do a sort/purge on the local table */
  ierr = _PTableSortPurgeValues(p);CHKERRQ(ierr);

  if (p->sync) {
    ierr = p->sync(p);CHKERRQ(ierr);
  }
  p->issynchronized = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/* collective */
PetscErrorCode PTableSetup(PTable p)
{
  PetscErrorCode ierr;

  if (p->issetup) PetscFunctionReturn(0);
  switch (p->type) {
    case PTABLE_UNINIT:
      SETERRQ(p->comm,PETSC_ERR_ORDER,"Table type has not been set");
      break;

    case PTABLE_DENSE:
      p->size = p->end - p->start;
      // set methods
      p->create    = PTableCreate_DENSE;
      p->setvalue  = PTableSetValue_DENSE;
      p->getvalues = PTableGetValues_DENSE;
      p->hasvalues = PTableHasValues_DENSE;
      p->destroy   = PTableDestroy_DENSE;
      break;

    case PTABLE_SPARSE:
      p->size = 0;
      // set methods
      p->create    = PTableCreate_SPARSE;
      p->setvalue  = PTableSetValue_SPARSE;
      p->getvalues = PTableGetValues_SPARSE;
      p->hasvalues = PTableHasValues_SPARSE;
      p->destroy   = PTableDestroy_SPARSE;
      break;

    default:
      SETERRQ(p->comm,PETSC_ERR_ORDER,"Table type has not been set");
      break;
  }
  p->sync = NULL;
  if (p->commsize > 1) {
    p->sync = PTableSync_MPI;
  }

  // create, setup context
  ierr = p->create(p);CHKERRQ(ierr);
  p->issetup = PETSC_TRUE;

  /* collect min,max ranges */
  if (p->cache) {
    long int g_start,g_end,k;

    ierr = MPI_Allreduce(&p->start,&g_start,1,MPI_LONG,MPI_MIN,p->comm);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&p->end,  &g_end,  1,MPI_LONG,MPI_MAX,p->comm);CHKERRQ(ierr);
    ierr = MPI_Allgather(&p->start,1,MPI_LONG,p->start_all,1,MPI_LONG,p->comm);CHKERRQ(ierr);
    ierr = MPI_Allgather(&p->end,  1,MPI_LONG,p->end_all,  1,MPI_LONG,p->comm);CHKERRQ(ierr);

    // 0 15
    // 10 30
    for (k=1; k<p->commsize; k++) {
      //if (p->start_all[k] < p->start_all[k-1]) SETERRQ(p->comm,PETSC_ERR_SUP,"Start ranges must increase monotonically");
      if (p->start_all[k] >= p->start_all[k-1]) {
        if (p->start_all[k] < p->end_all[k-1]) {
          SETERRQ(p->comm,PETSC_ERR_SUP,"Start range must increase monotonically. Ranges cannot overlap");
        }
      }
      if (p->end_all[k] >= p->start_all[k-1]) {
        if (p->end_all[k] < p->end_all[k-1]) {
          SETERRQ(p->comm,PETSC_ERR_SUP,"End range must increase monotonically. Ranges cannot overlap");
        }
      }
    }


    ierr = PTableSetRange(p->cache,g_start,g_end);CHKERRQ(ierr);
    ierr = PTableSetup(p->cache);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/* collective */
PetscErrorCode PTableSetType(PTable p,PTableType type)
{
  if (p->issetup) SETERRQ(p->comm,PETSC_ERR_ORDER,"Cannot change table type PTableSetup() called");
  p->type = type;
  PetscFunctionReturn(0);
}

/* not collective */
PetscErrorCode PTableSetRange(PTable p,long int start,long int end)
{
  if (p->issetup) SETERRQ(p->comm,PETSC_ERR_ORDER,"Cannot resize table after PTableSetup() called");
  p->start = start;
  p->end   = end;
  PetscFunctionReturn(0);
}

PetscErrorCode PTableGetRange(PTable p,long int *start,long int *end)
{
  if (!p->issetup) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not setup. Call PTableSetup() first");
  if (start) { *start = p->start; }
  if (end) {   *end   = p->end; }
  PetscFunctionReturn(0);
}

PetscErrorCode PTableSetValue_SEQ(PTable p,long int row,long int val_j)
{
  PetscErrorCode ierr;

  p->issynchronized = PETSC_TRUE;
  if (row >= p->start && row < p->end) {
    ierr = p->setvalue(p,row,val_j);CHKERRQ(ierr);
  } else {
    SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"row %ld is outside table range [%ld,%ld)",row,p->start,p->end);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PTableSetValue_MPI(PTable p,long int row,long int val_j)
{
  PetscErrorCode ierr;

  if (row < p->cache->start || row >= p->cache->end) {
    SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"row %ld is outside table range [%ld,%ld)",row,p->cache->start,p->cache->end);
  }

  if (row >= p->start && row < p->end) {
    p->issynchronized = PETSC_FALSE;
    ierr = p->setvalue(p,row,val_j);CHKERRQ(ierr);
  } else {
    p->issynchronized = PETSC_FALSE;
    p->cache->issynchronized = PETSC_TRUE; /* always true as SEQ */
    ierr = p->cache->setvalue(p->cache,row,val_j);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


/* not collective */
PetscErrorCode PTableSetValue(PTable p,long int row,long int val_j)
{
  PetscErrorCode ierr;

  if (!p->issetup) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not setup. Call PTableSetup() first");
  if (p->commsize > 1) {
    ierr = PTableSetValue_MPI(p,row,val_j);CHKERRQ(ierr);
  } else {
    ierr = PTableSetValue_SEQ(p,row,val_j);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* not collective */
PetscErrorCode PTableGetValues(PTable p,long int row,long int *len,const long int *vals[])
{
  PetscErrorCode ierr;
  if (!p->issetup) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not setup. Call PTableSetup() first");
  if (!p->issynchronized) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not synchronized. Call PTableSynchronize() first");
  if (row >= p->start && row < p->end) {
    ierr = p->getvalues(p,row,len,vals);CHKERRQ(ierr);
  } else {
    SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"row %ld is outside table range [%ld,%ld)",row,p->start,p->end);
  }
  PetscFunctionReturn(0);
}

/* not collective */
PetscErrorCode PTableHasValues(PTable p,long int row,PetscBool *found)
{
  PetscErrorCode ierr;
  if (!p->issetup) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not setup. Call PTableSetup() first");
  if (!p->issynchronized) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not synchronized. Call PTableSynchronize() first");
  if (row >= p->start && row < p->end) {
    *found = PETSC_FALSE;
    ierr = p->hasvalues(p,row,found);CHKERRQ(ierr);
  } else {
    SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"row %ld is outside table range [%ld,%ld)",row,p->start,p->end);
  }
  PetscFunctionReturn(0);
}

/* not collective */
PetscErrorCode PTableGetNumberOfEntriesLocal(PTable p,long int *nentries)
{
  if (!p->issetup) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not setup. Call PTableSetup() first");
  if (!p->issynchronized) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not synchronized. Call PTableSynchronize()  first");

  switch (p->type) {
    case PTABLE_DENSE:
    {
      PTable_DENSE *dense = (PTable_DENSE*)p->data;
      *nentries = dense->nentries;
    }
      break;

    case PTABLE_SPARSE:
    {
      PTable_SPARSE *sparse = (PTable_SPARSE*)p->data;
      *nentries = sparse->nentries;
    }
      break;

    default:
      SETERRQ(p->comm,PETSC_ERR_ORDER,"Table type has not been set");
      break;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode _PTableFlattenIntoIS_SPARSE(PTable p,PetscBool ignore_empty_slots,IS *from_is,IS *to_is)
{
  PTable_SPARSE *sparse = (PTable_SPARSE*)p->data;
  long int k,length;
  PetscInt i;
  PetscInt *indices_f,*indices_t;
  long int start = 0,end = 0;
  PetscErrorCode ierr;

  ierr = PTableSparse_PrepareKeyList(sparse);CHKERRQ(ierr);

  if (ignore_empty_slots) {
    ierr = PTableGetNumberOfEntriesLocal(p,&length);CHKERRQ(ierr);
  } else {
    ierr = PTableGetRange(p,&start,&end);CHKERRQ(ierr);
    length = end - start;
  }
  ierr = PetscMalloc1(length,&indices_f);CHKERRQ(ierr);
  ierr = PetscMalloc1(length,&indices_t);CHKERRQ(ierr);

  for (i=0; i<length; i++) {
    indices_f[i] = -1;
    indices_t[i] = -1;
  }

  /* traverse */
  if (ignore_empty_slots) {
    for (k=0; k<sparse->key_list_len; k++) {
      long int     key,len_k;
      const long int *val_k;

      key = sparse->key_list[k];
      ierr = PTableGetValues(p,key,&len_k,&val_k);CHKERRQ(ierr);
      if (len_k > 0 && val_k != NULL) {
        indices_f[k] = (PetscInt)(start + key); /* take care - countering acting the offset */
        indices_t[k] = (PetscInt)val_k[0];
      }
    }
  } else {
    for (k=0; k<sparse->key_list_len; k++) {
      long int     key,len_k;
      const long int *val_k;

      key = sparse->key_list[k];
      if (key > length) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"key index out of range");

      ierr = PTableGetValues(p,key,&len_k,&val_k);CHKERRQ(ierr);
      if (len_k > 0 && val_k != NULL) {
        /* note that the first arg "key" is local to the rank so key <= length */
        indices_f[key] = (PetscInt)(start + key);  /* take care - countering acting the offset */
        indices_t[key] = (PetscInt)val_k[0];
      }
    }
  }

  if (from_is) {
    ierr = ISCreateGeneral(p->comm,(PetscInt)length,indices_f,PETSC_COPY_VALUES,from_is);CHKERRQ(ierr);
  }
  if (to_is) {
    ierr = ISCreateGeneral(p->comm,(PetscInt)length,indices_t,PETSC_COPY_VALUES,to_is);CHKERRQ(ierr);
  }

  ierr = PetscFree(indices_f);CHKERRQ(ierr);
  ierr = PetscFree(indices_t);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode _PTableFlattenIntoIS_DENSE(PTable p,PetscBool ignore_empty_slots,IS *from_is,IS *to_is)
{
  long int k,length;
  PetscInt i;
  PetscInt *indices_f,*indices_t;
  long int start = 0,end = 0;
  PetscErrorCode ierr;

  ierr = PTableGetNumberOfEntriesLocal(p,&length);CHKERRQ(ierr);
  ierr = PTableGetRange(p,&start,&end);CHKERRQ(ierr);

  ierr = PetscMalloc1(length,&indices_f);CHKERRQ(ierr);
  ierr = PetscMalloc1(length,&indices_t);CHKERRQ(ierr);

  for (i=0; i<length; i++) {
    indices_f[i] = -1;
    indices_t[i] = -1;
  }

  /* traverse */
  for (k=start; k<end; k++) {
    long int     key,len_k;
    const long int *val_k;

    key = k;
    ierr = PTableGetValues(p,key,&len_k,&val_k);CHKERRQ(ierr);
    if (len_k > 0 && val_k != NULL) {
      indices_f[k-start] = (PetscInt)key;
      indices_t[k-start] = (PetscInt)val_k[0];
    }
  }

  if (from_is) {
    ierr = ISCreateGeneral(p->comm,(PetscInt)length,indices_f,PETSC_COPY_VALUES,from_is);CHKERRQ(ierr);
  }
  if (to_is) {
    ierr = ISCreateGeneral(p->comm,(PetscInt)length,indices_t,PETSC_COPY_VALUES,to_is);CHKERRQ(ierr);
  }

  ierr = PetscFree(indices_f);CHKERRQ(ierr);
  ierr = PetscFree(indices_t);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
 IS[local_index] = global_index
 collective
*/
PetscErrorCode PTableFlattenIntoIS(PTable p,PetscBool ignore_empty_slots,IS *fis,IS *tis)
{
  PetscErrorCode ierr;

  if (!p->issetup) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not setup. Call PTableSetup() first");
  if (!p->issynchronized) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Table is not synchronized. Call PTableSynchronize()  first");

  if (p->type == PTABLE_SPARSE) {
    ierr = _PTableFlattenIntoIS_SPARSE(p,ignore_empty_slots,fis,tis);CHKERRQ(ierr);
  } else if (p->type == PTABLE_DENSE) {
    ierr = _PTableFlattenIntoIS_DENSE(p,ignore_empty_slots,fis,tis);CHKERRQ(ierr);
  } else {
    SETERRQ(p->comm,PETSC_ERR_ORDER,"Table type has not been set");
  }

  PetscFunctionReturn(0);
}

void hashdemo(void)
{
  khash_t(ptable64)  *hash;
  PTableValues *dummy;
  PTableValues *row;
  int ret;
  khiter_t iterator;
  long int key;
  long int size;
  khint_t k;

  hash = kh_init(ptable64);

  dummy = malloc(sizeof(PTableValues));
  dummy->len = 1;
  dummy->allocated_len = 1;
  dummy->values = NULL;

  /* insert key into slot 1 */
  iterator = kh_put(ptable64, hash, -1, &ret);
  row = malloc(sizeof(PTableValues));
  row->len = -1;
  row->allocated_len = -1;
  row->values = NULL;
  (kh_value(hash, iterator)) = row;


  key = 10;
  row = malloc(sizeof(PTableValues));
  row->len = 10;
  row->allocated_len = 10;
  row->values = NULL;
  printf("key %ld : row %p\n",key,row);

  iterator = kh_put(ptable64, hash, key, &ret);
  (kh_value(hash, iterator)) = row;

  key = 100000;
  row = malloc(sizeof(PTableValues));
  row->len = 100000;
  row->allocated_len = 100000;
  row->values = NULL;
  printf("key %ld : row %p\n",key,row);

  iterator = kh_put(ptable64, hash, key, &ret);
  (kh_value(hash, iterator)) = row;

  key = -2;
  row = malloc(sizeof(PTableValues));
  row->len = -2;
  row->allocated_len = -2;
  row->values = NULL;
  printf("key %ld : row %p\n",key,row);

  iterator = kh_put(ptable64, hash, key, &ret);
  (kh_value(hash, iterator)) = row;


  key = 100;
  row = malloc(sizeof(PTableValues));
  row->len = 100;
  row->allocated_len = 100;
  row->values = NULL;
  printf("key %ld : row %p\n",key,row);

  iterator = kh_put(ptable64, hash, key, &ret);
  (kh_value(hash, iterator)) = row;

  key = 15;
  row = malloc(sizeof(PTableValues));
  row->len = 15;
  row->allocated_len = 15;
  row->values = NULL;
  printf("key %ld : row %p\n",key,row);

  iterator = kh_put(ptable64, hash, key, &ret);
  (kh_value(hash, iterator)) = row;


  row = NULL;

  /* fetch */
  key = 10;
  iterator = kh_get(ptable64, hash, key);
  row = kh_value(hash, iterator);
  printf("key %ld --> row %p\n",key,row);
  printf("key %ld --> %ld %ld\n",key,row->len,row->allocated_len);


  size = kh_size(hash);
  printf("size %ld \n",size);

  // delete thing with key = 10
  free(row);
  kh_del(ptable64, hash, iterator);

  size = kh_size(hash);
  printf("size %ld \n",size);
  //


  khiter_t k0,k1;

  k0 = kh_begin(hash);
  k1 = kh_end(hash);

  for (k=k0; k<k1; k++) {
    key = 0;
    key = kh_key(hash, k);
    //iterator = kh_get(ptable64, hash, key);
    if (kh_exist(hash, k)) {
      printf("k %u : key %ld\n",k,key);

      iterator = kh_get(ptable64, hash, key);
      row = kh_value(hash, iterator);
      printf("   --> row %p\n",row);
    }

  }

  kh_destroy(ptable64, hash);
}

/* testing code */
#if 0

PetscErrorCode ex1(void)
{
  PTable table;
  PetscErrorCode ierr;

  ierr = PTableCreate(PETSC_COMM_SELF,&table);CHKERRQ(ierr);
  ierr = PTableSetRange(table,0,30);CHKERRQ(ierr);
  ierr = PTableSetType(table,PTABLE_SPARSE);CHKERRQ(ierr);
  ierr = PTableSetup(table);CHKERRQ(ierr);

  ierr = PTableSetValue(table,5,0);CHKERRQ(ierr);
  ierr = PTableSetValue(table,6,19);CHKERRQ(ierr);
  // ierr = PTableSetValue(table,7,101);CHKERRQ(ierr);
  //ierr = PTableSetValue(table,7,-1);CHKERRQ(ierr);
  //for (int i=0; i<128*3; i++) {
  for (int i=0; i<5000000; i++) {
    ierr = PTableSetValue(table,7,i);CHKERRQ(ierr);
  }

  ierr = PTableSynchronize(table);CHKERRQ(ierr);

  //ierr = PTableView(table);CHKERRQ(ierr);
  //ierr = PTableViewLite(table);CHKERRQ(ierr);
  ierr = PTableDestroy(&table);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ex2(void)
{
  PTable table;
  PetscErrorCode ierr;
  PetscMPIInt rank;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  ierr = PTableCreate(PETSC_COMM_WORLD,&table);CHKERRQ(ierr);
  ierr = PTableSetRange(table,rank*30,rank*30+30);CHKERRQ(ierr);
  ierr = PTableSetType(table,PTABLE_SPARSE);CHKERRQ(ierr);
  ierr = PTableSetup(table);CHKERRQ(ierr);

  ierr = PTableSetValue(table,5,8870);CHKERRQ(ierr);
  /*
  ierr = PTableSetValue(table,5,0);CHKERRQ(ierr);
  ierr = PTableSetValue(table,6,19);CHKERRQ(ierr);
  // ierr = PTableSetValue(table,7,101);CHKERRQ(ierr);
  //ierr = PTableSetValue(table,7,-1);CHKERRQ(ierr);
  for (int i=0; i<21; i++) {
    ierr = PTableSetValue(table,7,i);CHKERRQ(ierr);
  }


  */
  if (rank == 0) {
    ierr = PTableSetValue(table,55,501);CHKERRQ(ierr);
    ierr = PTableSetValue(table,55,502);CHKERRQ(ierr);

    ierr = PTableSetValue(table,50,2);CHKERRQ(ierr);
    ierr = PTableSetValue(table,50,22);CHKERRQ(ierr);
    ierr = PTableSetValue(table,50,222);CHKERRQ(ierr);
    for (int i=0; i<128*3; i++) {
      ierr = PTableSetValue(table,59,i+100000);CHKERRQ(ierr);
    }
  }

  ierr = PTableSynchronize(table);CHKERRQ(ierr);

  MPI_Barrier(table->comm);
  ierr = PTableView(table);CHKERRQ(ierr);
  if (rank == 0) {
    //  ierr = PTableView_Self(table->cache);CHKERRQ(ierr);
  }

  ierr = PTableDestroy(&table);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ex3(void)
{
  PTable table;
  IS fis,tis;
  PetscErrorCode ierr;

  ierr = PTableCreate(PETSC_COMM_SELF,&table);CHKERRQ(ierr);
  ierr = PTableSetRange(table,0,30);CHKERRQ(ierr);
  ierr = PTableSetType(table,PTABLE_SPARSE);CHKERRQ(ierr);
  ierr = PTableSetup(table);CHKERRQ(ierr);

  ierr = PTableSetValue(table,5,23);CHKERRQ(ierr);
  ierr = PTableSetValue(table,6,19);CHKERRQ(ierr);
  ierr = PTableSetValue(table,7,101);CHKERRQ(ierr);

  ierr = PTableSynchronize(table);CHKERRQ(ierr);

  ierr = PTableView(table);CHKERRQ(ierr);

  ierr = PTableFlattenIntoIS(table,PETSC_TRUE,&fis,&tis);CHKERRQ(ierr);
  //ierr = PTableFlattenIntoIS(table,PETSC_FALSE,&fis,&tis);CHKERRQ(ierr);
  ierr = ISView(fis,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = ISView(tis,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = PTableDestroy(&table);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/* overlapping ranges */
PetscErrorCode ex4(void)
{
  PTable table;
  PetscErrorCode ierr;
  PetscMPIInt rank,size;

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only for comm.size=2");
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  ierr = PTableCreate(PETSC_COMM_WORLD,&table);CHKERRQ(ierr);
  if (rank == 0) {
    ierr = PTableSetRange(table,0,15);CHKERRQ(ierr);
  } else {
    ierr = PTableSetRange(table,10,30);CHKERRQ(ierr);
  }
  ierr = PTableSetType(table,PTABLE_SPARSE);CHKERRQ(ierr);
  ierr = PTableSetup(table);CHKERRQ(ierr);

  if (rank == 0) {
    ierr = PTableSetValue(table,5,8870);CHKERRQ(ierr);
  }
  if (rank == 1) {
    ierr = PTableSetValue(table,12,8);CHKERRQ(ierr);
  }

  ierr = PTableSynchronize(table);CHKERRQ(ierr);

  MPI_Barrier(table->comm);
  ierr = PTableView(table);CHKERRQ(ierr);
  ierr = PTableDestroy(&table);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc,&argv,(char *)0,NULL);CHKERRQ(ierr);

  //hashdemo();
  //ex1(); /* test sparse/dense in serial */
  ex2(); /* test sparse/dense in parallel */
  //ex3(); /* test for flatten PTable into IS */
  //ex4(); /* test for overlapping ranges */

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

#endif


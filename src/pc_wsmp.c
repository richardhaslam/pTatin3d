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
 **    filename:   pc_wsmp.c
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
 Notes:
 [1] Compatible with PWSSMP version V 17.01.01
 [2] PWSSMP appears to fail when using C indices.
     wsmp->IPARM[5 -1] = 0 does not appear to function correctly
     I also note that some of the factorization steps alter the value of wsmp->IPARM[5 -1].
     e.g. If you set it to zero, the symbolic factorization changes its value to 1
     If this gets fixed, search this file for "fortran" and correct the +1 / -1 accordingly
 [3] If you run with MPI-ranks = 1, by default this PC implementation will try to use WSSMP.
     If WSSMP hasn't found, it will call PWSSMP.
     Support with WSSMP is activated by the compiler flag
       -DTATIN_HAVE_WSSMP
     Support for PWSSMP is activated via the compiler flag
       -DTATIN_HAVE_PWSSMP
     In parallel (MPI-ranks > 1), you must have compiled the code with support for PWSSMP
 [4] For consistency with the WSMP documentation, the parameter names, e.g.,
       -pc_wsmp_iparm31
     relate to the parameters identified in the manual as e.g. IPARM[31].
     That is, the pc_wsmp option strings use the fortran index style (starting from 1)
*/

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petsc/private/pcimpl.h>     /*I "petscpc.h" I*/
#include <petscksp.h>

/* prototypes for symbols in WSSMP and PWSSMP */
extern void wsetmaxthrds_(int *nthreads);

extern void wssmp_ ( int *n, int ia[], int ja[], double a[], double diag[], int perm[], int iperm[], double b[], int *ldb, int *nrhs, int aux[], int *naux, int mrp[], int iparam[], double dparam[] );
extern void wsmp_clear_(void);
extern void wsffree_(void);

extern void pwssmp_ ( int *n, int ia[], int ja[], double a[], double diag[], int perm[], int iperm[], double b[], int *ldb, int *nrhs, int aux[], int *naux, int mrp[], int iparam[], double dparam[] );
extern void pwsmp_clear_(void);
extern void wsetmpicomm_(int *fcomm);
extern void pwsffree_(void);

/*
 If support for PWSSMP is defined, we automatically have support for WSSMP
 as all the symbols are contained in libpwsmp.so
 We thus set TATIN_HAVE_WSSMP. The consequence of this is that when running in serial,
 wssmp_() will be called instead of pwssmp_()
*/
#ifdef TATIN_HAVE_PWSSMP
  #ifndef TATIN_HAVE_WSSMP
    #define TATIN_HAVE_WSSMP
  #endif
#endif

typedef struct {
  PetscBool sequential,symmetric;
  int    IPARM[64];
  double DPARM[64];
  int    Nlocal; /* local size on processor */
  int    Nglobal;
  int    nnz;
  int    *IA,*JA;
  double *AVALS;
  double *DIAG; /* size Nlocal */
  int    *PERM,*INVP; /* size Nglobal */
  int    NRHS;
  double *B; /* placeholder for rhs and solution */
  int    LDB;
  int    AUX,NAUX;
  int    *MRP;
} PC_WSMP;


PetscErrorCode WSMPSetFromOptions_Ordering(PC pc,PC_WSMP *wsmp);
PetscErrorCode WSMPSetFromOptions_SymbolicFactorization(PC pc,PC_WSMP *wsmp);
PetscErrorCode WSMPSetFromOptions_NumericFactorization(PC pc,PC_WSMP *wsmp);
PetscErrorCode WSMPSetFromOptions_BackSubstitution(PC pc,PC_WSMP *wsmp);
PetscErrorCode WSMPSetFromOptions_IterativeRefinement(PC pc,PC_WSMP *wsmp);

#undef __FUNCT__
#define __FUNCT__ "PCWSMP_CheckCSR"
PetscErrorCode PCWSMP_CheckCSR(Mat A,PC_WSMP *wsmp)
{
  PetscErrorCode ierr;
  MPI_Comm       comm;
  Mat            B;
  PetscInt       i,m,M,_Nlocal,_nnz;
  PetscInt       *_ia,*_ja;
  PetscScalar    *_avals;
  Vec            x,y,y1;
  PetscReal      nrm;
  PetscScalar    range[2];
  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
  ierr = MatGetSize(A,&M,NULL);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m,NULL);CHKERRQ(ierr);
  
  ierr = MatCreate(comm,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,m,m,M,M);CHKERRQ(ierr);
  ierr = MatSetType(B,MATSBAIJ);CHKERRQ(ierr);
  
  _Nlocal = (PetscInt)wsmp->Nlocal;
  _nnz    = (PetscInt)wsmp->nnz;
  ierr = PetscMalloc1(_Nlocal+1,&_ia);CHKERRQ(ierr);
  ierr = PetscMalloc1(_nnz,&_ja);CHKERRQ(ierr);
  ierr = PetscMalloc1(_nnz,&_avals);CHKERRQ(ierr);
  
  for (i=0; i<_Nlocal+1; i++) {
    _ia[i] = (PetscInt)wsmp->IA[i] - 1; /* convert the fortran indices to C indices starting at 0 */
  }
  for (i=0; i<_nnz; i++) {
    _ja[i] = (PetscInt)wsmp->JA[i] - 1; /* convert the fortran indices to C indices starting at 0 */
    _avals[i] = (PetscScalar)wsmp->AVALS[i];
  }
  
  ierr = MatSeqSBAIJSetPreallocationCSR(B,1,_ia,_ja,_avals);CHKERRQ(ierr);
  ierr = MatMPISBAIJSetPreallocationCSR(B,1,_ia,_ja,_avals);CHKERRQ(ierr);
  
  ierr = MatCreateVecs(A,&x,&y);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&y1);CHKERRQ(ierr);
  ierr = VecSetRandom(x,NULL);CHKERRQ(ierr);
  
  ierr = MatMult(A,x,y);CHKERRQ(ierr);
  ierr = VecMin(y,NULL,&range[0]);CHKERRQ(ierr);
  ierr = VecMax(y,NULL,&range[1]);CHKERRQ(ierr);
  PetscPrintf(comm,"{wsmp} min(A.x): %+1.12e \n",range[0]);
  PetscPrintf(comm,"{wsmp} max(A.x): %+1.12e \n",range[1]);
  
  ierr = MatMult(B,x,y1);CHKERRQ(ierr);
  ierr = VecMin(y1,NULL,&range[0]);CHKERRQ(ierr);
  ierr = VecMax(y1,NULL,&range[1]);CHKERRQ(ierr);
  PetscPrintf(comm,"{wsmp} min(B.x): %+1.12e \n",range[0]);
  PetscPrintf(comm,"{wsmp} max(B.x): %+1.12e \n",range[1]);
  
  ierr = VecAXPY(y,-1.0,y1);CHKERRQ(ierr);
  ierr = VecNorm(y,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = VecMin(y,NULL,&range[0]);CHKERRQ(ierr);
  ierr = VecMax(y,NULL,&range[1]);CHKERRQ(ierr);
  PetscPrintf(comm,"{wsmp} nrm: %+1.12e \n",nrm);
  PetscPrintf(comm,"{wsmp} min: %+1.12e \n",range[0]);
  PetscPrintf(comm,"{wsmp} max: %+1.12e \n",range[1]);
  
  ierr = PetscFree(_ia);CHKERRQ(ierr);
  ierr = PetscFree(_ja);CHKERRQ(ierr);
  ierr = PetscFree(_avals);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = VecDestroy(&y1);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCWSMP_VecView"
PetscErrorCode PCWSMP_VecView(const char name[],PC_WSMP *wsmp)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  FILE           *fp = NULL;
  char           fname[PETSC_MAX_PATH_LEN];
  PetscInt       i;
  
  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  
  PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"%s_rank%D.dat",name,rank);
  fp = fopen(fname,"w");
  if (!fp) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",fname);
  }
  
  PetscFPrintf(PETSC_COMM_SELF,fp,"%D\n",wsmp->Nlocal);
  for (i=0; i<wsmp->Nlocal; i++) {
    PetscFPrintf(PETSC_COMM_SELF,fp,"%D %+1.8e\n",i,wsmp->B[i]);
  }
  fclose(fp);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCWSMP_CSRView"
PetscErrorCode PCWSMP_CSRView(Mat A,int ia[],int ja[],double aij[])
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  FILE           *fp;
  char           name[PETSC_MAX_PATH_LEN];
  PetscInt       i,m,nnz;
  MPI_Comm       comm;
  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m,NULL);CHKERRQ(ierr);
  
  /* ia */
  if (ia) {
    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"ia_rank%D.dat",rank);
    fp = fopen(name,"w");
    
    PetscFPrintf(PETSC_COMM_SELF,fp,"%D\n",m+1);
    for (i=0; i<m; i++) {
      PetscFPrintf(PETSC_COMM_SELF,fp,"%D %D\n",i,ia[i]);
    }
    PetscFPrintf(PETSC_COMM_SELF,fp,"%D %D\n",m,ia[m]);
    fclose(fp);
  }
  
  /* ja */
  if (ja) {
    if (!ia) {
      SETERRQ(comm,PETSC_ERR_USER,"Cannot write ja[] as ia[] isn't provided <require nnz");
    }
    nnz = (PetscInt)ia[m];
    
    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"ja_rank%D.dat",rank);
    fp = fopen(name,"w");
    
    PetscFPrintf(PETSC_COMM_SELF,fp,"%D\n",nnz);
    for (i=0; i<nnz; i++) {
      PetscFPrintf(PETSC_COMM_SELF,fp,"%D %D\n",i,ja[i]);
    }
    fclose(fp);
  }
  
  /* aij */
  if (aij) {
    if (!ia) {
      SETERRQ(comm,PETSC_ERR_USER,"Cannot write aij[] as ia[] isn't provided <require nnz>");
    }
    nnz = (PetscInt)ia[m];
    
    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"aij_rank%D.dat",rank);
    fp = fopen(name,"w");
    
    PetscFPrintf(PETSC_COMM_SELF,fp,"%D\n",nnz);
    for (i=0; i<nnz; i++) {
      PetscFPrintf(PETSC_COMM_SELF,fp,"%D %+1.8e\n",i,aij[i]);
    }
    fclose(fp);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Default_MatIsSymmetric"
PetscErrorCode Default_MatIsSymmetric(Mat A,PetscReal tol,PetscBool *flg)
{
  PetscErrorCode ierr;
  Vec            x,y1,y2;
  PetscReal      min,max;
  
  PetscFunctionBegin;
  ierr = MatCreateVecs(A,&y2,&y1);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&x,NULL);CHKERRQ(ierr);
  
  ierr = VecSetRandom(x,NULL);CHKERRQ(ierr);
  ierr = MatMult(A,x,y1);CHKERRQ(ierr);
  ierr = MatMultTranspose(A,x,y2);CHKERRQ(ierr);
  
  ierr = VecAXPY(y1,-1.0,y2);CHKERRQ(ierr);
  ierr = VecAbs(y1);CHKERRQ(ierr);
  ierr = VecMin(y1,NULL,&min);CHKERRQ(ierr);
  ierr = VecMax(y1,NULL,&max);CHKERRQ(ierr);
  
  if (max < tol) {
    *flg = PETSC_TRUE;
  } else {
    *flg = PETSC_FALSE;
  }
  
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y1);CHKERRQ(ierr);
  ierr = VecDestroy(&y2);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCWSMP_MatIsSymmetric"
PetscErrorCode PCWSMP_MatIsSymmetric(Mat A,PetscReal tol,PetscBool *flg)
{
  PetscErrorCode ierr;
  PetscMPIInt    size;
  MPI_Comm       comm;
  
  PetscFunctionBegin;
  PetscObjectGetComm((PetscObject)A,&comm);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  if (size == 1) {
    ierr = MatIsSymmetric(A,tol,flg);CHKERRQ(ierr);
  } else {
    ierr = Default_MatIsSymmetric(A,tol,flg);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

/* helpers to get AIJ info from sequential and parallel matrices */
#undef __FUNCT__
#define __FUNCT__ "PCWSMP_ExtractUpperTriangular_MatSeqAIJ"
PetscErrorCode PCWSMP_ExtractUpperTriangular_MatSeqAIJ(Mat parent,Mat A,int *_nnz_ut,int **_ia_ut,int **_ja_ut,PetscBool reuse,double **_vals)
{
  PetscErrorCode ierr;
  const PetscInt *ia;
  const PetscInt *ja;
  PetscInt       m,nnz_i,cnt,idx,nnz_ut,i,j;
  PetscBool      done;
  PetscScalar    *array;
  int            *ia_ut,*ja_ut;
  double         *vals;
  PetscInt       start,end;
  MPI_Comm       comm;
  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(parent,&start,&end);CHKERRQ(ierr);
  
  ierr = MatSeqAIJGetArray(A,&array);CHKERRQ(ierr);
  ierr = MatGetRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
  if (!done) {
    SETERRQ(comm,PETSC_ERR_SUP,"MatGetRowIJ failed... aborting");
  }
  
  /* two pass - a) compute new nnz */
  cnt = 0;
  nnz_ut = 0;
  for (i=0; i<m; i++) {
    nnz_i = ia[i+1]-ia[i];
    for (j=cnt; j<cnt+nnz_i; j++) {
      if (ja[j] >= (i+start)) {
        nnz_ut++;
      }
    }
    cnt = cnt + nnz_i;
  }
  
  /* allocate and store */
  if (_ia_ut) { ierr = PetscMalloc1(m+1,&ia_ut);CHKERRQ(ierr); }
  if (_ja_ut) { ierr = PetscMalloc1(nnz_ut,&ja_ut);CHKERRQ(ierr); }
  if (_vals) {
    vals = *_vals;
    if (!reuse) {
      ierr = PetscMalloc1(nnz_ut,&vals);CHKERRQ(ierr);
    }
  }
  
  /* ia_ut */
  if (_ia_ut) {
    cnt = 0;
    idx = 0;
    for (i=0; i<m; i++) {
      PetscInt nnz_i_ut = 0;
      
      nnz_i = ia[i+1]-ia[i];
      for (j=cnt; j<cnt+nnz_i; j++) {
        if (ja[j] >= (i+start)) {
          nnz_i_ut++;
        }
      }
      ia_ut[i] = (int)idx + 1; /* fortran +1 */
      if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ia_ut[] failed");}
      
      idx = idx + nnz_i_ut;
      cnt = cnt + nnz_i;
    }
    ia_ut[m] = (int)idx + 1; /* fortran +1 */
  }
  
  /* ja_ut */
  if (_ja_ut) {
    idx = 0;
    cnt = 0;
    for (i=0; i<m; i++) {
      nnz_i = ia[i+1]-ia[i];
      for (j=cnt; j<cnt+nnz_i; j++) {
        if (ja[j] >= (i+start)) {
          
          ja_ut[idx] = (int)ja[j] + 1; /* fortran +1 */
          if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ja_ut[] failed");}
          idx++;
        }
      }
      cnt = cnt + nnz_i;
    }
  }
  
  /* aij */
  if (_vals) {
    idx = 0;
    cnt = 0;
    for (i=0; i<m; i++) {
      nnz_i = ia[i+1]-ia[i];
      for (j=cnt; j<cnt+nnz_i; j++) {
        if (ja[j] >= (i+start)) {
          vals[idx] = (double)array[j];
          if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"aij[] failed");}
          idx++;
        }
      }
      cnt = cnt + nnz_i;
    }
  }
  
  ierr = MatRestoreRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
  if (!done) SETERRQ(comm,PETSC_ERR_SUP,"MatRestoreRowIJ failed... aborting");
  ierr = MatSeqAIJRestoreArray(A,&array);CHKERRQ(ierr);
  
  if (_nnz_ut) { *_nnz_ut = (int)nnz_ut; }
  if (_ia_ut) { *_ia_ut  = ia_ut; }
  if (_ja_ut) { *_ja_ut  = ja_ut; }
  if (_vals) { *_vals   = vals; }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "xxx_PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ"
PetscErrorCode xxx_PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(Mat A,int *_nnz_ut,int **_ia_ut,int **_ja_ut)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PCWSMP_ExtractUpperTriangular_MatSeqAIJ(A,A,_nnz_ut,_ia_ut,_ja_ut,PETSC_FALSE,NULL);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "xxx_PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ"
PetscErrorCode xxx_PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(Mat A,int *_nnz_ut,int **_ia_ut,int **_ja_ut)
{
  PetscErrorCode ierr;
  PetscInt       m,n,M,N,start,end;
  IS             irow,icol;
  Mat            *smat;
  
  PetscFunctionBegin;
  ierr = MatGetSize(A,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
  
  ierr = ISCreateStride(PETSC_COMM_SELF,m,start,1,&irow);CHKERRQ(ierr);
  ierr = ISCreateStride(PETSC_COMM_SELF,N,0,1,&icol);CHKERRQ(ierr);
  
  ierr = MatCreateSubMatrices(A,1,&irow,&icol,MAT_INITIAL_MATRIX,&smat);CHKERRQ(ierr);
  
  ierr = PCWSMP_ExtractUpperTriangular_MatSeqAIJ(A,smat[0],_nnz_ut,_ia_ut,_ja_ut,PETSC_FALSE,NULL);CHKERRQ(ierr);
  
  ierr = MatDestroyMatrices(1,&smat);CHKERRQ(ierr);
  ierr = ISDestroy(&irow);CHKERRQ(ierr);
  ierr = ISDestroy(&icol);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "xxx_PCWSMP_ExtractUpperTriangularAIJ"
PetscErrorCode xxx_PCWSMP_ExtractUpperTriangularAIJ(Mat A,PetscBool reuse,int nnz_ut,double **_vals_ut)
{
  PetscErrorCode ierr;
  PetscInt       m,n,M,N,start,end;
  IS             irow,icol;
  Mat            *smat;
  int            nnz;
  
  PetscFunctionBegin;
  ierr = MatGetSize(A,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
  
  ierr = ISCreateStride(PETSC_COMM_SELF,m,start,1,&irow);CHKERRQ(ierr);
  ierr = ISCreateStride(PETSC_COMM_SELF,N,0,1,&icol);CHKERRQ(ierr);
  
  ierr = MatCreateSubMatrices(A,1,&irow,&icol,MAT_INITIAL_MATRIX,&smat);CHKERRQ(ierr);
  
  ierr = PCWSMP_ExtractUpperTriangular_MatSeqAIJ(A,smat[0],&nnz,NULL,NULL,reuse,_vals_ut);CHKERRQ(ierr);
  if (reuse) {
    if (nnz_ut != nnz) SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Number of non-zeros appears to have changed which is not permitted when re-using operator");
  }
  
  ierr = MatDestroyMatrices(1,&smat);CHKERRQ(ierr);
  ierr = ISDestroy(&irow);CHKERRQ(ierr);
  ierr = ISDestroy(&icol);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/* wrappers for WSMP to hide some ugly #if and if sequential type statements */
#undef __FUNCT__
#define __FUNCT__ "call_wsmp"
PetscErrorCode call_wsmp(MPI_Comm comm,PC_WSMP *wsmp)
{
  PetscFunctionBegin;

  if (wsmp->sequential) {
    #ifdef TATIN_HAVE_WSSMP
      wssmp_ ( &wsmp->Nlocal, wsmp->IA, wsmp->JA, wsmp->AVALS, wsmp->DIAG, wsmp->PERM, wsmp->INVP, wsmp->B, &wsmp->LDB, &wsmp->NRHS,
               &wsmp->AUX, &wsmp->NAUX, wsmp->MRP, wsmp->IPARM, wsmp->DPARM );
    
      if (wsmp->IPARM[64 -1] != 0) {
        SETERRQ1(comm,PETSC_ERR_USER,"[wsmp] WSSMP generated the following error code: %d",wsmp->IPARM[64 -1]);
      }
    #else
      #ifdef TATIN_HAVE_PWSSMP
        pwssmp_ ( &wsmp->Nlocal, wsmp->IA, wsmp->JA, wsmp->AVALS, wsmp->DIAG, wsmp->PERM, wsmp->INVP, wsmp->B, &wsmp->LDB, &wsmp->NRHS,
                  &wsmp->AUX, &wsmp->NAUX, wsmp->MRP, wsmp->IPARM, wsmp->DPARM );
        if (wsmp->IPARM[64 -1] != 0) {
          SETERRQ1(comm,PETSC_ERR_USER,"[wsmp] PWSSMP generated the following error code: %d",wsmp->IPARM[64 -1]);
        }
      #else
        SETERRQ(comm,PETSC_ERR_SUP,"[wsmp] Missing external package (WSSMP or PWSSMP) needed for type -pc_type \"wsmp\" when using 1 MPI rank");
      #endif
    #endif
  } else {
    #ifdef TATIN_HAVE_PWSSMP
      pwssmp_ ( &wsmp->Nlocal, wsmp->IA, wsmp->JA, wsmp->AVALS, wsmp->DIAG, wsmp->PERM, wsmp->INVP, wsmp->B, &wsmp->LDB, &wsmp->NRHS,
                &wsmp->AUX, &wsmp->NAUX, wsmp->MRP, wsmp->IPARM, wsmp->DPARM );
    
      if (wsmp->IPARM[64 -1] != 0) {
        if (wsmp->IPARM[64 -1] == -200) {
          PetscPrintf(comm,"[wsmp] PWSSMP generated error code indicating matrix is too small for a parallel solve - suggest using fewer MPI-ranks\n");
        }
        SETERRQ1(comm,PETSC_ERR_USER,"[wsmp] PWSSMP generated the following error code: %d",wsmp->IPARM[64 -1]);
      }
    #else
      SETERRQ(comm,PETSC_ERR_SUP,"[wsmp] Missing external package (PWSSMP) needed for type -pc_type \"wsmp\" when using > 1 MPI rank");
    #endif
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "call_wsffree"
PetscErrorCode call_wsffree(PC_WSMP *wsmp)
{
  PetscFunctionBegin;
  if (wsmp->sequential) {
#ifdef TATIN_HAVE_WSSMP
    wsffree_();
#endif
#ifdef TATIN_HAVE_PWSSMP
    pwsffree_();
#endif
  } else {
#ifdef TATIN_HAVE_PWSSMP
    pwsffree_();
#endif
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "call_wsmp_clear"
PetscErrorCode call_wsmp_clear(PC_WSMP *wsmp)
{
  PetscFunctionBegin;
  if (wsmp->sequential) {
#ifdef TATIN_HAVE_WSSMP
    wsmp_clear_();
#endif
#ifdef TATIN_HAVE_PWSSMP
    pwsmp_clear_();
#endif
  } else {
#ifdef TATIN_HAVE_PWSSMP
    pwsmp_clear_();
#endif
  }
  PetscFunctionReturn(0);
}

/*
 subroutine interface_wsetmpicomm(comm,commout)
 integer ierr,my_rank,size
 include 'mpif.h'
 
 call mpi_comm_dup(comm,fcomm,ierr)
 call wsetmpicomm_(fcomm)
 commout = fcomm
 
 
 #undef __FUNCT__
 #define __FUNCT__ "call_wsmp_wsetmpicomm"
 PetscErrorCode call_wsmp_wsetmpicomm(PC_WSMP *wsmp,MPI_Comm comm)
 {
 integer ierr,my_rank,size
 include 'mpif.h'
 
 call interface_wsetmpicomm(comm,&wsmp->fcomm)
 }
 */

#undef __FUNCT__
#define __FUNCT__ "call_wsmp_wsetmpicomm"
PetscErrorCode call_wsmp_wsetmpicomm(PC_WSMP *wsmp,MPI_Comm comm)
{
#ifdef TATIN_HAVE_PWSSMP
  MPI_Fint fcomm;
  
  PetscFunctionBegin;
  fcomm = MPI_Comm_c2f(comm);
  wsetmpicomm_(&fcomm);
#else
  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"[wsmp] Missing external package needed for type -pc_type \"wsmp\" (PWSSMP)");
#endif
  
  PetscFunctionReturn(0);
}

/* implementations for PC_WSMP */
#undef __FUNCT__
#define __FUNCT__ "PCSetUp_WSMP"
static PetscErrorCode PCSetUp_WSMP(PC pc)
{
  PetscErrorCode ierr;
  PC_WSMP        *wsmp = (PC_WSMP*)pc->data;
  Mat            A,B;
  PetscLogDouble t0,t1;
  MPI_Comm       comm;
  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
  ierr = PCGetOperators(pc,&A,&B);CHKERRQ(ierr);
  
  /* construction phase */
  if (!pc->setupcalled) {
    PetscInt M,N,m,n;
    
    /* Determine if matrix is symmetric */
    ierr = PCWSMP_MatIsSymmetric(B,1.0e-8,&wsmp->symmetric);CHKERRQ(ierr);
    if (!wsmp->symmetric) SETERRQ(comm,PETSC_ERR_SUP,"[wsmp] Only support for symmetric WSMP implemenatations exists");
    
    ierr = MatGetSize(B,&M,&N);CHKERRQ(ierr);
    ierr = MatGetLocalSize(B,&m,&n);CHKERRQ(ierr);
    
    /* get local size */
    if (m != n) {
      SETERRQ(comm,PETSC_ERR_SUP,"[wsmp] WSMP only supports square matrices");
    }
    
    wsmp->Nglobal = (int)M;
    wsmp->Nlocal = (int)m;
    wsmp->LDB = wsmp->Nlocal;
    wsmp->NRHS = 1;
    wsmp->NAUX = 0;
    
    /* allocate other local variables */
    wsmp->DIAG = NULL;
    ierr = PetscMalloc1(wsmp->Nglobal,&wsmp->PERM);CHKERRQ(ierr); ierr = PetscMemzero(wsmp->PERM,sizeof(int)   *wsmp->Nglobal);CHKERRQ(ierr);
    ierr = PetscMalloc1(wsmp->Nglobal,&wsmp->INVP);CHKERRQ(ierr); ierr = PetscMemzero(wsmp->INVP,sizeof(int)   *wsmp->Nglobal);CHKERRQ(ierr);
    ierr = PetscMalloc1(wsmp->Nlocal,&wsmp->B);CHKERRQ(ierr);     ierr = PetscMemzero(wsmp->B,   sizeof(double)*wsmp->Nlocal);CHKERRQ(ierr);
    ierr = PetscMalloc1(wsmp->Nlocal,&wsmp->MRP);CHKERRQ(ierr);   ierr = PetscMemzero(wsmp->MRP, sizeof(int)   *wsmp->Nlocal);CHKERRQ(ierr);
    
    /* Fetch ia, ja from matrix */
    PetscTime(&t0);
    if (wsmp->sequential) {
      ierr = xxx_PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
    } else {
      ierr = xxx_PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
    }
    PetscTime(&t1);
    PetscInfo1(pc,"CSR extraction --- %1.2e sec\n",t1-t0);
    
    {
      PetscBool view = PETSC_FALSE;
      PetscOptionsGetBool(NULL,NULL,"-pc_wsmp_debug",&view,NULL);
      if (view) {
        ierr = PCWSMP_CSRView(B,wsmp->IA,wsmp->JA,NULL);CHKERRQ(ierr);
      }
    }
    
    /* -- [wsmp] : ordering -- */
    PetscTime(&t0);
    wsmp->IPARM[2 -1] = 1;
    wsmp->IPARM[3 -1] = 1;
    ierr = WSMPSetFromOptions_Ordering(pc,wsmp);CHKERRQ(ierr);
    ierr = call_wsmp(comm,wsmp);CHKERRQ(ierr);
    PetscTime(&t1);
    PetscInfo1(pc,"Ordering --- %1.2e sec\n",t1-t0);
    
    /* -- [wsmp] : symbolic factorization -- */
    wsmp->IPARM[2 -1] = 2;
    wsmp->IPARM[3 -1] = 2;
    PetscTime(&t0);
    ierr = WSMPSetFromOptions_SymbolicFactorization(pc,wsmp);CHKERRQ(ierr);
    ierr = call_wsmp(comm,wsmp);CHKERRQ(ierr);
    PetscTime(&t1);
    PetscInfo1(pc,"Symbolic factorization --- %1.2e sec\n",t1-t0);
  } else {

    /* Determine if matrix is symmetric (either due to non-zero pattern or changing aij values */
    ierr = PCWSMP_MatIsSymmetric(B,1.0e-8,&wsmp->symmetric);CHKERRQ(ierr);
    if (!wsmp->symmetric) SETERRQ(comm,PETSC_ERR_SUP,"[wsmp] Only support for symmetric WSMP implemenatations exists");

    if (pc->flag != SAME_NONZERO_PATTERN) {
      /* Repeated solve but non-zero structure has changed, use new data */
      /* Re-construct ordering and symbolic factorization */
      wsmp->nnz = 0;
      ierr = PetscFree(wsmp->IA);CHKERRQ(ierr);
      ierr = PetscFree(wsmp->JA);CHKERRQ(ierr);
      
      PetscTime(&t0);
      if (wsmp->sequential) {
        ierr = xxx_PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
      } else {
        ierr = xxx_PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
      }
      PetscTime(&t1);
      PetscInfo1(pc,"CSR extraction --- %1.2e sec\n",t1-t0);
      
      /* -- [wsmp] : ordering -- */
      PetscTime(&t0);
      wsmp->IPARM[2 -1] = 1;
      wsmp->IPARM[3 -1] = 1;
      ierr = call_wsmp(comm,wsmp);CHKERRQ(ierr);
      PetscTime(&t1);
      PetscInfo1(pc,"Symmetric ordering --- %1.2e sec\n",t1-t0);
      
      /* -- [wsmp] : symbolic factorization -- */
      wsmp->IPARM[2 -1] = 2;
      wsmp->IPARM[3 -1] = 2;
      PetscTime(&t0);
      ierr = call_wsmp(comm,wsmp);CHKERRQ(ierr);
      PetscTime(&t1);
      PetscInfo1(pc,"Symbolic factorization --- %1.2e sec\n",t1-t0);
    }
  }
  
  /* Fetch matrix entries */
  PetscTime(&t0);
  if (!pc->setupcalled) {
    /* If first time we are in PCSetUp, use new data */
    ierr = xxx_PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_FALSE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
  } else {
    if (pc->flag != SAME_NONZERO_PATTERN) {
      /* Repeated solve but non-zero structure has changed, use new data */
      ierr = PetscFree(wsmp->AVALS);CHKERRQ(ierr);
      /* release internal memory for factors */
      ierr = call_wsffree(wsmp);CHKERRQ(ierr);
      
      ierr = xxx_PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_FALSE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
    } else {
      /* Repeated solve but non-zero structure is the same, re-use data */
      ierr = xxx_PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_TRUE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
    }
  }
  PetscTime(&t1);
  PetscInfo1(pc,"CSR extraction --- %1.2e sec\n",t1-t0);
  
  {
    PetscBool view = PETSC_FALSE;
    PetscOptionsGetBool(NULL,NULL,"-pc_wsmp_debug",&view,NULL);
    if (view) {
      ierr = PCWSMP_CSRView(B,NULL,NULL,wsmp->AVALS);CHKERRQ(ierr);
    }
  }
  {
    PetscBool check = PETSC_FALSE;
    PetscOptionsGetBool(NULL,NULL,"-pc_wsmp_check_csr",&check,NULL);
    if (check) {
      ierr = PCWSMP_CheckCSR(B,wsmp);CHKERRQ(ierr);
    }
  }
  
  /* -- [wsmp] : numeric factorization - Cholesky -- */
  PetscTime(&t0);
  wsmp->IPARM[2 -1] = 3;
  wsmp->IPARM[3 -1] = 3;
  
  ierr = WSMPSetFromOptions_NumericFactorization(pc,wsmp);CHKERRQ(ierr);
  ierr = call_wsmp(comm,wsmp);CHKERRQ(ierr);
  PetscTime(&t1);
  PetscInfo1(pc,"Numerical factorization --- %1.2e sec\n",t1-t0);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_WSMP"
static PetscErrorCode PCApply_WSMP(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;
  PC_WSMP        *wsmp = (PC_WSMP*)pc->data;
  const PetscScalar    *_val;
  PetscScalar    *_val_y;
  PetscInt       i,m;
  static int     beenhere = 0;
  PetscBool      view = PETSC_FALSE;
  int            rank;
  MPI_Comm       comm;
  
  PetscFunctionBegin;
  PetscOptionsGetBool(NULL,NULL,"-pc_wsmp_debug",&view,NULL);
  
  ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  ierr = VecGetLocalSize(x,&m);CHKERRQ(ierr);
  ierr = VecGetArrayRead(x,&_val);CHKERRQ(ierr);
  for (i=0; i<m; i++) {
    wsmp->B[i] = PetscRealPart(_val[i]);
  }
  ierr = VecRestoreArrayRead(x,&_val);CHKERRQ(ierr);
  
  if (view && (beenhere==0)) {
    ierr = PCWSMP_VecView("rhs",wsmp);CHKERRQ(ierr);
  }
  
  /* apply L^T L factors  - Back substitution */
  wsmp->IPARM[2 -1] = 4;
  wsmp->IPARM[3 -1] = 4;
  ierr = WSMPSetFromOptions_BackSubstitution(pc,wsmp);CHKERRQ(ierr);
  ierr = call_wsmp(comm,wsmp);CHKERRQ(ierr);
  
  /* iterative refinement */
  wsmp->IPARM[2 -1] = 5;
  wsmp->IPARM[3 -1] = 5;
  ierr = WSMPSetFromOptions_IterativeRefinement(pc,wsmp);CHKERRQ(ierr);
  ierr = call_wsmp(comm,wsmp);CHKERRQ(ierr);
  
  if (view && (beenhere==0)) {
    ierr = PCWSMP_VecView("sol",wsmp);CHKERRQ(ierr);
  }
  
  ierr = VecGetArray(y,&_val_y);CHKERRQ(ierr);
  for (i=0; i<m; i++) {
    _val_y[i] = (PetscScalar)wsmp->B[i];
  }
  ierr = VecRestoreArray(y,&_val_y);CHKERRQ(ierr);
  
  beenhere = 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCReset_WSMP"
static PetscErrorCode PCReset_WSMP(PC pc)
{
  PetscErrorCode   ierr;
  PC_WSMP          *wsmp = (PC_WSMP*)pc->data;
  
  PetscFunctionBegin;
  ierr = PetscFree(wsmp->IA);CHKERRQ(ierr);
  ierr = PetscFree(wsmp->JA);CHKERRQ(ierr);
  ierr = PetscFree(wsmp->AVALS);CHKERRQ(ierr);
  if (wsmp->DIAG) { ierr = PetscFree(wsmp->DIAG);CHKERRQ(ierr); }
  ierr = PetscFree(wsmp->PERM);CHKERRQ(ierr);
  ierr = PetscFree(wsmp->INVP);CHKERRQ(ierr);
  ierr = PetscFree(wsmp->B);CHKERRQ(ierr);
  ierr = PetscFree(wsmp->MRP);CHKERRQ(ierr);
  
  ierr = call_wsmp_clear(wsmp);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_WSMP"
static PetscErrorCode PCDestroy_WSMP(PC pc)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PCReset_WSMP(pc);CHKERRQ(ierr);
  
  ierr = PetscFree(pc->data);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WSMPSetFromOptions_Ordering"
PetscErrorCode WSMPSetFromOptions_Ordering(PC pc,PC_WSMP *wsmp)
{
  PetscErrorCode ierr;
  PetscBool      found;
  PetscInt       index,ival;
  PetscReal      dval;
  const char     *prefix;
  
  PetscFunctionBegin;
  prefix = ((PetscObject)pc)->prefix;
  
  /* recommended for finite elements */
  if (wsmp->sequential) {
    wsmp->IPARM[16 -1] = 1;
    wsmp->IPARM[17 -1] = 0;
    wsmp->IPARM[18 -1] = 1;
    wsmp->IPARM[19 -1] = 0;
    wsmp->IPARM[20 -1] = 2;
  } else {
    wsmp->IPARM[16 -1] = 1;
    wsmp->IPARM[17 -1] = 0;
    wsmp->IPARM[18 -1] = 0;
    wsmp->IPARM[19 -1] = 0;
    wsmp->IPARM[20 -1] = 2;
  }
  
  index = 10;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm10",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 15;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm15",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 16;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm16",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 17;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm17",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 18;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm18",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 19;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm19",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 20;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm20",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 27;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm27",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 31;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm31",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 10;  dval = (PetscInt)wsmp->DPARM[index-1];
  ierr = PetscOptionsGetReal(NULL,prefix,"-pc_wsmp_dparm10",&dval,&found);CHKERRQ(ierr); if (found) { wsmp->DPARM[index-1] = (double)dval; }
  
  index = 11;  dval = (PetscInt)wsmp->DPARM[index-1];
  ierr = PetscOptionsGetReal(NULL,prefix,"-pc_wsmp_dparm11",&dval,&found);CHKERRQ(ierr); if (found) { wsmp->DPARM[index-1] = (double)dval; }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WSMPSetFromOptions_SymbolicFactorization"
PetscErrorCode WSMPSetFromOptions_SymbolicFactorization(PC pc,PC_WSMP *wsmp)
{
  PetscErrorCode ierr;
  PetscInt       ival,index;
  PetscBool      found;
  const char     *prefix;
  
  PetscFunctionBegin;
  prefix = ((PetscObject)pc)->prefix;
  index = 35;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm35",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WSMPSetFromOptions_BackSubstitution"
PetscErrorCode WSMPSetFromOptions_BackSubstitution(PC pc,PC_WSMP *wsmp)
{
  PetscErrorCode ierr;
  PetscInt       ival,index;
  PetscBool      found;
  const char     *prefix;
  
  PetscFunctionBegin;
  prefix = ((PetscObject)pc)->prefix;
  index = 30;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm30",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WSMPSetFromOptions_IterativeRefinement"
PetscErrorCode WSMPSetFromOptions_IterativeRefinement(PC pc,PC_WSMP *wsmp)
{
  PetscErrorCode ierr;
  PetscInt       ival,index;
  PetscBool      found;
  PetscReal      dval;
  const char     *prefix;
  
  PetscFunctionBegin;
  prefix = ((PetscObject)pc)->prefix;
  index = 6;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm6",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 7;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm7",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 36;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm36",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 6;  dval = (PetscInt)wsmp->DPARM[index-1];
  ierr = PetscOptionsGetReal(NULL,prefix,"-pc_wsmp_dparm6",&dval,&found);CHKERRQ(ierr); if (found) { wsmp->DPARM[index-1] = (double)dval; }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WSMPSetFromOptions_NumericFactorization"
PetscErrorCode WSMPSetFromOptions_NumericFactorization(PC pc,PC_WSMP *wsmp)
{
  PetscErrorCode ierr;
  PetscInt       ival,index;
  PetscBool      found;
  PetscReal      dval;
  const char     *prefix;
  
  PetscFunctionBegin;
  prefix = ((PetscObject)pc)->prefix;
  index = 8;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm8",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 11;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm11",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 12;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm12",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 13;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm13",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 14;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm14",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 25;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm25",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 26;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm26",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 34;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm34",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 6;  dval = (PetscInt)wsmp->DPARM[index-1];
  ierr = PetscOptionsGetReal(NULL,prefix,"-pc_wsmp_dparm6",&dval,&found);CHKERRQ(ierr); if (found) { wsmp->DPARM[index-1] = (double)dval; }
  
  index = 12;  dval = (PetscInt)wsmp->DPARM[index-1];
  ierr = PetscOptionsGetReal(NULL,prefix,"-pc_wsmp_dparm12",&dval,&found);CHKERRQ(ierr); if (found) { wsmp->DPARM[index-1] = (double)dval; }
  
  index = 15;  dval = (PetscInt)wsmp->DPARM[index-1];
  ierr = PetscOptionsGetReal(NULL,prefix,"-pc_wsmp_dparm15",&dval,&found);CHKERRQ(ierr); if (found) { wsmp->DPARM[index-1] = (double)dval; }
  
  index = 21;  dval = (PetscInt)wsmp->DPARM[index-1];
  ierr = PetscOptionsGetReal(NULL,prefix,"-pc_wsmp_dparm21",&dval,&found);CHKERRQ(ierr); if (found) { wsmp->DPARM[index-1] = (double)dval; }
  
  index = 22;  dval = (PetscInt)wsmp->DPARM[index-1];
  ierr = PetscOptionsGetReal(NULL,prefix,"-pc_wsmp_dparm22",&dval,&found);CHKERRQ(ierr); if (found) { wsmp->DPARM[index-1] = (double)dval; }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetFromOptions_WSMP"
static PetscErrorCode PCSetFromOptions_WSMP(PetscOptionItems *PetscOptionsObject,PC pc)
{
  PetscErrorCode ierr;
  PC_WSMP        *wsmp = (PC_WSMP*)pc->data;
  PetscInt       ival,index;
  PetscBool      found;
  
  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"WSMP options [init]");CHKERRQ(ierr);
  
  index = 8;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsInt("-pc_wsmp_iparm8","","PCWSMPSetIParm",ival,&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 9;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsInt("-pc_wsmp_iparm9","","PCWSMPSetIParm",ival,&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 13;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsInt("-pc_wsmp_iparm13","","PCWSMPSetIParm",ival,&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 28;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsInt("-pc_wsmp_iparm28","","PCWSMPSetIParm",ival,&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 29;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsInt("-pc_wsmp_iparm29","","PCWSMPSetIParm",ival,&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  index = 32;  ival = (PetscInt)wsmp->IPARM[index-1];
  ierr = PetscOptionsInt("-pc_wsmp_iparm32","","PCWSMPSetIParm",ival,&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
  
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCView_WSMP"
static PetscErrorCode PCView_WSMP(PC pc,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PC_WSMP        *wsmp = (PC_WSMP*)pc->data;
  PetscBool      iascii,isstring;
  
  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERSTRING,&isstring);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  WSMP preconditioner\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(4)  = %d <matrix format>\n",wsmp->IPARM[4 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(5)  = %d <C or ForTran numbering>\n",wsmp->IPARM[5 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(6)  = %d <number of iterative refinement steps performed>\n",wsmp->IPARM[6 -1]);
    if ( wsmp->IPARM[7 -1] <= 3) {
      ierr = PetscViewerASCIIPrintf(viewer,"    iparm(7)  = %d <double precision>\n",wsmp->IPARM[7-1]);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"    iparm(7)  = %d <quadruple precision>\n",wsmp->IPARM[7-1]);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(8)  = %d <matrix permutation option>\n",wsmp->IPARM[8 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(9)  = %d <rhs permutation option>\n",wsmp->IPARM[9 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(10) = %d <matrix scaling option>\n",wsmp->IPARM[10 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(11) = %d <pivot type>\n",wsmp->IPARM[11 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(12) = %d <pivot activation flag>\n",wsmp->IPARM[12 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(13) = %d <num. diagonal entries perturbed>\n",wsmp->IPARM[13 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(14) = %d <memory reuse flag>\n",wsmp->IPARM[14 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(15) = %d <N1 cols. factored before remainder>\n",wsmp->IPARM[15 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(16) = %d <ordering type for perm. vectors>\n",wsmp->IPARM[16 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(17) = %d <max. num. nodes in subgraph>\n",wsmp->IPARM[17 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(18) = %d <force minimum local fill>\n",wsmp->IPARM[18 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(19) = %d <random seed for graph permutation>\n",wsmp->IPARM[19 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(20) = %d <matrix characteristics flag>\n",wsmp->IPARM[20 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(21) = %d <N (input matrix size) - M (rank)>\n",wsmp->IPARM[21 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(22) = %d <num. negative eigenvalues>\n",wsmp->IPARM[22 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(23) = %d <double words for factorisation>\n",wsmp->IPARM[23 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(24) = %d <nnz in triangular factor>\n",wsmp->IPARM[24 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(25) = %d <condition num. flag>\n",wsmp->IPARM[25 -1]);
    if (!wsmp->sequential) {
      ierr = PetscViewerASCIIPrintf(viewer,"    iparm(26) = %d <block size for dense matrix computations>\n",wsmp->IPARM[26 -1]);
      ierr = PetscViewerASCIIPrintf(viewer,"    iparm(27) = %d <load balance option>\n",wsmp->IPARM[27 -1]);
      ierr = PetscViewerASCIIPrintf(viewer,"    iparm(28) = %d <multiple RHS partitioning parameter>\n",wsmp->IPARM[28 -1]);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(29) = %d <garbage collection flag>\n",wsmp->IPARM[29 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(30) = %d <forward/backward solve flag>\n",wsmp->IPARM[30 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(31) = %d <factorization type (L.L^T=0, L.D.L^T=1 flag>\n",wsmp->IPARM[31 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(32) = %d <diagonal request flag>\n",wsmp->IPARM[32 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(33) = %d <num cpus used by process in SMP mode>\n",wsmp->IPARM[33 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(34) = %d <num of rows containing zeros>\n",wsmp->IPARM[34 -1]);
    if (!wsmp->sequential) {
      ierr = PetscViewerASCIIPrintf(viewer,"    iparm(35) = %d <distribute Cholesky work flag>\n",wsmp->IPARM[35 -1]);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"    iparm(36) = %d <out of core (OOC) flag>\n",wsmp->IPARM[36 -1]);
    
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(4)  = %1.4e <largest diagonal entry>\n",wsmp->DPARM[4 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(5)  = %1.4e <smallest diagonal entry>\n",wsmp->DPARM[5 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(6)  = %1.4e <relative tol for iterative refinement>\n",wsmp->DPARM[6 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(7)  = %1.4e <relative error for iterative refinement>\n",wsmp->DPARM[7 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(10) = %1.4e <lower threshold on diag without pivoting>\n",wsmp->DPARM[10 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(11) = %1.4e <Bunch-Kaufman pivoting threshold>\n",wsmp->DPARM[11 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(12) = %1.4e <pivoting control>\n",wsmp->DPARM[12 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(13) = %1.4e <num. supernodes>\n",wsmp->DPARM[13 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(21) = %1.4e <replacement pivot value if pivot below lower threshold : dparm(10)>\n",wsmp->DPARM[21 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(22) = %1.4e <replacement pivot value for pivots > dparm(10) but < dparm(12)>\n",wsmp->DPARM[22 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(23) = %1.4e <num. floating point operations>\n",wsmp->DPARM[23 -1]);
    ierr = PetscViewerASCIIPrintf(viewer,"    dparm(24) = %1.4e <expected num. floating point operations>\n",wsmp->DPARM[24 -1]);
    
  } else if (isstring) {
    ierr = PetscViewerStringSPrintf(viewer," WSMP preconditioner");CHKERRQ(ierr);
  } else {
    SETERRQ1(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"Viewer type %s not supported for PC WSMP",((PetscObject)viewer)->type_name);
  }
  
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCCreate_WSMP"
PetscErrorCode PCCreate_WSMP(PC pc)
{
  PetscErrorCode ierr;
  PC_WSMP        *wsmp;
  PetscMPIInt    size;
  
  PetscFunctionBegin;
  ierr = PetscNewLog(pc,&wsmp);CHKERRQ(ierr);
  pc->data = (void*)wsmp;
  
  /* determine if sequential call or parallel call required */
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pc),&size);CHKERRQ(ierr);
  if (size == 1) {
    wsmp->sequential = PETSC_TRUE;
  } else {
    wsmp->sequential = PETSC_FALSE;
  }
  //{
  //    int mt = 1;
  //    wsetmaxthrds_(&mt);
  //}
  /* -- [wsmp] : initialize with default parameters -- */
  wsmp->IPARM[1 -1] = 0;
  wsmp->IPARM[2 -1] = 0;
  wsmp->IPARM[3 -1] = 0;
  ierr = call_wsmp(PetscObjectComm((PetscObject)pc),wsmp);CHKERRQ(ierr);
  wsmp->IPARM[4 -1] = 0; /* CSR/CSC matrix format */
  wsmp->IPARM[5 -1] = 1; /* F style numbering */
  
  pc->ops->apply           = PCApply_WSMP;
  pc->ops->applytranspose  = 0;
  pc->ops->setup           = PCSetUp_WSMP;
  pc->ops->destroy         = PCDestroy_WSMP;
  pc->ops->reset           = PCReset_WSMP;
  pc->ops->setfromoptions  = PCSetFromOptions_WSMP;
  pc->ops->view            = PCView_WSMP;
  
  PetscFunctionReturn(0);
}
EXTERN_C_END

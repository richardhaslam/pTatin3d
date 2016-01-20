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

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petsc/private/pcimpl.h>     /*I "petscpc.h" I*/
#include <petscksp.h>

#define USE_COMPACT_FORM

//#define HAVE_WSSMP
extern void wssmp_ ( int *n, int ia[], int ja[], double a[], double diag[], int perm[], int iperm[], double b[], int *ldb, int *nrhs, int aux[], int *naux, int mrp[], int iparam[], double dparam[] );
extern void wsmp_clear_(void);
extern void wsffree_(void);

//#define HAVE_PWSSMP
extern void pwssmp_ ( int *n, int ia[], int ja[], double a[], double diag[], int perm[], int iperm[], double b[], int *ldb, int *nrhs, int aux[], int *naux, int mrp[], int iparam[], double dparam[] );
extern void pwsmp_clear_(void);
extern void wsetmpicomm_(int *fcomm);
extern void pwsffree_(void);


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
    int    AUX,NAUX; /* depreciated ununsed crap */
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
    PetscMalloc(sizeof(PetscInt)*(_Nlocal+1),&_ia);
    PetscMalloc(sizeof(PetscInt)*_nnz,&_ja);
    PetscMalloc(sizeof(PetscScalar)*_nnz,&_avals);
    
    for (i=0; i<_Nlocal+1; i++) {
        _ia[i] = (PetscInt)wsmp->IA[i];
    }
    for (i=0; i<_nnz; i++) {
        _ja[i] = (PetscInt)wsmp->JA[i];
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
    
    PetscFree(_ia);
    PetscFree(_ja);
    PetscFree(_avals);
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
    
    // ia
    if (ia) {
        PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"ia_rank%D.dat",rank);
        fp = fopen(name,"w");
        
        //PetscFPrintf(PETSC_COMM_SELF,fp,"# idx ia[] \n");
        PetscFPrintf(PETSC_COMM_SELF,fp,"%D\n",m+1);
        for (i=0; i<m; i++) {
            PetscFPrintf(PETSC_COMM_SELF,fp,"%D %D\n",i,ia[i]);
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"%D %D\n",m,ia[m]);
        fclose(fp);
    }
    
    // ja
    if (ja) {
        if (!ia) {
            SETERRQ(comm,PETSC_ERR_USER,"Cannot write ja[] as ia[] isn't provided <require nnz");
        }
        nnz = (PetscInt)ia[m];
        
        PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"ja_rank%D.dat",rank);
        fp = fopen(name,"w");
        
        //PetscFPrintf(PETSC_COMM_SELF,fp,"# idx ia[] \n");
        PetscFPrintf(PETSC_COMM_SELF,fp,"%D\n",nnz);
        for (i=0; i<nnz; i++) {
            PetscFPrintf(PETSC_COMM_SELF,fp,"%D %D\n",i,ja[i]);
        }
        fclose(fp);
    }
    
    // aij
    if (aij) {
        if (!ia) {
            SETERRQ(comm,PETSC_ERR_USER,"Cannot write aij[] as ia[] isn't provided <require nnz>");
        }
        nnz = (PetscInt)ia[m];
        
        PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"aij_rank%D.dat",rank);
        fp = fopen(name,"w");
        
        //PetscFPrintf(PETSC_COMM_SELF,fp,"# idx aij[] \n");
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
    PetscReal      drop_tol = 1.0e-12;
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
                PetscReal abs_value;
                
                abs_value = PetscAbsReal(array[j]);
                if (abs_value > drop_tol) {
                    nnz_ut++;
                }
            }
        }
        cnt = cnt + nnz_i;
    }
    //printf("nnz = %d \n",nnz_ut);
    
    /* allocate and store */
    if (_ia_ut) { PetscMalloc(sizeof(int)*(m+1),&ia_ut); }
    if (_ja_ut) { PetscMalloc(sizeof(int)*(nnz_ut),&ja_ut); }
    if (_vals) {
        vals = *_vals;
        if (!reuse) {
            PetscMalloc(sizeof(double)*(nnz_ut),&vals);
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
                    PetscReal abs_value;
                    
                    abs_value = PetscAbsReal(array[j]);
                    if (abs_value > drop_tol) {
                        nnz_i_ut++;
                    }
                }
            }
            ia_ut[i] = (int)idx + 0; /* fortran +1 */
            if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ia_ut[] failed");}
            
            idx = idx + nnz_i_ut;
            cnt = cnt + nnz_i;
        }
        ia_ut[m] = (int)idx + 0; /* fortran +1 */
    }
    
    /* ja_ut */
    if (_ja_ut) {
        idx = 0;
        cnt = 0;
        for (i=0; i<m; i++) {
            
            nnz_i = ia[i+1]-ia[i];
            for (j=cnt; j<cnt+nnz_i; j++) {
                if (ja[j] >= (i+start)) {
                    PetscReal abs_value;
                    
                    abs_value = PetscAbsReal(array[j]);
                    if (abs_value > drop_tol) {
                        
                        ja_ut[idx] = (int)ja[j] + 0; /* fortran +1 */
                        if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ja_ut[] failed");}
                        idx++;
                    }
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
                    PetscReal value,abs_value;
                    
                    value     = PetscRealPart(array[j]);
                    abs_value = PetscAbsReal(array[j]);
                    if (abs_value > drop_tol) {
                        vals[idx] = (double)array[j];
                        if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"aij[] failed");}
                        idx++;
                    }
                }
            }
            cnt = cnt + nnz_i;
        }
    }
    
    ierr = MatRestoreRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
    if (!done) {
        SETERRQ(comm,PETSC_ERR_SUP,"MatRestoreRowIJ failed... aborting");
    }
    ierr = MatSeqAIJRestoreArray(A,&array);CHKERRQ(ierr);
    
    if (_nnz_ut) { *_nnz_ut = (int)nnz_ut; }
    if (_ia_ut) { *_ia_ut  = ia_ut; }
    if (_ja_ut) { *_ja_ut  = ja_ut; }
    if (_vals) { *_vals   = vals; }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ"
PetscErrorCode PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(Mat A,int *_nnz_ut,int **_ia_ut,int **_ja_ut)
{
    PetscErrorCode ierr;
    const PetscInt *ia;
    const PetscInt *ja;
    PetscInt       m,nnz_i,cnt,idx,nnz_ut,i,j;
    PetscBool      done;
    int            *ia_ut,*ja_ut;
    MPI_Comm       comm;
    
    PetscFunctionBegin;
    ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
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
            if (ja[j] >= i) {
                nnz_ut++;
            }
        }
        cnt = cnt + nnz_i;
    }
    
    /* allocate and store */
    PetscMalloc(sizeof(int)*(m+1),&ia_ut);
    PetscMalloc(sizeof(int)*(nnz_ut),&ja_ut);
    
    /* ia_ut */
    cnt = 0;
    idx = 0;
    for (i=0; i<m; i++) {
        PetscInt nnz_i_ut = 0;
        
        nnz_i = ia[i+1]-ia[i];
        for (j=cnt; j<cnt+nnz_i; j++) {
            if (ja[j] >= i) {
                nnz_i_ut++;
            }
        }
        ia_ut[i] = (int)idx + 0; /* fortran +1 */
        if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ia_ut[] failed");}
        
        idx = idx + nnz_i_ut;
        cnt = cnt + nnz_i;
    }
    ia_ut[m] = (int)idx + 0; /* fortran +1 */
    if (idx != (nnz_ut)) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ia_ut[m] failed");}
    
    /* ja_ut */
    idx = 0;
    cnt = 0;
    for (i=0; i<m; i++) {
        
        nnz_i = ia[i+1]-ia[i];
        for (j=cnt; j<cnt+nnz_i; j++) {
            if (ja[j] >= i) {
                ja_ut[idx] = (int)ja[j] + 0; /* fortran +1 */
                if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ja_ut[] failed");}
                idx++;
            }
        }
        cnt = cnt + nnz_i;
    }
    
    ierr = MatRestoreRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
    if (!done) {
        SETERRQ(comm,PETSC_ERR_SUP,"MatRestoreRowIJ failed... aborting");
    }
    
    *_nnz_ut = (int)nnz_ut;
    *_ia_ut  = ia_ut;
    *_ja_ut  = ja_ut;
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ"
PetscErrorCode PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(Mat A,int *_nnz_ut,int **_ia_ut,int **_ja_ut)
{
    PetscErrorCode ierr;
    const PetscInt *ia;
    const PetscInt *ja;
    PetscInt       m,n,M,N,nnz_i,cnt,idx,nnz_ut,i,j,start,end;
    PetscBool      done;
    int            *ia_ut,*ja_ut;
    IS             irow,icol;
    Mat            *smat;
    MPI_Comm       comm;
    
    PetscFunctionBegin;
    ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
    ierr = MatGetSize(A,&M,&N);CHKERRQ(ierr);
    ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
    
    ierr = ISCreateStride(PETSC_COMM_SELF,m,start,1,&irow);CHKERRQ(ierr);
    ierr = ISCreateStride(PETSC_COMM_SELF,N,0,1,&icol);CHKERRQ(ierr);
    
    /* Use MatGetSubMatrices() to get a sequential matrix */
    ierr = MatGetSubMatrices(A,1,&irow,&icol,MAT_INITIAL_MATRIX,&smat);CHKERRQ(ierr);
    
    ierr = MatGetRowIJ(smat[0],0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
    if (!done) {
        SETERRQ(comm,PETSC_ERR_SUP,"MatGetRowIJ failed... aborting");
    }
    
    /* two pass - a) compute new nnz */
    cnt = 0;
    nnz_ut = 0;
    for (i=0; i<m; i++) {
        nnz_i = ia[i+1]-ia[i];
        for (j=cnt; j<cnt+nnz_i; j++) {
            if (ja[j] >= i+start) {
                nnz_ut++;
            }
        }
        cnt = cnt + nnz_i;
    }
    
    /* allocate and store */
    PetscMalloc(sizeof(int)*(m+1),&ia_ut);
    PetscMalloc(sizeof(int)*(nnz_ut),&ja_ut);
    
    /* ia_ut */
    cnt = 0;
    idx = 0;
    for (i=0; i<m; i++) {
        PetscInt nnz_i_ut = 0;
        
        nnz_i = ia[i+1]-ia[i];
        for (j=cnt; j<cnt+nnz_i; j++) {
            if (ja[j] >= i+start) {
                nnz_i_ut++;
            }
        }
        ia_ut[i] = (int)idx;
        if (idx >= nnz_ut) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ia_ut[] failed");}
        
        idx = idx + nnz_i_ut;
        cnt = cnt + nnz_i;
    }
    ia_ut[m] = (int)idx;
    if (idx != (nnz_ut)) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"ia_ut[m] failed");}
    
    
    /* ja_ut */
    idx = 0;
    cnt = 0;
    for (i=0; i<m; i++) {
        
        nnz_i = ia[i+1]-ia[i];
        for (j=cnt; j<cnt+nnz_i; j++) {
            if (ja[j] >= i+start) {
                ja_ut[idx] = (int)ja[j];
                idx++;
            }
        }
        cnt = cnt + nnz_i;
    }
    ierr = MatRestoreRowIJ(smat[0],0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
    if (!done) {
        SETERRQ(comm,PETSC_ERR_SUP,"MatRestoreRowIJ failed... aborting");
    }
    
    ierr = MatDestroyMatrices(1,&smat);CHKERRQ(ierr);
    ierr = ISDestroy(&irow);CHKERRQ(ierr);
    ierr = ISDestroy(&icol);CHKERRQ(ierr);
    
    *_nnz_ut = (int)nnz_ut;
    *_ia_ut  = ia_ut;
    *_ja_ut  = ja_ut;
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCWSMP_ExtractUpperTriangularAIJ"
PetscErrorCode PCWSMP_ExtractUpperTriangularAIJ(Mat A,PetscBool reuse,int nnz_ut,double **_vals_ut)
{
    PetscErrorCode    ierr;
    PetscInt          m,n,i,j,start,end;
    PetscInt          ncols,idx;
    const PetscInt    *cols;
    const PetscScalar *vals;
    double            *vals_ut;
    
    PetscFunctionBegin;
    if (!reuse) {
        PetscMalloc(sizeof(double)*nnz_ut,&vals_ut);
    } else {
        vals_ut = *_vals_ut;
    }
    
    ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
    ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
    
    idx = 0;
    for (i=start; i<start+m; i++) {
        ierr = MatGetRow(A,i,&ncols,&cols,&vals);CHKERRQ(ierr);
        
        for (j=0; j<ncols; j++) {
            if (cols[j] >= i) {
                vals_ut[idx] = PetscRealPart(vals[j]);
                idx++;
            }
        }
        ierr = MatRestoreRow(A,i,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
    if (!reuse) {
        *_vals_ut = vals_ut;
    }
    
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
    
    /* Use MatGetSubMatrices() to get a sequential matrix */
    ierr = MatGetSubMatrices(A,1,&irow,&icol,MAT_INITIAL_MATRIX,&smat);CHKERRQ(ierr);
    
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
    
    PetscFunctionBegin;
    ierr = MatGetSize(A,&M,&N);CHKERRQ(ierr);
    ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
    
    ierr = ISCreateStride(PETSC_COMM_SELF,m,start,1,&irow);CHKERRQ(ierr);
    ierr = ISCreateStride(PETSC_COMM_SELF,N,0,1,&icol);CHKERRQ(ierr);
    
    /* Use MatGetSubMatrices() to get a sequential matrix */
    ierr = MatGetSubMatrices(A,1,&irow,&icol,MAT_INITIAL_MATRIX,&smat);CHKERRQ(ierr);
    
    ierr = PCWSMP_ExtractUpperTriangular_MatSeqAIJ(A,smat[0],NULL,NULL,NULL,reuse,_vals_ut);CHKERRQ(ierr);
    
    ierr = MatDestroyMatrices(1,&smat);CHKERRQ(ierr);
    ierr = ISDestroy(&irow);CHKERRQ(ierr);
    ierr = ISDestroy(&icol);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

/* wrappers for WSMP to hide some ugly #if and if sequential type statements */
#undef __FUNCT__
#define __FUNCT__ "call_wsmp"
PetscErrorCode call_wsmp(PC_WSMP *wsmp)
{
    PetscFunctionBegin;
    if (wsmp->symmetric) {
        if (wsmp->sequential) {
            
#ifdef HAVE_WSSMP
            wssmp_ ( &wsmp->Nlocal, wsmp->IA, wsmp->JA, wsmp->AVALS, wsmp->DIAG, wsmp->PERM, wsmp->INVP, wsmp->B, &wsmp->LDB, &wsmp->NRHS,
                    &wsmp->AUX, &wsmp->NAUX, wsmp->MRP, wsmp->IPARM, wsmp->DPARM );
            
            if (wsmp->IPARM[64 -1] != 0) {
                SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"[wsmp] WSSMP generated the following error: %d",wsmp->IPARM[64 -1]);
            }
#else
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"[wsmp] Missing external package needed for type -pc_type \"wsmp\" (WSSMP)");
#endif
            
        } else {
#ifdef HAVE_PWSSMP
            pwssmp_ ( &wsmp->Nlocal, wsmp->IA, wsmp->JA, wsmp->AVALS, wsmp->DIAG, wsmp->PERM, wsmp->INVP, wsmp->B, &wsmp->LDB, &wsmp->NRHS,
                     &wsmp->AUX, &wsmp->NAUX, wsmp->MRP, wsmp->IPARM, wsmp->DPARM );
            
            if (wsmp->IPARM[64 -1] != 0) {
                SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"[wsmp] PWSSMP generated the following error: %d",wsmp->IPARM[64 -1]);
            }
#else
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"[wsmp] Missing external package needed for type -pc_type \"wsmp\" (PWSSMP)");
#endif
        }
    } else {
        
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "call_wsffree"
PetscErrorCode call_wsffree(PC_WSMP *wsmp)
{
    PetscFunctionBegin;
    if (wsmp->sequential) {
#ifdef HAVE_WSSMP
        wsffree_();
#endif
    } else {
#ifdef HAVE_PWSSMP
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
#ifdef HAVE_WSSMP
        wsmp_clear_();
#endif
    } else {
#ifdef HAVE_PWSSMP
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
#ifdef HAVE_PWSSMP
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
        if (!wsmp->symmetric) {
            SETERRQ(comm,PETSC_ERR_SUP,"[wsmp] Only support for symmetric WSMP implemenatations exists");
        }
        
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
        PetscMalloc(sizeof(int)   *wsmp->Nglobal,&wsmp->PERM); PetscMemzero(wsmp->PERM,sizeof(int)   *wsmp->Nglobal);
        PetscMalloc(sizeof(int)   *wsmp->Nglobal,&wsmp->INVP); PetscMemzero(wsmp->INVP,sizeof(int)   *wsmp->Nglobal);
        PetscMalloc(sizeof(double)*wsmp->Nlocal,&wsmp->B);     PetscMemzero(wsmp->B,   sizeof(double)*wsmp->Nlocal);
        PetscMalloc(sizeof(int)   *wsmp->Nlocal,&wsmp->MRP);   PetscMemzero(wsmp->MRP, sizeof(int)   *wsmp->Nlocal);
        
        /* Fetch ia, ja from matrix */
        PetscTime(&t0);
        if (wsmp->sequential) {
#ifndef USE_COMPACT_FORM
            ierr = PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
#else
            ierr = xxx_PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
#endif
        } else {
#ifndef USE_COMPACT_FORM
            ierr = PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
#else
            ierr = xxx_PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
#endif
        }
        PetscTime(&t1);
        PetscPrintf(comm," --- wsmp: done IJ --- %1.2e sec\n",t1-t0);
        
        /* -- [wsmp] : ordering -- */
        PetscTime(&t0);
        wsmp->IPARM[2 -1] = 1;
        wsmp->IPARM[3 -1] = 1;
        ierr = WSMPSetFromOptions_Ordering(pc,wsmp);CHKERRQ(ierr);
        ierr = call_wsmp(wsmp);CHKERRQ(ierr);
        PetscTime(&t1);
        PetscPrintf(comm," --- wsmp: done ordering --- %1.2e sec\n",t1-t0);
        
        /* -- [wsmp] : symbolic factorization -- */
        wsmp->IPARM[2 -1] = 2;
        wsmp->IPARM[3 -1] = 2;
        PetscTime(&t0);
        ierr = WSMPSetFromOptions_SymbolicFactorization(pc,wsmp);CHKERRQ(ierr);
        ierr = call_wsmp(wsmp);CHKERRQ(ierr);
        PetscTime(&t1);
        PetscPrintf(comm," --- wsmp: done sym fac --- %1.2e sec\n",t1-t0);
        
        PetscPrintf(comm,"[wsmp][sym. fact.] Num. nonzeros in factors = %d \n",wsmp->IPARM[24 -1]);
        //PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Estimated memory usage for factors = 1000 X %d \n",wsmp->IPARM[23 -1]);
        
        PetscPrintf(comm,"[wsmp][sym. fact.] Actual number of FLOPS in factorization =  %1.4e \n",wsmp->DPARM[23 -1]);
        PetscPrintf(comm,"[wsmp][sym. fact.] Factorization MegaFlops = %1.4e \n",wsmp->DPARM[23 -1]*1.0e-6 / (t1-t0) );
        
    } else {
        
        if (pc->flag != SAME_NONZERO_PATTERN) {
            /* Repeated solve but non-zero structure has changed, use new data */
            /* Re-constuct ordering and symbolic factorization */
            wsmp->nnz = 0;
            PetscFree(wsmp->IA);
            PetscFree(wsmp->JA);
            
            PetscTime(&t0);
            if (wsmp->sequential) {
#ifndef USE_COMPACT_FORM
                ierr = PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
#else
                ierr = xxx_PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
#endif
            } else {
#ifndef USE_COMPACT_FORM
                ierr = PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
#else
                ierr = xxx_PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
#endif
            }
            PetscTime(&t1);
            PetscPrintf(comm," --- wsmp: done IJ --- %1.2e sec\n",t1-t0);
            
            /* -- [wsmp] : ordering -- */
            PetscTime(&t0);
            wsmp->IPARM[2 -1] = 1;
            wsmp->IPARM[3 -1] = 1;
            ierr = call_wsmp(wsmp);CHKERRQ(ierr);
            PetscTime(&t1);
            PetscPrintf(comm," --- wsmp: done sym ordering --- %1.2e sec\n",t1-t0);
            
            /* -- [wsmp] : symbolic factorization -- */
            wsmp->IPARM[2 -1] = 2;
            wsmp->IPARM[3 -1] = 2;
            PetscTime(&t0);
            ierr = call_wsmp(wsmp);CHKERRQ(ierr);
            PetscTime(&t1);
            PetscPrintf(comm," --- wsmp: done sym fac --- %1.2e sec\n",t1-t0);
            
            PetscPrintf(comm,"[wsmp][sym. fact.] Num. nonzeros in factors = %d \n",wsmp->IPARM[24 -1]);
            //PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Estimated memory usage for factors = 1000 X %d \n",wsmp->IPARM[23 -1]);
            
            PetscPrintf(comm,"[wsmp][sym. fact.] Actual number of FLOPS in factorization =  %1.4e \n",wsmp->DPARM[23 -1]);
            PetscPrintf(comm,"[wsmp][sym. fact.] Factorization MegaFlops = %1.4e \n",wsmp->DPARM[23 -1]*1.0e-6 / (t1-t0) );
        }
    }
    
    /* Fetch matrix entries */
    PetscTime(&t0);
    if (!pc->setupcalled) {
        /* If first time we are in PCSetUp, use new data */
#ifndef USE_COMPACT_FORM
        ierr = PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_FALSE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
#else
        ierr = xxx_PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_FALSE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
#endif
    } else {
        if (pc->flag != SAME_NONZERO_PATTERN) {
            /* Repeated solve but non-zero structure has changed, use new data */
            PetscFree(wsmp->AVALS);
            /* release internal memory for factors */
            ierr = call_wsffree(wsmp);CHKERRQ(ierr);
            
#ifndef USE_COMPACT_FORM
            ierr = PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_FALSE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
#else
            ierr = xxx_PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_FALSE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
#endif
        } else {
            /* Repeated solve but non-zero structure is the same, re-use data */
#ifndef USE_COMPACT_FORM
            ierr = PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_TRUE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
#else
            ierr = xxx_PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_TRUE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
#endif
        }
    }
    PetscTime(&t1);
    PetscPrintf(comm," --- wsmp: done AIJ --- %1.2e sec\n",t1-t0);
    
    {
        PetscBool view = PETSC_FALSE;
        PetscOptionsGetBool(NULL,NULL,"-pc_wsmp_debug",&view,NULL);
        if (view) {
            ierr = PCWSMP_CSRView(B,wsmp->IA,wsmp->JA,wsmp->AVALS);CHKERRQ(ierr);
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
    ierr = call_wsmp(wsmp);CHKERRQ(ierr);
    PetscTime(&t1);
    PetscPrintf(comm," --- wsmp: done num fac --- %1.2e sec\n",t1-t0);
    //PetscPrintf(PETSC_COMM_SELF,"[wsmp][num. fact.] Actual memory usage for factors = 1000 X %d \n",wsmp->IPARM[23 -1]);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_WSMP"
static PetscErrorCode PCApply_WSMP(PC pc,Vec x,Vec y)
{
    PetscErrorCode ierr;
    PC_WSMP        *wsmp = (PC_WSMP*)pc->data;
    PetscScalar    *_val;
    PetscInt       i,m;
    static int     beenhere = 0;
    PetscBool      view = PETSC_FALSE;
    int            rank;
    MPI_Comm       comm;
    
    PetscFunctionBegin;
    PetscOptionsGetBool(NULL,NULL,"-pc_wsmp_debug",&view,NULL);
    
    ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
    MPI_Comm_rank(comm,&rank);
    //if (rank == 0) {
    //for (i=0; i<wsmp->Nglobal; i++) {
    //	printf("perm %d : iperm %d \n",wsmp->PERM[i],wsmp->INVP[i]);
    //}
    //}
    
    ierr = VecGetLocalSize(x,&m);CHKERRQ(ierr);
    ierr = VecGetArray(x,&_val);CHKERRQ(ierr);
    for (i=0; i<m; i++) {
        wsmp->B[i] = PetscRealPart(_val[i]);
    }
    ierr = VecRestoreArray(x,&_val);CHKERRQ(ierr);
    
    if (view && (beenhere==0)) {
        ierr = PCWSMP_VecView("rhs",wsmp);CHKERRQ(ierr);
    }
    
    /* apply L^T L factors  - Back substitution */
    wsmp->IPARM[2 -1] = 4;
    wsmp->IPARM[3 -1] = 4;
    ierr = WSMPSetFromOptions_BackSubstitution(pc,wsmp);CHKERRQ(ierr);
    ierr = call_wsmp(wsmp);CHKERRQ(ierr);
    
    /* iterative refinement */
    wsmp->IPARM[2 -1] = 5;
    wsmp->IPARM[3 -1] = 5;
    ierr = WSMPSetFromOptions_IterativeRefinement(pc,wsmp);CHKERRQ(ierr);
    ierr = call_wsmp(wsmp);CHKERRQ(ierr);
    
    if (view && (beenhere==0)) {
        ierr = PCWSMP_VecView("sol",wsmp);CHKERRQ(ierr);
    }
    
    ierr = VecGetArray(y,&_val);CHKERRQ(ierr);
    for (i=0; i<m; i++) {
        _val[i] = (PetscScalar)wsmp->B[i];
    }
    ierr = VecRestoreArray(y,&_val);CHKERRQ(ierr);
    
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
    PetscFree(wsmp->IA);
    PetscFree(wsmp->JA);
    PetscFree(wsmp->AVALS);
    if (wsmp->DIAG) { PetscFree(wsmp->DIAG); }
    PetscFree(wsmp->PERM);
    PetscFree(wsmp->INVP);
    PetscFree(wsmp->B);
    PetscFree(wsmp->MRP);
    
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
    
    ierr = PetscOptionsGetInt(NULL,prefix,"-pc_wsmp_iparm12",&ival,&found);CHKERRQ(ierr); if (found) { PetscMPIIntCast(ival,&wsmp->IPARM[index-1]); }
    
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
    
    /* -- [wsmp] : initialize with default parameters -- */
    wsmp->IPARM[1 -1] = 0;
    wsmp->IPARM[2 -1] = 0;
    wsmp->IPARM[3 -1] = 0;
    ierr = call_wsmp(wsmp);CHKERRQ(ierr);
    wsmp->IPARM[4 -1] = 0; /* CSR/CSC matrix format */
    wsmp->IPARM[5 -1] = 0; /* C style numbering */
    
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

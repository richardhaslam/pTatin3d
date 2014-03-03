
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <private/pcimpl.h>     /*I "petscpc.h" I*/
#include <petscksp.h>


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


PetscErrorCode WSMPSetFromOptions_Ordering(PC_WSMP *wsmp);
PetscErrorCode WSMPSetFromOptions_SymbolicFactorization(PC_WSMP *wsmp);
PetscErrorCode WSMPSetFromOptions_NumericFactorization(PC_WSMP *wsmp);


#undef __FUNCT__
#define __FUNCT__ "Default_MatIsSymmetric"
PetscErrorCode Default_MatIsSymmetric(Mat A,PetscReal tol,PetscBool *flg)
{
	PetscErrorCode ierr;
	Vec x,y1,y2;
	PetscReal min,max;
	
	ierr = MatGetVecs(A,&y2,&y1);CHKERRQ(ierr);
	ierr = MatGetVecs(A,&x,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = VecSetRandom(x,PETSC_NULL);CHKERRQ(ierr);
	ierr = MatMult(A,x,y1);CHKERRQ(ierr);
	ierr = MatMultTranspose(A,x,y2);CHKERRQ(ierr);
	
	ierr = VecAXPY(y1,-1.0,y2);CHKERRQ(ierr);
	ierr = VecAbs(y1);CHKERRQ(ierr);
	ierr = VecMin(y1,PETSC_NULL,&min);CHKERRQ(ierr);
	ierr = VecMax(y1,PETSC_NULL,&max);CHKERRQ(ierr);
	
	if (max < tol) {
		*flg = PETSC_TRUE;
	} else {
		*flg = PETSC_FALSE;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCWSMP_MatIsSymmetric"
PetscErrorCode PCWSMP_MatIsSymmetric(Mat A,PetscReal tol,PetscBool *flg)
{
	PetscErrorCode ierr;
	PetscMPIInt size;
	MPI_Comm comm;
	
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
#define __FUNCT__ "PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ"
PetscErrorCode PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(Mat A,int *_nnz_ut,int **_ia_ut,int **_ja_ut)
{
	PetscErrorCode ierr;
	PetscInt *ia;
	PetscInt *ja;
	PetscInt m,nnz_i,cnt,idx,nnz_full,nnz_ut,i,j;
	PetscBool done;
	int *ia_ut,*ja_ut;
	
	
	ierr = MatGetRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
	if (!done) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MatGetRowIJ failed... aborting");
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
	//PetscPrintf(PETSC_COMM_WORLD,"nnz_ut = %d \n",nnz_ut);
	
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
		
		idx = idx + nnz_i_ut;
		cnt = cnt + nnz_i;
	}
	ia_ut[m] = (int)idx + 0; /* fortran +1 */
	
	/* ja_ut */
	idx = 0;
	cnt = 0;
	for (i=0; i<m; i++) {
		
		nnz_i = ia[i+1]-ia[i];
		for (j=cnt; j<cnt+nnz_i; j++) {
			if (ja[j] >= i) {
				ja_ut[idx] = (int)ja[j] + 0; /* fortran +1 */
				idx++;
			}
		}
		
		cnt = cnt + nnz_i;
	}
	
	ierr = MatRestoreRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
	if (!done) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MatRestoreRowIJ failed... aborting");
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
	PetscInt *ia;
	PetscInt *ja;
	PetscInt m,n,M,N,nnz_i,cnt,idx,nnz_ut,i,j,start,end;
	PetscBool done;
	int *ia_ut,*ja_ut;
	IS irow,icol;
	Mat *smat;
	PetscBool view = PETSC_FALSE;
	
	
	ierr = MatGetSize(A,&M,&N);CHKERRQ(ierr);
	ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
	
	ierr = ISCreateStride(PETSC_COMM_SELF,m,start,1,&irow);CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_SELF,N,0,1,&icol);CHKERRQ(ierr);
	
	/* Use MatGetSubMatrices() to get a sequential matrix */
	ierr = MatGetSubMatrices(A,1,&irow,&icol,MAT_INITIAL_MATRIX,&smat);CHKERRQ(ierr);
	
	ierr = MatGetRowIJ(smat[0],0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
	if (!done) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MatGetRowIJ failed... aborting");
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
	//PetscPrintf(PETSC_COMM_SELF,"nnz_ut = %d \n",nnz_ut);
	
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
		
		idx = idx + nnz_i_ut;
		cnt = cnt + nnz_i;
	}
	ia_ut[m] = (int)idx;
	
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
	
	if (view) {
		PetscMPIInt rank;
		FILE *fp;
		char name[1000];
		PetscInt nnz_i_ut = 0;
		
		// ia
		MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		sprintf(name,"ia_rank%d.dat",rank);
		fp = fopen(name,"w");
        
		cnt = 0;
		idx = 0;
		for (i=0; i<m; i++) {
			nnz_i_ut = 0;
            
			nnz_i = ia[i+1]-ia[i];
			for (j=cnt; j<cnt+nnz_i; j++) {
				if (ja[j] >= i+start) {
					nnz_i_ut++;
				}
			}
			fprintf(fp,"[%d] %d \n",i,(int)idx);
			
			idx = idx + nnz_i_ut;
			cnt = cnt + nnz_i;
		}
		fprintf(fp,"[%d] %d \n",m,(int)idx);
		
		fclose(fp);
        
		// ja
		sprintf(name,"ja_rank%d.dat",rank);
		fp = fopen(name,"w");
        
		idx = 0;
		cnt = 0;
		for (i=0; i<m; i++) {
			
			nnz_i = ia[i+1]-ia[i];
			for (j=cnt; j<cnt+nnz_i; j++) {
				if (ja[j] >= i+start) {
					fprintf(fp,"[%d] %d %d \n",idx,i,(int)ja[j]);
					idx++;
				}
			}
			cnt = cnt + nnz_i;
		}
		
		fclose(fp);
		
	}
	
	ierr = MatRestoreRowIJ(smat[0],0,PETSC_FALSE,PETSC_FALSE,&m,&ia,&ja,&done);CHKERRQ(ierr);
	if (!done) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MatRestoreRowIJ failed... aborting");
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
	PetscErrorCode ierr;
	PetscInt m,n,i,j,start,end;
	PetscInt ncols,idx;
	const PetscInt *cols;
	const PetscScalar *vals;
	double *vals_ut;
	PetscBool view = PETSC_FALSE;
	
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
    
	if (view) {
		PetscMPIInt rank;
		FILE *fp;
		char name[1000];
		
		//ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_COMMON);CHKERRQ(ierr);
		//ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		
		MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		sprintf(name,"Aut_rank%d.dat",rank);
		fp = fopen(name,"w");
		for (i=start; i<start+m; i++) {
			ierr = MatGetRow(A,i,&ncols,&cols,&vals);CHKERRQ(ierr);
			for (j=0; j<ncols; j++) {
				if (cols[j] >= i) {
					fprintf(fp,"[%d,%d] %1.4e: ",i,cols[j],vals[j]);
				}
			} fprintf(fp,"\n");
			ierr = MatRestoreRow(A,i,&ncols,&cols,&vals);CHKERRQ(ierr);
		}
		fclose(fp);
	}
	
	
	if (!reuse) {
		*_vals_ut = vals_ut;
	}
	
	PetscFunctionReturn(0);
}

/* wrappers for WSMP to hide some ugly #if and if sequential type statements */
#undef __FUNCT__
#define __FUNCT__ "call_wsmp"
PetscErrorCode call_wsmp(PC_WSMP *wsmp)
{
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
	MPI_Fint fcomm;
	
	fcomm = MPI_Comm_c2f(comm);
#ifdef HAVE_PWSSMP
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
    PetscErrorCode   ierr;
    PC_WSMP          *wsmp = (PC_WSMP*)pc->data;
	Mat A,B;
	PetscLogDouble   t0,t1;
	
    PetscFunctionBegin;
	
	ierr = PCGetOperators(pc,&A,&B,PETSC_NULL);CHKERRQ(ierr);
	
	if (!pc->setfromoptionscalled) {
		/* -- [wsmp] : initialize with default parameters -- */
		wsmp->IPARM[1 -1] = 0;
		wsmp->IPARM[2 -1] = 0;
		wsmp->IPARM[3 -1] = 0;
		ierr = call_wsmp(wsmp);CHKERRQ(ierr);
		wsmp->IPARM[4 -1] = 0; /* CSR/CSC matrix format */
		wsmp->IPARM[5 -1] = 0; /* C style numbering */
	}
	
	/* construction phase */
    if (!pc->setupcalled) {
		PetscInt M,N,m,n;
		
		/* Determine if matrix is symmetric */
		ierr = PCWSMP_MatIsSymmetric(B,1.0e-8,&wsmp->symmetric);CHKERRQ(ierr);
		if (!wsmp->symmetric) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"[wsmp] Only support for symmetric WSMP implemenatations exists");
		}
		
		ierr = MatGetSize(B,&M,&N);CHKERRQ(ierr);
		ierr = MatGetLocalSize(B,&m,&n);CHKERRQ(ierr);
		
		/* get local size */
		if (m != n) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"[wsmp] WSMP only supports square matrices");
		}
		
		wsmp->Nglobal = (int)M;
		wsmp->Nlocal = (int)m;
		wsmp->LDB = wsmp->Nlocal;
		wsmp->NRHS = 1;
		wsmp->NAUX = 0;
		
		/* allocate other local variables */
		wsmp->DIAG = NULL;
		PetscMalloc(sizeof(int)*wsmp->Nglobal,&wsmp->PERM);
		PetscMalloc(sizeof(int)*wsmp->Nglobal,&wsmp->INVP);
		PetscMalloc(sizeof(double)*wsmp->Nlocal,&wsmp->B);
		PetscMalloc(sizeof(int)*wsmp->Nlocal,&wsmp->MRP);
		
		/* Fetch ia, ja from matrix */
		if (wsmp->sequential) {
			ierr = PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
		} else {
			ierr = PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
		}
		
		
		/* -- [wsmp] : ordering -- */
		wsmp->IPARM[2 -1] = 1;
		wsmp->IPARM[3 -1] = 1;
		ierr = WSMPSetFromOptions_Ordering(wsmp);CHKERRQ(ierr);
		ierr = call_wsmp(wsmp);CHKERRQ(ierr);
		
		/* -- [wsmp] : symbolic factorization -- */
		wsmp->IPARM[2 -1] = 2;
		wsmp->IPARM[3 -1] = 2;
		PetscGetTime(&t0);
		ierr = WSMPSetFromOptions_SymbolicFactorization(wsmp);CHKERRQ(ierr);
		ierr = call_wsmp(wsmp);CHKERRQ(ierr);
		PetscGetTime(&t1);
		
		PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Num. nonzeros in factors = %d \n",wsmp->IPARM[24 -1]);
		PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Estimated memory usage for factors = 1000 X %d \n",wsmp->IPARM[23 -1]);
        
		PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Actual number of FLOPS in factorization =  %1.4e \n",wsmp->DPARM[23 -1]);
		PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Factorization MegaFlops = %1.4e \n",wsmp->DPARM[23 -1]*1.0e-6 / (t1-t0) );
        
    } else {
		
		if (pc->flag != SAME_NONZERO_PATTERN) {
			/* Repeated solve but non-zero structure has changed, use new data */
			/* Re-constuct ordering and symbolic factorization */
			wsmp->nnz = 0;
			PetscFree(wsmp->IA);
			PetscFree(wsmp->JA);
			
			if (wsmp->sequential) {
				ierr = PCWSMP_ExtractUpperTriangularIJ_MatSeqAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
			} else {
				ierr = PCWSMP_ExtractUpperTriangularIJ_MatMPIAIJ(B,&wsmp->nnz,&wsmp->IA,&wsmp->JA);CHKERRQ(ierr);
			}
			
			/* -- [wsmp] : ordering -- */
			wsmp->IPARM[2 -1] = 1;
			wsmp->IPARM[3 -1] = 1;
			ierr = WSMPSetFromOptions_Ordering(wsmp);CHKERRQ(ierr);
			ierr = call_wsmp(wsmp);CHKERRQ(ierr);
			
			/* -- [wsmp] : symbolic factorization -- */
			wsmp->IPARM[2 -1] = 2;
			wsmp->IPARM[3 -1] = 2;
            PetscGetTime(&t0);
			ierr = WSMPSetFromOptions_SymbolicFactorization(wsmp);CHKERRQ(ierr);
			ierr = call_wsmp(wsmp);CHKERRQ(ierr);
            PetscGetTime(&t1);
            
			PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Num. nonzeros in factors = %d \n",wsmp->IPARM[24 -1]);
			PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Estimated memory usage for factors = 1000 X %d \n",wsmp->IPARM[23 -1]);
			
			PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Actual number of FLOPS in factorization =  %1.4e \n",wsmp->DPARM[23 -1]);
			PetscPrintf(PETSC_COMM_SELF,"[wsmp][sym. fact.] Factorization MegaFlops = %1.4e \n",wsmp->DPARM[23 -1]*1.0e-6 / (t1-t0) );
		}
		
	}
	
	/* Fetch matrix entries */
	
    if (!pc->setupcalled) {
		/* If first time we are in PCSetUp, use new data */
		ierr = PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_FALSE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
	} else {
		if (pc->flag != SAME_NONZERO_PATTERN) {
			/* Repeated solve but non-zero structure has changed, use new data */
			PetscFree(wsmp->AVALS);
			/* release internal memory for factors */
			ierr = call_wsffree(wsmp);CHKERRQ(ierr);
			
			ierr = PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_FALSE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
		} else {
			/* Repeated solve but non-zero structure is the same, re-use data */
			ierr = PCWSMP_ExtractUpperTriangularAIJ(B,PETSC_TRUE,wsmp->nnz,&wsmp->AVALS);CHKERRQ(ierr);
		}
	}
	
	/* -- [wsmp] : numeric factorization - Cholesky -- */
	wsmp->IPARM[2 -1] = 3;
	wsmp->IPARM[3 -1] = 3;
	ierr = WSMPSetFromOptions_NumericFactorization(wsmp);CHKERRQ(ierr);
	ierr = call_wsmp(wsmp);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_SELF,"[wsmp][num. fact.] Actual memory usage for factors = 1000 X %d \n",wsmp->IPARM[23 -1]);
	
	//if (params%iproc==0) write(816,'(a,es15.4)') 'Factorization MegaFlops          = ',(system%dparm(23) * 1.d-6) / (tt2-tt1)
	
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PCApply_WSMP"
static PetscErrorCode PCApply_WSMP(PC pc,Vec x,Vec y)
{
    PetscErrorCode   ierr;
    PC_WSMP          *wsmp = (PC_WSMP*)pc->data;
	PetscScalar      *_val;
	PetscInt         i,m;
	
	
	ierr = VecGetLocalSize(x,&m);CHKERRQ(ierr);
	ierr = VecGetArray(x,&_val);CHKERRQ(ierr);
	for (i=0; i<m; i++) {
		wsmp->B[i] = PetscRealPart(_val[i]);
	}
	ierr = VecRestoreArray(x,&_val);CHKERRQ(ierr);
	
	/* apply L^T L factors  - Back substitution */
	/* 1/ set args for wsmp call */
	wsmp->IPARM[2 -1] = 4;
	wsmp->IPARM[3 -1] = 4;
	ierr = call_wsmp(wsmp);CHKERRQ(ierr);
	
	ierr = VecGetArray(y,&_val);CHKERRQ(ierr);
	for (i=0; i<m; i++) {
		_val[i] = (PetscScalar)wsmp->B[i];
	}
	ierr = VecRestoreArray(y,&_val);CHKERRQ(ierr);
	
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
    PetscErrorCode   ierr;
    PC_WSMP          *wsmp = (PC_WSMP*)pc->data;
	
	
    PetscFunctionBegin;
    ierr = PCReset_WSMP(pc);CHKERRQ(ierr);
	
    ierr = PetscFree(pc->data);CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WSMPSetFromOptions_Ordering"
PetscErrorCode WSMPSetFromOptions_Ordering(PC_WSMP *wsmp)
{
	PetscErrorCode ierr;
    
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
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WSMPSetFromOptions_SymbolicFactorization"
PetscErrorCode WSMPSetFromOptions_SymbolicFactorization(PC_WSMP *wsmp)
{
	PetscErrorCode ierr;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WSMPSetFromOptions_NumericFactorization"
PetscErrorCode WSMPSetFromOptions_NumericFactorization(PC_WSMP *wsmp)
{
	PetscErrorCode ierr;
    
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetFromOptions_WSMP"
static PetscErrorCode PCSetFromOptions_WSMP(PC pc)
{
    PetscErrorCode   ierr;
    PC_WSMP          *wsmp = (PC_WSMP*)pc->data;
	
	
    PetscFunctionBegin;
    ierr = PetscOptionsHead("WSMP options");CHKERRQ(ierr);
	
	/* -- [wsmp] : initialize with default parameters -- */
	wsmp->IPARM[1 -1] = 0;
	wsmp->IPARM[2 -1] = 0;
	wsmp->IPARM[3 -1] = 0;
	ierr = call_wsmp(wsmp);CHKERRQ(ierr);
	wsmp->IPARM[4 -1] = 0; /* CSR/CSC matrix format */
	wsmp->IPARM[5 -1] = 0; /* C style numbering */
	
	/* put option calls here to overide defaults */
	// ierr = PetscOptionsInt("-pc_wsmp","Factor to reduce parent communication size by","PCSemiRedundantSetFactor",red->nsubcomm_factor,&red->nsubcomm_factor,0);CHKERRQ(ierr);
	
	// iterative refinement steps
	//wsmp->IPARM[6 -1] = X;
	
	// pivot styles IPARM[11 -1] = 0, 1, 2
	// pivot values DPARM[10 -1], DPARM[64 -1], DPARM[21 -1]
    
	ierr = WSMPSetFromOptions_Ordering(wsmp);CHKERRQ(ierr);
	ierr = WSMPSetFromOptions_SymbolicFactorization(wsmp);CHKERRQ(ierr);
	ierr = WSMPSetFromOptions_NumericFactorization(wsmp);CHKERRQ(ierr);
	
	
	ierr = PetscOptionsTail();CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCView_WSMP"
static PetscErrorCode PCView_WSMP(PC pc,PetscViewer viewer)
{
    PetscErrorCode   ierr;
    PC_WSMP          *wsmp = (PC_WSMP*)pc->data;
    PetscBool        iascii,isstring;
	
	
    PetscFunctionBegin;
    ierr = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
    ierr = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERSTRING,&isstring);CHKERRQ(ierr);
    if (iascii) {
		
		ierr = PetscViewerASCIIPrintf(viewer,"  WSMP preconditioner\n");CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(4)  = %d <matrix format>\n",wsmp->IPARM[4 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(5)  = %d <C or ForTran numbering>\n",wsmp->IPARM[5 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(6)  = %d <number of iterative refinement steps performed>\n",wsmp->IPARM[6 -1]);
		if ( wsmp->IPARM[7 -1] <= 3) {
			ierr = PetscViewerASCIIPrintf(viewer,"    iparm(7)  = 0,1,2,3 <double precision>\n");
		} else {
			ierr = PetscViewerASCIIPrintf(viewer,"    iparm(7)  = 4,5,6,7 <quadruple precision>\n");
		}
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(8)  = %d <matrix permutation option>\n",wsmp->IPARM[8 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(9)  = %d <rhs permutation option>\n",wsmp->IPARM[9 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(10) = %d <matrix scaling option>\n",wsmp->IPARM[10 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(11) = %d <pivot type>\n",wsmp->IPARM[11 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(16) = %d <ordering type for perm. vectors>\n",wsmp->IPARM[16 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(17) = %d <max. num. nodes in subgraph>\n",wsmp->IPARM[17 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(18) = %d <force minimum local fill>\n",wsmp->IPARM[18 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(19) = %d <random seed for graph permutation>\n",wsmp->IPARM[19 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(20) = %d <matrix characteristics flag\n",wsmp->IPARM[20 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(22) = %d <num. negative eigenvalues>\n",wsmp->IPARM[22 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(23) = %d <double words for factorisation>\n",wsmp->IPARM[23 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(24) = %d <nnz in triangular factor>\n",wsmp->IPARM[24 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    iparm(33) = %d <num cpus>\n",wsmp->IPARM[33 -1]);
		
		if (!wsmp->sequential) {
			
		}
		
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(4)  = %1.4e <largest diagonal entry>\n",wsmp->DPARM[4 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(5)  = %1.4e <smallest diagonal entry>\n",wsmp->DPARM[5 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(6)  = %1.4e <relative tol for iterative refinement>\n",wsmp->DPARM[6 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(7)  = %1.4e <relative error for iterative refinement>\n",wsmp->DPARM[7 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(10) = %1.4e <lower threshold on diag without pivoting>\n",wsmp->DPARM[10 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(11) = %1.4e <Bunch-Kaufman pivoting threshold>\n",wsmp->DPARM[11 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(12) = %1.4e <pivoting control>\n",wsmp->DPARM[12 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(13) = %1.4e <num. supernodes>\n",wsmp->DPARM[13 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(23) = %1.4e <num. floating point operations>\n",wsmp->DPARM[23 -1]);
		ierr = PetscViewerASCIIPrintf(viewer,"    dparm(24) = %1.4e <expected num. floating point operations>\n",wsmp->DPARM[24 -1]);
		
    } else if (isstring) {
		ierr = PetscViewerStringSPrintf(viewer," WSMP preconditioner");CHKERRQ(ierr);
    } else {
        SETERRQ1(((PetscObject)pc)->comm,PETSC_ERR_SUP,"Viewer type %s not supported for PC WSMP",((PetscObject)viewer)->type_name);
    }
	
    PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCCreate_WSMP"
PetscErrorCode PCCreate_WSMP(PC pc)
{
    PetscErrorCode   ierr;
    PC_WSMP          *wsmp;
    PetscMPIInt      size;
    
	
    PetscFunctionBegin;
    ierr = PetscNewLog(pc,PC_WSMP,&wsmp);CHKERRQ(ierr);
    pc->data            = (void*)wsmp;
	
	/* determine if sequential call or parallel call required */
	ierr = MPI_Comm_size(((PetscObject)pc)->comm,&size);CHKERRQ(ierr);
	if (size == 1) {
		wsmp->sequential = PETSC_TRUE;
	} else {
		wsmp->sequential = PETSC_FALSE;
	}
	
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

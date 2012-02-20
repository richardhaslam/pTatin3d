/*
 
 Notes about fieldsplit - matnest
 1) MatNest requires MatMultAdd() be defined on the sub-matrice
 
 How is the best way to configure the sub matrices

 When using MatGetSubMatrix from the MF operator for the complete stokes operator, do this
 -stokes_Amf_A11_mf (XXX A11_mf)
 
 -stokes_Amf_A11_Auu_mf (XXX A11_Auu_mf)
 -stokes_Amf_A11_Avv_mf
 -stokes_Amf_A11_Aww_mf

 When using the matrix free A11 operator, do this
 -stokes_A_A11_Auu_mf  (XXX Auu_mf)
 -stokes_A_A11_Avv_mf
 -stokes_A_A11_Aww_mf

*/


#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
#include "private/matimpl.h" /*I   "petscmat.h"   I*/


#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "swarm_fields.h"

#include "dmda_duplicate.h"
#include "dmda_bcs.h"
#include "quadrature.h"
#include "stokes_operators.h"

struct _p_MatStokesMF {
	PetscInt    mu,mp,Mu,Mp;
	PetscInt    level;
	PetscInt    ii;
	DM          stokes_pack,daUVW,dap;
	BCList      u_bclist,p_bclist;
	Quadrature  volQ;
	DM          daU;
	IS          isUVW,isU,isV,isW,isP; /* Need the IS's for GetSubMatrix */
	PetscInt    refcnt;
};

struct _p_MatA11MF {
	PetscInt    mu,Mu;
	DM          daUVW;
	BCList      u_bclist;
	Quadrature  volQ;
	DM          daU; /* Optionally need this */
	IS          isUVW; /* Needed for full column space */
	IS          isU,isV,isW; /* Optionally: Need the IS's for GetSubMatrix */
	PetscInt    refcnt;
	/* Not sure I need this at all */
	PetscInt    level;
	PetscInt    ii;
};


#undef __FUNCT__  
#define __FUNCT__ "MatStokesMFCreate"
PetscErrorCode MatStokesMFCreate(MatStokesMF *B)
{
	PetscErrorCode ierr;
	MatStokesMF Stk;
	PetscFunctionBegin;
	
	ierr = PetscMalloc(sizeof(struct _p_MatStokesMF),&Stk);CHKERRQ(ierr);
	ierr = PetscMemzero(Stk,sizeof(struct _p_MatStokesMF));CHKERRQ(ierr);
	
	*B = Stk;
	PetscFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "MatA11MFCreate"
PetscErrorCode MatA11MFCreate(MatA11MF *B)
{
	PetscErrorCode ierr;
	MatA11MF A11;
	PetscFunctionBegin;
	
	ierr = PetscMalloc(sizeof(struct _p_MatA11MF),&A11);CHKERRQ(ierr);
	ierr = PetscMemzero(A11,sizeof(struct _p_MatA11MF));CHKERRQ(ierr);
	
	*B = A11;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatStokesMFSetup"
PetscErrorCode MatStokesMFSetup(MatStokesMF StkCtx,PhysCompStokes user)
{
	PetscErrorCode ierr;
	Vec X,u,p;
	PetscInt mu,mp,Mu,Mp;
	DM dau,dap,pack;
	IS             *is;
	PetscInt n,start,offset;
	PetscBool same;
	
	PetscFunctionBegin;
	
	StkCtx->stokes_pack = user->stokes_pack;   ierr = PetscObjectReference((PetscObject)user->stokes_pack);CHKERRQ(ierr);
	StkCtx->daUVW       = user->dav;           ierr = PetscObjectReference((PetscObject)user->dav);CHKERRQ(ierr); 
	StkCtx->dap         = user->dap;           ierr = PetscObjectReference((PetscObject)user->dap);CHKERRQ(ierr);
	StkCtx->volQ        = user->volQ;
	StkCtx->u_bclist    = user->u_bclist;
	StkCtx->p_bclist    = user->p_bclist;
	
	pack = user->stokes_pack;
	
	/* is composite */
	same = PETSC_FALSE;
	ierr = PetscTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
  /* Fetch the DA's */
	dau = user->dav;
	dap = user->dap;
	
	/* Sizes */
	ierr = DMGetGlobalVector(pack,&X);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(pack,X,&u,&p);CHKERRQ(ierr);
	ierr = VecGetSize(u,&Mu);CHKERRQ(ierr);
	ierr = VecGetLocalSize(u,&mu);CHKERRQ(ierr);
	ierr = VecGetSize(p,&Mp);CHKERRQ(ierr);
	ierr = VecGetLocalSize(p,&mp);CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(u,&start,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(pack,X,&u,&p);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(pack,&X);CHKERRQ(ierr);
	
	StkCtx->mu = mu;
	StkCtx->mp = mp;
	StkCtx->Mu = Mu;
	StkCtx->Mp = Mp;
	
	ierr = DMCompositeGetGlobalISs(pack,&is);CHKERRQ(ierr);
	StkCtx->isUVW = is[0];
	StkCtx->isP   = is[1];
	ierr = PetscFree(is);CHKERRQ(ierr);
	
	n = (mu/3);
	offset = start + 0;
	ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&StkCtx->isU);CHKERRQ(ierr);
	offset = start + 1;
	ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&StkCtx->isV);CHKERRQ(ierr);
	offset = start + 2;
	ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&StkCtx->isW);CHKERRQ(ierr);
	
	ierr = DMDADuplicateLayout(dau,1,2,DMDA_STENCIL_BOX,&StkCtx->daU);CHKERRQ(ierr);
	
	
	StkCtx->refcnt = 1;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatA11MFSetup"
PetscErrorCode MatA11MFSetup(MatA11MF A11Ctx,DM dav,Quadrature volQ,BCList u_bclist)
{
	PetscErrorCode ierr;
	Vec X;
	PetscInt mu,Mu;
	DM dau,dap,pack;
	IS  *is;
	PetscInt n,start,offset;
	PetscBool same;
	
	PetscFunctionBegin;
	
	A11Ctx->daUVW       = dav;           ierr = PetscObjectReference((PetscObject)dav);CHKERRQ(ierr); 
	A11Ctx->volQ        = volQ;
	A11Ctx->u_bclist    = u_bclist;
	
  /* Fetch the DA's */
	dau = dav;
	
	/* Sizes */
	ierr = DMGetGlobalVector(dau,&X);CHKERRQ(ierr);
	ierr = VecGetSize(X,&Mu);CHKERRQ(ierr);
	ierr = VecGetLocalSize(X,&mu);CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(X,&start,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(dau,&X);CHKERRQ(ierr);
	
	A11Ctx->mu = mu;
	A11Ctx->Mu = Mu;
	
	n = (mu/3);
	offset = start + 0;
	ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&A11Ctx->isU);CHKERRQ(ierr);
	offset = start + 1;
	ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&A11Ctx->isV);CHKERRQ(ierr);
	offset = start + 2;
	ierr = ISCreateStride(PETSC_COMM_WORLD, n,offset,3,&A11Ctx->isW);CHKERRQ(ierr);
	
	ierr = DMDADuplicateLayout(dau,1,2,DMDA_STENCIL_BOX,&A11Ctx->daU);CHKERRQ(ierr);
	
	A11Ctx->refcnt = 1;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatStokesMFDestroy"
PetscErrorCode MatStokesMFDestroy(MatStokesMF *B)
{
	MatStokesMF A;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if(!B) { PetscFunctionReturn(0); }
	A = *B;
	
	A->refcnt--;
	if (A->refcnt==0) {
		ierr = DMDestroy(&A->daUVW);CHKERRQ(ierr);
		ierr = DMDestroy(&A->dap);CHKERRQ(ierr);
		ierr = DMDestroy(&A->stokes_pack);CHKERRQ(ierr);
		
		ierr = ISDestroy(&A->isU);CHKERRQ(ierr);
		ierr = ISDestroy(&A->isV);CHKERRQ(ierr);
		ierr = ISDestroy(&A->isW);CHKERRQ(ierr);
		ierr = ISDestroy(&A->isP);CHKERRQ(ierr);
		ierr = ISDestroy(&A->isUVW);CHKERRQ(ierr);
		ierr = DMDestroy(&A->daU);CHKERRQ(ierr);
		ierr = PetscFree(A);CHKERRQ(ierr);

		*B = PETSC_NULL;
	}
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatA11MFDestroy"
PetscErrorCode MatA11MFDestroy(MatA11MF *B)
{
	MatA11MF A;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	
	if(!B) { PetscFunctionReturn(0); }
	A = *B;
	
	A->refcnt--;
	if (A->refcnt==0) {
		ierr = DMDestroy(&A->daUVW);CHKERRQ(ierr);
		ierr = ISDestroy(&A->isUVW);CHKERRQ(ierr);
		
		ierr = ISDestroy(&A->isU);CHKERRQ(ierr);
		ierr = ISDestroy(&A->isV);CHKERRQ(ierr);
		ierr = ISDestroy(&A->isW);CHKERRQ(ierr);
		ierr = DMDestroy(&A->daU);CHKERRQ(ierr);
		ierr = PetscFree(A);CHKERRQ(ierr);

		*B = PETSC_NULL;
	}
	
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatDestroy_MatStokesMF"
PetscErrorCode MatDestroy_MatStokesMF(Mat A)
{
	MatStokesMF       ctx;
	PetscErrorCode  ierr;
	
  PetscFunctionBegin;
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&ctx);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "MatDestroy_MatA11MF"
PetscErrorCode MatDestroy_MatA11MF(Mat A)
{
	MatA11MF       ctx;
	PetscErrorCode  ierr;
	
  PetscFunctionBegin;
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	ierr = MatA11MFDestroy(&ctx);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatStokesMFCopy"
PetscErrorCode MatStokesMFCopy(MatStokesMF A,MatStokesMF *B)
{
	PetscErrorCode ierr;
	MatStokesMF Stk;
	
	PetscFunctionBegin;

	ierr = MatStokesMFCreate(&Stk);CHKERRQ(ierr);
	
	Stk->mu    = A->mu;
	Stk->mp    = A->mp;
	Stk->Mu    = A->Mu;
	Stk->Mp    = A->Mp;
	Stk->level = A->level;
	Stk->ii    = A->ii;

	Stk->u_bclist = A->u_bclist;
	Stk->p_bclist = A->p_bclist;
	Stk->volQ     = A->volQ;
	
	Stk->stokes_pack  = A->stokes_pack;    ierr = PetscObjectReference((PetscObject)A->stokes_pack);CHKERRQ(ierr);
	
	Stk->isU   = A->isU;             ierr = PetscObjectReference((PetscObject)A->isU);CHKERRQ(ierr);
	Stk->isV   = A->isV;             ierr = PetscObjectReference((PetscObject)A->isV);CHKERRQ(ierr);
	Stk->isW   = A->isW;             ierr = PetscObjectReference((PetscObject)A->isW);CHKERRQ(ierr);
	Stk->isP   = A->isP;             ierr = PetscObjectReference((PetscObject)A->isP);CHKERRQ(ierr);
	Stk->isUVW = A->isUVW;           ierr = PetscObjectReference((PetscObject)A->isUVW);CHKERRQ(ierr);
	
	Stk->daUVW = A->daUVW;        ierr = PetscObjectReference((PetscObject)A->daUVW);CHKERRQ(ierr);
	Stk->daU   = A->daU;          ierr = PetscObjectReference((PetscObject)A->daU);CHKERRQ(ierr);
	Stk->dap   = A->dap;          ierr = PetscObjectReference((PetscObject)A->dap);CHKERRQ(ierr);
	
	Stk->refcnt = 1;
	
	*B = Stk;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatA11MFCopy"
PetscErrorCode MatA11MFCopy(MatA11MF A,MatA11MF *B)
{
	PetscErrorCode ierr;
	MatA11MF A11;
	
	PetscFunctionBegin;
	
	ierr = MatA11MFCreate(&A11);CHKERRQ(ierr);
	
	A11->mu    = A->mu;
	A11->Mu    = A->Mu;
	A11->level = A->level;
	A11->ii    = A->ii;
	
	A11->u_bclist = A->u_bclist;
	A11->volQ     = A->volQ;
	
	A11->isU   = A->isU;             ierr = PetscObjectReference((PetscObject)A->isU);CHKERRQ(ierr);
	A11->isV   = A->isV;             ierr = PetscObjectReference((PetscObject)A->isV);CHKERRQ(ierr);
	A11->isW   = A->isW;             ierr = PetscObjectReference((PetscObject)A->isW);CHKERRQ(ierr);
	
	A11->daUVW = A->daUVW;        ierr = PetscObjectReference((PetscObject)A->daUVW);CHKERRQ(ierr);
	A11->daU   = A->daU;          ierr = PetscObjectReference((PetscObject)A->daU);CHKERRQ(ierr);
	
	A11->refcnt = 1;
	
	*B = A11;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatCopy_StokesMF_A11MF"
PetscErrorCode MatCopy_StokesMF_A11MF(MatStokesMF A,MatA11MF *B)
{
	PetscErrorCode ierr;
	MatA11MF A11;
	
	PetscFunctionBegin;
	
	ierr = MatA11MFCreate(&A11);CHKERRQ(ierr);
	
	A11->mu    = A->mu;
	A11->Mu    = A->Mu;
	A11->level = A->level;
	A11->ii    = A->ii;
	
	A11->u_bclist = A->u_bclist;
	A11->volQ     = A->volQ;
	
	A11->isU   = A->isU;             ierr = PetscObjectReference((PetscObject)A->isU);CHKERRQ(ierr);
	A11->isV   = A->isV;             ierr = PetscObjectReference((PetscObject)A->isV);CHKERRQ(ierr);
	A11->isW   = A->isW;             ierr = PetscObjectReference((PetscObject)A->isW);CHKERRQ(ierr);
	
	A11->daUVW = A->daUVW;        ierr = PetscObjectReference((PetscObject)A->daUVW);CHKERRQ(ierr);
	A11->daU   = A->daU;          ierr = PetscObjectReference((PetscObject)A->daU);CHKERRQ(ierr);
	
	A11->refcnt = 1;
	
	*B = A11;
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "MatMultAdd_basic"
PetscErrorCode MatMultAdd_basic(Mat A,Vec v1,Vec v2,Vec v3)
{
	PetscErrorCode ierr;
	PetscScalar *LA_v2,*LA_v3;
	PetscInt i,n;
  PetscFunctionBegin;
	
	ierr = MatMult(A,v1,v3);CHKERRQ(ierr);
	ierr = VecGetLocalSize(v2,&n);CHKERRQ(ierr);
	ierr = VecGetArray(v2,&LA_v2);CHKERRQ(ierr);
	ierr = VecGetArray(v3,&LA_v3);CHKERRQ(ierr);
	for (i=0; i<n; i++) {
		LA_v3[i] = LA_v3[i] + LA_v2[i];
	}
	ierr = VecRestoreArray(v2,&LA_v2);CHKERRQ(ierr);
	ierr = VecRestoreArray(v3,&LA_v3);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetSubMatrix_MFStokes_A"
PetscErrorCode MatGetSubMatrix_MFStokes_A(Mat A,IS isr,IS isc,MatReuse cll,Mat *B)
{
	MatStokesMF ctx,copyA;
	MatA11MF copyA11;
	PetscBool f1,f2,isFullCol,same;
	PetscBool is_Auu_ii_mf = PETSC_FALSE;
	PetscBool is_Auu_mf = PETSC_FALSE;
	PetscInt ii,d;
	IS is_list[3];
	char *label[] = { "UU", "VV", "WW" };
	char *prefix = NULL;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	
	same = PETSC_FALSE;
	ierr = ISEqual(isr,isc,&same);CHKERRQ(ierr);

	ierr = MatGetOptionsPrefix(A,(const char**)&prefix);CHKERRQ(ierr);
	
	if (same) {
		/* either {u,v,w}, u,v,w,P */
		
		// [1]
		is_Auu_mf = PETSC_FALSE;
		ierr = PetscOptionsGetBool(prefix,"-A11_mf",&is_Auu_mf,0);CHKERRQ(ierr);
		
		f1 = f2 = PETSC_FALSE;
		ierr = ISEqual(isr,ctx->isUVW,&f1);CHKERRQ(ierr);
		ierr = ISEqual(isc,ctx->isUVW,&f2);CHKERRQ(ierr);
		if ((f1==PETSC_TRUE) && (f2==PETSC_TRUE)) {
			PetscPrintf(PETSC_COMM_WORLD,"Fetching UVW block = A(1,1)\n");
			if (!is_Auu_mf) {
				if (cll==MAT_INITIAL_MATRIX) {
					ierr = DMGetMatrix(ctx->daUVW,MATAIJ,B);CHKERRQ(ierr);
				} else {
					ierr = MatZeroEntries(*B);CHKERRQ(ierr);
				}
				PetscPrintf(PETSC_COMM_WORLD,"ierr = AssembleAUU_Stokes();CHKERRQ(ierr);  TODO \n");
				
			} else {
				if (cll==MAT_INITIAL_MATRIX) {
					PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
					ierr = MatCopy_StokesMF_A11MF(ctx,&copyA11);CHKERRQ(ierr);
					ierr = StokesQ2P1CreateMatrix_MFOperator_A11(copyA11,B);CHKERRQ(ierr);
					ierr = MatA11MFDestroy(&copyA11);CHKERRQ(ierr);
				} else {
					// to nothing
				}
			}
			PetscFunctionReturn(0);
		}
		
		// [2]
		is_Auu_ii_mf = PETSC_FALSE;
		//ierr = PetscOptionsGetBool(prefix,"-A11iimf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
		//ierr = PetscOptionsGetBool(PETSC_NULL,"-Auu_ii_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
		
		is_list[0] = ctx->isU;
		is_list[1] = ctx->isV;
		is_list[2] = ctx->isW;
		
		for (d=0; d<3; d++) {
			f1 = f2 = PETSC_FALSE;
			ierr = ISEqual(isr,is_list[d],&f1);CHKERRQ(ierr);
			ierr = ISEqual(isc,is_list[d],&f2);CHKERRQ(ierr);

			is_Auu_ii_mf = PETSC_FALSE;
			switch (d) {
				case 0:
					ierr = PetscOptionsGetBool(prefix,"-A11_Auu_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
				case 1:
					ierr = PetscOptionsGetBool(prefix,"-A11_Avv_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
				case 2:
					ierr = PetscOptionsGetBool(prefix,"-A11_Aww_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
			}			
			
			if ((f1==PETSC_TRUE) && (f2==PETSC_TRUE)) {
				PetscPrintf(PETSC_COMM_WORLD,"Fetching %s block = A(1,1)_%s\n",label[d],label[d]);
				
				if (!is_Auu_ii_mf) {
					if (cll==MAT_INITIAL_MATRIX) {
						ierr = DMGetMatrix(ctx->daU,MATAIJ,B);CHKERRQ(ierr);
					} else {
						ierr = MatZeroEntries(*B);CHKERRQ(ierr);
					}
					PetscPrintf(PETSC_COMM_WORLD,"ierr = AssembleAUiUi_Stokes();CHKERRQ(ierr);  TODO\n");
				} else {
					if (cll==MAT_INITIAL_MATRIX) {
						PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
						PetscPrintf(PETSC_COMM_WORLD,"ierr = MatCreateAUiUiStokesMF(A,d,B);CHKERRQ(ierr);  TODO \n");
					} else {
						// to nothing
					}
				}
				PetscFunctionReturn(0);
			}
		}		
		
		// [3]
		f1 = f2 = PETSC_FALSE;
		ierr = ISEqual(isr,ctx->isP,&f1);CHKERRQ(ierr);
		ierr = ISEqual(isc,ctx->isP,&f2);CHKERRQ(ierr);
		if ((f1==PETSC_TRUE) && (f2==PETSC_TRUE)) {
			
			//
			PetscPrintf(PETSC_COMM_WORLD,"Fetching PP stab block = A(2,2)\n");
			PetscPrintf(PETSC_COMM_WORLD,"  ***** !!!WARNING!!! Stokes operator A has been requested to fetch the NULL block A(2,2), returning the pressure mass matrix *****\n");
			
			if (cll==MAT_INITIAL_MATRIX) {
				ierr = DMGetMatrix(ctx->dap,MATAIJ,B);CHKERRQ(ierr);
			} else {
				ierr = MatZeroEntries(*B);CHKERRQ(ierr);
			}
			PetscPrintf(PETSC_COMM_WORLD,"ierr = AssembleAPP();CHKERRQ(ierr);  TODO\n");
			
			PetscFunctionReturn(0);
		}
		
	} else {
		/* either [U,P], [P,U] or [Ui,{U,P}] */
		
		// [1]
		f1 = f2 = PETSC_FALSE;
		ierr = ISEqual(isr,ctx->isUVW,&f1);CHKERRQ(ierr);
		ierr = ISEqual(isc,ctx->isP,&f2);CHKERRQ(ierr);
		if ( (f1==PETSC_TRUE) && (f2==PETSC_TRUE) ) {
			PetscPrintf(PETSC_COMM_WORLD,"Fetching UP block\n");
			
			if (cll==MAT_INITIAL_MATRIX) {
				PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
				ierr = MatStokesMFCopy(ctx,&copyA);CHKERRQ(ierr);
				ierr = StokesQ2P1CreateMatrix_MFOperator_A12(copyA,B);CHKERRQ(ierr);
				ierr = MatStokesMFDestroy(&copyA);CHKERRQ(ierr);
			} else {
				// to nothing
			}
			PetscFunctionReturn(0);
		}
		
		// [2]
		f1 = f2 = PETSC_FALSE;
		ierr = ISEqual(isr,ctx->isP,&f1);CHKERRQ(ierr);
		ierr = ISEqual(isc,ctx->isUVW,&f2);CHKERRQ(ierr);
		if ( (f1==PETSC_TRUE) && (f2==PETSC_TRUE) ) {
			PetscPrintf(PETSC_COMM_WORLD,"Fetching PU block\n");
			
			if (cll==MAT_INITIAL_MATRIX) {
				PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
				ierr = MatStokesMFCopy(ctx,&copyA);CHKERRQ(ierr);
				ierr = StokesQ2P1CreateMatrix_MFOperator_A21(copyA,B);CHKERRQ(ierr);
				ierr = MatStokesMFDestroy(&copyA);CHKERRQ(ierr);
			} else {
				// to nothing
			}
			PetscFunctionReturn(0);
		}
		
		// [3]
		/* check if full column space */
		isFullCol = PETSC_FALSE;
		same = PETSC_FALSE;
		ierr = PetscTypeCompare((PetscObject)isc,"stride",&same);CHKERRQ(ierr);
		if (same) {
			PetscInt n,first,step;
			ierr = ISStrideGetInfo(isc,&first,&step);CHKERRQ(ierr);
			ierr = ISGetLocalSize(isc,&n);CHKERRQ(ierr);
			if  ((A->cmap->n==n)
					 && (A->cmap->rstart==first)
					 && (step==1) ) {
				isFullCol = PETSC_TRUE;
				PetscPrintf(PETSC_COMM_WORLD,"Detected full column space\n");
				
				ii = -1;
				f1 = f2 = PETSC_FALSE;
				ierr = ISEqual(isr,ctx->isUVW,&f1);CHKERRQ(ierr);
				if (f1==PETSC_TRUE) { ii = 0; }
				ierr = ISEqual(isr,ctx->isP,&f2);CHKERRQ(ierr);
				if (f2==PETSC_TRUE) { ii = 1; }
				
				if (ii==-1) {
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Requested row which matrix-free operator doesn't define");
				}
				
				if (cll==MAT_INITIAL_MATRIX) {
					PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
					PetscPrintf(PETSC_COMM_WORLD,"ierr = MatCreateAiStokesMF(A,ii,B);CHKERRQ(ierr); TODO\n");
				} else {
					// to nothing
				}
				
				PetscFunctionReturn(0);
			}
		}
		
		
	}
	
	
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Something bad happened... I don't know how to retrieve the sub matrix you requested");
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetSubMatrix_MFStokes_A11"
PetscErrorCode MatGetSubMatrix_MFStokes_A11(Mat A,IS isr,IS isc,MatReuse cll,Mat *B)
{
	MatA11MF  ctx,copy;
	PetscBool f1,f2,isFullCol,same;
	PetscBool is_Auu_ii_mf = PETSC_FALSE;
	PetscBool is_Auu_mf = PETSC_FALSE;
	PetscInt ii,d;
	IS is_list[3];
	char *label[] = { "UU", "VV", "WW" };
	char *prefix = NULL;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	
	same = PETSC_FALSE;
	ierr = ISEqual(isr,isc,&same);CHKERRQ(ierr);
	
	ierr = MatGetOptionsPrefix(A,(const char**)&prefix);CHKERRQ(ierr);
	
	if (same) {
		/* either {u,v,w} */
		
		// [2]
		is_Auu_ii_mf = PETSC_FALSE;
		//ierr = PetscOptionsGetBool(prefix,"-A11iimf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
		//ierr = PetscOptionsGetBool(PETSC_NULL,"-Auu_ii_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
		
		is_list[0] = ctx->isU;
		is_list[1] = ctx->isV;
		is_list[2] = ctx->isW;
		
		for (d=0; d<3; d++) {
			f1 = f2 = PETSC_FALSE;
			ierr = ISEqual(isr,is_list[d],&f1);CHKERRQ(ierr);
			ierr = ISEqual(isc,is_list[d],&f2);CHKERRQ(ierr);
			
			is_Auu_ii_mf = PETSC_FALSE;
			switch (d) {
				case 0:
					ierr = PetscOptionsGetBool(prefix,"-Auu_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
				case 1:
					ierr = PetscOptionsGetBool(prefix,"-Avv_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
				case 2:
					ierr = PetscOptionsGetBool(prefix,"-Aww_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
			}			
			
			if ((f1==PETSC_TRUE) && (f2==PETSC_TRUE)) {
				PetscPrintf(PETSC_COMM_WORLD,"Fetching %s block = A(1,1)_%s\n",label[d],label[d]);
				
				if (!is_Auu_ii_mf) {
					if (cll==MAT_INITIAL_MATRIX) {
						ierr = DMGetMatrix(ctx->daU,MATAIJ,B);CHKERRQ(ierr);
					} else {
						ierr = MatZeroEntries(*B);CHKERRQ(ierr);
					}
					PetscPrintf(PETSC_COMM_WORLD,"ierr = AssembleAUiUi_Stokes();CHKERRQ(ierr);  TODO\n");
				} else {
					if (cll==MAT_INITIAL_MATRIX) {
						PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
						PetscPrintf(PETSC_COMM_WORLD,"ierr = MatCreateAUiUiStokesMF(A,d,B);CHKERRQ(ierr);  TODO \n");
					} else {
						// to nothing
					}
				}
				PetscFunctionReturn(0);
			}
		}		
		
	} else {
		// [3]
		/* check if full column space */
		isFullCol = PETSC_FALSE;
		same = PETSC_FALSE;
		ierr = PetscTypeCompare((PetscObject)isc,"stride",&same);CHKERRQ(ierr);
		if (same) {
			PetscInt n,first,step;
			ierr = ISStrideGetInfo(isc,&first,&step);CHKERRQ(ierr);
			ierr = ISGetLocalSize(isc,&n);CHKERRQ(ierr);
			if  ((A->cmap->n==n)
					 && (A->cmap->rstart==first)
					 && (step==1) ) {
				isFullCol = PETSC_TRUE;
				PetscPrintf(PETSC_COMM_WORLD,"Detected full column space\n");
				
				is_list[0] = ctx->isU;
				is_list[1] = ctx->isV;
				is_list[2] = ctx->isW;
				
				f1 = PETSC_FALSE;
				for (d=0; d<3; d++) {
					ierr = ISEqual(isr,is_list[d],&f1);CHKERRQ(ierr);
					if (f1==PETSC_TRUE){ break; }
				}
				
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Requested row which matrix-free operator doesn't define - need to determine which uu,vv,ww component requested");
				
				if (cll==MAT_INITIAL_MATRIX) {
					PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
					PetscPrintf(PETSC_COMM_WORLD,"ierr = MatCreateAiStokesMF(A,ii,B);CHKERRQ(ierr); TODO\n");
				} else {
					// to nothing
				}
				
				PetscFunctionReturn(0);
			}
		}
		
		
	}
	
	
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Something bad happened... I don't know how to retrieve the sub matrix you requested");
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A"
PetscErrorCode MatMult_MFStokes_A(Mat A,Vec X,Vec Y)
{
	MatStokesMF       ctx;
  PetscErrorCode    ierr;
  DM                stokes_pack,dau,dap;
  DMDALocalInfo     infou,infop;
  Vec               XUloc,XPloc,YUloc,YPloc;
	Vec               Xu,Xp,Yu,Yp;
  PetscScalar       *LA_XUloc,*LA_XPloc;
  PetscScalar       *LA_YUloc,*LA_YPloc;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	stokes_pack = ctx->stokes_pack;
	
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dau,&infou);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dap,&infop);CHKERRQ(ierr);
	
  ierr = DMCompositeGetLocalVectors(stokes_pack,&XUloc,&XPloc);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&YUloc,&YPloc);CHKERRQ(ierr);
	
	/* get the local (ghosted) entries for each physics */
	ierr = DMCompositeScatter(stokes_pack,X,XUloc,XPloc);CHKERRQ(ierr);
	
	/* Zero entries in local vectors corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out columns when the mat-mult is performed */
	ierr = BCListInsertLocalZero(ctx->u_bclist,XUloc);CHKERRQ(ierr);
	/* if we have pressure boundary conditions */
	/*
	 ierr = BCListInsertLocalZero(ctx->p_bclist,XPloc);CHKERRQ(ierr);
	 */
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	ierr = VecGetArray(XPloc,&LA_XPloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YUloc);CHKERRQ(ierr);
	ierr = VecZeroEntries(YPloc);CHKERRQ(ierr);
	ierr = VecGetArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	ierr = VecGetArray(YPloc,&LA_YPloc);CHKERRQ(ierr);
	
	/* momentum + continuity */
	//ierr = MF_Stokes_yAx(ctx,dau,LA_XUloc,dap,LA_XPloc,LA_YUloc,LA_YPloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(YPloc,&LA_YPloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XPloc,&LA_XPloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = VecZeroEntries(Y);CHKERRQ(ierr);
	ierr = DMCompositeGather(stokes_pack,Y,ADD_VALUES,YUloc,YPloc);CHKERRQ(ierr);
	
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&YUloc,&YPloc);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&XUloc,&XPloc);CHKERRQ(ierr);
	
	/* modify Y for the boundary conditions, y_k = scale_k(x_k) */
	ierr = DMCompositeGetAccess(stokes_pack,Y,&Yu,&Yp);CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
	
	/* Clobbering entries in global vector corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out rows when the mat-mult is performed */
	ierr = BCListInsertDirichlet_MatMult(ctx->u_bclist,Xu,Yu);CHKERRQ(ierr);
	/* if we have pressure boundary conditions */
	/*
	 ierr = BCListInsertDirichlet_MatMult(ctx->p_bclist,Xp,Yp);CHKERRQ(ierr);
	 */
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(stokes_pack,Y,&Yu,&Yp);CHKERRQ(ierr);
	
	{
		PetscPrintf(PETSC_COMM_WORLD,"%s: DUMMY MAT MULT FOR PLUMBING TESTING \n", __FUNCT__);
		ierr = VecCopy(X,Y);CHKERRQ(ierr);
	}
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A11"
PetscErrorCode MatMult_MFStokes_A11(Mat A,Vec X,Vec Y)
{
	MatA11MF          ctx;
  PetscErrorCode    ierr;
  DM                dau;
  DMDALocalInfo     infou;
  Vec               XUloc,YUloc;
  PetscScalar       *LA_XUloc;
  PetscScalar       *LA_YUloc;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	
  ierr = DMDAGetLocalInfo(dau,&infou);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&XUloc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dau,&YUloc);CHKERRQ(ierr);
	
	/* get the local (ghosted) entries for each physics */
	ierr = DMGlobalToLocalBegin(dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	
	/* Zero entries in local vectors corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out columns when the mat-mult is performed */
	ierr = BCListInsertLocalZero(ctx->u_bclist,XUloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YUloc);CHKERRQ(ierr);
	ierr = VecGetArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	
	/* momentum + continuity */
	//ierr = MF_Stokes_yAx(ctx,dau,LA_XUloc,LA_YUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = VecZeroEntries(Y);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(dau,YUloc,ADD_VALUES,Y);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (dau,YUloc,ADD_VALUES,Y);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(dau,&YUloc);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dau,&XUloc);CHKERRQ(ierr);
	
	/* modify Y for the boundary conditions, y_k = scale_k(x_k) */
	/* Clobbering entries in global vector corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out rows when the mat-mult is performed */
	ierr = BCListInsertDirichlet_MatMult(ctx->u_bclist,X,Y);CHKERRQ(ierr);
	
	{
		PetscPrintf(PETSC_COMM_WORLD,"%s: DUMMY MAT MULT FOR PLUMBING TESTING \n", __FUNCT__);
		ierr = VecCopy(X,Y);CHKERRQ(ierr);
	}
  PetscFunctionReturn(0);
}

/* 
 IN:  X - a pressure vector 
 OUT: Y - a velocity vector 
*/
#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A12"
PetscErrorCode MatMult_MFStokes_A12(Mat A,Vec X,Vec Y)
{
	MatStokesMF       ctx;
  PetscErrorCode    ierr;
  DM                dau,dap;
  DMDALocalInfo     infou,infop;
  Vec               XPloc,YUloc;
	Vec               Xp,Yu;
  PetscScalar       *LA_XPloc;
  PetscScalar       *LA_YUloc;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	dap = ctx->dap;
	
  ierr = DMDAGetLocalInfo(dau,&infou);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dap,&infop);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dap,&XPloc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dau,&YUloc);CHKERRQ(ierr);
	
	/* get the local (ghosted) entries for each physics */
	ierr = DMGlobalToLocalBegin(dap,X,INSERT_VALUES,XPloc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (dap,X,INSERT_VALUES,XPloc);CHKERRQ(ierr);
	
	/* Zero entries in local vectors corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out columns when the mat-mult is performed */
	/* if we have pressure boundary conditions */
	/*
	 ierr = BCListInsertLocalZero(ctx->p_bclist,XPloc);CHKERRQ(ierr);
	 */
	
	ierr = VecGetArray(XPloc,&LA_XPloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YUloc);CHKERRQ(ierr);
	ierr = VecGetArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	
	/* momentum + continuity */
	//ierr = MF_Stokes_yAx(ctx,dau,LA_XUloc,LA_YUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XPloc,&LA_XPloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = VecZeroEntries(Y);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(dau,YUloc,ADD_VALUES,Y);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (dau,YUloc,ADD_VALUES,Y);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(dap,&XPloc);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dau,&YUloc);CHKERRQ(ierr);
	
	/* modify Y for the boundary conditions, y_k = 0, 0 is inserted as this is an off-diagonal operator */
	/* Clobbering entries in global vector corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out rows when the mat-mult is performed */
	ierr = BCListInsertZero(ctx->u_bclist,Y);CHKERRQ(ierr);
	
	{
		PetscPrintf(PETSC_COMM_WORLD,"%s: DUMMY MAT MULT FOR PLUMBING TESTING \n", __FUNCT__);
		ierr = VecSet(Y,1.0);CHKERRQ(ierr);
	}
  PetscFunctionReturn(0);
}

/* 
 IN:  X - a velocity vector 
 OUT: Y - a pressure vector 
 */
#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A21"
PetscErrorCode MatMult_MFStokes_A21(Mat A,Vec X,Vec Y)
{
	MatStokesMF       ctx;
  PetscErrorCode    ierr;
  DM                dau,dap;
  DMDALocalInfo     infou,infop;
  Vec               XUloc,YPloc;
	Vec               Xu,Yp;
  PetscScalar       *LA_XUloc;
  PetscScalar       *LA_YPloc;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	dap = ctx->dap;
	
  ierr = DMDAGetLocalInfo(dau,&infou);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dap,&infop);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&XUloc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dap,&YPloc);CHKERRQ(ierr);
	
	/* get the local (ghosted) entries for each physics */
	ierr = DMGlobalToLocalBegin(dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	
	/* Zero entries in local vectors corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out columns when the mat-mult is performed */
	ierr = BCListInsertLocalZero(ctx->u_bclist,XUloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YPloc);CHKERRQ(ierr);
	ierr = VecGetArray(YPloc,&LA_YPloc);CHKERRQ(ierr);
	
	/* momentum + continuity */
	//ierr = MF_Stokes_yAx(ctx,dau,LA_XUloc,LA_YUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(YPloc,&LA_YPloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = VecZeroEntries(Y);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(dap,YPloc,ADD_VALUES,Y);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (dap,YPloc,ADD_VALUES,Y);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(dap,&YPloc);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dau,&XUloc);CHKERRQ(ierr);
	
	/* modify Y for the boundary conditions, y_k = 0, 0 is inserted as this is an off-diagonal operator */
	/* Clobbering entries in global vector corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out rows when the mat-mult is performed */
	/* if we have pressure boundary conditions */
	/*
	ierr = BCListInsertZero(ctx->p_bclist,Y);CHKERRQ(ierr);
	 */
	
	{
		PetscPrintf(PETSC_COMM_WORLD,"%s: DUMMY MAT MULT FOR PLUMBING TESTING \n", __FUNCT__);
		ierr = VecSet(Y,1.0);CHKERRQ(ierr);
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A(MatStokesMF Stk,Mat *A11)
{
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	Stk->refcnt++;
	ierr = MatCreateShell(PETSC_COMM_WORLD,Stk->mu+Stk->mp,Stk->mu+Stk->mp,Stk->Mu+Stk->Mp,Stk->Mu+Stk->Mp,(void*)Stk,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_GET_SUBMATRIX,(void(*)(void))MatGetSubMatrix_MFStokes_A);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatStokesMF);CHKERRQ(ierr);
	
	*A11 = B;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A11"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11(MatA11MF A11,Mat *A)
{
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	A11->refcnt++;
	ierr = MatCreateShell(PETSC_COMM_WORLD,A11->mu,A11->mu,A11->Mu,A11->Mu,(void*)A11,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A11);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_GET_SUBMATRIX,(void(*)(void))MatGetSubMatrix_MFStokes_A11);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatA11MF);CHKERRQ(ierr);
	ierr = MatSetBlockSize(B,3);CHKERRQ(ierr);
	
	*A = B;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A12"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A12(MatStokesMF Stk,Mat *A12)
{
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	Stk->refcnt++;
	ierr = MatCreateShell(PETSC_COMM_WORLD,Stk->mu,Stk->mp,Stk->Mu,Stk->Mp,(void*)Stk,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A12);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatStokesMF);CHKERRQ(ierr);

	*A12 = B;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A21"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A21(MatStokesMF Stk,Mat *A21)
{
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	Stk->refcnt++;
	ierr = MatCreateShell(PETSC_COMM_WORLD,Stk->mp,Stk->mu,Stk->Mp,Stk->Mu,(void*)Stk,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A21);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatStokesMF);CHKERRQ(ierr);
	
	*A21 = B;

	PetscFunctionReturn(0);
}

/* 
 Should be
 PetscErrorCode StokesQ2P1CreateMatrix_MFOperator(PhysCompStokes user,Mat *B)
 */
#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_Operator"
PetscErrorCode StokesQ2P1CreateMatrix_Operator(PhysCompStokes user,Mat *B)
{
	PetscBool      same;
  DM             dau,dap,pack;
	Mat            A;
	MPI_Comm       comm;
	PetscInt       mu,Mu,mp,Mp;
	Vec            X,u,p;
	MatStokesMF    StkCtx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	pack = user->stokes_pack;

	/* is composite */
	same = PETSC_FALSE;
	ierr = PetscTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
	/* Create submatrices */
	comm = ((PetscObject)pack)->comm;
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A(StkCtx,&A);CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(A,"stokes_Amf_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	
	*B = A;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrixNest_Operator"
PetscErrorCode StokesQ2P1CreateMatrixNest_Operator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B)
{
	PetscBool      same;
  DM             dau,dap,pack;
	Mat            A,Auu,Aup,Apu,bA[2][2];
	IS             *is;
	MPI_Comm       comm;
	PetscInt       mu,Mu,mp,Mp,i,j;
	Vec            X,u,p;
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);

	pack = user->stokes_pack;

	
	/* is composite */
	same = PETSC_FALSE;
	ierr = PetscTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
  /* Fetch the DA's */
  ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	/* Create submatrices */
	comm = ((PetscObject)pack)->comm;

	/* A11 */
	if (tA11==0) {
		ierr = DMGetMatrix(dau,MATAIJ,&Auu);CHKERRQ(ierr);
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
		ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
	}
	ierr = MatSetOptionsPrefix(Auu,"stokes_A_A11_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);
	
	/* A12 */
	if (tA12==0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stokes->A12 doesn't support assembled operator");
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A12(StkCtx,&Aup);CHKERRQ(ierr);
		ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	}
	ierr = MatSetOptionsPrefix(Aup,"stokes_A_A12_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Aup);CHKERRQ(ierr);

	/* A21 */
	if (tA21==0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stokes->A21 doesn't support assembled operator");
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A21(StkCtx,&Apu);CHKERRQ(ierr); 
		ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	}
	ierr = MatSetOptionsPrefix(Apu,"stokes_A_A21_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Apu);CHKERRQ(ierr);
	
	/* Create nest */
	ierr = DMCompositeGetGlobalISs(pack,&is);CHKERRQ(ierr);
	
	bA[0][0] = Auu; bA[0][1] = Aup;
	bA[1][0] = Apu; bA[1][1] = PETSC_NULL;
  ierr = MatCreateNest(PETSC_COMM_WORLD,2,is,2,is,&bA[0][0],&A);CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(A,"stokes_A_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	*B = A;
	
	/* tidy up */
	for (i=0; i<2; i++) {
		for (j=0; j<2; j++) {
			if (bA[i][j]) { ierr = MatDestroy(&bA[i][j]);CHKERRQ(ierr); }
		}
	}
	ierr = ISDestroy(&is[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is[1]);CHKERRQ(ierr);
	ierr = PetscFree(is);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrixNest_PCOperator"
PetscErrorCode StokesQ2P1CreateMatrixNest_PCOperator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B)
{
	PetscBool      same;
  DM             dau,dap,pack;
	Mat            A,Auu,Aup,Apu,Spp,bA[2][2];
	IS             *is;
	MPI_Comm       comm;
	PetscInt       mu,Mu,mp,Mp,i,j;
	MatStokesMF    StkCtx;
	MatA11MF       A11Ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	ierr = MatCopy_StokesMF_A11MF(StkCtx,&A11Ctx);CHKERRQ(ierr);
	pack = user->stokes_pack;
	
	/* is composite */
	same = PETSC_FALSE;
	ierr = PetscTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
  /* Fetch the DA's */
  ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	/* Create submatrices */
	comm = ((PetscObject)pack)->comm;
	
	/* A11 */
	if (tA11==0) {
		ierr = DMGetMatrix(dau,MATAIJ,&Auu);CHKERRQ(ierr);  ierr = MatSetOptionsPrefix(Auu,"Buu");CHKERRQ(ierr);
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
	}

	/* Schur complement */
	ierr = DMGetMatrix(dap,MATAIJ,&Spp);CHKERRQ(ierr);  ierr = MatSetOptionsPrefix(Spp,"S*");CHKERRQ(ierr);

	/* A12 */
	if (tA12==0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stokes->A12 doesn't support assembled operator");
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A12(StkCtx,&Aup);CHKERRQ(ierr); ierr = MatSetOptionsPrefix(Aup,"Bup");CHKERRQ(ierr);
	}
	
	/* A21 */
	if (tA21==0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stokes->A21 doesn't support assembled operator");
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A21(StkCtx,&Apu);CHKERRQ(ierr); ierr = MatSetOptionsPrefix(Apu,"Bpu");CHKERRQ(ierr);
	}
	
	/* Create nest */
	ierr = DMCompositeGetGlobalISs(pack,&is);CHKERRQ(ierr);
	
	bA[0][0] = Auu;        bA[0][1] = Aup;
	bA[1][0] = Apu;        bA[1][1] = Spp;
  ierr = MatCreateNest(PETSC_COMM_WORLD,2,is,2,is,&bA[0][0],&A);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	*B = A;
	
	/* tidy up */
	for (i=0; i<2; i++) {
		for (j=0; j<2; j++) {
			if (bA[i][j]) { ierr = MatDestroy(&bA[i][j]);CHKERRQ(ierr); }
		}
	}
	ierr = ISDestroy(&is[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is[1]);CHKERRQ(ierr);
	ierr = PetscFree(is);CHKERRQ(ierr);
	ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

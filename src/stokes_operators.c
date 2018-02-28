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
 **    filename:   stokes_operators.c
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
-

*/


#include "petsc/private/petscimpl.h"
#include "petsc.h" /*I   "petscmat.h"   I*/
#include "petscsys.h"
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "data_bucket.h"

#include "dmda_duplicate.h"
#include "dmda_bcs.h"
#include "quadrature.h"

#include "stokes_operators.h"
#include "stokes_operators_mf.h"
#include "stokes_assembly.h"

PetscLogEvent MAT_MultMFA11;
PetscLogEvent MAT_MultMFA11_stp;
PetscLogEvent MAT_MultMFA11_cto;
PetscLogEvent MAT_MultMFA11_ker;
PetscLogEvent MAT_MultMFA11_cfr;

PetscLogEvent MAT_MultMFA11_sub;
PetscLogEvent MAT_MultMFA11_rto;
PetscLogEvent MAT_MultMFA11_rfr;
PetscLogEvent MAT_MultMFA11_SUP;

PetscLogEvent MAT_MultMFA; 
PetscLogEvent MAT_MultMFA12; 
PetscLogEvent MAT_MultMFA21; 

PetscLogEvent MAT_MultMFA_QuasiNewtonX; 
PetscLogEvent MAT_MultMFA11_QuasiNewtonX; 
PetscLogEvent MAT_MultMFA12_QuasiNewtonX; 
PetscLogEvent MAT_MultMFA21_QuasiNewtonX; 

PetscLogEvent MAT_MultMFA11_gtl;
PetscLogEvent MAT_MultMFA11_ltg;
PetscLogEvent MAT_MultMFA11_yz;
PetscLogEvent MAT_MultMFA11_gv;

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
	PetscFunctionList MatMult_flist = NULL;
	PetscFunctionList SetUp_flist = NULL;
	PetscFunctionList Destroy_flist = NULL;
	PetscErrorCode ierr;
	MatA11MF A11;
	char optype[64] = "ref";
	PetscFunctionBegin;
	
	ierr = PetscMalloc(sizeof(struct _p_MatA11MF),&A11);CHKERRQ(ierr);
	ierr = PetscMemzero(A11,sizeof(struct _p_MatA11MF));CHKERRQ(ierr);

  A11->is_setup       = PETSC_FALSE;
  A11->ctx            = NULL;
  A11->SpMVOp_MatMult = NULL;
  A11->SpMVOp_SetUp   = NULL;
  A11->SpMVOp_Destroy = NULL;
  
	ierr = PetscFunctionListAdd(&MatMult_flist,"ref",MFStokesWrapper_A11);CHKERRQ(ierr);
	ierr = PetscFunctionListAdd(&MatMult_flist,"tensor",MFStokesWrapper_A11_Tensor);CHKERRQ(ierr);
#ifdef TATIN_HAVE_AVX
	ierr = PetscFunctionListAdd(&MatMult_flist,"avx",MFStokesWrapper_A11_AVX);CHKERRQ(ierr);
#endif
#ifdef TATIN_HAVE_CUDA
	ierr = PetscFunctionListAdd(&MatMult_flist,"cuda",MFStokesWrapper_A11_CUDA);CHKERRQ(ierr);
	ierr = PetscFunctionListAdd(&SetUp_flist,"cuda",MFA11SetUp_CUDA);CHKERRQ(ierr);
	ierr = PetscFunctionListAdd(&Destroy_flist,"cuda",MFA11Destroy_CUDA);CHKERRQ(ierr);
#endif
#ifdef TATIN_HAVE_OPENCL
	ierr = PetscFunctionListAdd(&MatMult_flist,"opencl",MFStokesWrapper_A11_OpenCL);CHKERRQ(ierr);
	ierr = PetscFunctionListAdd(&SetUp_flist,"opencl",MFA11SetUp_OpenCL);CHKERRQ(ierr);
	ierr = PetscFunctionListAdd(&Destroy_flist,"opencl",MFA11Destroy_OpenCL);CHKERRQ(ierr);
#endif
#if defined(TATIN_HAVE_AVX) && defined(TATIN_HAVE_CUDA)
	ierr = PetscFunctionListAdd(&MatMult_flist,"subrepart",MFStokesWrapper_A11_SubRepart);CHKERRQ(ierr);
	ierr = PetscFunctionListAdd(&SetUp_flist,"subrepart",MFA11SetUp_SubRepart);CHKERRQ(ierr);
	ierr = PetscFunctionListAdd(&Destroy_flist,"subrepart",MFA11Destroy_SubRepart);CHKERRQ(ierr);
#endif
	ierr = PetscOptionsGetString(NULL,NULL,"-a11_op",optype,sizeof optype,NULL);CHKERRQ(ierr);
	ierr = PetscFunctionListFind(MatMult_flist,optype,&A11->SpMVOp_MatMult);CHKERRQ(ierr);
	ierr = PetscFunctionListFind(SetUp_flist,optype,&A11->SpMVOp_SetUp);CHKERRQ(ierr);
	ierr = PetscFunctionListFind(Destroy_flist,optype,&A11->SpMVOp_Destroy);CHKERRQ(ierr);
	if (!A11->SpMVOp_MatMult) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"No -a11_op %s",optype);
	ierr = PetscFunctionListDestroy(&MatMult_flist);CHKERRQ(ierr);
	ierr = PetscFunctionListDestroy(&SetUp_flist);CHKERRQ(ierr);
	ierr = PetscFunctionListDestroy(&Destroy_flist);CHKERRQ(ierr);
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
	DM dau,pack;
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
	ierr = PetscObjectTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
    /* Fetch the DA's */
	dau = user->dav;
	//dap = user->dap;
	
	/* Sizes */
	ierr = DMGetGlobalVector(pack,&X);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(pack,X,&u,&p);CHKERRQ(ierr);
	ierr = VecGetSize(u,&Mu);CHKERRQ(ierr);
	ierr = VecGetLocalSize(u,&mu);CHKERRQ(ierr);
	ierr = VecGetSize(p,&Mp);CHKERRQ(ierr);
	ierr = VecGetLocalSize(p,&mp);CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(u,&start,NULL);CHKERRQ(ierr);
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
	ierr = ISCreateStride(PetscObjectComm((PetscObject)user->dav), n,offset,3,&StkCtx->isU);CHKERRQ(ierr);
	offset = start + 1;
	ierr = ISCreateStride(PetscObjectComm((PetscObject)user->dav), n,offset,3,&StkCtx->isV);CHKERRQ(ierr);
	offset = start + 2;
	ierr = ISCreateStride(PetscObjectComm((PetscObject)user->dav), n,offset,3,&StkCtx->isW);CHKERRQ(ierr);
	
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
	DM dau;
	PetscInt n,start,offset;
	
	PetscFunctionBegin;
	
  if (A11Ctx->is_setup) PetscFunctionReturn(0);
  
	A11Ctx->daUVW       = dav;           ierr = PetscObjectReference((PetscObject)dav);CHKERRQ(ierr); 
	A11Ctx->volQ        = volQ;
	A11Ctx->u_bclist    = u_bclist;
	
  /* Fetch the DA's */
	dau = dav;
	
	/* Sizes */
	ierr = DMGetGlobalVector(dau,&X);CHKERRQ(ierr);
	ierr = VecGetSize(X,&Mu);CHKERRQ(ierr);
	ierr = VecGetLocalSize(X,&mu);CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(X,&start,NULL);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(dau,&X);CHKERRQ(ierr);
	
	A11Ctx->mu = mu;
	A11Ctx->Mu = Mu;
	
	n = (mu/3);
	offset = start + 0;
	ierr = ISCreateStride(PetscObjectComm((PetscObject)dav), n,offset,3,&A11Ctx->isU);CHKERRQ(ierr);
	offset = start + 1;
	ierr = ISCreateStride(PetscObjectComm((PetscObject)dav), n,offset,3,&A11Ctx->isV);CHKERRQ(ierr);
	offset = start + 2;
	ierr = ISCreateStride(PetscObjectComm((PetscObject)dav), n,offset,3,&A11Ctx->isW);CHKERRQ(ierr);
	
	ierr = DMDADuplicateLayout(dau,1,2,DMDA_STENCIL_BOX,&A11Ctx->daU);CHKERRQ(ierr);
	
  if (A11Ctx->SpMVOp_SetUp) {
    ierr = A11Ctx->SpMVOp_SetUp(A11Ctx);CHKERRQ(ierr);
  }
  
	A11Ctx->refcnt = 1;
	A11Ctx->is_setup = PETSC_TRUE;
  
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
		if (A->daUVW) { ierr = DMDestroy(&A->daUVW);CHKERRQ(ierr); }
		if (A->dap) { ierr = DMDestroy(&A->dap);CHKERRQ(ierr); }
		if (A->stokes_pack) { ierr = DMDestroy(&A->stokes_pack);CHKERRQ(ierr); }
		
		if (A->isU) { ierr = ISDestroy(&A->isU);CHKERRQ(ierr); }
		if (A->isV) { ierr = ISDestroy(&A->isV);CHKERRQ(ierr); }
		if (A->isW) { ierr = ISDestroy(&A->isW);CHKERRQ(ierr); }
		if (A->isP) { ierr = ISDestroy(&A->isP);CHKERRQ(ierr); }
		if (A->isUVW) { ierr = ISDestroy(&A->isUVW);CHKERRQ(ierr); }
		if (A->daU) { ierr = DMDestroy(&A->daU);CHKERRQ(ierr); }
		ierr = PetscFree(A);CHKERRQ(ierr);

		*B = NULL;
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
		if (A->daUVW) { ierr = DMDestroy(&A->daUVW);CHKERRQ(ierr); }
		if (A->isUVW) { ierr = ISDestroy(&A->isUVW);CHKERRQ(ierr); }
		
		if (A->isU) { ierr = ISDestroy(&A->isU);CHKERRQ(ierr); }
		if (A->isV) { ierr = ISDestroy(&A->isV);CHKERRQ(ierr); }
		if (A->isW) { ierr = ISDestroy(&A->isW);CHKERRQ(ierr); }
		if (A->daU) { ierr = DMDestroy(&A->daU);CHKERRQ(ierr); }
    
    if (A->SpMVOp_Destroy) {
      ierr = A->SpMVOp_Destroy(A);CHKERRQ(ierr);
    }
    
		ierr = PetscFree(A);CHKERRQ(ierr);

		*B = NULL;
	}
	
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatDestroy_MatStokesMF"
PetscErrorCode MatDestroy_MatStokesMF(Mat A)
{
	MatStokesMF     ctx;
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
#define __FUNCT__ "MatDestroy_MatA11MF_QuasiNewtonX"
PetscErrorCode MatDestroy_MatA11MF_QuasiNewtonX(Mat A)
{
	MatA11MF        ctx;
	Vec             Xloc;
	PetscErrorCode  ierr;
	
  PetscFunctionBegin;
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);

	/* fetch shifted coordinate vector from Mat A */
	Xloc = NULL;
	ierr = PetscObjectQuery((PetscObject)A,"MatA11_QuasiNewtonX",(PetscObject*)&Xloc);CHKERRQ(ierr);
	if (!Xloc) { SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Require a shifted coordinate vector to have been set via PetscObjectCompose()"); }
	ierr = VecDestroy(&Xloc);CHKERRQ(ierr);
	
	ierr = MatA11MFDestroy(&ctx);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatDestroy_MatStokesMF_QuasiNewtonX"
PetscErrorCode MatDestroy_MatStokesMF_QuasiNewtonX(Mat A)
{
	MatStokesMF     ctx;
	Vec             Xloc;
	PetscErrorCode  ierr;
	
  PetscFunctionBegin;
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);

	/* fetch shifted coordinate vector from Mat A */
	Xloc = NULL;
	ierr = PetscObjectQuery((PetscObject)A,"MatA_QuasiNewtonX",(PetscObject*)&Xloc);CHKERRQ(ierr);
	if (!Xloc) { SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Require a shifted coordinate vector to have been set via PetscObjectCompose()"); }
	ierr = VecDestroy(&Xloc);CHKERRQ(ierr);
	
	ierr = MatStokesMFDestroy(&ctx);CHKERRQ(ierr);
	
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

	Stk->u_bclist = A->u_bclist;
	Stk->p_bclist = A->p_bclist;
	Stk->volQ     = A->volQ;
	
	Stk->stokes_pack  = A->stokes_pack;    ierr = PetscObjectReference((PetscObject)A->stokes_pack);CHKERRQ(ierr);
	
	Stk->isU   = A->isU;             if (A->isU) { ierr = PetscObjectReference((PetscObject)A->isU);CHKERRQ(ierr); }
	Stk->isV   = A->isV;             if (A->isV) { ierr = PetscObjectReference((PetscObject)A->isV);CHKERRQ(ierr); }
	Stk->isW   = A->isW;             if (A->isW) { ierr = PetscObjectReference((PetscObject)A->isW);CHKERRQ(ierr); }
	Stk->isP   = A->isP;             if (A->isP) { ierr = PetscObjectReference((PetscObject)A->isP);CHKERRQ(ierr); }
	Stk->isUVW = A->isUVW;           if (A->isUVW) { ierr = PetscObjectReference((PetscObject)A->isUVW);CHKERRQ(ierr); }
	
	Stk->daUVW = A->daUVW;        if (A->daUVW) { ierr = PetscObjectReference((PetscObject)A->daUVW);CHKERRQ(ierr); }
	Stk->daU   = A->daU;          if (A->daU) { ierr = PetscObjectReference((PetscObject)A->daU);CHKERRQ(ierr); }
	Stk->dap   = A->dap;          if (A->dap) { ierr = PetscObjectReference((PetscObject)A->dap);CHKERRQ(ierr); }
	
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
	
	A11->u_bclist = A->u_bclist;
	A11->volQ     = A->volQ;
	
	A11->isU   = A->isU;             if (A->isU) { ierr = PetscObjectReference((PetscObject)A->isU);CHKERRQ(ierr); }
	A11->isV   = A->isV;             if (A->isV) { ierr = PetscObjectReference((PetscObject)A->isV);CHKERRQ(ierr); }
	A11->isW   = A->isW;             if (A->isW) { ierr = PetscObjectReference((PetscObject)A->isW);CHKERRQ(ierr); }
	
	A11->daUVW = A->daUVW;        if (A->daUVW) { ierr = PetscObjectReference((PetscObject)A->daUVW);CHKERRQ(ierr); }
	A11->daU   = A->daU;          if (A->daU) { ierr = PetscObjectReference((PetscObject)A->daU);CHKERRQ(ierr); }
	
  if (A11->SpMVOp_SetUp) {
    ierr = A11->SpMVOp_SetUp(A11);CHKERRQ(ierr);
  }

	A11->refcnt = 1;
  A11->is_setup = PETSC_TRUE;

	*B = A11;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatCopy_StokesMF_A11MF"
PetscErrorCode MatCopy_StokesMF_A11MF(MatStokesMF A,MatA11MF *B)
{
	PetscErrorCode ierr;
	MatA11MF       A11;
	
	PetscFunctionBegin;
	ierr = MatA11MFCreate(&A11);CHKERRQ(ierr);

	A11->mu       = A->mu;
	A11->Mu       = A->Mu;
	A11->u_bclist = A->u_bclist;
	A11->volQ     = A->volQ;
	A11->isU      = A->isU;   if (A->isU)   { ierr = PetscObjectReference((PetscObject)A->isU);  CHKERRQ(ierr); }
	A11->isV      = A->isV;   if (A->isV)   { ierr = PetscObjectReference((PetscObject)A->isV);  CHKERRQ(ierr); }
	A11->isW      = A->isW;   if (A->isW)   { ierr = PetscObjectReference((PetscObject)A->isW);  CHKERRQ(ierr); }
	A11->daUVW    = A->daUVW; if (A->daUVW) { ierr = PetscObjectReference((PetscObject)A->daUVW);CHKERRQ(ierr); }
	A11->daU      = A->daU;   if (A->daU)   { ierr = PetscObjectReference((PetscObject)A->daU);  CHKERRQ(ierr); }
  if (A11->SpMVOp_SetUp) {
    ierr = A11->SpMVOp_SetUp(A11);CHKERRQ(ierr);
  }
	A11->refcnt   = 1;
	A11->is_setup = PETSC_TRUE;

	*B = A11;
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "MatMultAdd_basic"
PetscErrorCode MatMultAdd_basic(Mat A,Vec v1,Vec v2,Vec v3)
{
	PetscErrorCode ierr;
	PetscScalar *LA_v2,*LA_v3,*LA_tmp;
	PetscInt i,n;
	Vec tmp;
  PetscFunctionBegin;
	
	ierr = VecDuplicate(v2,&tmp);CHKERRQ(ierr);
	ierr = MatMult(A,v1,tmp);CHKERRQ(ierr);
	
	ierr = VecGetLocalSize(v2,&n);CHKERRQ(ierr);
	ierr = VecGetArray(v2,&LA_v2);CHKERRQ(ierr);
	ierr = VecGetArray(v3,&LA_v3);CHKERRQ(ierr);
	ierr = VecGetArray(tmp,&LA_tmp);CHKERRQ(ierr);
	for (i=0; i<n; i++) {
//		printf("  [%d] : v2 = %lf .. v3 = %lf \n", i,LA_tmp[i],LA_v2[i]);
		LA_v3[i] = LA_tmp[i] + LA_v2[i];
	}
	ierr = VecRestoreArray(v3,&LA_v3);CHKERRQ(ierr);
	ierr = VecRestoreArray(v2,&LA_v2);CHKERRQ(ierr);
	ierr = VecRestoreArray(tmp,&LA_tmp);CHKERRQ(ierr);
	ierr = VecDestroy(&tmp);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMultTransposeAdd_generic"
PetscErrorCode MatMultTransposeAdd_generic(Mat mat,Vec v1,Vec v2,Vec v3)
{
	Vec vt;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = VecDuplicate(v1,&vt);CHKERRQ(ierr);
	ierr = VecCopy(v1,vt);CHKERRQ(ierr);
	
	ierr = MatMult(mat,v1,vt);CHKERRQ(ierr);
	ierr = VecAXPY(vt,1.0,v2);
	ierr = VecCopy(vt,v3);CHKERRQ(ierr);
	ierr = VecDestroy(&vt);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
 
 
*/
#undef __FUNCT__  
#define __FUNCT__ "MatGetSubMatrix_MFStokes_A"
PetscErrorCode MatGetSubMatrix_MFStokes_A(Mat A,IS isr,IS isc,MatReuse cll,Mat *B)
{
	MatStokesMF ctx,copyA;
	MatA11MF copyA11;
	PetscBool f1,f2,same;
	PetscBool is_Auu_ii_mf = PETSC_FALSE;
	PetscBool is_Auu_mf = PETSC_FALSE;
	PetscInt ii,d;
	IS is_list[3];
	const char *label[] = { "UU", "VV", "WW" };
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
		ierr = PetscOptionsGetBool(NULL,prefix,"-A11_mf",&is_Auu_mf,0);CHKERRQ(ierr);
		
		f1 = f2 = PETSC_FALSE;
		ierr = ISEqual(isr,ctx->isUVW,&f1);CHKERRQ(ierr);
		ierr = ISEqual(isc,ctx->isUVW,&f2);CHKERRQ(ierr);
		if ((f1==PETSC_TRUE) && (f2==PETSC_TRUE)) {
			PetscPrintf(PetscObjectComm((PetscObject)A),"Fetching UVW block = A(1,1)\n");
			if (!is_Auu_mf) {
				if (cll==MAT_INITIAL_MATRIX) {
					ierr = DMSetMatType(ctx->daUVW,MATAIJ);CHKERRQ(ierr);
					ierr = DMCreateMatrix(ctx->daUVW,B);CHKERRQ(ierr);
				} else {
					ierr = MatZeroEntries(*B);CHKERRQ(ierr);
				}
				//PetscPrintf(PETSC_COMM_WORLD,"ierr = AssembleAUU_Stokes();CHKERRQ(ierr);  TODO \n");
				ierr = MatAssemble_StokesA_AUU(*B,ctx->daUVW,ctx->u_bclist,ctx->volQ);CHKERRQ(ierr);
			} else {
				if (cll==MAT_INITIAL_MATRIX) {
					PetscPrintf(PetscObjectComm((PetscObject)A),"  defining matrix free operator\n");
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
		//ierr = PetscOptionsGetBool(NULL,prefix,"-A11iimf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
		//ierr = PetscOptionsGetBool(NULL,NULL,"-Auu_ii_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
		
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
					ierr = PetscOptionsGetBool(NULL,prefix,"-A11_Auu_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
				case 1:
					ierr = PetscOptionsGetBool(NULL,prefix,"-A11_Avv_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
				case 2:
					ierr = PetscOptionsGetBool(NULL,prefix,"-A11_Aww_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
			}			
			
			if ((f1==PETSC_TRUE) && (f2==PETSC_TRUE)) {
				PetscPrintf(PetscObjectComm((PetscObject)A),"Fetching %s block = A(1,1)_%s\n",label[d],label[d]);
				
				if (!is_Auu_ii_mf) {
					if (cll==MAT_INITIAL_MATRIX) {
						ierr = DMSetMatType(ctx->daU,MATAIJ);CHKERRQ(ierr);
						ierr = DMCreateMatrix(ctx->daU,B);CHKERRQ(ierr);
					} else {
						ierr = MatZeroEntries(*B);CHKERRQ(ierr);
					}
					PetscPrintf(PetscObjectComm((PetscObject)A),"ierr = AssembleAUiUi_Stokes();CHKERRQ(ierr);  TODO\n");
				} else {
					if (cll==MAT_INITIAL_MATRIX) {
						PetscPrintf(PetscObjectComm((PetscObject)A),"  defining matrix free operator\n");
						PetscPrintf(PetscObjectComm((PetscObject)A),"ierr = MatCreateAUiUiStokesMF(A,d,B);CHKERRQ(ierr);  TODO \n");
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
			PetscPrintf(PetscObjectComm((PetscObject)A),"Fetching PP stab block = A(2,2)\n");
			PetscPrintf(PetscObjectComm((PetscObject)A),"  ***** !!!WARNING!!! Stokes operator A has been requested to fetch the NULL block A(2,2), returning the pressure mass matrix *****\n");
			
			if (cll==MAT_INITIAL_MATRIX) {
				ierr = DMSetMatType(ctx->dap,MATAIJ);CHKERRQ(ierr);
				ierr = DMCreateMatrix(ctx->dap,B);CHKERRQ(ierr);
			} else {
				ierr = MatZeroEntries(*B);CHKERRQ(ierr);
			}
			PetscPrintf(PetscObjectComm((PetscObject)A),"ierr = AssembleAPP();CHKERRQ(ierr);  TODO\n");
			
			PetscFunctionReturn(0);
		}
		
	} else {
		/* either [U,P], [P,U] or [Ui,{U,P}] */
		
		// [1]
		f1 = f2 = PETSC_FALSE;
		ierr = ISEqual(isr,ctx->isUVW,&f1);CHKERRQ(ierr);
		ierr = ISEqual(isc,ctx->isP,&f2);CHKERRQ(ierr);
		if ( (f1==PETSC_TRUE) && (f2==PETSC_TRUE) ) {
			//PetscPrintf(PETSC_COMM_WORLD,"Fetching UP block\n");
			
			if (cll==MAT_INITIAL_MATRIX) {
				//PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
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
			//PetscPrintf(PETSC_COMM_WORLD,"Fetching PU block\n");
			
			if (cll==MAT_INITIAL_MATRIX) {
				//PetscPrintf(PETSC_COMM_WORLD,"  defining matrix free operator\n");
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
		//isFullCol = PETSC_FALSE;
		same = PETSC_FALSE;
		ierr = PetscObjectTypeCompare((PetscObject)isc,"stride",&same);CHKERRQ(ierr);
		if (same) {
			PetscInt nc,n,cstart,first,step;
			ierr = ISStrideGetInfo(isc,&first,&step);CHKERRQ(ierr);
			ierr = ISGetLocalSize(isc,&n);CHKERRQ(ierr);
            ierr = MatGetLocalSize(A,NULL,&nc);CHKERRQ(ierr);
            ierr = MatGetOwnershipRange(A,&cstart,NULL);CHKERRQ(ierr);
			if  (    (nc == n)
					 && (cstart == first)
					 && (step==1) ) {
				//isFullCol = PETSC_TRUE;
				PetscPrintf(PetscObjectComm((PetscObject)A),"Detected full column space\n");
				
				ii = -1;
				f1 = f2 = PETSC_FALSE;
				ierr = ISEqual(isr,ctx->isUVW,&f1);CHKERRQ(ierr);
				if (f1==PETSC_TRUE) { ii = 0; }
				ierr = ISEqual(isr,ctx->isP,&f2);CHKERRQ(ierr);
				if (f2==PETSC_TRUE) { ii = 1; }
				
				if (ii==-1) {
					SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Requested row which matrix-free operator doesn't define");
				}
				
				if (cll==MAT_INITIAL_MATRIX) {
					PetscPrintf(PetscObjectComm((PetscObject)A),"  defining matrix free operator\n");
					PetscPrintf(PetscObjectComm((PetscObject)A),"ierr = MatCreateAiStokesMF(A,ii,B);CHKERRQ(ierr); TODO\n");
				} else {
					// to nothing
				}
				
				PetscFunctionReturn(0);
			}
		}
		
		
	}
	
	
	SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Something bad happened... I don't know how to retrieve the sub matrix you requested");
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetSubMatrix_MFStokes_A11"
PetscErrorCode MatGetSubMatrix_MFStokes_A11(Mat A,IS isr,IS isc,MatReuse cll,Mat *B)
{
	MatA11MF  ctx;
	PetscBool f1,f2,same;
	PetscBool is_Auu_ii_mf = PETSC_FALSE;
	PetscInt d;
	IS is_list[3];
	const char *label[] = { "UU", "VV", "WW" };
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
		//ierr = PetscOptionsGetBool(NULL,prefix,"-A11iimf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
		//ierr = PetscOptionsGetBool(NULL,NULL,"-Auu_ii_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
		
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
					ierr = PetscOptionsGetBool(NULL,prefix,"-Auu_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
				case 1:
					ierr = PetscOptionsGetBool(NULL,prefix,"-Avv_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
				case 2:
					ierr = PetscOptionsGetBool(NULL,prefix,"-Aww_mf",&is_Auu_ii_mf,0);CHKERRQ(ierr);
					break;
			}			
			
			if ((f1==PETSC_TRUE) && (f2==PETSC_TRUE)) {
				PetscPrintf(PetscObjectComm((PetscObject)A),"Fetching %s block = A(1,1)_%s\n",label[d],label[d]);
				
				if (!is_Auu_ii_mf) {
					if (cll==MAT_INITIAL_MATRIX) {
						ierr = DMSetMatType(ctx->daU,MATAIJ);CHKERRQ(ierr);
						ierr = DMCreateMatrix(ctx->daU,B);CHKERRQ(ierr);
					} else {
						ierr = MatZeroEntries(*B);CHKERRQ(ierr);
					}
					PetscPrintf(PetscObjectComm((PetscObject)A),"ierr = AssembleAUiUi_Stokes();CHKERRQ(ierr);  TODO\n");
				} else {
					if (cll==MAT_INITIAL_MATRIX) {
						PetscPrintf(PetscObjectComm((PetscObject)A),"  defining matrix free operator\n");
						PetscPrintf(PetscObjectComm((PetscObject)A),"ierr = MatCreateAUiUiStokesMF(A,d,B);CHKERRQ(ierr);  TODO \n");
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
		//isFullCol = PETSC_FALSE;
		same = PETSC_FALSE;
		ierr = PetscObjectTypeCompare((PetscObject)isc,"stride",&same);CHKERRQ(ierr);
		if (same) {
			PetscInt nc,n,cstart,first,step;
			ierr = ISStrideGetInfo(isc,&first,&step);CHKERRQ(ierr);
			ierr = ISGetLocalSize(isc,&n);CHKERRQ(ierr);
            ierr = MatGetLocalSize(A,NULL,&nc);CHKERRQ(ierr);
            ierr = MatGetOwnershipRange(A,&cstart,NULL);CHKERRQ(ierr);
			if  (   (nc == n)
					 && (cstart == first)
					 && (step == 1) ) {
				//isFullCol = PETSC_TRUE;
				PetscPrintf(PetscObjectComm((PetscObject)A),"Detected full column space\n");
				
				is_list[0] = ctx->isU;
				is_list[1] = ctx->isV;
				is_list[2] = ctx->isW;
				
				f1 = PETSC_FALSE;
				for (d=0; d<3; d++) {
					ierr = ISEqual(isr,is_list[d],&f1);CHKERRQ(ierr);
					if (f1==PETSC_TRUE){ break; }
				}
				
				SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Requested row which matrix-free operator doesn't define - need to determine which uu,vv,ww component requested");
				
				if (cll==MAT_INITIAL_MATRIX) {
					PetscPrintf(PetscObjectComm((PetscObject)A),"  defining matrix free operator\n");
					PetscPrintf(PetscObjectComm((PetscObject)A),"ierr = MatCreateAiStokesMF(A,ii,B);CHKERRQ(ierr); TODO\n");
				} else {
					// to nothing
				}
				
				PetscFunctionReturn(0);
			}
		}
		
		
	}
	
	
	SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Something bad happened... I don't know how to retrieve the sub matrix you requested");
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A"
PetscErrorCode MatMult_MFStokes_A(Mat A,Vec X,Vec Y)
{
	MatStokesMF       ctx;
  PetscErrorCode    ierr;
  DM                stokes_pack,dau,dap;
  Vec               XUloc,XPloc,YUloc,YPloc;
	Vec               Xu,Xp,Yu,Yp;
  PetscScalar       *LA_XUloc,*LA_XPloc;
  PetscScalar       *LA_YUloc,*LA_YPloc;
	
  PetscFunctionBegin;
  
	ierr = PetscLogEventBegin(MAT_MultMFA,A,X,Y,0);CHKERRQ(ierr);
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	stokes_pack = ctx->stokes_pack;
	
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	
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
#ifdef TATIN_HAVE_AVX
  ierr = MFStokesWrapper_A_AVX(ctx->volQ,dau,LA_XUloc,dap,LA_XPloc,LA_YUloc,LA_YPloc);CHKERRQ(ierr);
#else
	ierr = MFStokesWrapper_A(ctx->volQ,dau,LA_XUloc,dap,LA_XPloc,LA_YUloc,LA_YPloc);CHKERRQ(ierr);
#endif
  
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
	
//	{
//		PetscPrintf(PETSC_COMM_WORLD,"%s: DUMMY MAT MULT FOR PLUMBING TESTING \n", __FUNCT__);
//		ierr = VecCopy(X,Y);CHKERRQ(ierr);
//	}
	
	ierr = PetscLogEventEnd(MAT_MultMFA,A,X,Y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A_QuasiNewtonX"
PetscErrorCode MatMult_MFStokes_A_QuasiNewtonX(Mat A,Vec X,Vec Y)
{
	MatStokesMF       ctx;
  PetscErrorCode    ierr;
  DM                stokes_pack,dau,dap,dax;
  Vec               XUloc,XPloc,YUloc,YPloc;
	Vec               Xu,Xp,Yu,Yp;
  PetscScalar       *LA_XUloc,*LA_XPloc;
  PetscScalar       *LA_YUloc,*LA_YPloc;
	Vec               Xloc;
	PetscScalar       *LA_Xloc;
  PetscFunctionBegin;
  
	ierr = PetscLogEventBegin(MAT_MultMFA_QuasiNewtonX,A,X,Y,0);CHKERRQ(ierr);
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	stokes_pack = ctx->stokes_pack;
	
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dau,&dax);CHKERRQ(ierr);
	
  ierr = DMCompositeGetLocalVectors(stokes_pack,&XUloc,&XPloc);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&YUloc,&YPloc);CHKERRQ(ierr);
	/* fetch shifted coordinate vector from Mat A */
	Xloc = NULL;
	ierr = PetscObjectQuery((PetscObject)A,"MatA_QuasiNewtonX",(PetscObject*)&Xloc);CHKERRQ(ierr);
	if (!Xloc) { SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Require a shifted coordinate vector to have been set via PetscObjectCompose()"); }
	
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
	ierr = VecGetArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YUloc);CHKERRQ(ierr);
	ierr = VecZeroEntries(YPloc);CHKERRQ(ierr);
	ierr = VecGetArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	ierr = VecGetArray(YPloc,&LA_YPloc);CHKERRQ(ierr);
	
	/* momentum + continuity */
	ierr = MFStokesWrapper_A_UPX(ctx->volQ,dau,LA_XUloc,dap,LA_XPloc, dax,LA_Xloc, LA_YUloc,LA_YPloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
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
	
	//	{
	//		PetscPrintf(PETSC_COMM_WORLD,"%s: DUMMY MAT MULT FOR PLUMBING TESTING \n", __FUNCT__);
	//		ierr = VecCopy(X,Y);CHKERRQ(ierr);
	//	}
	
	ierr = PetscLogEventEnd(MAT_MultMFA_QuasiNewtonX,A,X,Y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetDiagonal_MFStokes_A11"
PetscErrorCode MatGetDiagonal_MFStokes_A11(Mat A,Vec X)
{
	MatA11MF          ctx;
  PetscErrorCode    ierr;
  DM                dau;
  Vec               XUloc;
  PetscScalar       *LA_XUloc;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	
  ierr = DMGetLocalVector(dau,&XUloc);CHKERRQ(ierr);

	/* Zero input X */
	ierr = VecZeroEntries(X);CHKERRQ(ierr);
	ierr = VecZeroEntries(XUloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* A11 from momentum */
	ierr = MFStokesWrapper_diagA11(ctx->volQ,dau,LA_XUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = DMLocalToGlobalBegin(dau,XUloc,ADD_VALUES,X);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (dau,XUloc,ADD_VALUES,X);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(dau,&XUloc);CHKERRQ(ierr);
	
	/* modify X for the boundary conditions, x_k = scale_k(x_k) */

	/* FOR THE MOMENT THE DIAGONAL IS ALWAYS 1 x_k = scale_k(1.0) */
	/* Clobbering entries in global vector corresponding to dirichlet boundary conditions */
	ierr = BCListInsertValueIntoDirichletSlot(ctx->u_bclist,1.0,X);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetDiagonal_MFStokes_A11LowOrder"
PetscErrorCode MatGetDiagonal_MFStokes_A11LowOrder(Mat A,Vec X)
{
	MatA11MF          ctx;
  PetscErrorCode    ierr;
  DM                dau;
  Vec               XUloc;
  PetscScalar       *LA_XUloc;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	
  ierr = DMGetLocalVector(dau,&XUloc);CHKERRQ(ierr);
	
	/* Zero input X */
	ierr = VecZeroEntries(X);CHKERRQ(ierr);
	ierr = VecZeroEntries(XUloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* A11 from momentum */
	ierr = MFStokesWrapper_diagA11LowOrder(ctx->volQ,dau,LA_XUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = DMLocalToGlobalBegin(dau,XUloc,ADD_VALUES,X);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (dau,XUloc,ADD_VALUES,X);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(dau,&XUloc);CHKERRQ(ierr);
	
	/* modify X for the boundary conditions, x_k = scale_k(x_k) */
	
	/* FOR THE MOMENT THE DIAGONAL IS ALWAYS 1 x_k = scale_k(1.0) */
	/* Clobbering entries in global vector corresponding to dirichlet boundary conditions */
	ierr = BCListInsertValueIntoDirichletSlot(ctx->u_bclist,1.0,X);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A11"
PetscErrorCode MatMult_MFStokes_A11(Mat A,Vec X,Vec Y)
{
	MatA11MF          ctx;
  PetscErrorCode    ierr;
  DM                dau;
  Vec               XUloc,YUloc;
  PetscScalar       *LA_XUloc;
  PetscScalar       *LA_YUloc;
//	PetscBool         use_low_order_geometry = PETSC_FALSE;
	PetscObjectState  state;
	
  PetscFunctionBegin;
  
	ierr = PetscLogEventBegin(MAT_MultMFA11,A,X,Y,0);CHKERRQ(ierr);
//	ierr = PetscOptionsGetBool(NULL,NULL,"-use_low_order_geometry",&use_low_order_geometry,NULL);CHKERRQ(ierr);
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  
  if (!ctx->is_setup) SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_USER,"MF operator not setup");
  
	ierr = PetscObjectStateGet((PetscObject)A,&state);CHKERRQ(ierr);
	ctx->state = state;
  
	dau = ctx->daUVW;
	
	ierr = PetscLogEventBegin(MAT_MultMFA11_gv,A,X,Y,0);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dau,&XUloc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dau,&YUloc);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(MAT_MultMFA11_gv,A,X,Y,0);CHKERRQ(ierr);
	
	/* get the local (ghosted) entries for each physics */
	ierr = PetscLogEventBegin(MAT_MultMFA11_gtl,A,X,Y,0);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(MAT_MultMFA11_gtl,A,X,Y,0);CHKERRQ(ierr);
	
	/* Zero entries in local vectors corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out columns when the mat-mult is performed */
	ierr = BCListInsertLocalZero(ctx->u_bclist,XUloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YUloc);CHKERRQ(ierr);
	ierr = VecGetArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	
	/* momentum */
	ierr = (*ctx->SpMVOp_MatMult)(ctx,ctx->volQ,dau,LA_XUloc,LA_YUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = PetscLogEventBegin(MAT_MultMFA11_yz,A,X,Y,0);CHKERRQ(ierr);
	ierr = VecZeroEntries(Y);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(MAT_MultMFA11_yz,A,X,Y,0);CHKERRQ(ierr);
	ierr = PetscLogEventBegin(MAT_MultMFA11_ltg,A,X,Y,0);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(dau,YUloc,ADD_VALUES,Y);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (dau,YUloc,ADD_VALUES,Y);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(dau,&YUloc);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dau,&XUloc);CHKERRQ(ierr);
	ierr = PetscLogEventEnd(MAT_MultMFA11_ltg,A,X,Y,0);CHKERRQ(ierr);
	
	/* modify Y for the boundary conditions, y_k = scale_k(x_k) */
	/* Clobbering entries in global vector corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out rows when the mat-mult is performed */
	ierr = BCListInsertDirichlet_MatMult(ctx->u_bclist,X,Y);CHKERRQ(ierr);
	
	ierr = PetscLogEventEnd(MAT_MultMFA11,A,X,Y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A11LowOrder"
PetscErrorCode MatMult_MFStokes_A11LowOrder(Mat A,Vec X,Vec Y)
{
	MatA11MF          ctx;
  PetscErrorCode    ierr;
  DM                dau;
  Vec               XUloc,YUloc;
  PetscScalar       *LA_XUloc;
  PetscScalar       *LA_YUloc;
	PetscInt          low_order_geometry_type = 1; /* 0 - none; 1 - 1gp on Jac; 2 - 2x2x2 quad; 3 - 1x1x1 quad */
	
  PetscFunctionBegin;
  
	//ierr = PetscOptionsGetInt(NULL,NULL,"-low_order_geometry_type",&low_order_geometry_type,NULL);CHKERRQ(ierr);
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	
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
	
	/* momentum */
	switch (low_order_geometry_type) {
		case 0:
			ierr = MFStokesWrapper_A11(ctx,ctx->volQ,dau,LA_XUloc,LA_YUloc);CHKERRQ(ierr);
			break;

		case 1:
			ierr = MFStokesWrapper_A11PC(ctx->volQ,dau,LA_XUloc,LA_YUloc);CHKERRQ(ierr);
			break;

		case 2:
			ierr = MFStokesWrapper_A11PC_2x2x2(ctx->volQ,dau,LA_XUloc,LA_YUloc);CHKERRQ(ierr);
			break;

		case 3:
			// totally useless - u solve requires loads of iterations
			ierr = MFStokesWrapper_A11PC_1x1x1(ctx->volQ,dau,LA_XUloc,LA_YUloc);CHKERRQ(ierr);
			break;
	}
	
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
  Vec               XPloc,YUloc;
  PetscScalar       *LA_XPloc;
  PetscScalar       *LA_YUloc;
	
  PetscFunctionBegin;
  
	ierr = PetscLogEventBegin(MAT_MultMFA12,A,X,Y,0);CHKERRQ(ierr);

	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	dap = ctx->dap;
	
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
	
	/* grad */
#ifdef TATIN_HAVE_AVX
  ierr = MFStokesWrapper_A12_AVX(ctx->volQ,dau,dap,LA_XPloc,LA_YUloc);CHKERRQ(ierr);
#else
  ierr = MFStokesWrapper_A12(ctx->volQ,dau,dap,LA_XPloc,LA_YUloc);CHKERRQ(ierr);
#endif
  
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
	
//	{
//		PetscPrintf(PETSC_COMM_WORLD,"%s: DUMMY MAT MULT FOR PLUMBING TESTING \n", __FUNCT__);
//		ierr = VecSet(Y,1.0);CHKERRQ(ierr);
//	}
	ierr = PetscLogEventEnd(MAT_MultMFA12,A,X,Y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A12_QuasiNewtonX"
PetscErrorCode MatMult_MFStokes_A12_QuasiNewtonX(Mat A,Vec X,Vec Y)
{
	MatStokesMF       ctx;
  PetscErrorCode    ierr;
  DM                dau,dap,dax;
  Vec               XPloc,YUloc,Xloc;
  PetscScalar       *LA_XPloc;
  PetscScalar       *LA_YUloc;
  PetscScalar       *LA_Xloc;
	
  PetscFunctionBegin;

	ierr = PetscLogEventBegin(MAT_MultMFA12_QuasiNewtonX,A,X,Y,0);CHKERRQ(ierr);
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	dap = ctx->dap;
	ierr = DMGetCoordinateDM(dau,&dax);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dap,&XPloc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dau,&YUloc);CHKERRQ(ierr);
	/* fetch shifted coordinate vector from Mat A */
	Xloc = NULL;
	ierr = PetscObjectQuery((PetscObject)A,"MatA_QuasiNewtonX",(PetscObject*)&Xloc);CHKERRQ(ierr);
	if (!Xloc) { SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Require a shifted coordinate vector to have been set via PetscObjectCompose()"); }
	
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
	ierr = VecGetArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YUloc);CHKERRQ(ierr);
	ierr = VecGetArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	
	/* grad */
	ierr = MFStokesWrapper_A12_UPX(ctx->volQ,dau,dap,LA_XPloc,dax,LA_Xloc,LA_YUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XPloc,&LA_XPloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
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
	
	ierr = PetscLogEventEnd(MAT_MultMFA12_QuasiNewtonX,A,X,Y,0);CHKERRQ(ierr);
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
  Vec               XUloc,YPloc;
  PetscScalar       *LA_XUloc;
  PetscScalar       *LA_YPloc;
	
  PetscFunctionBegin;

	ierr = PetscLogEventBegin(MAT_MultMFA21,A,X,Y,0);CHKERRQ(ierr);
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	dap = ctx->dap;
	
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
	
	/* div */
	ierr = MFStokesWrapper_A21(ctx->volQ,dau,dap,LA_XUloc,LA_YPloc);CHKERRQ(ierr);
	
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
	
//	{
//		PetscPrintf(PETSC_COMM_WORLD,"%s: DUMMY MAT MULT FOR PLUMBING TESTING \n", __FUNCT__);
//		ierr = VecSet(Y,1.0);CHKERRQ(ierr);
//	}
	ierr = PetscLogEventEnd(MAT_MultMFA21,A,X,Y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A21_QuasiNewtonX"
PetscErrorCode MatMult_MFStokes_A21_QuasiNewtonX(Mat A,Vec X,Vec Y)
{
	MatStokesMF       ctx;
  PetscErrorCode    ierr;
  DM                dau,dap,dax;
  Vec               XUloc,YPloc,Xloc;
  PetscScalar       *LA_XUloc;
  PetscScalar       *LA_YPloc;
  PetscScalar       *LA_Xloc;
	
  PetscFunctionBegin;
  
	ierr = PetscLogEventBegin(MAT_MultMFA21_QuasiNewtonX,A,X,Y,0);CHKERRQ(ierr);

	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	dap = ctx->dap;
	ierr = DMGetCoordinateDM(dau,&dax);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&XUloc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dap,&YPloc);CHKERRQ(ierr);
	/* fetch shifted coordinate vector from Mat A */
	Xloc = NULL;
	ierr = PetscObjectQuery((PetscObject)A,"MatA_QuasiNewtonX",(PetscObject*)&Xloc);CHKERRQ(ierr);
	if (!Xloc) { SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Require a shifted coordinate vector to have been set via PetscObjectCompose()"); }
	
	/* get the local (ghosted) entries for each physics */
	ierr = DMGlobalToLocalBegin(dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	
	/* Zero entries in local vectors corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out columns when the mat-mult is performed */
	ierr = BCListInsertLocalZero(ctx->u_bclist,XUloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	ierr = VecGetArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YPloc);CHKERRQ(ierr);
	ierr = VecGetArray(YPloc,&LA_YPloc);CHKERRQ(ierr);
	
	/* div */
	ierr = MFStokesWrapper_A21_UPX(ctx->volQ,dau,LA_XUloc,dap,dax,LA_Xloc,LA_YPloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(YPloc,&LA_YPloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
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
	
	ierr = PetscLogEventEnd(MAT_MultMFA21_QuasiNewtonX,A,X,Y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_MFStokes_A11_QuasiNewtonX"
PetscErrorCode MatMult_MFStokes_A11_QuasiNewtonX(Mat A,Vec X,Vec Y)
{
	MatA11MF          ctx;
  PetscErrorCode    ierr;
  DM                dau,dax;
  Vec               XUloc,YUloc,Xloc;
  PetscScalar       *LA_XUloc,*LA_Xloc;
  PetscScalar       *LA_YUloc;
	
  PetscFunctionBegin;

	ierr = PetscLogEventBegin(MAT_MultMFA11_QuasiNewtonX,A,X,Y,0);CHKERRQ(ierr);
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	ierr = DMGetCoordinateDM(dau,&dax);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&XUloc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dau,&YUloc);CHKERRQ(ierr);
	/* fetch shifted coordinate vector from Mat A */
	Xloc = NULL;
	ierr = PetscObjectQuery((PetscObject)A,"MatA11_QuasiNewtonX",(PetscObject*)&Xloc);CHKERRQ(ierr);
	if (!Xloc) { SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Require a shifted coordinate vector to have been set via PetscObjectCompose()"); }
	
	/* get the local (ghosted) entries for each physics */
	ierr = DMGlobalToLocalBegin(dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (dau,X,INSERT_VALUES,XUloc);CHKERRQ(ierr);
		
	/* Zero entries in local vectors corresponding to dirichlet boundary conditions */
	/* This has the affect of zeroing out columns when the mat-mult is performed */
	ierr = BCListInsertLocalZero(ctx->u_bclist,XUloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	ierr = VecGetArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
	/* compute Ax - b */
	ierr = VecZeroEntries(YUloc);CHKERRQ(ierr);
	ierr = VecGetArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	
	/* momentum */
	ierr = MFStokesWrapper_A11_UPX(ctx->volQ,dau,LA_XUloc,dax,LA_Xloc,LA_YUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(YUloc,&LA_YUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
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
	
	ierr = PetscLogEventEnd(MAT_MultMFA11_QuasiNewtonX,A,X,Y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetDiagonal_MFStokes_A11_QuasiNewtonX"
PetscErrorCode MatGetDiagonal_MFStokes_A11_QuasiNewtonX(Mat A,Vec X)
{
	MatA11MF          ctx;
  PetscErrorCode    ierr;
  DM                dau,dax;
  Vec               XUloc,Xloc;
  PetscScalar       *LA_XUloc,*LA_Xloc;
	
  PetscFunctionBegin;
  
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	dau = ctx->daUVW;
	ierr = DMGetCoordinateDM(dau,&dax);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&XUloc);CHKERRQ(ierr);
	
	/* fetch shifted coordinate vector from Mat A */
	Xloc = NULL;
	ierr = PetscObjectQuery((PetscObject)A,"MatA11_QuasiNewtonX",(PetscObject*)&Xloc);CHKERRQ(ierr);
	if (!Xloc) { SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Require a shifted coordinate vector to have been set via PetscObjectCompose()"); }
	
	/* Zero input X */
	ierr = VecZeroEntries(X);CHKERRQ(ierr);
	ierr = VecZeroEntries(XUloc);CHKERRQ(ierr);
	
	ierr = VecGetArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	ierr = VecGetArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
	/* A11 from momentum */
	ierr = MFStokesWrapper_diagA11_UPX(ctx->volQ,dau,dax,LA_Xloc,LA_XUloc);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(XUloc,&LA_XUloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Xloc,&LA_Xloc);CHKERRQ(ierr);
	
	/* do global fem summation */
	ierr = DMLocalToGlobalBegin(dau,XUloc,ADD_VALUES,X);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (dau,XUloc,ADD_VALUES,X);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(dau,&XUloc);CHKERRQ(ierr);
	
	/* modify X for the boundary conditions, x_k = scale_k(x_k) */
	
	/* FOR THE MOMENT THE DIAGONAL IS ALWAYS 1 x_k = scale_k(1.0) */
	/* Clobbering entries in global vector corresponding to dirichlet boundary conditions */
	ierr = BCListInsertValueIntoDirichletSlot(ctx->u_bclist,1.0,X);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A(MatStokesMF Stk,Mat *A)
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
	
	*A = B;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A_QuasiNewtonX"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A_QuasiNewtonX(MatStokesMF Stk,Mat *A)
{
	DM  dax;
	Vec Xloc;
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	Stk->refcnt++;
	ierr = MatCreateShell(PETSC_COMM_WORLD,Stk->mu+Stk->mp,Stk->mu+Stk->mp,Stk->Mu+Stk->Mp,Stk->Mu+Stk->Mp,(void*)Stk,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A_QuasiNewtonX);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_GET_SUBMATRIX,(void(*)(void))MatGetSubMatrix_MFStokes_A);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatStokesMF_QuasiNewtonX);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(Stk->daUVW,&dax);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(dax,&Xloc);CHKERRQ(ierr);
	ierr = PetscObjectCompose((PetscObject)B,"MatA_QuasiNewtonX",(PetscObject)Xloc);CHKERRQ(ierr);
	
	*A = B;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A11"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11(MatA11MF A11,Mat *A)
{
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

  if ((A11->refcnt > 1) && (A11->ctx)) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Not clear how to safely share MF-SpMV context");
  }

	ierr = MatCreateShell(PETSC_COMM_WORLD,A11->mu,A11->mu,A11->Mu,A11->Mu,(void*)A11,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A11);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_GET_SUBMATRIX,(void(*)(void))MatGetSubMatrix_MFStokes_A11);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatA11MF);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_MFStokes_A11);CHKERRQ(ierr);
	ierr = MatSetBlockSize(B,3);CHKERRQ(ierr);
	A11->refcnt++;	
	*A = B;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A11_QuasiNewtonX"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11_QuasiNewtonX(MatA11MF A11,Mat *A)
{
	DM  dax;
	Vec Xloc;
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
  if ((A11->refcnt > 1) && (A11->ctx)) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Not clear how to safely share MF-SpMV context");
  }
	ierr = MatCreateShell(PETSC_COMM_WORLD,A11->mu,A11->mu,A11->Mu,A11->Mu,(void*)A11,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A11_QuasiNewtonX);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_GET_SUBMATRIX,(void(*)(void))MatGetSubMatrix_MFStokes_A11);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatA11MF_QuasiNewtonX);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_MFStokes_A11_QuasiNewtonX);CHKERRQ(ierr);
	ierr = MatSetBlockSize(B,3);CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(A11->daUVW,&dax);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(dax,&Xloc);CHKERRQ(ierr);
	ierr = PetscObjectCompose((PetscObject)B,"MatA11_QuasiNewtonX",(PetscObject)Xloc);CHKERRQ(ierr);
	A11->refcnt++;	
	*A = B;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A11LowOrder"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A11LowOrder(MatA11MF A11,Mat *A)
{
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
  if ((A11->refcnt > 1) && (A11->ctx)) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Not clear how to safely share MF-SpMV context");
  }

	ierr = MatCreateShell(PETSC_COMM_WORLD,A11->mu,A11->mu,A11->Mu,A11->Mu,(void*)A11,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A11LowOrder);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
//	ierr = MatShellSetOperation(B,MATOP_GET_SUBMATRIX,(void(*)(void))MatGetSubMatrix_MFStokes_A11LowOrder);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatA11MF);CHKERRQ(ierr);
	// does the true diagonal work better??
	ierr = MatShellSetOperation(B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_MFStokes_A11);CHKERRQ(ierr);
	//ierr = MatShellSetOperation(B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_MFStokes_A11LowOrder);CHKERRQ(ierr);
	ierr = MatSetBlockSize(B,3);CHKERRQ(ierr);

	 A11->refcnt++;	
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
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A12_QuasiNewtonX"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A12_QuasiNewtonX(MatStokesMF Stk,Mat *A12)
{
	DM  dax;
	Vec Xloc;
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	Stk->refcnt++;
	ierr = MatCreateShell(PETSC_COMM_WORLD,Stk->mu,Stk->mp,Stk->Mu,Stk->Mp,(void*)Stk,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A12_QuasiNewtonX);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatStokesMF_QuasiNewtonX);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(Stk->daUVW,&dax);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(dax,&Xloc);CHKERRQ(ierr);
	ierr = PetscObjectCompose((PetscObject)B,"MatA_QuasiNewtonX",(PetscObject)Xloc);CHKERRQ(ierr);
	
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

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_A21_QuasiNewtonX"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_A21_QuasiNewtonX(MatStokesMF Stk,Mat *A21)
{
	DM  dax;
	Vec Xloc;
	Mat B;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	Stk->refcnt++;
	ierr = MatCreateShell(PETSC_COMM_WORLD,Stk->mp,Stk->mu,Stk->Mp,Stk->Mu,(void*)Stk,&B);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_MFStokes_A21_QuasiNewtonX);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_basic);CHKERRQ(ierr);
	ierr = MatShellSetOperation(B,MATOP_DESTROY,(void(*)(void))MatDestroy_MatStokesMF_QuasiNewtonX);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(Stk->daUVW,&dax);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(dax,&Xloc);CHKERRQ(ierr);
	ierr = PetscObjectCompose((PetscObject)B,"MatA_QuasiNewtonX",(PetscObject)Xloc);CHKERRQ(ierr);
	
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
	DM             pack;
	Mat            A;
	MatStokesMF    StkCtx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	pack = user->stokes_pack;

	/* is composite */
	same = PETSC_FALSE;
	ierr = PetscObjectTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
	/* Create submatrices */
	//comm = PetscObjectComm((PetscObject)pack);
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A(StkCtx,&A);CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(A,"stokes_Amf_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
	
	*B = A;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_MFOperator_QuasiNewtonX"
PetscErrorCode StokesQ2P1CreateMatrix_MFOperator_QuasiNewtonX(PhysCompStokes user,Mat *B)
{
	PetscBool      same;
    DM             pack;
	Mat            A;
	MatStokesMF    StkCtx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = MatStokesMFCreate(&StkCtx);CHKERRQ(ierr);
	ierr = MatStokesMFSetup(StkCtx,user);CHKERRQ(ierr);
	pack = user->stokes_pack;
	
	/* is composite */
	same = PETSC_FALSE;
	ierr = PetscObjectTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
	/* Create submatrices */
	//comm = ((PetscObject)pack)->comm;
	
	ierr = StokesQ2P1CreateMatrix_MFOperator_A_QuasiNewtonX(StkCtx,&A);CHKERRQ(ierr);
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
	PetscInt       i,j;
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
	ierr = PetscObjectTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
	/* Fetch the DA's */
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	/* Create submatrices */
	//comm = PetscObjectComm((PetscObject)pack);

	/* A11 */
	if (tA11==0) {
		ierr = DMSetMatType(dau,MATAIJ);CHKERRQ(ierr);
		ierr = DMCreateMatrix(dau,&Auu);CHKERRQ(ierr);
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
	}
	ierr = MatSetOptionsPrefix(Auu,"stokes_A_A11_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);
	
	/* A12 */
	if (tA12==0) {
		SETERRQ(PetscObjectComm((PetscObject)dau),PETSC_ERR_SUP,"Stokes->A12 doesn't support assembled operator");
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A12(StkCtx,&Aup);CHKERRQ(ierr);
	}
	ierr = MatSetOptionsPrefix(Aup,"stokes_A_A12_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Aup);CHKERRQ(ierr);

	/* A21 */
	if (tA21==0) {
		SETERRQ(PetscObjectComm((PetscObject)dau),PETSC_ERR_SUP,"Stokes->A21 doesn't support assembled operator");
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A21(StkCtx,&Apu);CHKERRQ(ierr); 
	}
	ierr = MatSetOptionsPrefix(Apu,"stokes_A_A21_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Apu);CHKERRQ(ierr);
	
	/* Create nest */
	ierr = DMCompositeGetGlobalISs(pack,&is);CHKERRQ(ierr);
	
	bA[0][0] = Auu; bA[0][1] = Aup;
	bA[1][0] = Apu; bA[1][1] = NULL;
	ierr = MatCreateNest(PetscObjectComm((PetscObject)dau),2,is,2,is,&bA[0][0],&A);CHKERRQ(ierr);
	ierr = MatSetOptionsPrefix(A,"stokes_A_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	*B = A;
	
	/* tidy up */
	ierr = MatStokesMFDestroy(&StkCtx);CHKERRQ(ierr);
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

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrixNest_PCOperator"
PetscErrorCode StokesQ2P1CreateMatrixNest_PCOperator(PhysCompStokes user,PetscInt tA11,PetscInt tA12,PetscInt tA21,Mat *B)
{
	PetscBool      same;
  DM             dau,dap,pack;
	Mat            A,Auu,Aup,Apu,Spp,bA[2][2];
	IS             *is;
	PetscInt       i,j;
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
	ierr = PetscObjectTypeCompare((PetscObject)pack,DMCOMPOSITE,&same);CHKERRQ(ierr);
	if (!same) PetscFunctionReturn(0);
	
  /* Fetch the DA's */
  ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	/* Create submatrices */
	//comm = PetscObjectComm((PetscObject)pack);
	
	/* A11 */
	if (tA11==0) {
//		ierr = DMCreateMatrix(dau,MATAIJ,&Auu);CHKERRQ(ierr);
//
		ierr = DMSetMatType(dau,MATSBAIJ);CHKERRQ(ierr);
		ierr = DMCreateMatrix(dau,&Auu);CHKERRQ(ierr);
		ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
// 
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
	}
	ierr = MatSetOptionsPrefix(Auu,"Buu_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);
	
	/* Schur complement */
//	ierr = DMCreateMatrix(dap,MATAIJ,&Spp);CHKERRQ(ierr);

	ierr = DMSetMatType(dap,MATSBAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(dap,&Spp);CHKERRQ(ierr);
	ierr = MatSetOption(Spp,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);

	ierr = MatSetOptionsPrefix(Spp,"S*_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Spp);CHKERRQ(ierr);

	/* A12 */
	if (tA12==0) {
		SETERRQ(PetscObjectComm((PetscObject)dau),PETSC_ERR_SUP,"Stokes->A12 doesn't support assembled operator");
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A12(StkCtx,&Aup);CHKERRQ(ierr);
	}
	ierr = MatSetOptionsPrefix(Aup,"Bup_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Aup);CHKERRQ(ierr);
	
	/* A21 */
	if (tA21==0) {
		SETERRQ(PetscObjectComm((PetscObject)dau),PETSC_ERR_SUP,"Stokes->A21 doesn't support assembled operator");
	} else {
		ierr = StokesQ2P1CreateMatrix_MFOperator_A21(StkCtx,&Apu);CHKERRQ(ierr);
	}
	ierr = MatSetOptionsPrefix(Apu,"Bpu_");CHKERRQ(ierr);
	ierr = MatSetFromOptions(Apu);CHKERRQ(ierr);
	
	/* Create nest */
	ierr = DMCompositeGetGlobalISs(pack,&is);CHKERRQ(ierr);
	
	bA[0][0] = Auu;        bA[0][1] = Aup;
	bA[1][0] = Apu;        bA[1][1] = Spp;
  ierr = MatCreateNest(PetscObjectComm((PetscObject)dau),2,is,2,is,&bA[0][0],&A);CHKERRQ(ierr);
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

#undef __FUNCT__  
#define __FUNCT__ "StokesA12Preallocation_basic"
PetscErrorCode StokesA12Preallocation_basic(Mat mat,DM dav,DM dap)
{
	PetscInt nnz,onnz;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* n_pressure_basis * neighbour_cells = 4 x 8 */
	nnz = 32;
	ierr = PetscOptionsGetInt(NULL,NULL,"-A12_preallocation_nnz",&nnz,&flg);CHKERRQ(ierr);

        /*
        I've tuned this parameter for speed, not memory usage
        BRUTUS:
        8 cpus  mx^3 = 60 onnz =  4  assembly > 5 mins
                          onnz =  8  assembly > 5 mins
                          onnz = 16  assembly = 4.5953e+01 sec
	                  onnz = 24  assembly = 4.1227e+00 sec

        64 cpus mx^3 = 120 onnz =  4 assembly = sec
                           onnz =  8 assembly = sec
                           onnz = 16 assembly = 8.5164e+01 sec
	                   onnz = 24 assembly = 5.3824e+00 sec
        */
	onnz = 16;
	ierr = PetscOptionsGetInt(NULL,NULL,"-A12_preallocation_onnz",&onnz,&flg);CHKERRQ(ierr);
	
	PetscPrintf(PetscObjectComm((PetscObject)mat),"StokesA12Preallocation_basic: using nnz = %D and onnz = %D \n", nnz,onnz );
	
	ierr = MatSeqAIJSetPreallocation(mat,nnz,NULL);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(mat,nnz,NULL,onnz,NULL);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesA21Preallocation_basic"
PetscErrorCode StokesA21Preallocation_basic(Mat mat,DM dav,DM dap)
{
	PetscInt nnz,onnz;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	/* Each pressure dof is connected to all vel dofs in a single cell, 27 * 3 */
	nnz = 27 * 3;
	ierr = PetscOptionsGetInt(NULL,NULL,"-A21_preallocation_nnz",&nnz,&flg);CHKERRQ(ierr);

        /*
        I've tuned this parameter for speed, not memory usage
        BRUTUS:
        8 cpus  mx^3 = 60 onnz = 4   assembly = 2.2707e+02 sec
        8 cpus  mx^3 = 60 onnz = 10  assembly > 5 mins
	                  onnz = 20  assembly > 5 mins
	                  onnz = 30  assembly = 1.9846e+01 sec
	                  onnz = 40  assembly = 3.0207e+01 sec

        64 cpus mx^3 = 120 onnz = 4 assembly = 3.6358e+02 sec

	                   onnz = 30 assembly = 2.7790e+01 sec 
	                   onnz = 40 assembly = 4.0755e+01 sec
	                   onnz = 50 assembly = 3.3143e+00 sec
	                   onnz = 60 assembly = 2.9067e+00 sec
        */
	onnz = 30;
	ierr = PetscOptionsGetInt(NULL,NULL,"-A21_preallocation_onnz",&onnz,&flg);CHKERRQ(ierr);
	
	PetscPrintf(PetscObjectComm((PetscObject)mat),"StokesA21Preallocation_basic: using nnz = %D and onnz = %D \n", nnz,onnz );
	
	ierr = MatSeqAIJSetPreallocation(mat,nnz,NULL);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(mat,nnz,NULL,onnz,NULL);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_A12"
PetscErrorCode StokesQ2P1CreateMatrix_A12(PhysCompStokes user,Mat *mat)
{
  DM             pack,dav,dap;
	Mat            Aup;
	MPI_Comm       comm;
	PetscInt       mu,Mu,mp,Mp;
	Vec            X,Xu,Xp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
  /* Fetch the DA's */
	pack = user->stokes_pack;
  ierr = DMCompositeGetEntries(pack,&dav,&dap);CHKERRQ(ierr);
	
	/* Fetch sizes */
	ierr = DMGetGlobalVector(pack,&X);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(pack,X,&Xu,&Xp);CHKERRQ(ierr);
	ierr = VecGetSize(Xu,&Mu);CHKERRQ(ierr);
	ierr = VecGetLocalSize(Xu,&mu);CHKERRQ(ierr);
	ierr = VecGetSize(Xp,&Mp);CHKERRQ(ierr);
	ierr = VecGetLocalSize(Xp,&mp);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(pack,X,&Xu,&Xp);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(pack,&X);CHKERRQ(ierr);
	
	/* Create matrix - this is not meant to be efficient non-zero allocation */
	comm = PetscObjectComm((PetscObject)pack);
  ierr = MatCreate(comm,&Aup);CHKERRQ(ierr);
	ierr = MatSetType(Aup,MATAIJ);CHKERRQ(ierr);
	ierr = MatSetSizes(Aup,mu,mp,Mu,Mp);CHKERRQ(ierr);
	
	/* Do some preallocation */
	ierr = StokesA12Preallocation_basic(Aup,dav,dap);CHKERRQ(ierr);
	
	ierr = MatSetOption(Aup,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE);CHKERRQ(ierr);
	
	ierr = MatSetFromOptions(Aup);CHKERRQ(ierr);
	
	*mat = Aup;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "StokesQ2P1CreateMatrix_A21"
PetscErrorCode StokesQ2P1CreateMatrix_A21(PhysCompStokes user,Mat *mat)
{
  DM             pack,dav,dap;
	Mat            Apu;
	MPI_Comm       comm;
	PetscInt       mu,Mu,mp,Mp;
	Vec            X,Xu,Xp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
  /* Fetch the DA's */
	pack = user->stokes_pack;
  ierr = DMCompositeGetEntries(pack,&dav,&dap);CHKERRQ(ierr);
	
	/* Fetch sizes */
	ierr = DMGetGlobalVector(pack,&X);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(pack,X,&Xu,&Xp);CHKERRQ(ierr);
	ierr = VecGetSize(Xu,&Mu);CHKERRQ(ierr);
	ierr = VecGetLocalSize(Xu,&mu);CHKERRQ(ierr);
	ierr = VecGetSize(Xp,&Mp);CHKERRQ(ierr);
	ierr = VecGetLocalSize(Xp,&mp);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(pack,X,&Xu,&Xp);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(pack,&X);CHKERRQ(ierr);
	
	/* Create matrix - this is not meant to be efficient non-zero allocation */
	comm = PetscObjectComm((PetscObject)pack);
  ierr = MatCreate(comm,&Apu);CHKERRQ(ierr);
	ierr = MatSetType(Apu,MATAIJ);CHKERRQ(ierr);
	ierr = MatSetSizes(Apu,mp,mu,Mp,Mu);CHKERRQ(ierr);
	
	/* Do some preallocation */
	ierr = StokesA21Preallocation_basic(Apu,dav,dap);CHKERRQ(ierr);
	
	ierr = MatSetOption(Apu,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE);CHKERRQ(ierr);
	
	ierr = MatSetFromOptions(Apu);CHKERRQ(ierr);
	
	*mat = Apu;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatCreate_StokesA11_asm"
PetscErrorCode MatCreate_StokesA11_asm(PhysCompStokes user,const char prefix[],Mat *mat)
{
  DM             pack,dav,dap;
	PetscBool      same1,same2,same3;
	Mat            Auu;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	pack = user->stokes_pack;
  ierr = DMCompositeGetEntries(pack,&dav,&dap);CHKERRQ(ierr);

	ierr = DMSetMatType(dav,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(dav,&Auu);CHKERRQ(ierr);
	//ierr = DMCreateMatrix(dav,MATSBAIJ,&Auu);CHKERRQ(ierr);

	if (prefix) {
		ierr = MatSetOptionsPrefix(Auu,prefix);CHKERRQ(ierr);
	}
	ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);

	ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSBAIJ,&same1);CHKERRQ(ierr);
	ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSEQSBAIJ,&same2);CHKERRQ(ierr);
	ierr = PetscObjectTypeCompare((PetscObject)Auu,MATMPISBAIJ,&same3);CHKERRQ(ierr);
	if (same1||same2||same3) {
		ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
	}
	
	*mat = Auu;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatShellGetMatStokesMF"
PetscErrorCode MatShellGetMatStokesMF(Mat A,MatStokesMF *mf)
{
	MatStokesMF     ctx;
	PetscErrorCode  ierr;
	
  PetscFunctionBegin;
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	*mf = ctx;
	
  PetscFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "MatShellGetMatA11MF"
PetscErrorCode MatShellGetMatA11MF(Mat A,MatA11MF *mf)
{
	MatA11MF     ctx;
	PetscErrorCode  ierr;
	
  PetscFunctionBegin;
	
	ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
	*mf = ctx;
	
  PetscFunctionReturn(0);
}




#include <private/kspimpl.h>                    /*I "petscksp.h" I*/

/*  
 Private data structure for Chebychev Iteration 
*/

typedef struct {
  PetscReal emin,emax;   /* store user provided estimates of extreme eigenvalues */
  KSP kspest;            /* KSP used to estimate eigenvalues */
  PC  pcnone;            /* Dummy PC to drop in so PCSetFromOptions doesn't get called extra times */
  PetscReal tform[4];    /* transform from Krylov estimates to Chebychev bounds */
  PetscBool estimate_current;

	PetscBool    use_red_norm;
  VecScatter   scatterin; /* scatter used to move all values to each processor group (subcommunicator) */
  PetscSubcomm psubcomm;          
  PetscInt     nsubcomm;  /* num of data structure PetscSubcomm */
	Vec          xdup,xsub;
	
} KSP_ChebychevRN;


#undef __FUNCT__  
#define __FUNCT__ "KSPVecNorm_ChebychevRN"
static PetscErrorCode KSPVecNorm_ChebychevRN(KSP ksp,Vec x,NormType type,PetscReal *val)
{
  KSP_ChebychevRN  *cheb = (KSP_ChebychevRN*)ksp->data;
  PetscErrorCode   ierr;
  PetscScalar      *array;
	
  PetscFunctionBegin;
  /* scatter x to xdup */
  ierr = VecScatterBegin(cheb->scatterin,x,cheb->xdup,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd  (cheb->scatterin,x,cheb->xdup,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  
  /* place xdup's local array into xsub */
  ierr = VecGetArray(cheb->xdup,&array);CHKERRQ(ierr);
  ierr = VecPlaceArray(cheb->xsub,(const PetscScalar*)array);CHKERRQ(ierr);
	
  /* compute VecNorm on each processor subset */
  ierr = VecNorm(cheb->xsub,type,val);CHKERRQ(ierr);

  ierr = VecResetArray(cheb->xsub);CHKERRQ(ierr);
  ierr = VecRestoreArray(cheb->xdup,&array);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSPSetUp_ChebychevRN"
PetscErrorCode KSPSetUp_ChebychevRN(KSP ksp)
{
  KSP_ChebychevRN  *cheb = (KSP_ChebychevRN*)ksp->data;
  MPI_Comm         comm = ((PetscObject)ksp)->comm,subcomm;
  PetscErrorCode   ierr;
	
  PetscFunctionBegin;
  ierr = KSPDefaultGetWork(ksp,3);CHKERRQ(ierr);
  cheb->estimate_current = PETSC_FALSE;
	
	if (cheb->use_red_norm) {
    if (!cheb->psubcomm) {
      ierr = PetscSubcommCreate(comm,&cheb->psubcomm);CHKERRQ(ierr);
      ierr = PetscSubcommSetNumber(cheb->psubcomm,cheb->nsubcomm);CHKERRQ(ierr);
      ierr = PetscSubcommSetType(cheb->psubcomm,PETSC_SUBCOMM_INTERLACED);CHKERRQ(ierr);
      ierr = PetscLogObjectMemory(ksp,sizeof(PetscSubcomm));CHKERRQ(ierr);
			
      subcomm = cheb->psubcomm->comm;
    } else {
			subcomm = cheb->psubcomm->comm;
    }
	}
	
	if (cheb->use_red_norm) {
    if (!cheb->scatterin) {
			Mat             mat;
			Vec             vec;
			PetscInt        mstart,mend,mlocal,m,mlocal_sub,rstart_sub,rend_sub,mloc_sub;
			PetscMPIInt     size;
			PetscMPIInt     subsize,subrank;
			const PetscInt  *range;
			
			ierr = KSPGetOperators(ksp,&mat,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			ierr = MatGetVecs(mat,PETSC_NULL,&vec);CHKERRQ(ierr);

			/* create working vectors xsub */
			ierr = VecGetSize(vec,&m);CHKERRQ(ierr);
			ierr = VecGetLocalSize(vec,&mlocal);CHKERRQ(ierr);  
			ierr = VecGetOwnershipRange(vec,&mstart,&mend);CHKERRQ(ierr);
			
			/* get local size of xsub/ysub */    
			ierr = MPI_Comm_size(subcomm,&subsize);CHKERRQ(ierr);
			ierr = MPI_Comm_rank(subcomm,&subrank);CHKERRQ(ierr);
			ierr = VecGetOwnershipRanges(vec,&range);CHKERRQ(ierr);
			rstart_sub = range[cheb->psubcomm->n*subrank]; /* rstart in xsub */    
			if (subrank+1 < subsize){
				rend_sub = range[cheb->psubcomm->n*(subrank+1)];
			} else {
				rend_sub = m; 
			}
			mloc_sub = rend_sub - rstart_sub;
			printf("mloc_sub = %d \n",mloc_sub);
			printf("nsubcomm %d \n",cheb->nsubcomm);
			MPI_Barrier(PETSC_COMM_WORLD);

			/* create xsub with empty local arrays, because xdup's arrays will be placed into it */
			ierr = VecCreateMPIWithArray(subcomm,mloc_sub,PETSC_DECIDE,PETSC_NULL,&cheb->xsub);CHKERRQ(ierr);
			/* create xdup. Note: we use communicator dupcomm, not ((PetscObject)pc)->comm! */      
			ierr = VecCreateMPI(cheb->psubcomm->dupparent,mloc_sub,PETSC_DECIDE,&cheb->xdup);CHKERRQ(ierr);
			
			/* create vec scatters */
			{
				IS       is1,is2;
				PetscInt *idx1,*idx2,i,j,k;
				
				ierr = PetscMalloc2(cheb->psubcomm->n*mlocal,PetscInt,&idx1,cheb->psubcomm->n*mlocal,PetscInt,&idx2);CHKERRQ(ierr);
				j = 0;
				for (k=0; k<cheb->psubcomm->n; k++){
					for (i=mstart; i<mend; i++){
						idx1[j]   = i;
						idx2[j++] = i + m*k;
					}
				}
				ierr = ISCreateGeneral(comm,cheb->psubcomm->n*mlocal,idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
				ierr = ISCreateGeneral(comm,cheb->psubcomm->n*mlocal,idx2,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);      
				ierr = VecScatterCreate(vec,is1,cheb->xdup,is2,&cheb->scatterin);CHKERRQ(ierr);

				ierr = ISDestroy(&is1);CHKERRQ(ierr);
				ierr = ISDestroy(&is2);CHKERRQ(ierr);
				ierr = PetscFree2(idx1,idx2);CHKERRQ(ierr);
			}
			ierr = VecDestroy(&vec);CHKERRQ(ierr);
		}
		
	}
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "KSPChebychevSetEigenvalues_ChebychevRN"
PetscErrorCode  KSPChebychevSetEigenvalues_ChebychevRN(KSP ksp,PetscReal emax,PetscReal emin)
{
  KSP_ChebychevRN *chebychevP = (KSP_ChebychevRN*)ksp->data;
	
  PetscFunctionBegin;
  if (emax <= emin) SETERRQ2(((PetscObject)ksp)->comm,PETSC_ERR_ARG_INCOMP,"Maximum eigenvalue must be larger than minimum: max %g min %G",emax,emin);
  if (emax*emin <= 0.0) SETERRQ2(((PetscObject)ksp)->comm,PETSC_ERR_ARG_INCOMP,"Both eigenvalues must be of the same sign: max %G min %G",emax,emin);
  chebychevP->emax = emax;
  chebychevP->emin = emin;
  PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "KSPChebychevSetEstimateEigenvalues_ChebychevRN"
PetscErrorCode  KSPChebychevSetEstimateEigenvalues_ChebychevRN(KSP ksp,PetscReal a,PetscReal b,PetscReal c,PetscReal d)
{
  KSP_ChebychevRN  *cheb = (KSP_ChebychevRN*)ksp->data;
  PetscErrorCode   ierr;
	
  PetscFunctionBegin;
  if (a != 0.0 || b != 0.0 || c != 0.0 || d != 0.0) {
    if (!cheb->kspest) {
      PetscBool nonzero;
			
      ierr = KSPCreate(((PetscObject)ksp)->comm,&cheb->kspest);CHKERRQ(ierr);
			
      ierr = KSPGetPC(cheb->kspest,&cheb->pcnone);CHKERRQ(ierr);
      ierr = PetscObjectReference((PetscObject)cheb->pcnone);CHKERRQ(ierr);
      ierr = PCSetType(cheb->pcnone,PCNONE);CHKERRQ(ierr);
      ierr = KSPSetPC(cheb->kspest,ksp->pc);CHKERRQ(ierr);
			
      ierr = KSPGetInitialGuessNonzero(ksp,&nonzero);CHKERRQ(ierr);
      ierr = KSPSetInitialGuessNonzero(cheb->kspest,nonzero);CHKERRQ(ierr);
      ierr = KSPSetComputeSingularValues(cheb->kspest,PETSC_TRUE);CHKERRQ(ierr);
			
      /* Estimate with a fixed number of iterations */
      ierr = KSPSetConvergenceTest(cheb->kspest,KSPSkipConverged,0,0);CHKERRQ(ierr);
      ierr = KSPSetNormType(cheb->kspest,KSP_NORM_NONE);CHKERRQ(ierr);
      ierr = KSPSetTolerances(cheb->kspest,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,5);CHKERRQ(ierr);
			
      if (a >= 0) cheb->tform[0] = a;
      if (b >= 0) cheb->tform[1] = b;
      if (c >= 0) cheb->tform[2] = c;
      if (d >= 0) cheb->tform[3] = d;
    }
  } else {
    ierr = KSPDestroy(&cheb->kspest);CHKERRQ(ierr);
    ierr = PCDestroy(&cheb->pcnone);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "KSPChebychevRNSetEigenvalues"
/*@
 KSPChebychevSetEigenvalues - Sets estimates for the extreme eigenvalues
 of the preconditioned problem.
 
 Logically Collective on KSP
 
 Input Parameters:
 +  ksp - the Krylov space context
 -  emax, emin - the eigenvalue estimates
 
 Options Database:
 .  -ksp_chebychev_eigenvalues emin,emax
 
 Note: If you run with the Krylov method of KSPCG with the option -ksp_monitor_singular_value it will 
 for that given matrix and preconditioner an estimate of the extreme eigenvalues.
 
 Level: intermediate
 
 .keywords: KSP, Chebyshev, set, eigenvalues
 @*/
PetscErrorCode  KSPChebychevRNSetEigenvalues(KSP ksp,PetscReal emax,PetscReal emin)
{
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  PetscValidLogicalCollectiveReal(ksp,emax,2);
  PetscValidLogicalCollectiveReal(ksp,emin,3);
  ierr = PetscTryMethod(ksp,"KSPChebychevRNSetEigenvalues_C",(KSP,PetscReal,PetscReal),(ksp,emax,emin));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSPChebychevRNSetEstimateEigenvalues"
/*@
 KSPChebychevSetEstimateEigenvalues - Automatically estimate the eigenvalues to use for Chebychev
 
 Logically Collective on KSP
 
 Input Parameters:
 +  ksp - the Krylov space context
 .  a - multiple of min eigenvalue estimate to use for min Chebychev bound (or PETSC_DECIDE)
 .  b - multiple of max eigenvalue estimate to use for min Chebychev bound (or PETSC_DECIDE)
 .  c - multiple of min eigenvalue estimate to use for max Chebychev bound (or PETSC_DECIDE)
 -  d - multiple of max eigenvalue estimate to use for max Chebychev bound (or PETSC_DECIDE)
 
 Options Database:
 .  -ksp_chebychev_estimate_eigenvalues a,b,c,d
 
 Notes:
 The Chebychev bounds are estimated using
 .vb
 minbound = a*minest + b*maxest
 maxbound = c*minest + d*maxest
 .ve
 The default configuration targets the upper part of the spectrum for use as a multigrid smoother, so only the maximum eigenvalue estimate is used.
 The minimum eigenvalue estimate obtained by Krylov iteration is typically not accurate until the method has converged.
 
 If 0.0 is passed for all transform arguments (a,b,c,d), eigenvalue estimation is disabled.
 
 Level: intermediate
 
 .keywords: KSP, Chebyshev, set, eigenvalues
 @*/
PetscErrorCode KSPChebychevRNSetEstimateEigenvalues(KSP ksp,PetscReal a,PetscReal b,PetscReal c,PetscReal d)
{
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,1);
  PetscValidLogicalCollectiveReal(ksp,a,2);
  PetscValidLogicalCollectiveReal(ksp,b,3);
  PetscValidLogicalCollectiveReal(ksp,c,4);
  PetscValidLogicalCollectiveReal(ksp,d,5);
  ierr = PetscTryMethod(ksp,"KSPChebychevRNSetEstimateEigenvalues_C",(KSP,PetscReal,PetscReal,PetscReal,PetscReal),(ksp,a,b,c,d));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSPSetFromOptions_ChebychevRN"
PetscErrorCode KSPSetFromOptions_ChebychevRN(KSP ksp)
{
  KSP_ChebychevRN  *cheb = (KSP_ChebychevRN*)ksp->data;
  PetscErrorCode   ierr;
  PetscInt         two = 2,four = 4;
  PetscReal        tform[4] = {PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE};
  PetscBool        flg;
	PetscMPIInt      size;
	
  PetscFunctionBegin;
  ierr = PetscOptionsHead("KSP ChebychevRN Options");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-ksp_chebychev_global_reduction_redundant_number","Number of redundant vecs used in VecNorm/VecDot","KSPChebychevGlobalReductionRedundantSetNumber",cheb->nsubcomm,&cheb->nsubcomm,0);CHKERRQ(ierr);
	ierr = PetscOptionsRealArray("-ksp_chebychev_eigenvalues","extreme eigenvalues","KSPChebychevRNSetEigenvalues",&cheb->emin,&two,0);CHKERRQ(ierr);
  ierr = PetscOptionsRealArray("-ksp_chebychev_estimate_eigenvalues","estimate eigenvalues using a Krylov method, then use this transform for Chebychev eigenvalue bounds","KSPChebychevRNSetEstimateEigenvalues",tform,&four,&flg);CHKERRQ(ierr);
  if (flg) {ierr = KSPChebychevRNSetEstimateEigenvalues(ksp,tform[0],tform[1],tform[2],tform[3]);CHKERRQ(ierr);}
  if (cheb->kspest) {
    /* Mask the PC so that PCSetFromOptions does not do anything */
    ierr = KSPSetPC(cheb->kspest,cheb->pcnone);CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(cheb->kspest,((PetscObject)ksp)->prefix);CHKERRQ(ierr);
    ierr = KSPAppendOptionsPrefix(cheb->kspest,"est_");CHKERRQ(ierr);
    if (!((PetscObject)cheb->kspest)->type_name) {
      ierr = KSPSetType(cheb->kspest,KSPGMRES);CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(cheb->kspest);CHKERRQ(ierr);
    ierr = KSPSetPC(cheb->kspest,ksp->pc);CHKERRQ(ierr);
  }

  ierr = MPI_Comm_size(((PetscObject)ksp)->comm,&size);CHKERRQ(ierr);
	if (cheb->nsubcomm != size) {
		cheb->use_red_norm = PETSC_TRUE;
	}
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSolve_ChebychevRN"
PetscErrorCode KSPSolve_ChebychevRN(KSP ksp)
{
  KSP_ChebychevRN  *cheb = (KSP_ChebychevRN*)ksp->data;
  PetscErrorCode   ierr;
  PetscInt         k,kp1,km1,maxit,ktmp,i;
  PetscScalar      alpha,omegaprod,mu,omega,Gamma,c[3],scale;
  PetscReal        rnorm = 0.0;
  Vec              x,b,p[3],r;
  Mat              Amat,Pmat;
  MatStructure     pflag;
  PetscBool        diagonalscale;
	
  PetscFunctionBegin;
  ierr    = PCGetDiagonalScale(ksp->pc,&diagonalscale);CHKERRQ(ierr);
  if (diagonalscale) SETERRQ1(((PetscObject)ksp)->comm,PETSC_ERR_SUP,"Krylov method %s does not support diagonal scaling",((PetscObject)ksp)->type_name);
	
  if (cheb->kspest && !cheb->estimate_current) {
    PetscReal max,min;
    PetscBool nonzero;
    Vec X = ksp->vec_sol;
    ierr = KSPGetInitialGuessNonzero(ksp,&nonzero);CHKERRQ(ierr);
    if (nonzero) {ierr = VecDuplicate(ksp->vec_sol,&X);CHKERRQ(ierr);}
    ierr = KSPSolve(cheb->kspest,ksp->vec_rhs,X);CHKERRQ(ierr);
    if (nonzero) {ierr = VecDestroy(&X);CHKERRQ(ierr);}
    else {ierr = VecZeroEntries(X);CHKERRQ(ierr);}
    ierr = KSPComputeExtremeSingularValues(cheb->kspest,&max,&min);CHKERRQ(ierr);
    cheb->emin = cheb->tform[0]*min + cheb->tform[1]*max;
    cheb->emax = cheb->tform[2]*min + cheb->tform[3]*max;
    cheb->estimate_current = PETSC_TRUE;
  }
	
  ksp->its = 0;
  ierr     = PCGetOperators(ksp->pc,&Amat,&Pmat,&pflag);CHKERRQ(ierr);
  maxit    = ksp->max_it;
	
  /* These three point to the three active solutions, we
	 rotate these three at each solution update */
  km1    = 0; k = 1; kp1 = 2;
  x      = ksp->vec_sol;
  b      = ksp->vec_rhs;
  p[km1] = x;
  p[k]   = ksp->work[0];
  p[kp1] = ksp->work[1];
  r      = ksp->work[2];
	
  /* use scale*B as our preconditioner */
  scale  = 2.0/(cheb->emax + cheb->emin);
	
  /*   -alpha <=  scale*lambda(B^{-1}A) <= alpha   */
  alpha  = 1.0 - scale*(cheb->emin); ;
  Gamma  = 1.0;
  mu     = 1.0/alpha; 
  omegaprod = 2.0/alpha;
	
  c[km1] = 1.0;
  c[k]   = mu;
	
  if (!ksp->guess_zero) {
    ierr = KSP_MatMult(ksp,Amat,x,r);CHKERRQ(ierr);     /*  r = b - Ax     */
    ierr = VecAYPX(r,-1.0,b);CHKERRQ(ierr);
  } else {
    ierr = VecCopy(b,r);CHKERRQ(ierr);
  }
	
  ierr = KSP_PCApply(ksp,r,p[k]);CHKERRQ(ierr);  /* p[k] = scale B^{-1}r + x */
  ierr = VecAYPX(p[k],scale,x);CHKERRQ(ierr);
	
  for (i=0; i<maxit; i++) {
    ierr = PetscObjectTakeAccess(ksp);CHKERRQ(ierr);
    ksp->its++;
    ierr = PetscObjectGrantAccess(ksp);CHKERRQ(ierr);
    c[kp1] = 2.0*mu*c[k] - c[km1];
    omega = omegaprod*c[k]/c[kp1];
		
    ierr = KSP_MatMult(ksp,Amat,p[k],r);CHKERRQ(ierr);                 /*  r = b - Ap[k]    */
    ierr = VecAYPX(r,-1.0,b);CHKERRQ(ierr);                       
    ierr = KSP_PCApply(ksp,r,p[kp1]);CHKERRQ(ierr);             /*  p[kp1] = B^{-1}z  */
		
    /* calculate residual norm if requested */
		if (ksp->lagnorm) {

			if ((i==0) || (i%ksp->chknorm == 0)) {
				
				if (ksp->normtype != KSP_NORM_NONE) {
					if (ksp->normtype == KSP_NORM_UNPRECONDITIONED) {
						if (cheb->use_red_norm) {
							ierr = KSPVecNorm_ChebychevRN(ksp,r,NORM_2,&rnorm);CHKERRQ(ierr);
						} else {
							ierr = VecNorm(r,NORM_2,&rnorm);CHKERRQ(ierr);
						}
					}
					else {
						if (cheb->use_red_norm) {
							ierr = KSPVecNorm_ChebychevRN(ksp,p[kp1],NORM_2,&rnorm);CHKERRQ(ierr);
						} else {
						ierr = VecNorm(p[kp1],NORM_2,&rnorm);CHKERRQ(ierr);
						}
					}
					ierr = PetscObjectTakeAccess(ksp);CHKERRQ(ierr);
					ksp->rnorm                              = rnorm;
					ierr = PetscObjectGrantAccess(ksp);CHKERRQ(ierr);
					ksp->vec_sol = p[k]; 
					KSPLogResidualHistory(ksp,rnorm);
					ierr = KSPMonitor(ksp,i,rnorm);CHKERRQ(ierr);
					ierr = (*ksp->converged)(ksp,i,rnorm,&ksp->reason,ksp->cnvP);CHKERRQ(ierr);
					if (ksp->reason) break;
				}
			} else {
				/* lagging */
				ierr = PetscObjectTakeAccess(ksp);CHKERRQ(ierr);
				ksp->rnorm                              = rnorm;
				ierr = PetscObjectGrantAccess(ksp);CHKERRQ(ierr);
				KSPLogResidualHistory(ksp,rnorm);
				ierr = KSPMonitor(ksp,i,rnorm);CHKERRQ(ierr);
				ierr = (*ksp->converged)(ksp,i,rnorm,&ksp->reason,ksp->cnvP);CHKERRQ(ierr);
				if (ksp->reason) break;
			}
			
		} else {
			
			if (ksp->normtype != KSP_NORM_NONE) {
				if (ksp->normtype == KSP_NORM_UNPRECONDITIONED) {
					if (cheb->use_red_norm) {
						ierr = KSPVecNorm_ChebychevRN(ksp,r,NORM_2,&rnorm);CHKERRQ(ierr);
					} else {
						ierr = VecNorm(r,NORM_2,&rnorm);CHKERRQ(ierr);
					}
				}
				else {
					if (cheb->use_red_norm) {
						ierr = KSPVecNorm_ChebychevRN(ksp,p[kp1],NORM_2,&rnorm);CHKERRQ(ierr);
					} else {
						ierr = VecNorm(p[kp1],NORM_2,&rnorm);CHKERRQ(ierr);
					}
				}
				ierr = PetscObjectTakeAccess(ksp);CHKERRQ(ierr);
				ksp->rnorm                              = rnorm;
				ierr = PetscObjectGrantAccess(ksp);CHKERRQ(ierr);
				ksp->vec_sol = p[k]; 
				KSPLogResidualHistory(ksp,rnorm);
				ierr = KSPMonitor(ksp,i,rnorm);CHKERRQ(ierr);
				ierr = (*ksp->converged)(ksp,i,rnorm,&ksp->reason,ksp->cnvP);CHKERRQ(ierr);
				if (ksp->reason) break;
			}
		}
			
		
		
    /* y^{k+1} = omega(y^{k} - y^{k-1} + Gamma*r^{k}) + y^{k-1} */
    ierr = VecScale(p[kp1],omega*Gamma*scale);CHKERRQ(ierr);
    ierr = VecAXPY(p[kp1],1.0-omega,p[km1]);CHKERRQ(ierr);
    ierr = VecAXPY(p[kp1],omega,p[k]);CHKERRQ(ierr);
		
    ktmp = km1;
    km1  = k;
    k    = kp1;
    kp1  = ktmp;
  }
	
	
  if (!ksp->reason) {
    if (ksp->normtype != KSP_NORM_NONE) {
      ierr = KSP_MatMult(ksp,Amat,p[k],r);CHKERRQ(ierr);       /*  r = b - Ap[k]    */
      ierr = VecAYPX(r,-1.0,b);CHKERRQ(ierr);
      if (ksp->normtype == KSP_NORM_UNPRECONDITIONED) {
				if (cheb->use_red_norm) {
					ierr = KSPVecNorm_ChebychevRN(ksp,r,NORM_2,&rnorm);CHKERRQ(ierr);
				} else {
					ierr = VecNorm(r,NORM_2,&rnorm);CHKERRQ(ierr);
				}
      } else {
				ierr = KSP_PCApply(ksp,r,p[kp1]);CHKERRQ(ierr); /* p[kp1] = B^{-1}z */
				if (cheb->use_red_norm) {
					ierr = KSPVecNorm_ChebychevRN(ksp,p[kp1],NORM_2,&rnorm);CHKERRQ(ierr);
				} else {
					ierr = VecNorm(p[kp1],NORM_2,&rnorm);CHKERRQ(ierr);
				}
      }
      ierr = PetscObjectTakeAccess(ksp);CHKERRQ(ierr);
      ksp->rnorm = rnorm;
      ierr = PetscObjectGrantAccess(ksp);CHKERRQ(ierr);
      ksp->vec_sol = p[k]; 
      KSPLogResidualHistory(ksp,rnorm);
      ierr = KSPMonitor(ksp,i,rnorm);CHKERRQ(ierr);
    }
		
    if (ksp->its >= ksp->max_it) {
      if (ksp->normtype != KSP_NORM_NONE) {
				ierr = (*ksp->converged)(ksp,i,rnorm,&ksp->reason,ksp->cnvP);CHKERRQ(ierr);
				if (!ksp->reason) ksp->reason = KSP_DIVERGED_ITS;
      } else { 
				ksp->reason = KSP_CONVERGED_ITS;
      }
    }
		
		
  }
	
  /* make sure solution is in vector x */
  ksp->vec_sol = x;
  if (k) {
    ierr = VecCopy(p[k],x);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSPView_ChebychevRN" 
PetscErrorCode KSPView_ChebychevRN(KSP ksp,PetscViewer viewer)
{
  KSP_ChebychevRN  *cheb = (KSP_ChebychevRN*)ksp->data;
  PetscErrorCode   ierr;
  PetscBool        iascii;
	
  PetscFunctionBegin;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  ChebychevRN: eigenvalue estimates:  min = %G, max = %G\n",cheb->emin,cheb->emax);CHKERRQ(ierr);
    if (cheb->kspest) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = KSPView(cheb->kspest,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    if (!cheb->psubcomm) {
      ierr = PetscViewerASCIIPrintf(viewer,"  Redundant VecNorm: Not activated\n");CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"  Redundant VecNorm: First (color=0) of %D processor subsets follows\n",cheb->nsubcomm);CHKERRQ(ierr);
		}			
		
  } else {
    SETERRQ1(((PetscObject)ksp)->comm,PETSC_ERR_SUP,"Viewer type %s not supported for KSP Chebychev",((PetscObject)viewer)->type_name);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSPDestroy_ChebychevRN"
PetscErrorCode KSPDestroy_ChebychevRN(KSP ksp)
{
  KSP_ChebychevRN  *cheb = (KSP_ChebychevRN*)ksp->data;
  PetscErrorCode   ierr;
	
  PetscFunctionBegin;
	if (cheb->use_red_norm) {
		ierr = VecDestroy(&cheb->xsub);CHKERRQ(ierr);
		ierr = VecDestroy(&cheb->xdup);CHKERRQ(ierr);
		ierr = VecScatterDestroy(&cheb->scatterin);CHKERRQ(ierr);
		ierr = PetscSubcommDestroy(&cheb->psubcomm);CHKERRQ(ierr);
	}
  ierr = KSPDestroy(&cheb->kspest);CHKERRQ(ierr);
  ierr = PCDestroy(&cheb->pcnone);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPChebychevRNSetEigenvalues_C","",PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPChebychevRNSetEstimateEigenvalues_C","",PETSC_NULL);CHKERRQ(ierr);
  ierr = KSPDefaultDestroy(ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*MC
 KSPCHEBYCHEV - The preconditioned Chebychev iterative method
 
 Options Database Keys:
 +   -ksp_chebychev_eigenvalues <emin,emax> - set approximations to the smallest and largest eigenvalues
 of the preconditioned operator. If these are accurate you will get much faster convergence.
 -   -ksp_chebychev_estimate_eigenvalues <a,b,c,d> - estimate eigenvalues using a Krylov method, then use this
 transform for Chebychev eigenvalue bounds (KSPChebychevSetEstimateEigenvalues)
 
 
 Level: beginner
 
 Notes: The Chebychev method requires both the matrix and preconditioner to 
 be symmetric positive (semi) definite.
 Only support for left preconditioning.
 
 .seealso:  KSPCreate(), KSPSetType(), KSPType (for list of available types), KSP,
 KSPChebychevSetEigenvalues(), KSPRICHARDSON, KSPCG
 
 M*/

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "KSPCreate_ChebychevRN"
PetscErrorCode  KSPCreate_ChebychevRN(KSP ksp)
{
  PetscErrorCode   ierr;
  KSP_ChebychevRN  *chebychevP;
  PetscMPIInt      size;
	
  PetscFunctionBegin;
  ierr = PetscNewLog(ksp,KSP_ChebychevRN,&chebychevP);CHKERRQ(ierr);
	
  ksp->data                      = (void*)chebychevP;
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_LEFT,2);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_LEFT,1);CHKERRQ(ierr);
	
  chebychevP->emin               = 1.e-2;
  chebychevP->emax               = 1.e+2;
	
  chebychevP->tform[0]           = 0.0;
  chebychevP->tform[1]           = 0.02;
  chebychevP->tform[1]           = 0;
  chebychevP->tform[2]           = 1.1;
	
	chebychevP->use_red_norm       = PETSC_FALSE;
  ierr = MPI_Comm_size(((PetscObject)ksp)->comm,&size);CHKERRQ(ierr);
  chebychevP->nsubcomm           = size;
	
  ksp->ops->setup                = KSPSetUp_ChebychevRN;
  ksp->ops->solve                = KSPSolve_ChebychevRN;
  ksp->ops->destroy              = KSPDestroy_ChebychevRN;
  ksp->ops->buildsolution        = KSPDefaultBuildSolution;
  ksp->ops->buildresidual        = KSPDefaultBuildResidual;
  ksp->ops->setfromoptions       = KSPSetFromOptions_ChebychevRN;
  ksp->ops->view                 = KSPView_ChebychevRN;
	
  ierr = PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPChebychevSetEigenvalues_C",
																					 "KSPChebychevSetEigenvalues_ChebychevRN",
																					 KSPChebychevSetEigenvalues_ChebychevRN);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunctionDynamic((PetscObject)ksp,"KSPChebychevSetEstimateEigenvalues_C",
																					 "KSPChebychevSetEstimateEigenvalues_ChebychevRN",
																					 KSPChebychevSetEstimateEigenvalues_ChebychevRN);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
EXTERN_C_END


#undef __FUNCT__  
#define __FUNCT__ "pTatinKSPRegister"
PetscErrorCode pTatinKSPRegister(void)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = KSPRegisterDynamic("chebychevrn","./","KSPCreate_ChebychevRN",KSPCreate_ChebychevRN);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

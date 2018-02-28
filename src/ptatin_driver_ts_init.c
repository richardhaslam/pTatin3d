

static const char help[] = "pTatin3d 1.0.0: Generate an initial condition \n\n";

#include "petsc/private/dmdaimpl.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "ptatin_log.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "dmda_checkpoint.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1.h"

#include "dmda_project_coords.h"
#include "stokes_operators.h"
#include "stokes_assembly.h"
#include "stokes_form_function.h"
#include "monitors.h"

#define MAX_MG_LEVELS 20

typedef enum { OP_TYPE_REDISC_ASM=0, OP_TYPE_REDISC_MF, OP_TYPE_GALERKIN } OperatorType;

typedef struct {
  PetscInt     nlevels;
  OperatorType *level_type;
  Mat          *operatorA11;
  Mat          *operatorB11;
  DM           *dav_hierarchy;
  Mat          *interpolation_v;
  Mat          *interpolation_eta;
  Quadrature   *volQ;
  BCList       *u_bclist;
  IS           *is_stokes_field;
} AuuMultiLevelCtx;


#undef __FUNCT__
#define __FUNCT__ "pTatin3dStokesBuildMeshHierarchy"
PetscErrorCode pTatin3dStokesBuildMeshHierarchy(DM dav,PetscInt nlevels,DM dav_hierarchy[])
{
  PetscErrorCode ierr;
  DM *coarsened_list;
  PetscInt k;
  
  PetscFunctionBegin;
  
  /* set up mg */
  dav->ops->coarsenhierarchy = DMCoarsenHierarchy2_DA;
  
  dav_hierarchy[ nlevels-1 ] = dav;
  ierr = PetscObjectReference((PetscObject)dav);CHKERRQ(ierr);
  
  /* Coarsen nlevels - 1 times, and add levels into list so that level 0 is the coarsest */
  ierr = PetscMalloc1(nlevels-1,&coarsened_list);CHKERRQ(ierr);
  ierr = DMCoarsenHierarchy(dav,nlevels-1,coarsened_list);CHKERRQ(ierr);
  for (k=0; k<nlevels-1; k++) {
    dav_hierarchy[ nlevels-2-k ] = coarsened_list[k];
  }
  ierr = PetscFree(coarsened_list);CHKERRQ(ierr);
  
  /* Set all dav's to be of type Q2 */
  for (k=0; k<nlevels-1; k++) {
    ierr = PetscObjectSetOptionsPrefix((PetscObject)dav_hierarchy[k],"stk_velocity_");CHKERRQ(ierr);
    ierr = DMDASetElementType_Q2(dav_hierarchy[k]);CHKERRQ(ierr);
  }
  
  /* inject coordinates */
  ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin3dStokesReportMeshHierarchy"
PetscErrorCode pTatin3dStokesReportMeshHierarchy(PetscInt nlevels,DM dav_hierarchy[])
{
  PetscErrorCode ierr;
  PetscInt       k,lmx,lmy,lmz;
  PetscMPIInt    rank,size;
  
  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  /* Report mesh sizes */
  for (k=0; k<nlevels; k++) {
    ierr = DMDAGetSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"         level [%2D]: global Q2 elements (%D x %D x %D) \n", k,lmx,lmy,lmz );
  }
  
  for (k=0; k<nlevels; k++) {
    PetscInt mp,np,pp,*_mx,*_my,*_mz,ii,jj,kk;
    
    ierr = DMDAGetOwnershipRangesElementQ2(dav_hierarchy[k],&mp,&np,&pp,NULL,NULL,NULL,&_mx,&_my,&_mz);CHKERRQ(ierr);
    
    PetscPrintf(PETSC_COMM_WORLD,"level [%2D]: [total cores %4D]: np-I [%4D]: element range I [ ", k,size,mp );
    for (ii=0; ii<mp; ii++) {
      PetscPrintf(PETSC_COMM_WORLD,"%4D", _mx[ii] );
      if (ii != mp-1) { PetscPrintf(PETSC_COMM_WORLD,", "); }
    }PetscPrintf(PETSC_COMM_WORLD," ]\n");
    
    PetscPrintf(PETSC_COMM_WORLD,"                                np-J [%4D]: element range J [ ",np);
    for (jj=0; jj<np; jj++) {
      PetscPrintf(PETSC_COMM_WORLD,"%4D", _my[jj] );
      if (jj != np-1) { PetscPrintf(PETSC_COMM_WORLD,", "); }
    }PetscPrintf(PETSC_COMM_WORLD," ]\n");
    
    PetscPrintf(PETSC_COMM_WORLD,"                                np-K [%4D]: element range K [ ",pp);
    for (kk=0; kk<pp; kk++) {
      PetscPrintf(PETSC_COMM_WORLD,"%4D", _mz[kk] );
      if (kk != pp-1) { PetscPrintf(PETSC_COMM_WORLD,", "); }
    }PetscPrintf(PETSC_COMM_WORLD," ]\n");
    
    ierr = PetscFree(_mx);CHKERRQ(ierr);
    ierr = PetscFree(_my);CHKERRQ(ierr);
    ierr = PetscFree(_mz);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin3dCreateStokesOperators"
PetscErrorCode pTatin3dCreateStokesOperators(PhysCompStokes stokes_ctx,IS is_stokes_field[],
                                             PetscInt nlevels,DM dav_hierarchy[],Mat interpolation_v[],
                                             BCList u_bclist[],Quadrature volQ[],
                                             OperatorType level_type[],
                                             Mat *_A,Mat operatorA11[],Mat *_B,Mat operatorB11[])
{
  Mat            A,B;
  DM             dap;
  PetscInt       k,max;
  PetscBool      flg,smart_defaults_activated;
  PetscInt       _level_type[MAX_MG_LEVELS];
  static int     been_here = 0;
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  dap = stokes_ctx->dap;
  
  /* A operator */
  ierr = StokesQ2P1CreateMatrix_Operator(stokes_ctx,&A);CHKERRQ(ierr);
  /* memory saving - only need daU IF you want to split A11 into A11uu,A11vv,A11ww */
  {
    MatStokesMF mf;
    
    ierr = MatShellGetMatStokesMF(A,&mf);CHKERRQ(ierr);
    ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
    mf->daU = NULL;
    ierr = ISDestroy(&mf->isU);CHKERRQ(ierr);
    ierr = ISDestroy(&mf->isV);CHKERRQ(ierr);
    ierr = ISDestroy(&mf->isW);CHKERRQ(ierr);
  }
  
  /* B operator */
  {
    Mat         Aup,Apu,Spp,bA[2][2];
    MatStokesMF StkCtx;
    
    ierr = MatShellGetMatStokesMF(A,&StkCtx);CHKERRQ(ierr);
    
    /* Schur complement */
    //ierr = DMSetMatType(dap,MATSBAIJ);CHKERRQ(ierr);
    ierr = DMCreateMatrix(dap,&Spp);CHKERRQ(ierr);
    ierr = MatSetOptionsPrefix(Spp,"S*_");CHKERRQ(ierr);
    ierr = MatSetOption(Spp,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatSetFromOptions(Spp);CHKERRQ(ierr);
    
    /* A12 */
    ierr = StokesQ2P1CreateMatrix_MFOperator_A12(StkCtx,&Aup);CHKERRQ(ierr);
    ierr = MatSetOptionsPrefix(Aup,"Bup_");CHKERRQ(ierr);
    ierr = MatSetFromOptions(Aup);CHKERRQ(ierr);
    
    /* A21 */
    ierr = StokesQ2P1CreateMatrix_MFOperator_A21(StkCtx,&Apu);CHKERRQ(ierr);
    ierr = MatSetOptionsPrefix(Apu,"Bpu_");CHKERRQ(ierr);
    ierr = MatSetFromOptions(Apu);CHKERRQ(ierr);
    
    /* nest */
    bA[0][0] = NULL; bA[0][1] = Aup;
    bA[1][0] = Apu;        bA[1][1] = Spp;
    
    ierr = MatCreateNest(PETSC_COMM_WORLD,2,is_stokes_field,2,is_stokes_field,&bA[0][0],&B);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    /* tidy up - hand back destruction to B */
    ierr = MatDestroy(&Aup);CHKERRQ(ierr);
    ierr = MatDestroy(&Apu);CHKERRQ(ierr);
    ierr = MatDestroy(&Spp);CHKERRQ(ierr);
  }
  
  /* A11 operator */
  /* defaults */
  smart_defaults_activated = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-auto_operator_type",&smart_defaults_activated,&flg);CHKERRQ(ierr);
  
  if (!smart_defaults_activated) {
    _level_type[0] = (PetscInt)OP_TYPE_REDISC_ASM;
    for (k=1; k<nlevels; k++) {
      _level_type[k] = (PetscInt)OP_TYPE_REDISC_MF;
    }
  } else {
    PetscInt b,lmx,lmy,lmz;
    
    // rule 1 // all grids > 20^3 per rank -> use MF, otherwise assemble
    // rule 2 // if rule invoked, use ASM on next coarsest, followed by Galerkin

    for (k=0; k<nlevels-1; k++) {
      _level_type[k] = (PetscInt)OP_TYPE_GALERKIN;
    }
    _level_type[nlevels-1] = (PetscInt)OP_TYPE_REDISC_ASM;

    // rule 1 //
    b = nlevels;
    for (k=1; k<nlevels; k++) {
      ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
      if ((lmx > 20) || (lmy > 20) || (lmz > 20)) {
        _level_type[k-1] = (PetscInt)OP_TYPE_REDISC_ASM;
        _level_type[k] = (PetscInt)OP_TYPE_REDISC_MF;
        b = k;
        break;
      }
    }
    for (k=b; k<nlevels; k++) {
      _level_type[k] = (PetscInt)OP_TYPE_REDISC_MF;
    }
  }
  
  max = nlevels;
  ierr = PetscOptionsGetIntArray(NULL,NULL,"-A11_operator_type",_level_type,&max,&flg);CHKERRQ(ierr);
  if (flg) {
    if (max != nlevels) SETERRQ2(PetscObjectComm((PetscObject)dap),PETSC_ERR_USER,"Incorrect number of values provided to option -A11_operator_type. Expected %D values, found %D values",nlevels,max);
  }
  for (k=nlevels-1; k>=0; k--) {
    level_type[k] = (OperatorType)_level_type[k];
  }
  for (k=nlevels-1; k>=0; k--) {
    
    switch (level_type[k]) {
        
      case OP_TYPE_REDISC_ASM:
      {
        Mat Auu;
        PetscBool same1 = PETSC_FALSE,same2 = PETSC_FALSE,same3 = PETSC_FALSE;
        Vec X;
        MatNullSpace nullsp;
        
        /* use -stk_velocity_da_mat_type sbaij or -Buu_da_mat_type sbaij */
        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Re-discretisation :: assembled operator \n", k);
        //ierr = DMSetMatType(dav_hierarchy[k],MATSBAIJ);CHKERRQ(ierr);
        ierr = DMCreateMatrix(dav_hierarchy[k],&Auu);CHKERRQ(ierr);
        ierr = MatSetOptionsPrefix(Auu,"Buu_");CHKERRQ(ierr);
        ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSBAIJ,&same1);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSEQSBAIJ,&same2);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATMPISBAIJ,&same3);CHKERRQ(ierr);
        if (same1||same2||same3) {
          ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
        }
        /* should move assembly into jacobian */
        ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
        ierr = MatAssemble_StokesA_AUU(Auu,dav_hierarchy[k],u_bclist[k],volQ[k]);CHKERRQ(ierr);
        
        operatorA11[k] = Auu;
        operatorB11[k] = Auu;
        ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
        ierr = DMGetCoordinates(dav_hierarchy[k],&X);CHKERRQ(ierr);
        ierr = MatNullSpaceCreateRigidBody(X,&nullsp);CHKERRQ(ierr);
        ierr = MatSetNearNullSpace(Auu,nullsp);CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(&nullsp);CHKERRQ(ierr);
      }
        break;
        
      case OP_TYPE_REDISC_MF:
      {
        Mat Auu;
        MatA11MF mf,A11Ctx;
        
        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Re-discretisation :: matrix free operator \n", k);
        ierr = MatA11MFCreate(&A11Ctx);CHKERRQ(ierr);
        ierr = MatA11MFSetup(A11Ctx,dav_hierarchy[k],volQ[k],u_bclist[k]);CHKERRQ(ierr);
        
        ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
        /* memory saving - only need daU IF you want to split A11 into A11uu,A11vv,A11ww */
        ierr = MatShellGetMatA11MF(Auu,&mf);CHKERRQ(ierr);
        ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
        mf->daU = NULL;
        ierr = ISDestroy(&mf->isU);CHKERRQ(ierr);
        ierr = ISDestroy(&mf->isV);CHKERRQ(ierr);
        ierr = ISDestroy(&mf->isW);CHKERRQ(ierr);
        /* --- */
        operatorA11[k] = Auu;
        
        operatorB11[k] = Auu;
        ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
        
        ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
      }
        break;
        
      case OP_TYPE_GALERKIN:
      {
        Mat Auu;
        Vec X;
        MatNullSpace nullsp;
        
        if (k==nlevels-1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Galerkin :: assembled operator \n", k);
        
        /* should move coarse grid assembly into jacobian */
        ierr = MatPtAP(operatorA11[k+1],interpolation_v[k+1],MAT_INITIAL_MATRIX,1.0,&Auu);CHKERRQ(ierr);
        
        operatorA11[k] = Auu;
        operatorB11[k] = Auu;
        ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
        ierr = DMGetCoordinates(dav_hierarchy[k],&X);CHKERRQ(ierr);
        ierr = MatNullSpaceCreateRigidBody(X,&nullsp);CHKERRQ(ierr);
        ierr = MatSetBlockSize(Auu,3);CHKERRQ(ierr);
        ierr = MatSetNearNullSpace(Auu,nullsp);CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(&nullsp);CHKERRQ(ierr);
      }
        break;
        
      default:
        break;
    }
  }	
  
  /* Set fine A11 into nest */
  ierr = MatNestSetSubMat(B,0,0,operatorA11[nlevels-1]);CHKERRQ(ierr);
  
  *_A = A;
  *_B = B;
  
  been_here = 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HMG_SetUp"
PetscErrorCode HMG_SetUp(AuuMultiLevelCtx *mlctx, pTatinCtx user)
{
  DM             dav_hierarchy[MAX_MG_LEVELS];
  Mat            interpolation_v[MAX_MG_LEVELS],interpolation_eta[MAX_MG_LEVELS];
  PetscInt       k,nlevels;
  Quadrature     volQ[MAX_MG_LEVELS];
  BCList         u_bclist[MAX_MG_LEVELS];
  DM             dmstokes,dav;
  IS             *is_stokes_field;
  PhysCompStokes stokes = NULL;
  pTatinModel     model = NULL;
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dav,NULL);CHKERRQ(ierr);
  
  ierr = DMCompositeGetGlobalISs(dmstokes,&is_stokes_field);CHKERRQ(ierr);

  nlevels = 1;
  PetscOptionsGetInt(NULL,NULL,"-dau_nlevels",&nlevels,NULL);
  PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%D x %D x %D) : MG levels %D  \n", user->mx,user->my,user->mz,nlevels );
  ierr = pTatin3dStokesBuildMeshHierarchy(dav,nlevels,dav_hierarchy);CHKERRQ(ierr);
  ierr = pTatin3dStokesReportMeshHierarchy(nlevels,dav_hierarchy);CHKERRQ(ierr);
  ierr = pTatinLogNote(user,"  [Velocity multi-grid hierarchy]");CHKERRQ(ierr);
  for (k=nlevels-1; k>=0; k--) {
    char name[PETSC_MAX_PATH_LEN];
    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"vel_dmda_Lv%D",k);
    ierr = pTatinLogBasicDMDA(user,name,dav_hierarchy[k]);CHKERRQ(ierr);
  }
  
  /* Coarse grid setup: Define interpolation operators for velocity space */
  interpolation_v[0] = NULL;
  for (k=0; k<nlevels-1; k++) {
    ierr = DMCreateInterpolation(dav_hierarchy[k],dav_hierarchy[k+1],&interpolation_v[k+1],NULL);CHKERRQ(ierr);
  }
  
  /* Coarse grid setup: Define interpolation operators for scalar space */
  interpolation_eta[0] = NULL;
  for (k=1; k<nlevels; k++) {
    ierr = MatMAIJRedimension(interpolation_v[k],1,&interpolation_eta[k]);CHKERRQ(ierr);
  }
  PetscPrintf(PETSC_COMM_WORLD,"Generated velocity mesh hierarchy --> ");
  pTatinGetRangeCurrentMemoryUsage(NULL);
  
  /* Coarse grid setup: Define material properties on gauss points */
  for (k=0; k<nlevels-1; k++) {
    PetscInt ncells,lmx,lmy,lmz;
    PetscInt np_per_dim;
    
    np_per_dim = 3;
    ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
    ncells = lmx * lmy * lmz;
    ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&volQ[k]);CHKERRQ(ierr);
  }
  volQ[nlevels-1] = stokes->volQ;
  PetscPrintf(PETSC_COMM_WORLD,"Generated quadrature point hierarchy --> ");
  pTatinGetRangeCurrentMemoryUsage(NULL);
  
  /* Coarse grid setup: Define boundary conditions */
  for (k=0; k<nlevels-1; k++) {
    ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);
  }
  u_bclist[nlevels-1] = stokes->u_bclist;
  
  /* Coarse grid setup: Configure boundary conditions */
  ierr = pTatinModel_ApplyBoundaryConditionMG(nlevels,u_bclist,dav_hierarchy,model,user);CHKERRQ(ierr);
  
  /* set all pointers into mg context */
  mlctx->is_stokes_field     = is_stokes_field;
  
  ierr = PetscMalloc1(nlevels,&mlctx->dav_hierarchy);CHKERRQ(ierr);
  ierr = PetscMalloc1(nlevels,&mlctx->interpolation_v);CHKERRQ(ierr);
  ierr = PetscMalloc1(nlevels,&mlctx->interpolation_eta);CHKERRQ(ierr);
  ierr = PetscMalloc1(nlevels,&mlctx->volQ);CHKERRQ(ierr);
  ierr = PetscMalloc1(nlevels,&mlctx->u_bclist);CHKERRQ(ierr);
  mlctx->nlevels                  = nlevels;
  for (k=0; k<nlevels; k++) {
    mlctx->dav_hierarchy[k]       = dav_hierarchy[k];
    mlctx->interpolation_v[k]     = interpolation_v[k];
    mlctx->interpolation_eta[k]   = interpolation_eta[k];
    mlctx->volQ[k]                = volQ[k];
    mlctx->u_bclist[k]            = u_bclist[k];
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HMGOperator_SetUp"
PetscErrorCode HMGOperator_SetUp(AuuMultiLevelCtx *mlctx,pTatinCtx user,Mat *A,Mat *B)
{
  OperatorType   level_type[MAX_MG_LEVELS];
  Mat            operatorA11[MAX_MG_LEVELS],operatorB11[MAX_MG_LEVELS];
  PhysCompStokes stokes = NULL;
  PetscInt       k,nlevels;
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);

  /* configure stokes opertors */
  ierr = pTatin3dCreateStokesOperators(stokes,mlctx->is_stokes_field,
                                       mlctx->nlevels,
                                       mlctx->dav_hierarchy,
                                       mlctx->interpolation_v,
                                       mlctx->u_bclist,
                                       mlctx->volQ,
                                       level_type,
                                       A,operatorA11,B,operatorB11);CHKERRQ(ierr);

  nlevels = mlctx->nlevels;
  ierr = PetscMalloc1(nlevels,&mlctx->level_type);CHKERRQ(ierr);
  ierr = PetscMalloc1(nlevels,&mlctx->operatorA11);CHKERRQ(ierr);
  ierr = PetscMalloc1(nlevels,&mlctx->operatorB11);CHKERRQ(ierr);
  for (k=0; k<nlevels; k++) {
    mlctx->level_type[k]  = level_type[k];
    mlctx->operatorA11[k] = operatorA11[k];
    mlctx->operatorB11[k] = operatorB11[k];
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESComposeWithMGCtx"
PetscErrorCode SNESComposeWithMGCtx(SNES snes,AuuMultiLevelCtx *mgctx)
{
  PetscErrorCode ierr;
  PetscContainer container;
  
  PetscFunctionBegin;
  ierr = PetscContainerCreate(PetscObjectComm((PetscObject)snes),&container);CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(container,(void*)mgctx);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)snes,"AuuMultiLevelCtx",(PetscObject)container);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESDestroyMGCtx"
PetscErrorCode SNESDestroyMGCtx(SNES snes)
{
  PetscErrorCode ierr;
  PetscContainer container;
  
  PetscFunctionBegin;
  container = NULL;
  ierr = PetscObjectQuery((PetscObject)snes,"AuuMultiLevelCtx",(PetscObject*)&container);CHKERRQ(ierr);
  if (container) {
    ierr = PetscContainerDestroy(&container);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin3dStokesKSPConfigureFSGMG"
PetscErrorCode pTatin3dStokesKSPConfigureFSGMG(KSP ksp,PetscInt nlevels,Mat operatorA11[],Mat operatorB11[],Mat interpolation_v[],DM dav_hierarchy[])
{
  PetscInt k,nsplits;
  PC       pc,pc_i;
  KSP      *sub_ksp,ksp_coarse,ksp_smoother;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
  
  ierr = KSPSetDM(sub_ksp[0],dav_hierarchy[nlevels-1]);CHKERRQ(ierr);
  ierr = KSPSetDMActive(sub_ksp[0],PETSC_FALSE);CHKERRQ(ierr);
  
  ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
  ierr = PCSetType(pc_i,PCMG);CHKERRQ(ierr);
  ierr = PCMGSetLevels(pc_i,nlevels,NULL);CHKERRQ(ierr);
  ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
  ierr = PCMGSetGalerkin(pc_i,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PCSetDM(pc_i,NULL);CHKERRQ(ierr);
  
  for( k=1; k<nlevels; k++ ){
    ierr = PCMGSetInterpolation(pc_i,k,interpolation_v[k]);CHKERRQ(ierr);
  }
  
  /* drop the operators in - i presume this will also need to be performed inside the jacobian each time the operators are modified */
  /* No - it looks like PCSetUp_MG will call set operators on all levels if the SetOperators was called on the finest, which should/is done by the SNES */
  ierr = PCMGGetCoarseSolve(pc_i,&ksp_coarse);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp_coarse,operatorA11[0],operatorA11[0]);CHKERRQ(ierr);
  
  ierr = KSPSetDM(ksp_coarse,dav_hierarchy[0]);CHKERRQ(ierr);
  ierr = KSPSetDMActive(ksp_coarse,PETSC_FALSE);CHKERRQ(ierr);
  
  for( k=1; k<nlevels; k++ ){
    ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);
    
    // use A for smoother, B for residual
    ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k]);CHKERRQ(ierr);
    ierr = KSPSetDM(ksp_smoother,dav_hierarchy[k]);CHKERRQ(ierr);
    ierr = KSPSetDMActive(ksp_smoother,PETSC_FALSE);CHKERRQ(ierr);
  }
  PetscFree(sub_ksp);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormJacobian_StokesMGAuu"
PetscErrorCode FormJacobian_StokesMGAuu(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  pTatinCtx         user;
  AuuMultiLevelCtx  *mlctx;
  DM                stokes_pack,dau,dap;
  PhysCompStokes    stokes;
  Vec               Uloc,Ploc;
  PetscScalar       *LA_Uloc,*LA_Ploc;
  PetscBool         is_mffd = PETSC_FALSE;
  PetscBool         is_nest = PETSC_FALSE;
  PetscBool         is_shell = PETSC_FALSE;
  PetscContainer    container;
  PetscInt          k;
  PetscErrorCode    ierr;
  
  PetscFunctionBegin;
  
  user = (pTatinCtx)ctx;
  
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;
  
  ierr = PetscObjectQuery((PetscObject)snes,"AuuMultiLevelCtx",(PetscObject*)&container);CHKERRQ(ierr);
  if (!container) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"No data with name \"AuuMultiLevelCtx\" was composed with SNES");
  ierr = PetscContainerGetPointer(container,(void*)&mlctx);CHKERRQ(ierr);
  
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
  
  ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
  ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
  
  /* nonlinearitiers: markers => quad points */
  ierr = pTatin_EvaluateRheologyNonlinearities(user,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
  
  /* interpolate coefficients */
  {
    int               npoints;
    DataField         PField_std;
    DataField         PField_stokes;
    MPntStd           *mp_std;
    MPntPStokes       *mp_stokes;
    
    DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
    DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
    
    DataBucketGetSizes(user->materialpoint_db,&npoints,NULL,NULL);
    mp_std    = PField_std->data; /* should write a function to do this */
    mp_stokes = PField_stokes->data; /* should write a function to do this */
    
    ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(user->coefficient_projection_type,npoints,mp_std,mp_stokes,mlctx->nlevels,mlctx->interpolation_eta,mlctx->dav_hierarchy,mlctx->volQ);CHKERRQ(ierr);
  }
  
  /* operator */
  ierr = PetscObjectTypeCompare((PetscObject)A,MATMFFD, &is_mffd);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)A,MATNEST, &is_nest);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)A,MATSHELL,&is_shell);CHKERRQ(ierr);
  
  if (is_nest) {
    Mat Auu;
    
    ierr = MatGetSubMatrix(A,mlctx->is_stokes_field[0],mlctx->is_stokes_field[0],MAT_INITIAL_MATRIX,&Auu);CHKERRQ(ierr);
    
    is_shell = PETSC_FALSE;
    ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSHELL,&is_shell);CHKERRQ(ierr);
    if (!is_shell) {
      ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
      ierr = MatAssemble_StokesA_AUU(Auu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
    }
    
    ierr = MatDestroy(&Auu);CHKERRQ(ierr);
  }
  /* If shell, do nothing */
  /* If mffd,  do nothing */
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  /* preconditioner operator for Jacobian */
  {
    Mat Buu,Bpp;
    
    ierr = MatGetSubMatrix(B,mlctx->is_stokes_field[0],mlctx->is_stokes_field[0],MAT_INITIAL_MATRIX,&Buu);CHKERRQ(ierr);
    ierr = MatGetSubMatrix(B,mlctx->is_stokes_field[1],mlctx->is_stokes_field[1],MAT_INITIAL_MATRIX,&Bpp);CHKERRQ(ierr);
    
    is_shell = PETSC_FALSE;
    ierr = PetscObjectTypeCompare((PetscObject)Buu,MATSHELL,&is_shell);CHKERRQ(ierr);
    if (!is_shell) {
      ierr = MatZeroEntries(Buu);CHKERRQ(ierr);
      ierr = MatAssemble_StokesA_AUU(Buu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
    }
    
    is_shell = PETSC_FALSE;
    ierr = PetscObjectTypeCompare((PetscObject)Bpp,MATSHELL,&is_shell);CHKERRQ(ierr);
    if (!is_shell) {
      ierr = MatZeroEntries(Bpp);CHKERRQ(ierr);
      ierr = MatAssemble_StokesPC_ScaledMassMatrix(Bpp,dau,dap,user->stokes_ctx->p_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
    }
    
    ierr = MatDestroy(&Buu);CHKERRQ(ierr);
    ierr = MatDestroy(&Bpp);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  /* Buu preconditioner for all other levels in the hierarchy */
  {
    KSP       ksp,*sub_ksp,ksp_smoother;
    PC        pc,pc_i;
    PetscInt  nsplits;
    
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
    
    for (k=mlctx->nlevels-2; k>=0; k--) {
      /* fetch smoother */
      if (k == 0) { ierr = PCMGGetCoarseSolve(pc_i,&ksp_smoother);CHKERRQ(ierr); }
      else {        ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr); }
      
      switch (mlctx->level_type[k]) {
          
        case OP_TYPE_REDISC_ASM:
          ierr = MatZeroEntries(mlctx->operatorB11[k]);CHKERRQ(ierr);
          ierr = MatAssemble_StokesA_AUU(mlctx->operatorB11[k],mlctx->dav_hierarchy[k],mlctx->u_bclist[k],mlctx->volQ[k]);CHKERRQ(ierr);
          
          ierr = KSPSetOperators(ksp_smoother,mlctx->operatorB11[k],mlctx->operatorB11[k]);CHKERRQ(ierr);
          /* hack for nested coarse solver */
          {
            KSP ksp_nested;
            PC pc_smoother;
            PetscBool is_nested_ksp;
            
            ierr = KSPGetPC(ksp_smoother,&pc_smoother);CHKERRQ(ierr);
            is_nested_ksp = PETSC_FALSE;
            ierr = PetscObjectTypeCompare((PetscObject)pc_smoother,PCKSP,&is_nested_ksp);CHKERRQ(ierr);
            if (is_nested_ksp) {
              ierr = PCKSPGetKSP(pc_smoother,&ksp_nested);CHKERRQ(ierr);
              ierr = KSPSetOperators(ksp_nested,mlctx->operatorB11[k],mlctx->operatorB11[k]);CHKERRQ(ierr);
            }
          }
          break;
          
        case OP_TYPE_REDISC_MF:
      //	ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k]);CHKERRQ(ierr);
          ierr = KSPSetOperators(ksp_smoother,mlctx->operatorA11[k],mlctx->operatorA11[k]);CHKERRQ(ierr);
          break;
          
        case OP_TYPE_GALERKIN:
        {
          Mat Auu_k;
          
          if (k == mlctx->nlevels-1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
          PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Galerkin :: assembled operator \n", k);
          
          /*
           ierr = MatPtAP(mlctx->operatorA11[k+1],mlctx->interpolation_v[k+1],MAT_INITIAL_MATRIX,1.0,&Auu_k);CHKERRQ(ierr);
           ierr = KSPSetOperators(ksp_smoother,Auu_k,Auu_k);CHKERRQ(ierr);
           mlctx->operatorA11[k] = Auu_k;
           mlctx->operatorB11[k] = Auu_k;
           ierr = PetscObjectReference((PetscObject)Auu_k);CHKERRQ(ierr);
           */
          Auu_k = mlctx->operatorA11[k];
          ierr = MatPtAP(mlctx->operatorA11[k+1],mlctx->interpolation_v[k+1],MAT_REUSE_MATRIX,1.0,&Auu_k);CHKERRQ(ierr);
          ierr = KSPSetOperators(ksp_smoother,Auu_k,Auu_k);CHKERRQ(ierr);
          mlctx->operatorB11[k] = Auu_k;
        }
          break;
          
        default:
          break;
      }
    }
    PetscFree(sub_ksp);
    
    /* push operators */
    for (k=mlctx->nlevels-1; k>=0; k--) {
      
      ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
      ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
      ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
      
      if (k == 0) {  ierr = PCMGGetCoarseSolve(pc_i,&ksp_smoother);CHKERRQ(ierr); }
      else {         ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr); }
      
      switch (mlctx->level_type[k]) {
          
        case OP_TYPE_REDISC_ASM:
          ierr = KSPSetOperators(ksp_smoother,mlctx->operatorB11[k],mlctx->operatorB11[k]);CHKERRQ(ierr);
          /* hack for nested coarse solver */
          {
            KSP ksp_nested;
            PC pc_smoother;
            PetscBool is_nested_ksp;
            
            ierr = KSPGetPC(ksp_smoother,&pc_smoother);CHKERRQ(ierr);
            is_nested_ksp = PETSC_FALSE;
            ierr = PetscObjectTypeCompare((PetscObject)pc_smoother,PCKSP,&is_nested_ksp);CHKERRQ(ierr);
            if (is_nested_ksp) {
              ierr = PCKSPGetKSP(pc_smoother,&ksp_nested);CHKERRQ(ierr);
              ierr = KSPSetOperators(ksp_nested,mlctx->operatorB11[k],mlctx->operatorB11[k]);CHKERRQ(ierr);
            }
          }
          break;
          
        case OP_TYPE_REDISC_MF:
          ierr = KSPSetOperators(ksp_smoother,mlctx->operatorA11[k],mlctx->operatorA11[k]);CHKERRQ(ierr);
          break;
          
        case OP_TYPE_GALERKIN:
          ierr = KSPSetOperators(ksp_smoother,mlctx->operatorA11[k],mlctx->operatorB11[k]);CHKERRQ(ierr);
          break;
      }
      PetscFree(sub_ksp);
    }
				
  }
  
  {
    PetscBool mg_dump_coarse = PETSC_FALSE;
    char filename[PETSC_MAX_PATH_LEN];
    PetscInt snes_it;
    PetscViewer viewer;
    
    PetscOptionsGetBool(NULL,NULL,"-ptatin_mg_dump_coarse_operator",&mg_dump_coarse,0);
    SNESGetIterationNumber(snes,&snes_it);
    
    if (mg_dump_coarse) {
      if (mlctx->level_type[0] != OP_TYPE_REDISC_MF) {
        
        
        PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/mg_coarse_operatorA_step%D_snes%D.mat",user->outputpath,user->step,snes_it);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
        MatView(mlctx->operatorA11[0],viewer);
        PetscViewerDestroy(&viewer);
        if (mlctx->operatorA11[0] != mlctx->operatorB11[0]) {
          PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/mg_coarse_operatorB_step%D_snes%D.mat",user->outputpath,user->step,snes_it);
          PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
          MatView(mlctx->operatorB11[0],viewer);
          PetscViewerDestroy(&viewer);
        }
      }
    }
  }
  
  
  /* clean up */
  ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
  
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinNonlinearStokesSolveCreate(pTatinCtx user,Mat A,Mat B,Vec F,AuuMultiLevelCtx *mgctx,SNES *s)
{
  SNES snes;
  KSP  ksp;
  PC pc;
  PetscErrorCode ierr;
  
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,A,B,FormJacobian_StokesMGAuu,user);CHKERRQ(ierr);
  
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = SNESComposeWithMGCtx(snes,mgctx);CHKERRQ(ierr);
  ierr = pTatin_Stokes_ActivateMonitors(user,snes);CHKERRQ(ierr);
  
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCFIELDSPLIT);CHKERRQ(ierr);
  ierr = PCFieldSplitSetIS(pc,"u",mgctx->is_stokes_field[0]);CHKERRQ(ierr);
  ierr = PCFieldSplitSetIS(pc,"p",mgctx->is_stokes_field[1]);CHKERRQ(ierr);
  
  /* configure uu split for galerkin multi-grid */
  ierr = pTatin3dStokesKSPConfigureFSGMG(ksp,mgctx->nlevels,mgctx->operatorA11,mgctx->operatorB11,mgctx->interpolation_v,mgctx->dav_hierarchy);CHKERRQ(ierr);
  ierr = pTatinLogBasic(user);CHKERRQ(ierr);

  *s = snes;
  
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinNonlinearStokesSolve(pTatinCtx user,SNES snes,Vec X,const char stagename[])
{
  PetscLogDouble time[2];
  PhysCompStokes stokes = NULL;
  DM             dmstokes;
  char           title[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"   --------- Stokes[%s] ---------\n",stagename);
  PetscTime(&time[0]);
  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
  PetscTime(&time[1]);
  PetscSNPrintf(title,PETSC_MAX_PATH_LEN-1,"Stokes[%s]",stagename);
  ierr = pTatinLogBasicSNES(user,   title,snes);CHKERRQ(ierr);
  ierr = pTatinLogBasicCPUtime(user,title,time[1]-time[0]);CHKERRQ(ierr);
  ierr = pTatinLogBasicStokesSolution(user,dmstokes,X);CHKERRQ(ierr);
  ierr = pTatinLogBasicStokesSolutionResiduals(user,snes,dmstokes,X);CHKERRQ(ierr);
  ierr = pTatinLogPetscLog(user,    title);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GenerateICStateFromModelDefinition"
PetscErrorCode GenerateICStateFromModelDefinition(void)
{
  pTatinCtx       user;
  pTatinModel     model = NULL;
  PhysCompStokes  stokes = NULL;
  PhysCompEnergy  energy = NULL;
  DM              dmstokes,dmv,dmp,dmenergy = NULL;
  Vec             X_s,X_e = NULL;
  Vec             velocity,pressure;
  PetscBool       activate_energy = PETSC_FALSE;
  DataBucket      materialpoint_db,material_constants_db;
  PetscErrorCode  ierr;
  
  PetscFunctionBegin;
  ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
  ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);

  /* driver specific options parsed here */

  /* Register all models */
  ierr = pTatinModelLoad(user);CHKERRQ(ierr);
  ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);

  ierr = pTatinModel_Initialize(model,user);CHKERRQ(ierr);
  PetscOptionsGetBool(NULL,NULL,"-activate_energy",&activate_energy,NULL);

  /* Create Stokes context */
  ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,&dmp);CHKERRQ(ierr);
  
  { /* IF I DON'T DO THIS, THE IS's OBTAINED FROM DMCompositeGetGlobalISs() are wrong !! */
    Vec X;
    ierr = DMGetGlobalVector(dmstokes,&X);CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(dmstokes,&X);CHKERRQ(ierr);
  }
  /* Pack all physics together */
  ierr = PetscObjectReference((PetscObject)dmstokes);CHKERRQ(ierr);
  user->pack = dmstokes;
  
  /* Create material points */
  ierr = pTatin3dCreateMaterialPoints(user,dmv);CHKERRQ(ierr);
  ierr = pTatinGetMaterialPoints(user,&materialpoint_db,NULL);CHKERRQ(ierr);

  /* Create energy context */
  /* NOTE - Calling pTatinPhysCompActivate_Energy() after pTatin3dCreateMaterialPoints() is essential */
  if (activate_energy) {
    ierr = pTatinPhysCompActivate_Energy(user,PETSC_TRUE);CHKERRQ(ierr);
    ierr = pTatinGetContext_Energy(user,&energy);CHKERRQ(ierr);
    dmenergy = energy->daT;
  }
  
  /* mesh geometry */
  ierr = pTatinModel_ApplyInitialMeshGeometry(model,user);CHKERRQ(ierr);
  if (activate_energy) {
    ierr = DMDAProjectCoordinatesQ2toQ1(dmv,dmenergy,energy->energy_mesh_type);CHKERRQ(ierr);
  }
  
  /* interpolate material point coordinates (needed if mesh was modified) */
  ierr = MaterialPointCoordinateSetUp(user,dmv);CHKERRQ(ierr);
  
  /* material geometry */
  ierr = pTatinModel_ApplyInitialMaterialGeometry(model,user);CHKERRQ(ierr);

  /* work vector for solution */
  ierr = DMCreateGlobalVector(dmstokes,&X_s);CHKERRQ(ierr);
  if (activate_energy) {
		ierr = DMCreateGlobalVector(dmenergy,&X_e);CHKERRQ(ierr);
		ierr = pTatinPhysCompAttachData_Energy(user,X_e,NULL);CHKERRQ(ierr);
  }
  
  /* initial condition */
  ierr = pTatinModel_ApplyInitialSolution(model,user,X_s);CHKERRQ(ierr);
  
  /* initial viscosity  */
  ierr = pTatinModel_ApplyInitialStokesVariableMarkers(model,user,X_s);CHKERRQ(ierr);
  
  /* boundary conditions */
  ierr = pTatinModel_ApplyBoundaryCondition(model,user);CHKERRQ(ierr);

  /* insert boundary conditions into solution vector */
  {
    ierr = DMCompositeGetAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);
    ierr = BCListInsert(stokes->u_bclist,velocity);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);
    
    if (activate_energy) {
      ierr = BCListInsert(energy->T_bclist,X_e);CHKERRQ(ierr);
    }
  }
  
  ierr = pTatinGetMaterialConstants(user,&material_constants_db);CHKERRQ(ierr);
  
  /* Configure for the initial condition */
  user->step = 0;
  user->time = 0.0;
  user->dt = 1.0e-10;
  
#if 0
  {
    AuuMultiLevelCtx mgctx;
    Mat A,B;
    Vec F;
    SNES snes;
    KSP  ksp;
    PC pc;
    PetscInt snes_its;
    PetscLogDouble time[2];
    
    ierr = HMG_SetUp(&mgctx,user);CHKERRQ(ierr);
    ierr = HMGOperator_SetUp(&mgctx,user,&A,&B);CHKERRQ(ierr);
    
    ierr = VecDuplicate(X_s,&F);CHKERRQ(ierr);
    
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
    ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,A,B,FormJacobian_StokesMGAuu,user);CHKERRQ(ierr);

    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
    ierr = SNESComposeWithMGCtx(snes,&mgctx);CHKERRQ(ierr);
    ierr = pTatin_Stokes_ActivateMonitors(user,snes);CHKERRQ(ierr);

    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCFIELDSPLIT);CHKERRQ(ierr);
    ierr = PCFieldSplitSetIS(pc,"u",mgctx.is_stokes_field[0]);CHKERRQ(ierr);
    ierr = PCFieldSplitSetIS(pc,"p",mgctx.is_stokes_field[1]);CHKERRQ(ierr);
    
    /* configure uu split for galerkin multi-grid */
    ierr = pTatin3dStokesKSPConfigureFSGMG(ksp,mgctx.nlevels,mgctx.operatorA11,mgctx.operatorB11,mgctx.interpolation_v,mgctx.dav_hierarchy);CHKERRQ(ierr);
    ierr = pTatinLogBasic(user);CHKERRQ(ierr);

    ierr = SNESGetTolerances(snes,0,0,0,&snes_its,0);CHKERRQ(ierr);
    ierr = SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
    
    PetscPrintf(PETSC_COMM_WORLD,"   --------- LINEAR STAGE ---------\n");
    PetscTime(&time[0]);
    ierr = SNESSolve(snes,NULL,X_s);CHKERRQ(ierr);
    PetscTime(&time[1]);
    ierr = pTatinLogBasicSNES(user,   "Stokes[LinearStage]",snes);CHKERRQ(ierr);
    ierr = pTatinLogBasicCPUtime(user,"Stokes[LinearStage]",time[1]-time[0]);CHKERRQ(ierr);
    ierr = pTatinLogBasicStokesSolution(user,dmstokes,X_s);CHKERRQ(ierr);
    ierr = pTatinLogBasicStokesSolutionResiduals(user,snes,dmstokes,X_s);CHKERRQ(ierr);
    ierr = pTatinLogPetscLog(user,    "Stokes[LinearStage]");CHKERRQ(ierr);
    
    ierr = VecDestroy(&F);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  }
#endif
  
  
#if 1
  {
    AuuMultiLevelCtx mgctx;
    Mat A,B;
    Vec F;
    SNES snes;
    SNESLineSearch linesearch;
    KSP  ksp;
    
    ierr = HMG_SetUp(&mgctx,user);CHKERRQ(ierr);
    ierr = HMGOperator_SetUp(&mgctx,user,&A,&B);CHKERRQ(ierr);
    
    ierr = VecDuplicate(X_s,&F);CHKERRQ(ierr);

    ierr = pTatinNonlinearStokesSolveCreate(user,A,B,F,&mgctx,&snes);CHKERRQ(ierr);
    
    /* configure as a linear solve */
    //ierr = SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = SNESSetType(snes,SNESNEWTONLS);CHKERRQ(ierr);
    ierr = SNESGetLineSearch(snes,&linesearch);CHKERRQ(ierr);
    ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC);CHKERRQ(ierr);
    ierr = SNESSetType(snes,SNESKSPONLY);CHKERRQ(ierr);
    //ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    //ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);

    ierr = pTatinNonlinearStokesSolve(user,snes,X_s,"Linear Stage");CHKERRQ(ierr);
  }
#endif
  
  
  
  {
    char output_path[PETSC_MAX_PATH_LEN];
    char output_path_ic[PETSC_MAX_PATH_LEN];
    
    ierr = PetscSNPrintf(output_path,PETSC_MAX_PATH_LEN-1,"%s",user->outputpath);CHKERRQ(ierr);
    ierr = PetscSNPrintf(output_path_ic,PETSC_MAX_PATH_LEN-1,"%s/step%d",user->outputpath,user->step);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(output_path_ic);CHKERRQ(ierr);
    ierr = PetscSNPrintf(user->outputpath,PETSC_MAX_PATH_LEN-1,"%s",output_path_ic);CHKERRQ(ierr);

    ierr = pTatinModel_Output(model,user,X_s,NULL);CHKERRQ(ierr);

    ierr = PetscSNPrintf(user->outputpath,PETSC_MAX_PATH_LEN-1,"%s",output_path);CHKERRQ(ierr);
  }

#if 0
  /* probably this should be in a stand alone function */
  {
    char g_checkpoint_path[PETSC_MAX_PATH_LEN];
    char checkpoint_path[PETSC_MAX_PATH_LEN];
    char checkpoint_prefix[PETSC_MAX_PATH_LEN];
    
    ierr = PetscSNPrintf(g_checkpoint_path,PETSC_MAX_PATH_LEN-1,"%s/checkpoints",user->outputpath);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(g_checkpoint_path);CHKERRQ(ierr);

    ierr = PetscSNPrintf(checkpoint_path,PETSC_MAX_PATH_LEN-1,"%s/checkpoints/intitial_condition",user->outputpath);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(checkpoint_path);CHKERRQ(ierr);

    ierr = PetscSNPrintf(checkpoint_prefix,PETSC_MAX_PATH_LEN-1,"%s/stokes_v",checkpoint_path);CHKERRQ(ierr);
    ierr = DMDACheckpointWrite(dmv,checkpoint_prefix);CHKERRQ(ierr);

    ierr = DMCompositeGetAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);
    ierr = PetscSNPrintf(checkpoint_prefix,PETSC_MAX_PATH_LEN-1,"%s/ptatinstate_stokes_Xv.pbvec",checkpoint_path);CHKERRQ(ierr);
    ierr = DMDAWriteVectorToFile(velocity,checkpoint_prefix,PETSC_FALSE);CHKERRQ(ierr);
    ierr = PetscSNPrintf(checkpoint_prefix,PETSC_MAX_PATH_LEN-1,"%s/ptatinstate_stokes_Xp.pbvec",checkpoint_path);CHKERRQ(ierr);
    ierr = DMDAWriteVectorToFile(pressure,checkpoint_prefix,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);

    
    if (activate_energy) {
      ierr = PhysCompCheckpointWrite_Energy(energy,PETSC_TRUE,checkpoint_path,NULL);CHKERRQ(ierr);

      ierr = PetscSNPrintf(checkpoint_prefix,PETSC_MAX_PATH_LEN-1,"%s/ptatinstate_energy_Xt.pbvec",checkpoint_path);CHKERRQ(ierr);
      ierr = DMDAWriteVectorToFile(X_e,checkpoint_prefix,PETSC_FALSE);CHKERRQ(ierr);
    }
    
    ierr = PetscSNPrintf(checkpoint_prefix,PETSC_MAX_PATH_LEN-1,"%s/materialpoint",checkpoint_path);CHKERRQ(ierr);
    DataBucketView(PETSC_COMM_WORLD,materialpoint_db,checkpoint_prefix,DATABUCKET_VIEW_NATIVE);

    ierr = PetscSNPrintf(checkpoint_prefix,PETSC_MAX_PATH_LEN-1,"%s/materialconstants",checkpoint_path);CHKERRQ(ierr);
    
    /* material_constants_db is a redundant object, e.g. it is identical on all ranks */
    /* Hence, we let only 1 rank write out the data file during checkpoint.write() */
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
    if (rank == 0) {
      DataBucketView(PETSC_COMM_SELF,material_constants_db,checkpoint_prefix,DATABUCKET_VIEW_NATIVE);
    }
    
    /* For checkpoint.read() we do a redundant read */
    DataBucket newdb;
    DataBucketLoadRedundant_NATIVE(PETSC_COMM_WORLD,checkpoint_prefix,&newdb);
    DataBucketView(PETSC_COMM_SELF,newdb,"mat-constants-red",DATABUCKET_VIEW_STDOUT);
  }
#endif

  
  
  
  
  /* last thing we do */
  {
    char checkpoints_path[PETSC_MAX_PATH_LEN];
    char checkpoint_path[PETSC_MAX_PATH_LEN];
    
    ierr = PetscSNPrintf(checkpoints_path,PETSC_MAX_PATH_LEN-1,"%s/checkpoints",user->outputpath);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(checkpoints_path);CHKERRQ(ierr);
    
    ierr = PetscSNPrintf(checkpoint_path,PETSC_MAX_PATH_LEN-1,"%s/intitial_condition",checkpoints_path);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(checkpoint_path);CHKERRQ(ierr);
    
    ierr = pTatinCtxCheckpointWrite(user,checkpoint_path,NULL,
                                    dmstokes,dmenergy,
                                    0,NULL,NULL,
                                    X_s,X_e,NULL,NULL);CHKERRQ(ierr);
    
  }
  
  ierr = VecDestroy(&X_e);CHKERRQ(ierr);
  ierr = VecDestroy(&X_s);CHKERRQ(ierr);
  ierr = pTatin3dDestroyContext(&user);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LoadICStateFromModelDefinition"
PetscErrorCode LoadICStateFromModelDefinition(pTatinCtx *pctx,Vec *v1,Vec *v2,PetscBool write_checkpoint)
{
  pTatinCtx       user;
  pTatinModel     model = NULL;
  PhysCompStokes  stokes = NULL;
  PhysCompEnergy  energy = NULL;
  DM              dmstokes,dmv,dmp,dmenergy = NULL;
  Vec             X_s,X_e = NULL;
  PetscBool       activate_energy = PETSC_FALSE;
  DataBucket      materialpoint_db,material_constants_db;
  PetscErrorCode  ierr;
  
  PetscFunctionBegin;
  ierr = pTatin3dLoadContext_FromFile(&user);CHKERRQ(ierr);
  ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
  
  /* driver specific options parsed here */
  
  /* Register all models */
  ierr = pTatinModelLoad(user);CHKERRQ(ierr);
  ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);
  
  ierr = pTatinModel_Initialize(model,user);CHKERRQ(ierr);
  PetscOptionsGetBool(NULL,NULL,"-activate_energy",&activate_energy,NULL);
  
  /* Create Stokes context */
  ierr = pTatin3d_PhysCompStokesLoad_FromFile(user);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,&dmp);CHKERRQ(ierr);

  { /* IF I DON'T DO THIS, THE IS's OBTAINED FROM DMCompositeGetGlobalISs() are wrong !! */
    Vec X;
    ierr = DMGetGlobalVector(dmstokes,&X);CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(dmstokes,&X);CHKERRQ(ierr);
  }
  /* Pack all physics together */
  ierr = PetscObjectReference((PetscObject)dmstokes);CHKERRQ(ierr);
  user->pack = dmstokes;
  
  /* Create material points */
  // material point load //ierr = pTatin3dCreateMaterialPoints(user,dmv);CHKERRQ(ierr);
  ierr = pTatin3dLoadMaterialPoints_FromFile(user,dmv);CHKERRQ(ierr);
  ierr = pTatinGetMaterialPoints(user,&materialpoint_db,NULL);CHKERRQ(ierr);

  /* Create energy context */
  /* NOTE - Calling pTatinPhysCompActivate_Energy() after pTatin3dCreateMaterialPoints() is not essential when restarting */
  /* The reason for this is pTatinPhysCompActivate_Energy_FromFile() does not register new fields, rather DataBucketLoad() does */
  if (activate_energy) {
    ierr = pTatinPhysCompActivate_Energy_FromFile(user);CHKERRQ(ierr);
    ierr = pTatinGetContext_Energy(user,&energy);CHKERRQ(ierr);
    dmenergy = energy->daT;
  }
  
  /* mesh geometry */
  // stokes dmda load // ierr = pTatinModel_ApplyInitialMeshGeometry(model,user);CHKERRQ(ierr);
  if (activate_energy) {
    // probably undeeded as well // ierr = DMDAProjectCoordinatesQ2toQ1(dmv,dmenergy,energy->energy_mesh_type);CHKERRQ(ierr);
  }
  
  /* interpolate material point coordinates (needed if mesh was modified) */
  // material point load //ierr = MaterialPointCoordinateSetUp(user,dmv);CHKERRQ(ierr);
  
  /* material geometry */
  // material point load //ierr = pTatinModel_ApplyInitialMaterialGeometry(model,user);CHKERRQ(ierr);
  
  /* work vector for solution */
  ierr = DMCreateGlobalVector(dmstokes,&X_s);CHKERRQ(ierr);
  if (activate_energy) {
    ierr = DMCreateGlobalVector(dmenergy,&X_e);CHKERRQ(ierr);
    ierr = pTatinPhysCompAttachData_Energy(user,X_e,NULL);CHKERRQ(ierr);
  }
  
  /* initial condition - call user method, then clobber */
  ierr = pTatinModel_ApplyInitialSolution(model,user,X_s);CHKERRQ(ierr);
  ierr = pTatin3dLoadState_FromFile(user,dmstokes,dmenergy,X_s,X_e);CHKERRQ(ierr);
  
  /* initial viscosity  */
  // material point load //ierr = pTatinModel_ApplyInitialStokesVariableMarkers(model,user,X_s);CHKERRQ(ierr);

  /* boundary conditions */
  ierr = pTatinModel_ApplyBoundaryCondition(model,user);CHKERRQ(ierr);
  
  /* insert boundary conditions into solution vector */
  // pTatin3dLoadState_FromFile() //
  
  ierr = pTatinGetMaterialConstants(user,&material_constants_db);CHKERRQ(ierr);
  
  {
    char output_path[PETSC_MAX_PATH_LEN];
    char output_path_ic[PETSC_MAX_PATH_LEN];
    
    ierr = PetscSNPrintf(output_path,PETSC_MAX_PATH_LEN-1,"%s",user->outputpath);CHKERRQ(ierr);
    ierr = PetscSNPrintf(output_path_ic,PETSC_MAX_PATH_LEN-1,"%s/fromfile",user->outputpath);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(output_path_ic);CHKERRQ(ierr);
    ierr = PetscSNPrintf(user->outputpath,PETSC_MAX_PATH_LEN-1,"%s",output_path_ic);CHKERRQ(ierr);
    
    ierr = pTatinModel_Output(model,user,X_s,"icbc");CHKERRQ(ierr);
    
    ierr = PetscSNPrintf(user->outputpath,PETSC_MAX_PATH_LEN-1,"%s",output_path);CHKERRQ(ierr);
  }
  
  if (write_checkpoint) {
    char checkpoints_path[PETSC_MAX_PATH_LEN];
    char checkpoint_path[PETSC_MAX_PATH_LEN];
    
    ierr = PetscSNPrintf(checkpoints_path,PETSC_MAX_PATH_LEN-1,"%s/checkpoints",user->outputpath);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(checkpoints_path);CHKERRQ(ierr);
    
    ierr = PetscSNPrintf(checkpoint_path,PETSC_MAX_PATH_LEN-1,"%s/intitial_condition",checkpoints_path);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(checkpoint_path);CHKERRQ(ierr);
    
    ierr = pTatinCtxCheckpointWrite(user,checkpoint_path,NULL,
                                    dmstokes,dmenergy,
                                    0,NULL,NULL,
                                    X_s,X_e,NULL,NULL);CHKERRQ(ierr);
  }

  if (v1) { *v2 = X_e; }
  else    { ierr = VecDestroy(&X_e);CHKERRQ(ierr); }
  if (v1) { *v1 = X_s; }
  else    { ierr = VecDestroy(&X_s);CHKERRQ(ierr); }

  *pctx = user;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DummyRun"
PetscErrorCode DummyRun(pTatinCtx pctx,Vec v1,Vec v2)
{
  PetscInt k,step0;
  Vec Xs = NULL,Xe = NULL;
  PetscBool energy_activated;
  PhysCompStokes  stokes = NULL;
  PhysCompEnergy  energy = NULL;
  DM dmstokes,dmenergy = NULL;
  PetscErrorCode ierr;
  pTatinModel     model = NULL;

  ierr = pTatinGetModel(pctx,&model);CHKERRQ(ierr);

  ierr = pTatinGetStokesContext(pctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);

  if (v1) {
    Xs = v1;
  } else {
    ierr = DMCreateGlobalVector(dmstokes,&Xs);CHKERRQ(ierr);
  }

  ierr = pTatinContextValid_Energy(pctx,&energy_activated);CHKERRQ(ierr);
  if (energy_activated) {
    ierr = pTatinGetContext_Energy(pctx,&energy);CHKERRQ(ierr);
    dmenergy = energy->daT;
    ierr = DMCreateGlobalVector(dmenergy,&Xe);CHKERRQ(ierr);
    
    if (v2) {
      Xe = v2;
    } else {
      ierr = DMCreateGlobalVector(dmenergy,&Xe);CHKERRQ(ierr);
    }
  }
  
  step0 = pctx->step + 1;
  for (k=step0; k<=pctx->nsteps; k++) {
    PetscReal dt_curr,dt_next;
    
    /* use last saved time step for any calculations, e.g. update particles */
    dt_curr = pctx->dt;
    
    printf("  [executing time step %3d] : curr time %+1.4e : advancing time by dt %+1.8e \n",k,pctx->time,dt_curr);
    
    /* solve mechanics + energy */
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /* compute new time step */
    dt_next = 0.0;
    //dt_next = 200.0 * rand()/((PetscReal)RAND_MAX);
    dt_next = 0.001245 * ((PetscReal)k) + 1.0e-4;
    
    if (dt_next > pctx->dt_max) { dt_next = pctx->dt_max; }
    if (dt_next < pctx->dt_min) { dt_next = pctx->dt_min; }
    
    pctx->step = k;
    pctx->time += dt_curr;
    pctx->dt = dt_next;
    if (energy_activated) {
      energy->time = pctx->time;
      energy->dt = pctx->dt;
    }
    ierr = pTatinLogBasic(pctx);CHKERRQ(ierr);
    
    {
      char io_path[PETSC_MAX_PATH_LEN],op[PETSC_MAX_PATH_LEN];
      
      ierr = PetscSNPrintf(op,PETSC_MAX_PATH_LEN-1,"%s",pctx->outputpath);CHKERRQ(ierr);
      
      ierr = PetscSNPrintf(io_path,PETSC_MAX_PATH_LEN-1,"%s/step%d",pctx->outputpath,k);CHKERRQ(ierr);
      ierr = pTatinCreateDirectory(io_path);CHKERRQ(ierr);

      ierr = PetscSNPrintf(pctx->outputpath,PETSC_MAX_PATH_LEN-1,"%s",io_path);CHKERRQ(ierr);
      ierr = pTatinModel_Output(model,pctx,Xs,NULL);CHKERRQ(ierr);

      ierr = PetscSNPrintf(pctx->outputpath,PETSC_MAX_PATH_LEN-1,"%s",op);CHKERRQ(ierr);
    }
    
    {
      char checkpoints_path[PETSC_MAX_PATH_LEN];
      char checkpoint_path[PETSC_MAX_PATH_LEN];
      
      ierr = PetscSNPrintf(checkpoints_path,PETSC_MAX_PATH_LEN-1,"%s/checkpoints",pctx->outputpath);CHKERRQ(ierr);
      ierr = pTatinCreateDirectory(checkpoints_path);CHKERRQ(ierr);
      
      ierr = PetscSNPrintf(checkpoint_path,PETSC_MAX_PATH_LEN-1,"%s/step%d",checkpoints_path,k);CHKERRQ(ierr);
      ierr = pTatinCreateDirectory(checkpoint_path);CHKERRQ(ierr);
      
      ierr = pTatinCtxCheckpointWrite(pctx,checkpoint_path,NULL,
                                      dmstokes,dmenergy,0,NULL,NULL,Xs,Xe,NULL,NULL);CHKERRQ(ierr);
    }
    
    if (pctx->time > pctx->time_max) break;
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char *argv[])
{
  PetscErrorCode ierr;
  PetscBool init = PETSC_FALSE,load = PETSC_FALSE, run = PETSC_FALSE;
  pTatinCtx pctx = NULL;
  
  ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
  ierr = pTatinModelRegisterAll();CHKERRQ(ierr);

  // linear only
  // non-linear with picard
  // isoastasy

  PetscOptionsGetBool(NULL,NULL,"-init",&init,NULL);
  if (init) {
    ierr = GenerateICStateFromModelDefinition();CHKERRQ(ierr);
  }
  
  PetscOptionsGetBool(NULL,NULL,"-load",&load,NULL);
  if (load) {
    ierr = LoadICStateFromModelDefinition(&pctx,NULL,NULL,PETSC_TRUE);CHKERRQ(ierr);
  }

  PetscOptionsGetBool(NULL,NULL,"-run",&run,NULL);
  if (run) {
    Vec Xup,Xt;
    ierr = LoadICStateFromModelDefinition(&pctx,&Xup,&Xt,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DummyRun(pctx,Xup,Xt);CHKERRQ(ierr);
    ierr = VecDestroy(&Xup);CHKERRQ(ierr);
    ierr = VecDestroy(&Xt);CHKERRQ(ierr);
  }
  
  ///ierr = linearSolve();CHKERRQ(ierr);

  //ierr = picardSolve();CHKERRQ(ierr);

  //ierr = dumpState();CHKERRQ(ierr);

  if (pctx) { ierr = pTatin3dDestroyContext(&pctx); }
  
  ierr = pTatinModelDeRegisterAll();CHKERRQ(ierr);
  ierr = pTatinFinalize();CHKERRQ(ierr);
  return 0;
}



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
#include "energy_assembly.h"
#include "dmda_checkpoint.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1.h"

#include "dmda_project_coords.h"
#include "stokes_operators.h"
#include "stokes_assembly.h"
#include "stokes_form_function.h"
#include "monitors.h"
#include "mp_advection.h"
#include "mesh_update.h"
#include "material_point_popcontrol.h"

#define MAX_MG_LEVELS 20

typedef enum { OP_TYPE_REDISC_ASM=0, OP_TYPE_REDISC_MF, OP_TYPE_GALERKIN, OP_TYPE_UDEF } OperatorType;

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


PetscErrorCode pTatin3dStokesBuildMeshHierarchy(DM dav,PetscInt nlevels,DM dav_hierarchy[])
{
  PetscErrorCode ierr;
  DM             *coarsened_list;
  PetscInt       k;

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
      if (ii != (mp-1)) { PetscPrintf(PETSC_COMM_WORLD,", "); }
    } PetscPrintf(PETSC_COMM_WORLD," ]\n");

    PetscPrintf(PETSC_COMM_WORLD,"                                np-J [%4D]: element range J [ ",np);
    for (jj=0; jj<np; jj++) {
      PetscPrintf(PETSC_COMM_WORLD,"%4D", _my[jj] );
      if (jj != (np-1)) { PetscPrintf(PETSC_COMM_WORLD,", "); }
    } PetscPrintf(PETSC_COMM_WORLD," ]\n");

    PetscPrintf(PETSC_COMM_WORLD,"                                np-K [%4D]: element range K [ ",pp);
    for (kk=0; kk<pp; kk++) {
      PetscPrintf(PETSC_COMM_WORLD,"%4D", _mz[kk] );
      if (kk != (pp-1)) { PetscPrintf(PETSC_COMM_WORLD,", "); }
    } PetscPrintf(PETSC_COMM_WORLD," ]\n");

    ierr = PetscFree(_mx);CHKERRQ(ierr);
    ierr = PetscFree(_my);CHKERRQ(ierr);
    ierr = PetscFree(_mz);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

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
    bA[1][0] = Apu;  bA[1][1] = Spp;

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
        Mat          Auu;
        PetscBool    same1 = PETSC_FALSE,same2 = PETSC_FALSE,same3 = PETSC_FALSE;
        Vec          X;
        MatNullSpace nullsp;

        /* use -stk_velocity_da_mat_type sbaij or -Buu_da_mat_type sbaij */
        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Re-discretisation :: assembled operator \n",k);
        //ierr = DMSetMatType(dav_hierarchy[k],MATSBAIJ);CHKERRQ(ierr);
        ierr = DMCreateMatrix(dav_hierarchy[k],&Auu);CHKERRQ(ierr);
        ierr = MatSetOptionsPrefix(Auu,"Buu_");CHKERRQ(ierr);
        ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSBAIJ,&same1);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSEQSBAIJ,&same2);CHKERRQ(ierr);
        ierr = PetscObjectTypeCompare((PetscObject)Auu,MATMPISBAIJ,&same3);CHKERRQ(ierr);
        if (same1 || same2 || same3) {
          ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
        }
        /* should move assembly into jacobian */
        //ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
        //ierr = MatAssemble_StokesA_AUU(Auu,dav_hierarchy[k],u_bclist[k],volQ[k]);CHKERRQ(ierr);

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
        Mat      Auu;
        MatA11MF mf,A11Ctx;

        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Re-discretisation :: matrix free operator \n",k);
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
        Mat          Auu;
        Vec          X;
        MatNullSpace nullsp;

        if (k == (nlevels-1)) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
        if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%D]: Coarse grid type :: Galerkin :: assembled operator \n",k);

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

      case OP_TYPE_UDEF:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported operator type OP_TYPE_UDEF");

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
  pTatinModel    model = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dav,NULL);CHKERRQ(ierr);

  ierr = DMCompositeGetGlobalISs(dmstokes,&is_stokes_field);CHKERRQ(ierr);

  nlevels = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-dau_nlevels",&nlevels,NULL);CHKERRQ(ierr);
  if (nlevels >= MAX_MG_LEVELS) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Maximum number of multi-grid levels is set by #define MAX_MG_LEVELS %D",MAX_MG_LEVELS);

  PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%D x %D x %D) : MG levels %D  \n",user->mx,user->my,user->mz,nlevels);
  ierr = pTatin3dStokesBuildMeshHierarchy(dav,nlevels,dav_hierarchy);CHKERRQ(ierr);
  ierr = pTatin3dStokesReportMeshHierarchy(nlevels,dav_hierarchy);CHKERRQ(ierr);
  ierr = pTatinLogNote(user,"  [Velocity multi-grid hierarchy]");CHKERRQ(ierr);
  for (k=nlevels-1; k>=0; k--) {
    char name[PETSC_MAX_PATH_LEN];

    ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"vel_dmda_Lv%D",k);CHKERRQ(ierr);
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
  /*
  PetscPrintf(PETSC_COMM_WORLD,"Generated velocity mesh hierarchy --> ");
  pTatinGetRangeCurrentMemoryUsage(NULL);
  */

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
  /*
  PetscPrintf(PETSC_COMM_WORLD,"Generated quadrature point hierarchy --> ");
  pTatinGetRangeCurrentMemoryUsage(NULL);
  */

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

  ierr = PetscMalloc1(nlevels,&mlctx->level_type);CHKERRQ(ierr);
  ierr = PetscMalloc1(nlevels,&mlctx->operatorA11);CHKERRQ(ierr);
  ierr = PetscMalloc1(nlevels,&mlctx->operatorB11);CHKERRQ(ierr);

  for (k=0; k<nlevels; k++) {
    mlctx->dav_hierarchy[k]       = NULL;
    mlctx->interpolation_v[k]     = NULL;
    mlctx->interpolation_eta[k]   = NULL;
    mlctx->volQ[k]                = NULL;
    mlctx->u_bclist[k]            = NULL;

    mlctx->level_type[k]          = OP_TYPE_UDEF;
    mlctx->operatorA11[k]         = NULL;
    mlctx->operatorB11[k]         = NULL;
  }

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

PetscErrorCode HMG_Destroy(AuuMultiLevelCtx *mlctx)
{
  PetscInt       k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (k=0; k<mlctx->nlevels-1; k++) {
    ierr = BCListDestroy(&mlctx->u_bclist[k]);CHKERRQ(ierr);
    ierr = QuadratureDestroy(&mlctx->volQ[k]);CHKERRQ(ierr);
  }
  for (k=0; k<mlctx->nlevels; k++) {
    ierr = MatDestroy(&mlctx->operatorA11[k]);CHKERRQ(ierr);
    ierr = MatDestroy(&mlctx->operatorB11[k]);CHKERRQ(ierr);
    ierr = MatDestroy(&mlctx->interpolation_v[k]);CHKERRQ(ierr);
    ierr = MatDestroy(&mlctx->interpolation_eta[k]);CHKERRQ(ierr);
    ierr = DMDestroy(&mlctx->dav_hierarchy[k]);CHKERRQ(ierr);
  }

  ierr = ISDestroy(&mlctx->is_stokes_field[0]);CHKERRQ(ierr);
  ierr = ISDestroy(&mlctx->is_stokes_field[1]);CHKERRQ(ierr);
  ierr = PetscFree(mlctx->is_stokes_field);CHKERRQ(ierr);

  ierr = PetscFree(mlctx->level_type);CHKERRQ(ierr);
  ierr = PetscFree(mlctx->operatorB11);CHKERRQ(ierr);
  ierr = PetscFree(mlctx->operatorA11);CHKERRQ(ierr);
  ierr = PetscFree(mlctx->dav_hierarchy);CHKERRQ(ierr);
  ierr = PetscFree(mlctx->interpolation_v);CHKERRQ(ierr);
  ierr = PetscFree(mlctx->interpolation_eta);CHKERRQ(ierr);
  ierr = PetscFree(mlctx->volQ);CHKERRQ(ierr);
  ierr = PetscFree(mlctx->u_bclist);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode HMGOperator_Destroy(AuuMultiLevelCtx *mlctx)
{
  PetscInt       k,nlevels;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  nlevels = mlctx->nlevels;
  for (k=0; k<nlevels; k++) {
    ierr = MatDestroy(&mlctx->operatorA11[k]);CHKERRQ(ierr);
    ierr = MatDestroy(&mlctx->operatorB11[k]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

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
  for (k=0; k<nlevels; k++) {
    mlctx->level_type[k]  = level_type[k];
    mlctx->operatorA11[k] = operatorA11[k];
    mlctx->operatorB11[k] = operatorB11[k];
  }
  PetscFunctionReturn(0);
}

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

PetscErrorCode pTatin3dStokesKSPConfigureFSGMG(KSP ksp,PetscInt nlevels,Mat operatorA11[],Mat operatorB11[],Mat interpolation_v[],DM dav_hierarchy[])
{
  PetscInt       k,nsplits;
  PC             pc,pc_i;
  KSP            *sub_ksp,ksp_coarse,ksp_smoother;
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
  ierr = PCMGSetGalerkin(pc_i,PC_MG_GALERKIN_NONE);CHKERRQ(ierr);
  ierr = PCSetDM(pc_i,NULL);CHKERRQ(ierr);

  for (k=1; k<nlevels; k++) {
    ierr = PCMGSetInterpolation(pc_i,k,interpolation_v[k]);CHKERRQ(ierr);
  }

  /* drop the operators in - i presume this will also need to be performed inside the jacobian each time the operators are modified */
  /* No - it looks like PCSetUp_MG will call set operators on all levels if the SetOperators was called on the finest, which should/is done by the SNES */
  ierr = PCMGGetCoarseSolve(pc_i,&ksp_coarse);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp_coarse,operatorA11[0],operatorA11[0]);CHKERRQ(ierr);

  ierr = KSPSetDM(ksp_coarse,dav_hierarchy[0]);CHKERRQ(ierr);
  ierr = KSPSetDMActive(ksp_coarse,PETSC_FALSE);CHKERRQ(ierr);

  for (k=1; k<nlevels; k++) {
    ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);

    // use A for smoother, B for residual
    ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k]);CHKERRQ(ierr);
    ierr = KSPSetDM(ksp_smoother,dav_hierarchy[k]);CHKERRQ(ierr);
    ierr = KSPSetDMActive(ksp_smoother,PETSC_FALSE);CHKERRQ(ierr);
  }
  ierr = PetscFree(sub_ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

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

    ierr = MatCreateSubMatrix(A,mlctx->is_stokes_field[0],mlctx->is_stokes_field[0],MAT_INITIAL_MATRIX,&Auu);CHKERRQ(ierr);

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

    ierr = MatCreateSubMatrix(B,mlctx->is_stokes_field[0],mlctx->is_stokes_field[0],MAT_INITIAL_MATRIX,&Buu);CHKERRQ(ierr);
    ierr = MatCreateSubMatrix(B,mlctx->is_stokes_field[1],mlctx->is_stokes_field[1],MAT_INITIAL_MATRIX,&Bpp);CHKERRQ(ierr);

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
            KSP       ksp_nested;
            PC        pc_smoother;
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
      //  ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k]);CHKERRQ(ierr);
          ierr = KSPSetOperators(ksp_smoother,mlctx->operatorA11[k],mlctx->operatorA11[k]);CHKERRQ(ierr);
          break;

        case OP_TYPE_GALERKIN:
        {
          Mat Auu_k;

          if (k == (mlctx->nlevels-1)) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
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

        case OP_TYPE_UDEF:
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported operator type OP_TYPE_UDEF");

        default:
          break;
      }
    }
    ierr = PetscFree(sub_ksp);CHKERRQ(ierr);

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
            KSP       ksp_nested;
            PC        pc_smoother;
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

        case OP_TYPE_UDEF:
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported operator type OP_TYPE_UDEF");

        case OP_TYPE_GALERKIN:
          ierr = KSPSetOperators(ksp_smoother,mlctx->operatorA11[k],mlctx->operatorB11[k]);CHKERRQ(ierr);
          break;
      }
      ierr = PetscFree(sub_ksp);CHKERRQ(ierr);
    }

  }

  {
    PetscBool   mg_dump_coarse = PETSC_FALSE;
    char        filename[PETSC_MAX_PATH_LEN];
    PetscInt    snes_it;
    PetscViewer viewer;

    ierr = PetscOptionsGetBool(NULL,NULL,"-ptatin_mg_dump_coarse_operator",&mg_dump_coarse,0);CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(snes,&snes_it);CHKERRQ(ierr);

    if (mg_dump_coarse) {
      if (mlctx->level_type[0] != OP_TYPE_REDISC_MF) {
        ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/mg_coarse_operatorA_step%D_snes%D.mat",user->outputpath,user->step,snes_it);CHKERRQ(ierr);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
        ierr = MatView(mlctx->operatorA11[0],viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
        if (mlctx->operatorA11[0] != mlctx->operatorB11[0]) {
          ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/mg_coarse_operatorB_step%D_snes%D.mat",user->outputpath,user->step,snes_it);CHKERRQ(ierr);
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
          ierr = MatView(mlctx->operatorB11[0],viewer);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
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
  SNES           snes;
  KSP            ksp;
  PC             pc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
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

  PetscFunctionBegin;
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);

  /* insert boundary conditions into solution vector */
  {
    Vec Xu,Xp;
    
    ierr = DMCompositeGetAccess(dmstokes,X,&Xu,&Xp);CHKERRQ(ierr);
    ierr = BCListInsert(stokes->u_bclist,Xu);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(dmstokes,X,&Xu,&Xp);CHKERRQ(ierr);
  }

  if (stagename) {
    PetscPrintf(PETSC_COMM_WORLD,"   --------- Stokes[%s] ---------\n",stagename);
  }
  PetscTime(&time[0]);
  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
  PetscTime(&time[1]);
  if (stagename) {
    ierr = PetscSNPrintf(title,PETSC_MAX_PATH_LEN-1,"Stokes[%s]",stagename);CHKERRQ(ierr);
  } else {
    ierr = PetscSNPrintf(title,PETSC_MAX_PATH_LEN-1,"Stokes");CHKERRQ(ierr);
  }
  ierr = pTatinLogBasicSNES(user,   title,snes);CHKERRQ(ierr);
  ierr = pTatinLogBasicCPUtime(user,title,time[1]-time[0]);CHKERRQ(ierr);
  ierr = pTatinLogBasicStokesSolution(user,dmstokes,X);CHKERRQ(ierr);
  ierr = pTatinLogBasicStokesSolutionResiduals(user,snes,dmstokes,X);CHKERRQ(ierr);
  //ierr = pTatinLogPetscLog(user,    title);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode GenerateICStateFromModelDefinition(pTatinCtx *pctx)
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
  PetscReal       surface_displacement_max = 1.0e32;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
  ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);

  /* driver specific options parsed here */
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt_max_surface_displacement",&surface_displacement_max,NULL);CHKERRQ(ierr);

  /* Register all models */
  ierr = pTatinModelLoad(user);CHKERRQ(ierr);
  ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);

  ierr = pTatinModel_Initialize(model,user);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-activate_energy",&activate_energy,NULL);CHKERRQ(ierr);

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
    ierr = pTatinLogBasicDMDA(user,"energy_dmda",dmenergy);CHKERRQ(ierr);
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
  ierr = DMCompositeGetAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);
  ierr = BCListInsert(stokes->u_bclist,velocity);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);

  if (activate_energy) {
    ierr = BCListInsert(energy->T_bclist,X_e);CHKERRQ(ierr);
  }

  ierr = pTatinGetMaterialConstants(user,&material_constants_db);CHKERRQ(ierr);

  /* Configure for the initial condition */
  user->step = 0;
  user->time = 0.0;
  user->dt = 1.0e-10;

  {
    AuuMultiLevelCtx mgctx;
    Mat              A,B;
    Vec              F;
    SNES             snes;
    SNESLineSearch   linesearch;
    //KSP              ksp;
    RheologyType     init_rheology_type;

    ierr = VecDuplicate(X_s,&F);CHKERRQ(ierr);
    ierr = HMG_SetUp(&mgctx,user);CHKERRQ(ierr);

    /* linear stage */
    ierr = HMGOperator_SetUp(&mgctx,user,&A,&B);CHKERRQ(ierr);
    ierr = pTatinNonlinearStokesSolveCreate(user,A,B,F,&mgctx,&snes);CHKERRQ(ierr);

    /* configure as a linear solve */
    init_rheology_type = user->rheology_constants.rheology_type;
    user->rheology_constants.rheology_type = RHEOLOGY_VISCOUS;
    ierr = SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);CHKERRQ(ierr);

    ierr = SNESSetType(snes,SNESNEWTONLS);CHKERRQ(ierr);
    ierr = SNESGetLineSearch(snes,&linesearch);CHKERRQ(ierr);
    ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC);CHKERRQ(ierr);
    //ierr = SNESSetType(snes,SNESKSPONLY);CHKERRQ(ierr);
    //ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    //ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    //ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);

    ierr = pTatinNonlinearStokesSolve(user,snes,X_s,"Linear Stage");CHKERRQ(ierr);

    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = MatDestroy(&B);CHKERRQ(ierr);
    ierr = SNESDestroyMGCtx(snes);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);
    ierr = HMGOperator_Destroy(&mgctx);CHKERRQ(ierr);

    /* picard stage */
    ierr = HMGOperator_SetUp(&mgctx,user,&A,&B);CHKERRQ(ierr);
    ierr = pTatinNonlinearStokesSolveCreate(user,A,B,F,&mgctx,&snes);CHKERRQ(ierr);

    user->rheology_constants.rheology_type = init_rheology_type;
    ierr = pTatinNonlinearStokesSolve(user,snes,X_s,"Picard Stage");CHKERRQ(ierr);

    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = MatDestroy(&B);CHKERRQ(ierr);
    ierr = SNESDestroyMGCtx(snes);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);
    ierr = HMGOperator_Destroy(&mgctx);CHKERRQ(ierr);

    ierr = HMG_Destroy(&mgctx);CHKERRQ(ierr);
    ierr = VecDestroy(&F);CHKERRQ(ierr);
  }

  /* compute timestep */
  user->dt = 1.0e32;
  {
    Vec       velocity,pressure;
    PetscReal timestep;

    ierr = DMCompositeGetAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);
    ierr = SwarmUpdatePosition_ComputeCourantStep(dmv,velocity,&timestep);CHKERRQ(ierr);
    ierr = pTatin_SetTimestep(user,"StkCourant",timestep);CHKERRQ(ierr);

    ierr = UpdateMeshGeometry_ComputeSurfaceCourantTimestep(dmv,velocity,surface_displacement_max,&timestep);CHKERRQ(ierr);
    ierr = pTatin_SetTimestep(user,"StkSurfaceCourant",timestep);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"  timestep[stokes] dt_courant = %1.4e \n", user->dt );
  }

  /* first time step, enforce to be super small */
  user->dt = user->dt * 1.0e-10;

  /* initialise the energy solver */
  if (activate_energy) {
    PetscReal timestep;

    ierr = pTatinPhysCompEnergy_Initialise(energy,X_e);CHKERRQ(ierr);

    /* first time this is called we REQUIRE that a valid time step is chosen */
    energy->dt = user->dt;
    ierr = pTatinPhysCompEnergy_UpdateALEVelocity(stokes,X_s,energy,energy->dt);CHKERRQ(ierr);
    ierr = pTatinPhysCompEnergy_ComputeTimestep(energy,energy->Told,&timestep);CHKERRQ(ierr);

    /*
     Note - we cannot use the time step for energy equation here.
     It seems silly, but to compute the adf-diff time step, we need to the ALE velocity,
     however to compute the ALE velocity we need to know the timestep.
     */
    PetscPrintf(PETSC_COMM_WORLD,"  timestep[adv-diff] dt_courant = %1.4e \n", timestep );
    energy->dt   = user->dt;
  }

  {
    char prefix[PETSC_MAX_PATH_LEN];

    ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step%D",user->step);CHKERRQ(ierr);
    ierr = pTatinModel_Output(model,user,X_s,prefix);CHKERRQ(ierr);
  }

  /* last thing we do */
  {
    char           checkpoints_path[PETSC_MAX_PATH_LEN];
    char           checkpoint_path[PETSC_MAX_PATH_LEN];
    char           filename[PETSC_MAX_PATH_LEN];
    PetscLogDouble time[2];

    ierr = PetscSNPrintf(checkpoints_path,PETSC_MAX_PATH_LEN-1,"%s/checkpoints",user->outputpath);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(checkpoints_path);CHKERRQ(ierr);

    ierr = PetscSNPrintf(checkpoint_path,PETSC_MAX_PATH_LEN-1,"%s/initial_condition",checkpoints_path);CHKERRQ(ierr);
    ierr = pTatinCreateDirectory(checkpoint_path);CHKERRQ(ierr);

    PetscTime(&time[0]);
    ierr = pTatinCtxCheckpointWrite(user,checkpoint_path,NULL,
                                    dmstokes,dmenergy,
                                    0,NULL,NULL,
                                    X_s,X_e,NULL,NULL);CHKERRQ(ierr);
    PetscTime(&time[1]);
    ierr = pTatinLogBasicCPUtime(user,"Checkpoint.write()",time[1]-time[0]);CHKERRQ(ierr);

    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/initial_condition/ptatin.options",checkpoints_path);CHKERRQ(ierr);
    ierr = pTatinWriteOptionsFile(filename);CHKERRQ(ierr);
  }

  /* write out a default string for restarting the job */
  {
    char        restartfile[PETSC_MAX_PATH_LEN];
    char        restartstring[PETSC_MAX_PATH_LEN];
    PetscMPIInt rank;

    ierr = PetscSNPrintf(restartfile,PETSC_MAX_PATH_LEN-1,"%s/restart.default",user->outputpath);CHKERRQ(ierr);
    ierr = PetscSNPrintf(restartstring,PETSC_MAX_PATH_LEN-1,"-restart_directory %s/checkpoints/initial_condition",user->outputpath);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    if (rank == 0) {
      FILE *fp;

      fp = fopen(restartfile,"w");
      fprintf(fp,"%s",restartstring);
      fclose(fp);
    }
  }

  ierr = VecDestroy(&X_e);CHKERRQ(ierr);
  ierr = VecDestroy(&X_s);CHKERRQ(ierr);

  *pctx = user;
  PetscFunctionReturn(0);
}

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
  PetscLogDouble  time[2];
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscTime(&time[0]);
  ierr = pTatin3dLoadContext_FromFile(&user);CHKERRQ(ierr);
  PetscTime(&time[1]);
  ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
  ierr = pTatinLogNote(user,"  [ptatin_driver.Load]");CHKERRQ(ierr);
  ierr = pTatinLogBasicCPUtime(user,"Checkpoint.read()",time[1]-time[0]);CHKERRQ(ierr);

  /* driver specific options parsed here */

  /* Register all models */
  ierr = pTatinModelLoad(user);CHKERRQ(ierr);
  ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);

  ierr = pTatinModel_Initialize(model,user);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-activate_energy",&activate_energy,NULL);CHKERRQ(ierr);

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
    ierr = PetscSNPrintf(output_path_ic,PETSC_MAX_PATH_LEN-1,"%s/fromfile",output_path);CHKERRQ(ierr);
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

    ierr = PetscSNPrintf(checkpoint_path,PETSC_MAX_PATH_LEN-1,"%s/initial_condition",checkpoints_path);CHKERRQ(ierr);
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

PetscErrorCode DummyRun(pTatinCtx pctx,Vec v1,Vec v2)
{
  PetscInt         k,step0;
  Vec              X_s = NULL,X_e = NULL, velocity,pressure;
  Vec              F_s = NULL,F_e = NULL;
  PetscBool        energy_activated;
  PhysCompStokes   stokes = NULL;
  PhysCompEnergy   energy = NULL;
  DM               dmstokes,dmv,dmp,dmenergy = NULL;
  PetscErrorCode   ierr;
  pTatinModel      model = NULL;
  PetscReal        surface_displacement_max = 1.0e32;
  PetscLogDouble   time[4];
  AuuMultiLevelCtx mgctx;
  PetscReal        timestep,dt_factor = 10.0;
  Mat              J_e = NULL;
  PetscMPIInt      rank;

  PetscFunctionBegin;
  ierr = pTatinLogNote(pctx,"  [ptatin_driver.Execute]");CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  /* driver specific options parsed here */
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt_max_surface_displacement",&surface_displacement_max,NULL);CHKERRQ(ierr);

  ierr = pTatinGetModel(pctx,&model);CHKERRQ(ierr);
  ierr = pTatinGetStokesContext(pctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,&dmp);CHKERRQ(ierr);

  if (v1) {
    X_s = v1;
  } else {
    ierr = DMCreateGlobalVector(dmstokes,&X_s);CHKERRQ(ierr);
  }
  ierr = VecDuplicate(X_s,&F_s);CHKERRQ(ierr);

  ierr = pTatinContextValid_Energy(pctx,&energy_activated);CHKERRQ(ierr);
  if (energy_activated) {
    ierr = pTatinGetContext_Energy(pctx,&energy);CHKERRQ(ierr);
    dmenergy = energy->daT;
    ierr = pTatinLogBasicDMDA(pctx,"energy_dmda",dmenergy);CHKERRQ(ierr);
    if (v2) {
      X_e = v2;
    } else {
      ierr = DMCreateGlobalVector(dmenergy,&X_e);CHKERRQ(ierr);
      ierr = pTatinPhysCompAttachData_Energy(pctx,X_e,NULL);CHKERRQ(ierr);
    }
    ierr = VecDuplicate(X_e,&F_e);CHKERRQ(ierr);
    ierr = DMSetMatType(dmenergy,MATAIJ);CHKERRQ(ierr);
    ierr = DMCreateMatrix(dmenergy,&J_e);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J_e);CHKERRQ(ierr);
  }

  ierr = HMG_SetUp(&mgctx,pctx);CHKERRQ(ierr);

  step0 = pctx->step + 1;
  for (k=step0; k<=pctx->nsteps; k++) {
    PetscReal dt_curr,dt_next;

    /* use last saved time step for any calculations, e.g. update particles */
    dt_curr = pctx->dt;

    PetscPrintf(PETSC_COMM_WORLD,"<<----------------------------------------------------------------------------------------------->>\n");
    PetscPrintf(PETSC_COMM_WORLD,"  [[ Executing time step %D ]] : current time %+1.4e : advancing time by dt %+1.4e \n",k,pctx->time,dt_curr);

    pctx->step = k;
    pctx->time += dt_curr;
    if (energy_activated) {
      energy->time = pctx->time;
      energy->dt = pctx->dt;
    }

    /* <<<------------------------------------------------------- >>> */
    /* <<<------------------------------------------------------- >>> */
    /* <<<------------------------------------------------------- >>> */

    ierr = pTatinLogBasic(pctx);CHKERRQ(ierr);

    /* update marker time dependent terms */
    /* e.g. e_plastic^1 = e_plastic^0 + dt * [ strain_rate_inv(u^0) ] */
    /*
     NOTE: for a consistent forward difference time integration we evaluate u^0 at x^0
     - thus this update is performed BEFORE we advect the markers
     */
    ierr = pTatin_UpdateCoefficientTemporalDependence_Stokes(pctx,X_s);CHKERRQ(ierr);

    /* update marker positions */
    ierr = DMCompositeGetAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);
    ierr = MaterialPointStd_UpdateGlobalCoordinates(pctx->materialpoint_db,dmv,velocity,pctx->dt);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);

    /* update mesh */
    ierr = pTatinModel_UpdateMeshGeometry(model,pctx,X_s);CHKERRQ(ierr);

    /* update mesh coordinate hierarchy */
    ierr = DMDARestrictCoordinatesHierarchy(mgctx.dav_hierarchy,mgctx.nlevels);CHKERRQ(ierr);

    /* 3 Update local coordinates and communicate */
    ierr = MaterialPointStd_UpdateCoordinates(pctx->materialpoint_db,dmv,pctx->materialpoint_ex);CHKERRQ(ierr);

    /* 3a - Add material */
    ierr = pTatinModel_ApplyMaterialBoundaryCondition(model,pctx);CHKERRQ(ierr);

    /* add / remove points if cells are over populated or depleted of points */
    ierr = MaterialPointPopulationControl_v1(pctx);CHKERRQ(ierr);

    /* update markers = >> gauss points */
    {
      int               npoints;
      DataField         PField_std;
      DataField         PField_stokes;
      MPntStd           *mp_std;
      MPntPStokes       *mp_stokes;

      DataBucketGetSizes(pctx->materialpoint_db,&npoints,NULL,NULL);
      DataBucketGetDataFieldByName(pctx->materialpoint_db, MPntStd_classname     , &PField_std);
      DataBucketGetDataFieldByName(pctx->materialpoint_db, MPntPStokes_classname , &PField_stokes);
      DataFieldGetEntries(PField_std,(void**)&mp_std);
      DataFieldGetEntries(PField_stokes,(void**)&mp_stokes);

      ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(pctx->coefficient_projection_type,npoints,mp_std,mp_stokes,mgctx.nlevels,mgctx.interpolation_eta,mgctx.dav_hierarchy,mgctx.volQ);CHKERRQ(ierr);

      DataFieldRestoreEntries(PField_stokes,(void**)&mp_stokes);
      DataFieldRestoreEntries(PField_std,(void**)&mp_std);
    }
    if (energy_activated) {
      /* copy current (undeformed) energy mesh coords, update energy mesh geometry */
      ierr = pTatinPhysCompEnergy_Update(energy,dmv,X_e);CHKERRQ(ierr);

      /* update v-V using new mesh coords and the previous mesh coords */
      ierr = pTatinPhysCompEnergy_UpdateALEVelocity(stokes,X_s,energy,energy->dt);CHKERRQ(ierr);

      /* update marker props on new mesh configuration */
      ierr = pTatinPhysCompEnergy_MPProjectionQ1(pctx);CHKERRQ(ierr);
    }

    /* Update boundary conditions */
    /* Fine level setup */
    ierr = pTatinModel_ApplyBoundaryCondition(model,pctx);CHKERRQ(ierr);
    /* Coarse grid setup: Configure boundary conditions */
    ierr = pTatinModel_ApplyBoundaryConditionMG(mgctx.nlevels,mgctx.u_bclist,mgctx.dav_hierarchy,model,pctx);CHKERRQ(ierr);

    /* solve energy + mechanics */
    /* (i) solve energy equation */
    if (energy_activated) {
      SNES snesT;

      ierr = SNESCreate(PETSC_COMM_WORLD,&snesT);CHKERRQ(ierr);
      ierr = SNESSetOptionsPrefix(snesT,"T_");CHKERRQ(ierr);
      ierr = SNESSetFunction(snesT,F_e,  SNES_FormFunctionEnergy,(void*)pctx);CHKERRQ(ierr);
      ierr = SNESSetJacobian(snesT,J_e,J_e,SNES_FormJacobianEnergy,(void*)pctx);CHKERRQ(ierr);
      ierr = SNESSetType(snesT,SNESKSPONLY);
      ierr = SNESSetFromOptions(snesT);CHKERRQ(ierr);

      PetscPrintf(PETSC_COMM_WORLD,"   [[ COMPUTING THERMAL FIELD FOR STEP : %D ]]\n",k);

      /* insert boundary conditions into solution vector */
      ierr = BCListInsert(energy->T_bclist,X_e);CHKERRQ(ierr);

      PetscTime(&time[0]);
      ierr = SNESSolve(snesT,NULL,X_e);CHKERRQ(ierr);
      PetscTime(&time[1]);
      ierr = pTatinLogBasicSNES(pctx,"Energy",snesT);CHKERRQ(ierr);
      ierr = pTatinLogBasicCPUtime(pctx,"Energy",time[1]-time[0]);CHKERRQ(ierr);
      ierr = SNESDestroy(&snesT);CHKERRQ(ierr);
    }

    /* (ii) solve stokes */
    {
      SNES snes;
      Mat  A,B;

      ierr = HMGOperator_SetUp(&mgctx,pctx,&A,&B);CHKERRQ(ierr);
      ierr = pTatinNonlinearStokesSolveCreate(pctx,A,B,F_s,&mgctx,&snes);CHKERRQ(ierr);

      PetscPrintf(PETSC_COMM_WORLD,"   [[ COMPUTING FLOW FIELD FOR STEP : %D ]]\n",k);
      ierr = pTatinNonlinearStokesSolve(pctx,snes,X_s,NULL);CHKERRQ(ierr);

      ierr = MatDestroy(&A);CHKERRQ(ierr);
      ierr = MatDestroy(&B);CHKERRQ(ierr);
      ierr = SNESDestroyMGCtx(snes);CHKERRQ(ierr);
      ierr = SNESDestroy(&snes);CHKERRQ(ierr);
      ierr = HMGOperator_Destroy(&mgctx);CHKERRQ(ierr);
    }

    /* compute new time step */
    dt_next = 0.0;
    //dt_next = 200.0 * rand()/((PetscReal)RAND_MAX);
    dt_next = 0.001245 * ((PetscReal)k) + 1.0e-4;

    /* compute timestep */
    pctx->dt = 1.0e32;
    ierr = DMCompositeGetAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);
    ierr = SwarmUpdatePosition_ComputeCourantStep(dmv,velocity,&timestep);CHKERRQ(ierr);
    timestep = timestep/dt_factor;
    ierr = pTatin_SetTimestep(pctx,"StkCourant",timestep);CHKERRQ(ierr);
    dt_next = pctx->dt;

    ierr = UpdateMeshGeometry_ComputeSurfaceCourantTimestep(dmv,velocity,surface_displacement_max,&timestep);CHKERRQ(ierr);
    ierr = pTatin_SetTimestep(pctx,"StkSurfaceCourant",timestep);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(dmstokes,X_s,&velocity,&pressure);CHKERRQ(ierr);
    dt_next = pctx->dt;

    //PetscPrintf(PETSC_COMM_WORLD,"  timestep_stokes[%D] dt_courant = %1.4e \n", k,pctx->dt );
    if (energy_activated) {
      ierr = pTatinPhysCompEnergy_ComputeTimestep(energy,X_e,&timestep);CHKERRQ(ierr);
      ierr = pTatin_SetTimestep(pctx,"AdvDiffCourant",timestep);CHKERRQ(ierr);
      dt_next = pctx->dt;
      //PetscPrintf(PETSC_COMM_WORLD,"  timestep_advdiff[%D] dt_courant = %1.4e \n", k,pctx->dt );
    }

    /* <<<------------------------------------------------------- >>> */
    /* <<<------------------------------------------------------- >>> */
    /* <<<------------------------------------------------------- >>> */

    // probably dont need this
    if (dt_next > pctx->dt_max) { dt_next = pctx->dt_max; }
    if (dt_next < pctx->dt_min) { dt_next = pctx->dt_min; }

    pctx->dt = dt_next;

    if (k%pctx->output_frequency == 0){
      char prefix[PETSC_MAX_PATH_LEN];

      ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step%D",k);CHKERRQ(ierr);
      PetscTime(&time[0]);
      ierr = pTatinModel_Output(model,pctx,X_s,prefix);CHKERRQ(ierr);
      PetscTime(&time[1]);
    }

    ierr = pTatin3dCheckpointManager(pctx,X_s);CHKERRQ(ierr);

    if (pctx->time > pctx->time_max) break;
  }

  ierr = MatDestroy(&J_e);CHKERRQ(ierr);
  ierr = VecDestroy(&F_s);CHKERRQ(ierr);
  ierr = VecDestroy(&F_e);CHKERRQ(ierr);
  ierr = HMG_Destroy(&mgctx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  PetscErrorCode ierr;
  PetscBool      init = PETSC_FALSE,load = PETSC_FALSE, run = PETSC_FALSE;
  pTatinCtx      pctx = NULL;

  ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
  ierr = pTatinModelRegisterAll();CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,NULL,"-init",&init,NULL);CHKERRQ(ierr);
  if (init) {
    ierr = GenerateICStateFromModelDefinition(&pctx);CHKERRQ(ierr);
    if (pctx) { ierr = pTatin3dDestroyContext(&pctx); }
    pctx = NULL;
  }

  ierr = PetscOptionsGetBool(NULL,NULL,"-load",&load,NULL);CHKERRQ(ierr);
  if (load) {
    ierr = LoadICStateFromModelDefinition(&pctx,NULL,NULL,PETSC_TRUE);CHKERRQ(ierr);
    /* do something */
    if (pctx) { ierr = pTatin3dDestroyContext(&pctx); }
    pctx = NULL;
  }

  ierr = PetscOptionsGetBool(NULL,NULL,"-run",&run,NULL);CHKERRQ(ierr);
  if (run || (!init && !load)) {
    Vec       Xup,Xt;
    PetscBool restart_string_found = PETSC_FALSE,flg = PETSC_FALSE;
    char      outputpath[PETSC_MAX_PATH_LEN];

    /* look for a default restart file */
    ierr = PetscOptionsGetString(NULL,NULL,"-output_path",outputpath,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (flg) {
      char fname[PETSC_MAX_PATH_LEN];

      ierr = PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"%s/restart.default",outputpath);CHKERRQ(ierr);
      ierr = pTatinTestFile(fname,'r',&restart_string_found);CHKERRQ(ierr);
      if (restart_string_found) {
        PetscPrintf(PETSC_COMM_WORLD,"[pTatin] Detected default restart option file helper: %s\n",fname);
        //ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD,NULL,fname,PETSC_TRUE);CHKERRQ(ierr);
        ierr = PetscOptionsInsert(NULL,&argc,&argv,fname);CHKERRQ(ierr);
      }
    }

    ierr = LoadICStateFromModelDefinition(&pctx,&Xup,&Xt,PETSC_FALSE);CHKERRQ(ierr);

    ierr = DummyRun(pctx,Xup,Xt);CHKERRQ(ierr);

    ierr = VecDestroy(&Xup);CHKERRQ(ierr);
    ierr = VecDestroy(&Xt);CHKERRQ(ierr);
  }

  if (pctx) { ierr = pTatin3dDestroyContext(&pctx); }

  ierr = pTatinModelDeRegisterAll();CHKERRQ(ierr);
  ierr = pTatinFinalize();CHKERRQ(ierr);
  return(0);
}

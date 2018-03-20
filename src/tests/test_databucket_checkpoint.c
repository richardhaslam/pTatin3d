
#include <ptatin3d.h>
#include <private/ptatin_impl.h>
#include <ptatin_init.h>
#include <ptatin_models.h>
#include <data_bucket.h>
#include <dmda_checkpoint.h>
#include <MPntStd_def.h>

/*
 stats[0,1,2] = min-x,y,z
 stats[3,4,5] = max-x,y,z
 stats[6,7,8] = average x,y,z
*/
PetscErrorCode MaterialPointComputeCoordStats(DataBucket db,PetscReal stats[])
{
  int d,p,npoints;
  long int npoints_,npoints_g_;
  DataField PField_std;
  double stats_[9],stats_g_[9];
  PetscErrorCode ierr;
  
  for (p=0; p<3; p++) {
    stats_[p] = 1.0e32;
  }
  for (p=3; p<6; p++) {
    stats_[p] = -1.0e32;
  }
  for (p=6; p<9; p++) {
    stats_[p] = 0.0;
  }
  
  DataBucketGetSizes(db,&npoints,0,0);
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);
  for (p=0; p<npoints; p++) {
    MPntStd *material_point;
    double  *position;
    
    DataFieldAccessPoint(PField_std,p,(void**)&material_point);
    MPntStdGetField_global_coord(material_point,&position);
    
    for (d=0; d<3; d++) {
      stats_[0+d] = PetscMin(stats_[0+d],position[d]);
    }

    for (d=0; d<3; d++) {
      stats_[3+d] = PetscMax(stats_[3+d],position[d]);
    }
    for (d=0; d<3; d++) {
      stats_[6+d] += position[d];
    }
  }
  DataFieldRestoreAccess(PField_std);

  ierr = MPI_Allreduce(&stats_[0],&stats_g_[0],3, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&stats_[3],&stats_g_[3],3, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&stats_[6],&stats_g_[6],3, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);

  npoints_ = (long int)npoints;
  ierr = MPI_Allreduce(&npoints_,&npoints_g_,1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  for (d=0; d<3; d++) {
    stats_g_[6+d] = stats_g_[6+d] / ((double)npoints_g_);
  }
  
  for (d=0; d<9; d++) {
    stats[d] = (PetscReal)stats_g_[d];
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode ptatin_db_checkpoint(void)
{
  PetscErrorCode ierr;
  DM             dav;
  DataBucket     material_points;
  pTatinCtx      user;
  pTatinModel    model;
  PetscLogDouble time[4];
  PetscReal      statsOrig[9],statsCP[9];
  PetscInt       k;
  char           filename[PETSC_MAX_PATH_LEN];
  
  PetscFunctionBegin;
  
  ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
  ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
  
  /* Register all models */
  ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
  /* Load model, call an initialization routines */
  ierr = pTatinModelLoad(user);CHKERRQ(ierr);
  
  ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);
  ierr = pTatinModel_Initialize(model,user);CHKERRQ(ierr);
  
  /* Generate physics modules */
  ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
  
  /* Pack all physics together */
  /* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
  ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
  user->pack = user->stokes_ctx->stokes_pack;
  
  /* fetch some local variables */
  dav = user->stokes_ctx->dav;
  
  ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
  
  /* mesh geometry */
  ierr = pTatinModel_ApplyInitialMeshGeometry(model,user);CHKERRQ(ierr);
  
  /* interpolate material point coordinates (needed if mesh was modified) */
  ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
  
  /* define material region/parameters */
  ierr = pTatinModel_ApplyInitialMaterialGeometry(model,user);CHKERRQ(ierr);
  
  ierr = pTatinGetMaterialPoints(user,&material_points,NULL);CHKERRQ(ierr);
  
  ierr = MaterialPointComputeCoordStats(material_points,statsOrig);CHKERRQ(ierr);
  
  {
    DM dms;
    
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/stokes",user->outputpath);CHKERRQ(ierr);
    ierr = DMDACheckpointWrite(dav,filename);CHKERRQ(ierr);

    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/stokes_dmda.json",user->outputpath);CHKERRQ(ierr);
    ierr = DMDACheckpointLoad(PETSC_COMM_WORLD,filename,&dms);CHKERRQ(ierr);
    ierr = DMView(dms,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    //ierr = DMDACheckpointLoad(PETSC_COMM_SELF,"./cptest/stokes_dmda.json",&dms);CHKERRQ(ierr);
    //ierr = DMView(dms,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    ierr = DMDestroy(&dms);CHKERRQ(ierr);
  }
  
  // Write databucket
  //DataBucketView(PETSC_COMM_WORLD,material_points,"mp-ascii",DATABUCKET_VIEW_STDOUT);
  //DataBucketView(PETSC_COMM_WORLD,material_points,"mp-binary",DATABUCKET_VIEW_BINARY);
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/mp-native",user->outputpath);CHKERRQ(ierr);
  PetscTime(&time[0]);
  DataBucketView(PETSC_COMM_WORLD,material_points,filename,DATABUCKET_VIEW_NATIVE);
  PetscTime(&time[1]);
  
  {
    DataBucket dbload;
    PetscMPIInt rank;
    
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/mp-native_db.json",user->outputpath);CHKERRQ(ierr);
    PetscTime(&time[2]);
    DataBucketLoad_NATIVE(PETSC_COMM_WORLD,filename,&dbload);
    PetscTime(&time[3]);
    
    DataBucketView(PETSC_COMM_WORLD,dbload,"mp-ascii",DATABUCKET_VIEW_STDOUT);

    // Better way to load redundant data - relies on the network rather than the file system
    //DataBucketLoadRedundant_NATIVE(PETSC_COMM_WORLD,"./cptest/mp-native_db.json",&dbload);
    //DataBucketView(PETSC_COMM_SELF,dbload,"mp-ascii[red]",DATABUCKET_VIEW_STDOUT);
    
    // Code below makes sense in two scenarios:
    //   (i) if you run once with N ranks then write; launch again on 1 rank and ONLY read
    //   (ii) if you are fine having every rank open and read the file
    //DataBucketLoad_NATIVE(PETSC_COMM_SELF,"cptest/mp-native_db.json",&dbload);
    //DataBucketView(PETSC_COMM_SELF,dbload,"mp-ascii",DATABUCKET_VIEW_STDOUT);

    /*
    // as a sanity check - write out the loaded data bucket data into a PV file
    {
      const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
      
      ierr = SwarmViewGeneric_ParaView(dbload,3,mp_prop_list,"./cptest","cpmprop");CHKERRQ(ierr);
    //ierr = SwarmOutputParaView_MPntStd(dbload,"./","mpstd");CHKERRQ(ierr);
    }
    */
     
    /*
    if (rank == 0) {
      const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
      
      ierr = SwarmViewGeneric_ParaView(dbload,3,mp_prop_list,"./","seqmp");CHKERRQ(ierr);

    //ierr = SwarmOutputParaView_MPntStd(dbload,"./","seqmp");CHKERRQ(ierr);
    }
    */
    
    ierr = MaterialPointComputeCoordStats(dbload,statsCP);CHKERRQ(ierr);

    DataBucketDestroy(&dbload);
  }
  
  PetscPrintf(PETSC_COMM_WORLD,"DataBucket checkpoint write %+1.4e (sec)\n",time[1]-time[0]);
  PetscPrintf(PETSC_COMM_WORLD,"DataBucket checkpoint load  %+1.4e (sec)\n",time[3]-time[2]);
  
  PetscPrintf(PETSC_COMM_WORLD,"  [Orig] x-coor range %+1.12e , %+1.12e\n",statsOrig[0],statsOrig[3]);
  PetscPrintf(PETSC_COMM_WORLD,"  [Orig] y-coor range %+1.12e , %+1.12e\n",statsOrig[1],statsOrig[4]);
  PetscPrintf(PETSC_COMM_WORLD,"  [Orig] z-coor range %+1.12e , %+1.12e\n",statsOrig[2],statsOrig[5]);
  PetscPrintf(PETSC_COMM_WORLD,"  [Orig] coor avg     %+1.12e , %+1.12e , %+1.12e\n",statsOrig[6],statsOrig[7],statsOrig[8]);

  PetscPrintf(PETSC_COMM_WORLD,"  [  CP] x-coor range %+1.12e , %+1.12e\n",statsCP[0],statsCP[3]);
  PetscPrintf(PETSC_COMM_WORLD,"  [  CP] y-coor range %+1.12e , %+1.12e\n",statsCP[1],statsCP[4]);
  PetscPrintf(PETSC_COMM_WORLD,"  [  CP] z-coor range %+1.12e , %+1.12e\n",statsCP[2],statsCP[5]);
  PetscPrintf(PETSC_COMM_WORLD,"  [  CP] coor avg     %+1.12e , %+1.12e , %+1.12e\n",statsCP[6],statsCP[7],statsCP[8]);

  for (k=0; k<9; k++) {
    PetscReal diff;
    diff = PetscAbsReal(statsOrig[k]-statsCP[k]);
    if (diff < 1.0e-8) {
      PetscPrintf(PETSC_COMM_WORLD,"  coord statistic %D PASSED\n",k);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,"  coord statistic %D FAILED\n",k);
    }
  }
  
  /* define boundary conditions */
  ierr = pTatinModel_ApplyBoundaryCondition(model,user);CHKERRQ(ierr);
  
  /* write out vts file */
  //ierr = pTatinModel_Output(model,user,X,"test");CHKERRQ(ierr);
  
  ierr = pTatin3dDestroyContext(&user);
  
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  
  ierr = pTatinInitialize(&argc,&argv,0,NULL);CHKERRQ(ierr);
  
  ierr = ptatin_db_checkpoint();CHKERRQ(ierr);
  
  ierr = pTatinFinalize();CHKERRQ(ierr);
  return 0;
}

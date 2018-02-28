
#include <ptatin3d.h>
#include <private/ptatin_impl.h>
#include <ptatin_init.h>
#include <ptatin_models.h>
#include <data_bucket.h>
#include <dmda_checkpoint.h>

#undef __FUNCT__
#define __FUNCT__ "ptatin_db_checkpoint"
PetscErrorCode ptatin_db_checkpoint(void)
{
  PetscErrorCode ierr;
  DM             multipys_pack,dav;
  DataBucket     material_points;
  pTatinCtx      user;
  Vec            X;
  pTatinModel    model;
  
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
  multipys_pack = user->pack;
  dav           = user->stokes_ctx->dav;
  
  ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
  
  ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
  
  /* mesh geometry */
  ierr = pTatinModel_ApplyInitialMeshGeometry(model,user);CHKERRQ(ierr);
  
  /* interpolate material point coordinates (needed if mesh was modified) */
  ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
  
  /* define material region/parameters */
  ierr = pTatinModel_ApplyInitialMaterialGeometry(model,user);CHKERRQ(ierr);
  
  ierr = pTatinGetMaterialPoints(user,&material_points,NULL);CHKERRQ(ierr);
  
  {
    DM dms;
    
    ierr = DMDACheckpointWrite(dav,"./cpdump/stokes");CHKERRQ(ierr);

    ierr = DMDACheckpointLoad(PETSC_COMM_WORLD,"./cpdump/stokes_dmda.json",&dms);CHKERRQ(ierr);
    ierr = DMView(dms,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    //ierr = DMDACheckpointLoad(PETSC_COMM_SELF,"./cpdump/stokes_dmda.json",&dms);CHKERRQ(ierr);
    //ierr = DMView(dms,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  }
  
  // Write to file
  //	ierr = MaterialPointDataBasicLoadIntoListFromFile(user->materialpoint_db,dav,PETSC_FALSE,"filters/coords_markers.dat","filters/phase_markers.dat");CHKERRQ(ierr);
  //DataBucketView(PETSC_COMM_WORLD,material_points,"mp-ascii",DATABUCKET_VIEW_STDOUT);
  //DataBucketView(PETSC_COMM_WORLD,material_points,"mp-binary",DATABUCKET_VIEW_BINARY);
  
  //DataBucketView(PETSC_COMM_WORLD,material_points,"cpdump/mp-native",DATABUCKET_VIEW_NATIVE);
  
  {
    DataBucket dbload;
    int rank;
    
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
    DataBucketLoad_NATIVE(PETSC_COMM_WORLD,"cpdump/mp-native_db.json",&dbload);
    DataBucketView(PETSC_COMM_SELF,dbload,"mp-ascii",DATABUCKET_VIEW_STDOUT);

    // Better way to load redundant data - relies on the network rather than the file system
    //DataBucketLoadRedundant_NATIVE(PETSC_COMM_WORLD,"cpdump/mp-native_db.json",&dbload);
    //DataBucketView(PETSC_COMM_SELF,dbload,"mp-ascii[red]",DATABUCKET_VIEW_STDOUT);
    
    // Code below makes sense in two scenarios:
    //   (i) if you run once with N ranks then write; launch again on 1 rank and ONLY read
    //   (ii) if you are fine having every rank open and read the file
    //DataBucketLoad_NATIVE(PETSC_COMM_SELF,"cpdump/mp-native_db.json",&dbload);
    //DataBucketView(PETSC_COMM_SELF,dbload,"mp-ascii",DATABUCKET_VIEW_STDOUT);

    {
      const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
      
      ierr = SwarmViewGeneric_ParaView(dbload,3,mp_prop_list,"./","allmprop");CHKERRQ(ierr);
    //ierr = SwarmOutputParaView_MPntStd(dbload,"./","mpstd");CHKERRQ(ierr);
    }

    /*
    if (rank == 0) {
      const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl };
      
      ierr = SwarmViewGeneric_ParaView(dbload,3,mp_prop_list,"./","seqmp");CHKERRQ(ierr);

    //ierr = SwarmOutputParaView_MPntStd(dbload,"./","seqmp");CHKERRQ(ierr);
    }
    */
  }
  
  
  /* define boundary conditions */
  ierr = pTatinModel_ApplyBoundaryCondition(model,user);CHKERRQ(ierr);
  
  /* write out vts file */
  ierr = pTatinModel_Output(model,user,X,"test");CHKERRQ(ierr);
  
  ierr = pTatin3dDestroyContext(&user);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  
  ierr = pTatinInitialize(&argc,&argv,0,NULL);CHKERRQ(ierr);
  
  ierr = ptatin_db_checkpoint();CHKERRQ(ierr);
  
  ierr = pTatinFinalize();CHKERRQ(ierr);
  return 0;
}

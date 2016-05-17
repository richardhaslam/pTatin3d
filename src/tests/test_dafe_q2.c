

#include <ptatin3d.h>
#include <ptatin_init.h>
#include <dafe.h>
#include <mesh_deformation.h>

#undef __FUNCT__
#define __FUNCT__ "test1_create"
PetscErrorCode test1_create(void)
{
  PetscErrorCode ierr;
  DM          dm;
  PetscViewer viewer;
  
  ierr = DMDAFECreate3d(PETSC_COMM_WORLD,3,DAFE_Q2,12,14,16,&dm);CHKERRQ(ierr);
  
  ierr = DMDAFESetUniformCoordinates(dm,0.0,1.0,-0.5,0.5,0.2,2.2);CHKERRQ(ierr);
  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERVTK);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,"dafe.vts");CHKERRQ(ierr);
  ierr = DMView(dm,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test2_createsubP1"
PetscErrorCode test2_createsubP1(void)
{
  PetscErrorCode ierr;
  DM          dm,dmsub;
  PetscInt    nref;
  
  ierr = DMDAFECreate3d(PETSC_COMM_WORLD,3,DAFE_Q2,12,14,16,&dm);CHKERRQ(ierr);
  ierr = DMDAFESetUniformCoordinates(dm,0.0,1.0,-0.5,0.5,0.2,2.2);CHKERRQ(ierr);
  
  nref = 2;
  PetscOptionsGetInt(NULL,NULL,"-ref",&nref,NULL);
  ierr = DMDAFECreateSubDMDAFE(dm,1,nref,DAFE_P1sub,&dmsub);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"---- Original DAFE ----\n");
  ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"---- Sub-divided DAFE ----\n");
  ierr = DMView(dmsub,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = DMDestroy(&dmsub);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test3_createsubQ1"
PetscErrorCode test3_createsubQ1(void)
{
  PetscErrorCode ierr;
  DM          dm,dmsub;
  PetscInt    nref,ndof;
  PetscViewer viewer;
  
  ierr = DMDAFECreate3d(PETSC_COMM_WORLD,3,DAFE_Q2,12,14,16,&dm);CHKERRQ(ierr);
  ierr = DMDAFESetUniformCoordinates(dm,0.0,1.0,-0.5,0.5,0.2,2.2);CHKERRQ(ierr);
  
  ndof = 1;
  nref = 1;
  PetscOptionsGetInt(NULL,NULL,"-ref",&nref,NULL);
  ierr = DMDAFECreateSubDMDAFE(dm,ndof,nref,DAFE_Q1sub,&dmsub);CHKERRQ(ierr);
  ierr = DMDAFESetUniformCoordinates(dmsub,0.0,1.0,-0.5,0.5,0.2,2.2);CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"---- Original DAFE ----\n");
  ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"---- Sub-divided DAFE ----\n");
  ierr = DMView(dmsub,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERVTK);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,"dafeQ2.vts");CHKERRQ(ierr);
  ierr = DMView(dm,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERVTK);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,"dafeQ1.vts");CHKERRQ(ierr);
  ierr = DMView(dmsub,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = DMDestroy(&dmsub);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "test4_projectQ1"
PetscErrorCode test4_projectQ1(void)
{
  PetscErrorCode ierr;
  DM          dm,dmsub,dmda;
  PetscInt    nref,ndof,mx;
  PetscViewer viewer;
  PetscLogDouble t0,t1;
  
  mx = 12;
  ierr = PetscOptionsGetInt(NULL,NULL,"-mx",&mx,NULL);CHKERRQ(ierr);
  ierr = DMDAFECreate3d(PETSC_COMM_WORLD,3,DAFE_Q2,mx,mx,mx,&dm);CHKERRQ(ierr);
  ierr = DMDAFESetUniformCoordinates(dm,-1.0,1.0,-1.0,1.0,-1.0,1.0);CHKERRQ(ierr);

  ierr = DMDAFEGetDA(dm,&dmda);CHKERRQ(ierr);
  ierr = MeshDeformation_GaussianBump_YMAX(dmda,-0.3,-5.6);CHKERRQ(ierr);
  
  ndof = 1;
  nref = 1;
  PetscOptionsGetInt(NULL,NULL,"-ref",&nref,NULL);
  ierr = DMDAFECreateSubDMDAFE(dm,ndof,nref,DAFE_Q1sub,&dmsub);CHKERRQ(ierr);

  PetscTime(&t0);
  ierr = DMDAFEProjectCoordinates(dmsub);CHKERRQ(ierr);
  PetscTime(&t1);
  PetscPrintf(PETSC_COMM_WORLD,"Project coordinates: %1.4e (sec)\n",t1-t0);
  
  PetscPrintf(PETSC_COMM_WORLD,"---- Original DAFE ----\n");
  ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"---- Sub-divided DAFE ----\n");
  ierr = DMView(dmsub,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERVTK);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,"dafeQ2.vts");CHKERRQ(ierr);
  ierr = DMView(dm,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERVTK);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,"dafeQ1.vts");CHKERRQ(ierr);
  ierr = DMView(dmsub,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = DMDestroy(&dmsub);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
  
	ierr = pTatinInitialize(&argc,&argv,0,NULL);CHKERRQ(ierr);
  
  //ierr = test1_create();CHKERRQ(ierr);
  //ierr = test2_createsubP1();CHKERRQ(ierr);
  //ierr = test3_createsubQ1();CHKERRQ(ierr);
  ierr = test4_projectQ1();CHKERRQ(ierr);
  
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

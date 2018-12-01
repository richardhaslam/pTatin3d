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
 **    filename:   ptatin3d.c
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

#include "petsc/private/dmdaimpl.h"
#include "petsc.h"

#include "ptatin3d.h"
#include "ptatin3d_defs.h"

#include "element_type_Q2.h"
#include "dmda_element_q2p1.h"
#include "dmda_view_petscvtk.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPStokesPl_def.h"
#include "material_point_std_utils.h"
#include "ptatin_utils.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "phys_comp_energy.h"
#include "output_paraview.h"
#include "ptatin_models.h"
#include "dmda_checkpoint.h"
#include "ptatin_log.h"
#include "dmda_project_coords.h"
#include "cJSON.h"
#include "cjson_utils.h"

#include "private/ptatin_impl.h"


/* PHYSICS MESHES */

PetscErrorCode pTatin3d_PhysCompStokesNew(pTatinCtx user)
{
  PetscErrorCode  ierr;
  PhysCompStokes  stokes;

  PetscFunctionBegin;
  ierr = PhysCompCreate_Stokes(&stokes);CHKERRQ(ierr);

  stokes->mx = user->mx;
  stokes->my = user->my;
  stokes->mz = user->mz;
  stokes->use_mf_stokes = user->use_mf_stokes;

  ierr = PhysCompCreateMesh_Stokes3d(stokes->mx,stokes->my,stokes->mz,stokes);CHKERRQ(ierr);
  ierr = PhysCompCreateBoundaryList_Stokes(stokes);CHKERRQ(ierr);
  ierr = PhysCompCreateVolumeQuadrature_Stokes(stokes);CHKERRQ(ierr);
  ierr = PhysCompCreateSurfaceQuadrature_Stokes(stokes);CHKERRQ(ierr);

  user->stokes_ctx = stokes;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_PhysCompStokesCreate(pTatinCtx user)
{
  PetscErrorCode  ierr;
  PhysCompStokes  stokes;
  PetscReal       grav[3];
  PetscInt        ncomponents;
  PetscBool       flg;

  PetscFunctionBegin;

  ierr = pTatin3d_PhysCompStokesNew(user);CHKERRQ(ierr);

  /* Default action - set gravity vector to be 0,1,0 to not break existing models which set -rho.g on the material points */
  /* Model initialize function can overload the gravity value */
  ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
  grav[0] = 0.0;
  grav[1] = 1.0;
  grav[2] = 0.0;
  ncomponents = 3;
  flg = PETSC_FALSE;
  PetscOptionsGetRealArray(NULL,NULL,"-stokes_gravity_vector",grav,&ncomponents,&flg);
  ierr = PhysCompStokesSetGravityVector(stokes,grav);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode default_pTatin3d_ModelOutput_VelocityPressure_Stokes(pTatinCtx ctx,Vec X,const char prefix[])
{
  PetscErrorCode ierr;
  char           *name;
  DM             stokes_pack;
  Vec            UP;
  PetscLogDouble t0,t1;
  static PetscBool beenhere=PETSC_FALSE;
  PetscFunctionBegin;

  PetscTime(&t0);
  // PVD
  {
    char *pvdfilename;
    char *vtkfilename;

    if (asprintf(&pvdfilename,"%s/timeseries_vp.pvd",ctx->outputpath) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
    if (prefix) {
      if (asprintf(&vtkfilename, "%s_vp.pvts",prefix) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
    } else {
      if (asprintf(&vtkfilename, "vp.pvts") < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
    }

    if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename ); }
    ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time, vtkfilename, "");CHKERRQ(ierr);
    beenhere = PETSC_TRUE;
    free(pvdfilename);
    free(vtkfilename);
  }

  // PVTS + VTS
  if (prefix) {
    if (asprintf(&name,"%s_vp",prefix) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
  } else {
    if (asprintf(&name,"vp") < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
  }

  //PetscPrintf(PETSC_COMM_WORLD,"[[DESIGN FLAW]] %s: require better physics modularity to extract (u,p) <---| (X) \n", PETSC_FUNCTION_NAME );

  stokes_pack = ctx->stokes_ctx->stokes_pack;
  UP = X;
  ierr = pTatinOutputParaViewMeshVelocityPressure(stokes_pack,UP,ctx->outputpath,name);CHKERRQ(ierr);
  free(name);
  PetscTime(&t1);
  /*PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_vp.(pvd,pvts,vts): CPU time %1.2e (sec) \n", PETSC_FUNCTION_NAME,prefix,t1-t0);*/

  PetscFunctionReturn(0);
}

PetscErrorCode _pTatin3d_ModelOutput_VelocityPressure_Stokes(pTatinCtx ctx,Vec X,const char root[],const char snap[],const char fileprefix[])
{
  PetscErrorCode ierr;
  DM             stokes_pack;
  Vec            UP;
  static PetscBool beenhere=PETSC_FALSE;
  char pvdfilename[PETSC_MAX_PATH_LEN],vtkfilename[PETSC_MAX_PATH_LEN],name[PETSC_MAX_PATH_LEN],pvtsoutputdir[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  if (snap) { PetscSNPrintf(pvtsoutputdir,PETSC_MAX_PATH_LEN-1,"%s/%s",root,snap);
  } else {    PetscSNPrintf(pvtsoutputdir,PETSC_MAX_PATH_LEN-1,"%s",root);         }

  // PVD
  PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_vp.pvd",root);

  if (fileprefix) { PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "%s_vp.pvts",fileprefix);
  } else {          PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "vp.pvts"); }

  if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename ); }
  ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time, vtkfilename, snap);CHKERRQ(ierr);
  beenhere = PETSC_TRUE;

  // PVTS + VTS
  if (fileprefix) { PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1,"%s_vp",fileprefix);
  } else {          PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1,"vp",fileprefix);    }

  stokes_pack = ctx->stokes_ctx->stokes_pack;
  UP = X;
  ierr = pTatinOutputParaViewMeshVelocityPressure(stokes_pack,UP,pvtsoutputdir,name);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_ModelOutput_VelocityPressure_Stokes_per_dir(pTatinCtx ctx,Vec X)
{
  char name[PETSC_MAX_PATH_LEN];
  PetscBool found;
  PetscErrorCode ierr;

  // create a snapshot directory based on step
  ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/step%D",ctx->outputpath,ctx->step);CHKERRQ(ierr);
  ierr = pTatinTestDirectory(name,'w',&found);CHKERRQ(ierr);
  if (!found) {
    ierr = pTatinCreateDirectory(name);CHKERRQ(ierr);
  }
  ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"step%D",ctx->step);CHKERRQ(ierr);
  ierr = _pTatin3d_ModelOutput_VelocityPressure_Stokes(ctx,X,ctx->outputpath,name,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_ModelOutput_VelocityPressure_Stokes_flat(pTatinCtx ctx,Vec X)
{
  char fileprefix[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  // create a file prefix based on step
  ierr = PetscSNPrintf(fileprefix,PETSC_MAX_PATH_LEN-1,"step%1.6D",ctx->step);CHKERRQ(ierr);
  ierr = _pTatin3d_ModelOutput_VelocityPressure_Stokes(ctx,X,ctx->outputpath,NULL,fileprefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_ModelOutput_VelocityPressure_Stokes(pTatinCtx ctx,Vec X,const char prefix[])
{
  char name[PETSC_MAX_PATH_LEN];
  PetscBool found;
  PetscErrorCode ierr;

  // create a snapshot directory based on step
  ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/step%D",ctx->outputpath,ctx->step);CHKERRQ(ierr);
  ierr = pTatinTestDirectory(name,'w',&found);CHKERRQ(ierr);
  if (!found) {
    ierr = pTatinCreateDirectory(name);CHKERRQ(ierr);
  }
  ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"step%D",ctx->step);CHKERRQ(ierr);
  ierr = _pTatin3d_ModelOutput_VelocityPressure_Stokes(ctx,X,ctx->outputpath,name,prefix);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_ModelOutputLite_Velocity_Stokes(pTatinCtx ctx,Vec X,const char prefix[])
{
  PetscErrorCode ierr;
  DM             stokes_pack;
  Vec            UP;
  PetscLogDouble t0,t1;
  static PetscBool beenhere=PETSC_FALSE;
  char name[PETSC_MAX_PATH_LEN],pvdfilename[PETSC_MAX_PATH_LEN],vtkfilename[PETSC_MAX_PATH_LEN],pvoutputdir[PETSC_MAX_PATH_LEN],root[PETSC_MAX_PATH_LEN];
  char stepprefix[PETSC_MAX_PATH_LEN];
  PetscBool found;

  PetscFunctionBegin;
  ierr = PetscSNPrintf(root,PETSC_MAX_PATH_LEN-1,"%s",ctx->outputpath);CHKERRQ(ierr);

  ierr = PetscSNPrintf(pvoutputdir,PETSC_MAX_PATH_LEN-1,"%s/step%D",root,ctx->step);CHKERRQ(ierr);
  ierr = pTatinTestDirectory(pvoutputdir,'w',&found);CHKERRQ(ierr);
  if (!found) { ierr = pTatinCreateDirectory(pvoutputdir);CHKERRQ(ierr); }

  PetscTime(&t0);
  // PVD
  PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_v.pvd",root);
  if (prefix) { PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "%s_v.pvts",prefix);
  } else {      PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "v.pvts");           }

  if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename ); }
  PetscSNPrintf(stepprefix,PETSC_MAX_PATH_LEN-1,"step%D",ctx->step);
  ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time, vtkfilename, stepprefix);CHKERRQ(ierr);
  beenhere = PETSC_TRUE;

  // PVTS + VTS
  if (prefix) { PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1,"%s_v",prefix);
  } else {      PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1,"v",prefix);    }

  stokes_pack = ctx->stokes_ctx->stokes_pack;
  UP = X;
  ierr = pTatinOutputLiteParaViewMeshVelocity(stokes_pack,UP,pvoutputdir,name);CHKERRQ(ierr);

  PetscTime(&t1);
  /*PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_v.(pvd,pvts,vts): CPU time %1.2e (sec) \n", PETSC_FUNCTION_NAME,prefix,t1-t0);*/

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(pTatinCtx ctx,Vec X,const char prefix[])
{
  PetscErrorCode ierr;
  DM             stokes_pack;
  PetscLogDouble t0,t1;
  static PetscBool beenhere=PETSC_FALSE;
  char           f1[PETSC_MAX_PATH_LEN];
  char           f3[PETSC_MAX_PATH_LEN];
  char pvdfilename[PETSC_MAX_PATH_LEN],vtkfilename[PETSC_MAX_PATH_LEN],pvoutputdir[PETSC_MAX_PATH_LEN],root[PETSC_MAX_PATH_LEN];
  char stepprefix[PETSC_MAX_PATH_LEN];
  PetscBool found;

  PetscFunctionBegin;
  PetscTime(&t0);
  // PVD
  ierr = PetscSNPrintf(root,PETSC_MAX_PATH_LEN-1,"%s",ctx->outputpath);CHKERRQ(ierr);

  ierr = PetscSNPrintf(pvoutputdir,PETSC_MAX_PATH_LEN-1,"%s/step%D",root,ctx->step);CHKERRQ(ierr);
  ierr = pTatinTestDirectory(pvoutputdir,'w',&found);CHKERRQ(ierr);
  if (!found) { ierr = pTatinCreateDirectory(pvoutputdir);CHKERRQ(ierr); }

  PetscTime(&t0);
  // PVD
  PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_X.pvd",root);
  if (prefix) { PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "%s_X.pvts",prefix);
  } else {      PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "X.pvts");           }

  if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename ); }
  PetscSNPrintf(stepprefix,PETSC_MAX_PATH_LEN-1,"step%D",ctx->step);
  ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time, vtkfilename, stepprefix);CHKERRQ(ierr);
  beenhere = PETSC_TRUE;

  stokes_pack = ctx->stokes_ctx->stokes_pack;

  /* dump the dmda's */
  /* dav,dap */
  if (prefix) {
    PetscSNPrintf(f1,PETSC_MAX_PATH_LEN-1,"%s/%s.dmda-velocity",pvoutputdir,prefix);
    PetscSNPrintf(f3,PETSC_MAX_PATH_LEN-1,"%s/%s.dmda-pressure",pvoutputdir,prefix);
  } else {
    PetscSNPrintf(f1,PETSC_MAX_PATH_LEN-1,"%s/dmda-velocity",pvoutputdir);
    PetscSNPrintf(f3,PETSC_MAX_PATH_LEN-1,"%s/dmda-pressure",pvoutputdir);
  }

  ierr = DMDACheckpointWrite(ctx->stokes_ctx->dav,f1);CHKERRQ(ierr);
  ierr = DMDACheckpointWrite(ctx->stokes_ctx->dap,f3);CHKERRQ(ierr);

  /* dump the vectors */
  {
    PetscViewer viewer;
    Vec         Xu,Xp;

    ierr = DMCompositeGetAccess(stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);

    if (prefix) { PetscSNPrintf(f1,PETSC_MAX_PATH_LEN-1,"%s/%s.dmda-Xu",pvoutputdir,prefix); }
    else {        PetscSNPrintf(f1,PETSC_MAX_PATH_LEN-1,"%s/dmda-Xu",pvoutputdir);           }
    PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", f1 );
    ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD,f1,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(Xu,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    if (prefix) { PetscSNPrintf(f1,PETSC_MAX_PATH_LEN-1,"%s/%s.dmda-Xp",pvoutputdir,prefix); }
    else {        PetscSNPrintf(f1,PETSC_MAX_PATH_LEN-1,"%s/dmda-Xp",pvoutputdir);           }
    PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", f1 );
    ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD,f1,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(Xp,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(ctx->stokes_ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
  }

  PetscTime(&t1);
  /*PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_{Xu,Xp}: CPU time %1.2e (sec) \n", PETSC_FUNCTION_NAME,prefix,t1-t0);*/

  PetscFunctionReturn(0);
}

/* MATERIAL POINTS */
PetscErrorCode SwarmDMDA3dDataExchangerCreate(DM da,DataEx *_de)
{
  DataEx de;
  const PetscMPIInt *neighborranks;
  PetscMPIInt neighborranks2[27],neighborcount;
  PetscInt i;
  PetscLogDouble t0,t1;
  PetscMPIInt rank;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)da),&rank);CHKERRQ(ierr);
  ierr = DMDAGetNeighbors(da,&neighborranks);CHKERRQ(ierr);

  neighborcount = 0;
  for (i=0; i<27; i++) {
    neighborranks2[i] = -1;
    if ( (neighborranks[i]>=0) && (neighborranks[i]!=rank) ) {
      neighborranks2[neighborcount] = neighborranks[i];
      neighborcount++;
    }
  }

  PetscTime(&t0);
  de = DataExCreate(PetscObjectComm((PetscObject)da),0);
  //  de = DataExCreate(PETSC_COMM_WORLD,0);
  ierr = DataExTopologyInitialize(de);CHKERRQ(ierr);
  for (i=0; i<neighborcount; i++) {
    ierr = DataExTopologyAddNeighbour(de,neighborranks2[i]);CHKERRQ(ierr);
  }
  ierr = DataExTopologyFinalize(de);CHKERRQ(ierr);
  PetscTime(&t1);
  PetscPrintf(de->comm,"[[SwarmDMDA3dDataExchangerCreate: time = %1.4e (sec)]]\n",t1-t0);

  *_de = de;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dCreateMaterialPoints(pTatinCtx ctx,DM dav)
{
  DataBucket     db;
  DataEx         ex;
  PetscLogDouble t0,t1;
  PetscInt       lmx,lmy,lmz;
  PetscBool      flg;
  int            npoints;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* register marker structures here */
  PetscTime(&t0);
  DataBucketCreate(&db);
  DataBucketRegisterField(db,MPntStd_classname,    sizeof(MPntStd),NULL);
  DataBucketRegisterField(db,MPntPStokes_classname,sizeof(MPntPStokes),NULL);
  DataBucketFinalize(db);
    DataBucketRegisterField(db,MPntPStokesPl_classname,sizeof(MPntPStokesPl),NULL);
  DataBucketFinalize(db);

  /* Choose type of projection (for eta and rho) */
  ctx->coefficient_projection_type = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-coefficient_projection_type",&ctx->coefficient_projection_type,&flg);CHKERRQ(ierr);
  switch (ctx->coefficient_projection_type) {
    case -1:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using null projection\n");
      break;

        /* P0 variants */
    case 0:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection [arithmetic avg]\n");
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 [arithmetic avg] not implemented");
      break;
    case 10:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection [harmonic avg]\n");
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 [harmonic avg] not implemented");
      break;
    case 20:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection [geometric avg]\n");
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 [geometric avg] not implemented");
      break;
    case 30:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection [dominant phase]\n");
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 [dominant phase] not implemented");
      break;

    case 1:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using Q1 projection\n");
      break;
    case 2:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using Q2 projection\n");
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = Q2 not implemented");
      break;
    case 3:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P1 projection\n");
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P1 not implemented");
      break;
    case 4:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using one2one projection\n");
      break;
    case 5:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using Q1 projection <sort ctx>\n");
      break;
      
    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER," -coefficient_projection_type = {0,1,2,4} implying {P0,Q1,Q2,one2on}");
      break;
  }

  /* Choose coordinate layout for material points, -mp_layout 0,1 */
  /* set initial size */
  ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
  DataBucketSetInitialSizes(db,lmx*lmy*lmz*2*2*2,1000);
  PetscTime(&t1);
  PetscPrintf(PETSC_COMM_WORLD,"[[Swarm initialization: %1.4lf (sec)]]\n", t1-t0);

  PetscTime(&t0);
  {
    /* defaults for lattice layout */
    PetscInt   Nxp[] = {2,2,2}; /* change with -lattice_layout_N{x,y,z} */
    PetscReal  perturb = 0.1;   /* change with -lattice_layout_perturb */
    /* defaults for random layout */
    PetscInt   nPerCell = 2 * 2 * 2; /* change with -random_layout_Np */

    PetscInt mplayout = 0;

    PetscOptionsGetInt(NULL,NULL,"-mp_layout",&mplayout,NULL);
    switch (mplayout) {
      case 0:
        ierr = SwarmMPntStd_CoordAssignment_LatticeLayout3d(dav,Nxp,perturb,db);CHKERRQ(ierr);
        break;
      case 1:
        ierr = SwarmMPntStd_CoordAssignment_RandomLayout3d(dav,nPerCell,db);CHKERRQ(ierr);
        break;
      case 2:
        ierr = SwarmMPntStd_CoordAssignment_GaussLayout3d(dav,db);CHKERRQ(ierr);
        break;
      default:
        ierr = SwarmMPntStd_CoordAssignment_LatticeLayout3d(dav,Nxp,perturb,db);CHKERRQ(ierr);
        break;
    }
  }
  PetscTime(&t1);
  DataBucketGetSizes(db,&npoints,NULL,NULL);
  PetscPrintf(PETSC_COMM_WORLD,"[[Swarm->coordinate assignment: %d points : %1.4lf (sec)]]\n", npoints,t1-t0);

  /* create the data exchanger need for parallel particle movement */
  ierr = SwarmDMDA3dDataExchangerCreate(dav,&ex);CHKERRQ(ierr);

  ctx->materialpoint_db = db;
  ctx->materialpoint_ex = ex;

  PetscFunctionReturn(0);
}

PetscErrorCode MaterialPointCoordinateSetUp(pTatinCtx ctx,DM da)
{
  PetscErrorCode ierr;
  DataBucket db;
  DM cda;
  Vec gcoords;
  PetscScalar *LA_gcoords;
  PetscInt nel,nen,i;
  int e,p,n_mp_points;
  const PetscInt *elnidx;
  PetscScalar Ni_p[Q2_NODES_PER_EL_3D];
  PetscScalar elcoords[3*Q2_NODES_PER_EL_3D];
  DataField      PField_std;


  PetscFunctionBegin;
  db = ctx->materialpoint_db;

  /* setup for coords */
  ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

  ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);

  DataBucketGetSizes(db,&n_mp_points,0,0);
  DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
  DataFieldGetAccess(PField_std);

  for (p=0; p<n_mp_points; p++) {
    MPntStd     *material_point;
    double      xp[3];
    double      *xi;

    DataFieldAccessPoint(PField_std,p,(void**)&material_point);

    MPntStdGetField_local_element_index(material_point,&e);
    ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);

    MPntStdGetField_local_coord(material_point,&xi);
    pTatin_ConstructNi_Q2_3D(xi,Ni_p);

    xp[0] = 0.0;
    xp[1] = 0.0;
    xp[2] = 0.0;
    for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
      xp[0] = xp[0] + Ni_p[i] * elcoords[NSD*i  ];
      xp[1] = xp[1] + Ni_p[i] * elcoords[NSD*i+1];
      xp[2] = xp[2] + Ni_p[i] * elcoords[NSD*i+2];
    }
    MPntStdSetField_global_coord(material_point,xp);
  }

  DataFieldRestoreAccess(PField_std);

  ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_ModelOutput_MPntStd(pTatinCtx ctx,const char prefix[])
{
  PetscErrorCode ierr;
  PetscLogDouble t0,t1;
  static PetscBool beenhere=PETSC_FALSE;
  char name[PETSC_MAX_PATH_LEN],pvdfilename[PETSC_MAX_PATH_LEN],vtkfilename[PETSC_MAX_PATH_LEN],pvoutputdir[PETSC_MAX_PATH_LEN],root[PETSC_MAX_PATH_LEN];
  char stepprefix[PETSC_MAX_PATH_LEN];
  PetscBool found;

  PetscFunctionBegin;
  ierr = PetscSNPrintf(root,PETSC_MAX_PATH_LEN-1,"%s",ctx->outputpath);CHKERRQ(ierr);

  ierr = PetscSNPrintf(pvoutputdir,PETSC_MAX_PATH_LEN-1,"%s/step%D",root,ctx->step);CHKERRQ(ierr);
  ierr = pTatinTestDirectory(pvoutputdir,'w',&found);CHKERRQ(ierr);
  if (!found) { ierr = pTatinCreateDirectory(pvoutputdir);CHKERRQ(ierr); }

  PetscTime(&t0);
  // PVD
  PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_mpoints_std.pvd",root);
  if (prefix) { PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "%s_mpoints_std.pvtu",prefix);
  } else {      PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "mpoints_std.pvtu");           }

  if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename ); }
  PetscSNPrintf(stepprefix,PETSC_MAX_PATH_LEN-1,"step%D",ctx->step);
  ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time, vtkfilename, stepprefix);CHKERRQ(ierr);
  beenhere = PETSC_TRUE;

  // PVTS + VTS
  if (prefix) { PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1,"%s_mpoints_std",prefix);
  } else {      PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1,"mpoints_std",prefix);    }

  ierr = SwarmOutputParaView_MPntStd(ctx->materialpoint_db,pvoutputdir,name);CHKERRQ(ierr);

  PetscTime(&t1);
  /*PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_mpoints_std.(pvd,pvtu,vtu): CPU time %1.2e (sec) \n", PETSC_FUNCTION_NAME,prefix,t1-t0);*/

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dCreateContext(pTatinCtx *ctx)
{
  pTatinCtx      user;
  PetscMPIInt    rank;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(&user);CHKERRQ(ierr);

  /* init */
  user->stokes_ctx = NULL;
  user->energy_ctx = NULL;
//  user->coords_ctx = NULL;

  user->pack     = NULL; /* DM composite for velocity and pressure */

  /* set defaults */
  user->restart_from_file         = PETSC_FALSE;
  user->checkpoint_every          = 1000;
  user->checkpoint_every_nsteps   = 1000;
  user->checkpoint_every_ncpumins = 90.0;
  user->checkpoint_disable        = PETSC_FALSE;

  user->mx               = 4;
  user->my               = 4;
  user->mz               = 4;
  user->use_mf_stokes    = PETSC_FALSE;
  user->solverstatistics = PETSC_FALSE;

  user->continuation_m   = 1;
  user->continuation_M   = 1;
  user->coefficient_projection_type = 1; /* Q1 */

  /* time step control */
  user->nsteps           = 1;
  user->dt_max           = 1.0e30;
  user->dt_min           = 0.0;
  user->dt               = 0.0;
  user->output_frequency = 1;
  user->time_max         = 1.0e32;
  user->time             = 0.0;
  user->step             = 0;
  user->dt_adv           = 0.0;
  user->use_constant_dt  = PETSC_FALSE;

  ierr = RheologyConstantsInitialise(&user->rheology_constants);CHKERRQ(ierr);
  ierr = MaterialConstantsCreate(&user->material_constants);CHKERRQ(ierr);
  ierr = MaterialConstantsInitialize(user->material_constants);CHKERRQ(ierr);
  ierr = PetscContainerCreate(PETSC_COMM_WORLD,&user->model_data);CHKERRQ(ierr);

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  if (rank == 0) {
    pTatinGenerateFormattedTimestamp(user->formatted_timestamp);
  }
  ierr = MPI_Bcast(user->formatted_timestamp,PETSC_MAX_PATH_LEN,MPI_CHAR,0,PETSC_COMM_WORLD);CHKERRQ(ierr);

  *ctx = user;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dDestroyContext(pTatinCtx *ctx)
{
  pTatinCtx user = *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (user->materialpoint_ex) { /* ierr = DataExView(user->materialpoint_ex);CHKERRQ(ierr); */ ierr = DataExDestroy(user->materialpoint_ex);CHKERRQ(ierr); }
  if (user->materialpoint_db) { DataBucketDestroy(&user->materialpoint_db); }

  if (user->energy_ctx) { ierr = PhysCompDestroy_Energy(&user->energy_ctx);CHKERRQ(ierr); }
  if (user->stokes_ctx) { ierr = PhysCompDestroy_Stokes(&user->stokes_ctx);CHKERRQ(ierr); }
  if (user->pack) {       ierr = DMDestroy(&user->pack);CHKERRQ(ierr); }

  /*
   if (user->Q) { ierr = QuadratureStokesDestroy(&user->Q);CHKERRQ(ierr); }
   for (e=0; e<QUAD_EDGES; e++) {
   if (user->surfQ[e]) { ierr = SurfaceQuadratureStokesDestroy(&user->surfQ[e]);CHKERRQ(ierr); }
   }
   if (user->p_bclist) { ierr = BCListDestroy(&user->p_bclist);CHKERRQ(ierr); }
   if (user->u_bclist) { ierr = BCListDestroy(&user->u_bclist);CHKERRQ(ierr); }
   if (user->dap) { ierr = DMDestroy(&user->dap);CHKERRQ(ierr); }
   if (user->dav) { ierr = DMDestroy(&user->dav);CHKERRQ(ierr); }
   */

  if (user->material_constants) { DataBucketDestroy(&user->material_constants); }

  ierr = PetscContainerDestroy(&user->model_data);CHKERRQ(ierr);

  {
    char  logfile[PETSC_MAX_PATH_LEN];
    PetscViewer viewer;

    sprintf(logfile,"%s/ptatin.petsc.log_summary-%s",user->outputpath,user->formatted_timestamp);

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,logfile,&viewer);CHKERRQ(ierr);
    ierr = PetscLogView(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  ierr = PetscLogView(user->log);CHKERRQ(ierr);
  ierr = pTatinLogCloseFile(user);CHKERRQ(ierr);

  ierr = PetscFree(user);CHKERRQ(ierr);

  *ctx = NULL;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinCtxGetModelData(pTatinCtx ctx,const char name[],void **data)
{
  PetscErrorCode ierr;
  PetscContainer container;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)ctx->model_data,name,(PetscObject*)&container);CHKERRQ(ierr);
  if (!container) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"No data with name \"%s\" was composed with ctx->model_data",name);
  ierr = PetscContainerGetPointer(container,data);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinCtxAttachModelData(pTatinCtx ctx,const char name[],void *data)
{
  PetscContainer container;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscContainerCreate(PETSC_COMM_WORLD,&container);CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(container,(void*)data);CHKERRQ(ierr);

  ierr = PetscObjectCompose((PetscObject)ctx->model_data,name,(PetscObject)container);CHKERRQ(ierr);
  ierr = PetscContainerDestroy(&container);CHKERRQ(ierr); /* decrement ref counter */

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dSetFromOptions(pTatinCtx ctx)
{
  char           optionsfile[PETSC_MAX_PATH_LEN];
  PetscInt       mx3 = 4;
  PetscBool      flg;
  PetscErrorCode ierr;

  /* parse options */
  ierr = PetscOptionsGetInt(NULL,NULL,"-mx",&ctx->mx,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-my",&ctx->my,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-mz",&ctx->mz,&flg);CHKERRQ(ierr);
  flg = PETSC_FALSE; ierr = PetscOptionsGetInt(NULL,NULL,"-mx3",&mx3,&flg);CHKERRQ(ierr);
  if (flg) {
    ctx->mx = mx3;
    ctx->my = mx3;
    ctx->mz = mx3;
  }

  ierr = PetscOptionsGetBool(NULL,NULL,"-use_mf_stokes",&ctx->use_mf_stokes,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-with_statistics",&ctx->solverstatistics,&flg);CHKERRQ(ierr);

  flg = PETSC_FALSE;
  ierr = PetscOptionsGetString(NULL,NULL,"-output_path",ctx->outputpath,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (flg == PETSC_FALSE) {
    sprintf(ctx->outputpath,"./output");
  }
  ierr = pTatinCreateDirectory(ctx->outputpath);CHKERRQ(ierr);

  /* checkpointing */
  ierr = PetscOptionsGetInt(NULL,NULL,"-checkpoint_every",&ctx->checkpoint_every,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-checkpoint_every_nsteps",&ctx->checkpoint_every_nsteps,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-checkpoint_every_ncpumins",&ctx->checkpoint_every_ncpumins,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-checkpoint_disable",&ctx->checkpoint_disable,&flg);CHKERRQ(ierr);

  /* time stepping */
  ierr = PetscOptionsGetInt(NULL,NULL,"-nsteps",&ctx->nsteps,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt_min",&ctx->dt_min,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt_max",&ctx->dt_max,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-time_max",&ctx->time_max,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-output_frequency",&ctx->output_frequency,&flg);CHKERRQ(ierr);
  {
    PetscReal constant_dt;

    ierr = PetscOptionsGetReal(NULL,NULL,"-constant_dt",&constant_dt,&flg);CHKERRQ(ierr);
    if (flg) {
      ctx->use_constant_dt = PETSC_TRUE;
      ctx->constant_dt     = constant_dt;
    }
  }

  /* open log file */
  ierr = pTatinLogOpenFile(ctx);CHKERRQ(ierr);
  ierr = pTatinLogHeader(ctx);CHKERRQ(ierr);

  sprintf(optionsfile,"%s/ptatin.options-%s",ctx->outputpath,ctx->formatted_timestamp);
  ierr = pTatinWriteOptionsFile(optionsfile);CHKERRQ(ierr);

  sprintf(optionsfile,"%s/ptatin.options",ctx->outputpath);
  ierr = pTatinWriteOptionsFile(optionsfile);CHKERRQ(ierr);

//  ierr = pTatinModelLoad(ctx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinModelLoad(pTatinCtx ctx)
{
  pTatinModel model;
  PetscBool flgname;
  char modelname[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  flgname = PETSC_FALSE;
  ierr = PetscOptionsGetString(NULL,NULL,"-ptatin_model",modelname,PETSC_MAX_PATH_LEN-1,&flgname);CHKERRQ(ierr);
  if (flgname) {
    ierr = pTatinModelGetByName(modelname,&model);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"  [pTatinModel]: -ptatin_model \"%s\" was detected\n",model->model_name);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"  [pTatinModel]: -ptatin_model wasn't specified - running boring \"template\" model\n");
    ierr = pTatinModelGetByName("template",&model);CHKERRQ(ierr);
  }

  model->ptat_ctx = ctx;
  ctx->model = model;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinGetTime(pTatinCtx ctx,PetscReal *time)
{
  if (time) { *time = ctx->time; }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinGetTimestep(pTatinCtx ctx,PetscReal *dt)
{
  if (dt) { *dt = ctx->dt; }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinGetMaterialPoints(pTatinCtx ctx,DataBucket *db,DataEx *de)
{
  if (db) { *db = ctx->materialpoint_db; }
  if (de) { *de = ctx->materialpoint_ex; }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinGetMaterialConstants(pTatinCtx ctx,DataBucket *db)
{
  if (db) { *db = ctx->material_constants; }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinGetModel(pTatinCtx ctx,pTatinModel *m)
{
  if (m) { *m = ctx->model; }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinGetRheology(pTatinCtx ctx,RheologyConstants **r)
{
  if (r) { *r = &ctx->rheology_constants; }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinGetStokesContext(pTatinCtx ctx,PhysCompStokes *s)
{
  if (s) { *s = ctx->stokes_ctx; }
  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dCheckpoint(pTatinCtx ctx,Vec X,const char prefix[])
{
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatin3dCheckpoint is deprecated");
  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dCheckpointManager(pTatinCtx ctx,Vec Xs)
{
  PetscErrorCode   ierr;
  PetscInt         checkpoint_every;
  PetscInt         checkpoint_every_nsteps,step;
  double           checkpoint_every_ncpumins, max_current_cpu_time, current_cpu_time;
  static double    last_cpu_time = 0.0;
  /*PetscBool      skip_existence_test = PETSC_TRUE;*/
  PhysCompStokes   stokes = NULL;
  PetscBool        energy_activated;
  PhysCompEnergy   energy = NULL;
  Vec              Xe = NULL;
  DM               dmv,dmp,dmstokes = NULL,dmenergy = NULL;
  PetscBool        exists,write_step_checkpoint;
  char             checkpoints_basedir[PETSC_MAX_PATH_LEN];
  char             test_dir[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  ierr = pTatinGetStokesContext(ctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,&dmp);CHKERRQ(ierr);
  ierr = pTatinContextValid_Energy(ctx,&energy_activated);CHKERRQ(ierr);
  if (energy_activated) {
    ierr = pTatinGetContext_Energy(ctx,&energy);CHKERRQ(ierr);
    dmenergy = energy->daT;
    ierr = pTatinPhysCompGetData_Energy(ctx,&Xe,NULL);CHKERRQ(ierr);
  }

  ierr = PetscSNPrintf(checkpoints_basedir,PETSC_MAX_PATH_LEN-1,"%s/checkpoints",ctx->outputpath);CHKERRQ(ierr);
  ierr = pTatinTestDirectory(checkpoints_basedir,'w',&exists);CHKERRQ(ierr);
  if (!exists) {
    ierr = pTatinCreateDirectory(checkpoints_basedir);CHKERRQ(ierr);
  }

  step                      = ctx->step;
  checkpoint_every          = ctx->checkpoint_every;
  checkpoint_every_nsteps   = ctx->checkpoint_every_nsteps;
  checkpoint_every_ncpumins = ctx->checkpoint_every_ncpumins;

  /* -------------------------------------- */
  /* check one - this file has a fixed name */
  if (step%checkpoint_every == 0) {

    ierr = PetscSNPrintf(test_dir,PETSC_MAX_PATH_LEN-1,"%s/checkpoints/default",ctx->outputpath);CHKERRQ(ierr);
    ierr = pTatinTestDirectory(test_dir,'w',&exists);CHKERRQ(ierr);
    if (!exists) {
      ierr = pTatinCreateDirectory(test_dir);CHKERRQ(ierr);
    }
    PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager[every]: Writing to dir %s\n",test_dir);
    // call checkpoint routine //
    //ierr = pTatin3dCheckpoint(ctx,X,NULL);CHKERRQ(ierr);
    ierr = pTatinCtxCheckpointWrite(ctx,test_dir,NULL,dmstokes,dmenergy,0,NULL,NULL,Xs,Xe,NULL,NULL);CHKERRQ(ierr);
  }

  write_step_checkpoint = PETSC_FALSE;

  /* -------------------------------------------------------------------- */
  /* check three - look at cpu time and decide if we need to write or not */
  PetscTime(&current_cpu_time);
  ierr = MPI_Allreduce(&current_cpu_time,&max_current_cpu_time,1,MPIU_PETSCLOGDOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
  max_current_cpu_time = max_current_cpu_time/60.0; /* convert sec to mins */

  if (max_current_cpu_time > last_cpu_time + checkpoint_every_ncpumins) {

    PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager[checkpoint_every_ncpumins]: Activated\n");
    write_step_checkpoint = PETSC_TRUE;

    last_cpu_time = max_current_cpu_time;
  }

  /* ----------------------------------------------------------------- */
  /* check two - these files have a file name related to the time step */
  if (step%checkpoint_every_nsteps == 0) {

    PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager[checkpoint_every_nsteps]: Activated\n");
    write_step_checkpoint = PETSC_TRUE;
  }

  /* if either the cpu based or step based checks returned true, write a checkpoint file */
  if (write_step_checkpoint) {
    PetscLogDouble time[2];
    char           restartfile[PETSC_MAX_PATH_LEN];
    char           restartstring[PETSC_MAX_PATH_LEN];
    PetscMPIInt    rank;

    ierr = PetscSNPrintf(test_dir,PETSC_MAX_PATH_LEN-1,"%s/step%d",checkpoints_basedir,step);CHKERRQ(ierr);

    ierr = pTatinTestDirectory(test_dir,'w',&exists);CHKERRQ(ierr);
    if (!exists) {
      ierr = pTatinCreateDirectory(test_dir);CHKERRQ(ierr);
    }
    PetscPrintf(PETSC_COMM_WORLD,"CheckpointManager: Writing to dir %s\n",test_dir);
    PetscTime(&time[0]);
    ierr = pTatinCtxCheckpointWrite(ctx,test_dir,NULL,dmstokes,dmenergy,0,NULL,NULL,Xs,Xe,NULL,NULL);CHKERRQ(ierr);
    PetscTime(&time[1]);
    ierr = pTatinLogBasicCPUtime(ctx,"Checkpoint.write()",time[1]-time[0]);CHKERRQ(ierr);

    /* write out a default string for restarting the job */
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = PetscSNPrintf(restartfile,PETSC_MAX_PATH_LEN-1,"%s/restart.default",ctx->outputpath);CHKERRQ(ierr);
    ierr = PetscSNPrintf(restartstring,PETSC_MAX_PATH_LEN-1,"-restart_directory %s/checkpoints/step%d",ctx->outputpath,step);CHKERRQ(ierr);
    if (rank == 0) {
      FILE *fp;
      fp = fopen(restartfile,"w");
      fprintf(fp,"%s",restartstring);
      fclose(fp);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinRestart_Initialize(pTatinCtx ctx,void *data)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinRestart_ApplyInitialMeshGeometry(pTatinCtx ctx,void *data)
{
  PetscFunctionBegin;
  /* load coordinates of the velocity DMDA */
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatinRestart_ApplyInitialMeshGeometry is deprecated");
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinRestart_ApplyInitialMaterialGeometry(pTatinCtx ctx,void *data)
{
  PetscFunctionBegin;
  /* load material points from file */
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatinRestart_ApplyInitialMaterialGeometry is deprecated");
  PetscFunctionReturn(0);
}

PetscErrorCode pTatinRestart_ApplyInitialSolution(pTatinCtx ctx,Vec X,void *data)
{
  PetscFunctionBegin;
  /* load state vectors (u,p,T) from file */
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatinRestart_ApplyInitialSolution is deprecated");
  PetscFunctionReturn(0);
}

PetscErrorCode  DMCoarsenHierarchy2_DA(DM da,PetscInt nlevels,DM dac[])
{
  PetscErrorCode ierr;
  PetscInt       i,n,*refx,*refy,*refz;
  PetscBool      view;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(da,DM_CLASSID,1);
  if (nlevels < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_OUTOFRANGE,"nlevels cannot be negative");
  if (nlevels == 0) PetscFunctionReturn(0);
  PetscValidPointer(dac,3);

  /* Get refinement factors, defaults taken from the coarse DMDA */
  ierr = PetscMalloc3(nlevels,&refx,nlevels,&refy,nlevels,&refz);CHKERRQ(ierr);
  for (i=0; i<nlevels; i++) {
    ierr = DMDAGetRefinementFactor(da,&refx[i],&refy[i],&refz[i]);CHKERRQ(ierr);
  }
  n = nlevels;
  ierr = PetscOptionsGetIntArray(NULL,((PetscObject)da)->prefix,"-da_refine_hierarchy_x",refx,&n,NULL);CHKERRQ(ierr);
  n = nlevels;
  ierr = PetscOptionsGetIntArray(NULL,((PetscObject)da)->prefix,"-da_refine_hierarchy_y",refy,&n,NULL);CHKERRQ(ierr);
  n = nlevels;
  ierr = PetscOptionsGetIntArray(NULL,((PetscObject)da)->prefix,"-da_refine_hierarchy_z",refz,&n,NULL);CHKERRQ(ierr);

  ierr = DMDASetCoarseningFactor(da,refx[nlevels-1],refy[nlevels-1],refz[nlevels-1]);CHKERRQ(ierr);
  ierr = DMDASetRefinementFactor(da,refx[nlevels-1],refy[nlevels-1],refz[nlevels-1]);CHKERRQ(ierr);
  ierr = DMCoarsen(da,PetscObjectComm((PetscObject)da),&dac[0]);CHKERRQ(ierr);
  for (i=1; i<nlevels; i++) {
    ierr = DMDASetCoarseningFactor(dac[i-1],refx[nlevels-1-i],refy[nlevels-1-i],refz[nlevels-1-i]);CHKERRQ(ierr);
    ierr = DMDASetRefinementFactor(dac[i-1],refx[nlevels-1-i],refy[nlevels-1-i],refz[nlevels-1-i]);CHKERRQ(ierr);
    ierr = DMCoarsen(dac[i-1],PetscObjectComm((PetscObject)da),&dac[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree3(refx,refy,refz);CHKERRQ(ierr);

  view = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,((PetscObject)da)->prefix,"-da_view_hierarchy",&view,NULL);CHKERRQ(ierr);
  if (view) {
    char levelname[PETSC_MAX_PATH_LEN];

    ierr = DMDAViewPetscVTK(da,NULL,"dav_fine.vtk");CHKERRQ(ierr);
    for (i=0; i<nlevels; i++) {
      PetscSNPrintf(levelname,PETSC_MAX_PATH_LEN-1,"dav_level%D.vtk",nlevels-1-i); /* do shift to make 0 named as the corsest */
      ierr = DMDAViewPetscVTK(dac[i],NULL,levelname);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/*
 For multiple, decoupled physics, several time scales are likely to be present.
 In general, we are unlikely to be able to know, from each physics, what all the time scales are going to be.

 Probably what we will do is this

 Solve physics_1,
 Compute appropriate time step for physics_1, dt_1

 Solve physics_2,
 Compute appropriate time step for physics_2, dt_2

 Solve physics_3,
 Compute appropriate time step for physics_3, dt_3

 dt = min( dt_1, dt_2, dt_3 ) subject to some global min/max cut offs.
 */
PetscErrorCode pTatin_SetTimestep(pTatinCtx ctx,const char timescale_name[],PetscReal dt_trial)
{
  PetscReal dt_current;

  PetscFunctionBegin;
  if (timescale_name) {
    PetscPrintf(PETSC_COMM_WORLD,"  TimeStep control(%.20s):",timescale_name);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"  TimeStep control:");
  }

  if (ctx->use_constant_dt) {
    ctx->dt = ctx->constant_dt;
    PetscPrintf(PETSC_COMM_WORLD," | using constant time step ==>> dt used = %1.4e |\n", ctx->dt );

    PetscFunctionReturn(0);
  }

  dt_current = ctx->dt;
  if (dt_trial < dt_current) {
    PetscPrintf(PETSC_COMM_WORLD," | current = %1.4e : trial = %1.4e [accepted]", dt_current, dt_trial);
    ctx->dt = dt_trial;
  }
  else {
    PetscPrintf(PETSC_COMM_WORLD," | current = %1.4e : trial = %1.4e", dt_current, dt_trial);
  }
  dt_current = ctx->dt;

  /* apply limiters */
  if (dt_current < ctx->dt_min) {
    dt_current = ctx->dt_min;
    PetscPrintf(PETSC_COMM_WORLD," | dt < dt_min : limited to = %1.4e", dt_current );
  }
  if (dt_current > ctx->dt_max) {
    dt_current = ctx->dt_max;
    PetscPrintf(PETSC_COMM_WORLD," | dt > dt_max : restricted to = %1.4e", dt_current );
  }
  ctx->dt = dt_current;
  PetscPrintf(PETSC_COMM_WORLD," | ==>> dt used = %1.4e |\n", ctx->dt );

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinCtxCheckpointWrite(pTatinCtx ctx,const char path[],const char prefix[],
                                        DM dms,DM dme,
                                        PetscInt nfields,const char *dmnames[],DM dmlist[],
                                        Vec Xs,Vec Xe,const char *fieldnames[],Vec veclist[])
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscMPIInt commsize,commrank;
  char jfilename[PETSC_MAX_PATH_LEN];
  char vfilename[3][PETSC_MAX_PATH_LEN],checkpoint_prefix[3][PETSC_MAX_PATH_LEN];
  PetscBool energy_activated;
  PhysCompStokes stokes = NULL;
  PhysCompEnergy energy = NULL;
  DM dmv,dmp;
  DataBucket materialpoint_db = NULL,material_constants_db = NULL;


  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)ctx->pack,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);

  if (path && prefix) SETERRQ(comm,PETSC_ERR_SUP,"Intended support for either ${path}/FILENAME or ${prefix}_FILENAME");
  if (prefix) SETERRQ(comm,PETSC_ERR_SUP,"Current implementation only supports ${path}/FILENAME format");
  if (nfields > 0) SETERRQ(comm,PETSC_ERR_SUP,"Only support for stokes + energy");

  if (path) {
    PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s/ptatin3dctx.json",path);

    PetscSNPrintf(vfilename[0],PETSC_MAX_PATH_LEN-1,"%s/ptatinstate_stokes_Xv.pbvec",path);
    PetscSNPrintf(vfilename[1],PETSC_MAX_PATH_LEN-1,"%s/ptatinstate_stokes_Xp.pbvec",path);
    PetscSNPrintf(vfilename[2],PETSC_MAX_PATH_LEN-1,"%s/ptatinstate_energy_Xt.pbvec",path);

    PetscSNPrintf(checkpoint_prefix[0],PETSC_MAX_PATH_LEN-1,"%s/stokes_v",path);
    PetscSNPrintf(checkpoint_prefix[1],PETSC_MAX_PATH_LEN-1,"%s/materialpoint",path);
    PetscSNPrintf(checkpoint_prefix[2],PETSC_MAX_PATH_LEN-1,"%s/materialconstants",path);
  } else {
    SETERRQ(comm,PETSC_ERR_SUP,"Current implementation only supports ${path}/FILENAME format");
  }

  /* stokes */
  ierr = pTatinGetStokesContext(ctx,&stokes);CHKERRQ(ierr);
  ierr = PhysCompStokesGetDMs(stokes,&dmv,&dmp);CHKERRQ(ierr);

  ierr = DMDACheckpointWrite(dmv,checkpoint_prefix[0]);CHKERRQ(ierr);
  {
    Vec velocity,pressure;

    ierr = DMCompositeGetAccess(dms,Xs,&velocity,&pressure);CHKERRQ(ierr);

    ierr = DMDAWriteVectorToFile(velocity,vfilename[0],PETSC_FALSE);CHKERRQ(ierr);
    ierr = DMDAWriteVectorToFile(pressure,vfilename[1],PETSC_FALSE);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(dms,Xs,&velocity,&pressure);CHKERRQ(ierr);
  }

  /* energy */
  ierr = pTatinContextValid_Energy(ctx,&energy_activated);CHKERRQ(ierr);
  if (energy_activated) {
    ierr = pTatinGetContext_Energy(ctx,&energy);CHKERRQ(ierr);
    ierr = PhysCompCheckpointWrite_Energy(energy,PETSC_FALSE,path,NULL);CHKERRQ(ierr);

    ierr = DMDAWriteVectorToFile(Xe,vfilename[2],PETSC_FALSE);CHKERRQ(ierr);
  }

  /* material points */
  ierr = pTatinGetMaterialPoints(ctx,&materialpoint_db,NULL);CHKERRQ(ierr);
  DataBucketView(PETSC_COMM_WORLD,materialpoint_db,checkpoint_prefix[1],DATABUCKET_VIEW_NATIVE);

  /* material constants */
  /* material_constants_db is a redundant object, e.g. it is identical on all ranks */
  /* Hence, we let only 1 rank write out the data file during checkpoint.write() */
  ierr = pTatinGetMaterialConstants(ctx,&material_constants_db);CHKERRQ(ierr);
  if (commrank == 0) {
    DataBucketView(PETSC_COMM_SELF,material_constants_db,checkpoint_prefix[2],DATABUCKET_VIEW_NATIVE);
  }

  /* user state */
  /*
  for (f=0; f<nfields; f++) {
    char uprefix[PETSC_MAX_PATH_LEN];
    char fprefix[PETSC_MAX_PATH_LEN];

    PetscSNPrintf(uprefix,PETSC_MAX_PATH_LEN-1,"%s/%s",path,dmnames[f]);
    ierr = DMDACheckpointWrite(dmlist[f],uprefix);CHKERRQ(ierr);

    PetscSNPrintf(fprefix,PETSC_MAX_PATH_LEN-1,"%s/ptatinstate_%s_%s.pbvec",path,dmnames[f],fieldnames[f]);
    ierr = DMDAWriteVectorToFile(veclist[f],fprefix,PETSC_FALSE);CHKERRQ(ierr);
  }
  */

  if (commrank == 0) {
    cJSON *jso_file = NULL,*jso_ptat = NULL,*jso_state,*jso_object,*content;
    char relpathtofile[PETSC_MAX_PATH_LEN];

    /* create json meta data file */
    jso_file = cJSON_CreateObject();

    jso_ptat = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_file,"pTatinCtx",jso_ptat);

    content = cJSON_CreateInt((int)commsize);  cJSON_AddItemToObject(jso_ptat,"commSize",content);
    content = cJSON_CreateInt((int)ctx->mx);  cJSON_AddItemToObject(jso_ptat,"mx",content);
    content = cJSON_CreateInt((int)ctx->my);  cJSON_AddItemToObject(jso_ptat,"my",content);
    content = cJSON_CreateInt((int)ctx->mz);  cJSON_AddItemToObject(jso_ptat,"mz",content);

    content = cJSON_CreateBool((int)ctx->restart_from_file);        cJSON_AddItemToObject(jso_ptat,"restartFromFile",content);
    content = cJSON_CreateString(ctx->restart_dir);                 cJSON_AddItemToObject(jso_ptat,"restartPath",content);
    content = cJSON_CreateInt((int)ctx->checkpoint_every);          cJSON_AddItemToObject(jso_ptat,"checkpointEvery",content);
    content = cJSON_CreateInt((int)ctx->checkpoint_every_nsteps);   cJSON_AddItemToObject(jso_ptat,"checkpointEveryNSteps",content);
    content = cJSON_CreateDouble(ctx->checkpoint_every_ncpumins);   cJSON_AddItemToObject(jso_ptat,"checkpointEveryNCPUMins",content);

    content = cJSON_CreateBool((int)ctx->use_mf_stokes);  cJSON_AddItemToObject(jso_ptat,"useMFStokes",content);

    content = cJSON_CreateString(ctx->formatted_timestamp);            cJSON_AddItemToObject(jso_ptat,"formattedTimestamp",content);
    content = cJSON_CreateString(ctx->outputpath);                     cJSON_AddItemToObject(jso_ptat,"outputPath",content);
    content = cJSON_CreateInt((int)ctx->coefficient_projection_type);  cJSON_AddItemToObject(jso_ptat,"coefficientProjectionType",content);

    content = cJSON_CreateBool((int)ctx->solverstatistics);  cJSON_AddItemToObject(jso_ptat,"solverStatistics",content);
    content = cJSON_CreateInt((int)ctx->continuation_m);     cJSON_AddItemToObject(jso_ptat,"continuation_m",content);
    content = cJSON_CreateInt((int)ctx->continuation_M);     cJSON_AddItemToObject(jso_ptat,"continuation_M",content);

    content = cJSON_CreateInt((int)ctx->step);     cJSON_AddItemToObject(jso_ptat,"timeStep",content);
    content = cJSON_CreateInt((int)ctx->nsteps);   cJSON_AddItemToObject(jso_ptat,"timeStepMax",content);

    content = cJSON_CreateDouble(ctx->dt);          cJSON_AddItemToObject(jso_ptat,"timeStepSize",content);
    content = cJSON_CreateDouble(ctx->dt_max);      cJSON_AddItemToObject(jso_ptat,"timeStepSizeMax",content);
    content = cJSON_CreateDouble(ctx->dt_min);      cJSON_AddItemToObject(jso_ptat,"timeStepSizeMin",content);
    content = cJSON_CreateDouble(ctx->dt_adv);      cJSON_AddItemToObject(jso_ptat,"timeStepSizeAdv",content);
    content = cJSON_CreateDouble(ctx->constant_dt); cJSON_AddItemToObject(jso_ptat,"constantTimeStepSize",content);
    content = cJSON_CreateBool((int)ctx->use_constant_dt);  cJSON_AddItemToObject(jso_ptat,"useConstantTimeStepSize",content);

    content = cJSON_CreateDouble(ctx->time);                cJSON_AddItemToObject(jso_ptat,"time",content);
    content = cJSON_CreateDouble(ctx->time_max);            cJSON_AddItemToObject(jso_ptat,"timeMax",content);

    content = cJSON_CreateInt((int)ctx->output_frequency);  cJSON_AddItemToObject(jso_ptat,"outputFrequency",content);

    /* references to data files */
    /* PhysCompStokes */
    jso_object = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_ptat,"stokes->gravity_vector",jso_object);
    content = cJSON_CreateString("double");          cJSON_AddItemToObject(jso_object,"ctype",content);
    content = cJSON_CreateInt((int)3);               cJSON_AddItemToObject(jso_object,"length",content);
    /* todo - change this to base64 encoded */
    content = cJSON_CreateString("ascii");           cJSON_AddItemToObject(jso_object,"dataFormat",content);
    content = cJSON_CreateDoubleArray(stokes->gravity_vector,3);   cJSON_AddItemToObject(jso_object,"data",content);
    content = cJSON_CreateDoubleArray(stokes->gravity_vector,3);   cJSON_AddItemToObject(jso_object,"data[ascii]",content);

    /* DMDA for velocity */
    jso_object = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_ptat,"stokes->dmv",jso_object);
    content = cJSON_CreateString("DMDA");                cJSON_AddItemToObject(jso_object,"ctype",content);
    content = cJSON_CreateString("json-meta");           cJSON_AddItemToObject(jso_object,"dataFormat",content);
    //content = cJSON_CreateString(checkpoint_prefix[0]);  cJSON_AddItemToObject(jso_object,"prefix",content);
    PetscSNPrintf(relpathtofile,PETSC_MAX_PATH_LEN-1,"%s_dmda.json",checkpoint_prefix[0]);
    content = cJSON_CreateString(relpathtofile);         cJSON_AddItemToObject(jso_object,"fileName",content);

    /* PhysCompEnergy */
    if (energy_activated) {
      jso_object = cJSON_CreateObject();
      cJSON_AddItemToObject(jso_ptat,"energy",jso_object);
      content = cJSON_CreateString("PhysCompEnergy");      cJSON_AddItemToObject(jso_object,"ctype",content);
      content = cJSON_CreateString("json-meta");           cJSON_AddItemToObject(jso_object,"dataFormat",content);
      //content = cJSON_CreateString(path);  cJSON_AddItemToObject(jso_object,"prefix",content);
      PetscSNPrintf(relpathtofile,PETSC_MAX_PATH_LEN-1,"%s/physcomp_energy.json",path);
      content = cJSON_CreateString(relpathtofile);         cJSON_AddItemToObject(jso_object,"fileName",content);
    }

    /* Material points */
    jso_object = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_ptat,"materialpoint_db",jso_object);
    content = cJSON_CreateString("DataBucket");          cJSON_AddItemToObject(jso_object,"ctype",content);
    content = cJSON_CreateString("native");              cJSON_AddItemToObject(jso_object,"dataFormat",content);
    //content = cJSON_CreateString(checkpoint_prefix[1]);  cJSON_AddItemToObject(jso_object,"prefix",content);
    PetscSNPrintf(relpathtofile,PETSC_MAX_PATH_LEN-1,"%s_db.json",checkpoint_prefix[1]);
    content = cJSON_CreateString(relpathtofile);         cJSON_AddItemToObject(jso_object,"fileName",content);

    /* Material constants */
    jso_object = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_ptat,"material_constants",jso_object);
    content = cJSON_CreateString("DataBucket");          cJSON_AddItemToObject(jso_object,"ctype",content);
    content = cJSON_CreateString("native");              cJSON_AddItemToObject(jso_object,"dataFormat",content);
    //content = cJSON_CreateString(checkpoint_prefix[2]);  cJSON_AddItemToObject(jso_object,"prefix",content);
    PetscSNPrintf(relpathtofile,PETSC_MAX_PATH_LEN-1,"%s_db.json",checkpoint_prefix[2]);
    content = cJSON_CreateString(relpathtofile);         cJSON_AddItemToObject(jso_object,"fileName",content);

    /* State vectors */
    jso_state = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_ptat,"stokes.Xs->Xv",jso_state);
    content = cJSON_CreateString("Vec");                      cJSON_AddItemToObject(jso_state,"ctype",content);
    content = cJSON_CreateString("petsc-binary");             cJSON_AddItemToObject(jso_state,"dataFormat",content);
    content = cJSON_CreateString(vfilename[0]);               cJSON_AddItemToObject(jso_state,"fileName",content);

    jso_state = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_ptat,"stokes.Xs->Xp",jso_state);
    content = cJSON_CreateString("Vec");                      cJSON_AddItemToObject(jso_state,"ctype",content);
    content = cJSON_CreateString("petsc-binary");             cJSON_AddItemToObject(jso_state,"dataFormat",content);
    content = cJSON_CreateString(vfilename[1]);               cJSON_AddItemToObject(jso_state,"fileName",content);

    if (energy_activated) {
      jso_state = cJSON_CreateObject();
      cJSON_AddItemToObject(jso_ptat,"energy.Xt",jso_state);
      content = cJSON_CreateString("Vec");                    cJSON_AddItemToObject(jso_state,"ctype",content);
      content = cJSON_CreateString("petsc-binary");           cJSON_AddItemToObject(jso_state,"dataFormat",content);
      content = cJSON_CreateString(vfilename[2]);             cJSON_AddItemToObject(jso_state,"fileName",content);
    }

    /*
    jso_state = cJSON_CreateArray();
    cJSON_AddItemToObject(jso_ptat,"userstate",jso_state);

    jso_object = cJSON_CreateObject();
    content = cJSON_CreateString();            cJSON_AddItemToObject(jso_object,"points",content);
    cJSON_AddItemToArray(jso_state,jso_object);
    */

    /* write json meta data file */
    {
      FILE *fp;
      char *jbuff = cJSON_Print(jso_file);

      fp = fopen(jfilename,"w");
      if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to open file %s",jfilename);
      fprintf(fp,"%s\n",jbuff);
      fclose(fp);
      free(jbuff);
    }

    cJSON_Delete(jso_file);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dLoadContext_FromFile(pTatinCtx *_ctx)
{
  pTatinCtx ctx;
  PetscErrorCode  ierr;
  char restart_dir[PETSC_MAX_PATH_LEN];
  PetscBool flg,found;
  PetscMPIInt commrank;
  char jfilename[PETSC_MAX_PATH_LEN];
  cJSON *jfile = NULL,*jptat = NULL,*jobj = NULL;
  MPI_Comm comm;
  PetscInt csize;
  char field_string[PETSC_MAX_PATH_LEN];


  PetscFunctionBegin;
  comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);

  ierr = pTatin3dCreateContext(&ctx);CHKERRQ(ierr);

  restart_dir[0] = '\0';
  ierr = PetscOptionsGetString(NULL,NULL,"-restart_directory",restart_dir,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Failed to locate essential option -restart_directory");
  PetscPrintf(PETSC_COMM_WORLD,"[pTatin] Found -restart_directory: %s\n",restart_dir);

  /* populate context with content from JSON file */
  PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s/ptatin3dctx.json",restart_dir);
  PetscPrintf(PETSC_COMM_WORLD,"[pTatin] Using checkpoint file: %s\n",jfilename);
  if (commrank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open JSON file \"%s\"",jfilename);
    jptat = cJSON_GetObjectItem(jfile,"pTatinCtx");
  }

  ierr = cJSONGetPetscInt(comm,jptat,"commSize",&csize,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscInt(comm,jptat,"mx",&ctx->mx,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscInt(comm,jptat,"my",&ctx->my,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscInt(comm,jptat,"mz",&ctx->mz,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscBool(comm,jptat,"restartFromFile",&ctx->restart_from_file,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscString(comm,jptat,"restartPath",ctx->restart_dir,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscInt(comm,jptat,"checkpointEvery",&ctx->checkpoint_every,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscInt(comm,jptat,"checkpointEveryNSteps",&ctx->checkpoint_every_nsteps,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscReal(comm,jptat,"checkpointEveryNCPUMins",&ctx->checkpoint_every_ncpumins,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscBool(comm,jptat,"useMFStokes",&ctx->use_mf_stokes,&found);CHKERRQ(ierr);

  /*ierr = cJSONGetPetscString(comm,jptat,"formattedTimestamp",ctx->formatted_timestamp,&found);CHKERRQ(ierr);*/
  ierr = cJSONGetPetscString(comm,jptat,"outputPath",ctx->outputpath,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscInt(comm,jptat,"coefficientProjectionType",&ctx->coefficient_projection_type,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscBool(comm,jptat,"solverStatistics",&ctx->solverstatistics,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscInt(comm,jptat,"continuation_m",&ctx->continuation_m,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscInt(comm,jptat,"continuation_M",&ctx->continuation_M,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscInt(comm,jptat,"timeStep",&ctx->step,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscInt(comm,jptat,"timeStepMax",&ctx->nsteps,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscReal(comm,jptat,"timeStepSize",&ctx->dt,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscReal(comm,jptat,"timeStepSizeMax",&ctx->dt_max,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscReal(comm,jptat,"timeStepSizeMin",&ctx->dt_min,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscReal(comm,jptat,"timeStepSizeAdv",&ctx->dt_adv,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscReal(comm,jptat,"constantTimeStepSize",&ctx->constant_dt,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscBool(comm,jptat,"useConstantTimeStepSize",&ctx->use_constant_dt,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscReal(comm,jptat,"time",&ctx->time,&found);CHKERRQ(ierr);
  ierr = cJSONGetPetscReal(comm,jptat,"timeMax",&ctx->time_max,&found);CHKERRQ(ierr);

  ierr = cJSONGetPetscInt(comm,jptat,"outputFrequency",&ctx->output_frequency,&found);CHKERRQ(ierr);

  /* load material constants from JSON file */
  if (ctx->material_constants) { DataBucketDestroy(&ctx->material_constants); }
  jobj = NULL;
  if (commrank == 0) {
    jobj = cJSON_GetObjectItem(jptat,"material_constants");
    if (!jobj) SETERRQ_JSONKEY(PETSC_COMM_SELF,"material_constants");
  }
  ierr = cJSONGetPetscString(comm,jobj,"fileName",field_string,&found);CHKERRQ(ierr);
  DataBucketLoadRedundant_NATIVE(comm,field_string,&ctx->material_constants);

  if (commrank == 0) {
    cJSON_Delete(jfile);
  }

  /* force these values */
  ctx->restart_from_file = PETSC_TRUE;
  PetscSNPrintf(ctx->restart_dir,PETSC_MAX_PATH_LEN-1,"%s",restart_dir);CHKERRQ(ierr);

  ierr = pTatinModelLoad(ctx);CHKERRQ(ierr);

  /* Over function pointers for loading model */
  //ctx->model->FP_pTatinModel_ApplyInitialSolution         = pTatinRestart_ApplyInitialSolution; /* could call users function and afterwards clobber */
  //ctx->model->FP_pTatinModel_ApplyInitialMeshGeometry     = pTatinRestart_ApplyInitialMeshGeometry; /* could call users function and afterwards clobber */
  //ctx->model->FP_pTatinModel_ApplyInitialMaterialGeometry = pTatinRestart_ApplyInitialMaterialGeometry; /* calling users function and afterwards clobbering would be expensive - avoid this */

  *_ctx = ctx;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dLoadState_FromFile(pTatinCtx ctx,DM dmstokes,DM dmenergy,Vec Xs,Vec Xt)
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscMPIInt commrank;
  char jfilename[PETSC_MAX_PATH_LEN];
  cJSON *jfile = NULL,*jptat = NULL,*jobj;
  PetscBool found,energy_activated;


  PetscFunctionBegin;
  comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s/ptatin3dctx.json",ctx->restart_dir);
  if (commrank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open JSON file \"%s\"",jfilename);
    jptat = cJSON_GetObjectItem(jfile,"pTatinCtx");
  }

  {
    char field_string_v[PETSC_MAX_PATH_LEN];
    char field_string_p[PETSC_MAX_PATH_LEN];
    Vec velocity,pressure;

    jobj = NULL;
    if (commrank == 0) { jobj = cJSON_GetObjectItem(jptat,"stokes.Xs->Xv"); if (!jobj) SETERRQ_JSONKEY(PETSC_COMM_SELF,"stokes.Xs->Xv"); }
    ierr = cJSONGetPetscString(comm,jobj,"fileName",field_string_v,&found);CHKERRQ(ierr);
    if (!found) SETERRQ_JSONKEY(comm,"fileName");

    jobj = NULL;
    if (commrank == 0) { jobj = cJSON_GetObjectItem(jptat,"stokes.Xs->Xp"); if (!jobj) SETERRQ_JSONKEY(PETSC_COMM_SELF,"stokes.Xs->Xp");}
    ierr = cJSONGetPetscString(comm,jobj,"fileName",field_string_p,&found);CHKERRQ(ierr);
    if (!found) SETERRQ_JSONKEY(comm,"fileName");

    ierr = DMCompositeGetAccess(dmstokes,Xs,&velocity,&pressure);CHKERRQ(ierr);

    ierr = VecLoadFromFile(velocity,field_string_v);CHKERRQ(ierr);
    ierr = VecLoadFromFile(pressure,field_string_p);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(dmstokes,Xs,&velocity,&pressure);CHKERRQ(ierr);
  }

  ierr = pTatinContextValid_Energy(ctx,&energy_activated);CHKERRQ(ierr);
  if (energy_activated) {
    char field_string_t[PETSC_MAX_PATH_LEN];

    jobj = NULL;
    if (commrank == 0) { jobj = cJSON_GetObjectItem(jptat,"energy.Xt"); if (!jobj) SETERRQ_JSONKEY(PETSC_COMM_SELF,"energy.Xt");}
    ierr = cJSONGetPetscString(comm,jobj,"fileName",field_string_t,&found);CHKERRQ(ierr);
    if (!found) SETERRQ_JSONKEY(comm,"fileName");

    ierr = VecLoadFromFile(Xt,field_string_t);CHKERRQ(ierr);
  }

  if (commrank == 0) {
    cJSON_Delete(jfile);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_PhysCompStokesLoad_FromFile(pTatinCtx ctx)
{
  PetscErrorCode  ierr;
  PhysCompStokes  stokes;
  PetscReal       grav[3];
  PetscInt        ncomponents;
  PetscBool       found;
  DM              dmv;
  MPI_Comm        comm;
  PetscMPIInt     commrank;
  char            jfilename[PETSC_MAX_PATH_LEN];
  cJSON           *jfile = NULL,*jptat = NULL,*jobj;


  PetscFunctionBegin;
  comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s/ptatin3dctx.json",ctx->restart_dir);
  if (commrank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open JSON file \"%s\"",jfilename);
    jptat = cJSON_GetObjectItem(jfile,"pTatinCtx");
  }

  {
    char field_string[PETSC_MAX_PATH_LEN];

    jobj = NULL;
    if (commrank == 0) { jobj = cJSON_GetObjectItem(jptat,"stokes->dmv"); if (!jobj) SETERRQ_JSONKEY(PETSC_COMM_SELF,"stokes->dmv"); }
    ierr = cJSONGetPetscString(comm,jobj,"fileName",field_string,&found);CHKERRQ(ierr);
    if (!found) SETERRQ_JSONKEY(comm,"fileName");

    ierr = DMDACheckpointLoad(comm,field_string,&dmv);CHKERRQ(ierr);
  }

  /* Default action - set gravity vector to be 0,1,0 to not break existing models which set -rho.g on the material points */
  /* Model initialize function can overload the gravity value */
  grav[0] = 0.0;
  grav[1] = 1.0;
  grav[2] = 0.0;

  {
    char field_string[PETSC_MAX_PATH_LEN];
    PetscReal *gravity_vec;
    PetscInt length;
    PetscBool is_ascii,is_base64;

    jobj = NULL;
    if (commrank == 0) {
      jobj = cJSON_GetObjectItem(jptat,"stokes->gravity_vector");
      if (!jobj) SETERRQ_JSONKEY(PETSC_COMM_SELF,"stokes->gravity_vector");
    }

    ierr = cJSONGetPetscInt(comm,jobj,"length",&length,&found);CHKERRQ(ierr);
    if (!found) SETERRQ_JSONKEY(comm,"length");

    ierr = cJSONGetPetscString(comm,jobj,"dataFormat",field_string,&found);CHKERRQ(ierr);
    if (!found) SETERRQ_JSONKEY(comm,"dataFormat");

    /* check dataFormat type */
    found = PETSC_FALSE;
    ierr = PetscStrncmp(field_string, "ascii",    5, &is_ascii);CHKERRQ(ierr);
    ierr = PetscStrncmp(field_string, "base64",    5, &is_base64);CHKERRQ(ierr);
    if (is_ascii) {
      ierr = cJSONGetPetscRealArray(comm,jobj,"data",&length,&gravity_vec,&found);CHKERRQ(ierr);
      if (!found) SETERRQ_JSONKEY(comm,"stokes->gravity_vector::data");
    } else if (is_base64) {
      SETERRQ(comm,PETSC_ERR_SUP,"Only support for reading ascii data");
    } else {
      SETERRQ(comm,PETSC_ERR_SUP,"Only support for reading ascii or base64 encoded data");
    }

    grav[0] = gravity_vec[0];
    grav[1] = gravity_vec[1];
    grav[2] = gravity_vec[2];

    ierr = PetscFree(gravity_vec);CHKERRQ(ierr);
  }

  if (commrank == 0) {
    cJSON_Delete(jfile);
  }

  ierr = PhysCompCreate_Stokes(&stokes);CHKERRQ(ierr);
  stokes->use_mf_stokes = ctx->use_mf_stokes;

  ierr = PhysCompSetup_Stokes(stokes,dmv);CHKERRQ(ierr);
  ierr = PhysCompCreateBoundaryList_Stokes(stokes);CHKERRQ(ierr);
  ierr = PhysCompCreateVolumeQuadrature_Stokes(stokes);CHKERRQ(ierr);
  ierr = PhysCompCreateSurfaceQuadrature_Stokes(stokes);CHKERRQ(ierr);

  /* Default action - set gravity vector to be 0,1,0 to not break existing models which set -rho.g on the material points */
  /* Model initialize function can overload the gravity value */
  ierr = PhysCompStokesSetGravityVector(stokes,grav);CHKERRQ(ierr);

  ncomponents = 3;
  found = PETSC_FALSE;
  PetscOptionsGetRealArray(NULL,NULL,"-stokes_gravity_vector",grav,&ncomponents,&found);
  if (found) {
    ierr = PhysCompStokesSetGravityVector(stokes,grav);CHKERRQ(ierr);
  }

  ctx->stokes_ctx = stokes;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3dLoadMaterialPoints_FromFile(pTatinCtx ctx,DM dmv)
{
  DataBucket     db;
  DataEx         ex;
  PetscLogDouble t0,t1;
  PetscBool      found;
  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscMPIInt commrank;
  char jfilename[PETSC_MAX_PATH_LEN],field_string[PETSC_MAX_PATH_LEN];
  cJSON *jfile = NULL,*jptat = NULL,*jobj;


  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)dmv,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s/ptatin3dctx.json",ctx->restart_dir);
  if (commrank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open JSON file \"%s\"",jfilename);
    jptat = cJSON_GetObjectItem(jfile,"pTatinCtx");
  }

  jobj = NULL;
  if (commrank == 0) { jobj = cJSON_GetObjectItem(jptat,"materialpoint_db"); if (!jobj) SETERRQ_JSONKEY(PETSC_COMM_SELF,"materialpoint_db"); }
  ierr = cJSONGetPetscString(comm,jobj,"fileName",field_string,&found);CHKERRQ(ierr);
  if (!found) SETERRQ_JSONKEY(comm,"fileName");

  /* Choose type of projection (for eta and rho) */
  ierr = PetscOptionsGetInt(NULL,NULL,"-coefficient_projection_type",&ctx->coefficient_projection_type,&found);CHKERRQ(ierr);
  switch (ctx->coefficient_projection_type) {
    case -1:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using null projection\n");
      break;

      /* P0 variants */
    case 0:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection [arithmetic avg]\n");
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 [arithmetic avg] not implemented");
      break;
    case 10:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection [harmonic avg]\n");
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 [harmonic avg] not implemented");
      break;
    case 20:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection [geometric avg]\n");
      //SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 [geometric avg] not implemented");
      break;
    case 30:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P0 projection [dominant phase]\n");
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P0 [dominant phase] not implemented");
      break;

    case 1:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using Q1 projection\n");
      break;
    case 2:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using Q2 projection\n");
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = Q2 not implemented");
      break;
    case 3:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using P1 projection\n");
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP," -coefficient_projection_type = P1 not implemented");
      break;
    case 4:
      PetscPrintf(PETSC_COMM_WORLD,"  MaterialPointsStokes: Using one2one projection\n");
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER," -coefficient_projection_type = {0,1,2,4} implying {P0,Q1,Q2,one2on}");
      break;
  }

  if (commrank == 0) {
    cJSON_Delete(jfile);
  }

  /* register marker structures here */
  PetscTime(&t0);
  DataBucketLoad_NATIVE(comm,field_string,&db);
  PetscTime(&t1);
  PetscPrintf(PETSC_COMM_WORLD,"[[Swarm initialization from file: %1.4lf (sec)]]\n",t1-t0);

  /* create the data exchanger need for parallel particle movement */
  ierr = SwarmDMDA3dDataExchangerCreate(dmv,&ex);CHKERRQ(ierr);

  ctx->materialpoint_db = db;
  ctx->materialpoint_ex = ex;

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinPhysCompActivate_Energy_FromFile(pTatinCtx ctx)
{
  PetscErrorCode ierr;
  PhysCompStokes stokes;
  PhysCompEnergy energy;
  MPI_Comm       comm;
  PetscMPIInt    commrank;
  char           jfilename[PETSC_MAX_PATH_LEN],field_string[PETSC_MAX_PATH_LEN];
  cJSON          *jfile = NULL,*jptat = NULL,*jobj;
  PetscBool      found;

  PetscFunctionBegin;
  comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm,&commrank);CHKERRQ(ierr);
  PetscSNPrintf(jfilename,PETSC_MAX_PATH_LEN-1,"%s/ptatin3dctx.json",ctx->restart_dir);
  if (commrank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open JSON file \"%s\"",jfilename);
    jptat = cJSON_GetObjectItem(jfile,"pTatinCtx");
  }

  jobj = NULL;
  if (commrank == 0) { jobj = cJSON_GetObjectItem(jptat,"energy"); if (!jobj) SETERRQ_JSONKEY(PETSC_COMM_SELF,"energy"); }
  ierr = cJSONGetPetscString(comm,jobj,"fileName",field_string,&found);CHKERRQ(ierr);
  if (!found) SETERRQ_JSONKEY(comm,"fileName");

  stokes = ctx->stokes_ctx;
  ierr = PhysCompLoad2_Energy(stokes->dav,field_string,&energy);CHKERRQ(ierr);

  /* Since this method is only called when restarting a job,
     the material point coefficients will have been registered when the data bucket was loaded
  */
  /*
  if (user->restart_from_file) {
  } else {
    ierr = PhysCompAddMaterialPointCoefficients_Energy(ctx->materialpoint_db);CHKERRQ(ierr);
  }
  */

  if (commrank == 0) {
    cJSON_Delete(jfile);
  }

  ctx->energy_ctx = energy;

  PetscFunctionReturn(0);
}

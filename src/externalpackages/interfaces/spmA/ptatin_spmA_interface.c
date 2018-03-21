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
 **    filename:   ptatin_spmA_interface.c
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

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_redundant.h"
#include "dmda_remesh.h"

#ifdef PTATIN_HAVE_SPMA
#include "spmA.h"
#endif

/*
SEQ
 Copy data

MPI
 special: no decomp in J

 general:
 Find ranks containing jmax
 Gather sizes on output_rank
 Send data
 Reshuffle data into natural ordering
*/
PetscErrorCode ptatin3d_DMDAAllGatherCoorJMax(DM dm,PetscMPIInt output_rank,long int *_nx,long int *_nz,double *_x[],double *_y[],double *_z[])
{
  DM         da_surf;
  double     *x,*y,*z;
  PetscInt        nx,ny,nz,si,sj,sk,si_p,ei_p,sj_p,ej_p,sk_p,ek_p;
  PetscMPIInt     rank,nproc;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&nproc);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank);CHKERRQ(ierr);

  ierr = DMDAGetInfo( dm,0,&nx,&ny,&nz,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = DMDAGetCorners( dm, &si,&sj,&sk, 0,0,0 );CHKERRQ(ierr);

  /* general gather list of nodes to collect */
  si_p = 0;
  ei_p = nx+1;

  sk_p = 0;
  ek_p = nz+1;

  sj_p = ny;
  ej_p = ny+1;

  /* on all non-"output_rank" ranks, just fetch a local piece of the domain so no uneccessary communication occurs */
  if (rank != output_rank) {
    si_p = si;
    ei_p = si+1;

    sk_p = sk;
    ek_p = sk+1;

    sj_p = sj;
    ej_p = sj+1;
  }
  ierr = DMDACreate3dRedundant( dm, si_p,ei_p, sj_p,ej_p, sk_p,ek_p, 1, &da_surf );CHKERRQ(ierr);

  x = y = z = NULL;
  if (rank == output_rank) {
    DM         cda_surf;
    Vec        coor_surf;
    DMDACoor3d ***LA_coor_surf;
    PetscInt   i,j,k;

    /* copy out coordinates */
    PetscMalloc(sizeof(double)*nx*nz,&x);
    PetscMalloc(sizeof(double)*nx*nz,&y);
    PetscMalloc(sizeof(double)*nx*nz,&z);

    ierr = DMGetCoordinates(da_surf,&coor_surf);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(da_surf,&cda_surf);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(cda_surf,coor_surf,&LA_coor_surf);CHKERRQ(ierr);
    j = 0;
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
        PetscInt idx;

        idx = i + k*nx;
        x[idx] = LA_coor_surf[k][j][i].x;

        y[idx] = LA_coor_surf[k][j][i].y;

        z[idx] = LA_coor_surf[k][j][i].z;
      }
    }
    ierr = DMDAVecRestoreArray(cda_surf,coor_surf,&LA_coor_surf);CHKERRQ(ierr);
  }

  ierr = DMDestroy(&da_surf);CHKERRQ(ierr);

  *_nx = -1;
  *_nz = -1;
  if (rank == output_rank) {
    *_nx = (long int)nx;
    *_nx = (long int)nz;
  }
  *_x  = x;
  *_y  = y;
  *_z  = z;

  PetscFunctionReturn(0);
}

/*
   Copy the chunk each rank requires into temporary array
   Send chunk
   Recieve and insert
 */
PetscErrorCode ptatin3d_DMDAAllScatterCoorJMax(PetscMPIInt intput_rank,double ymax[],DM dm)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

/*
   Interpolate mechanical model surface onto lem grid.
 */
PetscErrorCode ptatin3d_SEQLEMHelper_InterpolateM2L(DM m_surf,DM l_surf)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}


PetscErrorCode ptatin3d_SEQLEMHelper_InterpolateL2M(DM l_surf,DM m_surf)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

/*
 * ----------------------------------------------------------------------------------------------- *
 SPMA specfic functionality is embedded inside this #if statment
 */
#ifdef PTATIN_HAVE_SPMA

int spmA_InitialiseTopo_pTatin3d(SPMAData *spm,int nx,int nz,double x[],double z[],double h[])
{
  int i,j,ii;
  double xpos,ypos;

  for (i=0; i<spm->nx; i++) {
    for (j=0; j<spm->ny; j++) {
      int nid;
      double h_interp,sep,smax;

      nid = spmA_NID(i,j,spm->nx);
      xpos = spm->x[nid];
      ypos = spm->y[nid];

      /* interpolate height */
      smax = 1.0e32;
      for (ii=0; ii<nx*nz; ii++) {

        sep =  (x[ii]-xpos)*(x[ii]-xpos);
        sep += (z[ii]-ypos)*(z[ii]-ypos);
        if (sep < smax) {
          h_interp = h[ii];
        }
      }

      spm->h_old[nid] = h_interp;
    }
  }

  return 1;
}

PetscErrorCode _ptatin3d_ApplyLandscapeEvolutionModel_SPMA(pTatinCtx pctx,Vec X)
{
  PhysCompStokes  stokes;
  DM              stokes_pack,dau,dap;
  PetscReal       dt_mechanical;
  PetscMPIInt     rank,spm_rank,nproc;
  SPMAData *spm;
  PetscReal min[3],max[3];
  double Lx,Lz;
  int      ie;
  double dt,dt_final;
  long int nx,nz;
  double *x,*y,*z;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = pTatinGetStokesContext(pctx,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);

  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dau),&nproc);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dau),&rank);CHKERRQ(ierr);
  spm_rank = nproc - 1;

  ierr= DMDAGetBoundingBox(dau,min,max);CHKERRQ(ierr);
  Lx = (double)( max[0] - min[0] );
  Lz = (double)( max[2] - min[2] );

  /* fetch the mechincal model surface on the desired core */
  ierr = ptatin3d_DMDAAllGatherCoorJMax(dau,spm_rank,&nx,&nz,&x,&y,&z);CHKERRQ(ierr);

  /* set time step for spm based on mechanical timestep */
  ierr = pTatinGetTimestep(pctx,&dt_mechanical);CHKERRQ(ierr);

  dt_final = dt_mechanical;
  dt       = dt_mechanical / 40.0;

  if (rank == spm_rank) {
    /* --- ---------------------- --- */
    /* --- call spm functionality --- */
    ie = spmA_New(&spm);CHKERRQ((PetscErrorCode)ie);
    ie = spmA_Initialise(spm,(int)2*nx+1,(int)2*nz+1,(double)Lx,(double)Lz,1.0,0.0,0.0,dt,dt_final);CHKERRQ((PetscErrorCode)ie);

    ie = spmA_InitialiseTopo_pTatin3d(spm,nx,nz,x,z,y);CHKERRQ((PetscErrorCode)ie);
    //ie = spmA_InitialiseUplift_pTatin3d(spm);
    spm->output_frequency = 100;

    ie = spmA_OutputIC(spm,"pt3d2spma.dat");CHKERRQ((PetscErrorCode)ie);
    /*
    ie = spmA_Apply(spm);CHKERRQ((PetscErrorCode)ie);
    ie = spmA_Output(spm,"test.dat");CHKERRQ((PetscErrorCode)ie);
    ie = spmA_Destroy(&spm);CHKERRQ((PetscErrorCode)ie);
    */
    /* --- ---------------------- --- */
  }

#if 0
  /* push update suface value into sequential dmda coordinate vector */
  ierr = DMDAVecGetArray(cda_surf,coor_surf,&LA_coor_surf);CHKERRQ(ierr);
  j = 0;
  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) {
      PetscInt idx;

      idx = i + k*nx;
      LA_coor_surf[k][j][i].y = z[idx];
    }
  }
  ierr = DMDAVecRestoreArray(cda_surf,coor_surf,&LA_coor_surf);CHKERRQ(ierr);
#endif

  /* scatter new surface back to parallel mesh */

  if (x) { PetscFree(x); }
  if (y) { PetscFree(y); }
  if (z) { PetscFree(z); }

  PetscFunctionReturn(0);
}

#endif
/* ----------------------------------------------------------------------------------------------- */

PetscErrorCode ptatin3d_ApplyLandscapeEvolutionModel_SPMA(pTatinCtx pctx,Vec X)
{
#ifdef PTATIN_HAVE_SPMA
  PetscErrorCode ierr;
#endif

  PetscFunctionBegin;

#ifdef PTATIN_HAVE_SPMA
  ierr = _ptatin3d_ApplyLandscapeEvolutionModel_SPMA(pctx,X);CHKERRQ(ierr);
#else
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"pTatind3D must be compiled with external package <SPMA: A simple FD landscape evolution model>");
#endif

  PetscFunctionReturn(0);
}

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
 **    filename:   output_material_points.c
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


#include "petsc.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPStokesPl_def.h"
#include "MPntPEnergy_def.h"

#include "ptatin_utils.h"
#include "dmda_duplicate.h"
#include "dmda_element_q2p1.h"
#include "data_bucket.h"
#include "output_paraview.h"
#include "element_type_Q2.h"
#include "element_utils_q2.h"
#include "element_utils_q1.h"
#include "quadrature.h"

#include "output_material_points.h"


const char *MaterialPointVariableName[] =  {
  "region",
  "viscosity",
  "density",
  "viscous_strain",
  "plastic_strain",
  "yield_indicator",
  "diffusivity",
  "energy_source",
  0
};

const char *MaterialPointVariableParaviewDataType[] =  {
  "Int32",
  "Float64",
  "Float64",
  "Float32",
  "Float32",
  "Int16",
  "Float64",
  "Float64",
  0
};

PetscErrorCode _write_float(FILE *vtk_fp,const PetscInt mx,const PetscInt my,const PetscInt mz,float LA_cell[])
{
  PetscInt i,j,k;

  PetscFunctionBegin;
  for (k=0; k<mz; k++) {
    for (j=0; j<my; j++) {
      for (i=0; i<mx; i++) {
        int idx = i + j*(mx) + k*(mx)*(my);

        fprintf( vtk_fp,"      %1.6e \n", LA_cell[idx]);
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode _write_double(FILE *vtk_fp,const PetscInt mx,const PetscInt my,const PetscInt mz,double LA_cell[])
{
  PetscInt i,j,k;

  PetscFunctionBegin;
  for (k=0; k<mz; k++) {
    for (j=0; j<my; j++) {
      for (i=0; i<mx; i++) {
        int idx = i + j*(mx) + k*(mx)*(my);

        fprintf( vtk_fp,"      %1.6e \n", LA_cell[idx]);
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode _write_int(FILE *vtk_fp,const PetscInt mx,const PetscInt my,const PetscInt mz,int LA_cell[])
{
  PetscInt i,j,k;

  PetscFunctionBegin;
  for (k=0; k<mz; k++) {
    for (j=0; j<my; j++) {
      for (i=0; i<mx; i++) {
        int idx = i + j*(mx) + k*(mx)*(my);

        fprintf( vtk_fp,"      %d \n", LA_cell[idx]);
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode _write_short(FILE *vtk_fp,const PetscInt mx,const PetscInt my,const PetscInt mz,short LA_cell[])
{
  PetscInt i,j,k;

  PetscFunctionBegin;
  for (k=0; k<mz; k++) {
    for (j=0; j<my; j++) {
      for (i=0; i<mx; i++) {
        int idx = i + j*(mx) + k*(mx)*(my);

        fprintf( vtk_fp,"      %d \n", LA_cell[idx]);
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode _check_for_empty_cells_double(const PetscInt mx,const PetscInt my,const PetscInt mz,int cell_count[],double LA_cell[])
{
  int uei,uej,uek;
  int e,ei,ej,ek,eid2,ii,jj,kk;
  int constant_conversion_occurred;
  PetscFunctionBegin;

  constant_conversion_occurred = 0;
  for (e=0; e<mx*my*mz; e++) {
    if (cell_count[e] == 0) {
      double local_LA_cell;
      int    local_cell_count;

      constant_conversion_occurred = 1;

      /* convert e into q2 eidx */
      ek   = e/(mx*my);
      eid2 = e - ek * (mx*my);
      ej   = eid2/mx;
      ei   = eid2 - ej * mx;

      uei = ei/2;
      uej = ej/2;
      uek = ek/2;

      /* traverse the q2 cell and try a new average */
      local_LA_cell = 0.0;
      local_cell_count = 0;
      for (kk=0; kk<2; kk++) {
        for (jj=0; jj<2; jj++) {
          for (ii=0; ii<2; ii++) {
            int cidx,ci,cj,ck;

            ci = 2*uei + ii;
            cj = 2*uej + jj;
            ck = 2*uek + kk;

            cidx = ci + cj*mx + ck*mx*my;
            if (cidx >= mx*my*mz) {
              SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"cidx is too large");
            }

            local_LA_cell += LA_cell[cidx];
            local_cell_count += cell_count[cidx];
          }
        }
      }
      /* set the same values on the 8 sub cells */
      for (kk=0; kk<2; kk++) {
        for (jj=0; jj<2; jj++) {
          for (ii=0; ii<2; ii++) {
            int cidx,ci,cj,ck;

            ci = 2*uei + ii;
            cj = 2*uej + jj;
            ck = 2*uek + kk;

            cidx = ci + cj*mx + ck*mx*my;
            LA_cell[cidx]    = local_LA_cell;
            cell_count[cidx] = local_cell_count;
          }
        }
      }


    }
  }

  if (constant_conversion_occurred == 1) {
    for (e=0; e<mx*my*mz; e++) {
      if (cell_count[e] == 0) {
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cell contains zero markers");
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode _check_for_empty_cells_float(const PetscInt mx,const PetscInt my,const PetscInt mz,int cell_count[],float LA_cell[])
{
  int uei,uej,uek;
  int e,ei,ej,ek,eid2,ii,jj,kk;
  int constant_conversion_occurred;
  PetscFunctionBegin;

  constant_conversion_occurred = 0;
  for (e=0; e<mx*my*mz; e++) {
    if (cell_count[e] == 0) {
      float local_LA_cell;
      int   local_cell_count;

      constant_conversion_occurred = 1;

      /* convert e into q2 eidx */
      ek   = e/(mx*my);
      eid2 = e - ek * (mx*my);
      ej   = eid2/mx;
      ei   = eid2 - ej * mx;

      uei = ei/2;
      uej = ej/2;
      uek = ek/2;

      /* traverse the q2 cell and try a new average */
      local_LA_cell = 0.0;
      local_cell_count = 0;
      for (kk=0; kk<2; kk++) {
        for (jj=0; jj<2; jj++) {
          for (ii=0; ii<2; ii++) {
            int cidx,ci,cj,ck;

            ci = 2*uei + ii;
            cj = 2*uej + jj;
            ck = 2*uek + kk;

            cidx = ci + cj*mx + ck*mx*my;
            if (cidx >= mx*my*mz) {
              SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"cidx is too large");
            }

            local_LA_cell += LA_cell[cidx];
            local_cell_count += cell_count[cidx];
          }
        }
      }
      /* set the same values on the 8 sub cells */
      for (kk=0; kk<2; kk++) {
        for (jj=0; jj<2; jj++) {
          for (ii=0; ii<2; ii++) {
            int cidx,ci,cj,ck;

            ci = 2*uei + ii;
            cj = 2*uej + jj;
            ck = 2*uek + kk;

            cidx = ci + cj*mx + ck*mx*my;
            LA_cell[cidx]    = local_LA_cell;
            cell_count[cidx] = local_cell_count;
          }
        }
      }


    }
  }

  if (constant_conversion_occurred == 1) {
    for (e=0; e<mx*my*mz; e++) {
      if (cell_count[e] == 0) {
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cell contains zero markers");
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode _compute_cell_value_double(DataBucket db,MaterialPointVariable variable,const PetscInt mx,const PetscInt my,const PetscInt mz,double LA_cell[])
{
  int *cell_count;
  int e,ueid,ueid2,umx,umy,uei,uej,uek;
  double *xi_p;
  int ei,ej,ek,eidx;
  double var;
  int p,n_mp;
  MPAccess X;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  umx = mx/2;
  umy = my/2;

  ierr = PetscMalloc(sizeof(int)*mx*my*mz,&cell_count);CHKERRQ(ierr);
  ierr = PetscMemzero(cell_count,sizeof(int)*mx*my*mz);CHKERRQ(ierr);

  DataBucketGetSizes(db,&n_mp,NULL,NULL);

  ierr = MaterialPointGetAccess(db,&X);CHKERRQ(ierr);

  for (p=0; p<n_mp; p++) {
    ierr = MaterialPointGet_local_element_index(X,p,&ueid);CHKERRQ(ierr);
    ierr = MaterialPointGet_local_coord(X,p,&xi_p);CHKERRQ(ierr);

    uek   = ueid/(umx*umy);
    ueid2 = ueid - uek * (umx*umy);
    uej   = ueid2/umx;
    uei   = ueid2 - uej * umx;

    ei = 2*uei;
    ej = 2*uej;
    ek = 2*uek;

    if (xi_p[0] > 0.0) { ei++; }
    if (xi_p[1] > 0.0) { ej++; }
    if (xi_p[2] > 0.0) { ek++; }

    eidx = ei + ej*mx + ek*mx*my;
    if (eidx >= mx*my*mz) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"eidx is too large");
    }

    switch (variable) {
      case MPV_viscosity:
        ierr = MaterialPointGet_viscosity(X,p,&var);CHKERRQ(ierr);
        break;
      case MPV_density:
        ierr = MaterialPointGet_density(X,p,&var);CHKERRQ(ierr);
        break;
      case MPV_diffusivity:
        ierr = MaterialPointGet_diffusivity(X,p,&var);CHKERRQ(ierr);
        break;
      case MPV_heat_source:
        ierr = MaterialPointGet_heat_source(X,p,&var);CHKERRQ(ierr);
        break;

      case MPV_region:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_region is not of type \"double\"");
        break;
	  case MPV_viscous_strain:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_viscous_strain is not of type \"double\"");
        break;
      case MPV_plastic_strain:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_plastic_strain is not of type \"double\"");
        break;
      case MPV_yield_indicator:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_yield_indicator is not of type \"double\"");
        break;

      default:
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown material point variable type (MPV_type)");
        break;
    }

    LA_cell[eidx] += var;
    cell_count[eidx]++;
  }

  ierr = _check_for_empty_cells_double(mx,my,mz,cell_count,LA_cell);CHKERRQ(ierr);

  for (e=0; e<mx*my*mz; e++) {
    LA_cell[e] = LA_cell[e] / ( (double)(cell_count[e]) );
  }

  ierr = MaterialPointRestoreAccess(db,&X);CHKERRQ(ierr);
  ierr = PetscFree(cell_count);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode _compute_cell_value_float(DataBucket db,MaterialPointVariable variable,const PetscInt mx,const PetscInt my,const PetscInt mz,float LA_cell[])
{
  int *cell_count;
  int e,ueid,ueid2,umx,umy,uei,uej,uek;
  double *xi_p;
  int ei,ej,ek,eidx;
  float var;
  int p,n_mp;
  MPAccess X;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  umx = mx/2;
  umy = my/2;

  ierr = PetscMalloc(sizeof(int)*mx*my*mz,&cell_count);CHKERRQ(ierr);
  ierr = PetscMemzero(cell_count,sizeof(int)*mx*my*mz);CHKERRQ(ierr);

  DataBucketGetSizes(db,&n_mp,NULL,NULL);

  ierr = MaterialPointGetAccess(db,&X);CHKERRQ(ierr);

  for (p=0; p<n_mp; p++) {
    ierr = MaterialPointGet_local_element_index(X,p,&ueid);CHKERRQ(ierr);
    ierr = MaterialPointGet_local_coord(X,p,&xi_p);CHKERRQ(ierr);

    uek   = ueid/(umx*umy);
    ueid2 = ueid - uek * (umx*umy);
    uej   = ueid2/umx;
    uei   = ueid2 - uej * umx;

    ei = 2*uei;
    ej = 2*uej;
    ek = 2*uek;

    if (xi_p[0] > 0.0) { ei++; }
    if (xi_p[1] > 0.0) { ej++; }
    if (xi_p[2] > 0.0) { ek++; }

    eidx = ei + ej*mx + ek*mx*my;
    if (eidx >= mx*my*mz) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"eidx is too large");
    }

    switch (variable) {
      case MPV_plastic_strain:
        ierr = MaterialPointGet_plastic_strain(X,p,&var);CHKERRQ(ierr);
        break;
	  case MPV_viscous_strain:
        ierr = MaterialPointGet_viscous_strain(X,p,&var);CHKERRQ(ierr);
        break;

      case MPV_viscosity:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_viscosity is not of type \"float\"");
        break;
      case MPV_density:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_density is not of type \"float\"");
        break;
      case MPV_diffusivity:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_diffusivity is not of type \"float\"");
        break;
      case MPV_heat_source:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_heat_source is not of type \"float\"");
        break;

      case MPV_region:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_region is not of type \"float\"");
        break;
      case MPV_yield_indicator:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_yield_indicator is not of type \"float\"");
        break;

      default:
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown material point variable type (MPV_type)");
        break;
    }

    LA_cell[eidx] += var;
    cell_count[eidx]++;
  }

  ierr = _check_for_empty_cells_float(mx,my,mz,cell_count,LA_cell);CHKERRQ(ierr);

  for (e=0; e<mx*my*mz; e++) {
    LA_cell[e] = LA_cell[e] / ( (float)(cell_count[e]) );
  }

  ierr = MaterialPointRestoreAccess(db,&X);CHKERRQ(ierr);
  ierr = PetscFree(cell_count);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode _compute_cell_composition(DM dau,PetscScalar LA_gcoords[],DataBucket db,const PetscInt mx,const PetscInt my,const PetscInt mz,int LA_cell[])
{
  int *closest_point;
  int e,ueid,ueid2,umx,umy,uei,uej,uek,i,j,k,li,lj,lk,ii,jj,kk;
  double *xi_p,*x_p;
  int ei,ej,ek,eidx,idx;
  int p,n_mp;
  MPAccess X;
  PetscInt nel,nen;
  const PetscInt *elnidx;
  PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
  long int empty,gempty;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  umx = mx/2;
  umy = my/2;

  ierr = PetscMalloc(sizeof(int)*mx*my*mz,&closest_point);CHKERRQ(ierr);
  for (e=0; e<mx*my*mz; e++) {
    //LA_cell[e]       = -1;
    closest_point[e] = -1;
  }

  ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen,&elnidx);CHKERRQ(ierr);


  DataBucketGetSizes(db,&n_mp,NULL,NULL);
  ierr = MaterialPointGetAccess(db,&X);CHKERRQ(ierr);

  for (p=0; p<n_mp; p++) {
    ierr = MaterialPointGet_local_element_index(X,p,&ueid);CHKERRQ(ierr);
    ierr = MaterialPointGet_local_coord(X,p,&xi_p);CHKERRQ(ierr);
    ierr = MaterialPointGet_global_coord(X,p,&x_p);CHKERRQ(ierr);
    uek   = ueid/(umx*umy);

    ueid2 = ueid - uek * (umx*umy);
    uej   = ueid2/umx;
    uei   = ueid2 - uej * umx;

    ei = 2*uei;
    ej = 2*uej;
    ek = 2*uek;

    li = lj = lk = 0;
    if (xi_p[0] > 0.0) { ei++; li++; }
    if (xi_p[1] > 0.0) { ej++; lj++; }
    if (xi_p[2] > 0.0) { ek++; lk++; }

    eidx = ei + ej*mx + ek*mx*my;
    if (eidx >= mx*my*mz) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"eidx is too large");
    }


    if (closest_point[eidx] == -1) {
      /* cell unclaimed */
      closest_point[eidx] = p;
    } else {
      /* someone else claimed it */
      double *x_p_closest;
      int p_closest;
      double x_center[3];
      double xc[3*8];
      double dist_p_closest,dist_p;

      ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx[nen*ueid],LA_gcoords);CHKERRQ(ierr);
      for (k=lk; k<lk+2; k++) {
        for (j=lj; j<lj+2; j++) {
          for (i=li; i<li+2; i++) {
            int idx,gidx;

            gidx = i + j*3 + k*3*3;
            idx  = (i-li) + (j-lj)*2 + (k-lk)*2*2;

            xc[3*idx+0] = elcoords[3*gidx+0];
            xc[3*idx+1] = elcoords[3*gidx+1];
            xc[3*idx+2] = elcoords[3*gidx+2];
          }
        }
      }

      x_center[0] = x_center[1] = x_center[2] = 0.0;
      for (i=0; i<8; i++) {
        x_center[0] += 0.125 * xc[3*i+0];
        x_center[1] += 0.125 * xc[3*i+1];
        x_center[2] += 0.125 * xc[3*i+2];
      }

      p_closest = closest_point[eidx];
      ierr = MaterialPointGet_global_coord(X,p_closest,&x_p_closest);CHKERRQ(ierr);

      dist_p_closest = (x_center[0]-x_p_closest[0])*(x_center[0]-x_p_closest[0])
                     + (x_center[1]-x_p_closest[1])*(x_center[1]-x_p_closest[1])
                     + (x_center[2]-x_p_closest[2])*(x_center[2]-x_p_closest[2]);
      dist_p         = (x_center[0]-x_p[0])*(x_center[0]-x_p[0])
                     + (x_center[1]-x_p[1])*(x_center[1]-x_p[1])
                     + (x_center[2]-x_p[2])*(x_center[2]-x_p[2]);

      if (dist_p < dist_p_closest) {
        closest_point[eidx] = p;
      }
    }
  }

  /* check cells if empty */
  for (ek=0; ek<mz; ek++) {
    for (ej=0; ej<my; ej++) {
      for (ei=0; ei<mx; ei++) {
        int eidx;

        eidx = ei + ej*mx + ek*mx*my;

        if (closest_point[eidx] != -1) { continue; }

        //printf("** cell %d is empty **\n", eidx);

        /* else take action */
        {
          double x_center[3];
          double xc[3*8];
          double min_sep;
          int p_closest;
          double *x_p_closest;
          double dist_p_closest;

          /* compute cell centroid */
          uei = ei/2;
          uej = ej/2;
          uek = ek/2;
          ueid = uei + uej*(umx) + uek*(umx)*(umy);
          //printf("i,j,k %d %d %d : uei %d %d %d : ueid = %d : %d \n", ei,ej,ek,uei,uej,uej,ueid,umx*umy*umz);

          li = ei - 2*uei;
          lj = ej - 2*uej;
          lk = ek - 2*uek;
          //printf("  li = %d %d %d \n", li,lj,lk);

          ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx[nen*ueid],LA_gcoords);CHKERRQ(ierr);
          for (kk=lk; kk<lk+2; kk++) {
            for (jj=lj; jj<lj+2; jj++) {
              for (ii=li; ii<li+2; ii++) {
                int idx,gidx;

                gidx = ii + jj*3 + kk*3*3;
                idx  = (ii-li) + (jj-lj)*2 + (kk-lk)*2*2;

                xc[3*idx+0] = elcoords[3*gidx+0];
                xc[3*idx+1] = elcoords[3*gidx+1];
                xc[3*idx+2] = elcoords[3*gidx+2];
              }
            }
          }

          x_center[0] = x_center[1] = x_center[2] = 0.0;
          for (ii=0; ii<8; ii++) {
            x_center[0] += 0.125 * xc[3*ii+0];
            x_center[1] += 0.125 * xc[3*ii+1];
            x_center[2] += 0.125 * xc[3*ii+2];
          }
          //printf("  %1.4e %1.4e %1.4e \n", x_center[0],x_center[1],x_center[2]);

          /* scan neighbours */
          min_sep = 1.0e32;
          for (k=ek-2; k<=ek+2; k++) {
            if (k < 0)  { continue; }
            if (k >= mz) { continue; }

            for (j=ej-2; j<=ej+2; j++) {
              if (j < 0)  { continue; }
              if (j >= my) { continue; }

              for (i=ei-2; i<=ei+2; i++) {
                if (i < 0)  { continue; }
                if (i >= mx) { continue; }

                idx = i + j*mx + k*mx*my;
                //printf("  checking cell %d \n", idx );
                if ( idx >= mx*my*mz ) {
                  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"idx too large");
                }

                p_closest = closest_point[ idx ];
                if (p_closest == -1) { continue; } /* skip if neighbour cell is also empty */
                //printf("  found neighbour with owner\n");

                ierr = MaterialPointGet_global_coord(X,p_closest,&x_p_closest);CHKERRQ(ierr);

                /* compute dist from p1 to x_center */
                dist_p_closest = (x_center[0]-x_p_closest[0])*(x_center[0]-x_p_closest[0])
                               + (x_center[1]-x_p_closest[1])*(x_center[1]-x_p_closest[1])
                               + (x_center[2]-x_p_closest[2])*(x_center[2]-x_p_closest[2]);

                if (dist_p_closest < min_sep) {
                  min_sep = dist_p_closest;
                  closest_point[ eidx ] = p_closest;
                  //printf(" asssigning cell %d the value from cell %d which associated with marker %d \n", eidx,idx,p_closest);
                }

              }
            }
          }


        }

      }
    }
  }
  /* report cells if empty */
  empty = 0;
  for (e=0; e<mx*my*mz; e++) {
    if (closest_point[e] == -1) {
      empty++;
    }
  }
  ierr = MPI_Allreduce(&empty,&gempty,1,MPI_LONG,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  if (gempty != 0) {
    PetscPrintf(PETSC_COMM_WORLD,"WARNING(_compute_cell_composition): Detected %ld cells which could not assigned a composition\n",gempty);
  }

  /* assign phase based on nearest point */
  for (e=0; e<mx*my*mz; e++) {
    int pid,phase;

    pid = closest_point[e];
    if (pid == -1) {
      LA_cell[e] = -1;
    } else {
      ierr = MaterialPointGet_phase_index(X,pid,&phase);CHKERRQ(ierr);
      LA_cell[e] = phase;
    }
  }

  ierr = MaterialPointRestoreAccess(db,&X);CHKERRQ(ierr);
  ierr = PetscFree(closest_point);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

PetscErrorCode _compute_cell_nearest_point(DM dau,PetscScalar LA_gcoords[],DataBucket db,const PetscInt mx,const PetscInt my,const PetscInt mz,int closest_point[],PetscBool *empty_cells_detected)
{
  int e,ueid,ueid2,umx,umy,uei,uej,uek,i,j,k,li,lj,lk,ii,jj,kk;
  double *xi_p,*x_p;
  int ei,ej,ek,eidx,idx;
  int p,n_mp;
  MPAccess X;
  PetscInt nel,nen;
  const PetscInt *elnidx;
  PetscReal elcoords[3*Q2_NODES_PER_EL_3D];
  long int empty,gempty;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  umx = mx/2;
  umy = my/2;

  for (e=0; e<mx*my*mz; e++) {
    closest_point[e] = -1;
  }

  ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen,&elnidx);CHKERRQ(ierr);

  DataBucketGetSizes(db,&n_mp,NULL,NULL);
  ierr = MaterialPointGetAccess(db,&X);CHKERRQ(ierr);

  for (p=0; p<n_mp; p++) {
    ierr = MaterialPointGet_local_element_index(X,p,&ueid);CHKERRQ(ierr);
    ierr = MaterialPointGet_local_coord(X,p,&xi_p);CHKERRQ(ierr);
    ierr = MaterialPointGet_global_coord(X,p,&x_p);CHKERRQ(ierr);
    uek   = ueid/(umx*umy);

    ueid2 = ueid - uek * (umx*umy);
    uej   = ueid2/umx;
    uei   = ueid2 - uej * umx;

    ei = 2*uei;
    ej = 2*uej;
    ek = 2*uek;

    li = lj = lk = 0;
    if (xi_p[0] > 0.0) { ei++; li++; }
    if (xi_p[1] > 0.0) { ej++; lj++; }
    if (xi_p[2] > 0.0) { ek++; lk++; }

    eidx = ei + ej*mx + ek*mx*my;
    if (eidx >= mx*my*mz) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"eidx is too large");
    }


    if (closest_point[eidx] == -1) {
      /* cell unclaimed */
      closest_point[eidx] = p;
    } else {
      /* someone else claimed it */
      double *x_p_closest;
      int p_closest;
      double x_center[3];
      double xc[3*8];
      double dist_p_closest,dist_p;

      ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx[nen*ueid],LA_gcoords);CHKERRQ(ierr);
      for (k=lk; k<lk+2; k++) {
        for (j=lj; j<lj+2; j++) {
          for (i=li; i<li+2; i++) {
            int idx,gidx;

            gidx = i + j*3 + k*3*3;
            idx  = (i-li) + (j-lj)*2 + (k-lk)*2*2;

            xc[3*idx+0] = elcoords[3*gidx+0];
            xc[3*idx+1] = elcoords[3*gidx+1];
            xc[3*idx+2] = elcoords[3*gidx+2];
          }
        }
      }

      x_center[0] = x_center[1] = x_center[2] = 0.0;
      for (i=0; i<8; i++) {
        x_center[0] += 0.125 * xc[3*i+0];
        x_center[1] += 0.125 * xc[3*i+1];
        x_center[2] += 0.125 * xc[3*i+2];
      }

      p_closest = closest_point[eidx];
      ierr = MaterialPointGet_global_coord(X,p_closest,&x_p_closest);CHKERRQ(ierr);

      dist_p_closest = (x_center[0]-x_p_closest[0])*(x_center[0]-x_p_closest[0])
      + (x_center[1]-x_p_closest[1])*(x_center[1]-x_p_closest[1])
      + (x_center[2]-x_p_closest[2])*(x_center[2]-x_p_closest[2]);
      dist_p         = (x_center[0]-x_p[0])*(x_center[0]-x_p[0])
      + (x_center[1]-x_p[1])*(x_center[1]-x_p[1])
      + (x_center[2]-x_p[2])*(x_center[2]-x_p[2]);

      if (dist_p < dist_p_closest) {
        closest_point[eidx] = p;
      }
    }
  }

  /* check cells if empty */
  for (ek=0; ek<mz; ek++) {
    for (ej=0; ej<my; ej++) {
      for (ei=0; ei<mx; ei++) {
        int eidx;

        eidx = ei + ej*mx + ek*mx*my;

        if (closest_point[eidx] != -1) { continue; }

        /* else take action */
        {
          double x_center[3];
          double xc[3*8];
          double min_sep;
          int p_closest;
          double *x_p_closest;
          double dist_p_closest;

          /* compute cell centroid */
          uei = ei/2;
          uej = ej/2;
          uek = ek/2;
          ueid = uei + uej*(umx) + uek*(umx)*(umy);

          li = ei - 2*uei;
          lj = ej - 2*uej;
          lk = ek - 2*uek;

          ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx[nen*ueid],LA_gcoords);CHKERRQ(ierr);
          for (kk=lk; kk<lk+2; kk++) {
            for (jj=lj; jj<lj+2; jj++) {
              for (ii=li; ii<li+2; ii++) {
                int idx,gidx;

                gidx = ii + jj*3 + kk*3*3;
                idx  = (ii-li) + (jj-lj)*2 + (kk-lk)*2*2;

                xc[3*idx+0] = elcoords[3*gidx+0];
                xc[3*idx+1] = elcoords[3*gidx+1];
                xc[3*idx+2] = elcoords[3*gidx+2];
              }
            }
          }

          x_center[0] = x_center[1] = x_center[2] = 0.0;
          for (ii=0; ii<8; ii++) {
            x_center[0] += 0.125 * xc[3*ii+0];
            x_center[1] += 0.125 * xc[3*ii+1];
            x_center[2] += 0.125 * xc[3*ii+2];
          }

          /* scan neighbours */
          min_sep = 1.0e32;
          for (k=ek-2; k<=ek+2; k++) {
            if (k < 0)  { continue; }
            if (k >= mz) { continue; }

            for (j=ej-2; j<=ej+2; j++) {
              if (j < 0)  { continue; }
              if (j >= my) { continue; }

              for (i=ei-2; i<=ei+2; i++) {
                if (i < 0)  { continue; }
                if (i >= mx) { continue; }

                idx = i + j*mx + k*mx*my;
                if ( idx >= mx*my*mz ) {
                  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"idx too large");
                }

                p_closest = closest_point[ idx ];
                if (p_closest == -1) { continue; } /* skip if neighbour cell is also empty */

                ierr = MaterialPointGet_global_coord(X,p_closest,&x_p_closest);CHKERRQ(ierr);

                /* compute dist from p1 to x_center */
                dist_p_closest = (x_center[0]-x_p_closest[0])*(x_center[0]-x_p_closest[0])
                + (x_center[1]-x_p_closest[1])*(x_center[1]-x_p_closest[1])
                + (x_center[2]-x_p_closest[2])*(x_center[2]-x_p_closest[2]);

                if (dist_p_closest < min_sep) {
                  min_sep = dist_p_closest;
                  closest_point[ eidx ] = p_closest;
                }

              }
            }
          } /* end neighbour cell k */

        }

      } /* end cell k */
    }
  }
  /* report cells if empty */
  *empty_cells_detected = PETSC_FALSE;
  empty = 0;
  for (e=0; e<mx*my*mz; e++) {
    if (closest_point[e] == -1) {
      empty++;
    }
  }
  ierr = MPI_Allreduce(&empty,&gempty,1,MPI_LONG,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  if (gempty != 0) {
    PetscPrintf(PETSC_COMM_WORLD,"WARNING(_compute_cell_nearest_point): Detected %ld cells which could not assigned a composition\n",gempty);
    *empty_cells_detected = PETSC_TRUE;
  }

  ierr = MaterialPointRestoreAccess(db,&X);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode _compute_cell_value_short(DataBucket db,MaterialPointVariable variable,const PetscInt ncells,const int closest_point[],short LA_cell[])
{
  int e,pid,n_mp;
  MPAccess X;
  PetscErrorCode ierr;

  DataBucketGetSizes(db,&n_mp,NULL,NULL);
  ierr = MaterialPointGetAccess(db,&X);CHKERRQ(ierr);
  for (e=0; e<ncells; e++) {
    LA_cell[e] = -1;
  }

  switch (variable) {
    case MPV_region:
    {
      int var_i;

      for (e=0; e<ncells; e++) {
        pid = closest_point[e];
        if (pid != -1) {
          ierr = MaterialPointGet_phase_index(X,pid,&var_i);CHKERRQ(ierr);
          LA_cell[e] = (short)var_i;
        }
      }
    }
      break;

    case MPV_yield_indicator:
    {
      short var_s;

      for (e=0; e<ncells; e++) {
        pid = closest_point[e];
        if (pid != -1) {
          ierr = MaterialPointGet_yield_indicator(X,pid,&var_s);CHKERRQ(ierr);
          LA_cell[e] = var_s;
        }
      }
    }
      break;

    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot convert material point variable type %D to short - consult typedef enum {} MaterialPointVariable in output_material_point.h to understand what the quantity is and whether _compute_cell_value_short() should be updated",(PetscInt)variable);
      break;
  }
  ierr = MaterialPointRestoreAccess(db,&X);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinOutputParaViewMarkerFields_VTS(DM dau,DataBucket material_points,const int nvars,const MaterialPointVariable vars[],const char name[])
{
  PetscErrorCode ierr;
  DM cda;
  Vec gcoords;
  DMDACoor3d ***LA_gcoords;
  PetscScalar *LA_gc;
  PetscInt mx,my,mz;
  PetscInt i,j,k,esi,esj,esk;
  FILE* vtk_fp = NULL;
  PetscInt gsi,gsj,gsk,gm,gn,gp;
  int t;
  int    *i_LA_cell;
  short  *s_LA_cell;
  float  *f_LA_cell;
  double *d_LA_cell;

  PetscFunctionBegin;
  if ((vtk_fp = fopen ( name, "w")) == NULL)  {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
  }


  ierr = DMDAGetGhostCorners(dau,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
  ierr = DMDAGetCornersElementQ2(dau,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);

  ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);


  /* VTS HEADER - OPEN */
  fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
  PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "    <Piece Extent=\"%D %D %D %D %D %D\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);

  /* VTS COORD DATA */
  fprintf( vtk_fp, "    <Points>\n");
  fprintf( vtk_fp, "      <DataArray Name=\"coords\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (k=esk; k<esk+2*mz+1; k++) {
    for (j=esj; j<esj+2*my+1; j++) {
      for (i=esi; i<esi+2*mx+1; i++) {
        float xc,yc,zc;

        xc = (float)LA_gcoords[k][j][i].x; if (fabsf(xc) < 1.0e-12) { xc = 0.0; }
        yc = (float)LA_gcoords[k][j][i].y; if (fabsf(yc) < 1.0e-12) { yc = 0.0; }
        zc = (float)LA_gcoords[k][j][i].z; if (fabsf(zc) < 1.0e-12) { zc = 0.0; }


        fprintf( vtk_fp,"      %1.6e %1.6e %1.6e\n", xc,yc,zc );
      }
    }
  }
  fprintf( vtk_fp, "      </DataArray>\n");
  fprintf( vtk_fp, "    </Points>\n");
  ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);


  /* VTS CELL DATA */
  fprintf( vtk_fp, "    <CellData>\n");

  if (nvars == -1) {
    const char *mpv_name;

    t = 0;
    mpv_name = MaterialPointVariableName[t];
    while (mpv_name != NULL) {

      fprintf( vtk_fp, "      <DataArray Name=\"%s\" type=\"%s\" NumberOfComponents=\"1\" format=\"ascii\">\n",MaterialPointVariableName[t],MaterialPointVariableParaviewDataType[t]);

      switch (t) {

        case MPV_region:
          ierr = PetscMalloc(sizeof(int)*2*mx*2*my*2*mz,&i_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(i_LA_cell,sizeof(int)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = VecGetArray(gcoords,&LA_gc);CHKERRQ(ierr);
          ierr = _compute_cell_composition(dau,LA_gc,material_points,2*mx,2*my,2*mz,i_LA_cell);CHKERRQ(ierr);
          ierr = VecRestoreArray(gcoords,&LA_gc);CHKERRQ(ierr);

          ierr = _write_int(vtk_fp,2*mx,2*my,2*mz,i_LA_cell);CHKERRQ(ierr);


          ierr = PetscFree(i_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_viscosity:
          ierr = PetscMalloc(sizeof(double)*2*mx*2*my*2*mz,&d_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(d_LA_cell,sizeof(double)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_double(material_points,MPV_viscosity,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);
          ierr = _write_double(vtk_fp,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(d_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_density:
          ierr = PetscMalloc(sizeof(double)*2*mx*2*my*2*mz,&d_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(d_LA_cell,sizeof(double)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_double(material_points,MPV_density,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);
          ierr = _write_double(vtk_fp,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(d_LA_cell);CHKERRQ(ierr);
          break;
		  
		case MPV_viscous_strain:
          ierr = PetscMalloc(sizeof(float)*2*mx*2*my*2*mz,&f_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(f_LA_cell,sizeof(float)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_float(material_points,MPV_viscous_strain,2*mx,2*my,2*mz,f_LA_cell);CHKERRQ(ierr);
          ierr = _write_float(vtk_fp,2*mx,2*my,2*mz,f_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(f_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_plastic_strain:
          ierr = PetscMalloc(sizeof(float)*2*mx*2*my*2*mz,&f_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(f_LA_cell,sizeof(float)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_float(material_points,MPV_plastic_strain,2*mx,2*my,2*mz,f_LA_cell);CHKERRQ(ierr);
          ierr = _write_float(vtk_fp,2*mx,2*my,2*mz,f_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(f_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_yield_indicator:
        {
          int *closest_point;
          PetscBool empty_cells_detected;

          ierr = PetscMalloc1(2*mx*2*my*2*mz,&closest_point);CHKERRQ(ierr);

          ierr = PetscMalloc1(2*mx*2*my*2*mz,&s_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(s_LA_cell,sizeof(short)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = VecGetArray(gcoords,&LA_gc);CHKERRQ(ierr);
          ierr = _compute_cell_nearest_point(dau,LA_gc,material_points,2*mx,2*my,2*mz,closest_point,&empty_cells_detected);CHKERRQ(ierr);
          ierr = VecRestoreArray(gcoords,&LA_gc);CHKERRQ(ierr);

          ierr = _compute_cell_value_short(material_points,MPV_yield_indicator,2*mx*2*my*2*mz,(const int*)closest_point,s_LA_cell);CHKERRQ(ierr);

          ierr = _write_short(vtk_fp,2*mx,2*my,2*mz,s_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(s_LA_cell);CHKERRQ(ierr);
          ierr = PetscFree(closest_point);CHKERRQ(ierr);
        }
          break;

        case MPV_diffusivity:
          ierr = PetscMalloc(sizeof(double)*2*mx*2*my*2*mz,&d_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(d_LA_cell,sizeof(double)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_double(material_points,MPV_diffusivity,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);
          ierr = _write_double(vtk_fp,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(d_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_heat_source:
          ierr = PetscMalloc(sizeof(double)*2*mx*2*my*2*mz,&d_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(d_LA_cell,sizeof(double)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_double(material_points,MPV_heat_source,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);
          ierr = _write_double(vtk_fp,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(d_LA_cell);CHKERRQ(ierr);
          break;

                default:
                    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown material point variable type (MPV_type)");
                    break;
      }

      fprintf( vtk_fp, "      </DataArray>\n");

      t++;
      mpv_name = MaterialPointVariableName[t];
    }
  } else {
    for (t=0; t<nvars; t++) {
      MaterialPointVariable idx = vars[t];

      fprintf( vtk_fp, "      <DataArray Name=\"%s\" type=\"%s\" NumberOfComponents=\"1\" format=\"ascii\">\n",MaterialPointVariableName[idx],MaterialPointVariableParaviewDataType[idx]);

      switch (idx) {

        case MPV_region:
          ierr = PetscMalloc(sizeof(int)*2*mx*2*my*2*mz,&i_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(i_LA_cell,sizeof(int)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = VecGetArray(gcoords,&LA_gc);CHKERRQ(ierr);
          ierr = _compute_cell_composition(dau,LA_gc,material_points,2*mx,2*my,2*mz,i_LA_cell);CHKERRQ(ierr);
          ierr = VecRestoreArray(gcoords,&LA_gc);CHKERRQ(ierr);

          ierr = _write_int(vtk_fp,2*mx,2*my,2*mz,i_LA_cell);CHKERRQ(ierr);


          ierr = PetscFree(i_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_viscosity:
          ierr = PetscMalloc(sizeof(double)*2*mx*2*my*2*mz,&d_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(d_LA_cell,sizeof(double)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_double(material_points,MPV_viscosity,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);
          ierr = _write_double(vtk_fp,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(d_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_density:
          ierr = PetscMalloc(sizeof(double)*2*mx*2*my*2*mz,&d_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(d_LA_cell,sizeof(double)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_double(material_points,MPV_density,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);
          ierr = _write_double(vtk_fp,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(d_LA_cell);CHKERRQ(ierr);
          break;
		  
		case MPV_viscous_strain:
          ierr = PetscMalloc(sizeof(float)*2*mx*2*my*2*mz,&f_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(f_LA_cell,sizeof(float)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_float(material_points,MPV_viscous_strain,2*mx,2*my,2*mz,f_LA_cell);CHKERRQ(ierr);
          ierr = _write_float(vtk_fp,2*mx,2*my,2*mz,f_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(f_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_plastic_strain:
          ierr = PetscMalloc(sizeof(float)*2*mx*2*my*2*mz,&f_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(f_LA_cell,sizeof(float)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_float(material_points,MPV_plastic_strain,2*mx,2*my,2*mz,f_LA_cell);CHKERRQ(ierr);
          ierr = _write_float(vtk_fp,2*mx,2*my,2*mz,f_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(f_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_yield_indicator:
          ierr = PetscMalloc(sizeof(short)*2*mx*2*my*2*mz,&s_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(s_LA_cell,sizeof(short)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          PetscPrintf(PETSC_COMM_WORLD,"MPV_yield_indicator -> writer not yet completed\n");
          ierr = _write_short(vtk_fp,2*mx,2*my,2*mz,s_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(s_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_diffusivity:
          ierr = PetscMalloc(sizeof(double)*2*mx*2*my*2*mz,&d_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(d_LA_cell,sizeof(double)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_double(material_points,MPV_diffusivity,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);
          ierr = _write_double(vtk_fp,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(d_LA_cell);CHKERRQ(ierr);
          break;

        case MPV_heat_source:
          ierr = PetscMalloc(sizeof(double)*2*mx*2*my*2*mz,&d_LA_cell);CHKERRQ(ierr);
          ierr = PetscMemzero(d_LA_cell,sizeof(double)*2*mx*2*my*2*mz);CHKERRQ(ierr);

          ierr = _compute_cell_value_double(material_points,MPV_heat_source,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);
          ierr = _write_double(vtk_fp,2*mx,2*my,2*mz,d_LA_cell);CHKERRQ(ierr);

          ierr = PetscFree(d_LA_cell);CHKERRQ(ierr);
          break;

                default:
                    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown material point variable type (MPV_type)");
                    break;
      }


      fprintf( vtk_fp, "      </DataArray>\n");
    }
  }


  fprintf( vtk_fp, "    </CellData>\n");

  /* VTS NODAL DATA */
  fprintf( vtk_fp, "    <PointData>\n");
  fprintf( vtk_fp, "    </PointData>\n");

  /* VTS HEADER - CLOSE */
  fprintf( vtk_fp, "    </Piece>\n");
  fprintf( vtk_fp, "  </StructuredGrid>\n");
  fprintf( vtk_fp, "</VTKFile>\n");


  fclose( vtk_fp );

  PetscFunctionReturn(0);
}


PetscErrorCode pTatinOutputParaViewMarkerFields_PVTS(DM dau,const int nvars,const MaterialPointVariable vars[],const char prefix[],const char name[])
{
  PetscErrorCode ierr;
  FILE* vtk_fp = NULL;
  PetscInt M,N,P,swidth;
  PetscMPIInt rank;
  int t;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  vtk_fp = NULL;
  if (rank==0) {
    if ((vtk_fp = fopen ( name, "w")) == NULL)  {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
    }
  }


  /* VTS HEADER - OPEN */
  if(vtk_fp) fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
  if(vtk_fp) fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

  DMDAGetInfo( dau, 0, &M,&N,&P, 0,0,0, 0,&swidth, 0,0,0, 0 );
  if(vtk_fp) PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <PStructuredGrid GhostLevel=\"%D\" WholeExtent=\"%D %D %D %D %D %D\">\n", swidth, 0,M-1, 0,N-1, 0,P-1 ); /* note overlap = 1 for Q1 */

  /* VTS COORD DATA */
  if(vtk_fp) fprintf( vtk_fp, "    <PPoints>\n");
  if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"Float32\" Name=\"coords\" NumberOfComponents=\"3\"/>\n");
  if(vtk_fp) fprintf( vtk_fp, "    </PPoints>\n");


  /* VTS CELL DATA */
  if(vtk_fp) fprintf( vtk_fp, "    <PCellData>\n");
  if (nvars == -1) {
    const char *mpv_name;

    t = 0;
    mpv_name = MaterialPointVariableName[t];
    while (mpv_name != NULL) {
      if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",MaterialPointVariableParaviewDataType[t],MaterialPointVariableName[t]);
      t++;
      mpv_name = MaterialPointVariableName[t];
    }
  } else {
    for (t=0; t<nvars; t++) {
      MaterialPointVariable idx = vars[t];

      if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",MaterialPointVariableParaviewDataType[idx],MaterialPointVariableName[idx]);
    }
  }

  if(vtk_fp) fprintf( vtk_fp, "    </PCellData>\n");

  /* VTS NODAL DATA */
  if(vtk_fp) fprintf( vtk_fp, "    <PPointData>\n");
  if(vtk_fp) fprintf( vtk_fp, "    </PPointData>\n");

  /* write out the parallel information */
  ierr = DAQ2PieceExtendForGhostLevelZero(vtk_fp,2,dau,prefix);CHKERRQ(ierr);

  /* VTS HEADER - CLOSE */
  if(vtk_fp) fprintf( vtk_fp, "  </PStructuredGrid>\n");
  if(vtk_fp) fprintf( vtk_fp, "</VTKFile>\n");

  if(vtk_fp) fclose( vtk_fp );
  PetscFunctionReturn(0);
}


PetscErrorCode pTatinOutputParaViewMarkerFields(DM pack,DataBucket material_points,const int nvars,const MaterialPointVariable vars[],const char path[],const char prefix[])
{
  char *vtkfilename,*filename;
  PetscMPIInt rank;
  DM dau,dap;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  ierr = pTatinGenerateParallelVTKName(prefix,"vts",&vtkfilename);CHKERRQ(ierr);
  if (path) {
    if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
  } else {
    if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
  }

  ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);

  ierr = pTatinOutputParaViewMarkerFields_VTS(dau,material_points,nvars,vars,filename);CHKERRQ(ierr);

  free(filename);
  free(vtkfilename);

  ierr = pTatinGenerateVTKName(prefix,"pvts",&vtkfilename);CHKERRQ(ierr);
  if (path) {
    if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
  } else {
    if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
  }

  ierr = pTatinOutputParaViewMarkerFields_PVTS(dau,nvars,vars,prefix,filename);CHKERRQ(ierr);

  free(filename);
  free(vtkfilename);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatin3d_ModelOutput_MarkerCellFields(pTatinCtx ctx,const int nvars,const MaterialPointVariable vars[],const char prefix[])
{
  PetscErrorCode ierr;
  DM             stokes_pack;
  PetscLogDouble t0,t1;
  static PetscBool beenhere=PETSC_FALSE;
  DataBucket     material_points;
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
  PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_mpoints_cell.pvd",root);
  if (prefix) { PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "%s_mpoints_cell.pvts",prefix);
  } else {      PetscSNPrintf(vtkfilename, PETSC_MAX_PATH_LEN-1, "mpoints_cell.pvts");           }

  if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n", pvdfilename ); }
  PetscSNPrintf(stepprefix,PETSC_MAX_PATH_LEN-1,"step%D",ctx->step);
  ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time, vtkfilename, stepprefix);CHKERRQ(ierr);
  beenhere = PETSC_TRUE;

  // PVTS + VTS
  if (prefix) { PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1,"%s_mpoints_cell",prefix);
  } else {      PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1,"mpoints_cell",prefix);    }

  ierr = pTatinGetMaterialPoints(ctx,&material_points,NULL);CHKERRQ(ierr);
  stokes_pack = ctx->stokes_ctx->stokes_pack;
  ierr = pTatinOutputParaViewMarkerFields(stokes_pack,material_points,nvars,vars,pvoutputdir,name);CHKERRQ(ierr);

  PetscTime(&t1);
  /*PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_mpoints_cell.(pvd,pvts,vts): CPU time %1.2e (sec) \n", PETSC_FUNCTION_NAME,prefix,t1-t0);*/

  PetscFunctionReturn(0);
}


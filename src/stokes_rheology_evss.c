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
 **    filename:   stokes_rheology_evss.c
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

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPEnergy_def.h"
#include "MPntPStokesPl_def.h"
#include "MPntPEVSS_def.h"

#include "QPntVolCoefStokes_def.h"
#include "QPntVolCoefEnergy_def.h"
#include "QPntVolCoefSymTens_def.h"

#include "dmda_duplicate.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1.h"
#include "data_bucket.h"
#include "output_paraview.h"
#include "quadrature.h"
#include "element_type_Q2.h"
#include "material_point_utils.h"
#include "element_utils_q2.h"
#include "element_utils_q1.h"


#undef __FUNCT__
#define __FUNCT__ "MPntPEVSSProjection_P0"
PetscErrorCode MPntPEVSSProjection_P0(const int npoints,MPntStd mp_std[],MPntPEVSS mp_symtens[],DM dmu,Quadrature Q)
{
	PetscErrorCode     ierr;
  DataField          PField;
  PetscInt           nqp,q,p,e,nel,k;
  QPntVolCoefSymTens *quadrature_T,*cell_quadrature_T;
	
	PetscFunctionBegin;
  
  /* Get quadrature point data */
  nqp = Q->npoints;
  
  if (nqp < 2)  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"P0 projection requires at least two quadrature points per cell <lame/lazy design issue>");
  
  DataBucketGetDataFieldByName(Q->properties_db,QPntVolCoefSymTens_classname,&PField);
  DataFieldGetEntries(PField,(void**)&quadrature_T);
  
  ierr = DMDAGetElements_pTatinQ2P1(dmu,&nel,NULL,NULL);CHKERRQ(ierr);
  
	/* traverse elements and initialize zeroth component of stress on each quadrature point */
	for (e=0; e<nel; e++) {
    cell_quadrature_T = &quadrature_T[e*nqp];
    
    cell_quadrature_T[0].T[voigt_xx] = 0.0;
    cell_quadrature_T[0].T[voigt_yy] = 0.0;
    cell_quadrature_T[0].T[voigt_zz] = 0.0;
    
    cell_quadrature_T[0].T[voigt_xy] = 0.0;
    cell_quadrature_T[0].T[voigt_xz] = 0.0;
    cell_quadrature_T[0].T[voigt_yz] = 0.0;
    
    cell_quadrature_T[1].T[voigt_xx] = 0.0;
	}
  
  /* traverse points and sum */
  for (p=0; p<npoints; p++) {
    e = mp_std[p].wil;
    cell_quadrature_T = &quadrature_T[e*nqp];
    
    for (k=0; k<6; k++) {
      cell_quadrature_T[0].T[k] += mp_symtens[p].tau[k];
    }
    cell_quadrature_T[1].T[voigt_xx] += 1.0;
  }
  
  /* assign constant values to each cell */
  for (e=0; e<nel; e++) {
    PetscReal T0[6],np_per_el;
    
    cell_quadrature_T = &quadrature_T[e*nqp];
    
    np_per_el = cell_quadrature_T[1].T[voigt_xx];
    for (k=0; k<6; k++) {
      T0[k] = cell_quadrature_T[0].T[k] / np_per_el;
    }
    
    for (q=0; q<nqp; q++) {
      cell_quadrature_T[q].T[voigt_xx] = T0[voigt_xx];
      cell_quadrature_T[q].T[voigt_yy] = T0[voigt_yy];
      cell_quadrature_T[q].T[voigt_zz] = T0[voigt_zz];
      
      cell_quadrature_T[q].T[voigt_xy] = T0[voigt_xy];
      cell_quadrature_T[q].T[voigt_xz] = T0[voigt_xz];
      cell_quadrature_T[q].T[voigt_yz] = T0[voigt_yz];
    }
	}
	PetscFunctionReturn(0);
}


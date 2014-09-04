/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      stokes_rheology_viscous.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

#include "petsc.h"

#include "ptatin3d.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "swarm_fields.h"

#undef __FUNCT__
#define __FUNCT__ "EvaluateRheologyNonlinearitiesMarkers_Viscous"
PetscErrorCode EvaluateRheologyNonlinearitiesMarkers_Viscous(pTatinCtx user,DM dau,PetscScalar u[],DM dap,PetscScalar p[])
{
	PetscErrorCode ierr;
	
	int            pidx,n_mp_points;
	DataBucket     db;
	DataField      PField_std, PField_stokes;
	PetscScalar    min_eta,max_eta,min_eta_g,max_eta_g;
	PetscLogDouble t0,t1;
	double         eta_mp;
	
	PetscFunctionBegin;
	
	PetscTime(&t0);

	ierr = pTatinGetMaterialPoints(user,&db,NULL);CHKERRQ(ierr);
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));
	
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	/* marker loop */
	min_eta = 1.0e100;
	max_eta = 1.0e-100;
	for (pidx=0; pidx<n_mp_points; pidx++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		
		DataFieldAccessPoint(PField_std,   pidx,(void**)&material_point);
		DataFieldAccessPoint(PField_stokes,pidx,(void**)&mpprop_stokes);

		/* get viscosity on marker */
		MPntPStokesGetField_eta_effective(mpprop_stokes,&eta_mp);
		
		/* monitor bounds */
		if (eta_mp > max_eta) { max_eta = eta_mp; }
		if (eta_mp < min_eta) { min_eta = eta_mp; }
		
	}
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);

  ierr = MPI_Allreduce(&min_eta,&min_eta_g,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&max_eta,&max_eta_g,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscTime(&t1);
	
	PetscPrintf(PETSC_COMM_WORLD,"Update rheology (viscous) [mpoint]: (min,max)_eta %1.2e,%1.2e; log10(max/min) %1.2e; cpu time %1.2e (sec)\n",
              min_eta_g, max_eta_g, log10(max_eta_g/min_eta_g), t1-t0 );
	
	PetscFunctionReturn(0);
}


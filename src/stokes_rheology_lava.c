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
 **    Filename:      stokes_rheology_lava.c
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
#include "petscdm.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "data_bucket.h"
#include "element_type_Q2.h"
#include "dmda_element_q2p1.h"
#include "element_utils_q2.h"
#include "dmdae.h"
#include "dmda_element_q1.h"
#include "element_utils_q1.h"
#include "phys_comp_energy.h"
#include "ptatin3d_energy.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"

#include "material_constants.h"

#define ETA_SCALE 1.0e2

static inline void ComputeStressIsotropic3d(PetscReal eta,double D[NSD][NSD],double T[NSD][NSD])
{
	const double two_eta = 2.0 * eta;
	
	T[0][0] = two_eta * D[0][0];	T[0][1] = two_eta * D[0][1];		T[0][2] = two_eta * D[0][2];
	T[1][0] =           T[0][1];	T[1][1] = two_eta * D[1][1];		T[1][2] = two_eta * D[1][2];
	T[2][0] =           T[0][2];	T[2][1] =           T[1][2];		T[2][2] = two_eta * D[2][2];	
}

static inline void ComputeStrainRate3d(double ux[],double uy[],double uz[],double dNudx[],double dNudy[],double dNudz[],double D[NSD][NSD])
{
	int    k;
	double exx,eyy,ezz,exy,exz,eyz;
	
	exx=0.0;  eyy=0.0;  ezz=0.0;
	exy=0.0;  exz=0.0;  eyz=0.0;
	
	for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
		exx += dNudx[k] * ux[k];
		eyy += dNudy[k] * uy[k];
		ezz += dNudz[k] * uz[k];
		
		exy += dNudy[k] * ux[k] + dNudx[k] * uy[k];
		exz += dNudz[k] * ux[k] + dNudx[k] * uz[k];
		eyz += dNudz[k] * uy[k] + dNudy[k] * uz[k];
	}
	exy = 0.5 * exy;
	exz = 0.5 * exz;
	eyz = 0.5 * eyz;
	
	D[0][0] = exx;		D[0][1] = exy;		D[0][2] = exz;
	D[1][0] = exy;		D[1][1] = eyy;		D[1][2] = eyz;
	D[2][0] = exz;		D[2][1] = eyz;		D[2][2] = ezz;
}

static inline void ComputeDeformationGradient3d(double ux[],double uy[],double uz[],double dNudx[],double dNudy[],double dNudz[],double L[NSD][NSD])
{
	int i,j,k;
	
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			L[i][j] = 0.0;
		}
	}
	for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
		// du/dx_i
		L[0][0] += dNudx[k] * ux[k];
		L[0][1] += dNudy[k] * ux[k];
		L[0][2] += dNudz[k] * ux[k];
		// dv/dx_i
		L[1][0] += dNudx[k] * uy[k];
		L[1][1] += dNudy[k] * uy[k];
		L[1][2] += dNudz[k] * uy[k];
		// dw/dx_i
		L[2][0] += dNudx[k] * uz[k];
		L[2][1] += dNudy[k] * uz[k];
		L[2][2] += dNudz[k] * uz[k];
	}
}

static inline void ComputeSecondInvariant3d(double A[NSD][NSD],double *A2)
{
	int i,j;
	double sum = 0.0;
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			sum = sum + A[i][j]*A[i][j];
		}
	}
	*A2 = sqrt( 0.5 * sum );	
}

static inline void ComputeAverageTrace3d(double A[NSD][NSD],double *A2)
{
	const double one_third = 0.333333333333333;
	
	*A2 = one_third * ( A[0][0] + A[1][1] + A[2][2] );
	
}

#undef __FUNCT__
#define __FUNCT__ "private_EvaluateRheologyNonlinearitiesMarkers_LAVA"
PetscErrorCode private_EvaluateRheologyNonlinearitiesMarkers_LAVA(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],DM daT,PetscScalar Tfield[])
{
	PetscErrorCode ierr;
	
	int            pidx,n_mp_points;
	DataBucket     db,material_constants;
	DataField      PField_std,PField_stokes,PField_pls;
	PetscScalar    min_eta,max_eta,min_eta_g,max_eta_g;
	PetscLogDouble t0,t1;
	
	DM             cda;
	Vec            gcoords,gcoords_T;
	PetscReal      *LA_gcoords,*LA_gcoords_T;
	PetscInt       nel,nen_u,nen_p,eidx,k;
	PetscInt       nel_T,nen_T;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	const PetscInt *elnidx_T;
	PetscInt       vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscReal      elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal      elu[3*Q2_NODES_PER_EL_3D];
	PetscReal      elT[Q1_NODES_PER_EL_3D];
	PetscReal      ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	
	PetscReal NI_T[Q1_NODES_PER_EL_3D];
	PetscReal NI[Q2_NODES_PER_EL_3D],GNI[3][Q2_NODES_PER_EL_3D];
	PetscReal dNudx[Q2_NODES_PER_EL_3D],dNudy[Q2_NODES_PER_EL_3D],dNudz[Q2_NODES_PER_EL_3D];
	
	double         eta_mp;
	double         D_mp[NSD][NSD],Tpred_mp[NSD][NSD];
	double         inv2_D_mp,inv2_Tpred_mp;
	
	//DataField      PField_MatTypes;
	DataField      PField_ViscConst;
    /* structs or material constants */
	//MaterialConst_MaterialType      *MatType_data;
	MaterialConst_ViscosityConst    *ViscConst_data;
	long int       npoints_yielded,npoints_yielded_g;
	
	PetscFunctionBegin;
	
	PetscTime(&t0);
	
	/* access material point information */
	ierr = pTatinGetMaterialPoints(user,&db,NULL);CHKERRQ(ierr);
	/* PField_std global index marker, phase marker, ...*/
	/* PField_stokes contains: etaf, rhof */
	/* PField_pls contains: accumulated plastic strain, yield type */
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	//
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	//
	DataBucketGetDataFieldByName(db,MPntPStokesPl_classname,&PField_pls);
	DataFieldGetAccess(PField_pls);
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	/* get u,p element information */
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	if (daT) {
		/* access the coordinates for the temperature mesh */ /* THIS IS NOT ACTUALLY NEEDED */
		ierr = DMGetCoordinatesLocal(daT,&gcoords_T);CHKERRQ(ierr);
		ierr = VecGetArray(gcoords_T,&LA_gcoords_T);CHKERRQ(ierr);
		
		ierr = DMDAGetElementsQ1(daT,&nel_T,&nen_T,&elnidx_T);CHKERRQ(ierr);
		
		if (nel_T != nel) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Require code update to utilize nested Q1 mesh for temperature");
		}
	}
	
	/* access material constants */
	ierr = pTatinGetMaterialConstants(user,&material_constants);CHKERRQ(ierr);
	
	//DataBucketGetDataFieldByName(material_constants,MaterialConst_MaterialType_classname,  &PField_MatTypes);
	//MatType_data           = (MaterialConst_MaterialType*)PField_MatTypes->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_ViscosityConst_classname,&PField_ViscConst);
	ViscConst_data         = (MaterialConst_ViscosityConst*)PField_ViscConst->data;
	
	/* marker loop */
	min_eta = 1.0e100;
	max_eta = 1.0e-100;
	npoints_yielded = 0;
	
	for (pidx=0; pidx<n_mp_points; pidx++) {
		MPntStd       *mpprop_std;
		MPntPStokes   *mpprop_stokes;
		MPntPStokesPl *mpprop_pls;
		double        *xi_p;
		double        T_mp;
		int           region_idx;
		
		DataFieldAccessPoint(PField_std,   pidx,(void**)&mpprop_std);
		DataFieldAccessPoint(PField_stokes,pidx,(void**)&mpprop_stokes);
		DataFieldAccessPoint(PField_pls,   pidx,(void**)&mpprop_pls);
		
		/* Get marker types */
		region_idx = mpprop_std->phase;

		/* Get index of element containing this marker */
		eidx = mpprop_std->wil;
		
		/* Get element indices */
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*eidx]);CHKERRQ(ierr);
		
		/* Get element coordinates */
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*eidx],LA_gcoords);CHKERRQ(ierr);
		
		/* Get element velocity */
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*eidx],ufield);CHKERRQ(ierr);
		
		if (daT) {
			/* Get element temperature */
			ierr = DMDAEQ1_GetScalarElementField_3D(elT,(PetscInt*)&elnidx_T[nen_T*eidx],Tfield);CHKERRQ(ierr);
		}
		
		/* Get local coordinate of marker */
		xi_p = mpprop_std->xi;
		
		/* Prepare basis functions */
		/* grad.Ni */
		P3D_ConstructGNi_Q2_3D(xi_p,GNI);
		
		/* Get shape function derivatives */
		P3D_evaluate_global_derivatives_Q2(elcoords,GNI,dNudx,dNudy,dNudz);
		
		/* get velocity components */
		for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
		
		pTatin_ConstructNi_Q2_3D( xi_p, NI );
		
		
		T_mp = 0.0;
		if (daT) {
			/* Interpolate the temperature */
			/* NOTE: scaling is requred of xi_p if nested mesh is used */
			P3D_ConstructNi_Q1_3D(xi_p,NI_T);
            
			for (k=0; k<Q1_NODES_PER_EL_3D; k++) {
				T_mp += NI_T[k] * elT[k];
			}
		}
		
		/* get viscosity on marker */
		/*
		 Reference:
		 
		 Ross W. Griffiths. 
		 The dynamics of lava flows. 
		 Annual review of fluid mechanics, 
		 32(1), pp. 477-518, (2000)
		*/
		if (region_idx != 0) {
			const PetscReal phi_max = 0.68;
			const PetscReal gamma   = -0.04;
			const PetscReal T_e     = 1100.0; /* temperature in degrees C */
			const PetscReal T_s     = 600.0; /* temperature in degrees C */
			const PetscReal phi_0   = 0.0;
			PetscReal eta0,ratio,t_dep,phi_f,phi;
			
			eta0 = ViscConst_data[ region_idx ].eta0;
			phi_f = phi_max - phi_0;
			phi = phi_0 + phi_f * (T_e - T_mp) / (T_e - T_s);
			/*
			if (phi > phi_max) {
				phi = phi_max;
				ratio = 1.0e32;
			} else {
				ratio = pow( 1.0 - phi / phi_max, -2.5 );
			}
			*/
			if (phi >= phi_max) {
				phi = phi_max - 1.0e-10;
			}
			
			ratio = pow( 1.0 - phi / phi_max, -2.5 );
			t_dep = exp( -gamma * (T_e - T_mp) );
			
			eta_mp = eta0 * ratio * t_dep;

			//{
			//	double rad = sqrt(mpprop_std->coor[0]*mpprop_std->coor[0] + mpprop_std->coor[1]*mpprop_std->coor[1]);
			//	printf("rad %1.5e : T %1.5e : phi %1.5e : eta %1.5e \n",rad,T_mp,phi,eta_mp);
			//}
			
			if (eta_mp > 1.0e5/ETA_SCALE) {
				eta_mp = 1.0e5/ETA_SCALE;
			}
			
			
		} else {
			eta_mp = ViscConst_data[ region_idx ].eta0;
		}
		
		/* apply stress limiters to all markers not considered "air" */
		if (region_idx != 0) {
			double tau_yield_mp;
			
			/* Compute yield surface */
			/* 
			 Reference:
			 
			 Hideaki Miyamoto and Sho Sasaki. 
			 Numerical simulations of flood basalt lava flows: Roles of parameters on lava flow morphologies. 
			 Journal of Geophysical Research, 
			 103(B11), pp. 27489-27502, (1998).
			*/
			tau_yield_mp = pow( 10.0, 11.59 - 0.0089 * T_mp );
			tau_yield_mp = tau_yield_mp / ETA_SCALE;
			
			/* strain rate */
			ComputeStrainRate3d(ux,uy,uz,dNudx,dNudy,dNudz,D_mp);
			/* stress */
			ComputeStressIsotropic3d(eta_mp,D_mp,Tpred_mp);
			/* second inv stress */
			ComputeSecondInvariant3d(Tpred_mp,&inv2_Tpred_mp);
			
			MPntPStokesPlSetField_yield_indicator(mpprop_pls,0);
			
			if (inv2_Tpred_mp > tau_yield_mp) {
				ComputeSecondInvariant3d(D_mp,&inv2_D_mp);
				
				eta_mp = 0.5 * tau_yield_mp / inv2_D_mp;
				MPntPStokesPlSetField_yield_indicator(mpprop_pls,1);
				npoints_yielded++;
			}
			
		}		
		
		/* update viscosity on marker */
		MPntPStokesSetField_eta_effective(mpprop_stokes,eta_mp);

		/* monitor bounds */
		if (eta_mp > max_eta) { max_eta = eta_mp; }
		if (eta_mp < min_eta) { min_eta = eta_mp; }
	}  
	
	DataFieldRestoreAccess(PField_pls);
	DataFieldRestoreAccess(PField_stokes);
	DataFieldRestoreAccess(PField_std);
    
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	if (daT) {
		ierr = VecRestoreArray(gcoords_T,&LA_gcoords_T);CHKERRQ(ierr);
	}
	
	
	
	ierr = MPI_Allreduce(&min_eta,&min_eta_g,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&max_eta,&max_eta_g,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&npoints_yielded,&npoints_yielded_g,1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	PetscTime(&t1);
	
	PetscPrintf(PETSC_COMM_WORLD,"Update non-linearities (LAVA) [mpoint]: (min,max)_eta %1.2e,%1.2e; log10(max/min) %1.2e; npoints_yielded %ld; cpu time %1.2e (sec)\n",
                min_eta_g, max_eta_g, log10(max_eta_g/min_eta_g), npoints_yielded_g, t1-t0 );
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateRheologyNonlinearitiesMarkers_LAVA"
PetscErrorCode EvaluateRheologyNonlinearitiesMarkers_LAVA(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[])
{
	PhysCompEnergy energy;
	PetscBool found;
	DM daT;
	Vec temperature,temperature_l;
	PetscScalar *LA_temperature_l;
	PetscErrorCode ierr;
    
	PetscFunctionBegin;
	
	ierr = pTatinContextValid_Energy(user,&found);CHKERRQ(ierr);
	if (found) {
		ierr = pTatinGetContext_Energy(user,&energy);CHKERRQ(ierr);
		ierr = pTatinPhysCompGetData_Energy(user,&temperature,NULL);CHKERRQ(ierr);
		daT  = energy->daT;
        
		ierr = DMGetLocalVector(daT,&temperature_l);CHKERRQ(ierr);
		ierr = VecZeroEntries(temperature_l);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(daT,temperature,INSERT_VALUES,temperature_l);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  (daT,temperature,INSERT_VALUES,temperature_l);CHKERRQ(ierr);
		ierr = VecGetArray(temperature_l,&LA_temperature_l);CHKERRQ(ierr);
        
		ierr = private_EvaluateRheologyNonlinearitiesMarkers_LAVA(user,dau,ufield,dap,pfield,daT,LA_temperature_l);CHKERRQ(ierr);
		
		ierr = VecRestoreArray(temperature_l,&LA_temperature_l);CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(daT,&temperature_l);CHKERRQ(ierr);
		
	} else {
		ierr = private_EvaluateRheologyNonlinearitiesMarkers_LAVA(user,dau,ufield,dap,pfield,NULL,NULL);CHKERRQ(ierr);
	}
	
	
	PetscFunctionReturn(0);
}

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
 **    filename:   stokes_rheology_vp_std.c
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
#include "MPntPStokesPl_def.h"

#include "material_constants.h"

typedef enum { YTYPE_NONE=0, YTYPE_MISES=1, YTYPE_DP=2, YTYPE_TENSILE_FAILURE=3 } YieldTypeDefinition;

static inline void ComputeLinearSoft(float eplast,PetscReal emin,PetscReal emax, PetscReal X0, PetscReal Xinf, PetscReal *Xeff)
{
	*Xeff = X0;
	if (eplast > emin) {
		if (eplast > emax) {
			*Xeff = Xinf;
		} else {
			*Xeff  = X0 - (eplast-emin)/(emax-emin)*(X0-Xinf);
		}
	}
}

static inline void ComputeExponentialSoft(float eplast,PetscReal emin,PetscReal efold, PetscReal X0, PetscReal Xinf, PetscReal *Xeff)
{
	*Xeff = X0;
	if (eplast > emin) {
		*Xeff  = Xinf + (X0-Xinf) * exp(-(eplast-emin)/efold);
	}
}



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

#if 0
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
#endif

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

#if 0
static inline void ComputeAverageTrace3d(double A[NSD][NSD],double *A2)
{
	const double one_third = 0.333333333333333;
	
	*A2 = one_third * ( A[0][0] + A[1][1] + A[2][2] );
	
}
#endif

#undef __FUNCT__
#define __FUNCT__ "private_EvaluateRheologyNonlinearitiesMarkers_VPSTD"
PetscErrorCode private_EvaluateRheologyNonlinearitiesMarkers_VPSTD(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[],DM daT,PetscScalar Tfield[])
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
	PetscReal      elu[3*Q2_NODES_PER_EL_3D],elp[P_BASIS_FUNCTIONS];
	PetscReal      elT[Q1_NODES_PER_EL_3D];
	PetscReal      ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	
	PetscReal NI_T[Q1_NODES_PER_EL_3D];
	PetscReal NI[Q2_NODES_PER_EL_3D],GNI[3][Q2_NODES_PER_EL_3D],NIp[P_BASIS_FUNCTIONS];
	PetscReal dNudx[Q2_NODES_PER_EL_3D],dNudy[Q2_NODES_PER_EL_3D],dNudz[Q2_NODES_PER_EL_3D];
	
	double         eta_mp;
	double         D_mp[NSD][NSD],Tpred_mp[NSD][NSD];
	double         inv2_D_mp,inv2_Tpred_mp;
	
	DataField      PField_MatTypes;
	DataField      PField_DensityConst,PField_DensityBoussinesq;
	DataField      PField_ViscConst,PField_ViscZ,PField_ViscFK,PField_ViscArrh;
	DataField      PField_PlasticMises,PField_PlasticDP;
	DataField      PField_SoftLin,PField_SoftExpo;
    /* structs or material constants */
	MaterialConst_MaterialType      *MatType_data;
	MaterialConst_DensityConst      *DensityConst_data;
	MaterialConst_DensityBoussinesq *DensityBoussinesq_data;
	MaterialConst_ViscosityConst    *ViscConst_data;
	MaterialConst_ViscosityZ        *ViscZ_data;
	MaterialConst_ViscosityFK       *ViscFK_data;
	MaterialConst_ViscosityArrh     *ViscArrh_data;
	MaterialConst_PlasticMises      *PlasticMises_data;
	MaterialConst_PlasticDP         *PlasticDP_data;
	MaterialConst_SoftLin           *SoftLin_data;
	MaterialConst_SoftExpo          *SoftExpo_data;
	
	int            viscous_type,plastic_type,softening_type,density_type;
	long int       npoints_yielded,npoints_yielded_g;
	
	
	PetscFunctionBegin;
	
	PetscTime(&t0);
	
	/* access material point information */
	ierr = pTatinGetMaterialPoints(user,&db,NULL);CHKERRQ(ierr);
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	
    /* PField_std global index marker, phase marker, ...*/
    /* PField_stokes contains: etaf, rhof */
    /* PField_pls contains: accumulated plastic strain, yield type */
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	
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
		/* access the energy variables stored on the markers */
		//DataBucketGetDataFieldByName(db,MPntPEnergy_classname,&PField_energy);
		//DataFieldGetAccess(PField_energy);
		
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
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_MaterialType_classname,  &PField_MatTypes);
	MatType_data           = (MaterialConst_MaterialType*)PField_MatTypes->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_DensityConst_classname,  &PField_DensityConst);
	DensityConst_data      = (MaterialConst_DensityConst*)PField_DensityConst->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_DensityBoussinesq_classname,  &PField_DensityBoussinesq);
	DensityBoussinesq_data = (MaterialConst_DensityBoussinesq*)PField_DensityBoussinesq->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_ViscosityConst_classname,&PField_ViscConst);
	ViscConst_data         = (MaterialConst_ViscosityConst*)PField_ViscConst->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_ViscosityArrh_classname,&PField_ViscArrh);
	ViscArrh_data          = (MaterialConst_ViscosityArrh*)PField_ViscArrh->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_ViscosityFK_classname,&PField_ViscFK);
	ViscFK_data            = (MaterialConst_ViscosityFK*)PField_ViscFK->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_ViscosityZ_classname,&PField_ViscZ);
	ViscZ_data             = (MaterialConst_ViscosityZ*)PField_ViscZ->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_PlasticMises_classname,  &PField_PlasticMises);
	PlasticMises_data      = (MaterialConst_PlasticMises*)  PField_PlasticMises->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_PlasticDP_classname,  &PField_PlasticDP);
	PlasticDP_data         = (MaterialConst_PlasticDP*)  PField_PlasticDP->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_SoftLin_classname,  &PField_SoftLin);
	SoftLin_data           = (MaterialConst_SoftLin*)  PField_SoftLin->data;
	
	DataBucketGetDataFieldByName(material_constants,MaterialConst_SoftExpo_classname,  &PField_SoftExpo);
	SoftExpo_data          = (MaterialConst_SoftExpo*)  PField_SoftExpo->data;
	
	
	/* marker loop */
	min_eta = 1.0e100;
	max_eta = 1.0e-100;
	npoints_yielded = 0;
	for (pidx=0; pidx<n_mp_points; pidx++) {
		MPntStd       *mpprop_std;
		MPntPStokes   *mpprop_stokes;
		MPntPStokesPl *mpprop_pls;
		double        *xi_p;
		double        pressure_mp,y_mp,T_mp;
		int           region_idx;
		
		DataFieldAccessPoint(PField_std,   pidx,(void**)&mpprop_std);
		DataFieldAccessPoint(PField_stokes,pidx,(void**)&mpprop_stokes);
		DataFieldAccessPoint(PField_pls,   pidx,(void**)&mpprop_pls);
		
		/* Get marker types */
		region_idx = mpprop_std->phase;
		
		viscous_type   = MatType_data[ region_idx ].visc_type;
		plastic_type   = MatType_data[ region_idx ].plastic_type;
		density_type   = MatType_data[ region_idx ].density_type;
		softening_type = MatType_data[ region_idx ].softening_type;
		
		/* Get index of element containing this marker */
		eidx = mpprop_std->wil;
		
		/* Get element indices */
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*eidx]);CHKERRQ(ierr);
		
		/* Get element coordinates */
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*eidx],LA_gcoords);CHKERRQ(ierr);
		
		/* Get element velocity */
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*eidx],ufield);CHKERRQ(ierr);
		
		/* Get element pressure */
		ierr = DMDAGetScalarElementField(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*eidx],pfield);CHKERRQ(ierr);
		
		if (daT) {
			/* Get element temperature */
			ierr = DMDAEQ1_GetScalarElementField_3D(elT,(PetscInt*)&elnidx_T[nen_T*eidx],Tfield);CHKERRQ(ierr);
		}
		
		/* Get local coordinate of marker */
		xi_p = mpprop_std->xi;
		
		/* Prepare basis functions */
		/* grad.Ni */
		P3D_ConstructGNi_Q2_3D(xi_p,GNI);
		/* Mi */
		ConstructNi_pressure(xi_p,elcoords,NIp);
		
		/* Get shape function derivatives */
		P3D_evaluate_global_derivatives_Q2(elcoords,GNI,dNudx,dNudy,dNudz);
		
		/* Compute pressure at material point */
		pressure_mp = 0.0;
		for (k=0; k<P_BASIS_FUNCTIONS; k++) {
			pressure_mp += NIp[k] * elp[k];
		}
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
		//MPntPStokesGetField_eta_effective(mpprop_stokes,&eta_mp);		
		switch (viscous_type) {
				
			case VISCOUS_CONSTANT: {
				eta_mp = ViscConst_data[ region_idx ].eta0;
			}
				break;
				
			case VISCOUS_Z: {
				y_mp = mpprop_std->coor[1];
				eta_mp = ViscZ_data[ region_idx ].eta0*exp(-(ViscZ_data[ region_idx ].zref-y_mp)/ViscZ_data[ region_idx ].zeta);
				
			}
				break;
				
			case VISCOUS_FRANKK: {
				eta_mp  = ViscFK_data[ region_idx ].eta0*exp(-ViscFK_data[ region_idx ].theta*T_mp);
			}
				break;
				
			case VISCOUS_ARRHENIUS: {
				PetscScalar R       = 8.31440;
				PetscReal nexp      = ViscArrh_data[ region_idx ].nexp;
				PetscReal entalpy   = ViscArrh_data[ region_idx ].entalpy;
				PetscReal preexpA   = ViscArrh_data[ region_idx ].preexpA;
				PetscReal Vmol      = ViscArrh_data[ region_idx ].Vmol;
				PetscReal Tref      = ViscArrh_data[ region_idx ].Tref;
				PetscReal Ascale    = ViscArrh_data[ region_idx ].Ascale;
				PetscReal T_arrh    = T_mp + Tref ;
				PetscReal sr, eta, pressure; 
				
				ComputeStrainRate3d(ux,uy,uz,dNudx,dNudy,dNudz,D_mp);
                ComputeSecondInvariant3d(D_mp,&inv2_D_mp);
                
                sr = inv2_D_mp/ViscArrh_data[ region_idx ].Eta_scale*ViscArrh_data[ region_idx ].P_scale;
				if (sr < 1.0e-17) {
					sr = 1.0e-17;
				}
                
				pressure = ViscArrh_data[ region_idx ].P_scale*pressure_mp;
                
                if (pressure >= 0.0) {
                    entalpy = entalpy + pressure*Vmol;
                }

				eta  = Ascale*0.25*pow(sr,1.0/nexp - 1.0)*pow(0.75*preexpA,-1.0/nexp)*exp(entalpy/(nexp*R*T_arrh));
				eta_mp = eta/ViscArrh_data[ region_idx ].Eta_scale;
			}
				break;

			case VISCOUS_ARRHENIUS_2: {
				PetscScalar R       = 8.31440;
				PetscReal nexp      = ViscArrh_data[ region_idx ].nexp;
				PetscReal entalpy   = ViscArrh_data[ region_idx ].entalpy;
				PetscReal preexpA   = ViscArrh_data[ region_idx ].preexpA;
				PetscReal Vmol      = ViscArrh_data[ region_idx ].Vmol;
				PetscReal Tref      = ViscArrh_data[ region_idx ].Tref;
				PetscReal Ascale    = ViscArrh_data[ region_idx ].Ascale;
				PetscReal T_arrh    = T_mp + Tref ;
				PetscReal sr, eta, pressure;
				
				ComputeStrainRate3d(ux,uy,uz,dNudx,dNudy,dNudz,D_mp);
                ComputeSecondInvariant3d(D_mp,&inv2_D_mp);
                
                sr = inv2_D_mp/ViscArrh_data[ region_idx ].Eta_scale*ViscArrh_data[ region_idx ].P_scale;
				if (sr < 1.0e-17) {
					sr = 1.0e-17;
				}
                
				pressure = ViscArrh_data[ region_idx ].P_scale*pressure_mp;
                
                if (pressure >= 0.0) {
                    entalpy = entalpy + pressure*Vmol;
                }

				eta  = Ascale*pow(sr,1.0/nexp - 1.0)*pow(preexpA,-1.0/nexp)*exp(entalpy/(nexp*R*T_arrh));
				eta_mp = eta/ViscArrh_data[ region_idx ].Eta_scale;
			}
				break;
                
            default:
                SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"No default ViscousType specified. Valid choices are VISCOUS_CONSTANT, VISCOUS_Z, VISCOUS_FRANKK, VISCOUS_ARRHENIUS, VISCOUS_ARRHENIUS_2");
                
		}
		
		switch (plastic_type) {
				
			case PLASTIC_NONE: {
			}
				break;
				
			case PLASTIC_MISES: {
				double tau_yield_mp  = PlasticMises_data[ region_idx ].tau_yield;
				double tau_yield_inf = PlasticMises_data[ region_idx ].tau_yield_inf;
				
				switch (MatType_data[ region_idx ].softening_type) {
					case SOFTENING_NONE: {
						
					}
						break;
                        
					case SOFTENING_LINEAR: {
						float     eplastic;
						PetscReal emin = SoftLin_data[ region_idx ].eps_min;
						PetscReal emax = SoftLin_data[ region_idx ].eps_max;
						
						MPntPStokesPlGetField_plastic_strain(mpprop_pls,&eplastic);
						ComputeLinearSoft(eplastic,emin,emax,tau_yield_mp, tau_yield_inf, &tau_yield_mp);
					}
						break;    
						
					case SOFTENING_EXPONENTIAL: {
						float eplastic;
						PetscReal emin     = SoftExpo_data[ region_idx ].eps_min;
						PetscReal efold    = SoftExpo_data[ region_idx ].eps_fold; 
						MPntPStokesPlGetField_plastic_strain(mpprop_pls,&eplastic);
						ComputeExponentialSoft(eplastic,emin,efold,tau_yield_mp, tau_yield_inf, &tau_yield_mp);
					}
						break;
                    default:
                        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"No default SofteningType set. Valid choices are SOFTENING_NONE, SOFTENING_LINEAR, SOFTENING_EXPONENTIAL");
				}
				
				/* mark all markers as not yielding */
                MPntPStokesPlSetField_yield_indicator(mpprop_pls,YTYPE_NONE);
				
				/* strain rate */
				ComputeStrainRate3d(ux,uy,uz,dNudx,dNudy,dNudz,D_mp);
				/* stress */
				ComputeStressIsotropic3d(eta_mp,D_mp,Tpred_mp);
				/* second inv stress */
				ComputeSecondInvariant3d(Tpred_mp,&inv2_Tpred_mp);
				
				if (inv2_Tpred_mp > tau_yield_mp) {
					ComputeSecondInvariant3d(D_mp,&inv2_D_mp);
					
					eta_mp = 0.5 * tau_yield_mp / inv2_D_mp;
					npoints_yielded++;
					MPntPStokesPlSetField_yield_indicator(mpprop_pls,YTYPE_MISES);
				}
			}
				break;
				
			case PLASTIC_MISES_H: {
				double tau_yield_mp  = PlasticMises_data[ region_idx ].tau_yield;
				double tau_yield_inf = PlasticMises_data[ region_idx ].tau_yield;
				double eta_flow_mp,eta_yield_mp;
                
				switch (MatType_data[ region_idx ].softening_type) {
					case SOFTENING_NONE: {
						
					}
						break;
                        
					case SOFTENING_LINEAR: {
						float     eplastic;
						PetscReal emin = SoftLin_data[ region_idx ].eps_min;
						PetscReal emax = SoftLin_data[ region_idx ].eps_max;
						
						MPntPStokesPlGetField_plastic_strain(mpprop_pls,&eplastic);
						ComputeLinearSoft(eplastic,emin,emax,tau_yield_mp,tau_yield_inf,&tau_yield_mp);
					}
						break;
						
					case SOFTENING_EXPONENTIAL: {
						float eplastic;
						PetscReal emin     = SoftExpo_data[ region_idx ].eps_min;
						PetscReal efold    = SoftExpo_data[ region_idx ].eps_fold;
						MPntPStokesPlGetField_plastic_strain(mpprop_pls,&eplastic);
						ComputeExponentialSoft(eplastic,emin,efold,tau_yield_mp,tau_yield_inf,&tau_yield_mp);
					}
						break;
                    default:
                        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"No default SofteningType set. Valid choices are SOFTENING_NONE, SOFTENING_LINEAR, SOFTENING_EXPONENTIAL");
				}
				
                eta_flow_mp = eta_mp;
                
				/* strain rate */
				ComputeStrainRate3d(ux,uy,uz,dNudx,dNudy,dNudz,D_mp);
                ComputeSecondInvariant3d(D_mp,&inv2_D_mp);
                
                eta_yield_mp = 0.5 * tau_yield_mp / inv2_D_mp;

                eta_mp = 1.0/ (  1.0/eta_flow_mp + 1.0/eta_yield_mp );
                
                npoints_yielded++;
			}
				break;

			case PLASTIC_DP: {
				double    tau_yield_mp;
				short     yield_type;
				PetscReal phi     = PlasticDP_data[ region_idx ].phi;
				PetscReal Co      = PlasticDP_data[ region_idx ].Co;
				PetscReal phi_inf = PlasticDP_data[ region_idx ].phi_inf;
				PetscReal Co_inf  = PlasticDP_data[ region_idx ].Co_inf;         
				
				
				switch (MatType_data[ region_idx ].softening_type) {
                        
					case SOFTENING_NONE: {
						
					}
						break;
						
					case SOFTENING_LINEAR: {
						float     eplastic;
						PetscReal emin = SoftLin_data[ region_idx ].eps_min;
						PetscReal emax = SoftLin_data[ region_idx ].eps_max; 
						
						MPntPStokesPlGetField_plastic_strain(mpprop_pls,&eplastic);
						ComputeLinearSoft(eplastic,emin,emax,Co , Co_inf, &Co);
						ComputeLinearSoft(eplastic,emin,emax,phi, phi_inf, &phi);
					}
						break;    
						
					case SOFTENING_EXPONENTIAL: {
						float     eplastic;
						PetscReal emin  = SoftExpo_data[ region_idx ].eps_min;
						PetscReal efold = SoftExpo_data[ region_idx ].eps_fold; 
						
						MPntPStokesPlGetField_plastic_strain(mpprop_pls,&eplastic);
						ComputeExponentialSoft(eplastic,emin,efold,Co , Co_inf, &Co);
						ComputeExponentialSoft(eplastic,emin,efold,phi, phi_inf, &phi);
					}
						break;
                    default:
                        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"No default SofteningType set. Valid choices are SOFTENING_NONE, SOFTENING_LINEAR, SOFTENING_EXPONENTIAL");
				}
				
				/* mark all markers as not yielding */
				MPntPStokesPlSetField_yield_indicator(mpprop_pls,YTYPE_NONE);
				
				/* compute yield surface */
				tau_yield_mp = sin(phi) * pressure_mp + Co * cos(phi);
				
				/* identify yield type */
				yield_type = YTYPE_MISES;
				
				if ( tau_yield_mp < PlasticDP_data[region_idx].tens_cutoff) {
					/* failure in tension cutoff */
					tau_yield_mp = PlasticDP_data[region_idx].tens_cutoff;
					yield_type = 2;
				} else if (tau_yield_mp > PlasticDP_data[region_idx].hst_cutoff) {   
					/* failure at High stress cut off à la boris */
					tau_yield_mp = PlasticDP_data[region_idx].hst_cutoff;
					yield_type = 3;
				} 
				
				/* strain rate */
				ComputeStrainRate3d(ux,uy,uz,dNudx,dNudy,dNudz,D_mp);
				/* stress */
				ComputeStressIsotropic3d(eta_mp,D_mp,Tpred_mp);
				/* second inv stress */
				ComputeSecondInvariant3d(Tpred_mp,&inv2_Tpred_mp);
				
				if (inv2_Tpred_mp > tau_yield_mp) {
					ComputeSecondInvariant3d(D_mp,&inv2_D_mp);
					
					eta_mp = 0.5 * tau_yield_mp / inv2_D_mp;
					npoints_yielded++;
					MPntPStokesPlSetField_yield_indicator(mpprop_pls,yield_type);
				}
			}
				break;

            case PLASTIC_DP_H: {

            }
                
            default:
                SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"No default PlasticType set. Valid choices are PLASTIC_NONE, PLASTIC_MISES, PLASTIC_DP");
		}
		
		
		/* update viscosity on marker */
		MPntPStokesSetField_eta_effective(mpprop_stokes,eta_mp);
		/* monitor bounds */
		if (eta_mp > max_eta) { max_eta = eta_mp; }
		if (eta_mp < min_eta) { min_eta = eta_mp; }
		
		
		
		switch (density_type) {
                
			case DENSITY_CONSTANT: {
				PetscReal rho_mp;
				
				rho_mp = DensityConst_data[region_idx].density;
				MPntPStokesSetField_density(mpprop_stokes,rho_mp);
			}
				break;
				
			case DENSITY_BOUSSINESQ: {	
				PetscReal rho_mp;
				PetscReal rho0  = DensityBoussinesq_data[region_idx].density;
				PetscReal alpha = DensityBoussinesq_data[region_idx].alpha;
				PetscReal beta  = DensityBoussinesq_data[region_idx].beta;
				
				rho_mp = rho0*(1-alpha*T_mp+beta*pressure_mp);
				MPntPStokesSetField_density(mpprop_stokes,rho_mp);         
			}
				break;
            default:
                SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"No default DensityType set. Valid choices are DENSITY_CONSTANT, DENSITY_BOUSSINESQ");
				
		}
		
		
		
		/* Global cutoffs for viscosity */
		/* Here should I store these? */
		
		
	}  
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	DataFieldRestoreAccess(PField_pls);
    
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	if (daT) {
		//DataFieldRestoreAccess(PField_energy);
		
		ierr = VecRestoreArray(gcoords_T,&LA_gcoords_T);CHKERRQ(ierr);
	}
	
	
	
	ierr = MPI_Allreduce(&min_eta,&min_eta_g,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&max_eta,&max_eta_g,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&npoints_yielded,&npoints_yielded_g,1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	PetscTime(&t1);
	
	PetscPrintf(PETSC_COMM_WORLD,"Update non-linearities (VPSTD) [mpoint]: (min,max)_eta %1.2e,%1.2e; log10(max/min) %1.2e; npoints_yielded %ld; cpu time %1.2e (sec)\n",
                min_eta_g, max_eta_g, log10(max_eta_g/min_eta_g), npoints_yielded_g, t1-t0 );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateRheologyNonlinearitiesMarkers_VPSTD"
PetscErrorCode EvaluateRheologyNonlinearitiesMarkers_VPSTD(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[])
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
        
		ierr = private_EvaluateRheologyNonlinearitiesMarkers_VPSTD(user,dau,ufield,dap,pfield,daT,LA_temperature_l);CHKERRQ(ierr);
		
		ierr = VecRestoreArray(temperature_l,&LA_temperature_l);CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(daT,&temperature_l);CHKERRQ(ierr);
		
	} else {
		ierr = private_EvaluateRheologyNonlinearitiesMarkers_VPSTD(user,dau,ufield,dap,pfield,NULL,NULL);CHKERRQ(ierr);
	}
	
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ApplyViscosityCutOffMarkers_VPSTD"
PetscErrorCode ApplyViscosityCutOffMarkers_VPSTD(pTatinCtx user)
{
	PetscErrorCode    ierr;
	int               pidx,n_mp_points;
	DataBucket        db;
	DataField         PField_std,PField_stokes;
	double            min_eta,max_eta,min_eta_g,max_eta_g,min_cutoff,max_cutoff;
	PetscLogDouble    t0,t1;
	RheologyConstants *rheology;
	double            eta_mp;
	long int          npoints_cutoff[2],npoints_cutoff_g[2];
	
	PetscFunctionBegin;
	
	PetscTime(&t0);
	
	/* access material point information */
	ierr = pTatinGetMaterialPoints(user,&db,NULL);CHKERRQ(ierr);
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	rheology = &user->rheology_constants;
	/* marker loop */
	min_eta = 1.0e100;
	max_eta = 1.0e-100;
	npoints_cutoff[0] = 0;
	npoints_cutoff[1] = 0;
	
	min_cutoff = rheology->eta_lower_cutoff_global;
	max_cutoff = rheology->eta_upper_cutoff_global;
	
	for (pidx=0; pidx<n_mp_points; pidx++) {
		MPntStd     *mpprop_std;
		MPntPStokes *mpprop_stokes;
		PetscScalar min_cutoff_l,max_cutoff_l;
		int         region_idx;
		
		DataFieldAccessPoint(PField_std,   pidx,(void**)&mpprop_std);
		DataFieldAccessPoint(PField_stokes,pidx,(void**)&mpprop_stokes);
		
		/* Get marker types */
		region_idx = mpprop_std->phase;
		
		min_cutoff_l = rheology->eta_lower_cutoff[region_idx];
		if (min_cutoff_l < min_cutoff) { min_cutoff_l = min_cutoff;}
		
		max_cutoff_l = rheology->eta_upper_cutoff[region_idx];
		if (max_cutoff_l > max_cutoff) { max_cutoff_l = max_cutoff;}
		
		MPntPStokesGetField_eta_effective(mpprop_stokes,&eta_mp);
		
		if (eta_mp > max_cutoff_l) { eta_mp = max_cutoff_l; npoints_cutoff[1]++; }
		if (eta_mp < min_cutoff_l) { eta_mp = min_cutoff_l; npoints_cutoff[0]++; }
		
		/* update viscosity on marker */
		MPntPStokesSetField_eta_effective(mpprop_stokes,eta_mp);
		
		/* monitor bounds */
		if (eta_mp > max_eta) { max_eta = eta_mp; }
		if (eta_mp < min_eta) { min_eta = eta_mp; }
	}  
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	
	ierr = MPI_Allreduce(&min_eta,&min_eta_g,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&max_eta,&max_eta_g,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(npoints_cutoff,npoints_cutoff_g,2, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	PetscTime(&t1);
	
  if ((npoints_cutoff_g[0] > 0) || (npoints_cutoff_g[1] > 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Apply viscosity cutoff (VPSTD) [mpoint]: (min,max)_eta %1.2e,%1.2e; log10(max/min) %1.2e; npoints_cutoff (%ld,%ld); cpu time %1.2e (sec)\n", min_eta_g, max_eta_g, log10(max_eta_g/min_eta_g), npoints_cutoff_g[0],npoints_cutoff_g[1], t1-t0 );
  }
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "StokesCoefficient_UpdateTimeDependentQuantities_VPSTD"
PetscErrorCode StokesCoefficient_UpdateTimeDependentQuantities_VPSTD(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[])
{
	PetscErrorCode ierr;
	DM             cda;
	Vec            gcoords;
	PetscScalar    *LA_gcoords;
	int            pidx,n_mp_points;
	DataBucket     db;
	DataField      PField_std,PField;
	float          strain_mp;
	PetscInt       nel,nen_u,k;
	const PetscInt *elnidx_u;
	PetscReal      elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal      elu[3*Q2_NODES_PER_EL_3D];
	PetscReal      ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscReal      GNI[3][Q2_NODES_PER_EL_3D];
	PetscReal      dNudx[Q2_NODES_PER_EL_3D],dNudy[Q2_NODES_PER_EL_3D],dNudz[Q2_NODES_PER_EL_3D];
	int            eidx_mp;
	double         *xi_mp;
	double         D_mp[NSD][NSD];
	double         inv2_D_mp;
    short          yield_type;
	PetscReal      dt;
	PetscFunctionBegin;
	
	/* access current time step */
	ierr = pTatinGetTimestep(user,&dt);CHKERRQ(ierr);
	
	/* access material point information */
	ierr = pTatinGetMaterialPoints(user,&db,NULL);CHKERRQ(ierr);
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataBucketGetDataFieldByName(db,MPntPStokesPl_classname,&PField);
	DataFieldGetAccess(PField);
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	/* setup for elements */
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	/* marker loop */
	for (pidx=0; pidx<n_mp_points; pidx++) {
		MPntStd       *mpprop_std;
		MPntPStokesPl *mpprop;
		
		DataFieldAccessPoint(PField_std, pidx,(void**)&mpprop_std);
		DataFieldAccessPoint(PField,     pidx,(void**)&mpprop);
		
        
        MPntPStokesPlGetField_yield_indicator(mpprop,&yield_type);
        if (yield_type > 0) {
            eidx_mp = mpprop_std->wil;
            
            /* Get element coordinates */
            ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*eidx_mp],LA_gcoords);CHKERRQ(ierr);
            /* Get element velocity */
            ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*eidx_mp],ufield);CHKERRQ(ierr);
            /* get velocity components */
            for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
                ux[k] = elu[3*k  ];
                uy[k] = elu[3*k+1];
                uz[k] = elu[3*k+2];
            }
            
            /* Get local coordinate of marker */
            xi_mp = mpprop_std->xi;
            
            /* Prepare basis functions */
            /* grad.Ni */
            P3D_ConstructGNi_Q2_3D(xi_mp,GNI);
            /* Get shape function derivatives */
            P3D_evaluate_global_derivatives_Q2(elcoords,GNI,dNudx,dNudy,dNudz);
            
            /* strain rate */
            ComputeStrainRate3d(ux,uy,uz,dNudx,dNudy,dNudz,D_mp);
            /* second inv stress */
            ComputeSecondInvariant3d(D_mp,&inv2_D_mp);
            
            MPntPStokesPlGetField_plastic_strain(mpprop,&strain_mp);
            strain_mp = strain_mp + dt * inv2_D_mp;
            MPntPStokesPlSetField_plastic_strain(mpprop,strain_mp);
        }
	}  
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField);
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}



#include "petsc.h"
#include "petscdm.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "swarm_fields.h"
#include "element_type_Q2.h"
#include "dmda_element_q2p1.h"
#include "element_utils_q2.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"

#include "material_constants.h"

typedef enum { YTYPE_NONE=0, YTYPE_MISES=1, YTYPE_DP=2, YTYPE_TENSILE_FAILURE=3 } YieldTypeDefinition;



static inline void ComputeStressIsotropic3d(PetscReal eta,double D[NSD][NSD],double T[NSD][NSD])
{
	const double two_eta = 2.0 * eta;
	
	T[0][0] = two_eta * D[0][0];	T[0][1] = two_eta * D[0][1];		T[0][2] = two_eta * D[0][2];
	T[1][0] = T[0][1];						T[1][1] = two_eta * D[1][1];		T[1][2] = two_eta * D[1][2];
	T[2][0] = T[0][2];						T[2][1] = T[1][2];							T[2][2] = two_eta * D[2][2];	
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
	D[1][0] = exy;		D[1][1] = eyy;		D[1][2] = exy;
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

void P3D_evaluate_global_derivatives_Q2(PetscReal el_coords[NPE*3],PetscReal GNI[3][NPE],
																		 PetscReal dNudx[NPE],
																		 PetscReal dNudy[NPE],
																		 PetscReal dNudz[NPE] )
{
	PetscInt k,p;
	PetscReal t4, t6, t8, t10, t12, t14, t17;
	PetscReal J[3][3],iJ[3][3];
	
	J[0][0] = J[0][1] = J[0][2] = 0.0;
	J[1][0] = J[1][1] = J[1][2] = 0.0;
	J[2][0] = J[2][1] = J[2][2] = 0.0;

	for (k=0; k<NPE; k++) {
		PetscReal xc = el_coords[3*k+0];
		PetscReal yc = el_coords[3*k+1];
		PetscReal zc = el_coords[3*k+2];
		
		J[0][0] += GNI[0][k] * xc;
		J[0][1] += GNI[0][k] * yc;
		J[0][2] += GNI[0][k] * zc;
		
		J[1][0] += GNI[1][k] * xc;
		J[1][1] += GNI[1][k] * yc;
		J[1][2] += GNI[1][k] * zc;
		
		J[2][0] += GNI[2][k] * xc;
		J[2][1] += GNI[2][k] * yc;
		J[2][2] += GNI[2][k] * zc;
	}
	/* flops = [NPE] * 18 */
	
	t4  = J[2][0] * J[0][1];
	t6  = J[2][0] * J[0][2];
	t8  = J[1][0] * J[0][1];
	t10 = J[1][0] * J[0][2];
	t12 = J[0][0] * J[1][1];
	t14 = J[0][0] * J[1][2]; // 6
	t17 = 0.1e1 / (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1]);  // 12
	
	iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17;  // 4
	iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17; // 5
	iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17;  // 4
	iJ[1][0] = -(-J[2][0] * J[1][2] + J[1][0] * J[2][2]) * t17;// 6
	iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17;                // 4
	iJ[1][2] = -(-t10 + t14) * t17;                            // 4
	iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17; // 5
	iJ[2][1] = -(-t4 + J[0][0] * J[2][1]) * t17;               // 5
	iJ[2][2] = (-t8 + t12) * t17;                              // 3
	/* flops = [NQP] * 58 */
	
	/* shape function derivatives */
	for (k=0; k<NPE; k++) {
		dNudx[k] = iJ[0][0]*GNI[0][k] + iJ[0][1]*GNI[1][k] + iJ[0][2]*GNI[2][k];
		
		dNudy[k] = iJ[1][0]*GNI[0][k] + iJ[1][1]*GNI[1][k] + iJ[1][2]*GNI[2][k];
		
		dNudz[k] = iJ[2][0]*GNI[0][k] + iJ[2][1]*GNI[1][k] + iJ[2][2]*GNI[2][k];
	}
	/* flops = [NPE] * 15 */
	
	// TOTAL = 18 + 58 + 15 = 105
}


#undef __FUNCT__
#define __FUNCT__ "EvaluateRheologyNonlinearitiesMarkers_VPSTD"
PetscErrorCode EvaluateRheologyNonlinearitiesMarkers_VPSTD(pTatinCtx user,DM dau,PetscScalar ufield[],DM dap,PetscScalar pfield[])
{
	PetscErrorCode ierr;
	
	int            pidx,n_mp_points;
	DataBucket     db,material_constants;
	DataField      PField_std,PField_stokes;
	PetscScalar    min_eta,max_eta,min_eta_g,max_eta_g;
	PetscLogDouble t0,t1;

	DM             cda;
	Vec            gcoords;
	PetscReal      *LA_gcoords;
	PetscInt       nel,nen_u,nen_p,eidx,k;
	const PetscInt *elnidx_u;
	const PetscInt *elnidx_p;
	PetscInt       vel_el_lidx[3*U_BASIS_FUNCTIONS];
	PetscReal      elcoords[3*Q2_NODES_PER_EL_3D];
	PetscReal      elu[3*Q2_NODES_PER_EL_3D],elp[P_BASIS_FUNCTIONS];
	PetscReal      ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
	PetscInt       *gidx;

	PetscReal detJ,NI[Q2_NODES_PER_EL_3D],GNI[3][Q2_NODES_PER_EL_3D],NIp[P_BASIS_FUNCTIONS];
	PetscReal dNudx[Q2_NODES_PER_EL_3D],dNudy[Q2_NODES_PER_EL_3D],dNudz[Q2_NODES_PER_EL_3D];
	
	double         eta_mp;
	double         D_mp[NSD][NSD],Tpred_mp[NSD][NSD];
	double         inv2_D_mp,inv2_Tpred_mp;
	
	DataField      PField_MatTypes;
	DataField      PField_ViscConst;
	DataField      PField_PlasticMises,PField_PlasticDP;
	MaterialConst_MaterialType   *MatType_data;
	MaterialConst_ViscosityConst *ViscConst_data;
	MaterialConst_PlasticMises   *PlasticMises_data;
	int            viscous_type,plastic_type,softening_type;
	int            npoints_yielded,npoints_yielded_g;
	
	
	PetscFunctionBegin;
	
	PetscGetTime(&t0);

	/* access material point information */
	ierr = pTatinGetMaterialPoints(user,&db,PETSC_NULL);CHKERRQ(ierr);
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dau,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dau,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

	/* get u,p element information */
	ierr = DMDAGetGlobalIndices(dau,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
	
	/* access material constants */
	ierr = pTatinGetMaterialConstants(user,&material_constants);CHKERRQ(ierr);

	DataBucketGetDataFieldByName(material_constants,MaterialConst_MaterialType_classname,  &PField_MatTypes);
	DataBucketGetDataFieldByName(material_constants,MaterialConst_ViscosityConst_classname,&PField_ViscConst);
	DataBucketGetDataFieldByName(material_constants,MaterialConst_PlasticMises_classname,  &PField_PlasticMises);
	MatType_data      = (MaterialConst_MaterialType*)  PField_MatTypes->data;
	ViscConst_data    = (MaterialConst_ViscosityConst*)PField_ViscConst->data;
	PlasticMises_data = (MaterialConst_PlasticMises*)  PField_PlasticMises->data;
	
	
	/* marker loop */
	min_eta = 1.0e100;
	max_eta = 1.0e-100;
	npoints_yielded = 0;
	for (pidx=0; pidx<n_mp_points; pidx++) {
		MPntStd     *mpprop_std;
		MPntPStokes *mpprop_stokes;
		double      *xi_p;
		double      pressure_mp;
		int         region_idx;
		
		DataFieldAccessPoint(PField_std,   pidx,(void**)&mpprop_std);
		DataFieldAccessPoint(PField_stokes,pidx,(void**)&mpprop_stokes);

		/* Get marker types */
		region_idx = mpprop_std->phase;
		
		viscous_type = MatType_data[ region_idx ].visc_type;
		plastic_type = MatType_data[ region_idx ].plastic_type;
		
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
				
		/* get viscosity on marker */
		//MPntPStokesGetField_eta_effective(mpprop_stokes,&eta_mp);		
		switch (viscous_type) {
				
			case VISCOUS_CONSTANT: {
				eta_mp = ViscConst_data[ region_idx ].eta0;
			}
				break;
				
			case VISCOUS_FRANKK: {

			}
				break;
				
			case VISCOUS_ARRHENIUS: {

			}
				break;
		}
		
		switch (plastic_type) {
				
			case PLASTIC_NONE: {
			}
				break;
				
			case PLASTIC_MISES: {
				double tau_yield_mp;
				
				tau_yield_mp = PlasticMises_data[ region_idx ].tau_yield;
				
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
				}
				break;
			}
				
			case PLASTIC_DP: {
				
			}
				break;
		}
		
		
		
		/* Global cutoffs for viscosity */
		/* Here should I store these? */
		
		
		/* update viscosity on marker */
		MPntPStokesSetField_eta_effective(mpprop_stokes,eta_mp);

		
		
		/* monitor bounds */
		if (eta_mp > max_eta) { max_eta = eta_mp; }
		if (eta_mp < min_eta) { min_eta = eta_mp; }
		
	}
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

  ierr = MPI_Allreduce(&min_eta,&min_eta_g,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&max_eta,&max_eta_g,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&npoints_yielded,&npoints_yielded_g,1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);

	PetscGetTime(&t1);
/*	
	PetscPrintf(PETSC_COMM_WORLD,"Update non-linearities (VPSTD) [mpoint]: (min,max)_eta %1.2e,%1.2e; log10(max/min) %1.2e; npoints_yielded %d; cpu time %1.2e (sec)\n",
              min_eta_g, max_eta_g, log10(max_eta_g/min_eta_g), npoints_yielded_g, t1-t0 );
*/
	
	PetscFunctionReturn(0);
}


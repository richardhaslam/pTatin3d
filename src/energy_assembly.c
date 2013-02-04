
#include "petsc.h"
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_bcs.h"
#include "element_utils_q1.h"
#include "dmda_element_q1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"

#include "QPntVolCoefEnergy_def.h"
#include "phys_comp_energy.h"
#include "ptatin3d_energy.h"

#undef __FUNCT__
#define __FUNCT__ "AdvDiff_GetElementLocalIndices_Q1"
PetscErrorCode AdvDiff_GetElementLocalIndices_Q1(PetscInt el_localIndices[],PetscInt elnid[])
{
	PetscInt n;
	PetscFunctionBegin;
	for (n=0; n<8; n++) {
		el_localIndices[n] = elnid[n];
	}
	PetscFunctionReturn(0);
}


/* SUPG business */
/** AdvectionDiffusion_UpwindXiExact - Brooks, Hughes 1982 equation 2.4.2
 * \bar \xi = coth( \alpha ) - \frac{1}{\alpha} */
double AdvDiffResidualForceTerm_UpwindXiExact(double pecletNumber)
{
  if (fabs(pecletNumber) < 1.0e-8 ) {
    return 0.33333333333333 * pecletNumber;
  } else if (pecletNumber < -20.0) {
    return -1.0 - 1.0/pecletNumber;
  } else if (pecletNumber > 20.0) {
    return +1.0 - 1.0/pecletNumber;
  }
  return cosh( pecletNumber )/sinh( pecletNumber ) - 1.0/pecletNumber;
}

/** AdvectionDiffusion_UpwindXiDoublyAsymptoticAssumption - Brooks, Hughes 1982 equation 3.3.1
 * Simplification of \bar \xi = coth( \alpha ) - \frac{1}{\alpha} from Brooks, Hughes 1982 equation 2.4.2
 * 
 *            { -1               for \alpha <= -3
 * \bar \xi ~ { \frac{/alpha}{3} for -3 < \alpha <= 3
 *            { +1               for \alpha > +3 */
double AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption(double pecletNumber)
{
  if (pecletNumber <= -3.0) {
    return -1;
  } else if (pecletNumber <= 3.0) {
    return 0.33333333333333 * pecletNumber;
  } else {
    return 1.0;
  }
}

/** AdvectionDiffusion_UpwindXiCriticalAssumption - Brooks, Hughes 1982 equation 3.3.2
 * Simplification of \bar \xi = coth( \alpha ) - \frac{1}{\alpha} from Brooks, Hughes 1982 equation 2.4.2
 * 
 *            { -1 - \frac{1}{\alpha}   for \alpha <= -1
 * \bar \xi ~ { 0                       for -1 < \alpha <= +1
 *            { +1 - \frac{1}{\alpha}   for \alpha > +1             */
double AdvDiffResidualForceTerm_UpwindXiCriticalAssumption(double pecletNumber)
{
  if (pecletNumber <= -1.0) {
    return -1.0 - 1.0/pecletNumber;
  } else if (pecletNumber <= 1.0) {
    return 0.0;
  } else {
    return 1.0 - 1.0/pecletNumber;
  }
}

#undef __FUNCT__
#define __FUNCT__ "DASUPG3dComputeAverageCellSize"
PetscErrorCode DASUPG3dComputeAverageCellSize(PetscScalar el_coords[],PetscScalar DX[])
{
	PetscInt d,k;
	PetscReal min_x[NSD],max_x[NSD];
	PetscReal xc,yc,zc;
	
	PetscFunctionBegin;
	for (d=0; d<NSD; d++) {
		min_x[d] = 1.0e32;
		max_x[d] = -1.0e32;
	}
	
	for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		xc = el_coords[NSD*k+0];
		yc = el_coords[NSD*k+1];
		zc = el_coords[NSD*k+2];
		
		if (xc < min_x[0]) { min_x[0] = xc; }
		if (yc < min_x[1]) { min_x[1] = yc; }
		if (zc < min_x[2]) { min_x[2] = zc; }
		
		if (xc > max_x[0]) { max_x[0] = xc; }
		if (yc > max_x[1]) { max_x[1] = yc; }
		if (zc > max_x[2]) { max_x[2] = zc; }
	}
	
	for (d=0; d<NSD; d++) {
		DX[d] = max_x[d] - min_x[d];
	}
	PetscFunctionReturn(0);
}

/* Eqn 4.3.7 */
#undef __FUNCT__
#define __FUNCT__ "DASUPG3dComputeElementPecletNumber_qp"
PetscErrorCode DASUPG3dComputeElementPecletNumber_qp( PetscScalar el_coords[],PetscScalar u[],
																										 PetscScalar kappa_el,
																										 PetscScalar *alpha)
{
	PetscInt    d,k;
	PetscScalar DX[NSD];
  PetscScalar u_xi[NSD],one_dxi2[NSD];
  PetscScalar _alpha,u_norm,dxi[NSD];
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
	
	ierr = DASUPG3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
	
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
  for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] = 0.125 *  u[NSD*k + 0];
		u_xi[1] = 0.125 *  u[NSD*k + 1];
		u_xi[2] = 0.125 *  u[NSD*k + 2];
	}
	
	for (d=0; d<NSD; d++) {
		dxi[d]      = DX[d];
		one_dxi2[d] = 1.0/dxi[d]/dxi[d];
	}
	
  u_norm = sqrt( u_xi[0]*u_xi[0] + u_xi[1]*u_xi[1] + u_xi[2]*u_xi[2]);
	
  _alpha = 0.5 * u_norm / ( (1.0e-32+kappa_el) * sqrt( one_dxi2[0] + one_dxi2[1] + one_dxi2[2]) );
  *alpha  = _alpha;  
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DASUPG3dComputeElementTimestep_qp"
PetscErrorCode DASUPG3dComputeElementTimestep_qp(PetscScalar el_coords[],PetscScalar u[],PetscScalar kappa_el,PetscScalar *dta,PetscScalar *dtd)
{
	PetscInt    d,k;
	PetscScalar DX[NSD];
  PetscScalar u_xi[NSD],one_dxi2[NSD],dxi[NSD];
  PetscScalar CrFAC,dt_optimal,alpha;
  PetscScalar U,H,dt_diffusive;
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
	
	ierr = DASUPG3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
  
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
  for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] = 0.125 *  u[NSD*k + 0];
		u_xi[1] = 0.125 *  u[NSD*k + 1];
		u_xi[2] = 0.125 *  u[NSD*k + 2];
	}
	
	for (d=0; d<NSD; d++) {
		dxi[d]      = DX[d];
		one_dxi2[d] = 1.0/dxi[d]/dxi[d];
	}
	
  U = sqrt( u_xi[0]*u_xi[0] + u_xi[1]*u_xi[1] + u_xi[2]*u_xi[2] );
  H = 1.0 / sqrt( one_dxi2[0] + one_dxi2[1] + one_dxi2[2] );
	
  ierr = DASUPG3dComputeElementPecletNumber_qp(el_coords,u,kappa_el,&alpha);CHKERRQ(ierr);
	
  if (U<1.0e-32) {
    dt_diffusive = 0.5 * H*H / kappa_el;
    *dtd = dt_diffusive;
    *dta = 1.0e10;
  } else {
    if (alpha >= 100.0 ) {
      CrFAC = 0.4;
    } else {
      CrFAC = PetscMin(1.0,alpha);
    }
    dt_optimal = CrFAC * H / (U+1.0e-32);
    *dta = dt_optimal;
    *dtd = 1.0e10;
  }
	/*
	 printf("Element Pe:              %1.4f \n", alpha );
	 printf("        dx/dv:           %1.4f \n", H/(U+1.0e-32) );
	 printf("        0.5.dx.dx/kappa: %1.4f \n", 0.5*H*H/(kappa_el+1.0e-32) );
	 */
	
	*dta = H/(U+1.0e-32);
	*dtd = 0.5*H*H/(kappa_el+1.0e-32);
	
  PetscFunctionReturn(0);
}

/*
 Eqns 3.3.11 (transient) + 3.3.4, 3.3.5, 3.3.6
 */
#undef __FUNCT__
#define __FUNCT__ "DASUPG3dComputeElementStreamlineDiffusion_qp"
PetscErrorCode DASUPG3dComputeElementStreamlineDiffusion_qp(PetscScalar el_coords[],PetscScalar u[],
																														PetscInt nqp, PetscScalar qp_detJ[], PetscScalar qp_w[],
																														PetscScalar qp_kappa[],
																														PetscScalar *khat)
{
	PetscInt p,k,d;
	PetscScalar DX[NSD];
  PetscScalar u_xi[NSD];
  PetscScalar kappa_el,alpha_xi[NSD],_khat,dxi[NSD],xi[NSD],vol_el;
	PetscErrorCode ierr;
	
  PetscFunctionBegin;
	
	/* average velocity - cell centre velocity */
	ierr = DASUPG3dComputeAverageCellSize(el_coords,DX);CHKERRQ(ierr);
	
	/* average velocity - cell centre velocity */
	u_xi[0] = u_xi[1] = u_xi[2] = 0.0;
  for (k=0; k<NODES_PER_EL_Q1_3D; k++) {
		u_xi[0] = 0.125 *  u[NSD*k + 0];
		u_xi[1] = 0.125 *  u[NSD*k + 1];
		u_xi[2] = 0.125 *  u[NSD*k + 2];
	}
	
	/* integral average of diffusivity */
	kappa_el = 0.0;
	vol_el   = 0.0;
	for (p=0; p<nqp; p++) {
		kappa_el += qp_w[p] * qp_kappa[p] * qp_detJ[p]; 
		vol_el   += qp_w[p] * qp_detJ[p];
	}
  kappa_el = kappa_el / vol_el;
	
	for (d=0; d<NSD; d++) {
		dxi[d]      = DX[d];
		alpha_xi[d] = 0.5 * u_xi[d]  * dxi[d]  / kappa_el;
		xi[d]       = AdvDiffResidualForceTerm_UpwindXiExact(alpha_xi[d]);
	}
	
	//  khat = 0.5*(xi*u_xi*dxi + eta*v_eta*deta); /* steady state */
  _khat = 1.0 * 0.258198889747161 * ( xi[0]*u_xi[0]*dxi[0] + xi[1]*u_xi[1]*dxi[1] + xi[2]*u_xi[2]*dxi[2] ); /* transient case, sqrt(1/15) */
  *khat = _khat;
	
  PetscFunctionReturn(0);
}

/* SUPG IMPLEMENTATION */
void ConstructNiSUPG_Q1_3D(PetscScalar Up[],PetscScalar kappa_hat,PetscScalar Ni[],PetscScalar GNx[NSD][NODES_PER_EL_Q1_3D],PetscScalar Ni_supg[])
{
  PetscScalar uhat[NSD],unorm;
  PetscInt i;
	
  unorm = PetscSqrtScalar(Up[0]*Up[0] + Up[1]*Up[1] + Up[2]*Up[2]);
  uhat[0] = Up[0]/unorm;
  uhat[1] = Up[1]/unorm;
  uhat[2] = Up[2]/unorm;
	
  if (kappa_hat < 1.0e-30) {
    for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
      Ni_supg[i] = Ni[i];
    }
  } else {
    for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
      Ni_supg[i] = Ni[i] + kappa_hat * ( uhat[0] * GNx[0][i] + uhat[1] * GNx[1][i] + uhat[2] * GNx[2][i] ) / unorm;
    }
  }
}

#undef __FUNCT__  
#define __FUNCT__ "AElement_FormJacobian_T"
PetscErrorCode AElement_FormJacobian_T( PetscScalar Re[],PetscReal dt,PetscScalar el_coords[],
																			 PetscScalar gp_kappa[],
																			 PetscScalar el_V[],
																			 PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] )
{
  PetscInt    p,i,j;
  PetscReal   Ni_p[NODES_PER_EL_Q1_3D],Ni_supg_p[NODES_PER_EL_Q1_3D];
  PetscReal   GNi_p[NSD][NODES_PER_EL_Q1_3D],GNx_p[NSD][NODES_PER_EL_Q1_3D],gp_detJ[27];
  PetscScalar J_p,fac;
  PetscScalar kappa_p,v_p[NSD];
  PetscScalar kappa_hat;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
  /* compute constants for the element */
  for (p = 0; p < ngp; p++) {
    P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
    P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&gp_detJ[p],&GNx_p[0],&GNx_p[1],&GNx_p[2]);
	}
	
	ierr = DASUPG3dComputeElementStreamlineDiffusion_qp(el_coords,el_V,
																							        ngp,gp_detJ,gp_weight,gp_kappa,
																							        &kappa_hat);CHKERRQ(ierr);	
  /* evaluate integral */
  for (p = 0; p < ngp; p++) {
    P3D_ConstructNi_Q1_3D(&gp_xi[NSD*p],Ni_p);
    P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
    P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&J_p,&GNx_p[0],&GNx_p[1],&GNx_p[2]);
		
    fac = gp_weight[p]*J_p;
		
		kappa_p = gp_kappa[p];
    v_p[0] = 0.0;
		v_p[1] = 0.0;
		v_p[2] = 0.0;
    for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
      v_p[0]     += Ni_p[j] * el_V[NSD*j+0];      /* compute vx on the particle */
      v_p[1]     += Ni_p[j] * el_V[NSD*j+1];      /* compute vy on the particle */
      v_p[2]     += Ni_p[j] * el_V[NSD*j+2];      /* compute vy on the particle */
    }
    ConstructNiSUPG_Q1_3D(v_p,kappa_hat,Ni_p,GNx_p,Ni_supg_p);
		
    /* R = f - m phi_dot - c phi*/
		for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
			for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
				Re[j+i*NODES_PER_EL_Q1_3D] += fac * ( 
																						 //Ni_p[i] * Ni_p[j]
																						 Ni_supg_p[i] * Ni_p[j]
																						 + dt * Ni_supg_p[i] * ( v_p[0] * GNx_p[0][j] + v_p[1] * GNx_p[1][j] + v_p[2] * GNx_p[2][j] )
																						 + dt * kappa_p * ( GNx_p[0][i] * GNx_p[0][j] + GNx_p[1][i] * GNx_p[1][j] + GNx_p[2][i] * GNx_p[2][j] )
																						 );
			}
		}
  }
	
	 printf("e=\n");
	 for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
	 for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
	 printf("%lf ", Re[j+i*NODES_PER_EL_Q1_3D]); 
	 }printf("\n");
	 }
	 
	PetscFunctionReturn(0);
}

/* 
 Computes M + dt.(L + A)
 */
#undef __FUNCT__  
#define __FUNCT__ "FormJacobianEnergy"
PetscErrorCode FormJacobianEnergy(PetscReal time,Vec X,PetscReal dt,Mat *A,Mat *B,MatStructure *mstr,void *ctx)
{
  PhysCompEnergy data = (PhysCompEnergy)ctx;
	PetscInt          nqp;
	PetscScalar       *qp_xi,*qp_weight;
	Quadrature        volQ;
	QPntVolCoefEnergy *all_quadpoints,*cell_quadpoints;
	PetscScalar       qp_kappa[27];
  DM            da,cda;
	Vec           gcoords;
	PetscScalar   *LA_gcoords;
	PetscScalar   ADe[NODES_PER_EL_Q1_3D*NODES_PER_EL_Q1_3D];
	PetscScalar   el_coords[NSD*NODES_PER_EL_Q1_3D];
	PetscScalar   el_V[NSD*NODES_PER_EL_Q1_3D];
	/**/
	PetscInt       nel,nen,e,n;
	const PetscInt *elnidx;
	BCList         bclist;
	PetscInt       NUM_GINDICES,*GINDICES,ge_eqnums[Q1_NODES_PER_EL_3D];
	Vec            V;
  Vec            local_V;
  PetscScalar    *LA_V;
	PetscErrorCode ierr;
	
	
  PetscFunctionBegin;
	da     = data->daT;
	V      = data->u_minus_V;
	bclist = data->T_bclist;
	volQ   = data->volQ;
	
	/* trash old entries */
  ierr = MatZeroEntries(*B);CHKERRQ(ierr);
	
	/* quadrature */
	volQ      = data->volQ;
	nqp       = volQ->npoints;
	qp_xi     = volQ->q_xi_coor;
	qp_weight = volQ->q_weight;
	
	ierr = VolumeQuadratureGetAllCellData_Energy(volQ,&all_quadpoints);CHKERRQ(ierr);
	
  
	/* setup for coords */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
  /* get acces to the vector V */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMGetLocalVector(cda,&local_V);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(cda,V,INSERT_VALUES,local_V);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(  cda,V,INSERT_VALUES,local_V);CHKERRQ(ierr);
  ierr = VecGetArray(local_V,&LA_V);CHKERRQ(ierr);
	
	/* stuff for eqnums */
	ierr = DMDAGetGlobalIndices(da,&NUM_GINDICES,&GINDICES);CHKERRQ(ierr);
	ierr = BCListApplyDirichletMask(NUM_GINDICES,GINDICES,bclist);CHKERRQ(ierr);
	
	ierr = DMDAGetElementsQ1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		/* get coords for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		
		/* get velocity for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_V,(PetscInt*)&elnidx[nen*e],LA_V);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Energy(volQ,all_quadpoints,e,&cell_quadpoints);CHKERRQ(ierr);
		
		/* copy the diffusivity */
		for (n=0; n<nqp; n++) {
			qp_kappa[n] = cell_quadpoints[n].diffusivity;
		}
		
		/* initialise element stiffness matrix */
		ierr = PetscMemzero(ADe,sizeof(PetscScalar)*NODES_PER_EL_Q1_3D*NODES_PER_EL_Q1_3D);CHKERRQ(ierr);
		
		/* form element stiffness matrix */
		ierr = AElement_FormJacobian_T( ADe,dt,el_coords, qp_kappa, el_V, nqp,qp_xi,qp_weight );CHKERRQ(ierr);
		
		/* insert element matrix into global matrix */
		ierr = DMDAEQ1_GetElementLocalIndicesDOF(ge_eqnums,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		ierr = MatSetValues(*B,NODES_PER_EL_Q1_3D,ge_eqnums, NODES_PER_EL_Q1_3D,ge_eqnums, ADe, ADD_VALUES );CHKERRQ(ierr);
  }
	/* tidy up */
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(local_V,&LA_V);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(cda,&local_V);CHKERRQ(ierr);
	
	/* partial assembly */
	ierr = MatAssemblyBegin(*B, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*B, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = BCListRemoveDirichletMask(NUM_GINDICES,GINDICES,bclist);CHKERRQ(ierr);
	ierr = BCListInsertScaling(*B,NUM_GINDICES,GINDICES,bclist);CHKERRQ(ierr);
	
	/* assemble */
	ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (*A != *B) {
    ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
	*mstr = SAME_NONZERO_PATTERN;
  PetscFunctionReturn(0);
}

void AElement_FormFunctionLocal_SUPG_T(
																			 PetscScalar Re[],
																			 PetscReal dt,
																			 PetscScalar el_coords[],PetscScalar el_V[],
																			 PetscScalar el_phi[],PetscScalar el_phi_last[],
																			 PetscScalar gp_kappa[],PetscScalar gp_Q[],
																			 PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] )
{
  PetscInt    p,i,j;
  PetscScalar Ni_p[NODES_PER_EL_Q1_3D],Ni_supg_p[NODES_PER_EL_Q1_3D];
  PetscScalar GNi_p[NSD][NODES_PER_EL_Q1_3D],GNx_p[NSD][NODES_PER_EL_Q1_3D];
  PetscScalar phi_p,phi_last_p,f_p,v_p[NSD],kappa_p, gradphi_p[NSD],gradphiold_p[NSD],M_dotT_p;
  PetscScalar J_p,fac,gp_detJ[27];
  PetscScalar kappa_hat;
	
  /* compute constants for the element */
  for (p = 0; p < ngp; p++) {
		P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
		P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&gp_detJ[p],&GNx_p[0],&GNx_p[1],&GNx_p[2]);
	}
	
	DASUPG3dComputeElementStreamlineDiffusion_qp(el_coords,el_V,
																							 ngp,gp_detJ,gp_weight,gp_kappa,
																							 &kappa_hat);
	
	
  /* evaluate integral */
  for (p = 0; p < ngp; p++) {
    P3D_ConstructNi_Q1_3D(&gp_xi[NSD*p],Ni_p);
		P3D_ConstructGNi_Q1_3D(&gp_xi[NSD*p],GNi_p);
		P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&J_p,&GNx_p[0],&GNx_p[1],&GNx_p[2]);
		
    fac = gp_weight[p]*J_p;
		
		kappa_p      = gp_kappa[p];
		f_p          = gp_Q[p];
    phi_p        = 0.0;
    phi_last_p   = 0.0;
    gradphi_p[0] = 0.0;
		gradphi_p[1] = 0.0;
		gradphi_p[2] = 0.0;
    gradphiold_p[0] = 0.0;
		gradphiold_p[1] = 0.0;
		gradphiold_p[2] = 0.0;
    v_p[0] = 0.0;
		v_p[1] = 0.0;
		v_p[2] = 0.0;
    for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
      phi_p        += Ni_p[j] * el_phi[j];      /* compute phi on the particle */
      phi_last_p   += Ni_p[j] * el_phi_last[j];  /* compute phi_dot on the particle */
      gradphiold_p[0] += GNx_p[0][j] * el_phi_last[j];
      gradphiold_p[1] += GNx_p[1][j] * el_phi_last[j];
      gradphiold_p[2] += GNx_p[2][j] * el_phi_last[j];
			
      gradphi_p[0] += GNx_p[0][j] * el_phi[j];
      gradphi_p[1] += GNx_p[1][j] * el_phi[j];
      gradphi_p[2] += GNx_p[2][j] * el_phi[j];
			
      v_p[0]     += Ni_p[j] * el_V[NSD*j+0];      /* compute vx on the particle */
      v_p[1]     += Ni_p[j] * el_V[NSD*j+1];      /* compute vy on the particle */
      v_p[2]     += Ni_p[j] * el_V[NSD*j+2];      /* compute vz on the particle */
    }
    ConstructNiSUPG_Q1_3D(v_p,kappa_hat,Ni_p,GNx_p,Ni_supg_p);
		
		
		M_dotT_p = 0.0;
		for (j=0; j<NODES_PER_EL_Q1_3D; j++) {
			M_dotT_p += Ni_supg_p[j] * el_phi_last[j];
		}
		
    /* R = f - m phi_dot - c phi*/
    for (i = 0; i < NODES_PER_EL_Q1_3D; i++) {
      Re[ i ] += fac * (
												// -f - M T^k
												- dt * Ni_p[i] * f_p 
												// (C+L) T^k+1
												+ dt * kappa_p * ( GNx_p[0][i]*gradphi_p[0] + GNx_p[1][i]*gradphi_p[1] + GNx_p[2][i]*gradphi_p[2] )
												+ dt * Ni_supg_p[i] * ( v_p[0]*gradphi_p[0] + v_p[1]*gradphi_p[1] + v_p[2]*gradphi_p[2] ) 
												// M T^k+1
                        + Ni_supg_p[i] * phi_p 
                        - Ni_supg_p[i] * phi_last_p 
                        //+ Ni_p[i] * phi_p 
                        //- Ni_p[i] * phi_last_p 
												);
    }
  }
	
}

#undef __FUNCT__  
#define __FUNCT__ "FormFunctionLocal_SUPG_T"
PetscErrorCode FormFunctionLocal_SUPG_T(
																				PhysCompEnergy data,
																				PetscReal dt,
																				DM da,
																				PetscScalar *LA_V,
																				PetscScalar *LA_phi,
																				PetscScalar *LA_philast,
																				PetscScalar *LA_R)
{
  DM                     cda;
  Vec                    gcoords;
  PetscScalar            *LA_gcoords;
  PetscScalar            Re[NODES_PER_EL_Q1_3D];
  PetscScalar            el_coords[NSD*NODES_PER_EL_Q1_3D];
	PetscScalar            el_V[NSD*NODES_PER_EL_Q1_3D];
	PetscScalar            el_phi[NODES_PER_EL_Q1_3D];
	PetscScalar            el_philast[NODES_PER_EL_Q1_3D];
	
	PetscInt       nel,nen,e,n;
	const PetscInt *elnidx;
	BCList         bclist;
	PetscInt       NUM_GINDICES,*GINDICES,ge_eqnums[NODES_PER_EL_Q1_3D];
	
	PetscInt          nqp;
	PetscScalar       *qp_xi,*qp_weight;
	Quadrature        volQ;
	QPntVolCoefEnergy *all_quadpoints,*cell_quadpoints;
	PetscScalar       qp_kappa[27],qp_Q[27];
	
  PetscScalar  Ni_p[NODES_PER_EL_Q1_3D];
  PetscScalar  GNi_p[NSD][NODES_PER_EL_Q1_3D],GNx_p[NSD][NODES_PER_EL_Q1_3D];
  PetscScalar  J_p,fac;
  PetscScalar  phi_p;
  PetscInt     p,i;
	PetscReal    cg,c = 0.0;
  PetscErrorCode   ierr;
	
  PetscFunctionBegin;
	//da     = data->daT;
	//V      = data->u_minus_V;
	bclist = data->T_bclist;
	volQ   = data->volQ;
	
	/* quadrature */
	nqp       = volQ->npoints;
	qp_xi     = volQ->q_xi_coor;
	qp_weight = volQ->q_weight;
	
	ierr = VolumeQuadratureGetAllCellData_Energy(volQ,&all_quadpoints);CHKERRQ(ierr);
	
	/* setup for coords */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	/* stuff for eqnums */
	ierr = DMDAGetGlobalIndices(da,&NUM_GINDICES,&GINDICES);CHKERRQ(ierr);
	
	ierr = DMDAGetElementsQ1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	for (e=0;e<nel;e++) {
		
		ierr = DMDAEQ1_GetElementLocalIndicesDOF(ge_eqnums,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		/* get coords for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAEQ1_GetScalarElementField_3D(el_phi,(PetscInt*)&elnidx[nen*e],LA_phi);CHKERRQ(ierr);
		ierr = DMDAEQ1_GetScalarElementField_3D(el_philast,(PetscInt*)&elnidx[nen*e],LA_philast);CHKERRQ(ierr);
		
		/* get velocity for the element */
		ierr = DMDAEQ1_GetVectorElementField_3D(el_V,(PetscInt*)&elnidx[nen*e],LA_V);CHKERRQ(ierr);
		
		ierr = VolumeQuadratureGetCellData_Energy(volQ,all_quadpoints,e,&cell_quadpoints);CHKERRQ(ierr);
		
		/* copy the diffusivity and force */
		for (n=0; n<nqp; n++) {
			qp_kappa[n] = cell_quadpoints[n].diffusivity;
			qp_Q[n]     = cell_quadpoints[n].heat_source;
		}
		
		/* initialise element stiffness matrix */
		ierr = PetscMemzero(Re,sizeof(PetscScalar)*NODES_PER_EL_Q1_3D);CHKERRQ(ierr);
		
		/* form element stiffness matrix */
		AElement_FormFunctionLocal_SUPG_T(Re,dt,el_coords,el_V,el_phi,el_philast,qp_kappa,qp_Q,nqp,qp_xi,qp_weight);
		
		ierr = DMDAEQ1_SetValuesLocalStencil_AddValues_DOF(LA_R,1,ge_eqnums,Re);CHKERRQ(ierr);
		
		/* diagnostics */
		for (p=0; p<nqp; p++) {
			P3D_ConstructNi_Q1_3D(&qp_xi[NSD*p],Ni_p);
			P3D_ConstructGNi_Q1_3D(&qp_xi[NSD*p],GNi_p);
			P3D_evaluate_geometry_elementQ1(1,el_coords,&GNi_p,&J_p,&GNx_p[0],&GNx_p[1],&GNx_p[2]);
			
			phi_p = 0.0;
			for (i=0; i<NODES_PER_EL_Q1_3D; i++) {
				phi_p = phi_p + Ni_p[i] * el_philast[i];
			}
			
			fac = qp_weight[p]*J_p;
			c = c + phi_p * fac;
		}
		
	}
  ierr = MPI_Allreduce(&c,&cg,1,MPIU_REAL,MPI_MIN,((PetscObject)da)->comm);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\\int \\phi dV = %1.12e \n", cg );
	
  /* tidy up local arrays (input) */
  ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "FormFunctionEnergy"
PetscErrorCode FormFunctionEnergy(PetscReal time,Vec X,PetscReal dt,Vec F,void *ctx)
{
  PhysCompEnergy data  = (PhysCompEnergy)ctx;
  DM             da,cda;
	Vec            philoc, philastloc, Fphiloc;
	Vec            Vloc;
	PetscScalar    *LA_philoc, *LA_philastloc, *LA_Fphiloc;
	PetscScalar    *LA_V;
  PetscErrorCode ierr;
	
  PetscFunctionBegin;
	da = data->daT;
	
	
	ierr = VecZeroEntries(F);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(da,&philoc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(da,&philastloc);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(da,&Fphiloc);CHKERRQ(ierr);
	
	/* get local solution and time derivative */
	ierr = VecZeroEntries(philoc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,X,ADD_VALUES,philoc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (da,X,ADD_VALUES,philoc);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(philastloc);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,data->Told,ADD_VALUES,philastloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (da,data->Told,ADD_VALUES,philastloc);CHKERRQ(ierr);
	
	/* insert boundary conditions into local vectors */
	ierr = BCListInsertLocal(data->T_bclist,philoc);CHKERRQ(ierr);
	
	/* what do i do with the phi_dot */
	
	/* init residual */
	ierr = VecZeroEntries(Fphiloc);CHKERRQ(ierr);
	
	/* get arrays */
	ierr = VecGetArray(philoc,    &LA_philoc);CHKERRQ(ierr);
	ierr = VecGetArray(philastloc,&LA_philastloc);CHKERRQ(ierr);
	ierr = VecGetArray(Fphiloc,   &LA_Fphiloc);CHKERRQ(ierr);
	
	/* ============= */
	/* FORM_FUNCTION */
	
  /* get acces to the vector V */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMGetLocalVector(cda,&Vloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(cda,data->u_minus_V,INSERT_VALUES,Vloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(  cda,data->u_minus_V,INSERT_VALUES,Vloc);CHKERRQ(ierr);
  ierr = VecGetArray(Vloc,&LA_V);CHKERRQ(ierr);
	
	ierr = FormFunctionLocal_SUPG_T(data,dt,da,LA_V, LA_philoc,LA_philastloc,LA_Fphiloc);CHKERRQ(ierr);
	
  ierr = VecRestoreArray(Vloc,&LA_V);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(cda,&Vloc);CHKERRQ(ierr);
	/* ============= */
	
	ierr = VecRestoreArray(Fphiloc,   &LA_Fphiloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(philastloc,&LA_philastloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(philoc,    &LA_philoc);CHKERRQ(ierr);
	
	
	/* do global fem summation */
	ierr = VecZeroEntries(F);CHKERRQ(ierr);
	
  ierr = DMLocalToGlobalBegin(da,Fphiloc,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd  (da,Fphiloc,ADD_VALUES,F);CHKERRQ(ierr);
	
  ierr = DMRestoreLocalVector(da,&Fphiloc);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&philoc);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&philastloc);CHKERRQ(ierr);
	
	/* modify F for the boundary conditions, F_k = scale_k(x_k - phi_k) */
	ierr = BCListResidualDirichlet(data->T_bclist,X,F);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}




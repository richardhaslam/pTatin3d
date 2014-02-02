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
 **    Filename:      ptatin_driver_linear_ts.c
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

static const char help[] = "Stokes solver using Q2-Pm1 mixed finite elements.\n"
"3D prototype of the (p)ragmatic version of Tatin. (pTatin3d_v0.0)\n\n";

#include "petsc-private/daimpl.h" 

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"
#include "ptatin_log.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "stokes_form_function.h"
#include "stokes_operators.h"
#include "stokes_operators_mf.h"
#include "stokes_assembly.h"
#include "dmda_element_q2p1.h"
#include "dmda_duplicate.h"
#include "dmda_redundant.h"
#include "dmda_project_coords.h"
#include "monitors.h"
#include "mp_advection.h"
#include "material_point_popcontrol.h"

#include "dmda_element_q2p1.h"
#include "element_utils_q2.h"
#include "quadrature.h"

typedef enum { OP_TYPE_REDISC_ASM=0, OP_TYPE_REDISC_MF, OP_TYPE_GALERKIN, OP_TYPE_MFGALERKIN } OperatorType;

#undef __FUNCT__  
#define __FUNCT__ "MatAssembleMFGalerkin"
PetscErrorCode MatAssembleMFGalerkin(DM dav_fine,BCList u_bclist_fine,Quadrature volQ_fine,DM dav_coarse,Mat Ac)
{
	PetscInt    refi,refj,refk,nni,nnj,nnk,n,sei,sej,sek;
	PetscInt    ic,jc,kc,iif,jjf,kkf;
	PetscInt    lmk,lmj,lmi;
	PetscInt    lmk_coarse,lmj_coarse,lmi_coarse;
	DM          daf,dac,cda;
	DMDABoundaryType wrap[] = { DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE };
	Mat         Acell,Ac_el,P;
	Vec              gcoords;
	PetscReal        *LA_gcoords;	
	PetscInt         q,nqp,ii,jj,kk;
	PetscReal         WEIGHT[NQP],XI[NQP][3];
	QPntVolCoefStokes *quadrature_points,*cell_quadrature_points;
	PetscReal      NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal      detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscInt       nel,nen_v;
	const PetscInt *elnidx_v;
	PetscInt       nel_cell,nen_v_cell;
	const PetscInt *elnidx_v_cell;
	PetscInt       nel_coarse,nen_v_coarse;
	const PetscInt *elnidx_v_coarse;
	PetscReal      el_coords[3*Q2_NODES_PER_EL_3D],el_eta[MAX_QUAD_PNTS];
	PetscReal      Ae[Q2_NODES_PER_EL_3D * Q2_NODES_PER_EL_3D * U_DOFS * U_DOFS];
	PetscReal      fac,D[NSTRESS][NSTRESS],diagD[NSTRESS],B[6][3*Q2_NODES_PER_EL_3D];
	PetscInt       NUM_GINDICES,*GINDICES,ge_eqnums[3*Q2_NODES_PER_EL_3D];
	PetscInt       NUM_GINDICES_cell,*GINDICES_cell,ge_eqnums_cell[3*Q2_NODES_PER_EL_3D];
	PetscInt       NUM_GINDICES_coarse,*GINDICES_coarse,ge_eqnums_coarse[3*Q2_NODES_PER_EL_3D];
	PetscScalar    *Ac_entries;
	PetscLogDouble t0,t1,t[6];
  PetscErrorCode ierr;

	
	ierr = MatZeroEntries(Ac);CHKERRQ(ierr);

	ierr = DMDAGetLocalSizeElementQ2(dav_fine,&lmi,&lmj,&lmk);CHKERRQ(ierr);
	ierr = DMDAGetLocalSizeElementQ2(dav_coarse,&lmi_coarse,&lmj_coarse,&lmk_coarse);CHKERRQ(ierr);
	
	ierr = DMDAGetRefinementFactor(dav_fine,&refi,&refj,&refk);CHKERRQ(ierr);

	nni = 2*refi + 1;
	nnj = 2*refj + 1;
	nnk = 2*refk + 1;
	
	PetscPrintf(PETSC_COMM_WORLD,"MatAssembleMFGalerkin:\n");
	PetscPrintf(PETSC_COMM_WORLD,"  Q2 cell problem contains: [%d x %d x %d] elements, [%d x %d x %d] nodes\n",refi,refj,refk,nni,nnj,nnk);
	
	ierr = x_DMDACreate3d(PETSC_COMM_SELF,wrap,DMDA_STENCIL_BOX,nni,nnj,nnk, 1,1,1, 3,2, 0,0,0,&daf);CHKERRQ(ierr);
	ierr = DMDASetElementType_Q2(daf);CHKERRQ(ierr);
	//ierr = DMDASetUniformCoordinates(daf, 0.0,1.0, 0.0,1.0, 0.0,1.0);CHKERRQ(ierr); /* not compatable with using -da_processors_x etc */
	
	ierr = x_DMDACreate3d(PETSC_COMM_SELF,wrap,DMDA_STENCIL_BOX,3,3,3, 1,1,1, 3,2, 0,0,0,&dac);CHKERRQ(ierr);
	ierr = DMDASetElementType_Q2(dac);CHKERRQ(ierr);
	//ierr = DMDASetUniformCoordinates(dac, 0.0,1.0, 0.0,1.0, 0.0,1.0);CHKERRQ(ierr);
	
	ierr = DMCreateInterpolation(dac,daf,&P,PETSC_NULL);CHKERRQ(ierr);
	
	/* check nlocal_elements_i is divisible by ref_i: I don't want to build any off proc elements needed for the cell problem */
	ierr = DMDAGetCornersElementQ2(dav_fine,&sei,&sej,&sek,&lmi,&lmj,&lmk);CHKERRQ(ierr);
	if (sei%refi != 0) { SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Start element (i) must be divisible by %D",refi); }
	if (sej%refj != 0) { SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Start element (j) must be divisible by %D",refj); }
	if (sek%refk != 0) { SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Start element (k) must be divisible by %D",refk); }

	if (lmi%refi != 0) { SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Local element size (i) must be divisible by %D",refi); }
	if (lmj%refj != 0) { SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Local element size (j) must be divisible by %D",refj); }
	if (lmk%refk != 0) { SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Local element size (k) must be divisible by %D",refk); }
	
	ierr = DMCreateMatrix(daf,MATSEQAIJ,&Acell);CHKERRQ(ierr);
	ierr = MatPtAPSymbolic(Acell,P,1.0,&Ac_el);CHKERRQ(ierr);
	
	/* quadrature */
	nqp = volQ_fine->npoints;
	P3D_prepare_elementQ2(nqp,WEIGHT,XI,NI,GNI);
	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ_fine,&quadrature_points);CHKERRQ(ierr);
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dav_fine,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dav_fine,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	/* indices */
	ierr = DMDAGetGlobalIndices(dav_fine,&NUM_GINDICES,&GINDICES);CHKERRQ(ierr);
	ierr = BCListApplyDirichletMask(NUM_GINDICES,GINDICES,u_bclist_fine);CHKERRQ(ierr);

	ierr = DMDAGetGlobalIndices(daf,&NUM_GINDICES_cell,&GINDICES_cell);CHKERRQ(ierr);

	ierr = DMDAGetGlobalIndices(dav_coarse,&NUM_GINDICES_coarse,&GINDICES_coarse);CHKERRQ(ierr);
	
	/* loop and assemble */
	ierr = DMDAGetElements_pTatinQ2P1(dav_fine,&nel,&nen_v,&elnidx_v);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(daf,&nel_cell,&nen_v_cell,&elnidx_v_cell);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dav_coarse,&nel_coarse,&nen_v_coarse,&elnidx_v_coarse);CHKERRQ(ierr);
	
	/* loop over cells on coarse grid */
	t[0] = t[1] = t[2] = t[3] = t[4] = t[5] = 0.0;
	PetscGetTime(&t[0]);
	for (kc=0; kc<lmk/refk; kc++) {
		for (jc=0; jc<lmj/refj; jc++) {
			for (ic=0; ic<lmi/refi; ic++) {
				PetscInt cidx_coarse;
				
				cidx_coarse = ic + jc*lmi_coarse + kc*lmi_coarse*lmj_coarse;
				
				
				ierr = MatZeroEntries(Acell);CHKERRQ(ierr);
				ierr = MatZeroEntries(Ac_el);CHKERRQ(ierr);
				
				/* look over fine cells contained in coarse (i call this a cell problem) */
				PetscGetTime(&t0);
				for (kkf=refk*kc; kkf<refk*kc+refk; kkf++) {
					for (jjf=refj*jc; jjf<refj*jc+refj; jjf++) {
						for (iif=refi*ic; iif<refi*ic+refi; iif++) {
							PetscInt cidx;
							PetscInt cidx_cell;
					
							cidx      = iif + jjf*lmi + kkf*lmi*lmj;
							cidx_cell = (iif-refi*ic) + (jjf-refj*jc)*refi + (kkf-refk*kc)*refi*refj;
							
							/* get global indices */
							for (n=0; n<NPE; n++) {
								int NID;
								
								/* global indices of FE problem */
								NID = elnidx_v[NPE*cidx + n];
								ge_eqnums[3*n  ] = GINDICES[ 3*NID   ];
								ge_eqnums[3*n+1] = GINDICES[ 3*NID+1 ];
								ge_eqnums[3*n+2] = GINDICES[ 3*NID+2 ];
								
								/* local indices of FE problem relative to cell problem */
								NID = elnidx_v_cell[NPE*cidx_cell + n];
								ge_eqnums_cell[3*n  ] = GINDICES_cell[ 3*NID   ];
								ge_eqnums_cell[3*n+1] = GINDICES_cell[ 3*NID+1 ];
								ge_eqnums_cell[3*n+2] = GINDICES_cell[ 3*NID+2 ];
							}
					
							ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx_v[nen_v*cidx],LA_gcoords);CHKERRQ(ierr);
							
							ierr = VolumeQuadratureGetCellData_Stokes(volQ_fine,quadrature_points,cidx,&cell_quadrature_points);CHKERRQ(ierr);
							
							/* initialise element stiffness matrix */
							PetscMemzero( Ae, sizeof(PetscScalar)* Q2_NODES_PER_EL_3D * Q2_NODES_PER_EL_3D * U_DOFS * U_DOFS );
							
							P3D_evaluate_geometry_elementQ2(nqp,el_coords,GNI, detJ,dNudx,dNudy,dNudz);

							/* evaluate the viscosity */
							for (q=0; q<nqp; q++) {
								el_eta[q] = cell_quadrature_points[q].eta;
							}
							
							/* assemble */
							for (q=0; q<nqp; q++) {
								fac = WEIGHT[q] * detJ[q];
								
								for (n = 0; n<NPE; n++) {
									PetscScalar d_dx_i = dNudx[q ][n];
									PetscScalar d_dy_i = dNudy[q ][n];
									PetscScalar d_dz_i = dNudz[q ][n];
									
									B[0][3*n  ] = d_dx_i; B[0][3*n+1] = 0.0;     B[0][3*n+2] = 0.0;
									B[1][3*n  ] = 0.0;    B[1][3*n+1] = d_dy_i;  B[1][3*n+2] = 0.0;
									B[2][3*n  ] = 0.0;    B[2][3*n+1] = 0.0;     B[2][3*n+2] = d_dz_i;
									
									B[3][3*n] = d_dy_i;   B[3][3*n+1] = d_dx_i;  B[3][3*n+2] = 0.0;   /* e_xy */
									B[4][3*n] = d_dz_i;   B[4][3*n+1] = 0.0;     B[4][3*n+2] = d_dx_i;/* e_xz */
									B[5][3*n] = 0.0;      B[5][3*n+1] = d_dz_i;  B[5][3*n+2] = d_dy_i;/* e_yz */
								}
								diagD[0] = 2.0*fac*el_eta[q ];
								diagD[1] = 2.0*fac*el_eta[q ];
								diagD[2] = 2.0*fac*el_eta[q ];
								
								diagD[3] =     fac*el_eta[q ];
								diagD[4] =     fac*el_eta[q ];
								diagD[5] =     fac*el_eta[q ];
								
								/* form Bt tildeD B */
								for (ii = 0; ii<81; ii++) {
									for (jj = ii; jj<81; jj++) {
										for (kk = 0; kk<6; kk++) {
											Ae[ii*81+jj] += B[kk][ii]*diagD[kk]*B[kk][jj];
										}
									}
								}
							}
							/* fill lower triangular part */
							for (ii = 0; ii < 81; ii++) {
								for (jj = ii; jj < 81; jj++) {
									Ae[jj*81+ii] = Ae[ii*81+jj];
								}
							}
							
							/* mask out any row/cols associated with boundary conditions */
							for (n=0; n<3*NPE; n++) {
								if (ge_eqnums[n] < 0) {
									
									ii = n;
									for (jj=0; jj<81; jj++) {
										Ae[ii*81+jj] = 0.0;
									}
									
									jj = n;
									for (ii=0; ii<81; ii++) {
										Ae[ii*81+jj] = 0.0;
									}
								}
							}
							ierr = MatSetValues(Acell,Q2_NODES_PER_EL_3D * U_DOFS,ge_eqnums_cell,Q2_NODES_PER_EL_3D * U_DOFS,ge_eqnums_cell,Ae,ADD_VALUES);CHKERRQ(ierr);
							
						}
					}
				}
				ierr = MatAssemblyBegin(Acell,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
				ierr = MatAssemblyEnd (Acell,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
				
				/* For each cell problem, define the boundary conditions (kinda ugly hack as this doesnt use bc_list_fine */				
				PetscGetTime(&t1);
				t[2] += (t1-t0);
				
				t0 = t1;
				for (kkf=refk*kc; kkf<refk*kc+refk; kkf++) {
					for (jjf=refj*jc; jjf<refj*jc+refj; jjf++) {
						for (iif=refi*ic; iif<refi*ic+refi; iif++) {
							PetscInt cidx;
							PetscInt cidx_cell;
							
							cidx      = iif + jjf*lmi + kkf*lmi*lmj;
							cidx_cell = (iif-refi*ic) + (jjf-refj*jc)*refi + (kkf-refk*kc)*refi*refj;
							
							for (n=0; n<NPE; n++) {
								int NID;
								
								NID = elnidx_v[NPE*cidx + n];
								ge_eqnums[3*n  ] = GINDICES[ 3*NID   ];
								ge_eqnums[3*n+1] = GINDICES[ 3*NID+1 ];
								ge_eqnums[3*n+2] = GINDICES[ 3*NID+2 ];
								
								NID = elnidx_v_cell[NPE*cidx_cell + n];
								ge_eqnums_cell[3*n  ] = GINDICES_cell[ 3*NID   ];
								ge_eqnums_cell[3*n+1] = GINDICES_cell[ 3*NID+1 ];
								ge_eqnums_cell[3*n+2] = GINDICES_cell[ 3*NID+2 ];
							}

							for (n=0; n<3*NPE; n++) {
								if (ge_eqnums[n] < 0) {
									ii = ge_eqnums_cell[n];
									ierr = MatSetValue(Acell,ii,ii,1.0,INSERT_VALUES);CHKERRQ(ierr);
								}
							}
						}
					}
				}
				
				ierr = MatAssemblyBegin(Acell,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
				ierr = MatAssemblyEnd (Acell,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
				PetscGetTime(&t1);
				t[3] += (t1-t0);

				t0 = t1;
				ierr = MatPtAPNumeric(Acell,P,Ac_el);CHKERRQ(ierr);
				PetscGetTime(&t1);
				t[4] += (t1-t0);

				/* assemble coarse grid operator */
				
				/* get global indices */
				t0 = t1;
				for (n=0; n<NPE; n++) {
					int NID;
					
					NID = elnidx_v_coarse[NPE*cidx_coarse + n];
					ge_eqnums_coarse[3*n  ] = GINDICES_coarse[ 3*NID   ];
					ge_eqnums_coarse[3*n+1] = GINDICES_coarse[ 3*NID+1 ];
					ge_eqnums_coarse[3*n+2] = GINDICES_coarse[ 3*NID+2 ];
				}
				ierr = MatGetArray(Ac_el,&Ac_entries);CHKERRQ(ierr);

				ierr = MatSetValues(Ac,Q2_NODES_PER_EL_3D*U_DOFS,ge_eqnums_coarse,Q2_NODES_PER_EL_3D*U_DOFS,ge_eqnums_coarse,Ac_entries,ADD_VALUES);CHKERRQ(ierr);
				
				ierr = MatRestoreArray(Ac_el,&Ac_entries);CHKERRQ(ierr);
				PetscGetTime(&t1);
				t[5] += (t1-t0);

			}
		}
	}
	ierr = MatAssemblyBegin(Ac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd (Ac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	PetscGetTime(&t[1]);

	PetscPrintf(PETSC_COMM_WORLD,"  Time[Assemble cell prb]   %1.4e (sec)\n",t[2]);
	PetscPrintf(PETSC_COMM_WORLD,"  Time[BC insertion]        %1.4e (sec)\n",t[3]);
	PetscPrintf(PETSC_COMM_WORLD,"  Time[Form PtAP]           %1.4e (sec)\n",t[4]);
	PetscPrintf(PETSC_COMM_WORLD,"  Time[Assemble coarse prb] %1.4e (sec)\n",t[5]);
	PetscPrintf(PETSC_COMM_WORLD,"  Time[Total]               %1.4e (sec)\n",t[1]-t[0]);
	
	ierr = BCListRemoveDirichletMask(NUM_GINDICES,GINDICES,u_bclist_fine);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = MatDestroy(&Acell);CHKERRQ(ierr);
	ierr = MatDestroy(&Ac_el);CHKERRQ(ierr);
	ierr = DMDestroy(&daf);CHKERRQ(ierr);
	ierr = DMDestroy(&dac);CHKERRQ(ierr);
	
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "FormJacobian_Stokes"
PetscErrorCode FormJacobian_Stokes(SNES snes,Vec X,Mat *A,Mat *B,MatStructure *mstr,void *ctx)
{
  pTatinCtx         user;
  DM                stokes_pack,dau,dap;
	IS                *is;
	PhysCompStokes    stokes;
  Vec               Uloc,Ploc;
  PetscScalar       *LA_Uloc,*LA_Ploc;
	PetscBool         is_mffd = PETSC_FALSE;
	PetscBool         is_nest = PETSC_FALSE;
	PetscBool         is_shell = PETSC_FALSE;
  PetscErrorCode    ierr;
	
  PetscFunctionBegin;
	
	user = (pTatinCtx)ctx;

	ierr = pTatinGetStokesContext(user,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;

  ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(stokes_pack,X,Uloc,Ploc);CHKERRQ(ierr);
	ierr = VecGetArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecGetArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	ierr = DMCompositeGetGlobalISs(stokes_pack,&is);CHKERRQ(ierr);
	
	/* Jacobian */
	ierr = pTatin_EvaluateRheologyNonlinearities(user,dau,LA_Uloc,dap,LA_Ploc);CHKERRQ(ierr);
	
	ierr = PetscObjectTypeCompare((PetscObject)(*A),MATMFFD, &is_mffd);CHKERRQ(ierr);
	ierr = PetscObjectTypeCompare((PetscObject)(*A),MATNEST, &is_nest);CHKERRQ(ierr);
	ierr = PetscObjectTypeCompare((PetscObject)(*A),MATSHELL,&is_shell);CHKERRQ(ierr);

	if (is_nest) {
		Mat Auu;
		
		ierr = MatGetSubMatrix(*A,is[0],is[0],MAT_INITIAL_MATRIX,&Auu);CHKERRQ(ierr);

		is_shell = PETSC_FALSE;
		ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
			ierr = MatAssemble_StokesA_AUU(Auu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		ierr = MatDestroy(&Auu);CHKERRQ(ierr);
	}
	/* If shell, do nothing */
	/* If mffd,  do nothing */
	
	ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	/* preconditioner for Jacobian */
	{
		Mat Buu,Bpp;
		
		ierr = MatGetSubMatrix(*B,is[0],is[0],MAT_INITIAL_MATRIX,&Buu);CHKERRQ(ierr);
		ierr = MatGetSubMatrix(*B,is[1],is[1],MAT_INITIAL_MATRIX,&Bpp);CHKERRQ(ierr);
		
		is_shell = PETSC_FALSE;
		ierr = PetscObjectTypeCompare((PetscObject)Buu,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Buu);CHKERRQ(ierr);
			ierr = MatAssemble_StokesA_AUU(Buu,dau,user->stokes_ctx->u_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		is_shell = PETSC_FALSE;
		ierr = PetscObjectTypeCompare((PetscObject)Bpp,MATSHELL,&is_shell);CHKERRQ(ierr);
		if (!is_shell) {
			ierr = MatZeroEntries(Bpp);CHKERRQ(ierr);
			ierr = MatAssemble_StokesPC_ScaledMassMatrix(Bpp,dau,dap,user->stokes_ctx->p_bclist,user->stokes_ctx->volQ);CHKERRQ(ierr);
		}
		
		ierr = MatDestroy(&Buu);CHKERRQ(ierr);
		ierr = MatDestroy(&Bpp);CHKERRQ(ierr);		
  }
  ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	*mstr = DIFFERENT_NONZERO_PATTERN;
	
	/* clean up */
	ierr = ISDestroy(&is[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is[1]);CHKERRQ(ierr);
	ierr = PetscFree(is);CHKERRQ(ierr);
	
	ierr = VecRestoreArray(Uloc,&LA_Uloc);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ploc,&LA_Ploc);CHKERRQ(ierr);
	
  ierr = DMCompositeRestoreLocalVectors(stokes_pack,&Uloc,&Ploc);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dStokesBuildMeshHierarchy"
PetscErrorCode pTatin3dStokesBuildMeshHierarchy(DM dav,PetscInt nlevels,DM dav_hierarchy[])
{
	PetscErrorCode ierr;
	DM *coarsened_list;
	PetscInt k;
	
	PetscFunctionBegin;
	
	/* set up mg */
	dav->ops->coarsenhierarchy = DMCoarsenHierarchy2_DA;
	
	dav_hierarchy[ nlevels-1 ] = dav;
	ierr = PetscObjectReference((PetscObject)dav);CHKERRQ(ierr);
	
	/* Coarsen nlevels - 1 times, and add levels into list so that level 0 is the coarsest */
	ierr = PetscMalloc(sizeof(DM)*(nlevels-1),&coarsened_list);CHKERRQ(ierr);
	ierr = DMCoarsenHierarchy(dav,nlevels-1,coarsened_list);CHKERRQ(ierr);
	for (k=0; k<nlevels-1; k++) {
		dav_hierarchy[ nlevels-2-k ] = coarsened_list[k];
	}
	PetscFree(coarsened_list);
	
	/* Set all dav's to be of type Q2 */
	for (k=0; k<nlevels-1; k++) {
		ierr = PetscObjectSetOptionsPrefix((PetscObject)dav_hierarchy[k],"stk_velocity_");CHKERRQ(ierr);
		ierr = DMDASetElementType_Q2(dav_hierarchy[k]);CHKERRQ(ierr);
	}
	
	/* inject coordinates */
	ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);
	
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dStokesReportMeshHierarchy"
PetscErrorCode pTatin3dStokesReportMeshHierarchy(PetscInt nlevels,DM dav_hierarchy[])
{
	PetscErrorCode ierr;
	PetscInt       k,lmx,lmy,lmz,si,sj,sk;
	PetscInt       nels,nen;
	const PetscInt *els;
	PetscMPIInt    rank;
	
	PetscFunctionBegin;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	/* Report mesh sizes */
	for (k=0; k<nlevels; k++) {
		ierr = DMDAGetSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"         level [%2D]: global Q2 elements (%D x %D x %D) \n", k,lmx,lmy,lmz );
	}		

	
	for (k=0; k<nlevels; k++) {
		
		ierr = DMDAGetElements_pTatinQ2P1(dav_hierarchy[k],&nels,&nen,&els);CHKERRQ(ierr);
		ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		if (rank<1) {
			PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: local Q2 elements  (%D x %D x %D) \n", rank, k,lmx,lmy,lmz );
		}
	}

	for (k=0; k<nlevels; k++) {
		ierr = DMDAGetElements_pTatinQ2P1(dav_hierarchy[k],&nels,&nen,&els);CHKERRQ(ierr);
		
		ierr = DMDAGetCornersElementQ2(dav_hierarchy[k],&si,&sj,&sk,&lmx,&lmy,&lmz);CHKERRQ(ierr);
		si = si/2;
		sj = sj/2;
		sk = sk/2;
		if (rank<1) {
			PetscPrintf(PETSC_COMM_SELF,"[r%4D]: level [%2D]: element range [%D - %D] x [%D - %D] x [%D - %D] \n", rank, k,si,si+lmx-1,sj,sj+lmy-1,sk,sk+lmz-1 );
		}
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3dCreateStokesOperators"
PetscErrorCode pTatin3dCreateStokesOperators(PhysCompStokes stokes_ctx,IS is_stokes_field[],
																						 PetscInt nlevels,DM dav_hierarchy[],Mat interpolation_v[],
																						 BCList u_bclist[],Quadrature volQ[],
																						 Mat *_A,Mat operatorA11[],Mat *_B,Mat operatorB11[])
{
	Mat            A,B;
	OperatorType   level_type[10];
	DM             dav,dap;
	PetscInt       k,max;
	PetscBool      flg;
	static int     been_here = 0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	dav = stokes_ctx->dav;
	dap = stokes_ctx->dap;
	
	/* A operator */
	ierr = StokesQ2P1CreateMatrix_Operator(stokes_ctx,&A);CHKERRQ(ierr);
	/* memory saving - only need daU IF you want to split A11 into A11uu,A11vv,A11ww */
	{
		MatStokesMF mf;
		
		ierr = MatShellGetMatStokesMF(A,&mf);CHKERRQ(ierr);
		ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
		mf->daU = PETSC_NULL;
	}
	
	/* B operator */
	{
		Mat         Aup,Apu,Spp,bA[2][2];
		MatStokesMF StkCtx;
		
		ierr = MatShellGetMatStokesMF(A,&StkCtx);CHKERRQ(ierr);
		
		/* Schur complement */
		ierr = DMCreateMatrix(dap,MATSBAIJ,&Spp);CHKERRQ(ierr);
		ierr = MatSetOptionsPrefix(Spp,"S*_");CHKERRQ(ierr);
		ierr = MatSetOption(Spp,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
		ierr = MatSetFromOptions(Spp);CHKERRQ(ierr);
		
		/* A12 */
		ierr = StokesQ2P1CreateMatrix_MFOperator_A12(StkCtx,&Aup);CHKERRQ(ierr);
		ierr = MatSetOptionsPrefix(Aup,"Bup_");CHKERRQ(ierr);
		ierr = MatSetFromOptions(Aup);CHKERRQ(ierr);
		
		/* A21 */
		ierr = StokesQ2P1CreateMatrix_MFOperator_A21(StkCtx,&Apu);CHKERRQ(ierr);
		ierr = MatSetOptionsPrefix(Apu,"Bpu_");CHKERRQ(ierr);
		ierr = MatSetFromOptions(Apu);CHKERRQ(ierr);
		
		/* nest */
		bA[0][0] = PETSC_NULL; bA[0][1] = Aup;
		bA[1][0] = Apu;        bA[1][1] = Spp;
		
		ierr = MatCreateNest(PETSC_COMM_WORLD,2,is_stokes_field,2,is_stokes_field,&bA[0][0],&B);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		
		/* tidy up - hand back destruction to B */
		ierr = MatDestroy(&Aup);CHKERRQ(ierr);
		ierr = MatDestroy(&Apu);CHKERRQ(ierr);
		ierr = MatDestroy(&Spp);CHKERRQ(ierr);
	}
	
	/* A11 operator */	
	/* defaults */
	level_type[0] = OP_TYPE_REDISC_ASM;
	for (k=1; k<nlevels; k++) {
		level_type[k] = OP_TYPE_REDISC_MF;
	}
	
	max = nlevels;
	ierr = PetscOptionsGetIntArray(PETSC_NULL,"-A11_operator_type",(PetscInt*)level_type,&max,&flg);CHKERRQ(ierr);
	for (k=nlevels-1; k>=0; k--) {
		
		switch (level_type[k]) {
				
			case OP_TYPE_REDISC_ASM:
			{
				Mat Auu;
				PetscBool same1 = PETSC_FALSE,same2 = PETSC_FALSE,same3 = PETSC_FALSE;
				
				/* use -stk_velocity_da_mat_type sbaij or -Buu_da_mat_type sbaij */
				if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Re-discretisation :: assembled operator \n", k);
				ierr = DMCreateMatrix(dav_hierarchy[k],MATSBAIJ,&Auu);CHKERRQ(ierr);
				ierr = MatSetOptionsPrefix(Auu,"Buu_");CHKERRQ(ierr);
				ierr = MatSetFromOptions(Auu);CHKERRQ(ierr);
				ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSBAIJ,&same1);CHKERRQ(ierr);
				ierr = PetscObjectTypeCompare((PetscObject)Auu,MATSEQSBAIJ,&same2);CHKERRQ(ierr);
				ierr = PetscObjectTypeCompare((PetscObject)Auu,MATMPISBAIJ,&same3);CHKERRQ(ierr);
				if (same1||same2||same3) {
					ierr = MatSetOption(Auu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE);CHKERRQ(ierr);
				}
				/* should move assembly into jacobian */
				ierr = MatZeroEntries(Auu);CHKERRQ(ierr);
				ierr = MatAssemble_StokesA_AUU(Auu,dav_hierarchy[k],u_bclist[k],volQ[k]);CHKERRQ(ierr);
				
				operatorA11[k] = Auu;
				operatorB11[k] = Auu;
				ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
				
			}
				break;
				
			case OP_TYPE_REDISC_MF:
			{
				Mat Auu;
				MatA11MF mf,A11Ctx;
				
				if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Re-discretisation :: matrix free operator \n", k);
				ierr = MatA11MFCreate(&A11Ctx);CHKERRQ(ierr);
				ierr = MatA11MFSetup(A11Ctx,dav_hierarchy[k],volQ[k],u_bclist[k]);CHKERRQ(ierr);
				
				ierr = StokesQ2P1CreateMatrix_MFOperator_A11(A11Ctx,&Auu);CHKERRQ(ierr);
				ierr = MatShellGetMatA11MF(Auu,&mf);CHKERRQ(ierr);
				ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
				mf->daU = PETSC_NULL;				
				operatorA11[k] = Auu;
				
				{
					PetscBool use_low_order_geometry = PETSC_FALSE;
					
					ierr = PetscOptionsGetBool(PETSC_NULL,"-use_low_order_geometry",&use_low_order_geometry,PETSC_NULL);CHKERRQ(ierr);
					if (use_low_order_geometry==PETSC_TRUE) {
						Mat Buu;
						
						if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"  Activiting low order A11 operator \n");
						ierr = StokesQ2P1CreateMatrix_MFOperator_A11LowOrder(A11Ctx,&Buu);CHKERRQ(ierr);
						ierr = MatShellGetMatA11MF(Buu,&mf);CHKERRQ(ierr);
						ierr = DMDestroy(&mf->daU);CHKERRQ(ierr);
						mf->daU = PETSC_NULL;				
						operatorB11[k] = Buu;
						
					} else {
						operatorB11[k] = Auu;
						ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
					}
				}
				
				
				ierr = MatA11MFDestroy(&A11Ctx);CHKERRQ(ierr);
			}
				break;
				
			case OP_TYPE_GALERKIN:
			{
				Mat Auu;
				
				if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: Galerkin :: assembled operator \n", k);
				if (k==nlevels-1) {
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid on the finest level");
				}	
				if (level_type[k+1] == OP_TYPE_REDISC_MF) {
					SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use galerkin coarse grid. Finest level above must be of type OP_TYPE_REDISC_ASM, OP_TYPE_GALERKIN or OP_TYPE_MFGALERKIN",k+1);
				}
				
				/* should move coarse grid assembly into jacobian */
				ierr = MatPtAP(operatorA11[k+1],interpolation_v[k+1],MAT_INITIAL_MATRIX,1.0,&Auu);CHKERRQ(ierr);
				
				operatorA11[k] = Auu;
				operatorB11[k] = Auu;
				ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
			}
				break;
				
			case OP_TYPE_MFGALERKIN:
			{
				Mat Auu;
				
				if (!been_here) PetscPrintf(PETSC_COMM_WORLD,"Level [%d]: Coarse grid type :: MFGalerkin :: assembled operator \n", k);
				if (k==nlevels-1) {
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use mf-galerkin coarse grid on the finest level");
				}	
				if (level_type[k+1] != OP_TYPE_REDISC_MF) {
					SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Cannot use mf-galerkin. Next finest level[%D] must be of type OP_TYPE_REDISC_MF",k+1);
				}	
				
				
				ierr = DMCreateMatrix(dav_hierarchy[k],MATAIJ,&Auu);CHKERRQ(ierr);
				ierr = MatAssembleMFGalerkin(dav_hierarchy[k+1],u_bclist[k+1],volQ[k+1],dav_hierarchy[k],Auu);CHKERRQ(ierr);
				
				operatorA11[k] = Auu;
				operatorB11[k] = Auu;
				ierr = PetscObjectReference((PetscObject)Auu);CHKERRQ(ierr);
				
			}
				break;
				
			default:
				break;
		}
	}	
	/*
	{
		Mat Ac;
		
		ierr = DMCreateMatrix(dav_hierarchy[nlevels-2],MATAIJ,&Ac);CHKERRQ(ierr);
		ierr = MatAssembleMFGalerkin(dav_hierarchy[nlevels-1],u_bclist[nlevels-1],volQ[nlevels-1],dav_hierarchy[nlevels-2],Ac);CHKERRQ(ierr);
	}
	*/
	
	/* Set fine A11 into nest */
	ierr = MatNestSetSubMat(B,0,0,operatorA11[nlevels-1]);CHKERRQ(ierr);

	*_A = A;
	*_B = B;

	been_here = 1;
	PetscFunctionReturn(0);
}
	
#undef __FUNCT__  
#define __FUNCT__ "pTatin3dStokesKSPConfigureFSGMG"
PetscErrorCode pTatin3dStokesKSPConfigureFSGMG(KSP ksp,PetscInt nlevels,Mat operatorA11[],Mat operatorB11[],Mat interpolation_v[])
{
	PetscInt k,nsplits;
	PC       pc,pc_i;
	KSP      *sub_ksp,ksp_coarse,ksp_smoother;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = KSPSetUp(ksp);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
	
	ierr = KSPGetPC(sub_ksp[0],&pc_i);CHKERRQ(ierr);
	ierr = PCSetType(pc_i,PCMG);CHKERRQ(ierr);
	ierr = PCMGSetLevels(pc_i,nlevels,PETSC_NULL);CHKERRQ(ierr);
	ierr = PCMGSetType(pc_i,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
	ierr = PCMGSetGalerkin(pc_i,PETSC_FALSE);CHKERRQ(ierr);
		ierr = PCSetDM(pc_i,PETSC_NULL);CHKERRQ(ierr);
	
	for( k=1; k<nlevels; k++ ){
		ierr = PCMGSetInterpolation(pc_i,k,interpolation_v[k]);CHKERRQ(ierr);
	}
	
	/* drop the operators in - i presume this will also need to be performed inside the jacobian each time the operators are modified */
	/* No - it looks like PCSetUp_MG will call set operators on all levels if the SetOperators was called on the finest, which should/is done by the SNES */
	ierr = PCMGGetCoarseSolve(pc_i,&ksp_coarse);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp_coarse,operatorA11[0],operatorA11[0],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	for( k=1; k<nlevels; k++ ){
		PetscBool use_low_order_geometry = PETSC_FALSE;
		
		ierr = PCMGGetSmoother(pc_i,k,&ksp_smoother);CHKERRQ(ierr);
		
		// use A for smoother, B for residual
		ierr = PetscOptionsGetBool(PETSC_NULL,"-use_low_order_geometry",&use_low_order_geometry,PETSC_NULL);CHKERRQ(ierr);
		if (use_low_order_geometry==PETSC_TRUE) {
			ierr = KSPSetOperators(ksp_smoother,operatorB11[k],operatorB11[k],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		//ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorB11[k],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		} else {
			// Use A for smoother, lo
			ierr = KSPSetOperators(ksp_smoother,operatorA11[k],operatorA11[k],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		}
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_linear_viscous_forward_model_driver"
PetscErrorCode pTatin3d_linear_viscous_forward_model_driver(int argc,char **argv)
{
	DM        multipys_pack,dav,dap;
	pTatinCtx user;
	Mat       A,B;
	Vec       X,F;
	IS        *is_stokes_field;
	SNES      snes;
	KSP       ksp;
	PC        pc;
	PetscInt       nlevels,k;
	Mat            operatorA11[10],operatorB11[10];
	DM             dav_hierarchy[10];
	Mat            interpolation_v[10],interpolation_eta[10];
	Quadrature     volQ[10];
	BCList         u_bclist[10];
	PetscInt       step;
  pTatinModel    model;
	PetscLogDouble time[2];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
	ierr = pTatinLogNote(user,"[ptatin_driver_linear_ts] -> new simulation");CHKERRQ(ierr);
	
	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);
	ierr = pTatinLogNote2(user,"[ptatin_driver_linear_ts] -> model loaded:",model->model_name);CHKERRQ(ierr);
	/* Check if model is being restarted from a checkpointed file */
	ierr = pTatin3dRestart(user);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
	
	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	user->pack = user->stokes_ctx->stokes_pack;
	
	/* fetch some local variables */
	multipys_pack = user->pack;
	dav           = user->stokes_ctx->dav;
	dap           = user->stokes_ctx->dap;
	
	/* IF I DON'T DO THIS, THE IS's OBTAINED FROM DMCompositeGetGlobalISs() are wrong !! */
	{
		ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = DMRestoreGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	}
 ierr = DMCompositeGetGlobalISs(multipys_pack,&is_stokes_field);CHKERRQ(ierr);	
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);

	ierr = pTatinLogBasicDMDA(user,"Velocity",dav);CHKERRQ(ierr);
	ierr = pTatinLogBasicDMDA(user,"Pressure",dap);CHKERRQ(ierr);
		
	/* interpolate point coordinates (needed if mesh was modified) */
	//for (e=0; e<QUAD_EDGES; e++) {
	//	ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
	//}

	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);


	/* setup mg */
	nlevels = 1;
	PetscOptionsGetInt(PETSC_NULL,"-dau_nlevels",&nlevels,0);
	PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%d x %d x %d) : MG levels %d  \n", user->mx,user->my,user->mz,nlevels );
	ierr = pTatin3dStokesBuildMeshHierarchy(dav,nlevels,dav_hierarchy);CHKERRQ(ierr);
	ierr = pTatin3dStokesReportMeshHierarchy(nlevels,dav_hierarchy);CHKERRQ(ierr);
	ierr = pTatinLogNote(user,"  [Velocity multi-grid hierarchy]");CHKERRQ(ierr);
	for (k=nlevels-1; k>=0; k--) {
		char name[128];
		sprintf(name,"vel_dmda_Lv%d",k);
		ierr = pTatinLogBasicDMDA(user,name,dav_hierarchy[k]);CHKERRQ(ierr);

		//sprintf(name,"vel_dmda_Lv%d.vtk",k);
		//ierr = DMDAViewPetscVTK(dav_hierarchy[k],PETSC_NULL,name);CHKERRQ(ierr);
	}
	
	/* Define interpolation operators for velocity space */
	interpolation_v[0] = PETSC_NULL;
	for (k=0; k<nlevels-1; k++) {
		ierr = DMCreateInterpolation(dav_hierarchy[k],dav_hierarchy[k+1],&interpolation_v[k+1],PETSC_NULL);CHKERRQ(ierr);
	}

	/* Define interpolation operators for scalar space */
	interpolation_eta[0] = PETSC_NULL;
	for (k=1; k<nlevels; k++) {
		ierr = MatMAIJRedimension(interpolation_v[k],1,&interpolation_eta[k]);CHKERRQ(ierr);
	}
	
	/* define material properties on gauss points on coarse grids */
	for (k=0; k<nlevels-1; k++) {
		PetscInt ncells,lmx,lmy,lmz;
		PetscInt np_per_dim;
		
		np_per_dim = 3;
		ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		ncells = lmx * lmy * lmz;
		ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&volQ[k]);CHKERRQ(ierr);
	}
	volQ[nlevels-1] = user->stokes_ctx->volQ;
	
	
	/* define bounary list on coarse grids */
	for (k=0; k<nlevels-1; k++) {
		ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);
	}
	u_bclist[nlevels-1] = user->stokes_ctx->u_bclist;

	/* update markers = >> gauss points */
	{
		int               npoints;
		DataField         PField_std;
		DataField         PField_stokes;
		MPntStd           *mp_std;
		MPntPStokes       *mp_stokes;
		
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
		DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
		
		DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
		mp_std    = PField_std->data; /* should write a function to do this */
		mp_stokes = PField_stokes->data; /* should write a function to do this */
		
		ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(user->coefficient_projection_type,npoints,mp_std,mp_stokes,nlevels,interpolation_eta,dav_hierarchy,volQ);CHKERRQ(ierr);
	}
	
	
	/* define bc's for hiearchy */
	ierr = pTatinModel_ApplyBoundaryConditionMG(nlevels,u_bclist,dav_hierarchy,user->model,user);CHKERRQ(ierr);
	

	/* configure stokes opertors */
	ierr = pTatin3dCreateStokesOperators(user->stokes_ctx,is_stokes_field,
																			 nlevels,dav_hierarchy,interpolation_v,u_bclist,volQ,
																			 &A,operatorA11,&B,operatorB11);CHKERRQ(ierr);
	
	/* work vector for solution and residual */
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

	/* initial condition */
	ierr = pTatinModel_ApplyInitialSolution(user->model,user,X);CHKERRQ(ierr);
		
	/* boundary condition */
	{
		Vec velocity,pressure;
		
		ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		ierr = BCListInsert(user->stokes_ctx->u_bclist,velocity);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	}

	/* write the initial fields */
	ierr = pTatinModel_Output(user->model,user,X,"icbc");CHKERRQ(ierr);

	
	ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
	ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
	ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);
	ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
	
	/* configure for fieldsplit */
	ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);

	ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_StdoutStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
//	ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_ParaviewStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
	
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCFIELDSPLIT);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc,"u",is_stokes_field[0]);CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(pc,"p",is_stokes_field[1]);CHKERRQ(ierr);

	/* configure uu split for galerkin multi-grid */
	ierr = pTatin3dStokesKSPConfigureFSGMG(ksp,nlevels,operatorA11,operatorB11,interpolation_v);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD,"   [[ COMPUTING FLOW FIELD FOR STEP : %D ]]\n", 0 );
	ierr = pTatinLogBasic(user);CHKERRQ(ierr);

	PetscGetTime(&time[0]);
	ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
	PetscGetTime(&time[1]);
	ierr = pTatinLogBasicSNES(user,"Stokes",snes);CHKERRQ(ierr);
	ierr = pTatinLogBasicCPUtime(user,"Stokes",time[1]-time[0]);CHKERRQ(ierr);
	ierr = pTatinLogPetscLog(user,"Stokes");CHKERRQ(ierr);

	ierr = pTatinViewBasicStokesSolutionResiduals(user,snes,multipys_pack,X);CHKERRQ(ierr);
	ierr = pTatinViewBasicStokesSolution(user,multipys_pack,X);CHKERRQ(ierr);
	
	/* dump */
	ierr = pTatinModel_Output(user->model,user,X,"step000000");CHKERRQ(ierr);

	/* compute timestep */
	user->dt = 1.0e30;
	{
		Vec velocity,pressure;
		PetscReal timestep;
		
		ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);		
		ierr = SwarmUpdatePosition_ComputeCourantStep(dav_hierarchy[nlevels-1],velocity,&timestep);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		
		ierr = pTatin_SetTimestep(user,"StkCourant",timestep);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  timestep[] dt_courant = %1.4e \n", user->dt );
	}
	/* checkpoint step 0 */
	ierr = pTatin3dCheckpoint(user,X,"step000000");CHKERRQ(ierr);
	
	/*
	{
		Vec velocity,pressure;
		PetscReal un,pn,Xn;
		VecNorm(X,NORM_2,&Xn); printf("Xn %1.4e\n",Xn);
		ierr = DMCompositeGetAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
		VecNorm(velocity,NORM_2,&un); printf("un %1.4e\n",un);
		VecNorm(pressure,NORM_2,&pn); printf("pn %1.4e\n",pn);
		ierr = DMCompositeRestoreAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
	}
	*/
	
	/* tidy up */
	for (k=0; k<nlevels; k++) {
		ierr = MatDestroy(&operatorA11[k]);CHKERRQ(ierr);
		ierr = MatDestroy(&operatorB11[k]);CHKERRQ(ierr);
	}
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	
	
	for (step=1; step <= user->nsteps; step++) {
		char stepname[128];
		Vec velocity,pressure;
		PetscReal timestep;

		/* update context time information */
		user->step = step;
		
		PetscPrintf(PETSC_COMM_WORLD,"<<----------------------------------------------------------------------------------------------->>\n");
		PetscPrintf(PETSC_COMM_WORLD,"   [[ EXECUTING TIME STEP : %D ]]\n", user->step );
		PetscPrintf(PETSC_COMM_WORLD,"     dt    : %1.4e \n", user->dt );
		PetscPrintf(PETSC_COMM_WORLD,"     time  : %1.4e \n", user->time+user->dt );
		/* update context time information */
		user->time = user->time + user->dt;

		ierr = pTatinLogBasic(user);CHKERRQ(ierr);
		
		
		/* update markers */
		{
			int npoints;
			MPntStd *mp_std;
			DataField PField;
			
			DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
			DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname ,&PField);
			mp_std = PField->data;
			
			ierr = DMCompositeGetAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
			ierr = SwarmUpdatePosition_MPntStd_Euler(dav_hierarchy[nlevels-1],velocity,user->dt,npoints,mp_std);CHKERRQ(ierr);
			ierr = DMCompositeRestoreAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
		}
		
		/* update mesh */
		ierr = pTatinModel_UpdateMeshGeometry(user->model,user,X);CHKERRQ(ierr);
		
		/* update mesh coordinate hierarchy */
		ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);
		
		/* 3 Update local coordinates and communicate */
		ierr = MaterialPointStd_UpdateCoordinates(user->materialpoint_db,dav_hierarchy[nlevels-1],user->materialpoint_ex);CHKERRQ(ierr);
		
		/* 3a - Add material */
		ierr = pTatinModel_ApplyMaterialBoundaryCondition(model,user);CHKERRQ(ierr);
		
		/* add / remove points if cells are over populated or depleted of points */
		ierr = MaterialPointPopulationControl_v1(user);CHKERRQ(ierr);
		
		
		/* update markers = >> gauss points */
		{
			int               npoints;
			DataField         PField_std;
			DataField         PField_stokes;
			MPntStd           *mp_std;
			MPntPStokes       *mp_stokes;
			
			DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
			DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
			
			DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
			mp_std    = PField_std->data; /* should write a function to do this */
			mp_stokes = PField_stokes->data; /* should write a function to do this */
			
			ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(user->coefficient_projection_type,npoints,mp_std,mp_stokes,nlevels,interpolation_eta,dav_hierarchy,volQ);CHKERRQ(ierr);
		}

		/* Update boundary conditions */
		/* Fine level setup */
		ierr = pTatinModel_ApplyBoundaryCondition(model,user);CHKERRQ(ierr);
		/* Coarse grid setup: Configure boundary conditions */
		ierr = pTatinModel_ApplyBoundaryConditionMG(nlevels,u_bclist,dav_hierarchy,model,user);CHKERRQ(ierr);
		
		/* solve */
		/* a) configure stokes opertors */
		ierr = pTatin3dCreateStokesOperators(user->stokes_ctx,is_stokes_field,
																				 nlevels,dav_hierarchy,interpolation_v,u_bclist,volQ,
																				 &A,operatorA11,&B,operatorB11);CHKERRQ(ierr);
		/* b) create solver */
		ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
		ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
		ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);
		ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
		
		/* c) configure for fieldsplit */
		ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
		ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_StdoutStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
		//	ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_ParaviewStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
		
		ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
		ierr = PCSetType(pc,PCFIELDSPLIT);CHKERRQ(ierr);
		ierr = PCFieldSplitSetIS(pc,"u",is_stokes_field[0]);CHKERRQ(ierr);
		ierr = PCFieldSplitSetIS(pc,"p",is_stokes_field[1]);CHKERRQ(ierr);
		
		ierr = pTatin3dStokesKSPConfigureFSGMG(ksp,nlevels,operatorA11,operatorB11,interpolation_v);CHKERRQ(ierr);

		/* e) solve */
		PetscPrintf(PETSC_COMM_WORLD,"   [[ COMPUTING FLOW FIELD FOR STEP : %D ]]\n", step );
		PetscGetTime(&time[0]);
		ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
		PetscGetTime(&time[1]);
		ierr = pTatinLogBasicSNES(user,"Stokes",snes);CHKERRQ(ierr);
		ierr = pTatinLogBasicCPUtime(user,"Stokes",time[1]-time[0]);CHKERRQ(ierr);

		/*
		{
			Vec velocity,pressure;
			PetscReal un,pn,Xn;
			VecNorm(X,NORM_2,&Xn); printf("Xn %1.4e\n",Xn);
			ierr = DMCompositeGetAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
			VecNorm(velocity,NORM_2,&un); printf("un %1.4e\n",un);
			VecNorm(pressure,NORM_2,&pn); printf("pn %1.4e\n",pn);
			ierr = DMCompositeRestoreAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
		}
		*/
		
		/* output */
		if ( (step%user->output_frequency == 0) || (step == 1) ) {
			sprintf(stepname,"step%1.6d",step);
			ierr = pTatinModel_Output(user->model,user,X,stepname);CHKERRQ(ierr);
		}
		
		/* compute timestep */
		user->dt = 1.0e32;
		ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		ierr = SwarmUpdatePosition_ComputeCourantStep(dav_hierarchy[nlevels-1],velocity,&timestep);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		ierr = pTatin_SetTimestep(user,"StkCourant",timestep);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  timestep[%d] dt_courant = %1.4e \n", step,user->dt );
		
		/* CHECKPOINT */
		ierr = pTatin3dCheckpointManager(user,X);CHKERRQ(ierr);

		/* tidy up */
		for (k=0; k<nlevels; k++) {
			ierr = MatDestroy(&operatorA11[k]);CHKERRQ(ierr);
			ierr = MatDestroy(&operatorB11[k]);CHKERRQ(ierr);
		}
		ierr = MatDestroy(&A);CHKERRQ(ierr);
		ierr = MatDestroy(&B);CHKERRQ(ierr);
		ierr = SNESDestroy(&snes);CHKERRQ(ierr);

	
		/* update context time information */
		//user->time = user->time + user->dt;
		//user->step = step;
		/* Terminate time stepping */
		if (user->time >= user->time_max) {
			break;
		}
		
	}
	
	
	/* Clean up */
	for (k=0; k<nlevels-1; k++) {
		ierr = BCListDestroy(&u_bclist[k]);CHKERRQ(ierr);
	}
	for (k=0; k<nlevels-1; k++) {
		ierr = QuadratureDestroy(&volQ[k]);CHKERRQ(ierr);
	}
	for (k=1; k<nlevels; k++) {
		ierr = MatDestroy(&interpolation_v[k]);CHKERRQ(ierr);
		ierr = MatDestroy(&interpolation_eta[k]);CHKERRQ(ierr);
	}
	for (k=0; k<nlevels; k++) {
		ierr = DMDestroy(&dav_hierarchy[k]);CHKERRQ(ierr);
	}
	
	ierr = ISDestroy(&is_stokes_field[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is_stokes_field[1]);CHKERRQ(ierr);
	ierr = PetscFree(is_stokes_field);CHKERRQ(ierr);
	
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = pTatin3dDestroyContext(&user);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_linear_viscous_forward_model_driver_RESTART"
PetscErrorCode pTatin3d_linear_viscous_forward_model_driver_RESTART(int argc,char **argv)
{
	DM        multipys_pack,dav,dap;
	pTatinCtx user;
	Mat       A,B;
	Vec       X,F;
	IS        *is_stokes_field;
	SNES      snes;
	KSP       ksp;
	PC        pc;
	PetscInt       nlevels,k;
	Mat            operatorA11[10],operatorB11[10];
	DM             dav_hierarchy[10];
	Mat            interpolation_v[10],interpolation_eta[10];
	Quadrature     volQ[10];
	BCList         u_bclist[10];
  pTatinModel    model;
	
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);
	ierr = pTatinLogNote(user,"[ptatin_driver_linear_ts] -> restarted simulation");CHKERRQ(ierr);
	
	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	ierr = pTatinGetModel(user,&model);CHKERRQ(ierr);
	ierr = pTatinLogNote2(user,"[ptatin_driver_linear_ts] -> model loaded:",model->model_name);CHKERRQ(ierr);
	/* Check if model is being restarted from a checkpointed file */
	ierr = pTatin3dRestart(user);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);
	
	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	user->pack = user->stokes_ctx->stokes_pack;
	
	/* fetch some local variables */
	multipys_pack = user->pack;
	dav           = user->stokes_ctx->dav;
	dap           = user->stokes_ctx->dap;
	
	/* IF I DON'T DO THIS, THE IS's OBTAINED FROM DMCompositeGetGlobalISs() are wrong !! */
	{
		ierr = DMGetGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
		ierr = DMRestoreGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	}
	ierr = DMCompositeGetGlobalISs(multipys_pack,&is_stokes_field);CHKERRQ(ierr);	
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	ierr = pTatinLogBasicDMDA(user,"Velocity",dav);CHKERRQ(ierr);
	ierr = pTatinLogBasicDMDA(user,"Pressure",dap);CHKERRQ(ierr);
		
	/* interpolate point coordinates (needed if mesh was modified) */
	//for (e=0; e<QUAD_EDGES; e++) {
	//	ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
	//}
	
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);
	
	
	/* setup mg */
	nlevels = 1;
	PetscOptionsGetInt(PETSC_NULL,"-dau_nlevels",&nlevels,0);
	PetscPrintf(PETSC_COMM_WORLD,"Mesh size (%d x %d x %d) : MG levels %d  \n", user->mx,user->my,user->mz,nlevels );
	ierr = pTatin3dStokesBuildMeshHierarchy(dav,nlevels,dav_hierarchy);CHKERRQ(ierr);
	ierr = pTatin3dStokesReportMeshHierarchy(nlevels,dav_hierarchy);CHKERRQ(ierr);
	ierr = pTatinLogNote(user,"  [Velocity multi-grid hierarchy]");CHKERRQ(ierr);
	for (k=nlevels-1; k>=0; k--) {
		char name[128];
		sprintf(name,"vel_dmda_Lv%d",k);
		ierr = pTatinLogBasicDMDA(user,name,dav_hierarchy[k]);CHKERRQ(ierr);
	}
	
	/* Define interpolation operators for velocity space */
	interpolation_v[0] = PETSC_NULL;
	for (k=0; k<nlevels-1; k++) {
		ierr = DMCreateInterpolation(dav_hierarchy[k],dav_hierarchy[k+1],&interpolation_v[k+1],PETSC_NULL);CHKERRQ(ierr);
	}
	
	/* Define interpolation operators for scalar space */
	interpolation_eta[0] = PETSC_NULL;
	for (k=1; k<nlevels; k++) {
		ierr = MatMAIJRedimension(interpolation_v[k],1,&interpolation_eta[k]);CHKERRQ(ierr);
	}
	
	/* define material properties on gauss points on coarse grids */
	for (k=0; k<nlevels-1; k++) {
		PetscInt ncells,lmx,lmy,lmz;
		PetscInt np_per_dim;
		
		np_per_dim = 3;
		ierr = DMDAGetLocalSizeElementQ2(dav_hierarchy[k],&lmx,&lmy,&lmz);CHKERRQ(ierr);
		ncells = lmx * lmy * lmz;
		ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&volQ[k]);CHKERRQ(ierr);
	}
	volQ[nlevels-1] = user->stokes_ctx->volQ;
	
	
	/* define bounary list on coarse grids */
	for (k=0; k<nlevels-1; k++) {
		ierr = DMDABCListCreate(dav_hierarchy[k],&u_bclist[k]);CHKERRQ(ierr);
	}
	u_bclist[nlevels-1] = user->stokes_ctx->u_bclist;
	
	/* define bc's for hiearchy */
	ierr = pTatinModel_ApplyBoundaryConditionMG(nlevels,u_bclist,dav_hierarchy,user->model,user);CHKERRQ(ierr);
	
	/* configure stokes opertors */
	ierr = pTatin3dCreateStokesOperators(user->stokes_ctx,is_stokes_field,
																			 nlevels,dav_hierarchy,interpolation_v,u_bclist,volQ,
																			 &A,operatorA11,&B,operatorB11);CHKERRQ(ierr);
	
	/* work vector for solution and residual */
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
	
	/* initial condition  */
	ierr = pTatinModel_ApplyInitialSolution(user->model,user,X);CHKERRQ(ierr);
	
	/* boundary condition */
	{
		Vec velocity,pressure;
		
		ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		ierr = BCListInsert(user->stokes_ctx->u_bclist,velocity);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	}

	/*
	{
		Vec velocity,pressure;
		PetscReal un,pn,Xn;
		VecNorm(X,NORM_2,&Xn); printf("Xn %1.4e\n",Xn);
		ierr = DMCompositeGetAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
		VecNorm(velocity,NORM_2,&un); printf("un %1.4e\n",un);
		VecNorm(pressure,NORM_2,&pn); printf("pn %1.4e\n",pn);
		ierr = DMCompositeRestoreAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
	}
	*/
	
	/* loaded step N, solving for next step */
	user->step = user->step + 1;
	while (user->step <= user->nsteps) {
		char      stepname[128];
		Vec       velocity,pressure;
		PetscReal timestep;
		
		
		PetscPrintf(PETSC_COMM_WORLD,"<<----------------------------------------------------------------------------------------------->>\n");
		PetscPrintf(PETSC_COMM_WORLD,"   [[ EXECUTING TIME STEP : %D ]]\n", user->step );
		PetscPrintf(PETSC_COMM_WORLD,"     dt    : %1.4e \n", user->dt );
		PetscPrintf(PETSC_COMM_WORLD,"     time  : %1.4e \n", user->time + user->dt);
		/* update context time information */
		user->time = user->time + user->dt;

		ierr = pTatinLogBasic(user);CHKERRQ(ierr);
		
		/* update markers */
		{
			int npoints;
			MPntStd *mp_std;
			DataField PField;
			
			DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
			DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname ,&PField);
			mp_std = PField->data;
			
			ierr = DMCompositeGetAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
			ierr = SwarmUpdatePosition_MPntStd_Euler(dav_hierarchy[nlevels-1],velocity,user->dt,npoints,mp_std);CHKERRQ(ierr);
			ierr = DMCompositeRestoreAccess(user->pack,X,&velocity,&pressure);CHKERRQ(ierr);
		}
		
		/* update mesh */
		ierr = pTatinModel_UpdateMeshGeometry(user->model,user,X);CHKERRQ(ierr);
		
		/* update mesh coordinate hierarchy */
		ierr = DMDARestrictCoordinatesHierarchy(dav_hierarchy,nlevels);CHKERRQ(ierr);
		
		/* 3 Update local coordinates and communicate */
		ierr = MaterialPointStd_UpdateCoordinates(user->materialpoint_db,dav_hierarchy[nlevels-1],user->materialpoint_ex);CHKERRQ(ierr);
		
		/* 3a - Add material */
		ierr = pTatinModel_ApplyMaterialBoundaryCondition(model,user);CHKERRQ(ierr);

		/* add / remove points if cells are over populated or depleted of points */
		ierr = MaterialPointPopulationControl_v1(user);CHKERRQ(ierr);
		
		
		/* update markers = >> gauss points */
		{
			int               npoints;
			DataField         PField_std;
			DataField         PField_stokes;
			MPntStd           *mp_std;
			MPntPStokes       *mp_stokes;
			
			DataBucketGetDataFieldByName(user->materialpoint_db, MPntStd_classname     , &PField_std);
			DataBucketGetDataFieldByName(user->materialpoint_db, MPntPStokes_classname , &PField_stokes);
			
			DataBucketGetSizes(user->materialpoint_db,&npoints,PETSC_NULL,PETSC_NULL);
			mp_std    = PField_std->data; /* should write a function to do this */
			mp_stokes = PField_stokes->data; /* should write a function to do this */
			
			ierr = SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(user->coefficient_projection_type,npoints,mp_std,mp_stokes,nlevels,interpolation_eta,dav_hierarchy,volQ);CHKERRQ(ierr);
		}

		/* Update boundary conditions */
		/* Fine level setup */
		ierr = pTatinModel_ApplyBoundaryCondition(model,user);CHKERRQ(ierr);
		/* Coarse grid setup: Configure boundary conditions */
		ierr = pTatinModel_ApplyBoundaryConditionMG(nlevels,u_bclist,dav_hierarchy,model,user);CHKERRQ(ierr);
		
		/* solve */
		/* a) configure stokes opertors */
		ierr = pTatin3dCreateStokesOperators(user->stokes_ctx,is_stokes_field,
																				 nlevels,dav_hierarchy,interpolation_v,u_bclist,volQ,
																				 &A,operatorA11,&B,operatorB11);CHKERRQ(ierr);
		/* b) create solver */
		ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
		ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
		ierr = SNESSetJacobian(snes,A,B,FormJacobian_Stokes,user);CHKERRQ(ierr);
		ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
		
		/* c) configure for fieldsplit */
		ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
		ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_StdoutStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
		//	ierr = KSPMonitorSet(ksp,pTatin_KSPMonitor_ParaviewStokesResiduals3d,(void*)user,PETSC_NULL);CHKERRQ(ierr);
		
		ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
		ierr = PCSetType(pc,PCFIELDSPLIT);CHKERRQ(ierr);
		ierr = PCFieldSplitSetIS(pc,"u",is_stokes_field[0]);CHKERRQ(ierr);
		ierr = PCFieldSplitSetIS(pc,"p",is_stokes_field[1]);CHKERRQ(ierr);
		
		ierr = pTatin3dStokesKSPConfigureFSGMG(ksp,nlevels,operatorA11,operatorB11,interpolation_v);CHKERRQ(ierr);
		
		/* e) solve */
		PetscPrintf(PETSC_COMM_WORLD,"   [[ COMPUTING FLOW FIELD FOR STEP : %D ]]\n", user->step );
		ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);
		ierr = pTatinLogBasicSNES(user,"Stokes",snes);CHKERRQ(ierr);
		
		
		/* output */
		if ( (user->step%user->output_frequency == 0) || (user->step == 1) ) {
			sprintf(stepname,"step%1.6d",user->step);
			ierr = pTatinModel_Output(user->model,user,X,stepname);CHKERRQ(ierr);
		}
		
		/* compute timestep */
		user->dt = 1.0e32;
		ierr = DMCompositeGetAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		ierr = SwarmUpdatePosition_ComputeCourantStep(dav_hierarchy[nlevels-1],velocity,&timestep);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&velocity,&pressure);CHKERRQ(ierr);
		ierr = pTatin_SetTimestep(user,"StkCourant",timestep);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"  timestep[%d] dt_courant = %1.4e \n", user->step,user->dt );

		/* CHECKPOINT */
		ierr = pTatin3dCheckpointManager(user,X);CHKERRQ(ierr);
		
		/* tidy up */
		for (k=0; k<nlevels; k++) {
			ierr = MatDestroy(&operatorA11[k]);CHKERRQ(ierr);
			ierr = MatDestroy(&operatorB11[k]);CHKERRQ(ierr);
		}
		ierr = MatDestroy(&A);CHKERRQ(ierr);
		ierr = MatDestroy(&B);CHKERRQ(ierr);
		ierr = SNESDestroy(&snes);CHKERRQ(ierr);

		
		/* update context time information */
		user->step++;
		/* Terminate time stepping */
		if (user->time >= user->time_max) {
			break;
		}
	}
	
	
	
	/* Clean up */
	for (k=0; k<nlevels-1; k++) {
		ierr = BCListDestroy(&u_bclist[k]);CHKERRQ(ierr);
	}
	for (k=0; k<nlevels-1; k++) {
		ierr = QuadratureDestroy(&volQ[k]);CHKERRQ(ierr);
	}
	for (k=1; k<nlevels; k++) {
		ierr = MatDestroy(&interpolation_v[k]);CHKERRQ(ierr);
		ierr = MatDestroy(&interpolation_eta[k]);CHKERRQ(ierr);
	}
	for (k=0; k<nlevels; k++) {
		ierr = DMDestroy(&dav_hierarchy[k]);CHKERRQ(ierr);
	}
	
	ierr = ISDestroy(&is_stokes_field[0]);CHKERRQ(ierr);
	ierr = ISDestroy(&is_stokes_field[1]);CHKERRQ(ierr);
	ierr = PetscFree(is_stokes_field);CHKERRQ(ierr);
	
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = pTatin3dDestroyContext(&user);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	PetscBool restart,flg;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	restart = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-restart",&restart,&flg);CHKERRQ(ierr);
	if (!restart) {
		ierr = pTatin3d_linear_viscous_forward_model_driver(argc,argv);CHKERRQ(ierr);
	} else {
		ierr = pTatin3d_linear_viscous_forward_model_driver_RESTART(argc,argv);CHKERRQ(ierr);
	}
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

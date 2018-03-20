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
 **    filename:   energy_assembly.h
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

#ifndef __energy_assembly_h__
#define __energy_assembly_h__

#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscsnes.h"
#include "ptatin3d_defs.h"
#include "element_utils_q2.h"
#include "element_utils_q1.h"
#include "QPntVolCoefEnergy_def.h"

double AdvDiffResidualForceTerm_UpwindXiExact(double pecletNumber);
double AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption(double pecletNumber);
double AdvDiffResidualForceTerm_UpwindXiCriticalAssumption(double pecletNumber);
PetscErrorCode AdvDiff3dComputeAverageCellSize( PetscScalar el_coords[], PetscScalar DX[]);
PetscErrorCode AdvDiff3dComputeElementPecletNumber_qp( PetscScalar el_coords[],PetscScalar u[],
																										 PetscScalar kappa_el,
																										 PetscScalar *alpha);
PetscErrorCode AdvDiff3dComputeElementTimestep_qp(PetscScalar el_coords[],PetscScalar u[],PetscReal kappa_el,PetscReal *dta,PetscReal *dtd);
PetscErrorCode DASUPG3dComputeElementTimestep_qp(PetscScalar el_coords[],PetscScalar u[],PetscReal kappa_el,PetscReal *dta,PetscReal *dtd);
PetscErrorCode DASUPG3dComputeElementStreamlineDiffusion_qp(PetscScalar el_coords[],PetscScalar u[],
																														PetscInt nqp, PetscScalar qp_detJ[], PetscScalar qp_w[],
																														PetscScalar qp_kappa[],
																														PetscScalar *khat);
void ConstructNiSUPG_Q1_3D(PetscScalar Up[],PetscScalar kappa_hat,PetscScalar Ni[],PetscScalar GNx[NSD][NODES_PER_EL_Q1_3D],PetscScalar Ni_supg[]);


PetscErrorCode AElement_FormJacobian_T_supg( PetscScalar Re[],PetscReal dt,PetscScalar el_coords[],
																			 PetscScalar gp_kappa[],
																			 PetscScalar el_V[],
																			 PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] );

PetscErrorCode AElement_FormFunction_T_supg(
																	PetscScalar Re[],
																	PetscReal dt,
																	PetscScalar el_coords[],
																	PetscScalar el_coords_old[],
																	PetscScalar el_V[],
																	PetscScalar el_phi[],PetscScalar el_phi_last[],
																	PetscScalar gp_kappa[],PetscScalar gp_Q[],
																	PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] );

PetscErrorCode AdvDiffComputeTau_BrooksHughes(PetscScalar el_coords[],PetscScalar el_vel[],PetscScalar kappa_el,PetscScalar *tau);
PetscErrorCode AdvDiffComputeTau_TezduyarOsawa(PetscScalar el_coords[],PetscScalar el_vel[],PetscScalar kappa_cell,PetscScalar theta,PetscScalar dt,PetscScalar *tau);
PetscErrorCode AdvDiffComputeTau_UserDefinedConstant(PetscScalar const_t,PetscScalar *tau);

PetscErrorCode AElement_FormJacobian_T_tau(
																					 PetscScalar Re[],PetscReal dt,
																					 PetscScalar tau,
																					 PetscScalar el_coords[],
																					 PetscScalar gp_kappa[],
																					 PetscScalar el_V[],
																					 PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] );

PetscErrorCode AElement_FormFunction_T_tau(
																 PetscScalar Re[],PetscReal dt,
																 PetscReal tau,
																 PetscScalar el_coords[],
																 PetscScalar el_coords_old[],
																 PetscScalar el_V[],
																 PetscScalar el_phi[],PetscScalar el_phi_old[],
																 PetscScalar gp_kappa[],PetscScalar gp_Q[],
																 PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] );

PetscErrorCode TS_FormJacobianEnergy(PetscReal time,Vec X,PetscReal dt,Mat A,Mat B,void *ctx);
PetscErrorCode TS_FormFunctionEnergy(PetscReal time,Vec X,PetscReal dt,Vec F,void *ctx);

PetscErrorCode SNES_FormJacobianEnergy(SNES snes,Vec X,Mat A,Mat B,void *ctx);
PetscErrorCode SNES_FormFunctionEnergy(SNES snes,Vec X,Vec F,void *ctx);

#endif

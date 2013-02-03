
#ifndef __ptatin3d_energy_h__
#define __ptatin3d_energy_h__

#include "QPntVolCoefEnergy_def.h"

PetscErrorCode AdvDiff_GetElementLocalIndices_Q1(PetscInt el_localIndices[],PetscInt elnid[]);

double AdvDiffResidualForceTerm_UpwindXiExact(double pecletNumber);
double AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption(double pecletNumber);
double AdvDiffResidualForceTerm_UpwindXiCriticalAssumption(double pecletNumber);
PetscErrorCode DASUPG3dComputeAverageCellSize( PetscScalar el_coords[], PetscScalar DX[]);
PetscErrorCode DASUPG3dComputeElementPecletNumber_qp( PetscScalar el_coords[],PetscScalar u[],
																										 PetscScalar kappa_el,
																										 PetscScalar *alpha);
PetscErrorCode DASUPG3dComputeElementTimestep_qp(PetscScalar el_coords[],PetscScalar u[],PetscScalar kappa_el,PetscScalar *dta,PetscScalar *dtd);
PetscErrorCode DASUPG3dComputeElementStreamlineDiffusion_qp(PetscScalar el_coords[],PetscScalar u[],
																														PetscInt nqp, PetscScalar qp_detJ[], PetscScalar qp_w[],
																														PetscScalar qp_kappa[],
																														PetscScalar *khat);
void ConstructNiSUPG_Q1_3D(PetscScalar Up[],PetscScalar kappa_hat,PetscScalar Ni[],PetscScalar GNx[NSD][NODES_PER_EL_Q1_3D],PetscScalar Ni_supg[]);
PetscErrorCode AElement_SUPG3d_qp( PetscScalar Re[],PetscReal dt,PetscScalar el_coords[],
																	PetscScalar gp_kappa[],
																	PetscScalar el_V[],
																	PetscInt ngp,PetscScalar gp_xi[],PetscScalar gp_weight[] );
PetscErrorCode SUPGFormJacobian_qp(PetscReal time,Vec X,PetscReal dt,Mat *A,Mat *B,MatStructure *mstr,void *ctx);

PetscErrorCode VolumeQuadratureCreate_GaussLegendreEnergy(PetscInt nsd,PetscInt np_per_dim,PetscInt ncells,Quadrature *quadrature);
PetscErrorCode VolumeQuadratureGetAllCellData_Energy(Quadrature Q,QPntVolCoefEnergy *coeffs[]);
PetscErrorCode VolumeQuadratureGetCellData_Energy(Quadrature Q,QPntVolCoefEnergy coeffs[],PetscInt cidx,QPntVolCoefEnergy *cell[]);

#endif

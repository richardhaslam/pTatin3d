
#ifndef __MATERIAL_POINT_STD_UTILS_H__
#define __MATERIAL_POINT_STD_UTILS_H__

#include "swarm_fields.h"
#include "data_exchanger.h"
#include "MPntStd_def.h"

typedef enum { _FALSE=0, _TRUE=1 } Truth;

PetscErrorCode SwarmMPntStd_AssignUniquePointIdentifiers(MPI_Comm comm,DataBucket db,int start_pid,int end_pid);
PetscErrorCode SwarmMPntStd_CoordAssignment_LatticeLayout3d(DM da,PetscInt Nxp[],PetscReal perturb,DataBucket db);
PetscErrorCode SwarmMPntStd_CoordAssignment_RandomLayout3d(DM da,PetscInt nPerCell,DataBucket db);
PetscErrorCode SwarmOutputParaView_MPntStd(DataBucket db,const char path[],const char prefix[]);

/*
PetscErrorCode SwarmUpdatePosition_MPntStd_Euler(DM da,Vec velocity,PetscReal step,int npoints,MPntStd marker[]);

void InverseMappingDomain_2dQ2( 
															 double tolerance, int max_its,
															 Truth use_nonzero_guess, 
															 Truth monitor, Truth log,
															 const double coords[], const int mx, const int my, const int element[],
															 int np, MPntStd marker[] );

PetscErrorCode SwarmUpdatePosition_ComputeCourantStep(DM da,Vec velocity,PetscReal *step);
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep_2(DM da,Vec velocity,PetscReal *step);
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep_3(DM da,Vec velocity,PetscReal *step);

PetscErrorCode SwarmUpdateProperties_MPntStd(DataBucket db,pTatinCtx ctx,Vec X);
*/

#endif

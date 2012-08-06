
#ifndef __MATERIAL_POINT_STD_UTILS_H__
#define __MATERIAL_POINT_STD_UTILS_H__

#include "swarm_fields.h"
#include "data_exchanger.h"
#include "MPntStd_def.h"


PetscErrorCode SwarmMPntStd_AssignUniquePointIdentifiers(MPI_Comm comm,DataBucket db,int start_pid,int end_pid);
PetscErrorCode SwarmMPntStd_CoordAssignment_LatticeLayout3d(DM da,PetscInt Nxp[],PetscReal perturb,DataBucket db);
PetscErrorCode SwarmMPntStd_CoordAssignment_RandomLayout3d(DM da,PetscInt nPerCell,DataBucket db);
PetscErrorCode SwarmOutputParaView_MPntStd(DataBucket db,const char path[],const char prefix[]);

#endif

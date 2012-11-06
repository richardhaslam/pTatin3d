
#ifndef __ptatin3d_mesh_update_h__
#define __ptatin3d_mesh_update_h__

PetscErrorCode DMDABilinearizeQ2Elements(DM dau);
PetscErrorCode UpdateMeshGeometry_FullLagrangian(DM dav,Vec velocity,PetscReal step);

#endif


#if !defined(_DAFEIMPL_H)
#define _DAFEIMPL_H

#include <petscvec.h> /*I "petscvec.h" I*/
#include <petscmat.h>       /*I      "petscmat.h"          I*/
#include <petscdm.h>
#include <petsc/private/dmimpl.h>
#include <dafe.h>

struct _p_DAFEOps {
  PetscErrorCode (*EvaluateBasis)(PetscReal*,PetscReal*,PetscReal*);
  PetscErrorCode (*GetElementCoords)(DM,PetscInt,PetscReal*);
  PetscErrorCode (*GetElementNodeMap)(DM,PetscInt,const PetscInt**);
  PetscErrorCode (*GetElementBasisMap)(DM,PetscInt,const PetscInt**);
  PetscErrorCode (*GetParentElementIndex)(DM,PetscInt,PetscInt*);
  PetscErrorCode (*ConvertChildLocalCoordinate)(DM,PetscInt,PetscReal*,PetscReal*);
  PetscErrorCode (*ProjectCoordinates)(DM,DM);
};

typedef struct {
  DAFEType dafetype;
  DM da,parent_dafe;

  PetscInt ncomponents; /* e.g. scalar = 1 or vector = 2,3 */
  
  PetscInt nodes_per_el,nbasis; /* number of nodes (for coords) per elemenet, number of basis functions per element */
  PetscInt mx,my,mz; /* global number of elements */
  PetscInt lmx,lmy,lmz; /* local number of elements */
  PetscInt *lmx_range,*lmy_range,*lmz_range; /* owenership range of local elements */
  PetscInt corner_imin,corner_jmin,corner_kmin; /* index span of elements in terms of DA points */
  //PetscInt corner_imax,corner_jmax,corner_kmax;
  PetscInt *corner_imin_range,*corner_jmin_range,*corner_kmin_range;
  PetscInt lnelements,nelements; /* local, global number of elements */
  PetscInt *element_node_map;
  PetscInt *element_basis_map;
  
  DAFEOps ops;
} DM_DAFE;

PetscErrorCode _DAFE_GetElementsPk_3D(DM dm,PetscInt *_npe,PetscInt *_nel,PetscInt **_eidx);
PetscErrorCode _DAFE_CreatePkFromQ2_3d(DM dm,PetscInt ref,PetscInt nbasis);

PetscErrorCode DAFE_CreateQ2_3d(DM dm,PetscInt mi,PetscInt mj,PetscInt mk);
PetscErrorCode DAFE_CreateP0FromQ2_3d(DM dm,PetscInt ref);
PetscErrorCode DAFE_CreateP1FromQ2_3d(DM dm,PetscInt ref);
PetscErrorCode DAFE_CreateQ1FromQ2_3d(DM dm,PetscInt ref);


#endif /* _DAFEIMPL_H */



#if !defined(_DAFE_H)
#define _DAFE_H


#include <petscdm.h>

typedef enum {
  DAFE_Q2=0,
  DAFE_Q1,
  DAFE_Q2sub, /* The XXXsub variants were created from either a Q2, or Q1 space */
  DAFE_Q1sub,
  DAFE_P1sub,
  DAFE_P0sub
} DAFEType;

extern const char *DAFETypeNames[];

typedef struct _p_DAFEOps *DAFEOps;

PetscErrorCode DMDAFESetType(DM dm,DAFEType type);
PetscErrorCode DMDAFESetNumberOfComponents(DM dm,PetscInt nc);
PetscErrorCode DMDAFESetUniformCoordinates(DM dm,PetscReal xmin,PetscReal xmax,PetscReal ymin,PetscReal ymax,PetscReal zmin,PetscReal zmax);

PetscErrorCode DMDAFEGetNumberOfComponents(DM dm,PetscInt *d);
PetscErrorCode DMDAFEGetElementSize(DM dm,PetscInt *ne,PetscInt r[]);
PetscErrorCode DMDAFEGetLocalElementSize(DM dm,PetscInt *ne,PetscInt r[]);
PetscErrorCode DMDAFEGetElements(DM dm,PetscInt *ne,PetscInt *np,const PetscInt **map);
PetscErrorCode DMDAFEGetBasis(DM dm,PetscInt *ne,PetscInt *np,const PetscInt **map);

PetscErrorCode DMDAFEGetDA(DM dm,DM *da);
PetscErrorCode DMDAFEGetParentDM(DM dm,DM *parent);
PetscErrorCode DMDAFEGetOps(DM dm,DAFEOps *ops);

PetscErrorCode DMDAFECreate3d(MPI_Comm comm,PetscInt ncomponents,DAFEType dafetype,PetscInt mi,PetscInt mj,PetscInt mk,DM *_dm);
PetscErrorCode DMDAFECreateSubDMDAFE(DM dm,PetscInt ncomponents,PetscInt refine_factor,DAFEType type,DM *_dmsub);

PetscErrorCode DMDAFEProjectCoordinates(DM sub);
PetscErrorCode DMDAFEInjectCoordinates(DM dmchild);


#endif

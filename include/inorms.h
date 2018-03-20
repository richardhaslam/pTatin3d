
#ifndef __inorms_h__
#define __inorms_h__

#include <petsc.h>
#include <petscvec.h>

PetscErrorCode Evaluate_uL2_uH1_pL2(DM dav,DM dap,Vec sv,Vec sp,Quadrature volQ,
                                    PetscErrorCode (*exact_v)(PetscReal*,PetscReal*,void*),
                                    PetscErrorCode (*exact_L)(PetscReal*,PetscReal*,void*),
                                    PetscErrorCode (*exact_p)(PetscReal*,PetscReal*,void*),
                                    void *data,
                                    PetscReal *ev_L2,PetscReal *ev_H1s,PetscReal *ev_H1,PetscReal *ep_L2);

PetscErrorCode Evaluate_uL2_symuH1_pL2(DM dav,DM dap,Vec sv,Vec sp,Quadrature volQ,
                                       PetscErrorCode (*exact_v)(PetscReal*,PetscReal*,void*),
                                       PetscErrorCode (*exact_E)(PetscReal*,PetscReal*,void*),
                                       PetscErrorCode (*exact_p)(PetscReal*,PetscReal*,void*),
                                       void *data,
                                       PetscReal *ev_L2,PetscReal *ev_H1s,PetscReal *ev_H1,PetscReal *ep_L2);


#endif

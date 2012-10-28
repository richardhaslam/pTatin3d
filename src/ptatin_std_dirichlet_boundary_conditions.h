
#ifndef _ptatin_std_dirichlet_boundary_conditions_h__
#define _ptatin_std_dirichlet_boundary_conditions_h__


typedef enum { FRONT_FACE=1, BACK_FACE, EAST_FACE, NORTH_FACE, SOUTH_FACE, WEST_FACE } BoundaryFaceType;

PetscErrorCode DirichletBC_FreeSlip(BCList list,DM dav,BoundaryFaceType face);
PetscErrorCode DirichletBC_ApplyNormalVelocity(BCList list,DM dav,BoundaryFaceType face,PetscReal v_normal);
PetscErrorCode DirichletBC_ApplyStrainRateExx(BCList list,DM dav,PetscReal exx_bc);
PetscErrorCode DirichletBC_ApplyStrainRateExy(BCList list,DM dav,PetscReal exy_bc);
PetscErrorCode DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningY(BCList list,DM dav,PetscReal vx_bc);
PetscErrorCode DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ(BCList list,DM dav,PetscReal vx_bc);
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionX(BCList list,DM dav,PetscReal vx_bc);
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(BCList list,DM dav,PetscReal beta,PetscReal vx_bc);


#endif


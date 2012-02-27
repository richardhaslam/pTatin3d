

#include "petsc.h"
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_element_q2p1.h"
#include "quadrature.h"



#undef __FUNCT__
#define __FUNCT__ "StokesVelocity_GetElementLocalIndices"
PetscErrorCode StokesVelocity_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<27; n++) {
		el_localIndices[3*n  ] = 3*elnid[n]  ;
		el_localIndices[3*n+1] = 3*elnid[n]+1;
		el_localIndices[3*n+2] = 3*elnid[n]+2;
	}
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "StokesPressure_GetElementLocalIndices"
PetscErrorCode StokesPressure_GetElementLocalIndices(PetscInt el_localIndices[],PetscInt elnid[])
{
	PetscInt n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	for (n=0; n<P_BASIS_FUNCTIONS; n++) {
		el_localIndices[n] = elnid[n];
	}
	PetscFunctionReturn(0);
}


/* physics component loader */
#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreate_Stokes"
PetscErrorCode PhysCompCreate_Stokes(PhysCompStokes *ctx)
{
	PetscErrorCode ierr;
	PhysCompStokes stokes;
	
	PetscFunctionBegin;
	ierr = PetscMalloc(sizeof(struct _p_PhysCompStokes),&stokes);CHKERRQ(ierr);
	*ctx = stokes;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompDestroy_Stokes"
PetscErrorCode PhysCompDestroy_Stokes(PhysCompStokes *ctx)
{
	PetscErrorCode ierr;
	PhysCompStokes user;
	PetscInt e;
	
	PetscFunctionBegin;
	
	if (!ctx) {PetscFunctionReturn(0);}
	user = *ctx;
	
//	for (e=0; e<QUAD_EDGES; e++) {
//		if (user->surfQ[e]) { ierr = SurfaceQuadratureStokesDestroy(&user->surfQ[e]);CHKERRQ(ierr); }
//	}
//	if (user->Q) { ierr = QuadratureStokesDestroy(&user->Q);CHKERRQ(ierr); }
	if (user->p_bclist) { ierr = BCListDestroy(&user->p_bclist);CHKERRQ(ierr); }
	if (user->u_bclist) { ierr = BCListDestroy(&user->u_bclist);CHKERRQ(ierr); }
  if (user->stokes_pack) { ierr = DMDestroy(&user->stokes_pack);CHKERRQ(ierr); }
	if (user->dap) { ierr = DMDestroy(&user->dap);CHKERRQ(ierr); }
  if (user->dav) { ierr = DMDestroy(&user->dav);CHKERRQ(ierr); }
	if (user) { ierr = PetscFree(user);CHKERRQ(ierr); }
	
	*ctx = PETSC_NULL;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateMesh_Stokes3d"
PetscErrorCode PhysCompCreateMesh_Stokes3d(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx)
{
	DM dav,dap,multipys_pack;
	PetscInt vbasis_dofs;
	PetscInt pbasis_dofs;
	const PetscInt *lxp,*lyp,*lzp;
	PetscInt MX,MY,MZ,p,Mp,Np,Pp,*lxv,*lyv,*lzv,i;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	MX = mx;
	MY = my;
	MZ = mz;
	
	/* pressure */
	pbasis_dofs = P_BASIS_FUNCTIONS;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX,MY,MZ, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, pbasis_dofs,0, 0,0,0, &dap );CHKERRQ(ierr);
	ierr = DMDASetElementType_P1(dap);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(dap,&lxp,&lyp,&lzp);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dap,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Mp+1),&lxv);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Np+1),&lyv);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(Pp+1),&lzv);CHKERRQ(ierr);
	for (p=0; p<Mp; p++) {
		lxv[p] = lxp[p] * 2;
	} lxv[Mp-1]++;
	for (p=0; p<Np; p++) {
		lyv[p] = lyp[p] * 2;
	} lyv[Np-1]++;
	for (p=0; p<Pp; p++) {
		lzv[p] = lzp[p] * 2;
	} lzv[Pp-1]++;
	
	/* velocity */
	vbasis_dofs = 3;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, 2*MX+1,2*MY+1,2*MZ+1, Mp,Np,Pp, vbasis_dofs,2, lxv,lyv,lzv, &dav );CHKERRQ(ierr);
	ierr = DMDASetElementType_Q2(dav);CHKERRQ(ierr);
	ierr = PetscFree(lxv);CHKERRQ(ierr);
	ierr = PetscFree(lyv);CHKERRQ(ierr);
	ierr = PetscFree(lzv);CHKERRQ(ierr);
	
	/* set an initial geometry */
	ierr = DMDASetUniformCoordinates(dav,0.0,1.0, 0.0,1.0, 0.0,1.0);CHKERRQ(ierr);
	
	/* stokes */
	ierr = DMCompositeCreate(PETSC_COMM_WORLD,&multipys_pack);CHKERRQ(ierr);
	ierr = DMCompositeAddDM(multipys_pack,dav);CHKERRQ(ierr);	
	ierr = DMCompositeAddDM(multipys_pack,dap);CHKERRQ(ierr);	
	
	ierr = DMDASetFieldName(dav,0,"ux");CHKERRQ(ierr);
	ierr = DMDASetFieldName(dav,1,"uy");CHKERRQ(ierr);
	ierr = DMDASetFieldName(dav,2,"uz");CHKERRQ(ierr);
	switch (P_BASIS_FUNCTIONS) {
		case 1:
			ierr = DMDASetFieldName(dap,0,"P0_p");CHKERRQ(ierr);
			break;
		case 4:
			ierr = DMDASetFieldName(dap,0,"P1_p");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,1,"P1_dpdx");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,2,"P1_dpdy");CHKERRQ(ierr);
			ierr = DMDASetFieldName(dap,3,"P1_dpdz");CHKERRQ(ierr);
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Pressure space may ONLY contain 1 or 4 basis functions");
			break;
	}
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dav,"stk_velocity_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)dap,"stk_pressure_");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)multipys_pack,"stk_pack_");CHKERRQ(ierr);
	
	ctx->dav  = dav;
	ctx->dap  = dap;
	ctx->stokes_pack = multipys_pack;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateBoundaryList_Stokes"
PetscErrorCode PhysCompCreateBoundaryList_Stokes(PhysCompStokes ctx)
{
	DM dav;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* vel bc's */
	dav = ctx->dav;
	if (!dav) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"dav must be set"); }
	
	ierr = DMDABCListCreate(dav,&ctx->u_bclist);CHKERRQ(ierr);
	
	/* pressure bc's */
	ctx->p_bclist = PETSC_NULL;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompCreateVolumeQuadrature_Stokes"
PetscErrorCode PhysCompCreateVolumeQuadrature_Stokes(PhysCompStokes ctx)
{
	DM dav;
	PetscInt lmx,lmy,lmz;
	PetscInt np_per_dim,ncells;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	dav = ctx->dav;

	np_per_dim = 3;
  ierr = DMDAGetLocalSizeElementQ2(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	ncells = lmx * lmy * lmz;
	ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&ctx->volQ);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetValuesLocalStencil_AddValues_Stokes_Velocity"
PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Velocity(PetscScalar *fields_F,PetscInt u_eqn[],PetscScalar Fe_u[])
{
  PetscInt n,idx;
	
  PetscFunctionBegin;
  for (n = 0; n<U_BASIS_FUNCTIONS; n++) {
		idx = u_eqn[3*n  ];
    fields_F[idx] += Fe_u[NSD*n  ];
		
		idx = u_eqn[3*n+1];
    fields_F[idx] += Fe_u[NSD*n+1];

		idx = u_eqn[3*n+2];
    fields_F[idx] += Fe_u[NSD*n+2];
}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDASetValuesLocalStencil_AddValues_Stokes_Pressure"
PetscErrorCode DMDASetValuesLocalStencil_AddValues_Stokes_Pressure(PetscScalar *fields_F,PetscInt p_eqn[],PetscScalar Fe_p[])
{
  PetscInt n,idx;
	
  PetscFunctionBegin;
  for (n = 0; n<P_BASIS_FUNCTIONS; n++) {
		idx = p_eqn[n];
    fields_F[idx] += Fe_p[n];
  }
  PetscFunctionReturn(0);
}



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
 **    filename:   ptatin3d_stokes_q1macrop1.c
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


#include "petsc.h"
#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1macrop1.h"
#include "quadrature.h"
#include "dmda_checkpoint.h"



PetscErrorCode PhysCompCreateMesh_Stokes3d_Q1MacroP1(const PetscInt mx,const PetscInt my,const PetscInt mz,PhysCompStokes ctx)
{
	DM dav,dap,multipys_pack;
	PetscInt vbasis_dofs;
	PetscInt pbasis_dofs;
	const PetscInt *lxp,*lyp,*lzp;
	PetscInt MX,MY,MZ,Mp,Np,Pp,*lxv,*lyv,*lzv;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	MX = mx;
	MY = my;
	MZ = mz;

	/* velocity */
	vbasis_dofs = 3;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, 2*MX+1,2*MY+1,2*MZ+1, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, vbasis_dofs,1, NULL,NULL,NULL, &dav );CHKERRQ(ierr);
  ierr = DMSetUp(dav);CHKERRQ(ierr);
	ierr = DMDAESetType_Q1Macro(dav);CHKERRQ(ierr);

	ierr = DMDAGetInfo(dav,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAEQ1Macro_MixedSpace_GetOwnershipRangesElement(dav, 0,0,0, 0,0,0, &lxv,&lyv,&lzv);CHKERRQ(ierr);

	/* pressure */
	pbasis_dofs = P_BASIS_FUNCTIONS;
	ierr = DMDACreate3d( PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX,MY,MZ, Mp,Np,Pp, pbasis_dofs,0, lxv,lyv,lzv, &dap );CHKERRQ(ierr);
  ierr = DMSetUp(dap);CHKERRQ(ierr);
	ierr = DMDASetElementType_P1(dap);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(dap,&lxp,&lyp,&lzp);CHKERRQ(ierr);

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

	ierr = PetscFree(lxv);CHKERRQ(ierr);
	ierr = PetscFree(lyv);CHKERRQ(ierr);
	ierr = PetscFree(lzv);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/*
 For efficiency, i need a special quadrature rule to evaluate the weak form
 associated with the mixed terms involving u,p.
 The special rule consists of mapping the 2x2x2 quadrature rule used on the
 velocity cells into the local coordinate system used by the macro element.

*/
PetscErrorCode PhysCompCreateVolumeQuadrature_Stokes_Q1MacroP1(PhysCompStokes ctx)
{
	DM dav;
	PetscInt lmx,lmy,lmz;
	PetscInt np_per_dim,ncells;
	PetscErrorCode ierr;

	PetscFunctionBegin;

	dav = ctx->dav;

	np_per_dim = 2;
  ierr = DMDAEQ1Macro_NaturalSpace_GetLocalSizeElement(dav,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	ncells = lmx * lmy * lmz;
	ierr = VolumeQuadratureCreate_GaussLegendreStokes(3,np_per_dim,ncells,&ctx->volQ);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/* re-usable functions */
/*
PetscErrorCode PhysCompCreateBoundaryList_Stokes(PhysCompStokes ctx);
*/


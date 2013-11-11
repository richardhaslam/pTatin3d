/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      iplus_ops.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id:$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/


/*
 Lx,Ly,Lz
 
 eta {mantle,plume,slab}
   83, 5, ?

 rho {mantle,plume,slab}
   1413, 1373, ?
 
 typical A0,r0 we should compare with
 
*/

#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "dmda_element_q2p1.h"
#include "mesh_update.h"
#include "ptatin_models.h"
#include "model_utils.h"
#include "ptatin_utils.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"

#include "iplus_ctx.h"


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_iPLUS"
PetscErrorCode ModelInitialize_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscBool        flg;
	char             logfile[PETSC_MAX_PATH_LEN];
	PetscErrorCode   ierr;
	
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	data->mantle_eta = 85.0;    data->mantle_rho = 1413.0;
	data->plume_eta  = 5.0;     data->plume_rho  = 1373.0;
	data->slab_eta   = 0.0;     data->slab_rho   = 0.0;
		
	data->plume_pos[0] = 0.15;
	data->plume_pos[2] = 0.15;
	data->plume_pos[1] = 0.0;
	
	data->plume_radius = 0.015625;
	data->plume_A0     = 1.0681;
	PetscOptionsGetReal(PETSC_NULL,"-iplus_plume_r",&data->plume_radius,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-iplus_plume_A0",&data->plume_A0,&flg);

	data->refinement_type = 0;
	
	data->np_plume_x = 10;
	data->np_plume_z = 10;
	PetscOptionsGetInt(PETSC_NULL,"-iplus_plume_npx",&data->np_plume_x,&flg);
	PetscOptionsGetInt(PETSC_NULL,"-iplus_plume_npz",&data->np_plume_z,&flg);
	
	sprintf(logfile,"%s/iplus.logfile",c->outputpath);
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,logfile,&data->logviewer);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUS_DetermineElementsContainingPlumeInlet"
PetscErrorCode iPLUS_DetermineElementsContainingPlumeInlet(DM dav,iPLUSCtx *data)
{
	MPI_Comm       comm;
	PetscInt       rank,M,N,P,pI,pJ,pK,rI,rJ,rK,rIJ;
	PetscMPIInt    _rank;
	const PetscInt *elnidx_u;
	PetscInt       nel,nen_u,e,k;
	DM             cda;
	Vec            gcoords;
	PetscReal      *LA_gcoords;
	PetscReal      elcoords[3*Q2_NODES_PER_EL_3D];
	PetscErrorCode ierr;

	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* determine ranks living on the base j=0 */
	PetscObjectGetComm((PetscObject)dav,&comm);
	ierr = MPI_Comm_rank(comm,&_rank);CHKERRQ(ierr);
	rank = (PetscInt)_rank;
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,&pI,&pJ,&pK,0,0,0,0,0,0);CHKERRQ(ierr);
	
	/* convert rank to rI,rJ,rK */
	rK  = rank / (pI*pJ);
	rIJ = rank - rK * pI*pJ;  
	rJ = rIJ / pI;
	rI = rIJ - rJ*pI;
	
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dav,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);

	data->nplume_elements = 0;
	ierr = PetscMalloc(sizeof(PetscInt)*nel,&data->plume_element);CHKERRQ(ierr);
	ierr = PetscMemzero(data->plume_element,sizeof(PetscInt)*nel);CHKERRQ(ierr);
	
	if (rJ == 0) {
		for (e=0; e<nel; e++) {
			PetscReal el_y_min;
			PetscBool plume_in_element = PETSC_FALSE;
			
			ierr = DMDAGetElementCoordinatesQ2_3D(elcoords,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);

			/* check the lowest nodes in element */
			el_y_min = elcoords[3*0+1];
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
				PetscReal sep_xz2,r02,*pos;
				
				pos      = &elcoords[3*k];
				sep_xz2  = (pos[0] - data->plume_pos[0])*(pos[0] - data->plume_pos[0]);
				sep_xz2 += (pos[2] - data->plume_pos[2])*(pos[2] - data->plume_pos[2]);
				r02      = data->plume_pos[0]*data->plume_pos[0] + data->plume_pos[1]*data->plume_pos[1];
				
				if (el_y_min < 1.0e-10) {
					if (sep_xz2 <= r02) {
						plume_in_element = PETSC_TRUE;
					}
				}
			}
			
			if (plume_in_element) {
				data->plume_element[ data->nplume_elements ] = e;
				data->nplume_elements++;
			}					
			
		}
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

	if (data->nplume_elements) {
		PetscPrintf(PETSC_COMM_SELF,"  [%D] Located %D elements within the plume inlet\n",rank,data->nplume_elements);
	}
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_iPLUS"
PetscErrorCode ModelApplyInitialMeshGeometry_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;
	
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMDASetUniformCoordinates(dav,0.0,0.3,0.0,0.3,0.0,0.3);CHKERRQ(ierr);

	/* refine? */
	switch (data->refinement_type) {
		case 1:
			PetscPrintf(PETSC_COMM_WORLD,"Mesh refinement type: 1 (Subtle)\n");
			break;
		case 2:
			PetscPrintf(PETSC_COMM_WORLD,"Mesh refinement type: 2 (Aggressive)\n");
			break;
	}
	
	/* determine elements located within the plume */
	ierr = iPLUS_DetermineElementsContainingPlumeInlet(dav,data);CHKERRQ(ierr);

	/* compute initial volume of the domain */
	ierr = DMDAComputeMeshVolume(dav,&data->intial_domain_volume);CHKERRQ(ierr);	
	
	PetscFunctionReturn(0);
}

/*
 -find procs on base
 -srand(rank)
 -sweep bottom elements
 -create cell layout for each element in plume
 -insert points inside, initializing with plume values
*/
#undef __FUNCT__
#define __FUNCT__ "iPLUS_InsertPlumeMaterial"
PetscErrorCode iPLUS_InsertPlumeMaterial(DM dav,DataBucket materialpoint_db,iPLUSCtx *data)
{
	MPAccess         mpX;
	PetscInt         p,n_mpoints_orig,n_mpoints,new_points_max,new_points;
	const PetscInt   *elnidx_u;
	PetscInt         nel,nen_u,e,k;
	DM               cda;
	Vec              gcoords;
	PetscReal        *LA_gcoords;
	PetscReal        elcoords[3*Q2_NODES_PER_EL_3D];
	MPI_Comm         comm;
	PetscInt         eidx,pi,pk,np_plume_x,np_plume_z;
	PetscReal        dxi,dzeta,sep2;
  PetscScalar      el_coords[Q2_NODES_PER_EL_3D*NSD];
	PetscBool        inside_plume;
	PetscMPIInt      rank;
	PetscErrorCode   ierr;

	
	PetscFunctionBegin;
	
	if (data->nplume_elements == 0) {
		PetscFunctionReturn(0);
	}
	
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dav,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);

	/* 
	 At most we will need to insert 
	   nplume_elements . plume_x . plume_z
	 new material points 
	*/
	np_plume_x = data->np_plume_x;
	np_plume_z = data->np_plume_z;
	dxi    = 2.0/(PetscReal)np_plume_x;
	dzeta  = 2.0/(PetscReal)np_plume_z;

	DataBucketGetSizes(materialpoint_db,&n_mpoints_orig,0,0);
	new_points_max = data->nplume_elements * np_plume_x * np_plume_z;
	DataBucketSetSizes(materialpoint_db,n_mpoints_orig+new_points_max,-1);
	
	
	PetscObjectGetComm((PetscObject)dav,&comm);
	ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	ptatin_RandomNumberSetSeedRank(comm);

	new_points = 0;
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (e=0; e<data->nplume_elements; e++) {

		eidx = data->plume_element[ e ];
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx_u[nen_u*eidx],LA_gcoords);CHKERRQ(ierr);

		for (pk=0; pk<np_plume_z; pk++) {
			for (pi=0; pi<np_plume_x; pi++) {
				double xip[3],xp[NSD],Ni[Q2_NODES_PER_EL_3D];
				
				xip[0] = -1.0 + dxi    * (pi + 0.5);
				xip[1] = -1.0;
				xip[2] = -1.0 + dzeta  * (pk + 0.5);
		
				pTatin_ConstructNi_Q2_3D(xip,Ni);
				
				xp[0] = xp[1] = xp[2] = 0.0;
				for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
					xp[0] += Ni[k] * el_coords[NSD*k+0];
					xp[1] += Ni[k] * el_coords[NSD*k+1];
					xp[2] += Ni[k] * el_coords[NSD*k+2];
				}

				inside_plume = PETSC_FALSE;
				sep2  = (xp[0]-data->plume_pos[0])*(xp[0]-data->plume_pos[0]);
				sep2 += (xp[2]-data->plume_pos[2])*(xp[2]-data->plume_pos[2]);
				if (sep2 < data->plume_radius * data->plume_radius) {
					inside_plume = PETSC_TRUE;
				}
				
				if (inside_plume) {
					PetscInt mpidx;
					
					mpidx = n_mpoints_orig + new_points;
					
					ierr = MaterialPointSet_global_coord(mpX,mpidx,xp);CHKERRQ(ierr);
					ierr = MaterialPointSet_local_coord(mpX,mpidx,xip);CHKERRQ(ierr);
					ierr = MaterialPointSet_local_element_index(mpX,mpidx,eidx);CHKERRQ(ierr);
					
					new_points++;
				}
				
			}
		}
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	DataBucketSetSizes(materialpoint_db,n_mpoints_orig+new_points,-1);
	
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=n_mpoints_orig; p<n_mpoints; p++) {
		ierr = MaterialPointSet_phase_index(mpX,p,iPLUSMatPlume);CHKERRQ(ierr);
		ierr = MaterialPointSet_viscosity(mpX,p,data->plume_eta);CHKERRQ(ierr);
		ierr = MaterialPointSet_density(mpX,p,-GRAVITY*data->plume_rho);CHKERRQ(ierr);
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_SELF,"  [%D] added %D new plume points; total = %D\n",rank,new_points,n_mpoints);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_iPLUS"
PetscErrorCode ModelApplyInitialMaterialGeometry_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	DataBucket       materialpoint_db;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;

	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointSet_phase_index(mpX,p,iPLUSMatMantle);CHKERRQ(ierr);
		ierr = MaterialPointSet_viscosity(mpX,p,data->mantle_eta);CHKERRQ(ierr);
		ierr = MaterialPointSet_density(mpX,p,-GRAVITY*data->mantle_rho);CHKERRQ(ierr);
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = iPLUS_InsertPlumeMaterial(dav,materialpoint_db,data);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialSolution_iPLUS"
PetscErrorCode ModelApplyInitialSolution_iPLUS(pTatinCtx c,Vec X,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUS_Inflow_BCListEvaluator"
PetscBool iPLUS_Inflow_BCListEvaluator(PetscScalar position[],PetscScalar *value,void *ctx) 
{
	PetscBool      impose_dirichlet = PETSC_TRUE;
	iPLUSCtx       *model_data;
	PetscScalar    vy;
	PetscBool      inside_plume;
	PetscReal      sep2;
	PetscErrorCode ierr;

	
	PetscFunctionBegin;

	model_data = (iPLUSCtx*)ctx;
	
	inside_plume = PETSC_FALSE;
	sep2  = (position[0]-model_data->plume_pos[0])*(position[0]-model_data->plume_pos[0]);
	sep2 += (position[2]-model_data->plume_pos[2])*(position[2]-model_data->plume_pos[2]);
	if (sep2 < model_data->plume_radius * model_data->plume_radius) {
		inside_plume = PETSC_TRUE;
	}
	
	if (inside_plume) {
		PetscReal r2  = sep2;
		PetscReal r02 = model_data->plume_radius * model_data->plume_radius;
		
		vy = model_data->plume_A0 * ( r02 - r2 );
	} else {
		vy = 0.0;
	}
	*value = vy;
	return impose_dirichlet;
}

#undef __FUNCT__
#define __FUNCT__ "iPLUS_VelocityBC"
PetscErrorCode iPLUS_VelocityBC(BCList bclist,DM dav,pTatinCtx c,iPLUSCtx *data)
{
	PetscScalar    zero;
	PetscInt       d;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	zero = 0.0;

	/* south face */
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,iPLUS_Inflow_BCListEvaluator,(void*)data);CHKERRQ(ierr);

	/* no slip bottom */
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

	/* no slip */
	/*
	for (d=0; d<3; d++) {
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,d,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,d,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,d,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
		ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,d,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
	}
	*/
	ierr = DirichletBC_FreeSlip(bclist,dav,FRONT_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,EAST_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,BACK_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,WEST_FACE);CHKERRQ(ierr);

	
	/* north face is free surface */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_iPLUS"
PetscErrorCode ModelApplyBoundaryCondition_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	PetscErrorCode   ierr;

	
	PetscFunctionBegin;

	/* Define velocity boundary conditions */
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	ierr = iPLUS_VelocityBC(stokes->u_bclist,dav,c,data);CHKERRQ(ierr);

	/* Define boundary conditions for any other physics */
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_iPLUS"
PetscErrorCode ModelApplyBoundaryConditionMG_iPLUS(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscInt         n;
	PetscErrorCode   ierr;
	
	
	PetscFunctionBegin;
	
	/* Define velocity boundary conditions on each level within the MG hierarchy */
	for (n=0; n<nl; n++) {
		ierr = iPLUS_VelocityBC(bclist[n],dav[n],c,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_iPLUS"
PetscErrorCode ModelApplyMaterialBoundaryCondition_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	DataBucket       materialpoint_db;
	PetscReal        r02,vy,element_dy,frac_dy;
	PetscInt         N;
	static PetscReal displacement = 0.0;
	PetscErrorCode   ierr;
	
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	/* compute velocity at inlet */
	r02 = data->plume_radius * data->plume_radius;
	vy = data->plume_A0 * r02;
	/* record displacement */
	displacement += vy * c->dt;

	/* if displacement is larger than frac.dy, inject more points */
	ierr = DMDAGetInfo(dav,0,0,&N,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	element_dy = 0.3/((PetscReal)(N-1)/2);
	frac_dy = 0.1 * element_dy;
	printf("*** displ = %1.4e frac.dy = %1.4e \n",displacement,frac_dy);
	if (displacement < frac_dy) {
		PetscFunctionReturn(0);
	}
	

	ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);

	ierr = iPLUS_InsertPlumeMaterial(dav,materialpoint_db,data);CHKERRQ(ierr);

	/* reset displacement measure */
	displacement = 0.0;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_iPLUS"
PetscErrorCode ModelApplyUpdateMeshGeometry_iPLUS(pTatinCtx c,Vec X,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscReal        step;
	PhysCompStokes   stokes;
	DM               stokes_pack,dav,dap;
	Vec              velocity,pressure;
	PetscErrorCode   ierr;

	
	PetscFunctionBegin;
	
	ierr = pTatinGetTimestep(c,&step);CHKERRQ(ierr);
	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	ierr = UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(dav,velocity,step);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUSOutput_ComputeVerticalRangeOfRegion"
PetscErrorCode iPLUSOutput_ComputeVerticalRangeOfRegion(DataBucket materialpoint_db,PetscInt region_idx,PetscReal range_yp[])
{
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	PetscReal        _range_yp[2];
	double           *pos_p;
	int              region_p;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointGet_phase_index(mpX,p,&region_p);CHKERRQ(ierr);
		ierr = MaterialPointGet_global_coord(mpX,p,&pos_p);CHKERRQ(ierr);
		
		if (region_p == region_idx) {
			if (pos_p[1] < _range_yp[0]) { _range_yp[0] = pos_p[1]; }
			if (pos_p[1] > _range_yp[1]) { _range_yp[1] = pos_p[1]; }
		}
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&_range_yp[0],&range_yp[0],1,MPIU_REAL,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&_range_yp[1],&range_yp[1],1,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "iPLUSOutput_ComputeDomainVolume"
PetscErrorCode iPLUSOutput_ComputeDomainVolume(DM dav,PetscReal *volume)
{
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = DMDAComputeMeshVolume(dav,volume);CHKERRQ(ierr);	

	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelOutput_iPLUS"
PetscErrorCode ModelOutput_iPLUS(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	static int       beenhere = 0;
	FILE             *fp;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* ---- Velocity-Pressure Mesh Output ---- */
	/* [1] Standard viewer: v,p written out as binary in double */
	//ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

	/* [2] Light weight viewer: Only v is written out. v and coords are expressed as floats */
	ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(c,X,prefix);CHKERRQ(ierr);
	
	/* [3] Write out v,p into PETSc Vec. These can be used to restart pTatin */
	/*
	ierr = pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
	*/

	
	/* ---- Material Point Output ---- */
	/* [1] Basic viewer: Only reports coords, regionid and other internal data */
	ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
	
	/* [2] Customized viewer: User defines specific fields they want to view - NOTE not .pvd file will be created */
	/*
	{
		DataBucket                materialpoint_db;
		const int                 nf = 4;
		const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
		char                      mp_file_prefix[256];
		
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
		sprintf(mp_file_prefix,"%s_mpoints",prefix);
		ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
	}
	*/
	/* [3] Customized marker->cell viewer: Marker data is projected onto the velocity mesh. User defines specific fields */
	/*
	{
		const int                    nf = 3;
		const MaterialPointVariable  mp_prop_list[] = { MPV_viscosity, MPV_density, MPV_plastic_strain }; 
		
		ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
	}	
	*/
	
	/* iPlus specific output */
	{
		PhysCompStokes   stokes;
		DM               stokes_pack,dav,dap;
		DataBucket       materialpoint_db;
		PetscReal        volume,range_yp[2];
		
		ierr = pTatinGetMaterialPoints(c,&materialpoint_db,PETSC_NULL);CHKERRQ(ierr);
		
		ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
		stokes_pack = stokes->stokes_pack;
		ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
		
		ierr = iPLUSOutput_ComputeVerticalRangeOfRegion(materialpoint_db,iPLUSMatPlume,range_yp);
		ierr = iPLUSOutput_ComputeDomainVolume(dav,&volume);CHKERRQ(ierr);
		
		if (beenhere == 0) {
			PetscViewerASCIIPrintf(data->logviewer,"# iPLUS logfile\n");
			PetscViewerASCIIPrintf(data->logviewer,"# step \t time \t Omega0 \t Omega \t plume (y_min y_max) \n");
			beenhere = 1;
		}
		PetscViewerASCIIPrintf(data->logviewer,"%D\t%1.4e\t%1.6e\t%1.6e\t%1.4e\t%1.4e\n",
								c->step,c->time, data->intial_domain_volume, volume,range_yp[0],range_yp[1]);
		
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_iPLUS"
PetscErrorCode ModelDestroy_iPLUS(pTatinCtx c,void *ctx)
{
	iPLUSCtx         *data = (iPLUSCtx*)ctx;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	if (data->plume_element) {
		ierr = PetscFree(data->plume_element);CHKERRQ(ierr);
	}
	if (data->logviewer) {
		ierr = PetscViewerDestroy(&data->logviewer);CHKERRQ(ierr);
	}
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_iPLUS"
PetscErrorCode pTatinModelRegister_iPLUS(void)
{
	iPLUSCtx       *data;
	pTatinModel    m,model;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(iPLUSCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(iPLUSCtx));CHKERRQ(ierr);
	
	/* set initial values for model parameters */
	data->modeltype = iPLUsModelA;
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"iplus");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_iPLUS);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_iPLUS);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

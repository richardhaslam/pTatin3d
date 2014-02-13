
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
	/* scale inflow velocity */
	*value = vy * (1.0/ model_data->vel_scale); /* make velocity non-dimensional */
	return impose_dirichlet;
}

#undef __FUNCT__
#define __FUNCT__ "iPLUS_ApplyMaterialBoundaryCondition_Plume"
PetscErrorCode iPLUS_ApplyMaterialBoundaryCondition_Plume(pTatinCtx c,iPLUSCtx *data)
{
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
	displacement += vy * c->dt * data->time_scale; /* dimensionalize time, t = T' hat{dt}, T' = L'/V' */
	
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

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
 **    filename:   ptatin_driver_ic_spm.c
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

static const char help[] = "Stokes solver using Q2-Pm1 mixed finite elements.\n"
"3D prototype of the (p)ragmatic version of Tatin. (pTatin3d_v0.0)\n\n";


#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"

#include "material_point_utils.h"
#include "material_point_std_utils.h"
#include "ptatin_models.h"
#include "ptatin_utils.h"
#include "stokes_form_function.h"
#include "stokes_operators.h"
#include "dmda_redundant.h"
#include "dmda_update_coords.h"
#include "dmda_element_q1.h"
#include "dmda_view_petscvtk.h"

/*
 assumed y is the surface
 */
PetscErrorCode pTatinSurfaceMeshCreate(DM dav, DM *da_spm,Vec *_height)
{
	PetscErrorCode ierr;
	PetscInt si,sj,sk,nx,ny,nz,M,N,P;
	PetscInt si2d,sj2d,nx2d,ny2d,i,k;
	MPI_Comm comm;
	DM da_red_spm,da_surf;
	DM da_red_spm_coords,da_surf_coords;
	Vec coords_surf, coords_red_spm;
	DMDACoor3d ***LA_coords_red_spm;
	DMDACoor2d **LA_coords_surf;
	PetscScalar **LA_height;
	int rank;
	Vec height;
	
	PetscFunctionBegin;
	
	comm = PetscObjectComm((PetscObject)dav);
	ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDACreate2d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,M,P, PETSC_DECIDE,PETSC_DECIDE, 1,1, 0,0,&da_surf);CHKERRQ(ierr);	
  ierr = DMSetUp(da_surf);CHKERRQ(ierr);
	ierr = DMDASetUniformCoordinates(da_surf, 0.0,1.0, 0.0,1.0, 0.,0.);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da_surf,&si2d,&sj2d,0,&nx2d,&ny2d,0);CHKERRQ(ierr);

	/* fetch coordinates i need from the mechanical domain */
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);

	//printf("rank %d: vol  si,ei=%.4d-%.4d: sj,ej=%.4d-%.4d: sk,ek=%.4d-%.4d \n", rank, si,si+nx,sj,sj+ny,sk,sk+nz);
	//printf("rank %d: surf si,ei=%.4d-%.4d:                  sj,ej=%.4d-%.4d \n", rank, si2d,si2d+nx2d,sj2d,sj2d+ny2d);
	
	ierr = DMDACreate3dRedundant(dav,si2d,si2d+nx2d,N-1,N,sj2d,sj2d+ny2d, 1, &da_red_spm);CHKERRQ(ierr);
	//sprintf(name,"surf%d.vtk",rank);
	//ierr = DMDAViewPetscVTK(da_red_spm,0,name);CHKERRQ(ierr);
	
	/* copy these values into my parallel surface mesh x,y,z (vol) => x,y (surf) */
	ierr = DMGetCoordinates(da_surf,&coords_surf);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_red_spm,&coords_red_spm);CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(da_surf,&da_surf_coords);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(da_red_spm,&da_red_spm_coords);CHKERRQ(ierr);
	
	
	ierr = DMDAVecGetArray(da_surf_coords,coords_surf,&LA_coords_surf);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_red_spm_coords,coords_red_spm,&LA_coords_red_spm);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(da_surf,&si,&sk,0,&nx,&nz,0);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (i=si; i<si+nx; i++) {
			LA_coords_surf[k][i].x = LA_coords_red_spm[k-sk][0][i-si].x;
			LA_coords_surf[k][i].y = LA_coords_red_spm[k-sk][0][i-si].z;
		}
	}
	
	ierr = DMCreateGlobalVector(da_surf,&height);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_surf,height,&LA_height);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (i=si; i<si+nx; i++) {
			LA_height[k][i] = LA_coords_red_spm[k-sk][0][i-si].y;
		}
	}
	ierr = DMDAVecRestoreArray(da_surf,height,&LA_height);CHKERRQ(ierr);
	
	
	ierr = DMDAVecRestoreArray(da_red_spm_coords,coords_red_spm,&LA_coords_red_spm);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da_surf_coords,coords_surf,&LA_coords_surf);CHKERRQ(ierr);

	ierr = DMDestroy(&da_red_spm);CHKERRQ(ierr);

	ierr = DMDAUpdateGhostedCoordinates(da_surf);CHKERRQ(ierr);
	
	//sprintf(name,"surf.vtk");
	//ierr = DMDAViewPetscVTK(da_surf,height,name);CHKERRQ(ierr);
	
	ierr = DMDestroy(&da_red_spm);CHKERRQ(ierr);
	
	*_height = height;
	*da_spm = da_surf;
	
	PetscFunctionReturn(0);
}

/*
 parallel(nproces) surface mesh
 => redundant surface mesh overlapping volume mesh
 => insert onto volume mesh
*/

PetscErrorCode pTatin_InjectSurfaceMeshOntoMechanicalDomain(DM da_surf,Vec height,DM da_vol)
{
	PetscInt M,N,P,M2d,P2d;
	PetscInt si,sj,sk,nx,ny,nz;
	PetscInt Ml2d,Pl2d;
	PetscInt *surf_indices;
	Vec coords_vol;
	IS is_surf,is_local;
	VecScatter sctx;
	PetscInt c,i,k,j;
	DM cda_vol;
	Vec height_NATURAL,height_self;
	PetscScalar *LA_height_self;
	DMDACoor3d ***LA_coords_vol;
	MPI_Comm comm;
	int rank;
	
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	comm = PetscObjectComm((PetscObject)da_vol);
	ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	

	ierr = DMDAGetInfo(da_vol,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);

	/* determin pieces i need on the volume mesh */
	ierr = DMDAGetCorners(da_vol,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);

	/* create natural vector in i,j ordering for height */
	ierr = DMDACreateNaturalVector(da_surf,&height_NATURAL);CHKERRQ(ierr);
	ierr = DMDAGlobalToNaturalBegin( da_surf,height,INSERT_VALUES, height_NATURAL );CHKERRQ(ierr);
	ierr = DMDAGlobalToNaturalEnd(   da_surf,height,INSERT_VALUES, height_NATURAL );CHKERRQ(ierr);
	
	/* surface */
	/* identify the nodes on the surface of the mechanical domain */
	ierr = DMDAGetInfo(da_surf,0,&M2d,&P2d,0,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	Ml2d = nx;
	Pl2d = nz;
	ierr = PetscMalloc( sizeof(PetscInt)*Ml2d*Pl2d*1, &surf_indices );CHKERRQ(ierr);
	
	/* it would be nice if didn't have to scatter the surface to all cpus... */
	j = N - 1;
	c = 0;
//	if ((j>=sj) && (j<sj+ny)) {
		for( k=sk; k<sk+nz; k++ ){
			for( i=si; i<si+nx; i++ ){
				PetscInt nidx;
				
				nidx = (i) + (k)*M2d;
				surf_indices[c] = nidx; 
				c++;
			}
		}
//	}

	/* create local comm_self vector to store height */
	ierr = VecCreate(PETSC_COMM_SELF,&height_self);CHKERRQ(ierr);
	ierr = VecSetSizes(height_self,PETSC_DECIDE,nx*nz);CHKERRQ(ierr);
	ierr = VecSetFromOptions(height_self);CHKERRQ(ierr);
	
//	if(c!=0){	
	ierr = ISCreateGeneral( PETSC_COMM_WORLD, c, surf_indices, PETSC_USE_POINTER, &is_surf );CHKERRQ(ierr);
	ierr = ISCreateStride(  PETSC_COMM_SELF,  c,0,1,&is_local);CHKERRQ(ierr);
	ierr = VecScatterCreate( height_NATURAL,is_surf, height_self,is_local, &sctx );CHKERRQ(ierr);
	
	ierr = VecScatterBegin( sctx, height_NATURAL,height_self,INSERT_VALUES, SCATTER_FORWARD );CHKERRQ(ierr);
	ierr = VecScatterEnd(   sctx, height_NATURAL,height_self,INSERT_VALUES, SCATTER_FORWARD );CHKERRQ(ierr);
//	}	
	
	
	/* update mesh */
	ierr = DMGetCoordinates(da_vol,&coords_vol);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(da_vol,&cda_vol);CHKERRQ(ierr);
	
	ierr = DMDAVecGetArray(cda_vol,coords_vol,&LA_coords_vol);CHKERRQ(ierr);
	ierr = VecGetArray(height_self,&LA_height_self);CHKERRQ(ierr);
	
	/* identify the nodes on the surface of the mechanical domain */
	ierr = DMDAGetInfo(da_vol,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da_vol,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);	
	
	j = N - 1;
	if ((j>=sj) && (j<sj+ny)) {
		for (k=sk; k<sk+nz; k++) {
			for (i=si; i<si+nx; i++) {
				LA_coords_vol[k][j][i].y = LA_height_self[(i-si)+(k-sk)*nz];
			}
		}
	}	

	ierr = VecRestoreArray(height_self,&LA_height_self);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda_vol,coords_vol,&LA_coords_vol);CHKERRQ(ierr);
	
	/* insert surface values */
	ierr = DMDAUpdateGhostedCoordinates(da_vol);CHKERRQ(ierr);
	
//	if (c!=0){
	ierr = VecScatterDestroy(&sctx);CHKERRQ(ierr);
	ierr = ISDestroy(&is_surf);CHKERRQ(ierr);
	ierr = ISDestroy(&is_local);CHKERRQ(ierr);
//	}
	ierr = VecDestroy(&height_self);CHKERRQ(ierr);
	ierr = PetscFree(surf_indices);CHKERRQ(ierr);
	ierr = VecDestroy(&height_NATURAL);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* 
 maps (vol) x,y,z => (surf) x,y,height
*/
PetscErrorCode pTatin_InjectMechanicalDomainSurfaceOntoSurfaceMesh(DM da_vol,DM da_surf,Vec height)
{
	PetscInt si,sk,nx,nz,M,N,P;
	PetscInt si2d,sj2d,nx2d,ny2d,i,k;
	Vec coords_surf, coords_red_vol;
	DMDACoor3d ***LA_coords_red_vol;
	DMDACoor2d **LA_coords_surf;
	PetscScalar **LA_height;
	DM da_red_vol,cda_surf,cda_red_vol;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(da_vol,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da_surf,&si2d,&sj2d,0,&nx2d,&ny2d,0);CHKERRQ(ierr);
	ierr = DMDACreate3dRedundant(da_vol,si2d,si2d+nx2d,N-1,N,sj2d,sj2d+ny2d, 1, &da_red_vol);CHKERRQ(ierr);
	
	
	ierr = DMGetCoordinates(da_surf,&coords_surf);CHKERRQ(ierr);
	ierr = DMGetCoordinates(da_red_vol,&coords_red_vol);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(da_surf,&cda_surf);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(da_red_vol,&cda_red_vol);CHKERRQ(ierr);
	
	
	/* do the copy */
	ierr = DMDAVecGetArray(cda_surf,coords_surf,&LA_coords_surf);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_red_vol,coords_red_vol,&LA_coords_red_vol);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(da_surf,&si,&sk,0,&nx,&nz,0);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (i=si; i<si+nx; i++) {
			LA_coords_surf[k][i].x = LA_coords_red_vol[k-sk][0][i-si].x;
			LA_coords_surf[k][i].y = LA_coords_red_vol[k-sk][0][i-si].z;
		}
	}
	
	ierr = DMDAVecGetArray(da_surf,height,&LA_height);CHKERRQ(ierr);
	for (k=sk; k<sk+nz; k++) {
		for (i=si; i<si+nx; i++) {
			LA_height[k][i] = LA_coords_red_vol[k-sk][0][i-si].y;
		}
	}
	ierr = DMDAVecRestoreArray(da_surf,height,&LA_height);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda_red_vol,coords_red_vol,&LA_coords_red_vol);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda_surf,coords_surf,&LA_coords_surf);CHKERRQ(ierr);
	
	ierr = DMDestroy(&da_red_vol);CHKERRQ(ierr);
	
	ierr = DMDAUpdateGhostedCoordinates(da_surf);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}


PetscErrorCode pTatin3d_material_points_check_ic(int argc,char **argv)
{
	DM              multipys_pack,dav;
	PetscErrorCode  ierr;
	pTatinCtx       user;
	Vec             X,F,height;
	DM              daspm;
	
	PetscFunctionBegin;
	
	ierr = pTatin3dCreateContext(&user);CHKERRQ(ierr);
	ierr = pTatin3dSetFromOptions(user);CHKERRQ(ierr);

	/* Register all models */
	ierr = pTatinModelRegisterAll();CHKERRQ(ierr);
	/* Load model, call an initialization routines */
	ierr = pTatinModelLoad(user);CHKERRQ(ierr);
	
	ierr = pTatinModel_Initialize(user->model,user);CHKERRQ(ierr);
	
	/* Generate physics modules */
	ierr = pTatin3d_PhysCompStokesCreate(user);CHKERRQ(ierr);

	/* Pack all physics together */
	/* Here it's simple, we don't need a DM for this, just assign the pack DM to be equal to the stokes DM */
	ierr = PetscObjectReference((PetscObject)user->stokes_ctx->stokes_pack);CHKERRQ(ierr);
	user->pack = user->stokes_ctx->stokes_pack;

	/* fetch some local variables */
	multipys_pack = user->pack;
	dav           = user->stokes_ctx->dav;
	
	ierr = DMCreateGlobalVector(multipys_pack,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);	
	
	ierr = pTatin3dCreateMaterialPoints(user,dav);CHKERRQ(ierr);
	
	/* mesh geometry */
	ierr = pTatinModel_ApplyInitialMeshGeometry(user->model,user);CHKERRQ(ierr);
	
	/* interpolate point coordinates (needed if mesh was modified) */
	//ierr = QuadratureStokesCoordinateSetUp(user->stokes_ctx->Q,dav);CHKERRQ(ierr);
	//for (e=0; e<QUAD_EDGES; e++) {
	//	ierr = SurfaceQuadratureStokesGeometrySetUp(user->stokes_ctx->surfQ[e],dav);CHKERRQ(ierr);
	//}
	/* interpolate material point coordinates (needed if mesh was modified) */
	ierr = MaterialPointCoordinateSetUp(user,dav);CHKERRQ(ierr);
	
	/* material geometry */
	ierr = pTatinModel_ApplyInitialMaterialGeometry(user->model,user);CHKERRQ(ierr);
	
	/* boundary conditions */
	ierr = pTatinModel_ApplyBoundaryCondition(user->model,user);CHKERRQ(ierr);

	{
		Vec Xu,Xp;
		ierr = DMCompositeGetAccess(multipys_pack,X,&Xu,&Xp);CHKERRQ(ierr);
		ierr = BCListInsert(user->stokes_ctx->u_bclist,Xu);CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(multipys_pack,X,&Xu,&Xp);CHKERRQ(ierr);
	}
	
	/* test form function */
	{
		SNES snes;
		
		
		ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
		ierr = SNESSetFunction(snes,F,FormFunction_Stokes,user);CHKERRQ(ierr);  
		ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
		ierr = SNESDestroy(&snes);CHKERRQ(ierr);
	}
	
	/* test data bucket viewer */
	DataBucketView(PetscObjectComm((PetscObject)multipys_pack), user->materialpoint_db,"materialpoint_stokes",DATABUCKET_VIEW_STDOUT);
	DataBucketView(PETSC_COMM_SELF, user->material_constants,"material_constants",DATABUCKET_VIEW_STDOUT);
	

	/* write out the initial condition */
	ierr = pTatinModel_Output(user->model,user,X,"icbc");CHKERRQ(ierr);
	
	/* test generic viewer */
	{
		const int nf = 1;
		const MaterialPointField mp_prop_list[] = { MPField_Std }; 
		ierr = SwarmViewGeneric_ParaView(user->materialpoint_db,nf,mp_prop_list,user->outputpath,"test_MPStd");CHKERRQ(ierr);
	}
	{
		const int nf = 1;
		const MaterialPointField mp_prop_list[] = { MPField_Stokes }; 
		ierr = SwarmViewGeneric_ParaView(user->materialpoint_db,nf,mp_prop_list,user->outputpath,"test_MPStokes");CHKERRQ(ierr);
	}
	
	{
		const int nf = 2;
		const MaterialPointField mp_prop_list[] = { MPField_Std, MPField_Stokes }; 
		ierr = SwarmViewGeneric_ParaView(user->materialpoint_db,nf,mp_prop_list,user->outputpath,"test_MPStd_MPStokes");CHKERRQ(ierr);
	}
	
	
	/* generate a thermal mesh */
	// overlapping form
	{
		DM daq1;
		Vec phi;
		
		ierr = DMDACreateOverlappingQ1FromQ2(dav,1,&daq1);CHKERRQ(ierr);
		//ierr = DMDACreateNestedQ1FromQ2(dav,1,&daq1);CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(daq1,&phi);CHKERRQ(ierr);
		ierr = DMDAViewPetscVTK(daq1,phi,"phi_overlapping_q1.vtk");CHKERRQ(ierr);
		
		ierr = VecDestroy(&phi);CHKERRQ(ierr);
		ierr = DMDestroy(&daq1);CHKERRQ(ierr);
	}
	
	
	/* generate a parallel surface mesh */
	ierr = pTatinSurfaceMeshCreate(dav,&daspm,&height);CHKERRQ(ierr);
	ierr = DMDAViewPetscVTK(daspm,height,"surf_pt1.vtk");CHKERRQ(ierr);

	{
		Vec rand;
		
		VecDuplicate(height,&rand);
		VecSetRandom(rand,NULL);
		VecScale(rand,0.5);
		VecAXPY(height,1.0,rand);
		ierr = VecDestroy(&rand);CHKERRQ(ierr);
	}
	ierr = DMDAViewPetscVTK(daspm,height,"surf_pt2.vtk");CHKERRQ(ierr);
	
	ierr = pTatin_InjectSurfaceMeshOntoMechanicalDomain(daspm,height,dav);CHKERRQ(ierr);
	ierr = pTatinModel_Output(user->model,user,X,"spm");CHKERRQ(ierr);
	
	
	ierr = DMDASetUniformCoordinates(dav, 0.0,4.0, 0.0,1.0, 0.0,6.0);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	ierr = pTatin_InjectMechanicalDomainSurfaceOntoSurfaceMesh(dav,daspm,height);CHKERRQ(ierr);
	ierr = DMDAViewPetscVTK(daspm,height,"surf_pt3.vtk");CHKERRQ(ierr);
	
	
	ierr = DMDestroy(&daspm);CHKERRQ(ierr);
	ierr = VecDestroy(&height);CHKERRQ(ierr);
	
	
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	ierr = VecDestroy(&F);CHKERRQ(ierr);
	ierr = pTatin3dDestroyContext(&user);

	PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
	
	ierr = pTatin3d_material_points_check_ic(argc,argv);CHKERRQ(ierr);
	
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}

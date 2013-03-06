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
 **    Filename:      model_ops_folding2d.c
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
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_models.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"
#include "dmda_update_coords.h"

#include "model_folding2d_ctx.h"

const double gravity = 9.8;

#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_Folding2d"
PetscErrorCode ModelInitialize_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscInt n_int,n;
	PetscBool flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* assign defaults */
	data->max_layers = 100;
	
	data->n_interfaces = 2;
	PetscOptionsGetInt(PETSC_NULL,"-model_folding2d_n_interfaces",&data->n_interfaces,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the number of interfaces including the top and bottom boundaries (-model_folding2d_n_interfaces)");
	}
	
	PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_Lx",&data->Lx,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the x direction (-model_folding2d_Lx)");
	}
	
	/*PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_Ly",&data->Ly,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the y direction (-model_folding2d_Ly)");
	}
	
	PetscOptionsGetInt(PETSC_NULL,"-model_folding2d_Lz",&data->Lz,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide the length along the z direction (-model_folding2d_Lz)");
	}*/

	n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_folding2d_interface_heights",data->interface_heights,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide interface heights relative from the base of the model including the top and bottom boundaries (-model_folding2d_interface_heights)");
	}
	if (n_int != data->n_interfaces) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d interface heights relative from the base of the model including the top and bottom boundaries (-model_folding2d_interface_heights)",data->n_interfaces);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetIntArray(PETSC_NULL,"-model_folding2d_layer_res_j",data->layer_res_j,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer resolution list (-model_folding2d_layer_res_j)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer resolutions (-model_folding2d_layer_res_j)",data->n_interfaces-1);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_folding2d_layer_eta",data->eta,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer viscosity list (-model_folding2d_layer_eta)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer viscosity (-model_folding2d_layer_eta)",data->n_interfaces-1);
	}
	
	n_int = data->max_layers;
	PetscOptionsGetRealArray(PETSC_NULL,"-model_folding2d_layer_rho",data->rho,&n_int,&flg);
	if (!flg) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide layer density list (-model_folding2d_layer_rho)");
	}
	if (n_int != data->n_interfaces-1) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"User must provide %d layer density (-model_folding2d_layer_rho)",data->n_interfaces-1);
	}
	
	/* define the mesh size the y-direction for the global problem */
	c->my = 0;
	for (n=0; n<data->n_interfaces-1; n++) {
		c->my += data->layer_res_j[n];
	}
	
	data->bc_type = 0; /* 0 use vx compression ; 1 use exx compression */
	data->exx             = -1.0e-3;
	data->vx_commpression = 1.0;
	
	/* parse from command line or input file */
	ierr = PetscOptionsGetInt(PETSC_NULL,"-model_folding2d_bc_type",&data->bc_type,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_exx",&data->exx,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,"-model_folding2d_vx",&data->vx_commpression,&flg);CHKERRQ(ierr);
	
	data->Ly = data->interface_heights[data->n_interfaces-1];
	
	PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "BoundaryCondition_Folding2d"
PetscErrorCode BoundaryCondition_Folding2d(DM dav,BCList bclist,pTatinCtx c,ModelFolding2dCtx *data)
{
	PetscReal         exx;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	exx = data->exx;
	
	if (data->bc_type == 0) {
		/* compression east/west in the x-direction (0) [east-west] using constant velocity */
		ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,EAST_FACE,-data->vx_commpression);CHKERRQ(ierr);
		ierr = DirichletBC_ApplyNormalVelocity(bclist,dav,WEST_FACE, data->vx_commpression);CHKERRQ(ierr);
	} else if (data->bc_type == 1) {
		/* compression east/west in the x-direction (0) [east-west] using constant strain rate */
		ierr = DirichletBC_ApplyDirectStrainRate(bclist,dav,exx,0);CHKERRQ(ierr);
	} else {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Unknonwn boundary condition type");
	}
	
	/* free slip south (base) */
	ierr = DirichletBC_FreeSlip(bclist,dav,SOUTH_FACE);CHKERRQ(ierr);
	
	/* free surface north */
	/* do nothing! */
	
	/* free slip front/back to mimic 2d behaviour */
	ierr = DirichletBC_FreeSlip(bclist,dav,FRONT_FACE);CHKERRQ(ierr);
	ierr = DirichletBC_FreeSlip(bclist,dav,BACK_FACE);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_Folding2d"
PetscErrorCode ModelApplyBoundaryCondition_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscReal         exx;
	BCList            bclist;
	DM                dav;
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	
	exx = data->exx;

	bclist = c->stokes_ctx->u_bclist;
	dav    = c->stokes_ctx->dav;
	ierr = BoundaryCondition_Folding2d(dav,bclist,c,data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_Folding2d"
PetscErrorCode ModelApplyBoundaryConditionMG_Folding2d(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx user,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscInt n;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	for (n=0; n<nl; n++) {
		/* Define boundary conditions for each level in the MG hierarchy */
		ierr = BoundaryCondition_Folding2d(dav[n],bclist[n],user,data);CHKERRQ(ierr);
	}	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyMaterialBoundaryCondition_Folding2d"
PetscErrorCode ModelApplyMaterialBoundaryCondition_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Folding2dSetPerturbedInterfaces"
PetscErrorCode Folding2dSetPerturbedInterfaces(DM dav, PetscScalar interface_heights[], PetscInt layer_res_j[], PetscInt n_interfaces)
{
	PetscErrorCode ierr;
	PetscInt i,k,si,sj,sk,nx,ny,nz,M,N,P, interf, jinter;
	PetscScalar *random, amp = 1./4., dy;
	DM cda;
	Vec coord;
	DMDACoor3d ***LA_coord;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	ierr = DMDAGetInfo(dav,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dav,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
	ierr = DMDAGetCoordinateDA(dav,&cda);CHKERRQ(ierr);
	ierr = DMDAGetCoordinates(dav,&coord);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	
	
	/*Perturbes the interface for cylindrical folding*/
	ierr = PetscMalloc(M*sizeof(PetscScalar), &random);CHKERRQ(ierr);
	jinter = 0;
	for(interf = 1; interf < n_interfaces-1; interf++){
		jinter += 2*layer_res_j[interf-1];
		PetscPrintf(PETSC_COMM_WORLD,"jinter = %d (max=%d)\n", jinter,N-1 );

		srand(interf+2);//The seed changes with the interface but we have the same seed for each process.
		for(i = 0; i<M; i++){
			random[i] = 2.0 * rand()/(RAND_MAX+1.0) - 1.0; 
		}

		if ( (jinter>=sj) && (jinter<sj+ny) ) {
			
			dy = 0.5*((interface_heights[interf+1] - interface_heights[interf])/(PetscScalar)(layer_res_j[interf]) + (interface_heights[interf] - interface_heights[interf-1])/(PetscScalar)(layer_res_j[interf-1]) );
			PetscPrintf(PETSC_COMM_SELF," interface %d: using dy computed from avg %1.4e->%1.4e / my=%d :: %1.4e->%1.4e / my=%d \n", interf,
									interface_heights[interf+1],interface_heights[interf],layer_res_j[interf],
									interface_heights[interf],interface_heights[interf-1],layer_res_j[interf-1] );
			//dy = ( (interface_heights[interf] - interface_heights[interf-1]) )/( (PetscScalar)(layer_res_j[interf-1]) );
			//PetscPrintf(PETSC_COMM_WORLD," interface %d: using dy computed %1.4e->%1.4e and my=%d \n", interf,interface_heights[interf],interface_heights[interf-1],layer_res_j[interf-1]);
			
			for(i = si; i<si+nx; i++) {
				//LA_coord[sk][jinter][i].y += amp * dy * random[i];
				for(k = sk; k<sk+nz; k++) {
					//LA_coord[k][jinter][i].y = LA_coord[sk][jinter][i].y;
					LA_coord[k][jinter][i].y += amp * dy * random[i];
				}
			}
			
		}
	}
	
	ierr = PetscFree(random);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,coord,&LA_coord);CHKERRQ(ierr);
	ierr = DMDAUpdateGhostedCoordinates(dav);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPntGetField_global_element_IJKindex"
PetscErrorCode MPntGetField_global_element_IJKindex(DM da, MPntStd *material_point, PetscInt *I, PetscInt *J, PetscInt *K)
{
	PetscInt    li, lj, lk,lmx, lmy, lmz, si, sj, sk, localeid;	
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	MPntStdGetField_local_element_index(material_point,&localeid);
	ierr = DMDAGetCornersElementQ2(da,&si,&sj,&sk,&lmx,&lmy,&lmz);CHKERRQ(ierr);

	si = si/2; 
	sj = sj/2;
	sk = sk/2;
//	lmx -= si;
//	lmy -= sj;
//	lmz -= sk;
	//global/localrank = mx*my*k + mx*j + i;
	lk = (PetscInt)localeid/(lmx*lmy);
	lj = (PetscInt)(localeid - lk*(lmx*lmy))/lmx;
	li = localeid - lk*(lmx*lmy) - lj*lmx;

	if ( (li < 0) || (li>=lmx) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"I computed incorrectly"); }
	if ( (lj < 0) || (lj>=lmy) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"J computed incorrectly"); }
	if ( (lk < 0) || (lk>=lmz) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"K computed incorrectly"); }
	//printf("li,lj,lk %d %d %d \n", li,lj,lk );
	
	*K = lk + sk;
	*J = lj + sj;
	*I = li + si;
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryMaterialPoints_Folding2d"
PetscErrorCode InitialMaterialGeometryMaterialPoints_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	PetscErrorCode ierr;
			
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
			
			
	/* define properties on material points */
	db = c->materialpoint_db;
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
			
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));
			
			
	DataBucketGetSizes(db,&n_mp_points,0,0);
			
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		//double      *position;
		PetscReal      eta,rho;
		PetscInt    phase;
		PetscInt    layer, jmaxlayer, jminlayer;
		PetscInt    I, J, K;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		/* Access using the getter function provided for you (recommeneded for beginner user) */
		//MPntStdGetField_global_coord(material_point,&position)

    MPntGetField_global_element_IJKindex(c->stokes_ctx->dav,material_point, &I, &J, &K);
		phase = -1;
		eta =  0.0;
		rho = 0.0;
		jmaxlayer = jminlayer = 0;
		layer = 0;
		// gets the global element index (i,j,k)
		//....
		
		//Set the properties
		while( (phase == -1) && (layer < data->n_interfaces-1) ){
			jmaxlayer += data->layer_res_j[layer];
			
			if( (J<jmaxlayer) && (J>=jminlayer) ){
				phase = layer + 1;
				eta = data->eta[layer];
				rho = data->rho[layer];

				rho = - rho * gravity;
			}
			jminlayer += data->layer_res_j[layer];
			layer++;
		}

		/* user the setters provided for you */
		MPntStdSetField_phase_index(material_point,phase);
		MPntPStokesSetField_eta_effective(mpprop_stokes,eta);
		MPntPStokesSetField_density(mpprop_stokes,rho);
	}
			
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
			
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitialMaterialGeometryQuadraturePoints_Folding2d"
PetscErrorCode InitialMaterialGeometryQuadraturePoints_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	PhysCompStokes         user;
	QPntVolCoefStokes      *all_gausspoints,*cell_gausspoints;
	PetscInt               nqp,qp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	
	/* define properties on material points */
	db = c->materialpoint_db;
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&PField_stokes);
	DataFieldGetAccess(PField_stokes);
	DataFieldVerifyAccess(PField_stokes,sizeof(MPntPStokes));
	
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	

	/* get the quadrature points */
	user = c->stokes_ctx;
	ierr = VolumeQuadratureGetAllCellData_Stokes(user->volQ,&all_gausspoints);CHKERRQ(ierr);
	nqp = user->volQ->npoints;
	
	for (p=0; p<n_mp_points; p++) {
		MPntStd     *material_point;
		MPntPStokes *mpprop_stokes;
		//double      *position;
		PetscReal      eta,rho;
		PetscInt    phase;
		PetscInt    layer, jmaxlayer, jminlayer, localeid_p;
		PetscInt    I, J, K;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		DataFieldAccessPoint(PField_stokes,p,(void**)&mpprop_stokes);
		
    MPntGetField_global_element_IJKindex(c->stokes_ctx->dav,material_point, &I, &J, &K);

		//Set the properties
		phase = -1;
		eta =  0.0;
		rho = 0.0;
		jmaxlayer = jminlayer = 0;
		layer = 0;
		while( (phase == -1) && (layer < data->n_interfaces-1) ){
			jmaxlayer += data->layer_res_j[layer];
			
			if( (J<jmaxlayer) && (J>=jminlayer) ){
				phase = layer + 1;
				eta = data->eta[layer];
				rho = data->rho[layer];
			}
			jminlayer += data->layer_res_j[layer];
			layer++;
		}

		
		MPntStdGetField_local_element_index(material_point,&localeid_p);
		ierr = VolumeQuadratureGetCellData_Stokes(user->volQ,all_gausspoints,localeid_p,&cell_gausspoints);CHKERRQ(ierr);
		
		for (qp=0; qp<nqp; qp++) {
			cell_gausspoints[qp].eta  = eta;
			cell_gausspoints[qp].rho  = rho;

			cell_gausspoints[qp].Fu[0] = 0.0;
			cell_gausspoints[qp].Fu[1] = -rho * gravity;
			cell_gausspoints[qp].Fu[2] = 0.0;

			cell_gausspoints[qp].Fp = 0.0;
		}		
		
	}
	
	DataFieldRestoreAccess(PField_std);
	DataFieldRestoreAccess(PField_stokes);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_Folding2d"
PetscErrorCode ModelApplyInitialMaterialGeometry_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	int                    p,n_mp_points;
	DataBucket             db;
	DataField              PField_std,PField_stokes;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	ierr = InitialMaterialGeometryMaterialPoints_Folding2d(c,ctx);CHKERRQ(ierr);
	ierr = InitialMaterialGeometryQuadraturePoints_Folding2d(c,ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

		
		
		
#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_Folding2d"
PetscErrorCode ModelApplyInitialMeshGeometry_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscReal Lx,Ly,dx,dy,dz,Lz;
	PetscInt mx,my,mz, itf;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);

	/* step 1 - create structured grid */
	Lx = data->Lx;
	Ly = data->Ly;
	
	/*
	 The length of the model in z-direction is determined by the grid spacing in x and y.
	 We do this so that the elements do not have a large aspect ratio. This would occur
	 if we hard coded Lz to a constant number which is independnet of the grid resolution in x and y.
	 We choose Lz to be mz * min(dx,dy).
	 */
	
	
	mx = c->mx; 
	my = c->my; 
	mz = c->mz; 
	
	my=0;
	for(itf = 0; itf<data->n_interfaces -1; ++itf){
		my += data->layer_res_j[itf];
	}

	
	dx = Lx / ((PetscReal)mx);
	dy = Ly / ((PetscReal)my);
	dz = dx;
	if (dz < dy) {
		dz = dy;
	}
	Lz = mz * dz;
	
	ierr = DMDASetUniformCoordinates(c->stokes_ctx->dav, 0.0,Lx, data->interface_heights[0],Ly, 0.0,Lz);CHKERRQ(ierr);

	
	/* step 2 - define two interfaces and perturb coords along the interface */
	ierr = Folding2dSetPerturbedInterfaces(c->stokes_ctx->dav, data->interface_heights, data->layer_res_j, data->n_interfaces);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ModelApplyUpdateMeshGeometry_Folding2d"
PetscErrorCode ModelApplyUpdateMeshGeometry_Folding2d(pTatinCtx c,Vec X,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	/*Remeshing*/
	/*remeshfunction(DM dav, PetscScalar interface_heights[], PetscInt layer_res_j[], PetscInt n_interfaces){
	jinter = sj;
	 for(interf = 0; interf < n_interfaces-1; ++interf){
	 dy = (interface_heights[interf+1] - interface_heights[interf])/(PetscScalar)(2*layer_res_j[interf]);
	 for(i = si; i<si+nx; i++) {
	 for(k = sk; k<sk+nz; k++) {
	 for(j = jinter + 2*layer_res_j[interf]-1; j >jinter; --j){
	 LA_coord[k][j][i].y = LA_coord[k][j+1][i].y - dy;
	 }	
	 }
	 }
	 jinter += 2*layer_res_j[interf];
	 }
	 
	 }*/	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialCondition_Foldinf2D"
PetscErrorCode ModelInitialCondition_Foldinf2D(pTatinCtx c,Vec X,void *ctx)
{
    /*
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	DM stokes_pack,dau,dap;
	Vec velocity,pressure;
	PetscReal rho0;
	DMDAVecTraverse3d_HydrostaticPressureCalcCtx HPctx;
	DMDAVecTraverse3d_InterpCtx IntpCtx;*/
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	/*
	
	stokes_pack = c->stokes_ctx->stokes_pack;
	
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    
	ierr = VecZeroEntries(velocity);CHKERRQ(ierr);
	ierr = VecZeroEntries(pressure);CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(stokes_pack,X,&velocity,&pressure);CHKERRQ(ierr);
    */
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_Folding2d"
PetscErrorCode ModelOutput_Folding2d(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_Folding2d"
PetscErrorCode ModelDestroy_Folding2d(pTatinCtx c,void *ctx)
{
	ModelFolding2dCtx *data = (ModelFolding2dCtx*)ctx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", __FUNCT__);
	
	/* Free contents of structure */
	
	/* Free structure */
	ierr = PetscFree(data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_Folding2d"
PetscErrorCode pTatinModelRegister_Folding2d(void)
{
	ModelFolding2dCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* Allocate memory for the data structure for this model */
	ierr = PetscMalloc(sizeof(ModelFolding2dCtx),&data);CHKERRQ(ierr);
	ierr = PetscMemzero(data,sizeof(ModelFolding2dCtx));CHKERRQ(ierr);
	
	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

	/* Set name, model select via -ptatin_model NAME */
	ierr = pTatinModelSetName(m,"folding2d");CHKERRQ(ierr);

	/* Set model data */
	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);
	
	/* Set function pointers */
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Folding2d);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Folding2d);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

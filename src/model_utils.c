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
 **    Filename:      model_utils.c
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

#include "ptatin3d.h"
#include "ptatin3d_defs.h"
#include "private/ptatin_impl.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "mesh_quality_metrics.h"
#include "mesh_update.h"
#include "mesh_deformation.h"

#include "dmda_bcs.h"
#include "dmda_iterator.h"
#include "dmda_redundant.h"
#include "dmda_remesh.h"
#include "dmda_update_coords.h"
#include "dmda_duplicate.h"
#include "dmdae.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1.h"
#include "element_type_Q2.h"
#include "element_utils_q2.h"
#include "element_utils_q1.h"

#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPStokesPl_def.h"
#include "MPntPEnergy_def.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"

#include "model_utils.h"


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
#define __FUNCT__ "DMDAConvertLocalElementIndex2GlobalIJK"
PetscErrorCode DMDAConvertLocalElementIndex2GlobalIJK(DM da,PetscInt localeid,PetscInt *I,PetscInt *J,PetscInt *K)
{
	PetscInt       li,lj,lk,lmx,lmy,lmz,si,sj,sk;	
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = DMDAGetCornersElementQ2(da,&si,&sj,&sk,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	
	si = si/2; 
	sj = sj/2;
	sk = sk/2;

	//global/localrank = mx*my*k + mx*j + i;
	lk = (PetscInt)localeid/(lmx*lmy);
	lj = (PetscInt)(localeid - lk*(lmx*lmy))/lmx;
	li = localeid - lk*(lmx*lmy) - lj*lmx;
	
	if ( (li < 0) || (li >= lmx) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"I computed incorrectly"); }
	if ( (lj < 0) || (lj >= lmy) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"J computed incorrectly"); }
	if ( (lk < 0) || (lk >= lmz) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"K computed incorrectly"); }
	*K = lk + sk;
	*J = lj + sj;
	*I = li + si;

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAConvertLocalNodeIndex2GlobalIJK"
PetscErrorCode DMDAConvertLocalNodeIndex2GlobalIJK(DM da,PetscInt localnid,PetscInt *I,PetscInt *J,PetscInt *K)
{
	PetscInt       li,lj,lk,lnx,lny,lnz,si,sj,sk;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = DMDAGetCorners(da,&si,&sj,&sk,&lnx,&lny,&lnz);CHKERRQ(ierr);
	
	//global/localrank = mx*my*k + mx*j + i;
	lk = (PetscInt)localnid/(lnx*lny);
	lj = (PetscInt)(localnid - lk*(lnx*lny))/lnx;
	li = localnid - lk*(lnx*lny) - lj*lnx;
	
	if ( (li < 0) || (li >= lnx) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"I computed incorrectly"); }
	if ( (lj < 0) || (lj >= lny) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"J computed incorrectly"); }
	if ( (lk < 0) || (lk >= lnz) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"K computed incorrectly"); }
	*K = lk + sk;
	*J = lj + sj;
	*I = li + si;
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAConvertLocalGhostNodeIndex2GlobalIJK"
PetscErrorCode DMDAConvertLocalGhostNodeIndex2GlobalIJK(DM da,PetscInt localnid,PetscInt *I,PetscInt *J,PetscInt *K)
{
	PetscInt       li,lj,lk,lnx,lny,lnz,si,sj,sk;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&lnx,&lny,&lnz);CHKERRQ(ierr);
	
	//global/localrank = mx*my*k + mx*j + i;
	lk = (PetscInt)localnid/(lnx*lny);
	lj = (PetscInt)(localnid - lk*(lnx*lny))/lnx;
	li = localnid - lk*(lnx*lny) - lj*lnx;
	
	if ( (li < 0) || (li >= lnx) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"I computed incorrectly"); }
	if ( (lj < 0) || (lj >= lny) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"J computed incorrectly"); }
	if ( (lk < 0) || (lk >= lnz) ) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"K computed incorrectly"); }
	*K = lk + sk;
	*J = lj + sj;
	*I = li + si;
    
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinModelGetOptionReal"
PetscErrorCode pTatinModelGetOptionReal(const char option[],PetscReal *val,
																				const char error[],
																				const char default_opt[],
																				PetscBool essential)
{
	PetscBool flg;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	ierr = PetscOptionsGetReal(NULL,option,val,&flg);CHKERRQ(ierr);
	if (essential) {
		if (!flg) {
			if (!default_opt) {
				SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"ModelOptionMissing(%s): \n\t\t%s ",option,error);
			} else {
				SETERRQ3(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"ModelOptionMissing(%s): \n\t\t%s : Suggested default values %s ",option,error,default_opt);
			}
		}
	}
	PetscFunctionReturn(0);
}

/* Absolute value for double in C - it's called fabs(double a) or PetscAbsReal(PetscRea a) */
#undef __FUNCT__
#define __FUNCT__ "absolute"
PetscReal absolute(PetscReal a)
{   
	if(a < 0) {
        return -1.0 * a;
	} else {
		return a;
	}
}

/*
 Remove the linear trend
 Assume a regular spacing and chose an appropriate coordinate system such
 sum x_i = 0
*/
#undef __FUNCT__
#define __FUNCT__ "detrend"
PetscErrorCode detrend(PetscReal array[],PetscInt n)
{
	PetscReal x,y,a,b;
	PetscReal sy = 0.0,sxy = 0.0,sxx = 0.0;
	PetscInt  i;
	
	for (i=0, x=(-n/2.0+0.5); i<n; i++, x+=1.0) {
		y = array[i];
		sy += y;
		sxy += x * y;
		sxx += x * x;
	}
	a = sxy/sxx;
	b = sy/(PetscReal)n;
	
	for (i=0, x=(-n/2.0+0.5); i<n; i++, x+=1.0) {
		array[i] -= (a*x+b);
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "rednoise"
PetscErrorCode rednoise(PetscReal rnoise[],PetscInt n,PetscInt seed)
{
	PetscInt       i;
	PetscReal      maxi;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	srand(seed);
	rnoise[0] = 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
	
	for (i=1; i<n; i++) {
		rnoise[i] = rnoise[i-1] + 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
	}
	
	detrend(rnoise,n);
	
	maxi = 0.0;
	for (i=1; i<n; i++) {
		
		if (abs(rnoise[i]) > maxi) {
			maxi = abs(rnoise[i]);
		}
	}
	for (i=1; i<n; i++) {
		rnoise[i] /= maxi;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAVecTraverse_InitialThermalField3D"
PetscBool DMDAVecTraverse_InitialThermalField3D(PetscScalar pos[],PetscScalar *val,void *ctx)
{
    DMDA_thermalfield_init_params *thermalparams;
    PetscInt i,klay,nlt;
    PetscReal y,yldep,dtemp;

		thermalparams = (DMDA_thermalfield_init_params*)ctx;
		y   = pos[1] * thermalparams->lscale;
    nlt = thermalparams->nlayers;
//  Which layer contains the current node?
    klay = 0;

    for (i=0; i <= nlt-1; i++) {
        if (y <= thermalparams->ytop[i] && y >= thermalparams->ytop[i+1]) {
            klay = i;
            break;
        }
    }
    yldep = thermalparams->ytop[klay] - y;
    dtemp = thermalparams->hp[klay] * yldep * (thermalparams->thick[klay]-yldep/2.0e0) + thermalparams->qbase[klay]*yldep;
    dtemp = dtemp / thermalparams->cond[klay];
    *val  = thermalparams->ttop[klay] + dtemp;
/*    if (pos[0]==0&&pos[2]==0) {
        printf("y=%1.4e, T=%1.4e \n",y,*val);
        printf("ytop0=%1.4e\n",thermalparams->ytop[0]);
        printf("ytop1=%1.4e\n",thermalparams->ytop[1]);
        printf("ytop1=%1.4e\n",thermalparams->ytop[2]);
        printf("klay=%d\n",klay);
    }
 */

    return PETSC_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "DMDAComputeMeshVolume"
PetscErrorCode DMDAComputeMeshVolume(DM dm,PetscReal *value)
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e,p;
	const PetscInt  *el_nidx;
	const PetscInt  *gidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	PetscInt        ngp;
	PetscReal       WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal       _value,detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* setup quadrature */
	ngp = 27;
	P3D_prepare_elementQ2(ngp,WEIGHT,XI,NI,GNI);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dm,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = DMDAGetGlobalIndices(dm,0,&gidx);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);

	_value = 0.0;
	for (e=0; e<nel; e++) {
		PetscReal el_vol;
		
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		P3D_evaluate_geometry_elementQ2(ngp,el_coords,GNI, detJ,dNudx,dNudy,dNudz);
		
		el_vol = 0.0;
		for (p=0; p<ngp; p++) {
			el_vol = el_vol + 1.0 * WEIGHT[p] * detJ[p];
		}
		_value += el_vol;
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&_value,value,1,MPIU_REAL,MPI_SUM,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin3d_DefineVelocityMeshQuasi2D"
PetscErrorCode pTatin3d_DefineVelocityMeshQuasi2D(pTatinCtx c)
{
	c->mz = 1;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin3d_DefineVelocityMeshGeometryQuasi2D"
PetscErrorCode pTatin3d_DefineVelocityMeshGeometryQuasi2D(pTatinCtx c)
{
	PhysCompStokes    stokes;
	DM                stokes_pack,dav,dap;
	PetscReal Lz,min_dl[3],max_dl[3];
	PetscBool geom_max,flg;
	PetscErrorCode ierr;

	ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	/* determine min/max dx,dy,dz for mesh */
	ierr = DMDAComputeQ2ElementBoundingBox(dav,min_dl,max_dl);CHKERRQ(ierr);
	
	geom_max = PETSC_FALSE;
	Lz = 1.0e32;
	Lz = PetscMin(Lz,min_dl[0]);
	Lz = PetscMin(Lz,min_dl[1]);

	PetscOptionsGetBool(NULL,"-ptatin_geometry_quasi_2d_max",&geom_max,NULL);
	if (geom_max) {
		Lz = 1.0e-32;
		Lz = PetscMax(Lz,max_dl[0]);
		Lz = PetscMax(Lz,max_dl[1]);
	}

	if (geom_max) {
		PetscPrintf(PETSC_COMM_WORLD,"[[pTatin3d_DefineVelocityMeshGeometryQuasi2D]] Using Lz = %1.4e from max(dx,dy) \n",Lz );
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"[[pTatin3d_DefineVelocityMeshGeometryQuasi2D]] Using Lz = %1.4e from min(dx,dy) \n",Lz );
	}
	
	ierr = DMDASetUniformCoordinates1D(dav,2,0.0,Lz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAComputeQ2ElementBoundingBox"
PetscErrorCode DMDAComputeQ2ElementBoundingBox(DM dm,PetscReal gmin[],PetscReal gmax[])
{
	DM              cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e;
	const PetscInt  *el_nidx;
	const PetscInt  *gidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	PetscReal       dx,dy,dz,dl_min[3],dl_max[3];
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dm,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dm,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetGlobalIndices(dm,0,&gidx);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dm,&nel,&nen,&el_nidx);CHKERRQ(ierr);
	
	dl_min[0] = dl_min[1] = dl_min[2] = 1.0e32;
	dl_max[0] = dl_max[1] = dl_max[2] = 1.0e-32;
	
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		
		dx = fabs( el_coords[3*Q2_FACE_NODE_EAST +0] - el_coords[3*Q2_FACE_NODE_WEST +0]  );
		dy = fabs( el_coords[3*Q2_FACE_NODE_NORTH+1] - el_coords[3*Q2_FACE_NODE_SOUTH+1] );
		dz = fabs( el_coords[3*Q2_FACE_NODE_FRONT+2] - el_coords[3*Q2_FACE_NODE_BACK +2]  );
		
		if (dx < dl_min[0]) { dl_min[0] = dx; }
		if (dy < dl_min[1]) { dl_min[1] = dy; }
		if (dz < dl_min[2]) { dl_min[2] = dz; }
		
		if (dx > dl_max[0]) { dl_max[0] = dx; }
		if (dy > dl_max[1]) { dl_max[1] = dy; }
		if (dz > dl_max[2]) { dl_max[2] = dz; }
		
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(dl_min,gmin,3,MPIU_REAL,MPI_MIN,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
	ierr = MPI_Allreduce(dl_max,gmax,3,MPIU_REAL,MPI_MAX,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
 This should only ever be used for debugging.
 We scatter the vector created from a DMDA into the natural i+j*nx+k*nx*ny ordering,
 then we scatter this i,j,k ordered vector onto rank 0 and write the contents out in ascii.
*/
#undef __FUNCT__
#define __FUNCT__ "DMDAFieldViewAscii"
PetscErrorCode DMDAFieldViewAscii(DM dm,Vec field,const char filename[])
{
	PetscErrorCode ierr;
	Vec natural_field,natural_field_red;
	VecScatter ctx;
	FILE *fp;
	PetscInt M,N,P,dofs;
	const char *oname = NULL;
	MPI_Comm comm;
	PetscMPIInt rank;
	PetscInt i,n;
	PetscScalar *LA_field;
	

	PetscFunctionBegin;
	
	ierr = DMDACreateNaturalVector(dm,&natural_field);CHKERRQ(ierr);
	ierr = DMDAGlobalToNaturalBegin(dm,field,INSERT_VALUES,natural_field);CHKERRQ(ierr);
	ierr = DMDAGlobalToNaturalEnd(dm,field,INSERT_VALUES,natural_field);CHKERRQ(ierr);

	ierr = VecScatterCreateToZero(natural_field,&ctx,&natural_field_red);CHKERRQ(ierr);
	ierr = VecScatterBegin(ctx,natural_field,natural_field_red,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx,natural_field,natural_field_red,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
	
	ierr = VecDestroy(&natural_field);CHKERRQ(ierr);
	
	/*
	 # DMDAFieldViewAscii: 
	 # DMDA Vec (name) 
	 # M N P x y z
	 # dofs x
	*/
	ierr = DMDAGetInfo(dm,0,&M,&N,&P, 0,0,0, &dofs,0, 0,0,0, 0);CHKERRQ(ierr);
	ierr = PetscObjectGetName((PetscObject)field,&oname);CHKERRQ(ierr);
	ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
	
	if (rank == 0) {
		
		fp = fopen(filename,"w");
		if (fp == NULL) { SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Unable to open file %s on rank 0",filename); }

		fprintf(fp,"# DMDAFieldViewAscii\n");
		if (oname) {
			fprintf(fp,"# DMDA Vec %s\n",oname);
		} else {
			fprintf(fp,"# DMDA Vec\n");
		}
		fprintf(fp,"# M N P %d %d %d\n",M,N,P);
		fprintf(fp,"# dofs %d\n",dofs);
		ierr = VecGetSize(natural_field_red,&n);CHKERRQ(ierr);
		ierr = VecGetArray(natural_field_red,&LA_field);CHKERRQ(ierr);
		for (i=0; i<n; i++) {
			fprintf(fp,"%1.6e\n",LA_field[i]);
		}
		ierr = VecRestoreArray(natural_field_red,&LA_field);CHKERRQ(ierr);
		
		fclose(fp);
	}
	ierr = VecDestroy(&natural_field_red);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelUtilsComputeAiryIsostaticHeights_SEQ"
PetscErrorCode ModelUtilsComputeAiryIsostaticHeights_SEQ(PhysCompStokes stokes)
{
	DM              stokes_pack,dav,dap,cda;
	Vec             gcoords;
	PetscReal       *LA_gcoords;
	PetscInt        nel,nen,e,q;
	const PetscInt  *el_nidx;
	const PetscInt  *gidx;
	PetscReal       el_coords[3*Q2_NODES_PER_EL_3D];
	PetscInt        nqp;
	PetscReal       WEIGHT[NQP],XI[NQP][3],NI[NQP][NPE],GNI[NQP][3][NPE];
	PetscReal       detJ[NQP],dNudx[NQP][NPE],dNudy[NQP][NPE],dNudz[NQP][NPE];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscReal         *vol_col,*rho_col;
	PetscInt          idx,I,J,K,MX,MY,MZ;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
	/* setup quadrature */
	nqp = 27;
	P3D_prepare_elementQ2(nqp,WEIGHT,XI,NI,GNI);
	
	/* get dav */
	stokes_pack = stokes->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dav,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = DMDAGetGlobalIndices(dav,0,&gidx);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen,&el_nidx);CHKERRQ(ierr);

	/* quadrature for all cells */
	ierr = VolumeQuadratureGetAllCellData_Stokes(stokes->volQ,&all_gausspoints);CHKERRQ(ierr);
	
	ierr = DMDAGetSizeElementQ2(dav,&MX,&MY,&MZ);CHKERRQ(ierr);
	
	PetscMalloc(sizeof(PetscReal)*MX*MZ,&vol_col);  PetscMemzero(vol_col,sizeof(PetscReal)*MX*MZ);
	PetscMalloc(sizeof(PetscReal)*MX*MZ,&rho_col);  PetscMemzero(rho_col,sizeof(PetscReal)*MX*MZ);
	
	for (e=0; e<nel; e++) {
		PetscReal el_rho,el_vol;
		
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&el_nidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		P3D_evaluate_geometry_elementQ2(nqp,el_coords,GNI, detJ,dNudx,dNudy,dNudz);

		ierr = VolumeQuadratureGetCellData_Stokes(stokes->volQ,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		el_rho = 0.0;
		el_vol = 0.0;
		for (q=0; q<nqp; q++) {
			el_rho = el_rho + cell_gausspoints[q].rho * WEIGHT[q] * detJ[q];
			el_vol = el_vol + 1.0 * WEIGHT[q] * detJ[q];
		}
		
		ierr = DMDAConvertLocalElementIndex2GlobalIJK(dav,e,&I,&J,&K);CHKERRQ(ierr);
		idx = I + K*MX;
		rho_col[idx] = rho_col[idx] + el_rho;
		vol_col[idx] = vol_col[idx] + el_vol;
		
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	for (e=0; e<MX*MZ; e++) {
		rho_col[e] = rho_col[e] / vol_col[e];
	}
	
	PetscFree(vol_col);
	PetscFree(rho_col);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelUtilsComputeAiryIsostaticHeights"
PetscErrorCode ModelUtilsComputeAiryIsostaticHeights(PhysCompStokes stokes)
{
	PetscErrorCode ierr;

	
	ierr = ModelUtilsComputeAiryIsostaticHeights_SEQ(stokes);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPntStdComputeBoundingBox"
PetscErrorCode MPntStdComputeBoundingBox(DataBucket materialpoint_db,PetscReal gmin[],PetscReal gmax[])
{
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	PetscReal        min[3],max[3];
	double           *pos_p;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	
	/* initialize */
	min[0] =  1.0e32;
	min[1] =  1.0e32;
	min[2] =  1.0e32;
	max[0] = -1.0e32;
	max[1] = -1.0e32;
	max[2] = -1.0e32;
	
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		ierr = MaterialPointGet_global_coord(mpX,p,&pos_p);CHKERRQ(ierr);

		min[0] = PetscMin(min[0],pos_p[0]);
		min[1] = PetscMin(min[1],pos_p[1]);
		min[2] = PetscMin(min[2],pos_p[2]);

		max[0] = PetscMax(max[0],pos_p[0]);
		max[1] = PetscMax(max[1],pos_p[1]);
		max[2] = PetscMax(max[2],pos_p[2]);
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(min,gmin,3,MPIU_REAL,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(max,gmax,3,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPntStdComputeBoundingBoxInRange"
PetscErrorCode MPntStdComputeBoundingBoxInRange(DataBucket materialpoint_db,PetscReal rmin[],PetscReal rmax[],PetscReal gmin[],PetscReal gmax[])
{
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	PetscReal        min[3],max[3];
	double           *pos_p;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	
	/* initialize */
	min[0] =  1.0e32;
	min[1] =  1.0e32;
	min[2] =  1.0e32;
	max[0] = -1.0e32;
	max[1] = -1.0e32;
	max[2] = -1.0e32;
	
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		PetscReal p_i, range_min, range_max;
		PetscInt  idx;
		
		ierr = MaterialPointGet_global_coord(mpX,p,&pos_p);CHKERRQ(ierr);
		
		idx = 0;
		range_min = -1.0e32; if (rmin) { range_min = rmin[idx]; }
		range_max =  1.0e32; if (rmax) { range_max = rmax[idx]; }
		p_i = pos_p[idx];
		if ( (p_i >= range_min) && (p_i <= range_max) ) {
			min[idx] = PetscMin(min[idx],p_i);
			max[idx] = PetscMax(max[idx],p_i);
		}

		idx = 1;
		range_min = -1.0e32; if (rmin) { range_min = rmin[idx]; }
		range_max =  1.0e32; if (rmax) { range_max = rmax[idx]; }
		p_i = pos_p[idx];
		if ( (p_i >= range_min) && (p_i <= range_max) ) {
			min[idx] = PetscMin(min[idx],p_i);
			max[idx] = PetscMax(max[idx],p_i);
		}
		
		idx = 2;
		range_min = -1.0e32; if (rmin) { range_min = rmin[idx]; }
		range_max =  1.0e32; if (rmax) { range_max = rmax[idx]; }
		p_i = pos_p[idx];
		if ( (p_i >= range_min) && (p_i <= range_max) ) {
			min[idx] = PetscMin(min[idx],p_i);
			max[idx] = PetscMax(max[idx],p_i);
		}
		
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(min,gmin,3,MPIU_REAL,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(max,gmax,3,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPntStdComputeBoundingBoxInRangeInRegion"
PetscErrorCode MPntStdComputeBoundingBoxInRangeInRegion(DataBucket materialpoint_db,PetscReal rmin[],PetscReal rmax[],PetscInt region_idx,PetscReal gmin[],PetscReal gmax[])
{
	MPAccess         mpX;
	PetscInt         p,n_mpoints;
	PetscReal        min[3],max[3];
	double           *pos_p;
	int              region_p;
	PetscInt         found_region = 0,found_region_g;
	PetscErrorCode   ierr;
	
	PetscFunctionBegin;
	
	/* initialize */
	min[0] =  1.0e32;
	min[1] =  1.0e32;
	min[2] =  1.0e32;
	max[0] = -1.0e32;
	max[1] = -1.0e32;
	max[2] = -1.0e32;
	
	DataBucketGetSizes(materialpoint_db,&n_mpoints,0,0);
	ierr = MaterialPointGetAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	for (p=0; p<n_mpoints; p++) {
		PetscReal p_i, range_min, range_max;
		PetscInt  idx;
		
		ierr = MaterialPointGet_global_coord(mpX,p,&pos_p);CHKERRQ(ierr);
		ierr = MaterialPointGet_phase_index(mpX,p,&region_p);CHKERRQ(ierr);
		
		if ( (region_idx == -1) && (region_p != region_idx) ) { continue; }

		idx = 0;
		range_min = -1.0e32; if (rmin) { range_min = rmin[idx]; }
		range_max =  1.0e32; if (rmax) { range_max = rmax[idx]; }
		p_i = pos_p[idx];
		if ( (p_i >= range_min) && (p_i <= range_max) ) {
			min[idx] = PetscMin(min[idx],p_i);
			max[idx] = PetscMax(max[idx],p_i);
		}
		
		idx = 1;
		range_min = -1.0e32; if (rmin) { range_min = rmin[idx]; }
		range_max =  1.0e32; if (rmax) { range_max = rmax[idx]; }
		p_i = pos_p[idx];
		if ( (p_i >= range_min) && (p_i <= range_max) ) {
			min[idx] = PetscMin(min[idx],p_i);
			max[idx] = PetscMax(max[idx],p_i);
		}
		
		idx = 2;
		range_min = -1.0e32; if (rmin) { range_min = rmin[idx]; }
		range_max =  1.0e32; if (rmax) { range_max = rmax[idx]; }
		p_i = pos_p[idx];
		if ( (p_i >= range_min) && (p_i <= range_max) ) {
			min[idx] = PetscMin(min[idx],p_i);
			max[idx] = PetscMax(max[idx],p_i);
		}
		
		found_region = 1;
	}
	ierr = MaterialPointRestoreAccess(materialpoint_db,&mpX);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(min,gmin,3,MPIU_REAL,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allreduce(max,gmax,3,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	ierr = MPI_Allreduce(&found_region,&found_region_g,1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	/* no particles of the desired region were found on any processors, set min/max accordingly */
	if (found_region_g == 0) {
		gmin[0] = gmin[1] = gmin[2] = NAN;
		gmax[0] = gmax[1] = gmax[2] = NAN;
	}
	
	PetscFunctionReturn(0);
}

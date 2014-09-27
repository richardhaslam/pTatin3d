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
 **    Filename:      mms_def.c
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


#include <petsc.h>

#include <ptatin3d.h>
#include <ptatin3d_defs.h>
#include <private/ptatin_impl.h>
#include <ptatin_models.h>
#include <model_utils.h>
#include <ptatin3d_stokes.h>
#include <dmda_bcs.h>
#include <element_utils_q2.h>
#include <dmda_element_q2p1.h>

#include "MMSStokes2d_ep.c"

struct _p_MMSCtx {
    PetscInt  mms_idx;
    PetscBool geometry_quasi_2d;
};

typedef struct _p_MMSCtx *MMSCtx;


#undef __FUNCT__
#define __FUNCT__ "ModelInitialize_MMS"
PetscErrorCode ModelInitialize_MMS(pTatinCtx c,void *ctx)
{
    PetscErrorCode ierr;
	MMSCtx         user = (MMSCtx)ctx;
    
    user->mms_idx = 1;
    user->geometry_quasi_2d = PETSC_TRUE;
    
    
    //if (user->geometry_quasi_2d) {
    //    ierr = pTatin3d_DefineVelocityMeshQuasi2D(c);CHKERRQ(ierr);
    //}
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelInitialCondition_MMS"
PetscErrorCode ModelInitialCondition_MMS(pTatinCtx c,Vec X,void *ctx)
{
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMeshGeometry_MMS"
PetscErrorCode ModelApplyInitialMeshGeometry_MMS(pTatinCtx c,void *ctx)
{
    PetscErrorCode ierr;
	MMSCtx         user = (MMSCtx)ctx;
    PhysCompStokes stokes;
    DM             dav;
    PetscReal      Lx,Ly,Lz,gravity[3];
    
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMs(stokes,&dav,NULL);CHKERRQ(ierr);

    Lx = Ly = Lz = 1.0;
	ierr = DMDASetUniformCoordinates(dav,0.0,Lx, 0.0,Ly, 0.0,Lz);CHKERRQ(ierr);
    //if (user->geometry_quasi_2d) {
    //    ierr = pTatin3d_DefineVelocityMeshGeometryQuasi2D(c);CHKERRQ(ierr);
    //}
    
    gravity[0] = 0.0;
    gravity[1] = 0.0;
    gravity[2] = 0.0;
    ierr = PhysCompStokesSetGravityUnitVector(stokes,gravity);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MMSStokes1_ModelApplyInitialMaterialGeometry"
PetscErrorCode MMSStokes1_ModelApplyInitialMaterialGeometry(pTatinCtx c,MMSCtx user)
{
    PetscErrorCode ierr;
    PetscInt e,q,k,nqp,nel,nen_u;
    QPntVolCoefStokes *quadpoints,*cell_quadpoints;
    Quadrature volQ;
    PhysCompStokes stokes;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][Q2_NODES_PER_EL_3D],GNI[NQP][3][Q2_NODES_PER_EL_3D];
	const PetscInt *elnidx_u;
	DM dav,cda;
	Vec gcoords;
	PetscReal *LA_gcoords;

    
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);

    ierr = PhysCompStokesGetDMs(stokes,&dav,NULL);CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dav,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);

    /* quadrature setup */
    ierr = PhysCompStokesGetVolumeQuadrature(stokes,&volQ);CHKERRQ(ierr);
    ierr = VolumeQuadratureGetInfo(volQ,&nqp,&nel);CHKERRQ(ierr);
    P3D_prepare_elementQ2(nqp,WEIGHT,XI,NI,GNI);
    ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&quadpoints);CHKERRQ(ierr);

	for (e=0; e<nel; e++) {
        PetscReal el_coord[3*Q2_NODES_PER_EL_3D];
        
        /* get element coordinates */
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coord,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);

		ierr = VolumeQuadratureGetCellData_Stokes(volQ,quadpoints,e,&cell_quadpoints);CHKERRQ(ierr);
		for (q=0; q<nqp; q ++) {
            PetscReal pos_q[3],Fu_q[3],Fp_q,eta_q;

            /* compute quadrature coordinates */
            pos_q[0] = pos_q[1] = pos_q[2] = 0.0;
            for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
                pos_q[0] += NI[q][k] * el_coord[3*k+0];
                pos_q[1] += NI[q][k] * el_coord[3*k+1];
            }
            
            /* evaluate */
            ierr = MMSStokes2D_Fu(pos_q[0],pos_q[1],Fu_q);CHKERRQ(ierr);
            ierr = MMSStokes2D_Fp(pos_q[0],pos_q[1],&Fp_q);CHKERRQ(ierr);
            
            cell_quadpoints[q].Fu[0] = -Fu_q[0];
            cell_quadpoints[q].Fu[1] = -Fu_q[1];
            cell_quadpoints[q].Fu[2] = 0.0;

            cell_quadpoints[q].Fp    = Fp_q;
            
            ierr = MMSStokes2D_eta(pos_q[0],pos_q[1],&eta_q);CHKERRQ(ierr);
            cell_quadpoints[q].eta = eta_q;
            
        }
    }
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyInitialMaterialGeometry_MMS"
PetscErrorCode ModelApplyInitialMaterialGeometry_MMS(pTatinCtx c,void *ctx)
{
    PetscErrorCode ierr;
	MMSCtx         user = (MMSCtx)ctx;
    
    switch (user->mms_idx) {
        case 1:
            ierr = MMSStokes1_ModelApplyInitialMaterialGeometry(c,user);CHKERRQ(ierr);
            break;
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_MMS_u"
PetscBool BCListEvaluator_MMS_u(PetscScalar position[],PetscScalar *value,void *ctx)
{
    PetscErrorCode ierr;
	MMSCtx         user;
    PetscScalar    u_vec[2];
    
    user = (MMSCtx)ctx;
    switch (user->mms_idx) {
        case 1:
            ierr = MMSStokes2D_u(position[0],position[1],u_vec);CHKERRQ(ierr);
            break;
    }
    *value = (PetscScalar)u_vec[0];
    
	return PETSC_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_MMS_v"
PetscBool BCListEvaluator_MMS_v(PetscScalar position[],PetscScalar *value,void *ctx)
{
    PetscErrorCode ierr;
	MMSCtx         user;
    PetscScalar    u_vec[2];
    
    user = (MMSCtx)ctx;
    switch (user->mms_idx) {
        case 1:
            ierr = MMSStokes2D_u(position[0],position[1],u_vec);CHKERRQ(ierr);
            break;
    }
    *value = (PetscScalar)u_vec[1];
    
	return PETSC_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "BCListEvaluator_MMS_w"
PetscBool BCListEvaluator_MMS_w(PetscScalar position[],PetscScalar *value,void *ctx)
{
	MMSCtx         user;
    PetscScalar    u_vec[3];
    
    user = (MMSCtx)ctx;
    switch (user->mms_idx) {
        case 1:
            u_vec[0] = u_vec[1] = u_vec[2] = 0.0;
            break;
    }
    *value = (PetscScalar)u_vec[2];
    
	return PETSC_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "MMS_BoundaryCondition"
PetscErrorCode MMS_BoundaryCondition(BCList bclist,DM dav,MMSCtx mms)
{
    PetscErrorCode ierr;

    if (mms->geometry_quasi_2d) {
        PetscReal zero;
    
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);

        zero = 0.0;
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        /*// block out vz on all other faces //
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
        */
    } else {
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,2,BCListEvaluator_MMS_w,(void*)mms);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,2,BCListEvaluator_MMS_w,(void*)mms);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,2,BCListEvaluator_MMS_w,(void*)mms);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,2,BCListEvaluator_MMS_w,(void*)mms);CHKERRQ(ierr);

        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_MMS_w,(void*)mms);CHKERRQ(ierr);
        
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,0,BCListEvaluator_MMS_u,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,1,BCListEvaluator_MMS_v,(void*)mms);CHKERRQ(ierr);
        ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_MMS_w,(void*)mms);CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryConditionMG_MMS"
PetscErrorCode ModelApplyBoundaryConditionMG_MMS(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
    PetscErrorCode ierr;
    PetscInt       n;
	MMSCtx         user = (MMSCtx)ctx;
    
    for (n=0; n<nl; n++) {
        ierr = MMS_BoundaryCondition(bclist[n],dav[n],user);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelApplyBoundaryCondition_MMS"
PetscErrorCode ModelApplyBoundaryCondition_MMS(pTatinCtx c,void *ctx)
{
    PetscErrorCode ierr;
    PhysCompStokes stokes;
    BCList         bclist;
    DM             dav;
	MMSCtx         user = (MMSCtx)ctx;

    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
    bclist = stokes->u_bclist;
    ierr = PhysCompStokesGetDMs(stokes,&dav,NULL);CHKERRQ(ierr);
    
    ierr = MMS_BoundaryCondition(bclist,dav,user);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MMSComputeErrors"
PetscErrorCode MMSComputeErrors(pTatinCtx c,Vec X,MMSCtx user,
                                PetscReal *_h,PetscReal *_uL2,PetscReal *_uH1,PetscReal *_pL2,PetscReal *_etaL2)
{
    PetscErrorCode ierr;
    DM stokes_pack,dav,dap,cda;
    Vec velocity,pressure,gcoords;
    PetscScalar *LA_ufield,*LA_pfield,*LA_gcoords;
    PetscInt k,e,q,nqp,nel,nen_u,nen_p;
    const PetscInt *elnidx_u;
    const PetscInt *elnidx_p;
    PetscReal l_data[7],g_data[7];
    QPntVolCoefStokes *quadpoints,*cell_quadpoints;
    Quadrature volQ;
    PhysCompStokes stokes;
	PetscReal WEIGHT[NQP],XI[NQP][3],NI[NQP][Q2_NODES_PER_EL_3D],NIp[NQP][P_BASIS_FUNCTIONS],GNI[NQP][3][Q2_NODES_PER_EL_3D];
    PetscReal GNx[NQP][Q2_NODES_PER_EL_3D],GNy[NQP][Q2_NODES_PER_EL_3D],GNz[NQP][Q2_NODES_PER_EL_3D],detJ[NQP];
    PetscReal el_coord[3*Q2_NODES_PER_EL_3D];
	PetscReal elu[3*Q2_NODES_PER_EL_3D],elp[P_BASIS_FUNCTIONS];
	PetscReal ux[Q2_NODES_PER_EL_3D],uy[Q2_NODES_PER_EL_3D],uz[Q2_NODES_PER_EL_3D];
    PetscReal l_uL2,l_uH1,l_pL2,l_etaL2,l_pL1,l_p_mean,l_vol;
    PetscReal min_dh[3],max_dh[3];
    PetscReal p_mean,p_corrected;
    
    
    ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
    
    ierr = PhysCompStokesGetDMs(stokes,&dav,&dap);CHKERRQ(ierr);

    /* velocity/pressure setup */
	stokes_pack = stokes->stokes_pack;
    ierr = DMCompositeGetLocalVectors(stokes_pack,&velocity,&pressure);CHKERRQ(ierr);
	ierr = DMCompositeScatter(stokes_pack,X,velocity,pressure);CHKERRQ(ierr);
	ierr = VecGetArray(velocity,&LA_ufield);CHKERRQ(ierr);
	ierr = VecGetArray(pressure,&LA_pfield);CHKERRQ(ierr);
    
	/* determine min/max dx,dy,dz for mesh */
	ierr = DMDAComputeQ2ElementBoundingBox(dav,min_dh,max_dh);CHKERRQ(ierr);
    
    /* coordinate setup */
	ierr = DMGetCoordinateDM(dav,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dav,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);

	ierr = DMDAGetElements_pTatinQ2P1(dav,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	ierr = DMDAGetElements_pTatinQ2P1(dap,&nel,&nen_p,&elnidx_p);CHKERRQ(ierr);
    
    /* quadrature setup */
    ierr = PhysCompStokesGetVolumeQuadrature(stokes,&volQ);CHKERRQ(ierr);
    ierr = VolumeQuadratureGetInfo(volQ,&nqp,&nel);CHKERRQ(ierr);
    P3D_prepare_elementQ2(nqp,WEIGHT,XI,NI,GNI);
    ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&quadpoints);CHKERRQ(ierr);
    
    l_uL2 = l_uH1 = l_pL2 = l_etaL2 = l_pL1 = 0.0;
    l_p_mean = l_vol = 0.0;
    for (e=0; e<nel; e++) {
        /* get element coordinates */
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coord,(PetscInt*)&elnidx_u[nen_u*e],LA_gcoords);CHKERRQ(ierr);
        /* get element velocity */
		ierr = DMDAGetVectorElementFieldQ2_3D(elu,(PetscInt*)&elnidx_u[nen_u*e],LA_ufield);CHKERRQ(ierr);
		for (k=0; k<Q2_NODES_PER_EL_3D; k++ ) {
			ux[k] = elu[3*k  ];
			uy[k] = elu[3*k+1];
			uz[k] = elu[3*k+2];
		}
        /* get element pressure */
		ierr = DMDAGetScalarElementField(elp,nen_p,(PetscInt*)&elnidx_p[nen_p*e],LA_pfield);CHKERRQ(ierr);
        /* get pressure basis functions */
		for (q=0; q<nqp; q++) {
			PetscScalar xip[] = { XI[q][0], XI[q][1], XI[q][2] };
			ConstructNi_pressure(xip,el_coord,NIp[q]);
		}
        /* get element quadrature */
		ierr = VolumeQuadratureGetCellData_Stokes(volQ,quadpoints,e,&cell_quadpoints);CHKERRQ(ierr);
		
        P3D_evaluate_geometry_elementQ2(nqp,el_coord,GNI,detJ,GNx,GNy,GNz);
        
        for (q=0; q<nqp; q ++) {
            PetscReal pos_q[3],u_q[3],gradu2d_q[4],gradu_q[9],p_q,eta_q;
            PetscReal u_fe_q[3],gradu_fe_q[9],p_fe_q,eta_fe_q;
            PetscReal err_gradu_L2,err_u_L2;
            PetscInt d;
            
            /* compute quadrature coordinates */
            pos_q[0] = pos_q[1] = pos_q[2] = 0.0;
            for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
                pos_q[0] += NI[q][k] * el_coord[3*k+0];
                pos_q[1] += NI[q][k] * el_coord[3*k+1];
                pos_q[2] += NI[q][k] * el_coord[3*k+2];
            }
            
            /* evaluate */
            u_q[0] = u_q[1] = u_q[2] = 0.0;
            ierr = MMSStokes2D_u(pos_q[0],pos_q[1],u_q);CHKERRQ(ierr);
            ierr = MMSStokes2D_gradu(pos_q[0],pos_q[1],gradu2d_q);CHKERRQ(ierr);
            ierr = MMSStokes2D_p(pos_q[0],pos_q[1],&p_q);CHKERRQ(ierr);
            ierr = MMSStokes2D_eta(pos_q[0],pos_q[1],&eta_q);CHKERRQ(ierr);
            gradu_q[0] = gradu_q[1] = gradu_q[2] = 0.0;
            gradu_q[3] = gradu_q[4] = gradu_q[5] = 0.0;
            gradu_q[6] = gradu_q[7] = gradu_q[8] = 0.0;
            gradu_q[0] = gradu2d_q[0]; //u_x
            gradu_q[1] = gradu2d_q[1]; //u_y
            gradu_q[3] = gradu2d_q[2]; //v_x
            gradu_q[4] = gradu2d_q[3]; //v_y
            
            /* interpolate velocity */
            u_fe_q[0] = u_fe_q[1] = u_fe_q[2] = 0.0;
            for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
                u_fe_q[0] += NI[q][k] * ux[k];
                u_fe_q[1] += NI[q][k] * uy[k];
                u_fe_q[2] += NI[q][k] * uz[k];
            }
            /* interpolate gradients velocity */
            for (d=0; d<9; d++) {
                gradu_fe_q[d] = 0.0;
            }
            for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
                gradu_fe_q[0] += GNx[q][k] * ux[k];
                gradu_fe_q[1] += GNy[q][k] * ux[k];
                gradu_fe_q[2] += GNz[q][k] * ux[k];
                
                gradu_fe_q[3] += GNx[q][k] * uy[k];
                gradu_fe_q[4] += GNy[q][k] * uy[k];
                gradu_fe_q[5] += GNz[q][k] * uy[k];

                gradu_fe_q[6] += GNx[q][k] * uz[k];
                gradu_fe_q[7] += GNy[q][k] * uz[k];
                gradu_fe_q[8] += GNz[q][k] * uz[k];
            }
            /* interpolate pressure */
            p_fe_q = 0.0;
			for (k=0; k<P_BASIS_FUNCTIONS; k++) {
				p_fe_q += NIp[q][k] * elp[k];
			}
            l_p_mean += p_fe_q * WEIGHT[q] * detJ[q];

            //p_fe_q = p_fe_q - ( -8.263367837e-01);
            //printf("p = %1.6e %1.6e [exact]\n",p_fe_q,p_q);
            /* get eta */
            eta_fe_q = cell_quadpoints[q].eta;

            err_u_L2 = 0.0;
            for (d=0; d<3; d++) {
                err_u_L2 += (u_fe_q[d] - u_q[d])*(u_fe_q[d] - u_q[d]);
            }
            err_gradu_L2 = 0.0;
            for (d=0; d<9; d++) {
                err_gradu_L2 += (gradu_fe_q[d] - gradu_q[d])*(gradu_fe_q[d] - gradu_q[d]);
            }
            
            l_uL2    += err_u_L2 * WEIGHT[q] * detJ[q];
            l_uH1    += err_gradu_L2 * WEIGHT[q] * detJ[q];
            l_pL2    += (p_fe_q - p_q)*(p_fe_q - p_q) * WEIGHT[q] * detJ[q];
            l_etaL2  += (eta_fe_q - eta_q)*(eta_fe_q - eta_q) * WEIGHT[q] * detJ[q];
            l_pL1    += (p_fe_q - p_q) * WEIGHT[q] * detJ[q];
            l_vol    += 1.0 * WEIGHT[q] * detJ[q];
        }
    }
    l_data[0] = l_uL2;
    l_data[1] = l_uH1;
    l_data[2] = l_pL2;
    l_data[3] = l_etaL2;
    l_data[4] = l_p_mean;
    l_data[5] = l_pL1;
    l_data[6] = l_vol;
    
    ierr = MPI_Allreduce(l_data,g_data,7,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
    
    {
        PetscPrintf(PETSC_COMM_WORLD,"(dx,dy,dz)_max       = %1.4e %1.4e %1.4e\n",max_dh[0],max_dh[1],max_dh[2]);
        PetscPrintf(PETSC_COMM_WORLD,"\\int 1 dv           = %+1.9e\n",g_data[6]);
        PetscPrintf(PETSC_COMM_WORLD,"\\int p dv           = %+1.9e\n",g_data[4]);
        PetscPrintf(PETSC_COMM_WORLD,"\\int p dv / Omega   = %+1.9e\n",g_data[4]/g_data[6]);
        PetscPrintf(PETSC_COMM_WORLD,"\\int (p^h - p) dv   = %+1.9e\n",g_data[5]);
        PetscPrintf(PETSC_COMM_WORLD,"\\int (p^h - p)^2 dv = %+1.9e\n",g_data[2]);
        
        p_mean = g_data[4]/g_data[6];
        
        /*
         \int (p(x)-a)(p(x)-a) dV
         = \int p(x)*p(x) dV - 2 * a * \int p(x) dV + \int a * a dV
        */
        p_corrected = g_data[2] - 2.0*p_mean*g_data[5] + g_data[5]*g_data[5]*g_data[6];
        p_corrected = PetscSqrtReal(p_corrected);
        PetscPrintf(PETSC_COMM_WORLD,"\\int (p^h - p)^2 dv [corrected] = %+1.9e\n",p_corrected);
        
    }
    
    
	*_h = -1.0e32;
	*_h = PetscMax(*_h,max_dh[0]); // x
	*_h = PetscMax(*_h,max_dh[1]); // y
    /*
    *_uL2   = PetscSqrtReal(g_data[0]/max_dh[2]);
    *_uH1   = PetscSqrtReal(g_data[1]/max_dh[2]);
    *_pL2   = PetscSqrtReal(g_data[2]/max_dh[2]);
    *_etaL2 = PetscSqrtReal(g_data[3]/max_dh[2]);
     */
    *_uL2   = PetscSqrtReal(g_data[0]);
    *_uH1   = PetscSqrtReal(g_data[1]);
    *_pL2   = p_corrected;
    *_etaL2 = PetscSqrtReal(g_data[3]);

	ierr = VecRestoreArray(pressure,&LA_pfield);CHKERRQ(ierr);
	ierr = VecRestoreArray(velocity,&LA_ufield);CHKERRQ(ierr);
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
    ierr = DMCompositeRestoreLocalVectors(stokes_pack,&velocity,&pressure);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelOutput_MMS"
PetscErrorCode ModelOutput_MMS(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
    PetscErrorCode ierr;
	MMSCtx         user = (MMSCtx)ctx;
    PetscReal      h,uL2,uH1,pL2,etaL2;
    
    ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"--------------- MMS summary ---------------\n");
    ierr = MMSComputeErrors(c,X,user,&h,&uL2,&uH1,&pL2,&etaL2);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"h %1.6e | uL2 %1.6e : uH1 %1.6e : pL2 %1.6e : etaL2 %1.6e\n",h,uL2,uH1,pL2,etaL2);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ModelDestroy_MMS"
PetscErrorCode ModelDestroy_MMS(pTatinCtx c,void *ctx)
{
    MMSCtx         user = (MMSCtx)ctx;
    PetscErrorCode ierr;

    ierr = PetscFree(user);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "pTatinModelRegister_MMS"
PetscErrorCode pTatinModelRegister_MMS(void)
{
    MMSCtx         user;
    pTatinModel    model;
    PetscErrorCode ierr;
    
    
    /* Allocate memory for the data structure for this model */
    ierr = PetscNew(&user);CHKERRQ(ierr);
    
    /* create model type and register type */
    ierr = pTatinModelCreate(&model);CHKERRQ(ierr);
    ierr = pTatinModelSetName(model,"mms");CHKERRQ(ierr);
    ierr = pTatinModelSetUserData(model,user);CHKERRQ(ierr);
    
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_MMS);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelInitialCondition_MMS);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_MMS);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_MMS);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_MMS);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_MMS);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_MMS);CHKERRQ(ierr);
    ierr = pTatinModelSetFunctionPointer(model,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_MMS);CHKERRQ(ierr);
	
	/* Insert model into list */
	ierr = pTatinModelRegister(model);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

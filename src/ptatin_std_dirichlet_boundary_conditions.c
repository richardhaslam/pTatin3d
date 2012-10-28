/*

 Collection of standard boundary conditions for hex meshes.

 Definitions:
 Regarding directions, we will adopt the following nomenclecture
 NX,NY,NZ are the number of nodes in x,y,z direction respectively

 NorthFace: Implies boundary defined by all nodes with ny = NY-1
 
 EastFace:  Implies boundary defined by all nodes with nx = NX-1
 
 SouthFace: Implies boundary defined by all nodes with ny = 0
 
 WestFace:  Implies boundary defined by all nodes with nx = 0

 FrontFace: Implies boundary defined by all nodes with nz = NZ-1
 
 BackFace:  Implies boundary defined by all nodes with nz = 0
 
*/

#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "dmda_bcs.h"
#include "ptatin_std_dirichlet_boundary_conditions.h"


#undef __FUNCT__
#define __FUNCT__ "DirichletBC_FreeSlip"
PetscErrorCode DirichletBC_FreeSlip(BCList bclist,DM dav,BoundaryFaceType face)
{
	PetscScalar    zero = 0.0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	switch (face) {

		case EAST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case WEST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case NORTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;
			
		case SOUTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case FRONT_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;

		case BACK_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&zero);CHKERRQ(ierr);
			break;
	}
		
	PetscFunctionReturn(0);	
}

/* 
 Note: This is not general and will not work on deformed boundaries.
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyNormalVelocity"
PetscErrorCode DirichletBC_ApplyNormalVelocity(BCList bclist,DM dav,BoundaryFaceType face,PetscReal v_normal)
{
	PetscScalar    value = v_normal;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	switch (face) {
			
		case EAST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case WEST_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case NORTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case SOUTH_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case FRONT_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMAX_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
			
		case BACK_FACE:
			ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_KMIN_LOC,2,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
			break;
	}
	
	PetscFunctionReturn(0);	
}

/* 
 Extend mesh in x direction with a specified strain rate, Exx.
 The applied velocity has the same magnitude in east and west directions.

 Compute size of mesh (Lx) in x direction.
 Exx = du/dx = (vx_east - vx_west)/Lx = 2 Vx_bc / Lx
 Vx_bc = 0.5 . Exx . Lx 
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyStrainRateExx"
PetscErrorCode DirichletBC_ApplyStrainRateExx(BCList bclist,DM dav,PetscReal exx_bc)
{
	PetscScalar value;
	PetscReal MeshMin[3],MeshMax[3];
	PetscReal Lx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	Lx = (MeshMax[0] - MeshMin[0]);
	value = 0.5 * exx_bc * Lx;

	
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);

	value = -value;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_IMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);

	PetscFunctionReturn(0);	
}

/* 
 Apply tangential velocity along the north/south faces via a specified shear strain rate, Exy.
 Normal velocity on north/south face is set to be zero
 The applied velocity has the same magnitude in east and west directions.
 
 Compute size of mesh (Ly) in y direction.
 Exy = 0.5.(du/dy + dv/dx)= 0.5.(vx_east - vx_west)/Lx = Vx_bc / Ly
 Vx_bc = Exy . Ly
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyStrainRateExy"
PetscErrorCode DirichletBC_ApplyStrainRateExy(BCList bclist,DM dav,PetscReal exy_bc)
{
	PetscScalar value;
	PetscReal MeshMin[3],MeshMax[3];
	PetscReal Ly;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	Ly = (MeshMax[1] - MeshMin[1]);

	value = exy_bc * Ly;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,1,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	value = -value;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,1,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	
	value = 0.0;
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMAX_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	ierr = DMDABCListTraverse3d(bclist,dav,DMDABCList_JMIN_LOC,0,BCListEvaluator_constant,(void*)&value);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

/* 
 Extend/compress in x direction, compress/extend in y direction.

 v_shear(x,y) = value.(x.i - y.j)  
 
 Apply x component of v_shear() along east-west boundary faces
 Apply y component of v_shear() along north-south boundary faces 
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningY"
PetscErrorCode DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningY(BCList list,DM dav,PetscReal vx_bc)
{
	PetscScalar v_normal;
	PetscReal MeshMin[3],MeshMax[3];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	
	
	v_normal = vx_bc * MeshMax[0];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,EAST_FACE,v_normal);CHKERRQ(ierr);
	v_normal = vx_bc * MeshMin[0];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,WEST_FACE,v_normal);CHKERRQ(ierr);

	v_normal = -vx_bc * MeshMax[1];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,NORTH_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc * MeshMin[1];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,SOUTH_FACE,v_normal);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

/* 
 Extend/compress in x direction, compress/extend in z direction.
 
 v_shear(x,z) = value.(x.i - z.k)  
 
 Apply x component of v_shear() along east-west boundary faces
 Apply z component of v_shear() along front-back boundary faces 
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ"
PetscErrorCode DirichletBC_ApplyConstantAreaSection_ExtensionX_ShorteningZ(BCList list,DM dav,PetscReal vx_bc)
{
	PetscScalar v_normal;
	PetscReal MeshMin[3],MeshMax[3];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	
	
	v_normal = vx_bc * MeshMax[0];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,EAST_FACE,v_normal);CHKERRQ(ierr);
	v_normal = vx_bc * MeshMin[0];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,WEST_FACE,v_normal);CHKERRQ(ierr);
	
	v_normal = -vx_bc * MeshMax[2];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,FRONT_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc * MeshMin[2];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,BACK_FACE,v_normal);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

/* 
 Extend/compress in x direction, compress/extend in y,z direction.
 
 v_shear(x,y,z) = value.(x.i - 0.5.y.j - 0.5.z.k)  
 
 Apply x component of v_shear() along east-west boundary faces
 Apply y component of v_shear() along north-south boundary faces 
 Apply z component of v_shear() along front-back boundary faces 
*/
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantVolumeDomain_ExtensionX"
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionX(BCList list,DM dav,PetscReal vx_bc)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(list,dav,0.5,vx_bc);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);	
}

/* 
 Extend/compress in x direction, compress/extend in y,z direction.
 
 v_shear(x,y,z) = value.(x.i - beta.y.j - (1-beta).z.k)  
 
 Apply x component of v_shear() along east-west boundary faces
 Apply y component of v_shear() along north-south boundary faces 
 Apply z component of v_shear() along front-back boundary faces 
 */
#undef __FUNCT__
#define __FUNCT__ "DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening"
PetscErrorCode DirichletBC_ApplyConstantVolumeDomain_ExtensionXFractionShortening(BCList list,DM dav,PetscReal beta,PetscReal vx_bc)
{
	PetscScalar v_normal;
	PetscReal MeshMin[3],MeshMax[3];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = DMDAGetBoundingBox(dav,MeshMin,MeshMax);CHKERRQ(ierr);
	
	v_normal = vx_bc * MeshMax[0];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,EAST_FACE,v_normal);CHKERRQ(ierr);
	v_normal = vx_bc * MeshMin[0];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,WEST_FACE,v_normal);CHKERRQ(ierr);
	
	v_normal = -vx_bc * beta * MeshMax[1];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,NORTH_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc * beta * MeshMin[1];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,SOUTH_FACE,v_normal);CHKERRQ(ierr);
	
	v_normal = -vx_bc * (1.0-beta) * MeshMax[2];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,FRONT_FACE,v_normal);CHKERRQ(ierr);
	v_normal = -vx_bc * (1.0-beta) * MeshMin[2];
	ierr = DirichletBC_ApplyNormalVelocity(list,dav,BACK_FACE,v_normal);CHKERRQ(ierr);

	PetscFunctionReturn(0);	
}




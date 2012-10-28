

#include "stdio.h"
#include "stdlib.h"
#include "petsc.h"
#include "petscdm.h"

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_element_q2p1.h"
#include "quadrature.h"
#include "element_utils_q2.h"
#include "dmda_update_coords.h"
#include "element_utils_q1.h"

#include "mesh_update.h"

#define Q1_NODES_PER_EL_3D 8

#undef __FUNCT__  
#define __FUNCT__ "DMDABilinearizeQ2Elements"
PetscErrorCode DMDABilinearizeQ2Elements(DM dau)
{
	DM cda;
	Vec gcoords;
	PetscScalar *LA_gcoords;
	PetscInt nel,nen,e,n,k,ii,jj,kk;
	const PetscInt *elnidx;
	PetscScalar elcoordsQ2[3*Q2_NODES_PER_EL_3D];
	PetscScalar elcoordsQ1[3*Q1_NODES_PER_EL_3D],Ni[Q1_NODES_PER_EL_3D];
	PetscInt cnt;
	PetscScalar xi_nodal_coordsQ2[3*Q2_NODES_PER_EL_3D],x_new[3];
	PetscInt vel_el_lidx[3*Q2_NODES_PER_EL_3D];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* define some xi coords to interpolate to - these correspond to the nodal location of the q2 basis functions */
	/* this loop should be ordered the same way as the natural nodal ordering of the Q2 nodes in element space */
	cnt = 0;
	for (kk=0; kk<3; kk++) {
		for (jj=0; jj<3; jj++) {
			for (ii=0; ii<3; ii++) {
				xi_nodal_coordsQ2[3*cnt+0] = -1.0 + (double)ii;
				xi_nodal_coordsQ2[3*cnt+1] = -1.0 + (double)jj;
				xi_nodal_coordsQ2[3*cnt+2] = -1.0 + (double)kk;
				cnt++;
			}
		}
	}	
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(dau,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dau,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen,&elnidx);CHKERRQ(ierr);
	for (e=0;e<nel;e++) {
		ierr = DMDAGetElementCoordinatesQ2_3D(elcoordsQ2,(PetscInt*)&elnidx[nen*e],LA_gcoords);CHKERRQ(ierr);
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		for (k=0; k<3; k++) {
			elcoordsQ1[3*0+k] = elcoordsQ2[3*0+k];
			elcoordsQ1[3*1+k] = elcoordsQ2[3*2+k];
			elcoordsQ1[3*2+k] = elcoordsQ2[3*6+k];
			elcoordsQ1[3*3+k] = elcoordsQ2[3*8+k];
			
			elcoordsQ1[3*4+k] = elcoordsQ2[3*18+k];
			elcoordsQ1[3*5+k] = elcoordsQ2[3*20+k];
			elcoordsQ1[3*6+k] = elcoordsQ2[3*24+k];
			elcoordsQ1[3*7+k] = elcoordsQ2[3*26+k];
		}
		
		/* for each interior point */
		for (k=0; k<27; k++) {
			P3D_ConstructNi_Q1_3D(&xi_nodal_coordsQ2[3*k],Ni);
			
			/* inpterpolate */
			x_new[0] = 0.0;
			x_new[1] = 0.0;
			x_new[2] = 0.0;
			for (n=0; n<Q1_NODES_PER_EL_3D; n++) {
				x_new[0] += elcoordsQ1[3*n+0] * Ni[n];
				x_new[1] += elcoordsQ1[3*n+1] * Ni[n];
				x_new[2] += elcoordsQ1[3*n+2] * Ni[n];
			}
			
			elcoordsQ2[3*k+0] = x_new[0];
			elcoordsQ2[3*k+1] = x_new[1];
			elcoordsQ2[3*k+2] = x_new[2];
			//printf("e-%d: %1.4e  %1.4e  %1.4e \n", e,x_new[0],x_new[1],x_new[2] );
		}
		
		/* push modification */
		ierr = DMDASetValuesLocalStencil_InsertValues_Stokes_Velocity(LA_gcoords, vel_el_lidx,elcoordsQ2);CHKERRQ(ierr);
		
	}
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	/* send ghostes values into global vector */
	ierr = DMDASetCoordinatesFromLocalVector(dau,gcoords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}




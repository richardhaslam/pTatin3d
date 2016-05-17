
#include <petscvec.h>
#include <petscmat.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petsc/private/dmimpl.h>
#include <petsc/private/dmdaimpl.h>
#include <private/dafeimpl.h>
#include <dmda_update_coords.h>

PetscErrorCode _DAFEProjectCoordinates_Q1ToQ2(DM dafeq2,DM dafeq1);


#undef __FUNCT__
#define __FUNCT__ "DAFE_EvaluateBasisQ2_3D"
PetscErrorCode DAFE_EvaluateBasisQ2_3D(PetscReal elcoor[],PetscReal _xi[],PetscReal Ni[])
{
	PetscInt i,j,k,d,cnt;
	PetscReal basis_NI[3][3];
	
	for (d=0; d<3; d++) {
		double xi = _xi[d];
		
		basis_NI[d][0] = 0.5 * xi * (xi-1.0); // 0.5 * ( xi^2 - xi )
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); // 1 - xi^2
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; // 0.5 * ( xi^2 + xi )
	}
	
	cnt = 0;
	for (k=0; k<3; k++) {
		for (j=0; j<3; j++) {
			for (i=0; i<3; i++) {
				Ni[cnt] = basis_NI[0][i] * basis_NI[1][j] * basis_NI[2][k];
				cnt++;
			}
		}
	}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetParentElementIndex_Q1sub2Q2"
PetscErrorCode GetParentElementIndex_Q1sub2Q2(DM sub,PetscInt e_sub,PetscInt *e_p)
{
  DM_DAFE        *feq1 = (DM_DAFE*)sub->data;
  DM_DAFE        *feq2 = (DM_DAFE*)feq1->parent_dafe->data;
  PetscInt       ref,e_sub_ijk[3],e2d,e_ijk[3];

  /* convert e_sub -> (i,j,k) */
  e_sub_ijk[2] = e_sub / (feq1->lmx * feq1->lmy);
  e2d = e_sub - e_sub_ijk[2] * (feq1->lmx * feq1->lmy);
  e_sub_ijk[1] = e2d / feq1->lmx;
  e_sub_ijk[0] = e2d - e_sub_ijk[1] * feq1->lmx;
  
  ref = feq1->mx / feq2->mx;
  
  e_ijk[0] = e_sub_ijk[0] / ref;
  e_ijk[1] = e_sub_ijk[1] / ref;
  e_ijk[2] = e_sub_ijk[2] / ref;
  
  *e_p = e_ijk[0] + e_ijk[1] * feq2->lmx + e_ijk[2] * (feq2->lmx * feq2->lmy);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ConvertChildLocalCoordinateToParentElement_Ref1"
PetscErrorCode ConvertChildLocalCoordinateToParentElement_Ref1(DM sub,PetscInt e_sub,PetscReal xic[],PetscReal xi[])
{
  xi[0] = xic[0];
  xi[1] = xic[1];
  xi[2] = xic[2];
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ConvertChildLocalCoordinateToParentElement_Ref2"
PetscErrorCode ConvertChildLocalCoordinateToParentElement_Ref2(DM sub,PetscInt e_sub,PetscReal xic[],PetscReal xi[])
{
  /*
   [xic - (-1)] / 1 = [xi - (-1)] / 2
   */
  xi[0] = 2.0 * (xic[0] + 1.0) -1.0;
  xi[1] = 2.0 * (xic[1] + 1.0) -1.0;
  xi[2] = 2.0 * (xic[2] + 1.0) -1.0;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ConvertChildLocalCoordinateToParentElement_Q1sub2Q2"
PetscErrorCode ConvertChildLocalCoordinateToParentElement_Q1sub2Q2(DM sub,PetscInt e_sub,PetscReal xic[],PetscReal xi[])
{
  DM_DAFE        *feq1 = (DM_DAFE*)sub->data;
  DM_DAFE        *feq2 = (DM_DAFE*)feq1->parent_dafe->data;
  PetscInt       k,ref,e_sub_ijk[3],e2d,e_ijk[3],sub_ijk[3];
  PetscReal      dxi,ref_points_x1_1d[21];
  
  
  /* convert e_sub -> (i,j,k) */
  e_sub_ijk[2] = e_sub / (feq1->lmx * feq1->lmy);
  e2d = e_sub - e_sub_ijk[2] * (feq1->lmx * feq1->lmy);
  e_sub_ijk[1] = e2d / feq1->lmx;
  e_sub_ijk[0] = e2d - e_sub_ijk[1] * feq1->lmx;
  
  ref = feq1->mx / feq2->mx;

  if (ref > 20) SETERRQ(PetscObjectComm((PetscObject)sub),PETSC_ERR_SUP,"Refinement factors > 20 are not supported");
  
  e_ijk[0] = e_sub_ijk[0] / ref;
  e_ijk[1] = e_sub_ijk[1] / ref;
  e_ijk[2] = e_sub_ijk[2] / ref;

  sub_ijk[0] = e_sub_ijk[0] - e_ijk[0] * ref;
  sub_ijk[1] = e_sub_ijk[1] - e_ijk[1] * ref;
  sub_ijk[2] = e_sub_ijk[2] - e_ijk[2] * ref;
  
  
  dxi = 2.0 / ((PetscReal)ref);
  for (k=0; k<ref+1; k++) { ref_points_x1_1d[k] = -1.0 + k * dxi; }
  
  /*
  [xic - (-1)] / 2 = [xi - (ref_points_x1_1d)] / dxi
  */
  
  xi[0] = 0.5*(xic[0]+1.0)*dxi + ref_points_x1_1d[sub_ijk[0]];
  xi[1] = 0.5*(xic[1]+1.0)*dxi + ref_points_x1_1d[sub_ijk[1]];
  xi[2] = 0.5*(xic[2]+1.0)*dxi + ref_points_x1_1d[sub_ijk[2]];
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAFE_GetElementCoords_3D"
PetscErrorCode DAFE_GetElementCoords_3D(DM dm,PetscInt e,PetscReal elcoor[])
{
  DM_DAFE  *fe = (DM_DAFE*)dm->data;
  PetscInt *element,k;
  Vec coor;
  const PetscScalar *LA_coor;
  PetscErrorCode ierr;
  
  element = &fe->element_node_map[e*fe->nodes_per_el];
  ierr = DMGetCoordinatesLocal(fe->da,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  for (k=0; k<fe->nodes_per_el; k++) {
    PetscInt nidx = element[k];
    
    elcoor[3*k+0] = LA_coor[3*nidx+0];
    elcoor[3*k+1] = LA_coor[3*nidx+1];
    elcoor[3*k+2] = LA_coor[3*nidx+2];
  }
  ierr = VecRestoreArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAFE_GetElementNodeMap"
PetscErrorCode DAFE_GetElementNodeMap(DM dm,PetscInt e,const PetscInt **map)
{
  DM_DAFE  *fe = (DM_DAFE*)dm->data;
  
  if (!fe->element_node_map) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"DAFE does not provide element_node_map[]");
  *map = &fe->element_node_map[e*fe->nodes_per_el];
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAFE_GetElementBasisMap"
PetscErrorCode DAFE_GetElementBasisMap(DM dm,PetscInt e,const PetscInt **map)
{
  DM_DAFE  *fe = (DM_DAFE*)dm->data;
  
  if (!fe->element_basis_map) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"DAFE does not provide element_basis_map[]");
  *map = &fe->element_basis_map[e*fe->nbasis];
  
  PetscFunctionReturn(0);
}


/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DAFE_GetLocalSizeElementQ2"
PetscErrorCode DAFE_GetLocalSizeElementQ2(DM da,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	PetscInt i,j,k,start;
	PetscInt cntx,cnty,cntz;
	PetscInt si,sj,sk,m,n,p,M,N,P,width;
	PetscInt sig,sjg,skg,mg,ng,pg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sig,&sjg,&skg,&mg,&ng,&pg);CHKERRQ(ierr);
	if (width != 2) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Stencil width must be 2 for Q2");

	cntx = cnty = cntz = 0;
  
	
	/* ======================================================================================== */
	// x
	start = -1;
	for (i=si; i<si+m; i++) {
		if (i%2 == 0 && i == si && i != 0) { continue; } /* reject first ghost if its's even */
		if (i%2 == 0) {
			start = i;
			break;
		}
	}
	if (start == -1) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot determine start index in I"); }
	while (start + 2 * cntx < si+m) {
		PetscInt n0,n2;
		
		n0 = start + 2 * cntx;
		n2 = n0 + 2;
		
		if (n2 < sig+mg) {
			cntx++;
			continue;
		}
		
		if (si+m-n2 > 1) {
			cntx++;
			continue;
		}
		
		if (si+m-n2 <= 1) {
			break;
		}
	}
	/* ======================================================================================== */
	// y
	start = -1;
	for (j=sj; j<sj+n; j++) {
		if (j%2 == 0 && j == sj && j != 0) { continue; } /* reject first ghost if its's even */
		if (j%2 == 0) {
			start = j;
			break;
		}
	}
	if (start == -1) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot determine start index in J"); }
	while (start + 2 * cnty < sj+n) {
		PetscInt n0,n2;
		
		n0 = start + 2 * cnty;
		n2 = n0 + 2;
		
		/* if start and end of element are inside global range - keep it */
		if (n2 < sjg+ng) {
			cnty++;
			continue;
		}
		
		if (sj+n-n2 > 1) {
			cnty++;
			continue;
		}
		
		if (sj+n-n2 <= 1) {
			break;
		}
	}
	
	/* ======================================================================================== */
	// z
	start = -1;
	for (k=sk; k<sk+p; k++) {
		if (k%2 == 0 && k == sk && k != 0) { continue; } /* reject first ghost if its's even */
		if (k%2 == 0) {
			start = k;
			break;
		}
	}
	if (start == -1) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot determine start index in K"); }
	while (start + 2 * cntz < p+sk) {
		PetscInt n0,n2;
		
		n0 = start + 2 * cntz;
		n2 = n0 + 2;
		
		/* if start and end of element are inside global range - keep it */
		if (n2 < skg+pg) {
			cntz++;
			continue;
		}
		
		if (sk+p-n2 > 1) {
			cntz++;
			continue;
		}
		
		if (sk+p-n2 <= 1) {
			break;
		}
	}
	
	if (mx) { *mx = cntx; }
	if (my) { *my = cnty; }
	if (mz) { *mz = cntz; }
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DAFE_GetCornersElementQ2"
PetscErrorCode DAFE_GetCornersElementQ2(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	PetscInt i,j,k;
	PetscInt si,sj,sk,m,n,p,M,N,P,width;
	PetscInt sig,sjg,skg,mg,ng,pg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sig,&sjg,&skg,&mg,&ng,&pg);CHKERRQ(ierr);
	if (width != 2) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Stencil width must be 2 for Q2");
	
	// x
	for (i=si; i<si+m; i++) {
		if (i%2 == 0 && i == si && i != 0) { continue; } /* reject first ghost if its's even */
		if (i%2 == 0) {
			*sei = i;
			break;
		}
	}
  
	// y
	for (j=sj; j<sj+n; j++) {
		if (j%2 == 0 && j == sj && j != 0) { continue; } /* reject first ghost if its's even */
		if (j%2 == 0) {
			*sej = j;
			break;
		}
	}
	
	// z
	for (k=sk; k<sk+p; k++) {
		if (k%2 == 0 && k == sk && k != 0) { continue; } /* reject first ghost if its's even */
		if (k%2 == 0) {
			*sek = k;
			break;
		}
	}
	ierr = DAFE_GetLocalSizeElementQ2(da,mx,my,mz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* DA Q2 1D,2D,3D */
#undef __FUNCT__
#define __FUNCT__ "DAFE_GetOwnershipRangesElementQ2"
PetscErrorCode DAFE_GetOwnershipRangesElementQ2(DM da,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt **si,PetscInt **sj,PetscInt **sk,PetscInt **_mx,PetscInt **_my,PetscInt **_mz)
{
	PetscMPIInt nproc;
	MPI_Comm comm;
	PetscInt M,N,P,pM,pN,pP;
	PetscInt i,j,k,dim,esi,esj,esk,mx,my,mz;
	PetscInt *olx,*oly,*olz;
	PetscInt *lmx,*lmy,*lmz,*tmp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
	ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo( da, &dim, &M,&N,&P, &pM,&pN,&pP, 0, 0, 0,0,0, 0 );CHKERRQ(ierr);
	ierr = DAFE_GetCornersElementQ2(da,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
	if (dim == 1) {
		pN = 1;
		pP = 1;
	}
	if (dim == 2) {
		pP = 1;
	}
	
	ierr = PetscMalloc( sizeof(PetscInt)*(nproc), &tmp );CHKERRQ(ierr);
	
	ierr = PetscMalloc( sizeof(PetscInt)*(pM+1), &olx );CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*(pN+1), &oly );CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*(pP+1), &olz );CHKERRQ(ierr);
	
	ierr = PetscMalloc( sizeof(PetscInt)*(pM+1), &lmx );CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*(pN+1), &lmy );CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*(pP+1), &lmz );CHKERRQ(ierr);
	
	if (dim >= 1) {
		ierr = MPI_Allgather ( &esi, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		j = k = 0;
		for (i=0; i<pM; i++) {
			PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
			olx[i] = tmp[procid];
		}
		
		ierr = MPI_Allgather ( &mx, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		j = k = 0;
		for (i=0; i<pM; i++) {
			PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
			lmx[i] = tmp[procid];
		}
	}
	
	if (dim >= 2) {
		ierr = MPI_Allgather ( &esj, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		i = k = 0;
		for (j=0; j<pN; j++) {
			PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
			oly[j] = tmp[procid];
		}
		
		ierr = MPI_Allgather ( &my, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		i = k = 0;
		for (j=0; j<pN; j++) {
			PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
			lmy[j] = tmp[procid];
		}
	}
	
	if (dim == 3) {
		ierr = MPI_Allgather ( &esk, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		i = j = 0;
		for (k=0; k<pP; k++) {
			PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
			olz[k] = tmp[procid];
		}
		
		ierr = MPI_Allgather ( &mz, 1, MPIU_INT, tmp, 1, MPIU_INT, comm );CHKERRQ(ierr);
		i = j = 0;
		for (k=0; k<pP; k++) {
			PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
			lmz[k] = tmp[procid];
		}
	}
	
	if(m) { *m = pM; }
	if(n) { *n = pN; }
	if(p) { *p = pP; }
	
	if(si) { *si = olx; } else { ierr = PetscFree(olx);CHKERRQ(ierr); }
	if(sj) { *sj = oly; } else { ierr = PetscFree(oly);CHKERRQ(ierr); }
	if(sk) { *sk = olz; } else { ierr = PetscFree(olz);CHKERRQ(ierr); }
	
	if(_mx) { *_mx = lmx; } else { ierr = PetscFree(lmx);CHKERRQ(ierr); }
	if(_my) { *_my = lmy; } else { ierr = PetscFree(lmy);CHKERRQ(ierr); }
	if(_mz) { *_mz = lmz; } else { ierr = PetscFree(lmz);CHKERRQ(ierr); }
	
	ierr = PetscFree(tmp);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAFE_GetElementsQ2_3D"
PetscErrorCode DAFE_GetElementsQ2_3D(DM dm,PetscInt *_npe,PetscInt *_nel,PetscInt **_eidx)
{
	const PetscInt order = 2;
	PetscErrorCode ierr;
	PetscInt mx,my,mz,npe,M,N,P,nel;
	PetscInt ei,ej,ek,i,j,k,elcnt,esi,esj,esk,gsi,gsj,gsk,nid[27],n,X,Y,Z,width;
	PetscInt *el,*eidx;
	PetscInt dof;
	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width != 2) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Stencil width must be 2 for Q2");
	
	npe = (order + 1)*(order + 1)*(order + 1);
  
  ierr = DAFE_GetCornersElementQ2(dm,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
  nel = mx*my*mz;

  ierr = PetscMalloc(sizeof(PetscInt)*(nel*npe+1),&eidx);CHKERRQ(ierr);
  ierr = DMDAGetGhostCorners(dm,&gsi,&gsj,&gsk, &X,&Y,&Z);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dm,0, 0,0,0, 0,0,0, &dof,0, 0,0,0, 0);CHKERRQ(ierr);
  
  elcnt = 0;
  
  for (ek=0; ek<mz; ek++) {
    k = esk-gsk + 2*ek;
    for (ej=0; ej<my; ej++) {
      j = esj-gsj + 2*ej;
      for (ei=0; ei<mx; ei++) {
        i = esi-gsi + 2*ei;
        
        el = &eidx[npe*elcnt];
        
        nid[ 0] = (i  ) + (j  ) *X  + (k  ) *X*Y;
        nid[ 1] = (i+1) + (j  ) *X  + (k  ) *X*Y;
        nid[ 2] = (i+2) + (j  ) *X  + (k  ) *X*Y;
        
        nid[ 3] = (i  ) + (j+1) *X  + (k  ) *X*Y;
        nid[ 4] = (i+1) + (j+1) *X  + (k  ) *X*Y;
        nid[ 5] = (i+2) + (j+1) *X  + (k  ) *X*Y;
        
        nid[ 6] = (i  ) + (j+2) *X  + (k  ) *X*Y;
        nid[ 7] = (i+1) + (j+2) *X  + (k  ) *X*Y;
        nid[ 8] = (i+2) + (j+2) *X  + (k  ) *X*Y;
        //
        nid[ 9] = (i  ) + (j  ) *X  + (k+1) *X*Y;
        nid[10] = (i+1) + (j  ) *X  + (k+1) *X*Y;
        nid[11] = (i+2) + (j  ) *X  + (k+1) *X*Y;
        
        nid[12] = (i  ) + (j+1) *X  + (k+1) *X*Y;
        nid[13] = (i+1) + (j+1) *X  + (k+1) *X*Y;
        nid[14] = (i+2) + (j+1) *X  + (k+1) *X*Y;
        
        nid[15] = (i  ) + (j+2) *X  + (k+1) *X*Y;
        nid[16] = (i+1) + (j+2) *X  + (k+1) *X*Y;
        nid[17] = (i+2) + (j+2) *X  + (k+1) *X*Y;
        //
        nid[18] = (i  ) + (j  ) *X  + (k+2) *X*Y;
        nid[19] = (i+1) + (j  ) *X  + (k+2) *X*Y;
        nid[20] = (i+2) + (j  ) *X  + (k+2) *X*Y;
        
        nid[21] = (i  ) + (j+1) *X  + (k+2) *X*Y;
        nid[22] = (i+1) + (j+1) *X  + (k+2) *X*Y;
        nid[23] = (i+2) + (j+1) *X  + (k+2) *X*Y;
        
        nid[24] = (i  ) + (j+2) *X  + (k+2) *X*Y;
        nid[25] = (i+1) + (j+2) *X  + (k+2) *X*Y;
        nid[26] = (i+2) + (j+2) *X  + (k+2) *X*Y;
        
        for (n=0; n<npe; n++) {
          if (nid[n] > M*N*P) {
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of global nodes");
          }
          if (nid[n] > X*Y*Z) {
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of local nodes");
          }
          el[n] = nid[n]; //gidx[dof*nid[n]+0]/dof;
        }
        
        elcnt++;
      }
    }
  }
  
  if (_nel)  { *_nel = nel; }
  if (_npe)  { *_npe = npe; }
  if (_eidx) { *_eidx = eidx; }
  else { PetscFree(_eidx); }
	
	PetscFunctionReturn(0);
}

/* operations */


/* constructors for 3d */
#undef __FUNCT__
#define __FUNCT__ "DAFE_CreateQ2_3d"
PetscErrorCode DAFE_CreateQ2_3d(DM dm,PetscInt mi,PetscInt mj,PetscInt mk)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  DM             da;
  PetscErrorCode ierr;
  
  fe->nodes_per_el = 27;
  fe->nbasis = 27;
  
	ierr = DMDACreate3d(PetscObjectComm((PetscObject)dm),DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,2*mi+1,2*mj+1,2*mk+1,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,fe->ncomponents,2,NULL,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr); /* force coordinates */
  fe->da = da;
  if (fe->ncomponents > 1) {
    ierr = DMSetMatType(da,MATBAIJ);CHKERRQ(ierr);
  } else {
    ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
  }

  /* sizes */
  ierr = DAFE_GetLocalSizeElementQ2(da,&fe->lmx,&fe->lmy,&fe->lmz);CHKERRQ(ierr);
  fe->lnelements = fe->lmx * fe->lmy * fe->lmz;
  
  /* corners */
  ierr = DAFE_GetCornersElementQ2(da,&fe->corner_imin,&fe->corner_jmin,&fe->corner_kmin,NULL,NULL,NULL);CHKERRQ(ierr);
  
  /* ranges */
  ierr = DAFE_GetOwnershipRangesElementQ2(da,NULL,NULL,NULL,
              &fe->corner_imin_range,&fe->corner_jmin_range,&fe->corner_kmin_range,
              &fe->lmx_range,&fe->lmy_range,&fe->lmz_range);CHKERRQ(ierr);
  
  /* elements -> basis map */
  ierr = DAFE_GetElementsQ2_3D(da,NULL,&fe->lnelements,&fe->element_node_map);CHKERRQ(ierr);
  fe->element_basis_map = fe->element_node_map;

  /* set operations */
  fe->ops->EvaluateBasis      = DAFE_EvaluateBasisQ2_3D;
  fe->ops->GetElementCoords   = DAFE_GetElementCoords_3D;
  fe->ops->GetElementNodeMap  = DAFE_GetElementNodeMap;
  fe->ops->GetElementBasisMap = DAFE_GetElementBasisMap;
  fe->ops->ConvertChildLocalCoordinate = NULL;
  fe->ops->ProjectCoordinates = NULL;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DAFE_GetElementsPk_3D"
PetscErrorCode _DAFE_GetElementsPk_3D(DM dm,PetscInt *_npe,PetscInt *_nel,PetscInt **_eidx)
{
	PetscErrorCode ierr;
	PetscInt mx,my,mz,npe;
	PetscInt nel,d,ei,ej,ek,elcnt,esi,esj,nid[100],width;
	PetscInt *el,*eidx,M,N,P,dof;
	PetscFunctionBegin;
	
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width != 0) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Stencil width must be 0 for Pk with multi-dofs");
	
	npe = 1;
	ierr = DMDAGetInfo(dm,0, 0,0,0, 0,0,0, &dof,0, 0,0,0, 0);CHKERRQ(ierr);
  if (dof > 100) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"dof must be <= 100");
	
  ierr = DMDAGetCorners(dm,&esi,&esj,0,&mx,&my,&mz);CHKERRQ(ierr);
  nel = mx * my * mz;
  ierr = PetscMalloc(sizeof(PetscInt)*(nel*npe*dof+1),&eidx);CHKERRQ(ierr); /* we add one extra space to allow for cases when a proc has zero elements, and we don't want to malloc 0 bytes */
  
  elcnt = 0;
  for (ek=0; ek<mz; ek++) {
    for (ej=0; ej<my; ej++) {
      for (ei=0; ei<mx; ei++) {
        el = &eidx[npe*dof*elcnt];
        
        for (d=0; d<npe*dof; d++) {
          nid[d] = (npe*dof) * elcnt + d;
        }
        
        for (d=0; d<npe*dof; d++) {
          if (nid[d] > M*N*P*dof) {
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of global dofs");
          }
          el[d] = nid[d];
        }
        
        elcnt++;
      }
    }
  }
  
	if (_npe)  *_npe  = npe * dof;
	if (_nel)  *_nel  = nel;
	if (_eidx) *_eidx = eidx;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DAFE_CreatePkFromQ2_3d"
PetscErrorCode _DAFE_CreatePkFromQ2_3d(DM dm,PetscInt ref,PetscInt nbasis)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  DM             da,parent_dm,parent_da;
  PetscInt       MX,MY,MZ,Mp,Np,Pp,k;
  PetscErrorCode ierr;
  
  fe->nbasis = nbasis;
  fe->nodes_per_el = 0;
  
  ierr = DMDAFEGetParentDM(dm,&parent_dm);CHKERRQ(ierr);
  ierr = DMDAFEGetDA(parent_dm,&parent_da);CHKERRQ(ierr);
  ierr = DMDAGetInfo(parent_da,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
  
  MX = fe->mx;
  MY = fe->my;
  MZ = fe->mz;
  
  ierr = DAFE_GetOwnershipRangesElementQ2(parent_da,NULL,NULL,NULL,
                                          NULL,NULL,NULL,
                                          &fe->lmx_range,&fe->lmy_range,&fe->lmz_range);CHKERRQ(ierr);
  /* scale ranges from Q2 */
  for (k=0; k<Mp; k++) { fe->lmx_range[k] = fe->lmx_range[k] * ref; }
  for (k=0; k<Np; k++) { fe->lmy_range[k] = fe->lmy_range[k] * ref; }
  for (k=0; k<Pp; k++) { fe->lmz_range[k] = fe->lmz_range[k] * ref; }
  
	ierr = DMDACreate3d( PetscObjectComm((PetscObject)dm), DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, MX,MY,MZ, Mp,Np,Pp, nbasis,0, fe->lmx_range,fe->lmy_range,fe->lmz_range, &da );CHKERRQ(ierr);
  fe->da = da;
  if (fe->ncomponents > 1) {
    ierr = DMSetMatType(da,MATBAIJ);CHKERRQ(ierr);
  } else {
    ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
  }
  
  /* elements */
  ierr = _DAFE_GetElementsPk_3D(da,NULL,NULL,&fe->element_basis_map);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAFE_CreateP0FromQ2_3d"
PetscErrorCode DAFE_CreateP0FromQ2_3d(DM dm,PetscInt ref)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = _DAFE_CreatePkFromQ2_3d(dm,ref,1);CHKERRQ(ierr);
  
  /* set operations */
  fe->ops->EvaluateBasis      = NULL; /* todo */
  fe->ops->GetElementCoords   = NULL;
  fe->ops->GetElementNodeMap  = NULL;
  fe->ops->GetElementBasisMap = DAFE_GetElementBasisMap;
  fe->ops->GetParentElementIndex = GetParentElementIndex_Q1sub2Q2;
  
  if (ref == 1) {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Ref1;
  } else if (ref == 2) {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Ref2;
  } else {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Q1sub2Q2;
  }
  fe->ops->ProjectCoordinates = NULL;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAFE_CreateP1FromQ2_3d"
PetscErrorCode DAFE_CreateP1FromQ2_3d(DM dm,PetscInt ref)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = _DAFE_CreatePkFromQ2_3d(dm,ref,4);CHKERRQ(ierr);

  /* set operations */
  fe->ops->EvaluateBasis      = NULL; /* todo */
  fe->ops->GetElementCoords   = NULL;
  fe->ops->GetElementNodeMap  = NULL;
  fe->ops->GetElementBasisMap = DAFE_GetElementBasisMap;
  fe->ops->GetParentElementIndex = GetParentElementIndex_Q1sub2Q2;
  
  if (ref == 1) {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Ref1;
  } else if (ref == 2) {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Ref2;
  } else {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Q1sub2Q2;
  }
  fe->ops->ProjectCoordinates = NULL;
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DAFE_GetSubElementsQ1_3D"
PetscErrorCode DAFE_GetSubElementsQ1_3D(DM dafe,PetscInt *_npe,PetscInt *_nel,PetscInt **_eidx)
{
  DM_DAFE        *fe = (DM_DAFE*)dafe->data;
  DM             dm;
	const PetscInt order = 1;
	PetscErrorCode ierr;
	PetscInt mx,my,mz,npe, M,N,P;
	PetscInt nel,ei,ej,ek,i,j,k,elcnt,esi,esj,esk,gsi,gsj,gsk,nid[8],n,X,Y,Z,width;
	PetscInt *el,*eidx;
	PetscFunctionBegin;
	
  ierr = DMDAFEGetDA(dafe,&dm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dm,0, &M,&N,&P, 0,0,0, 0,&width, 0,0,0, 0);CHKERRQ(ierr);
	if (width != 1) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Stencil width must be 1 for Q1");
	
	npe = (order + 1)*(order + 1)*(order + 1);

  mx = fe->lmx;
  my = fe->lmy;
  mz = fe->lmz;
  nel = mx * my * mz;
  
  esi = fe->corner_imin;
  esj = fe->corner_jmin;
  esk = fe->corner_kmin;
  
  ierr = PetscMalloc(sizeof(PetscInt)*(nel*npe+1),&eidx);CHKERRQ(ierr);
  ierr = DMDAGetGhostCorners(dm,&gsi,&gsj,&gsk, &X,&Y,&Z);CHKERRQ(ierr);
  
  elcnt = 0;
  for (ek=0; ek<mz; ek++) {
    k = (esk-gsk) + ek;
    
    for (ej=0; ej<my; ej++) {
      j = (esj-gsj) + ej;
      
      for (ei=0; ei<mx; ei++) {
        i = (esi-gsi) + ei;
        
        el = &eidx[npe*elcnt];
        
        nid[0] = (i  ) + (j  ) *X  + (k  ) *X*Y;
        nid[1] = (i+1) + (j  ) *X  + (k  ) *X*Y;
        
        nid[2] = (i  ) + (j+1) *X  + (k  ) *X*Y;
        nid[3] = (i+1) + (j+1) *X  + (k  ) *X*Y;
        
        nid[4] = (i  ) + (j  ) *X  + (k+1) *X*Y;
        nid[5] = (i+1) + (j  ) *X  + (k+1) *X*Y;
        
        nid[6] = (i  ) + (j+1) *X  + (k+1) *X*Y;
        nid[7] = (i+1) + (j+1) *X  + (k+1) *X*Y;
        
        //printf("e=%d: [%d %d %d %d] [%d %d %d %d] \n",elcnt,nid[0],nid[1],nid[2],nid[3],nid[0+4],nid[1+4],nid[2+4],nid[3+4] );
        
        for (n=0; n<npe; n++) {
          if (nid[n] > M*N*P) {
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of global nodes");
          }
          
          if (nid[n] > X*Y*Z) {
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Local indexing exceeds number of local nodes");
          }
          
          el[n] = nid[n]; //gidx[dof*nid[n]+0]/dof;
        }
        
        elcnt++;
      }
    }
  }

	if (_npe)  *_npe = npe;
	if (_nel)  *_nel = nel;
	if (_eidx) *_eidx = eidx;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAFE_CreateQ1FromQ2_3d"
PetscErrorCode DAFE_CreateQ1FromQ2_3d(DM dm,PetscInt ref)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  DM             da,parent,daq2;
  PetscInt       i,j,k,Mp,Np,Pp,*lsip,*lsjp,*lskp,*lmx,*lmy,*lmz,sei,sej,sek;
  PetscInt       NX,NY,NZ,*NXp,*NYp,*NZp;
  PetscErrorCode ierr;
  
  
  if (ref > 1) {
    if (ref%2 != 0) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Refinement factor must be an even number");
  }
  
  ierr = DMDAFEGetParentDM(dm,&parent);CHKERRQ(ierr);
  ierr = DMDAFEGetDA(parent,&daq2);CHKERRQ(ierr);
  
  /* sizes */
  fe->nodes_per_el = 8;
  fe->nbasis = 8;
  
  
	ierr = DMDAGetInfo(daq2,0, 0,0,0, &Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);
  ierr = DAFE_GetCornersElementQ2(daq2,&sei,&sej,&sek,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DAFE_GetOwnershipRangesElementQ2(daq2,0,0,0,&lsip,&lsjp,&lskp,&lmx,&lmy,&lmz);CHKERRQ(ierr);

  /* corners */
  sei = sei/2;
  sej = sej/2;
  sek = sek/2;

  sei = sei * ref;
  sej = sej * ref;
  sek = sek * ref;
  
  fe->corner_imin = sei;
  fe->corner_jmin = sej;
  fe->corner_kmin = sek;
  
  /* ranges */
  /* reset to overlapping Q1 */
	for (i=0; i<Mp; i++) { lsip[i] = lsip[i]/2; }
	for (j=0; j<Np; j++) { lsjp[j] = lsjp[j]/2; }
	for (k=0; k<Pp; k++) { lskp[k] = lskp[k]/2; }
  
  /* scale by refinement factor */
  for (i=0; i<Mp; i++) { lsip[i] = lsip[i] * ref; }
	for (j=0; j<Np; j++) { lsjp[j] = lsjp[j] * ref; }
	for (k=0; k<Pp; k++) { lskp[k] = lskp[k] * ref; }

  fe->corner_imin_range = lsip;
  fe->corner_jmin_range = lsjp;
  fe->corner_kmin_range = lskp;
  
  /* scale local element count by refinement factor */
  for (i=0; i<Mp; i++) { lmx[i] = lmx[i] * ref; }
	for (j=0; j<Np; j++) { lmy[j] = lmy[j] * ref; }
	for (k=0; k<Pp; k++) { lmz[k] = lmz[k] * ref; }

  fe->lmx_range = lmx;
  fe->lmy_range = lmy;
  fe->lmz_range = lmz;
  
  ierr = PetscMalloc(sizeof(PetscInt)*Mp,&NXp);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Np,&NYp);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*Pp,&NZp);CHKERRQ(ierr);

  /* points per rank */
  for (i=0; i<Mp; i++) { NXp[i] = lmx[i]; }
  NXp[Mp-1]++;

	for (j=0; j<Np; j++) { NYp[j] = lmy[j]; }
  NYp[Np-1]++;

	for (k=0; k<Pp; k++) { NZp[k] = lmz[k]; }
  NZp[Pp-1]++;

  
  /* total number of points */
  NX = fe->mx + 1;
  NY = fe->my + 1;
  NZ = fe->mz + 1;
  
  ierr = DMDACreate3d( PetscObjectComm((PetscObject)dm), DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                      NX,NY,NZ, Mp,Np,Pp,
                      fe->ncomponents,1,
                      NXp,NYp,NZp, &da );CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr); /* force coordinates */
  fe->da = da;
  if (fe->ncomponents > 1) {
    ierr = DMSetMatType(da,MATBAIJ);CHKERRQ(ierr);
  } else {
    ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
  }

  /* elements -> basis map */
  ierr = DAFE_GetSubElementsQ1_3D(dm,NULL,NULL,&fe->element_node_map);CHKERRQ(ierr);
  fe->element_basis_map = fe->element_node_map;
  
  /* set operations */
  fe->ops->EvaluateBasis      = NULL; /* todo */
  fe->ops->GetElementCoords   = DAFE_GetElementCoords_3D;
  fe->ops->GetElementNodeMap  = DAFE_GetElementNodeMap;
  fe->ops->GetElementBasisMap = DAFE_GetElementBasisMap;
  fe->ops->GetParentElementIndex = GetParentElementIndex_Q1sub2Q2;

  if (ref == 1) {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Ref1;
  } else if (ref == 2) {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Ref2;
  } else {
    fe->ops->ConvertChildLocalCoordinate = ConvertChildLocalCoordinateToParentElement_Q1sub2Q2;
  }

  fe->ops->ProjectCoordinates = _DAFEProjectCoordinates_Q1ToQ2;
  
  PetscFree(NZp);
  PetscFree(NYp);
  PetscFree(NXp);
  
  PetscFunctionReturn(0);
}

/*
 Two special cases:
 ref = 1 (lmx matches)
 ref = 2 (q2.lmx = 2 x q1.lmx)
*/
#undef __FUNCT__
#define __FUNCT__ "_DAFEProjectCoordinates_Q1ToQ2"
PetscErrorCode _DAFEProjectCoordinates_Q1ToQ2(DM dafeq2,DM dafeq1)
{
  DM_DAFE        *feq2 = (DM_DAFE*)dafeq2->data;
  DM_DAFE        *feq1 = (DM_DAFE*)dafeq1->data;
  DM daq2,daq1;
  PetscInt ref;
  PetscErrorCode ierr;
  
  ierr = DMDAFEGetDA(dafeq2,&daq2);CHKERRQ(ierr);
  ierr = DMDAFEGetDA(dafeq1,&daq1);CHKERRQ(ierr);
  
  ref = feq1->mx / feq2->mx;
  
  if (ref <= 0) {
    SETERRQ(PetscObjectComm((PetscObject)dafeq2),PETSC_ERR_USER,"Refinement factor ratio must be >= 1");
  } else if (ref == 1) {
    Vec coordinates[2];
    PetscInt i,j,k,si,sj,sk,m,n,p,si2,sj2,sk2,m2,n2,p2;
    const PetscScalar *LA_c0;
    PetscScalar *LA_c1;
    
    ierr = DMGetCoordinatesLocal(daq2,&coordinates[0]);CHKERRQ(ierr);
    ierr = DMGetCoordinates(daq1,&coordinates[1]);CHKERRQ(ierr);

    ierr = VecGetArrayRead(coordinates[0],&LA_c0);CHKERRQ(ierr);
    ierr = VecGetArray(coordinates[1],&LA_c1);CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners(daq2,&si2,&sj2,&sk2,&m2,&n2,&p2);CHKERRQ(ierr);
    //printf("[q2] %d %d %d \n",m2,n2,p2);
    ierr = DMDAGetCorners(daq1,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
    //printf("[q1] %d %d %d \n",m,n,p);
    
    for (k=0; k<p; k++) {
      for (j=0; j<n; j++) {
        for (i=0; i<m; i++) {
          PetscInt nid,nid2,i2,j2,k2;
          PetscInt q1i,q1j,q1k;
          PetscInt q2i,q2j,q2k;
          
          q1i = si + i;
          q1j = sj + j;
          q1k = sk + k;
          
          q2i = 2*q1i;
          q2j = 2*q1j;
          q2k = 2*q1k;
          
          i2 = q2i - si2;
          j2 = q2j - sj2;
          k2 = q2k - sk2;
          
          nid = i + j*m + k*m*n;
          nid2 = i2 + (j2)*m2 + (k2)*m2*n2;
          
          //printf("[Q1] i,j,k %d %d %d : [Q2] %d %d %d \n",i,j,k,2*i,2*j,2*k);
          LA_c1[3*nid]   = LA_c0[3*nid2];
          LA_c1[3*nid+1] = LA_c0[3*nid2+1];
          LA_c1[3*nid+2] = LA_c0[3*nid2+2];
        }
      }
    }
    ierr = VecRestoreArray(coordinates[1],&LA_c1);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(coordinates[0],&LA_c0);CHKERRQ(ierr);
    
    ierr = DMDAUpdateGhostedCoordinates(daq1);CHKERRQ(ierr);

  } else if (ref == 2) {
    Vec coordinates[2],coordinates_local[2];
    
    ierr = DMGetCoordinates(daq2,&coordinates[0]);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(daq2,&coordinates_local[0]);CHKERRQ(ierr);

    ierr = DMGetCoordinates(daq1,&coordinates[1]);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(daq1,&coordinates_local[1]);CHKERRQ(ierr);

    ierr = VecCopy(coordinates[0],coordinates[1]);CHKERRQ(ierr); /* q1 <- q2 */
    ierr = VecCopy(coordinates_local[0],coordinates_local[1]);CHKERRQ(ierr);
    
  } else { /* general case - interpolate */
    PetscInt e,ep,i,k;
    PetscReal xic[3],xi[3];
    Vec coorQ1,coorQ2;
    PetscScalar *LA_coorQ1;
    const PetscScalar *LA_coorQ2;
    
    ierr = DMGetCoordinatesLocal(daq2,&coorQ2);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coorQ2,&LA_coorQ2);CHKERRQ(ierr);

    ierr = DMGetCoordinatesLocal(daq1,&coorQ1);CHKERRQ(ierr);
    ierr = VecGetArray(coorQ1,&LA_coorQ1);CHKERRQ(ierr);

    for (e=0; e<feq1->lnelements; e++) {
      const PetscInt *element_q1;
      const PetscInt *element_q2;
      PetscReal elcoorQ2[3*27],Ni[27],xq1[3];
      PetscReal xic_xi[]   = { -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, 1.0 };
      PetscReal xic_eta[]  = { -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0 };
      PetscReal xic_zeta[] = { -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0, 1.0 };
      
      ierr = feq1->ops->GetParentElementIndex(dafeq1,e,&ep);CHKERRQ(ierr);

      ierr = feq2->ops->GetElementNodeMap(dafeq2,ep,&element_q2);CHKERRQ(ierr);

      //ierr = feq2->ops->GetElementCoords(dafeq2,ep,elcoorQ2);CHKERRQ(ierr);
      for (i=0; i<27; i++) {
        elcoorQ2[3*i+0] = LA_coorQ2[3*element_q2[i]+0];
        elcoorQ2[3*i+1] = LA_coorQ2[3*element_q2[i]+1];
        elcoorQ2[3*i+2] = LA_coorQ2[3*element_q2[i]+2];
      }
      
      //printf("sub[%d] --> parent[%d]\n",e,ep);
      
      ierr = feq1->ops->GetElementNodeMap(dafeq1,e,&element_q1);CHKERRQ(ierr);

      for (k=0; k<8; k++) {
        xic[0] = xic_xi[k];
        xic[1] = xic_eta[k];
        xic[2] = xic_zeta[k];
        ierr = feq1->ops->ConvertChildLocalCoordinate(dafeq1,e,xic,xi);CHKERRQ(ierr);
        //printf("%1.4e %1.4e %1.4e\n",xi[0],xi[1],xi[2]);
        //xi[0] = xi[1] = xi[2] = 0.0;
        
        ierr = feq2->ops->EvaluateBasis(NULL,xi,Ni);
        xq1[0]= xq1[1] = xq1[2] = 0.0;
        for (i=0; i<feq2->nodes_per_el; i++) {
          xq1[0] += Ni[i] * elcoorQ2[3*i];
          xq1[1] += Ni[i] * elcoorQ2[3*i+1];
          xq1[2] += Ni[i] * elcoorQ2[3*i+2];
        }
        LA_coorQ1[3*element_q1[k]+0] = xq1[0];
        LA_coorQ1[3*element_q1[k]+1] = xq1[1];
        LA_coorQ1[3*element_q1[k]+2] = xq1[2];
      }
    
    }
    
    ierr = VecRestoreArray(coorQ1,&LA_coorQ1);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(coorQ2,&LA_coorQ2);CHKERRQ(ierr);
    
    ierr = DMDASetCoordinatesFromLocalVector(daq1,coorQ1);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

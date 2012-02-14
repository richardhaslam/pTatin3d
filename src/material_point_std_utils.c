

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "pTatin3d.h"
#include "element_type_Q2.h"
#include "dmda_element_q2p1.h"
#include "swarm_fields.h"
#include "data_exchanger.h"
#include "MPntStd_def.h"
#include "ptatin3d_stokes.h"
#include "output_paraview.h"


#undef __FUNCT__
#define __FUNCT__ "SwarmMPntStd_AssignUniquePointIdentifiers"
PetscErrorCode SwarmMPntStd_AssignUniquePointIdentifiers(MPI_Comm comm,DataBucket db,int start_pid,int end_pid)
{
	DataField    PField;
	long int     np_local, np_global, max_local, max;
	int          rank,p,L;
	
	
	PetscFunctionBegin;
	
	MPI_Comm_rank(comm,&rank);
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	
	DataBucketGetSizes(db,&L,0,0);
	
	/* find max pid presently in the system */
	max_local = 0;
	for (p=0; p<L; p++) {
		MPntStd *marker;
		DataFieldAccessPoint(PField,p,(void**)&marker);
		
		if ( marker->pid > max_local ) {
			max_local = marker->pid;
		}
	}
	MPI_Allreduce( &max_local, &max, 1, MPI_LONG, MPI_MAX, comm );
	PetscPrintf(PETSC_COMM_WORLD,"SwarmMPntStd_AssignUniquePointIdentifiers : max_pid = %ld \n",max);
	max = max + 1;
	
	/* give particles a unique identifier */
	np_local = (end_pid-start_pid);

	MPI_Scan( &np_local, &np_global, 1, MPI_LONG, MPI_SUM, comm );
	//printf("rank %d : np_local = %ld, np_global = %ld \n",rank,np_local,np_global);
	for (p=start_pid; p<end_pid; p++) {
		MPntStd *marker;
		DataFieldAccessPoint(PField,p,(void**)&marker);
		
		//marker->pid = max + (np_global-1) - (np_local-1-p);
		marker->pid = max + (np_global-np_local) + (p-start_pid);
		//printf("assigning %d -> pid = %ld \n", p, marker->pid );
	}
	
	
	DataFieldRestoreAccess(PField);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmMPntStd_CoordAssignment_LatticeLayout3d"
PetscErrorCode SwarmMPntStd_CoordAssignment_LatticeLayout3d(DM da,PetscInt Nxp[],PetscReal perturb,DataBucket db)
{
	DataField    PField;
	PetscInt     e,mE,nE,pE;
  DM           cda;
  Vec          gcoords;
  PetscScalar  *LA_coords;
  PetscScalar  el_coords[Q2_NODES_PER_EL_3D*NSD];
	int          ncells,np_per_cell;
	PetscInt     nel,nen;
	const PetscInt     *elnidx;
	PetscInt     p,k,pi,pj,pk;
	PetscReal    dxi,deta,dzeta;
	long int     np_local, np_global;
	int          rank;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	PetscOptionsGetReal(PETSC_NULL,"-lattice_layout_perturb", &perturb, PETSC_NULL );
	PetscOptionsGetInt(PETSC_NULL,"-lattice_layout_Nx", &Nxp[0], PETSC_NULL );
	PetscOptionsGetInt(PETSC_NULL,"-lattice_layout_Ny", &Nxp[1], PETSC_NULL );
	PetscOptionsGetInt(PETSC_NULL,"-lattice_layout_Nz", &Nxp[2], PETSC_NULL );
	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	// re-size //
	ncells = nel;
	np_per_cell = Nxp[0] * Nxp[1] * Nxp[2];
	DataBucketSetSizes(db,np_per_cell*ncells,-1);
	
	if (perturb<0.0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cannot use a negative perturbation");
	}
	if (perturb>1.0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cannot use a perturbation greater than 1.0");
	}
	
  /* setup for coords */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	
	dxi   = 2.0/(PetscReal)Nxp[0];
	deta  = 2.0/(PetscReal)Nxp[1];
	dzeta  = 2.0/(PetscReal)Nxp[2];
	
	p = 0;
	for (e = 0; e < ncells; e++) {
		/* get coords for the element */
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_coords);CHKERRQ(ierr);
		
		for (pk=0; pk<Nxp[2]; pk++) {
			for (pj=0; pj<Nxp[1]; pj++) {
				for (pi=0; pi<Nxp[0]; pi++) {
					MPntStd *marker;
					double xip_rand[NSD],xp_rand[NSD],Ni[Q2_NODES_PER_EL_3D];
					
					/* random between -1 <= xi,eta,zeta <= 1 */
					xip_rand[0] = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
					xip_rand[1] = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
					xip_rand[2] = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
					
					xip_rand[0] = perturb * dxi    * xip_rand[0];
					xip_rand[1] = perturb * deta   * xip_rand[1];
					xip_rand[2] = perturb * dzeta  * xip_rand[2];
					
					xip_rand[0] += -1.0 + dxi    * (pi + 0.5);
					xip_rand[1] += -1.0 + deta   * (pj + 0.5);
					xip_rand[2] += -1.0 + dzeta  * (pk + 0.5);
					
					if (fabs(xip_rand[0]) > 1.0) {
						SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"fabs(x-point coord) greater than 1.0");
					}
					if (fabs(xip_rand[1]) > 1.0) {
						SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"fabs(y-point coord) greater than 1.0");
					}
					if (fabs(xip_rand[2]) > 1.0) {
						SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"fabs(z-point coord) greater than 1.0");
					}
					
					pTatin_ConstructNi_Q2_3D(xip_rand,Ni);
					
					xp_rand[0] = xp_rand[1] = xp_rand[2] = 0.0;
					for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
						xp_rand[0] += Ni[k] * el_coords[NSD*k+0];
						xp_rand[1] += Ni[k] * el_coords[NSD*k+1];
						xp_rand[2] += Ni[k] * el_coords[NSD*k+2];
					}
					
					DataFieldAccessPoint(PField,p,(void**)&marker);
					
					marker->coor[0] = xp_rand[0];
					marker->coor[1] = xp_rand[1];
					marker->coor[2] = xp_rand[2];
					
					marker->xi[0] = xip_rand[0];
					marker->xi[1] = xip_rand[1];
					marker->xi[2] = xip_rand[2];
					
					marker->wil    = e;
					marker->pid    = 0;
					
					p++;
				}
			}
		}		
	}
	
	DataFieldRestoreAccess(PField);
	
	np_local = np_per_cell * ncells;
	ierr = SwarmMPntStd_AssignUniquePointIdentifiers(((PetscObject)da)->comm,db,0,np_local);CHKERRQ(ierr);
	
  ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmMPntStd_CoordAssignment_RandomLayout3d"
PetscErrorCode SwarmMPntStd_CoordAssignment_RandomLayout3d(DM da,PetscInt nPerCell,DataBucket db)
{
	DataField    PField;
	PetscInt     e,mE,nE,pE;
  DM           cda;
  Vec          gcoords;
  PetscScalar  *LA_coords;
  PetscScalar  el_coords[Q2_NODES_PER_EL_3D*NSD];
	int          ncells,np_per_cell;
	PetscInt     nel,nen;
	const PetscInt *elnidx;
	PetscInt     p,k,pi;
	long int     np_local, np_global;
	int          rank;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	PetscOptionsGetInt(PETSC_NULL,"-random_layout_Np", &nPerCell, PETSC_NULL );
	
	// re-size //
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen,&elnidx);CHKERRQ(ierr);
	ncells = nel;
	np_per_cell = nPerCell;
	DataBucketSetSizes(db,np_per_cell*ncells,-1);
	
  /* setup for coords */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	
	p = 0;
	for (e = 0; e < ncells; e++) {
		/* get coords for the element */
		ierr = DMDAGetElementCoordinatesQ2_3D(el_coords,(PetscInt*)&elnidx[nen*e],LA_coords);CHKERRQ(ierr);
		
		for (pi=0; pi<np_per_cell; pi++) {
			MPntStd *marker;
			double xip_rand[NSD],xp_rand[NSD],Ni[Q2_NODES_PER_EL_3D];
			
			/* random between -1 <= xi,eta,zeta <= 1 */
			xip_rand[0] = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
			xip_rand[1] = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
			xip_rand[2] = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
			
			pTatin_ConstructNi_Q2_3D(xip_rand,Ni);
			
			xp_rand[0] = xp_rand[1] = xp_rand[2] = 0.0;
			for (k=0; k<Q2_NODES_PER_EL_3D; k++) {
				xp_rand[0] += Ni[k] * el_coords[NSD*k+0];
				xp_rand[1] += Ni[k] * el_coords[NSD*k+1];
				xp_rand[2] += Ni[k] * el_coords[NSD*k+2];
			}
			
			DataFieldAccessPoint(PField,p,(void**)&marker);
			
			marker->coor[0] = xp_rand[0];
			marker->coor[1] = xp_rand[1];
			marker->coor[2] = xp_rand[2];
			
			marker->xi[0] = xip_rand[0];
			marker->xi[1] = xip_rand[1];
			marker->xi[2] = xip_rand[2];
			
			marker->wil    = e;
			marker->pid    = p;
			
			p++;
		}
		
	}
	
	np_local = np_per_cell * ncells;
	ierr = SwarmMPntStd_AssignUniquePointIdentifiers(((PetscObject)da)->comm,db,0,np_local);CHKERRQ(ierr);
	
	DataFieldRestoreAccess(PField);
	ierr = DMDAVecRestoreArray(cda,gcoords,&LA_coords);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmView_MPntStd_VTKascii"
PetscErrorCode SwarmView_MPntStd_VTKascii(DataBucket db,const char name[])
{
	PetscMPIInt rank;
	FILE *vtk_fp;
	PetscInt k;
	int npoints;
	PetscLogDouble t0,t1;
	DataField PField;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscGetTime(&t0);CHKERRQ(ierr);
	
	if ((vtk_fp = fopen ( name, "w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
	}
	
	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
	
#ifdef WORDSIZE_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	fprintf( vtk_fp, "\t<UnstructuredGrid>\n" );
	
	DataBucketGetSizes(db,&npoints,PETSC_NULL,PETSC_NULL);
	fprintf( vtk_fp, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",npoints,npoints );
	
	
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<Cells>\n");
	
	// connectivity //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	fprintf( vtk_fp, "\t\t\t\t");
	for(k=0;k<npoints;k++) {
		fprintf( vtk_fp,"%d ",k);
	}
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");	
	
	// offsets //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	fprintf( vtk_fp, "\t\t\t\t");
	for(k=0;k<npoints;k++) {
		fprintf( vtk_fp,"%d ",k+1);
	}
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	
	// types //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
	fprintf( vtk_fp, "\t\t\t\t");
	for(k=0;k<npoints;k++) {
		fprintf( vtk_fp,"1 "); // VTK_VERTEX //
	}
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	
	fprintf( vtk_fp, "\t\t\t</Cells>\n");
	
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<CellData>\n");
	fprintf( vtk_fp, "\t\t\t</CellData>\n");
	fprintf( vtk_fp, "\n");
	
	
	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	
	
	/* point coordinates */
	fprintf( vtk_fp, "\t\t\t<Points>\n");
	
	/* copy coordinates */
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for(k=0;k<npoints;k++) {
		MPntStd *marker;
		double *coords;
		
		DataFieldAccessPoint(PField,k,(void**)&marker);
		
		
		/* extract coords from your data type */
		//coords = elasticParticle->pos;
		MPntStdGetField_global_coord( marker,&coords );
		
		fprintf( vtk_fp,"\t\t\t\t\t%lf %lf %lf \n",coords[0],coords[1],coords[2]);
	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	
	fprintf( vtk_fp, "\t\t\t</Points>\n");
	fprintf( vtk_fp, "\n");
	
	DataFieldRestoreAccess(PField);
	
	/* point data BEGIN */
	fprintf( vtk_fp, "\t\t\t<PointData>\n");
	
	/* auto generated shit goes here */
	{
		MPntStd *marker = PField->data; /* should write a function to do this */
		
		MPntStdVTKWriteAsciiAllFields(vtk_fp,(const int)npoints,(const MPntStd*)marker );
	}
	fprintf( vtk_fp, "\t\t\t</PointData>\n");
	fprintf( vtk_fp, "\n");
	/* point data END */
	
	
	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</UnstructuredGrid>\n");
	
	fprintf( vtk_fp, "</VTKFile>\n");
	
	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}
	
	ierr = PetscGetTime(&t1);CHKERRQ(ierr);
#ifdef PROFILE_TIMING
	PetscPrintf(PETSC_COMM_WORLD,"VTKWriter(%s): Time %1.4e sec\n",__FUNCT__,t1-t0);
#endif
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmView_MPntStd_VTKappended_binary"
PetscErrorCode SwarmView_MPntStd_VTKappended_binary(DataBucket db,const char name[])
{
	PetscMPIInt rank;
	FILE *vtk_fp;
	PetscInt k;
	int npoints;
	PetscLogDouble t0,t1;
	DataField PField;
	int byte_offset,length;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscGetTime(&t0);CHKERRQ(ierr);
	
	if ((vtk_fp = fopen ( name, "w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
	}

	DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField);
	
	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
	
#ifdef WORDSIZE_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	fprintf( vtk_fp, "\t<UnstructuredGrid>\n" );
	
	DataBucketGetSizes(db,&npoints,PETSC_NULL,PETSC_NULL);
	fprintf( vtk_fp, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",npoints,npoints );
	
	
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<Cells>\n");
	
	byte_offset = 0;
	
	// connectivity //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
  byte_offset = byte_offset + sizeof(int) + npoints * sizeof(int);
	
	// offsets //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
  byte_offset = byte_offset + sizeof(int) + npoints * sizeof(int);
	
	// types //
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
  byte_offset = byte_offset + sizeof(int) + npoints * sizeof(unsigned char);
	
	fprintf( vtk_fp, "\t\t\t</Cells>\n");
	
	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<CellData>\n");
	fprintf( vtk_fp, "\t\t\t</CellData>\n");
	fprintf( vtk_fp, "\n");
	
	fprintf( vtk_fp, "\t\t\t<Points>\n");
	
	/* coordinates */
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",byte_offset);
  byte_offset = byte_offset + sizeof(int) + npoints * 3 * sizeof(double);
	
	fprintf( vtk_fp, "\t\t\t</Points>\n");
	fprintf( vtk_fp, "\n");
	
	/* point data BEGIN */
	fprintf( vtk_fp, "\t\t\t<PointData>\n");
	/* auto generated shit for the header goes here */
	{
		MPntStd *marker = PField->data; /* should write a function to do this */

		MPntStdVTKWriteBinaryAppendedHeaderAllFields(vtk_fp,&byte_offset,npoints,marker);
	}
	fprintf( vtk_fp, "\t\t\t</PointData>\n");
	fprintf( vtk_fp, "\n");
	/* point data END */

	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</UnstructuredGrid>\n");

	/* WRITE APPENDED DATA HERE */
	fprintf( vtk_fp,"\t<AppendedData encoding=\"raw\">\n");
	fprintf( vtk_fp,"_");

	/* connectivity, offsets, types, coords */
	////////////////////////////////////////////////////////
	/* write connectivity */
	length = sizeof(int)*npoints;
	fwrite( &length,sizeof(int),1,vtk_fp);
	for (k=0; k<npoints; k++) {
		int idx = k;
		fwrite( &idx, sizeof(int),1, vtk_fp );
	}
	////////////////////////////////////////////////////////
	/* write offset */
	length = sizeof(int)*npoints;
	fwrite( &length,sizeof(int),1,vtk_fp);
	for (k=0; k<npoints; k++) {
		int idx = k+1;
		fwrite( &idx, sizeof(int),1, vtk_fp );
	}
	////////////////////////////////////////////////////////
	/* write types */
	length = sizeof(unsigned char)*npoints;
	fwrite( &length,sizeof(int),1,vtk_fp);
	for (k=0; k<npoints; k++) {
		unsigned char idx = 1; /* VTK_VERTEX */
		fwrite( &idx, sizeof(unsigned char),1, vtk_fp );
	}
	////////////////////////////////////////////////////////
	/* write coordinates */
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));

	length = sizeof(double)*npoints*3;
	fwrite( &length,sizeof(int),1,vtk_fp);
	for (k=0; k<npoints; k++) {
		MPntStd *marker;
		double  *coor;
		double  coords_k[] = {0.0, 0.0, 0.0};
		
		DataFieldAccessPoint(PField,k,(void**)&marker);
		MPntStdGetField_global_coord(marker,&coor);
		coords_k[0] = coor[0];
		coords_k[1] = coor[1];
		coords_k[2] = coor[2];
		
		fwrite( coords_k, sizeof(double), 3, vtk_fp );
	}
	DataFieldRestoreAccess(PField);
	
	/* auto generated shit for the marker data goes here */
	{
		MPntStd *marker = PField->data;
		MPntStdVTKWriteBinaryAppendedDataAllFields(vtk_fp,npoints,marker);
	}
		
	fprintf( vtk_fp,"\n\t</AppendedData>\n");
	
	fprintf( vtk_fp, "</VTKFile>\n");
	
	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}
	
	ierr = PetscGetTime(&t1);CHKERRQ(ierr);
#ifdef PROFILE_TIMING
	PetscPrintf(PETSC_COMM_WORLD,"VTKWriter(%s): Time %1.4e sec\n",__FUNCT__,t1-t0);
#endif
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "__SwarmView_MPntStd_PVTU"
PetscErrorCode __SwarmView_MPntStd_PVTU(const char prefix[],const char name[])
{
	PetscMPIInt nproc,rank;
	FILE *vtk_fp;
	PetscInt i;
	char *sourcename;
	
	PetscFunctionBegin;
	
	if ((vtk_fp = fopen ( name, "w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
	}
	
	/* (VTK) generate pvts header */
	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
	
#ifdef WORDSIZE_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	/* define size of the nodal mesh based on the cell DM */
	fprintf( vtk_fp, "  <PUnstructuredGrid GhostLevel=\"0\">\n" ); /* note overlap = 0 */
	
	/* DUMP THE CELL REFERENCES */
	fprintf( vtk_fp, "    <PCellData>\n");
	fprintf( vtk_fp, "    </PCellData>\n");
	
	///////////////
	fprintf( vtk_fp, "    <PPoints>\n");
	fprintf( vtk_fp, "      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	fprintf( vtk_fp, "    </PPoints>\n");
	///////////////
	
	///////////////
  fprintf(vtk_fp, "    <PPointData>\n");
	MPntStdPVTUWriteAllPPointDataFields(vtk_fp);
  fprintf(vtk_fp, "    </PPointData>\n");
	///////////////
	
	/* write out the parallel information */
	MPI_Comm_size(PETSC_COMM_WORLD,&nproc);
	for (i=0; i<nproc; i++) {
		asprintf( &sourcename, "%s-subdomain%1.5d.vtu", prefix, i );
		fprintf( vtk_fp, "    <Piece Source=\"%s\"/>\n",sourcename);
		free(sourcename);
	}
	
	/* close the file */
	fprintf( vtk_fp, "  </PUnstructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");
	
	if(vtk_fp!=NULL){
		fclose( vtk_fp );
		vtk_fp = NULL;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmOutputParaView_MPntStd"
PetscErrorCode SwarmOutputParaView_MPntStd(DataBucket db,const char path[],const char prefix[])
{ 
	char *vtkfilename,*filename;
	PetscMPIInt rank;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = pTatinGenerateParallelVTKName(prefix,"vtu",&vtkfilename);CHKERRQ(ierr);
	if (path) {
		asprintf(&filename,"%s/%s",path,vtkfilename);
	} else {
		asprintf(&filename,"./%s",vtkfilename);
	}

//#ifdef __VTK_ASCII__
//	ierr = SwarmView_MPntStd_VTKascii( db,filename );CHKERRQ(ierr);
//#endif
//#ifndef __VTK_ASCII__
	ierr = SwarmView_MPntStd_VTKappended_binary(db,filename);CHKERRQ(ierr);
//#endif
	free(filename);
	free(vtkfilename);
	
	ierr = pTatinGenerateVTKName(prefix,"pvtu",&vtkfilename);CHKERRQ(ierr);
	if (path) {
		asprintf(&filename,"%s/%s",path,vtkfilename);
	} else {
		asprintf(&filename,"./%s",vtkfilename);
	}
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	if (rank==0) {
		ierr = __SwarmView_MPntStd_PVTU( prefix, filename );CHKERRQ(ierr);
	}
	free(filename);
	free(vtkfilename);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_MPntStd_Euler"
PetscErrorCode SwarmUpdatePosition_MPntStd_Euler(DM da,Vec velocity,PetscReal step,int npoints,MPntStd marker[])
{
	Vec             Lvelocity;
	PetscScalar     *LA_velocity;
	PetscScalar     el_velocity[Q2_NODES_PER_EL_3D*NSD];
	PetscInt        e,i,p,wil;
	PetscScalar     Ni_p[Q2_NODES_PER_EL_3D],vel_p[NSD];
	PetscInt nel,nen_u;
	const PetscInt *elnidx_u;
	PetscInt vel_el_lidx[U_BASIS_FUNCTIONS*3];
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* advect */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	for (p=0; p<npoints; p++) {
		MPntStd *marker_p = &marker[p];
		
		wil   = marker_p->wil;
		e     = wil;
		if (wil<0) { SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Point[%d] has wil_e < 0", wil ); }
		
		ierr = StokesVelocity_GetElementLocalIndices(vel_el_lidx,(PetscInt*)&elnidx_u[nen_u*e]);CHKERRQ(ierr);
		ierr = DMDAGetVectorElementFieldQ2_3D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		pTatin_ConstructNi_Q2_3D(marker_p->xi,Ni_p);
		
		vel_p[0] = vel_p[1] = vel_p[2] = 0.0;
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			vel_p[0] += Ni_p[i] * el_velocity[NSD*i+0];
			vel_p[1] += Ni_p[i] * el_velocity[NSD*i+1];
			vel_p[2] += Ni_p[i] * el_velocity[NSD*i+2];
		}
		
		marker_p->coor[0] = marker_p->coor[0] + step * vel_p[0];
		marker_p->coor[1] = marker_p->coor[1] + step * vel_p[1];
		marker_p->coor[2] = marker_p->coor[2] + step * vel_p[2];
	}
	
  ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#if 0
void LSF2dQ2(double _xi[],double Ni[])
{
  double xi  = _xi[0];
  double eta = _xi[1];
  
  Ni[0] = 0.5*eta*(eta-1.0)       * 0.5*xi*(xi-1.0); /*0-0*/
  Ni[1] = 0.5*eta*(eta-1.0)       * (1.0+xi)*(1.0-xi); /*0-1*/
  Ni[2] = 0.5*eta*(eta-1.0)       * 0.5*(1.0+xi)*xi; /*0-2*/
  
  Ni[3] = (1.0+eta)*(1.0-eta) * 0.5*xi*(xi-1.0); /*1-0*/
  Ni[4] = (1.0+eta)*(1.0-eta) * (1.0+xi)*(1.0-xi); /*1-1*/
  Ni[5] = (1.0+eta)*(1.0-eta) * 0.5*(1.0+xi)*xi; /*1-2*/
  
  Ni[6] = 0.5*(1.0+eta)*eta  * 0.5*xi*(xi-1.0); /*2-0*/
  Ni[7] = 0.5*(1.0+eta)*eta  * (1.0+xi)*(1.0-xi); /*2-1*/
  Ni[8] = 0.5*(1.0+eta)*eta  * 0.5*(1.0+xi)*xi; /*2-2*/
}

void LSFDeriv2dQ2(double _xi[],double **GNi)
{
  double xi  = _xi[0];
  double eta = _xi[1];
  
  GNi[0][0] = 0.5*eta*(eta-1.0)   * (xi-0.5); /*0-0*/
  GNi[0][1] = 0.5*eta*(eta-1.0)   * ( - 2.0 * xi ); /*0-1*/
  GNi[0][2] = 0.5*eta*(eta-1.0)   * 0.5*( 1.0 + 2.0 * xi ); /*0-2*/
  
  GNi[0][3] = (1.0+eta)*(1.0-eta) * (xi-0.5); /*1-0*/
  GNi[0][4] = (1.0+eta)*(1.0-eta) * ( - 2.0 * xi ); /*1-1*/
  GNi[0][5] = (1.0+eta)*(1.0-eta) * 0.5*( 1.0 + 2.0 * xi ); /*1-2*/
  
  GNi[0][6] = 0.5*(1.0+eta)*eta  * (xi-0.5); /*2-0*/\
  GNi[0][7] = 0.5*(1.0+eta)*eta  * ( - 2.0 * xi ); /*2-1*/
  GNi[0][8] = 0.5*(1.0+eta)*eta  * 0.5*( 1.0 + 2.0 * xi ); /*2-2*/
  
  GNi[1][0] = (eta - 0.5) * 0.5*xi*(xi-1.0);   /*0-0*/
  GNi[1][1] = (eta - 0.5) * (1.0+xi)*(1.0-xi); /*0-1*/
  GNi[1][2] = (eta - 0.5) * 0.5*(1.0+xi)*xi;   /*0-2*/
  
  GNi[1][3] = (-2.0*eta) * 0.5*xi*(xi-1.0);   /*1-0*/
  GNi[1][4] = (-2.0*eta) * (1.0+xi)*(1.0-xi); /*1-1*/
  GNi[1][5] = (-2.0*eta) * 0.5*(1.0+xi)*xi;   /*1-2*/
  
  GNi[1][6] = 0.5*(1.0 + 2.0*eta) * 0.5*xi*(xi-1.0); /*2-0*/
  GNi[1][7] = 0.5*(1.0 + 2.0*eta) * (1.0+xi)*(1.0-xi); /*2-1*/
  GNi[1][8] = 0.5*(1.0 + 2.0*eta) * 0.5*(1.0+xi)*xi; /*2-2*/
}

void LSF2dQ2_Interpolate(double xi[],int L,double field[],double val[])
{
	double Ni[9],sum;
	int d,i;
	
	LSF2dQ2( xi, Ni );
	for( d=0; d<L; d++ ) {
		val[d] = 0.0;
		for( i=0; i<9; i++ ) {
			val[d] = val[d] + Ni[i] * field[i*L+d];
		}
	}
}
void LSF2dQ2_CheckPartitionOfUnity(double xi[],double *val)
{
	const double tol = 1.0e-6;
	double Ni[9],sum;
	int i;
	
	LSF2dQ2(xi,Ni);
	sum = 0.0;
	for (i=0; i<9; i++) {
		sum += Ni[i];
	}
	*val = sum;
	if (fabs(sum-1.0) > tol) {
		printf("**** sum( N_i(xi,eta) ) = %1.8e > 1.0 ==> partition of unity is not satisified, point is likely outside of the element ****\n",sum);
	}
}
void LSF2dQ2_CheckGlobalCoordinate(double element_coord[],double xi[],double xp[],double err[])
{
	const double tol = 1.0e-6;
	double Ni[9],xp_interp[2];
	int i;
	
	LSF2dQ2(xi,Ni);
	xp_interp[0] = 0.0;
	xp_interp[1] = 0.0;
	
	for (i=0; i<9; i++) {
		xp_interp[0] += Ni[i] * element_coord[2*i + 0];
		xp_interp[1] += Ni[i] * element_coord[2*i + 1];
	}
	err[0] = (xp[0] - xp_interp[0]);
	err[1] = (xp[1] - xp_interp[1]);
	
	if (fabs(err[0]) > tol) {
		printf("**** |xp - N_i(xi,eta).x_i| = %1.8e > %1.8e ==> x: coordinate interpolation error occurred ****\n",fabs(err[0]),tol);
	}
	if (fabs(err[1]) > tol) {
		printf("**** |yp - N_i(xi,eta).y_i| = %1.8e > %1.8e ==> y: coordinate interpolation error occurred ****\n",fabs(err[1]),tol);
	}
}

void _compute_deltaX( double J[2][2], double f[], double h[] )
{
	double D,d1,d2,one_on_D;
	
	// Cramers rule
	D =  J[0][0] * J[1][1] - J[0][1] * J[1][0];
	d1 = f[0]    * J[1][1] - J[0][1] * f[1];
	d2 = J[0][0] * f[1]    - f[0]    * J[1][0];
	one_on_D = 1.0/D;
	
	h[0] = d1 * one_on_D;
	h[1] = d2 * one_on_D;
}

void _compute_J_2dQ2(double xi[],double vertex[],double J[2][2])
{
	int i;
	double _GNi[2][9], *GNi[2];
	
	GNi[0] = _GNi[0];
	GNi[1] = _GNi[1];
	
	J[0][0] = J[0][1] = 0.0;
	J[1][0] = J[1][1] = 0.0;
	
	LSFDeriv2dQ2(xi,GNi);
	for (i=0; i<9; i++) {
		int i2 = i*2;
		double x = vertex[i2];
		double y = vertex[i2+1];
		
		J[0][0] += x * GNi[0][i];
		J[0][1] += x * GNi[1][i];
		
		J[1][0] += y * GNi[0][i];
		J[1][1] += y * GNi[1][i];
	}
}

void _compute_F_2dQ2(double xi[],double vertex[],double pos[],double f[])
{
	int i;
	double Ni[9];
	
	/* Update F for the next iteration */
	f[0] = f[1] = 0.0;
	
	LSF2dQ2(xi,Ni);
	for (i=0; i<9; i++) {
		int i2   = i*2;
		int i2p1 = i2+1;
		
		f[0] += vertex[i2]   * Ni[i];
		f[1] += vertex[i2p1] * Ni[i];
	}
	f[0] = - f[0] + pos[0];
	f[1] = - f[1] + pos[1];
}

void InverseMappingDomain_2dQ2( 
															 double tolerance, int max_its,
															 Truth use_nonzero_guess, 
															 Truth monitor, Truth log,
															 const double coords[], const int mx, const int my, const int element[],
															 int np, MPntStd marker[] )
{
	const int dim = 2;
	const int nodesPerEl = Q2_NODES_PER_EL_2D; 
	double h[dim];
	double Jacobian[dim][dim];
	double f[dim];
	int i;
	int its;
	double residual2,tolerance2,F2;
	
	int p;
	int mx_origin, my_origin, wil_origin;
	double cxip[2],Lxip[2],Gxip[2];
	double dxi,deta,xi0,eta0;
	int I,J,wil_IJ,eid,k;
	double vertex[dim * nodesPerEl];
	int n0, n1, n2, n3;
	Truth point_found;
	
	tolerance2 = tolerance * tolerance; /* Eliminates the need to do a sqrt in the convergence test */
	
	if(log)printf("Domain: ncells = %d x %d = %d \n", mx,my,mx*my );
	
	/* map domain to [-1,1]x[-1,1] domain */
	dxi  = 2.0/((double)mx);
	deta = 2.0/((double)my);
	if(log)printf("Domain: (dxi,eta) = (%1.8e,%1.8e)\n",dxi,deta );
	
	for( p=0; p<np; p++ ) {
		MPntStd *marker_p = &marker[p];
		
		/* copy these values */
		cxip[0] = marker_p->xi[0];
		cxip[1] = marker_p->xi[1];
		
		/* Check for an initial guess initial guess */
		if( use_nonzero_guess == _FALSE ) {
			Gxip[0] = 0.0;
			Gxip[1] = 0.0;
		}
		else {
			/* convert wil => IJ */
			wil_IJ = marker_p->wil;
			J = wil_IJ/mx;
			I = wil_IJ - J*mx;
			if(log)printf("init I,J = %d %d \n", I,J );
			/* convert Lxip => Gxip */
			xi0  = -1.0 + I*dxi;
			eta0 = -1.0 + J*deta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Gxip[0] = 2.0 * (cxip[0]-xi0)/dxi - 1.0;
			//			Gxip[1] = 2.0 * (cxip[1]-eta0)/deta - 1.0;
			// x*-x*0/dx = (x+1)/2
			Gxip[0] = dxi  * (cxip[0]+1.0)/2.0 + xi0;
			Gxip[1] = deta * (cxip[1]+1.0)/2.0 + eta0;
			if(log)printf("[Lxi-init] = %1.8e %1.8e \n", cxip[0], cxip[1] );
			if(log)printf("[Gxi-init] = %1.8e %1.8e \n", Gxip[0], Gxip[1] );
			
		}
		if(monitor)printf("point[%d]: pos = ( %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e ) \n", p, marker_p->coor[0],marker_p->coor[1], marker_p->xi[0], marker_p->xi[1] );
		
		point_found = _FALSE;
		
		its = 0;
		do {
			if(log)printf("iteration: %d\n",its);
			/* convert Gxi to IJ */
			I = (Gxip[0]+1.0)/dxi;
			J = (Gxip[1]+1.0)/deta;
			
			if( I==mx ) I--;
			if( J==my ) J--;
			
			if( (I<0) || (J<0) ) {
				if(log)printf("  I(%d),J(%d) negative Gxip %1.8e,%1.8e \n",I,J,Gxip[0],Gxip[1]);
				break;
			}
			if( I>=mx ) { 
				if(log)printf("  I too large \n");
				break;
			}
			if( J>=my ) {
				if(log)printf("  J too large \n");
				break;
			}
			
			
			/* Get coords of cell IJ */
			wil_IJ = I + J * mx;
			if(log)printf("  I,J=%d/%d : wil_IJ %d : nid = ", I,J,wil_IJ);
			for (k=0; k<Q2_NODES_PER_EL_2D; k++) {
				int nid = element[wil_IJ*Q2_NODES_PER_EL_2D+k];
				
				vertex[2*k+0] = coords[2*nid+0];
				vertex[2*k+1] = coords[2*nid+1];
				if(log)printf("%d ", nid);
			}
			if(log)printf("\n");
			
			
			if(log) {
				printf("  [vertex] ");
				for (k=0; k<Q2_NODES_PER_EL_2D; k++) {
					printf("(%1.8e , %1.8e) ",vertex[2*k+0],vertex[2*k+1] );
				}
				printf("\n");
			}
			
			/* convert global (domain) xi TO local (element) xi  */
			xi0  = -1.0 + I*dxi;
			eta0 = -1.0 + J*deta;
			
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			//			Lxip[0] = 0.5*(Gxip[0]+1.0)*dxi  + xi0;
			//			Lxip[1] = 0.5*(Gxip[1]+1.0)*deta + eta0;
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0)/dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0)/deta - 1.0;
			
			if(log)printf("  Lxi,Lxeta = %1.8e, %1.8e (%d,%d) \n", Lxip[0],Lxip[1],I,J );
			
			_compute_F_2dQ2( Lxip, vertex, marker_p->coor, f );
			if( monitor ) printf("%4d InverseMapping : F = ( %+1.8e, %+1.8e ) : xi = ( %+1.8e, %+1.8e ) \n", its, f[0],f[1], Lxip[0], Lxip[1] );
			
			/* Check for convergence */
			F2 = (f[0]*f[0]+f[1]*f[1]);
			if( F2 < tolerance2 ) {
				if( monitor ) printf("%4d InverseMapping : converged : Norm of F %1.8e \n", its, sqrt(F2) );
				point_found = _TRUE;
				break;
			}
			
			_compute_J_2dQ2( Lxip, vertex, Jacobian );
			
			/* compute update */
			_compute_deltaX( Jacobian, f, h );
			if(log)printf("  [delta] = %1.8e %1.8e \n", h[0],h[1] );
			
			/* update Lxip */
			Lxip[0] += 10.0e-1 *h[0];
			Lxip[1] += 10.0e-1 *h[1];
			if(log)printf("  [corrected] Lxi,Lxeta = %1.8e, %1.8e \n", Lxip[0],Lxip[1] );
			
			residual2 = ( h[0]*h[0] + h[1]*h[1] );
			if( residual2 < tolerance2 ) {
				if( monitor ) printf("%4d InverseMapping : converged : Norm of correction %1.8e \n", its, sqrt(residual2) );
				point_found = _TRUE;
				break;
			}
			
			/* convert Lxip => Gxip */
			//			xi0  = -1.0 + I*dxi;
			//			eta0 = -1.0 + J*deta;
			
			//			Gxip[0] = 2.0 * (Lxip[0]-xi0)/dxi   - 1.0;
			//			Gxip[1] = 2.0 * (Lxip[1]-eta0)/deta - 1.0;
			// x*-x*0/dx = (x+1)/2
			Gxip[0] = dxi  * (Lxip[0]+1.0)/2.0 + xi0;
			Gxip[1] = deta * (Lxip[1]+1.0)/2.0 + eta0;
			if(log)printf("  [Gxi] = %1.8e %1.8e \n", Gxip[0], Gxip[1] );
			
			if (Gxip[0]<-1.0) { 
				Gxip[0] = -1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			if (Gxip[1]<-1.0) {
				Gxip[1] = -1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			if (Gxip[0]>1.0) {
				Gxip[0] = 1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			if (Gxip[1]>1.0) {
				Gxip[1] = 1.0;
				if(log)printf("  correction outside box: correcting \n");
			}
			
			its++;
		} while(its<max_its);
		
		if( monitor && point_found==_FALSE ){
			if( its>=max_its ) {
				printf("%4d %s : Reached maximum iterations (%d) without converging. \n", its, __FUNCTION__, max_its );
			}
			else {
				printf("%4d %s : Newton broke down, diverged or stagnated after (%d) iterations without converging. \n", its, __FUNCTION__, its );
			}
		}
		
		/* if at the end of the solve, it still looks like the point is outside the mapped domain, mark point as not being found */
		if( fabs(Gxip[0]) > 1.0 ) { point_found =_FALSE; }
		if( fabs(Gxip[1]) > 1.0 ) { point_found =_FALSE; }
		
		/* update local variables */
		if( point_found==_FALSE ) {
			Lxip[0] = NAN;
			Lxip[1] = NAN;
			wil_IJ  = -1;
		}
		else {
			/* convert Gxi to IJ */
			I = (Gxip[0]+1.0)/dxi;
			J = (Gxip[1]+1.0)/deta;
			if( I==mx ) I--;
			if( J==my ) J--;
			
			if( I>=mx ) {
				if(log)printf("  I too large \n");
				break;
			}
			if( J>=my ) {
				if(log)printf("  J too large \n");
				break;
			}
			
			/* convert global (domain) xi TO local (element) xi  */
			/*  Gxi-(-1)/2 = (xp - x0)/dx */
			xi0  = -1.0 + I*dxi;
			eta0 = -1.0 + J*deta;
			
			// x*-x*0/dx = (x+1)/2
			Lxip[0] = 2.0*(Gxip[0]-xi0)/dxi   - 1.0;
			Lxip[1] = 2.0*(Gxip[1]-eta0)/deta - 1.0;
			
			wil_IJ = I + J * mx;
		}
		
		/* set into vector */
		marker_p->xi[0] = Lxip[0];
		marker_p->xi[1] = Lxip[1];
		marker_p->wil   = wil_IJ;
	}
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_ComputeCourantStep"
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep(DM da,Vec velocity,PetscReal *step)
{
	Vec             Lvelocity, gcoords;
	PetscScalar     *LA_velocity, *LA_coords;
	PetscScalar     el_coords[Q2_NODES_PER_EL_2D*NSD];
	PetscScalar     el_velocity[Q2_NODES_PER_EL_2D*NSD];
	PetscInt        e,i;
	PetscScalar     Ni_p[Q2_NODES_PER_EL_2D],xi_p[NSD];
	PetscInt        nel,nen_u;
	const PetscInt  *elnidx_u;
	DM              cda;
	double          dt_min_local, dt_min;
	MPI_Comm        comm;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;

  /* setup for coords */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);

	/* setup velocity */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);

	dt_min_local = 1.0e32;
	for (e=0; e<nel; e++) {
		double line_x1[2],line_x2[2];
		double dt_x1,dt_x2;
		
		/* get coords for the element */
		ierr = DMDAGetElementCoordinatesQ2_2D(el_coords,(PetscInt*)&elnidx_u[nen_u*e],LA_coords);CHKERRQ(ierr);
		/* get velocity */
		ierr = DMDAGetVectorElementFieldQ2_2D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		line_x1[0] = fabs(el_coords[NSD*5+0] - el_coords[NSD*3+0]);
		line_x1[1] = fabs(el_coords[NSD*5+1] - el_coords[NSD*3+1]);

		line_x2[0] = fabs(el_coords[NSD*7+0] - el_coords[NSD*1+0]);
		line_x2[1] = fabs(el_coords[NSD*7+1] - el_coords[NSD*1+1]);
		
		dt_x1 = ( line_x1[0]*line_x1[0] + line_x1[1]*line_x1[1] ) / fabs( el_velocity[NSD*4+0] );
		dt_x2 = ( line_x2[0]*line_x2[0] + line_x2[1]*line_x2[1] ) / fabs( el_velocity[NSD*4+1] );
		
		//printf("%lf %lf \n", dt_x1, dt_x2 );
		if (dt_x1 < dt_min_local) { dt_min_local = dt_x1; }
		if (dt_x2 < dt_min_local) { dt_min_local = dt_x2; }
	}
	dt_min_local = sqrt( dt_min_local );

  ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);

	ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&dt_min_local,&dt_min,1,MPI_DOUBLE,MPI_MIN,comm);CHKERRQ(ierr);
	
	*step = dt_min;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_ComputeCourantStep_2"
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep_2(DM da,Vec velocity,PetscReal *step)
{
	Vec             Lvelocity, gcoords;
	PetscScalar     *LA_velocity, *LA_coords;
	PetscScalar     el_coords[Q2_NODES_PER_EL_2D*NSD];
	PetscScalar     el_velocity[Q2_NODES_PER_EL_2D*NSD];
	PetscInt        e,i;
	PetscScalar     Ni_p[Q2_NODES_PER_EL_2D],xi_p[NSD];
	PetscInt        nel,nen_u;
	const PetscInt  *elnidx_u;
	DM              cda;
	double          dt_min_local, dt_min;
	MPI_Comm        comm;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
  /* setup for coords */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	/* setup velocity */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	dt_min_local = 1.0e32;
	for (e=0; e<nel; e++) {
		double line_x1[NSD],line_x2[NSD];
		double dt_x1,dt_x2;
		double vc[NSD];
		
		/* get coords for the element */
		ierr = DMDAGetElementCoordinatesQ2_2D(el_coords,(PetscInt*)&elnidx_u[nen_u*e],LA_coords);CHKERRQ(ierr);
		/* get velocity */
		ierr = DMDAGetVectorElementFieldQ2_2D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		vc[0] = el_velocity[NSD*4+0];
		vc[1] = el_velocity[NSD*4+1];
		
		line_x1[0] = fabs(el_coords[NSD*5+0] - el_coords[NSD*3+0]);
		line_x1[1] = fabs(el_coords[NSD*5+1] - el_coords[NSD*3+1]);
		
		line_x2[0] = fabs(el_coords[NSD*7+0] - el_coords[NSD*1+0]);
		line_x2[1] = fabs(el_coords[NSD*7+1] - el_coords[NSD*1+1]);
		
		dt_x1 = line_x1[0] / sqrt( vc[0]*vc[0] + vc[1]*vc[1] );
		dt_x2 = line_x2[1] / sqrt( vc[0]*vc[0] + vc[1]*vc[1] );
		
		if (dt_x1 < dt_min_local) { dt_min_local = dt_x1; }
		if (dt_x2 < dt_min_local) { dt_min_local = dt_x2; }
	}
  ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	
	ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&dt_min_local,&dt_min,1,MPI_DOUBLE,MPI_MIN,comm);CHKERRQ(ierr);
	
	*step = dt_min;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdatePosition_ComputeCourantStep_3"
PetscErrorCode SwarmUpdatePosition_ComputeCourantStep_3(DM da,Vec velocity,PetscReal *step)
{
	Vec             Lvelocity, gcoords;
	PetscScalar     *LA_velocity, *LA_coords;
	PetscScalar     el_coords[Q2_NODES_PER_EL_2D*NSD];
	PetscScalar     el_velocity[Q2_NODES_PER_EL_2D*NSD];
	PetscInt        e,i;
	PetscScalar     Ni_p[Q2_NODES_PER_EL_2D],xi_p[NSD];
	PetscInt        nel,nen_u;
	const PetscInt  *elnidx_u;
	DM              cda;
	double          dt_min_local, dt_min;
	MPI_Comm        comm;
	PetscErrorCode  ierr;
	
	PetscFunctionBegin;
	
  /* setup for coords */
  ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
  ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
  ierr = VecGetArray(gcoords,&LA_coords);CHKERRQ(ierr);
	
	/* setup velocity */
	ierr = DMGetLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,velocity,INSERT_VALUES,Lvelocity);CHKERRQ(ierr);
	ierr = VecGetArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	dt_min_local = 1.0e32;
	for (e=0; e<nel; e++) {
		double line_x1[NSD],line_x2[NSD];
		double dt_x1,dt_x2,dt_step;
		double vc[NSD];
		
		/* get coords for the element */
		ierr = DMDAGetElementCoordinatesQ2_2D(el_coords,(PetscInt*)&elnidx_u[nen_u*e],LA_coords);CHKERRQ(ierr);
		/* get velocity */
		ierr = DMDAGetVectorElementFieldQ2_2D(el_velocity,(PetscInt*)&elnidx_u[nen_u*e],LA_velocity);CHKERRQ(ierr);
		
		vc[0] = el_velocity[NSD*4+0];
		vc[1] = el_velocity[NSD*4+1];
		
		line_x1[0] = fabs(el_coords[NSD*5+0] - el_coords[NSD*3+0]);
		line_x1[1] = fabs(el_coords[NSD*5+1] - el_coords[NSD*3+1]);
		
		line_x2[0] = fabs(el_coords[NSD*7+0] - el_coords[NSD*1+0]);
		line_x2[1] = fabs(el_coords[NSD*7+1] - el_coords[NSD*1+1]);
		
		dt_x1 = line_x1[0] / fabs( vc[0] );
		dt_x2 = line_x2[1] / fabs( vc[1] );
		
		dt_step = ( dt_x1*dt_x1 + dt_x2*dt_x2 );
		
		if (dt_step < dt_min_local) { dt_min_local = dt_step; }
	}
	dt_min_local = sqrt( dt_min_local );
	
  ierr = VecRestoreArray(gcoords,&LA_coords);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lvelocity,&LA_velocity);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&Lvelocity);CHKERRQ(ierr);
	
	ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
	ierr = MPI_Allreduce(&dt_min_local,&dt_min,1,MPI_DOUBLE,MPI_MIN,comm);CHKERRQ(ierr);
	
	*step = dt_min;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_SwarmUpdatePosition_Communication_MPntStd"
PetscErrorCode _SwarmUpdatePosition_Communication_MPntStd(DM da,DataBucket db,DataEx de)
{
	DataField            PField;
	PetscInt p,npoints,npoints_global_init,npoints_global_fin;
	MPntStd *recv_data;
	int n,neighborcount, *neighborranks2;
	int recv_length,npoints_accepted;
	PetscMPIInt rank,size;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* communucate */
	MPI_Comm_size(((PetscObject)da)->comm,&size);
	if (size==1) {
		PetscFunctionReturn(0);
	}
	
	MPI_Comm_rank(((PetscObject)da)->comm,&rank);
	
	neighborcount  = de->n_neighbour_procs;
	neighborranks2 = de->neighbour_procs;
	
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField);
	DataFieldGetAccess(PField);
	DataFieldVerifyAccess( PField,sizeof(MPntStd));
	DataBucketGetSizes(db,&npoints,0,0);
	
	MPI_Allreduce(&npoints,&npoints_global_init,1,MPI_INT,MPI_SUM,de->comm);
	
	/* figure out how many points left processor */
	ierr = DataExInitializeSendCount(de);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		PetscBool onproc;
		MPntStd   *marker;
		
		DataFieldAccessPoint(PField,p,(void**)&marker);
		onproc = PETSC_TRUE;
		if (marker->wil==-1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc==PETSC_FALSE) {
			for (n=0; n<neighborcount; n++) {
				//	printf("  DataEx: rank %d sending %d to %d \n",rank,1,neighborranks2[n] );
				ierr = DataExAddToSendCount( de, neighborranks2[n], 1 );CHKERRQ(ierr);
			}
		}
	}
	ierr = DataExFinalizeSendCount(de);CHKERRQ(ierr);
	
	DataFieldRestoreAccess(PField);
	
	/* pack points which left processor */
	DataFieldGetAccess(PField);
	ierr = DataExPackInitialize(de,sizeof(MPntStd));CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		PetscBool onproc;
		MPntStd   *marker;
		
		DataFieldAccessPoint(PField,p,(void**)&marker);
		onproc = PETSC_TRUE;
		if (marker->wil==-1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc==PETSC_FALSE) {
			for (n=0; n<neighborcount; n++) {
				ierr = DataExPackData( de, neighborranks2[n], 1,(void*)marker );CHKERRQ(ierr);
			}
		}
	}		
	ierr = DataExPackFinalize(de);CHKERRQ(ierr);
	DataFieldRestoreAccess(PField);
	
	/* remove points which left processor */
	DataBucketGetSizes(db,&npoints,0,0);
	DataFieldGetAccess(PField);
	for (p=0; p<npoints; p++) {
		PetscBool onproc;
		MPntStd   *marker;
		
		DataFieldAccessPoint(PField,p,(void**)&marker);
		onproc = PETSC_TRUE;
		if (marker->wil==-1) {
			onproc = PETSC_FALSE;
		}
		
		if (onproc==PETSC_FALSE) { 
			/* kill point */
			DataBucketRemovePointAtIndex(db,p);
			DataBucketGetSizes(db,&npoints,0,0); /* you need to update npoints as the list size decreases! */
			p--; /* check replacement point */
		}
	}		
	DataFieldRestoreAccess(PField);
	
	// START communicate //
	ierr = DataExBegin(de);CHKERRQ(ierr);
	ierr = DataExEnd(de);CHKERRQ(ierr);
	// END communicate //
	
	// receive, if i own them, add new points to list //
	ierr = DataExGetRecvData( de, &recv_length, (void**)&recv_data );CHKERRQ(ierr);
	{
		int totalsent;
		MPI_Allreduce(&recv_length,&totalsent,1,MPI_INT,MPI_SUM,de->comm);
		PetscPrintf(PETSC_COMM_WORLD,"  DataEx: total points sent = %d \n", totalsent);
	}
	
	
	
	/* update the local coordinates and cell owner for all recieved points */
	{
		DM cda;
		Vec gcoords;
		PetscScalar *LA_gcoords;
		double tolerance;
		int max_its;
		Truth use_nonzero_guess, monitor, log;
		PetscInt lmx,lmy;
		PetscInt nel,nen_u;
		MPntStd *local_list;
		const PetscInt *elnidx_u;
		
		/* setup for coords */
		ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
		ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
		ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
		
		ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
		
		ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,PETSC_NULL);CHKERRQ(ierr);
		
		/* point location parameters */
		tolerance         = 1.0e-10;
		max_its           = 10;
		use_nonzero_guess = _TRUE;
		monitor           = _FALSE;
		log               = _FALSE;
		
		local_list = (MPntStd*)recv_data;
		InverseMappingDomain_2dQ2( 		 tolerance, max_its,
															use_nonzero_guess, 
															monitor, log,
															(const double*)LA_gcoords, (const int)lmx,(const int)lmy, (const int*)elnidx_u,
															recv_length, local_list );
		
		ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	}
	
	
	/* accepte all points living locally */
	npoints_accepted = 0;
	for (p=0; p<recv_length; p++) {
		PetscBool onproc;
		MPntStd   *marker;
		
		marker = &recv_data[p];
		
		onproc = PETSC_TRUE;
		if (marker->wil==-1) {
			onproc = PETSC_FALSE;
		}

		if (onproc == PETSC_TRUE) {
			int end;
			
			DataBucketAddPoint(db);
			DataBucketGetSizes(db,&end,0,0);
			end = end - 1;
			DataFieldInsertPoint(PField, end, (void*)marker );
			npoints_accepted++;
		}
	}	
	//printf("  DataEx: rank %d accepted %d new points \n",rank,npoints_accepted );
	
	
	DataBucketGetSizes(db,&npoints,0,0);
	MPI_Allreduce(&npoints,&npoints_global_fin,1,MPI_INT,MPI_SUM,de->comm);
	PetscPrintf(PETSC_COMM_WORLD,"  SwarmUpdatePosition(Communication): num. points global ( init. = %d : final = %d )\n", npoints_global_init,npoints_global_fin);
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdateProperties_MPntStd"
PetscErrorCode SwarmUpdateProperties_MPntStd(DataBucket db,pTatinCtx ctx,Vec X)
{
	BTruth found;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	DataBucketQueryDataFieldByName(db,MPntStd_classname,&found);
	if(found==BFALSE) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cannot find DataField with name %s \n", MPntStd_classname );
	}
	
	
	
	PetscFunctionReturn(0);
}
#endif


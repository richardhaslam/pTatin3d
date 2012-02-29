

#include "petsc.h"

#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "dmda_duplicate.h"
#include "dmda_element_q2p1.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "output_paraview.h"
#include "quadrature.h"
#include "element_type_Q2.h"
#include "material_point_utils.h"
#include "element_utils_q2.h"


#undef __FUNCT__
#define __FUNCT__ "MaterialPointGeneric_VTKWriteBinaryAppendedHeaderAllFields"
PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,DataBucket db,int *byte_offset,const int nfields,const MaterialPointField list[])
{
	PetscErrorCode ierr;
	int n,npoints;
	
	PetscFunctionBegin;

	DataBucketGetSizes(db,&npoints,PETSC_NULL,PETSC_NULL);
	for (n=0; n<nfields; n++) {
			switch (list[n]) {

				/* auto generated shit for the marker data goes here */
				case MPField_Std:
				{
					DataField PField_std;
					MPntStd   *marker_std;
					
					DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField_std);
					DataFieldGetAccess(PField_std);
					marker_std = PField_std->data;
					
					MPntStdVTKWriteBinaryAppendedHeaderAllFields(vtk_fp,byte_offset,(const int)npoints,(const MPntStd*)marker_std);
					DataFieldRestoreAccess(PField_std);
				}
					break;
				
				case MPField_Stokes:
				{
					DataField   PField_stokes;
					MPntPStokes *marker_stokes;
					
					DataBucketGetDataFieldByName(db, MPntPStokes_classname ,&PField_stokes);
					DataFieldGetAccess(PField_stokes);
					marker_stokes = PField_stokes->data;

					MPntPStokesVTKWriteBinaryAppendedHeaderAllFields(vtk_fp,byte_offset,(const int)npoints,(const MPntPStokes*)marker_stokes);
					DataFieldRestoreAccess(PField_stokes);
				}
					break;
					
				default:
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown material point field");
					break;
			}
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointGeneric_VTKWriteBinaryAppendedDataAllFields"
PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,DataBucket db,const int nfields,const MaterialPointField list[])
{
	PetscErrorCode ierr;
	int n,npoints;
	
	PetscFunctionBegin;
	
	DataBucketGetSizes(db,&npoints,PETSC_NULL,PETSC_NULL);
	for (n=0; n<nfields; n++) {
		switch (list[n]) {
				
			/* auto generated shit for the marker data goes here */
			case MPField_Std:
			{
				DataField PField_std;
				MPntStd   *marker_std;
				
				DataBucketGetDataFieldByName(db, MPntStd_classname ,&PField_std);
				DataFieldGetAccess(PField_std);
				marker_std = PField_std->data;
				
				MPntStdVTKWriteBinaryAppendedDataAllFields(vtk_fp,(const int)npoints,(const MPntStd*)marker_std);
				DataFieldRestoreAccess(PField_std);
			}
				break;
				
			case MPField_Stokes:
			{
				DataField   PField_stokes;
				MPntPStokes *marker_stokes;
				
				DataBucketGetDataFieldByName(db, MPntPStokes_classname ,&PField_stokes);
				DataFieldGetAccess(PField_stokes);
				marker_stokes = PField_stokes->data;
				
				MPntPStokesVTKWriteBinaryAppendedDataAllFields(vtk_fp,(const int)npoints,(const MPntPStokes*)marker_stokes);
				DataFieldRestoreAccess(PField_stokes);
			}
				break;
				
			default:
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown material point field");
				break;
		}
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointGeneric_PVTUWriteAllPPointDataFields"
PetscErrorCode MaterialPointGeneric_PVTUWriteAllPPointDataFields(FILE *vtk_fp,const int nfields,const MaterialPointField list[])
{
	PetscErrorCode ierr;
	int n;
	
	PetscFunctionBegin;
	
	for (n=0; n<nfields; n++) {
		switch (list[n]) {
				
				/* auto generated shit for the marker data goes here */
			case MPField_Std:
				MPntStdPVTUWriteAllPPointDataFields(vtk_fp);
				break;
				
			case MPField_Stokes:
				MPntPStokesPVTUWriteAllPPointDataFields(vtk_fp);
				break;
				
			default:
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown material point field");
				break;
		}
	}
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SwarmViewGeneric_VTUXML_binary_appended"
PetscErrorCode SwarmViewGeneric_VTUXML_binary_appended(DataBucket db,const int nmpfields,const MaterialPointField mplist[],const char name[])
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
		ierr = MaterialPointGeneric_VTKWriteBinaryAppendedHeaderAllFields(vtk_fp,db,&byte_offset,nmpfields,mplist);CHKERRQ(ierr);
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
		
		ierr = MaterialPointGeneric_VTKWriteBinaryAppendedDataAllFields(vtk_fp,db,nmpfields,mplist);CHKERRQ(ierr);
	
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
#define __FUNCT__ "SwarmViewGeneric_PVTUXML"
PetscErrorCode SwarmViewGeneric_PVTUXML(const int nfields,const MaterialPointField list[],const char prefix[],const char name[])
{
	PetscMPIInt nproc,rank;
	FILE *vtk_fp;
	PetscInt i;
	char *sourcename;
	PetscErrorCode ierr;
	
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
	{
		/* auto generated shit for the marker data goes here */
		ierr = MaterialPointGeneric_PVTUWriteAllPPointDataFields(vtk_fp,nfields,list);CHKERRQ(ierr);
	}
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
#define __FUNCT__ "SwarmViewGeneric_ParaView"
PetscErrorCode SwarmViewGeneric_ParaView(DataBucket db,const int nfields,const MaterialPointField list[],const char path[],const char prefix[])
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
	ierr = SwarmViewGeneric_VTUXML_binary_appended(db,nfields,list,filename);CHKERRQ(ierr);
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
		ierr = SwarmViewGeneric_PVTUXML(nfields,list,prefix,filename);CHKERRQ(ierr);
	}
	free(filename);
	free(vtkfilename);
	
	PetscFunctionReturn(0);
}

void pTatinConstructNI_Q1_3D(const double _xi[],double Ni[])
{
  PetscScalar xi   = _xi[0];
  PetscScalar eta  = _xi[1];
  PetscScalar zeta = _xi[2];
  
  Ni[0] = 0.125*(1.0-xi)*(1.0-eta)*(1.0-zeta); /*0-0*/
  Ni[1] = 0.125*(1.0+xi)*(1.0-eta)*(1.0-zeta); /*0-1*/
  
  Ni[2] = 0.125*(1.0-xi)*(1.0+eta)*(1.0-zeta); /*1-0*/
  Ni[3] = 0.125*(1.0+xi)*(1.0+eta)*(1.0-zeta); /*1-1*/

  Ni[4] = 0.125*(1.0-xi)*(1.0-eta)*(1.0+zeta); /*0-0+1*/
  Ni[5] = 0.125*(1.0+xi)*(1.0-eta)*(1.0+zeta); /*0-1+1*/
  
  Ni[6] = 0.125*(1.0-xi)*(1.0+eta)*(1.0+zeta); /*1-0+1*/
  Ni[7] = 0.125*(1.0+xi)*(1.0+eta)*(1.0+zeta); /*1-1+1*/
}

/*

 6[2]--7--8[3]
  |         |
  3     4   5
  |         |
 0[0]--1--2[1]

 
 15---16---17
 |          |
 12   13   14
 |          |
 9----10---11

 
 24[6]---25---26[7]
  |              |
  21     22     23
  |              |
 18[4]---19---20[5]
 
*/
void pTatinConstructNI_Q1_on_Q2_3D(const double _xi[],double Ni[])
{
	PetscScalar NiQ1[8];
  PetscScalar xi   = _xi[0];
  PetscScalar eta  = _xi[1];
  PetscScalar zeta = _xi[2];
  
  NiQ1[   0] = 0.125*(1.0-xi)*(1.0-eta)*(1.0-zeta); /*0-0*/
  NiQ1[   1] = 0.125*(1.0+xi)*(1.0-eta)*(1.0-zeta); /*0-1*/
  
  NiQ1[   2] = 0.125*(1.0-xi)*(1.0+eta)*(1.0-zeta); /*1-0*/
  NiQ1[   3] = 0.125*(1.0+xi)*(1.0+eta)*(1.0-zeta); /*1-1*/
	
  NiQ1[   4] = 0.125*(1.0-xi)*(1.0-eta)*(1.0+zeta); /*0-0+1*/
  NiQ1[   5] = 0.125*(1.0+xi)*(1.0-eta)*(1.0+zeta); /*0-1+1*/
  
  NiQ1[   6] = 0.125*(1.0-xi)*(1.0+eta)*(1.0+zeta); /*1-0+1*/
  NiQ1[   7] = 0.125*(1.0+xi)*(1.0+eta)*(1.0+zeta); /*1-1+1*/

	/* vertex guys */
	Ni[   0] = NiQ1[0];
	Ni[   2] = NiQ1[1];
	Ni[   6] = NiQ1[2];
	Ni[   8] = NiQ1[3];
	
	Ni[  18] = NiQ1[4];
	Ni[  20] = NiQ1[5];
	Ni[  24] = NiQ1[6];
	Ni[  26] = NiQ1[7];
	
	/* edge guys */
	Ni[ 1] = 0.5 * ( NiQ1[0] + NiQ1[1] );
	Ni[ 3] = 0.5 * ( NiQ1[0] + NiQ1[2] );
	Ni[ 5] = 0.5 * ( NiQ1[1] + NiQ1[3] );
	Ni[ 7] = 0.5 * ( NiQ1[2] + NiQ1[3] );

	Ni[ 9] = 0.5 * ( NiQ1[0] + NiQ1[4] );
	Ni[11] = 0.5 * ( NiQ1[1] + NiQ1[5] );
	
	Ni[15] = 0.5 * ( NiQ1[2] + NiQ1[6] );
	Ni[17] = 0.5 * ( NiQ1[3] + NiQ1[7] );
	
	Ni[19] = 0.5 * ( NiQ1[4] + NiQ1[5] );
	Ni[21] = 0.5 * ( NiQ1[4] + NiQ1[6] );
	Ni[23] = 0.5 * ( NiQ1[5] + NiQ1[7] );
	Ni[25] = 0.5 * ( NiQ1[6] + NiQ1[7] );
	
	/* face */
	Ni[ 4] = 0.5 * ( NiQ1[0] + NiQ1[1] + NiQ1[2] + NiQ1[3] );
	Ni[10] = 0.5 * ( NiQ1[0] + NiQ1[1] + NiQ1[4] + NiQ1[5] );
	Ni[12] = 0.5 * ( NiQ1[0] + NiQ1[2] + NiQ1[4] + NiQ1[6] );
	Ni[14] = 0.5 * ( NiQ1[1] + NiQ1[3] + NiQ1[5] + NiQ1[7] );
	Ni[16] = 0.5 * ( NiQ1[2] + NiQ1[3] + NiQ1[6] + NiQ1[7] );
	Ni[22] = 0.5 * ( NiQ1[4] + NiQ1[5] + NiQ1[6] + NiQ1[7] );
	
	/* center */
	Ni[13] = 0.125 * ( NiQ1[0] + NiQ1[1] + NiQ1[2] + NiQ1[3] + NiQ1[4] + NiQ1[5] + NiQ1[6] + NiQ1[7] );
	
}


#undef __FUNCT__
#define __FUNCT__ "_SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes"
PetscErrorCode _SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes(
																																					DM clone,Vec properties_A1,Vec properties_A2,Vec properties_B,
																																					const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],Quadrature Q) 
{
	PetscScalar Ni_p[Q2_NODES_PER_EL_3D];
	PetscScalar NiQ1_p[8];
	PetscScalar Ae1[Q2_NODES_PER_EL_3D], Ae2[Q2_NODES_PER_EL_3D], Be[Q2_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A1, Lproperties_A2, Lproperties_B;
	PetscScalar *LA_properties_A1, *LA_properties_A2, *LA_properties_B;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e;
	const PetscInt *elnidx;
	
	PetscInt ngp;
	PetscScalar *xi_mp;
	PetscScalar NIu[MAX_QUAD_PNTS][U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscScalar range_eta[2],range_rho[2];
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_A2);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A2);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_B);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_B, &LA_properties_B);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(Q,&all_gausspoints);CHKERRQ(ierr);

	PetscGetTime(&t0);
	for (p=0; p<npoints; p++) {
		double *xi_p  = &mp_std[p].xi[0];
		double eta_p  = mp_stokes[p].eta;
		double rho_p  = mp_stokes[p].rho;
		
		ierr = PetscMemzero(Ae1,sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		ierr = PetscMemzero(Ae2,sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		ierr = PetscMemzero(Be, sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		
		pTatinConstructNI_Q1_3D(xi_p,NiQ1_p);
		
		Ni_p[0] = NiQ1_p[0];
		Ni_p[2] = NiQ1_p[1];
		Ni_p[6] = NiQ1_p[2];
		Ni_p[8] = NiQ1_p[3];

		Ni_p[0+18] = NiQ1_p[4];
		Ni_p[2+18] = NiQ1_p[5];
		Ni_p[6+18] = NiQ1_p[6];
		Ni_p[8+18] = NiQ1_p[7];
		
		
		Ni_p[1] = Ni_p[7] = 1.0;
		Ni_p[3] = Ni_p[4] = Ni_p[5] = 1.0;

		Ni_p[ 9] = Ni_p[10] = Ni_p[11] = 1.0;
		Ni_p[12] = Ni_p[13] = Ni_p[14] = 1.0;
		Ni_p[15] = Ni_p[16] = Ni_p[17] = 1.0;
		
		Ni_p[1+18] = Ni_p[7+18] = 1.0;
		Ni_p[3+18] = Ni_p[4+18] = Ni_p[5+18] = 1.0;
		
		
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			Ae1[i] = Ni_p[i] * eta_p;
			Ae2[i] = Ni_p[i] * rho_p;
			Be[i]  = Ni_p[i];
		}
		
		
		
		
		/* sum into local vectors */
		e = mp_std[p].wil;
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_A1, 1, el_lidx,Ae1);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_A2, 1, el_lidx,Ae2);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_B,  1, el_lidx,Be);CHKERRQ(ierr);
		
	}
	PetscGetTime(&t1);
	//PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (summation): %1.4lf ]\n",t1-t0);
	
  ierr = VecRestoreArray(Lproperties_B,&LA_properties_B);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A2,ADD_VALUES,properties_A2);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A2,ADD_VALUES,properties_A2);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	
	/* scale */
	ierr = VecPointwiseDivide( properties_A1, properties_A1, properties_B );CHKERRQ(ierr);
	ierr = VecPointwiseDivide( properties_A2, properties_A2, properties_B );CHKERRQ(ierr);
	/* ========================================= */
	
	/* scatter result back to local array and do the interpolation onto the quadrature points */
	ngp       = Q->npoints;
	xi_mp     = Q->q_xi_coor;
	for (p=0; p<ngp; p++) {
		PetscScalar *xip = &xi_mp[3*p];
		
		ierr = PetscMemzero(NIu[p], sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		
		pTatinConstructNI_Q1_3D(xip,NiQ1_p);
		NIu[p][0] = NiQ1_p[0];
		NIu[p][2] = NiQ1_p[1];
		NIu[p][6] = NiQ1_p[2];
		NIu[p][8] = NiQ1_p[3];

		NIu[p][0+18] = NiQ1_p[4];
		NIu[p][2+18] = NiQ1_p[5];
		NIu[p][6+18] = NiQ1_p[6];
		NIu[p][8+18] = NiQ1_p[7];
	}
	
	
	PetscGetTime(&t0);
	ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = VecZeroEntries(Lproperties_A2);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A2,INSERT_VALUES,Lproperties_A2);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A2,INSERT_VALUES,Lproperties_A2);CHKERRQ(ierr);
	PetscGetTime(&t1);
//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (scatter): %1.4lf ]\n",t1-t0);
	
	PetscGetTime(&t0);
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	for (e=0;e<nel;e++) {
		ierr = VolumeQuadratureGetCellData_Stokes(Q,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		ierr = DMDAGetScalarElementField(Ae1,nen,(PetscInt*)&elnidx[nen*e],LA_properties_A1);CHKERRQ(ierr);
		ierr = DMDAGetScalarElementField(Ae2,nen,(PetscInt*)&elnidx[nen*e],LA_properties_A2);CHKERRQ(ierr);
		
		/* zero out the junk */
		Ae1[1] = Ae1[7] = 0.0;
		Ae1[3] = Ae1[4] = Ae1[5] = 0.0;
		
		Ae2[1] = Ae2[7] = 0.0;
		Ae2[3] = Ae2[4] = Ae2[5] = 0.0;
		
		Ni_p[ 9] = Ni_p[10] = Ni_p[11] = 1.0;
		Ni_p[12] = Ni_p[13] = Ni_p[14] = 1.0;
		Ni_p[15] = Ni_p[16] = Ni_p[17] = 1.0;
		
		Ni_p[1+18] = Ni_p[7+18] = 0.0;
		Ni_p[3+18] = Ni_p[4+18] = Ni_p[5+18] = 0.0;
		
		
		/* The Q2 interpolant tends to overshoot / undershoot when you have viscosity jumps */
		range_eta[0] = 1.0e32;  /* min */
		range_eta[1] = -1.0e32; /* max */
		range_rho[0] = 1.0e32;
		range_rho[1] = -1.0e32;
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			if (Ae1[i]<range_eta[0]) { range_eta[0] = Ae1[i]; }
			if (Ae1[i]>range_eta[1]) { range_eta[1] = Ae1[i]; }
			if (Ae2[i]<range_rho[0]) { range_rho[0] = Ae2[i]; }
			if (Ae2[i]>range_rho[1]) { range_rho[1] = Ae2[i]; }
		}
		
		for (p=0; p<ngp; p++) {
			cell_gausspoints[p].eta = 0.0;
			cell_gausspoints[p].Fu[0] = 0.0;
			cell_gausspoints[p].Fu[1] = 0.0;
			cell_gausspoints[p].Fu[2] = 0.0;
			for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
				cell_gausspoints[p].eta    += NIu[p][i] * Ae1[i];
				cell_gausspoints[p].Fu[1]  += NIu[p][i] * Ae2[i];
			}
			if (cell_gausspoints[p].eta < range_eta[0]) { cell_gausspoints[p].eta = range_eta[0]; }
			if (cell_gausspoints[p].eta > range_eta[1]) { cell_gausspoints[p].eta = range_eta[1]; }
			if (cell_gausspoints[p].Fu[1] < range_rho[0]) { cell_gausspoints[p].Fu[1] = range_rho[0]; }
			if (cell_gausspoints[p].Fu[1] > range_rho[1]) { cell_gausspoints[p].Fu[1] = range_rho[1]; }
			
		}
	}
	
  ierr = VecRestoreArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	PetscGetTime(&t1);
//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A2);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes"
PetscErrorCode SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],DM da,Quadrature Q)
{
	PetscInt  dof;
	DM        clone;
	Vec       properties_A1, properties_A2, properties_B;
	PetscBool view;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* setup */
	dof = 1;
	ierr = DMDADuplicateLayout(da,dof,2,DMDA_STENCIL_BOX,&clone);CHKERRQ(ierr); /* Q2 - but we'll fake it as a Q1 with cells the same size as the Q2 guys */
	
	ierr = DMGetGlobalVector(clone,&properties_A1);CHKERRQ(ierr);  ierr = PetscObjectSetName( (PetscObject)properties_A1, "LocalL2ProjQ1_nu");CHKERRQ(ierr);
	ierr = DMGetGlobalVector(clone,&properties_A2);CHKERRQ(ierr);  ierr = PetscObjectSetName( (PetscObject)properties_A2, "LocalL2ProjQ1_rho");CHKERRQ(ierr);
	ierr = DMGetGlobalVector(clone,&properties_B);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(properties_A1);CHKERRQ(ierr);
	ierr = VecZeroEntries(properties_A2);CHKERRQ(ierr);
	ierr = VecZeroEntries(properties_B);CHKERRQ(ierr);
	
	/* compute */
	ierr = _SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes(
								clone, properties_A1,properties_A2,properties_B, 																																	 
								npoints, mp_std,mp_stokes, Q );CHKERRQ(ierr);
	
	/* view */
	view = PETSC_FALSE;
	PetscOptionsGetBool(PETSC_NULL,"-view_projected_marker_fields",&view,PETSC_NULL);
	if (view) {
		PetscViewer viewer;
		
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "SwarmUpdateProperties_LocalL2Proj_Stokes.vtk", &viewer);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = DMView(clone, viewer);CHKERRQ(ierr);
		ierr = VecView(properties_A1, viewer);CHKERRQ(ierr);
		ierr = VecView(properties_A2, viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	}
	
	
	/* destroy */
	ierr = DMRestoreGlobalVector(clone,&properties_B);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(clone,&properties_A2);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(clone,&properties_A1);CHKERRQ(ierr);
	
	ierr = DMDestroy(&clone);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "_SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints"
PetscErrorCode _SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints(
																																					DM clone,Vec properties_A1,Vec properties_A2,
																																					Quadrature Q) 
{
	PetscScalar Ni_p[Q2_NODES_PER_EL_3D];
	PetscScalar Ae1[Q2_NODES_PER_EL_3D], Ae2[Q2_NODES_PER_EL_3D], Be[Q2_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A1, Lproperties_A2;
	PetscScalar *LA_properties_A1, *LA_properties_A2;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e;
	const PetscInt *elnidx;
	
	PetscInt ngp;
	PetscScalar *xi_mp;
	PetscScalar NIu[MAX_QUAD_PNTS][U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscScalar range_eta[2],range_rho[2];
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_A2);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A2);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(Q,&all_gausspoints);CHKERRQ(ierr);
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A2,ADD_VALUES,properties_A2);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A2,ADD_VALUES,properties_A2);CHKERRQ(ierr);
	
	/* ========================================= */
	
	/* scatter result back to local array and do the interpolation onto the quadrature points */
	ngp       = Q->npoints;
	xi_mp     = Q->q_xi_coor;
	for (p=0; p<ngp; p++) {
		PetscScalar *xip = &xi_mp[3*p];
		
		ierr = PetscMemzero(NIu[p], sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		
		P3D_ConstructNi_Q2_3D(xip,NIu[p]);
	}
	
	PetscGetTime(&t0);
	ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = VecZeroEntries(Lproperties_A2);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A2,INSERT_VALUES,Lproperties_A2);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A2,INSERT_VALUES,Lproperties_A2);CHKERRQ(ierr);
	PetscGetTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (scatter): %1.4lf ]\n",t1-t0);
	
	PetscGetTime(&t0);
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	for (e=0;e<nel;e++) {
		ierr = VolumeQuadratureGetCellData_Stokes(Q,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		ierr = DMDAGetScalarElementField(Ae1,nen,(PetscInt*)&elnidx[nen*e],LA_properties_A1);CHKERRQ(ierr);
		ierr = DMDAGetScalarElementField(Ae2,nen,(PetscInt*)&elnidx[nen*e],LA_properties_A2);CHKERRQ(ierr);
		
		/* The Q2 interpolant tends to overshoot / undershoot when you have viscosity jumps */
		range_eta[0] = 1.0e32;  /* min */
		range_eta[1] = -1.0e32; /* max */
		range_rho[0] = 1.0e32;
		range_rho[1] = -1.0e32;
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			if (Ae1[i]<range_eta[0]) { range_eta[0] = Ae1[i]; }
			if (Ae1[i]>range_eta[1]) { range_eta[1] = Ae1[i]; }
			if (Ae2[i]<range_rho[0]) { range_rho[0] = Ae2[i]; }
			if (Ae2[i]>range_rho[1]) { range_rho[1] = Ae2[i]; }
		}
		
		for (p=0; p<ngp; p++) {
			cell_gausspoints[p].eta = 0.0;
			cell_gausspoints[p].Fu[0] = 0.0;
			cell_gausspoints[p].Fu[1] = 0.0;
			cell_gausspoints[p].Fu[2] = 0.0;

			cell_gausspoints[p].Fp = 0.0;

			for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
				cell_gausspoints[p].eta    += NIu[p][i] * Ae1[i];
				cell_gausspoints[p].Fu[1]  += NIu[p][i] * Ae2[i];
			}
			if (cell_gausspoints[p].eta < range_eta[0]) { cell_gausspoints[p].eta = range_eta[0]; }
			if (cell_gausspoints[p].eta > range_eta[1]) { cell_gausspoints[p].eta = range_eta[1]; }
			if (cell_gausspoints[p].Fu[1] < range_rho[0]) { cell_gausspoints[p].Fu[1] = range_rho[0]; }
			if (cell_gausspoints[p].Fu[1] > range_rho[1]) { cell_gausspoints[p].Fu[1] = range_rho[1]; }
		}
	}
	
  ierr = VecRestoreArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	PetscGetTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_A2);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_FineGrid"
PetscErrorCode _SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_FineGrid(
																																					DM clone,Vec properties_A1,Vec properties_A2,Vec properties_B,
																																					const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[]) 
{
	PetscScalar Ni_p[Q2_NODES_PER_EL_3D];
	PetscScalar Ae1[Q2_NODES_PER_EL_3D], Ae2[Q2_NODES_PER_EL_3D], Be[Q2_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A1, Lproperties_A2, Lproperties_B;
	PetscScalar *LA_properties_A1, *LA_properties_A2, *LA_properties_B;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e;
	const PetscInt *elnidx;
	
	PetscInt ngp;
	PetscScalar *xi_mp;
	PetscScalar NIu[MAX_QUAD_PNTS][U_BASIS_FUNCTIONS];
	PetscScalar range_eta[2],range_rho[2];
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_A2);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A2);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_B);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_B, &LA_properties_B);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	PetscGetTime(&t0);
	for (p=0; p<npoints; p++) {
		double *xi_p  = &mp_std[p].xi[0];
		double eta_p  = mp_stokes[p].eta;
		double rho_p  = mp_stokes[p].rho;
		
		ierr = PetscMemzero(Ae1,sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		ierr = PetscMemzero(Ae2,sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		ierr = PetscMemzero(Be, sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		
		pTatinConstructNI_Q1_on_Q2_3D(xi_p,Ni_p);
		
		for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
			Ae1[i] = Ni_p[i] * eta_p;
			Ae2[i] = Ni_p[i] * rho_p;
			Be[i]  = Ni_p[i];
		}
		
		/* sum into local vectors */
		e = mp_std[p].wil;
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_A1, 1, el_lidx,Ae1);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_A2, 1, el_lidx,Ae2);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_B,  1, el_lidx,Be);CHKERRQ(ierr);
		
	}
	PetscGetTime(&t1);
	//PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (summation): %1.4lf ]\n",t1-t0);
	
  ierr = VecRestoreArray(Lproperties_B,&LA_properties_B);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A2,ADD_VALUES,properties_A2);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A2,ADD_VALUES,properties_A2);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	
	/* scale */
	ierr = VecPointwiseDivide( properties_A1, properties_A1, properties_B );CHKERRQ(ierr);
	ierr = VecPointwiseDivide( properties_A2, properties_A2, properties_B );CHKERRQ(ierr);
	/* ========================================= */
	
	
	PetscGetTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A2);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy"
PetscErrorCode SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,Mat R[],DM da[],Quadrature Q[])
{
	PetscInt  dof,k;
	DM        clone[100];
	Vec       properties_A1[100], properties_A2[100], properties_B;
	PetscInt  ptype;
	PetscBool view,flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	/* setup */
	dof = 1;
	for (k=0; k<nlevels; k++) {
		ierr = DMDADuplicateLayout(da[k],dof,2,DMDA_STENCIL_BOX,&clone[k]);CHKERRQ(ierr); /* Q2 - but we'll fake it as a Q1 with cells the same size as the Q2 guys */
	}
	
	for (k=0; k<nlevels; k++) {
		ierr = DMGetGlobalVector(clone[k],&properties_A1[k]);CHKERRQ(ierr);  ierr = PetscObjectSetName( (PetscObject)properties_A1[k], "LocalL2ProjQ1_nu");CHKERRQ(ierr);
		ierr = DMGetGlobalVector(clone[k],&properties_A2[k]);CHKERRQ(ierr);  ierr = PetscObjectSetName( (PetscObject)properties_A2[k], "LocalL2ProjQ1_rho");CHKERRQ(ierr);
		ierr = VecZeroEntries(properties_A1[k]);CHKERRQ(ierr);
		ierr = VecZeroEntries(properties_A2[k]);CHKERRQ(ierr);
	}
	ierr = DMGetGlobalVector(clone[nlevels-1],&properties_B);CHKERRQ(ierr);
	ierr = VecZeroEntries(properties_B);CHKERRQ(ierr);
	
	/* compute */
	ierr = _SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_FineGrid(
																																		clone[nlevels-1], properties_A1[nlevels-1],properties_A2[nlevels-1],properties_B, 																																	 
																																		npoints, mp_std,mp_stokes );CHKERRQ(ierr);	

	ierr = _SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints(
																																		clone[nlevels-1], properties_A1[nlevels-1],properties_A2[nlevels-1],
																																		Q[nlevels-1] );CHKERRQ(ierr);
	/* view */
	view = PETSC_FALSE;
	PetscOptionsGetBool(PETSC_NULL,"-view_projected_marker_fields",&view,&flg);
	if (view) {
		PetscViewer viewer;
		
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "SwarmUpdateProperties_LocalL2Proj_Stokes_fine.vtk", &viewer);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = DMView(clone[nlevels-1], viewer);CHKERRQ(ierr);
		ierr = VecView(properties_A1[nlevels-1], viewer);CHKERRQ(ierr);
		ierr = VecView(properties_A2[nlevels-1], viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	}

	ptype = 0;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-mp_hierarchy_projection_type",&ptype,&flg);CHKERRQ(ierr);
	for (k=nlevels-1; k>=1; k--) {
		
		switch (ptype) {
			
			case 0:
			{
				Vec scale;

				/* This introduces scaling effects I need to deal with */
				ierr = DMGetInterpolationScale(clone[k-1],clone[k],R[k],&scale);CHKERRQ(ierr);
				
				ierr = MatRestrict(R[k],properties_A1[k],properties_A1[k-1]);CHKERRQ(ierr);
				ierr = MatRestrict(R[k],properties_A2[k],properties_A2[k-1]);CHKERRQ(ierr);
				
				ierr = VecPointwiseMult(properties_A1[k-1],properties_A1[k-1],scale);CHKERRQ(ierr);
				ierr = VecPointwiseMult(properties_A2[k-1],properties_A2[k-1],scale);CHKERRQ(ierr);
				
				ierr = VecDestroy(&scale);CHKERRQ(ierr);
			}
				break;
				
			case 1:
			{
				VecScatter inject;

				ierr = DMGetInjection(clone[k-1],clone[k],&inject);CHKERRQ(ierr);
				
				ierr = VecScatterBegin(inject,properties_A1[k],properties_A1[k-1],INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
				ierr = VecScatterEnd(inject  ,properties_A1[k],properties_A1[k-1],INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
				
				ierr = VecScatterBegin(inject,properties_A2[k],properties_A2[k-1],INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
				ierr = VecScatterEnd(inject  ,properties_A2[k],properties_A2[k-1],INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
				
				ierr = VecScatterDestroy(&inject);CHKERRQ(ierr);
			}
				break;
				
			default:
				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Must use -mp_hierarchy_projection_type {0,1}"); 
				break;
		}
		
		ierr = _SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints(
																												clone[k-1], properties_A1[k-1],properties_A2[k-1],
																												Q[k-1] );CHKERRQ(ierr);

		if (view) {
			PetscViewer viewer;
			char name[512];
			
			sprintf(name,"SwarmUpdateProperties_LocalL2Proj_Stokes_%d.vtk",k-1);
			ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &viewer);CHKERRQ(ierr);
			ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
			ierr = DMView(clone[k-1], viewer);CHKERRQ(ierr);
			ierr = VecView(properties_A1[k-1], viewer);CHKERRQ(ierr);
			ierr = VecView(properties_A2[k-1], viewer);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
		}
		
	}
	
	
	/* destroy */
	ierr = DMRestoreGlobalVector(clone[nlevels-1],&properties_B);CHKERRQ(ierr);

	for (k=0; k<nlevels; k++) {
		ierr = DMRestoreGlobalVector(clone[k],&properties_A2[k]);CHKERRQ(ierr);
		ierr = DMRestoreGlobalVector(clone[k],&properties_A1[k]);CHKERRQ(ierr);
	}	
	for (k=0; k<nlevels; k++) {
		ierr = DMDestroy(&clone[k]);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


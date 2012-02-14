

#include "petsc.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "material_point_utils.h"
#include "output_paraview.h"


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

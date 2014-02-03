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
 **    Filename:      material_point_utils.c
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

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPEnergy_def.h"
#include "MPntPStokesPl_def.h"

#include "QPntVolCoefStokes_def.h"
#include "QPntVolCoefEnergy_def.h"

#include "dmda_duplicate.h"
#include "dmda_element_q2p1.h"
#include "dmda_element_q1.h"
#include "swarm_fields.h"
#include "output_paraview.h"
#include "quadrature.h"
#include "element_type_Q2.h"
#include "material_point_utils.h"
#include "element_utils_q2.h"
#include "element_utils_q1.h"


#undef __FUNCT__
#define __FUNCT__ "MaterialPointGeneric_VTKWriteBinaryAppendedHeaderAllFields"
PetscErrorCode MaterialPointGeneric_VTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,DataBucket db,int *byte_offset,const int nfields,const MaterialPointField list[])
{
	PetscErrorCode ierr;
	int n,npoints;
	
	PetscFunctionBegin;

	DataBucketGetSizes(db,&npoints,NULL,NULL);
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

				case MPField_Energy:
				{
					DataField   PField_energy;
					MPntPStokes *marker_energy;
					
					DataBucketGetDataFieldByName(db, MPntPEnergy_classname ,&PField_energy);
					DataFieldGetAccess(PField_energy);
					marker_energy = PField_energy->data;
					
					MPntPEnergyVTKWriteBinaryAppendedHeaderAllFields(vtk_fp,byte_offset,(const int)npoints,(const MPntPEnergy*)marker_energy);
					DataFieldRestoreAccess(PField_energy);
				}
					break;
				case MPField_StokesPl:
				{
					DataField     PField_mp_prop;
					MPntPStokesPl *marker_prop;
					
					DataBucketGetDataFieldByName(db, MPntPStokesPl_classname ,&PField_mp_prop);
					DataFieldGetAccess(PField_mp_prop);
					marker_prop = PField_mp_prop->data;
                    
					MPntPStokesPlVTKWriteBinaryAppendedHeaderAllFields(vtk_fp,byte_offset,(const int)npoints,(const MPntPStokesPl*)marker_prop);
					DataFieldRestoreAccess(PField_mp_prop);
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
	
	DataBucketGetSizes(db,&npoints,NULL,NULL);
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

			case MPField_Energy:
			{
				DataField   PField_energy;
				MPntPEnergy *marker_energy;
				
				DataBucketGetDataFieldByName(db, MPntPEnergy_classname ,&PField_energy);
				DataFieldGetAccess(PField_energy);
				marker_energy = PField_energy->data;
				
				MPntPEnergyVTKWriteBinaryAppendedDataAllFields(vtk_fp,(const int)npoints,(const MPntPEnergy*)marker_energy);
				DataFieldRestoreAccess(PField_energy);
			}
				break;
            
			case MPField_StokesPl:
			{
				DataField     PField_mp_prop;
				MPntPStokesPl *marker_prop;
				
				DataBucketGetDataFieldByName(db, MPntPStokesPl_classname ,&PField_mp_prop);
				DataFieldGetAccess(PField_mp_prop);
				marker_prop = PField_mp_prop->data;
				
				MPntPStokesPlVTKWriteBinaryAppendedDataAllFields(vtk_fp,(const int)npoints,(const MPntPStokesPl*)marker_prop);
				DataFieldRestoreAccess(PField_mp_prop);
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
				
			case MPField_Energy:
				MPntPEnergyPVTUWriteAllPPointDataFields(vtk_fp);
				break;
            
			case MPField_StokesPl:
				MPntPStokesPlPVTUWriteAllPPointDataFields(vtk_fp);
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
	ierr = PetscTime(&t0);CHKERRQ(ierr);
	
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
	
	DataBucketGetSizes(db,&npoints,NULL,NULL);
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
	
	ierr = PetscTime(&t1);CHKERRQ(ierr);
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
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&nproc);CHKERRQ(ierr);
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
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
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

	PetscTime(&t0);
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
	PetscTime(&t1);
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
	
	
	PetscTime(&t0);
	ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = VecZeroEntries(Lproperties_A2);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A2,INSERT_VALUES,Lproperties_A2);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A2,INSERT_VALUES,Lproperties_A2);CHKERRQ(ierr);
	PetscTime(&t1);
//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (scatter): %1.4lf ]\n",t1-t0);
	
	PetscTime(&t0);
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	//printf("_SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes NEL %d \n", nel );
	for (e=0;e<nel;e++) {
		ierr = VolumeQuadratureGetCellData_Stokes(Q,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		ierr = DMDAGetScalarElementField(Ae1,nen,(PetscInt*)&elnidx[nen*e],LA_properties_A1);CHKERRQ(ierr);
		ierr = DMDAGetScalarElementField(Ae2,nen,(PetscInt*)&elnidx[nen*e],LA_properties_A2);CHKERRQ(ierr);
				
		/* The Q2 interpolant tends to overshoot, Q1 shouldn't but we check anyway / undershoot when you have viscosity jumps */
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
			cell_gausspoints[p].rho = 0.0;

			cell_gausspoints[p].Fu[0] = 0.0;
			cell_gausspoints[p].Fu[1] = 0.0;
			cell_gausspoints[p].Fu[2] = 0.0;
			cell_gausspoints[p].Fp = 0.0;
			
			for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
				cell_gausspoints[p].eta += NIu[p][i] * Ae1[i];
				cell_gausspoints[p].rho += NIu[p][i] * Ae2[i];
				//cell_gausspoints[p].Fu[1]  += NIu[p][i] * Ae2[i];
			}
			if (cell_gausspoints[p].eta < range_eta[0]) { cell_gausspoints[p].eta = range_eta[0]; }
			if (cell_gausspoints[p].eta > range_eta[1]) { cell_gausspoints[p].eta = range_eta[1]; }
			if (cell_gausspoints[p].rho < range_rho[0]) { cell_gausspoints[p].rho = range_rho[0]; }
			if (cell_gausspoints[p].rho > range_rho[1]) { cell_gausspoints[p].rho = range_rho[1]; }
			
			//printf("e=%d: p=%d: eta = %1.4e: Fu1 = %1.4e \n", e, p, cell_gausspoints[p].eta, cell_gausspoints[p].Fu[1] );		
		}
	}
	
  ierr = VecRestoreArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	PetscTime(&t1);
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
	PetscOptionsGetBool(NULL,"-view_projected_marker_fields",&view,NULL);
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
	PetscScalar NiQ1_p[8];
	PetscScalar Ae1[Q2_NODES_PER_EL_3D], Ae2[Q2_NODES_PER_EL_3D];
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
		
		//P3D_ConstructNi_Q2_3D(xip,NIu[p]);
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
	
	PetscTime(&t0);
	ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = VecZeroEntries(Lproperties_A2);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A2,INSERT_VALUES,Lproperties_A2);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A2,INSERT_VALUES,Lproperties_A2);CHKERRQ(ierr);
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (scatter): %1.4lf ]\n",t1-t0);
	
	PetscTime(&t0);
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	//printf("_SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints NEL %d \n", nel );
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
			cell_gausspoints[p].rho = 0.0;
			
			cell_gausspoints[p].Fu[0] = 0.0;
			cell_gausspoints[p].Fu[1] = 0.0;
			cell_gausspoints[p].Fu[2] = 0.0;
			cell_gausspoints[p].Fp = 0.0;

			for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
				cell_gausspoints[p].eta += NIu[p][i] * Ae1[i];
				cell_gausspoints[p].rho += NIu[p][i] * Ae2[i];
				
				//cell_gausspoints[p].Fu[1]  += NIu[p][i] * Ae2[i];
			}
			if (cell_gausspoints[p].eta < range_eta[0]) { cell_gausspoints[p].eta = range_eta[0]; }
			if (cell_gausspoints[p].eta > range_eta[1]) { cell_gausspoints[p].eta = range_eta[1]; }
			if (cell_gausspoints[p].rho < range_rho[0]) { cell_gausspoints[p].rho = range_rho[0]; }
			if (cell_gausspoints[p].rho > range_rho[1]) { cell_gausspoints[p].rho = range_rho[1]; }

			//printf("e=%d: p=%d: eta = %1.4e: rho = %1.4e \n", e, p, cell_gausspoints[p].eta, cell_gausspoints[p].rho );		
		}
	}
	
  ierr = VecRestoreArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	PetscTime(&t1);
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
	PetscScalar NiQ1_p[8];
	PetscScalar Ni_p[Q2_NODES_PER_EL_3D];
	PetscScalar Ae1[Q2_NODES_PER_EL_3D], Ae2[Q2_NODES_PER_EL_3D], Be[Q2_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A1, Lproperties_A2, Lproperties_B;
	PetscScalar *LA_properties_A1, *LA_properties_A2, *LA_properties_B;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e;
	const PetscInt *elnidx;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_A2);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A2);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_B);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A2,&LA_properties_A2);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_B, &LA_properties_B);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	PetscTime(&t0);
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
	PetscTime(&t1);
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
	
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A2);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy"
PetscErrorCode SwarmUpdateGaussPropertiesLocalL2Projection_Q1_MPntPStokes_Hierarchy(PetscInt coefficient_projection_type,const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,Mat R[],DM da[],Quadrature Q[])
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

	/* 
	 If null projection is chosen, then we assume we have defined quadrature point values on the fine mesh already (by some other means), 
	 thus we do not interpolate the smoothed nodal viscosity onto the quadrature points 
	 */
	if (coefficient_projection_type != -1) {
		ierr = _SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints(
																																			clone[nlevels-1], properties_A1[nlevels-1],properties_A2[nlevels-1],
																																			Q[nlevels-1] );CHKERRQ(ierr);
	}
	
	/* view */
	view = PETSC_FALSE;
	PetscOptionsGetBool(NULL,"-view_projected_marker_fields",&view,&flg);
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
	ierr = PetscOptionsGetInt(NULL,"-mp_hierarchy_projection_type",&ptype,&flg);CHKERRQ(ierr);
	for (k=nlevels-1; k>=1; k--) {
		
		switch (ptype) {
			
			case 0:
			{
				Vec scale;

				/* This introduces scaling effects I need to deal with */
				ierr = DMCreateInterpolationScale(clone[k-1],clone[k],R[k],&scale);CHKERRQ(ierr);
				
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

				ierr = DMCreateInjection(clone[k-1],clone[k],&inject);CHKERRQ(ierr);
				
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

/* generic projection for stokes */
#undef __FUNCT__
#define __FUNCT__ "_compute_memory_offsets"
PetscErrorCode _compute_memory_offsets(void *ref,void *target,size_t *size)
{
	int i;
	size_t len;
	void *stride;
	
	PetscFunctionBegin;
	*size = -1;
	len = 0;
	for (i=0; i<64; i++) {
		stride = (char*)ref + len;
		if (stride == target) {
			*size = len;
		}
		len = len + sizeof(char);
	}
	if ( (*size) == -1 ) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Cannot determine memory offset");
	}	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPntPStokesComputeMemberOffsets"
PetscErrorCode MPntPStokesComputeMemberOffsets(size_t property_offsets[])
{
	MPntPStokes    stokes;
	size_t         s;
	int            i,N;
	static int     been_here = 0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	N = MPntPStokes_nmembers;
	PetscMemzero(property_offsets,sizeof(size_t)*N);
	
	ierr = _compute_memory_offsets(&stokes,&stokes.eta,&s); CHKERRQ(ierr);
	property_offsets[0] = s;
	ierr = _compute_memory_offsets(&stokes,&stokes.rho,&s); CHKERRQ(ierr);
	property_offsets[1] = s;

	if (!been_here) {
		for (i=0; i<N; i++) {
			PetscPrintf(PETSC_COMM_WORLD,"MPntPStokes field offset[%d] %zu \n", i,property_offsets[i]);
		}
		been_here++;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPntPStokesPlComputeMemberOffsets"
PetscErrorCode MPntPStokesPlComputeMemberOffsets(size_t property_offsets[])
{
	MPntPStokesPl  stokespl;
	size_t         s;
	int            i,N;
	static int     been_here = 0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	N = MPntPStokesPl_nmembers;
	PetscMemzero(property_offsets,sizeof(size_t)*N);

	ierr = _compute_memory_offsets(&stokespl,&stokespl.e_plastic,&s); CHKERRQ(ierr);
	property_offsets[0] = s;
	ierr = _compute_memory_offsets(&stokespl,&stokespl.is_yielding,&s); CHKERRQ(ierr);
	property_offsets[1] = s;
	
	if (!been_here) {
		for (i=0; i<N; i++) {
			PetscPrintf(PETSC_COMM_WORLD,"MPntPStokesPl field offset[%d] %zu \n", i,property_offsets[i]);
		}
		been_here++;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MPntPEnergyComputeMemberOffsets"
PetscErrorCode MPntPEnergyComputeMemberOffsets(size_t property_offsets[])
{
	MPntPEnergy    energy;
	size_t         s;
	int            i,N;
	static int     been_here = 0;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	N = MPntPEnergy_nmembers;
	PetscMemzero(property_offsets,sizeof(size_t)*N);

	ierr = _compute_memory_offsets(&energy,&energy.diffusivity,&s); CHKERRQ(ierr);
	property_offsets[0] = s;
	ierr = _compute_memory_offsets(&energy,&energy.heat_source,&s); CHKERRQ(ierr);
	property_offsets[1] = s;

	if (!been_here) {
		for (i=0; i<N; i++) {
			PetscPrintf(PETSC_COMM_WORLD,"MPntPEnergy field offset[%d] %zu \n", i,property_offsets[i]);
		}
		been_here++;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPntVolCoefStokesComputeMemberOffsets"
PetscErrorCode QPntVolCoefStokesComputeMemberOffsets(size_t property_offsets[])
{
	QPntVolCoefStokes stokes;
	size_t            s;
	int               i,N;
	static int        been_here = 0;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	
	N = QPntVolCoefStokes_nmembers;
	PetscMemzero(property_offsets,sizeof(size_t)*N);

	ierr = _compute_memory_offsets(&stokes,&stokes.eta,&s); CHKERRQ(ierr);
	property_offsets[0] = s;
	ierr = _compute_memory_offsets(&stokes,&stokes.rho,&s); CHKERRQ(ierr);
	property_offsets[1] = s;
	ierr = _compute_memory_offsets(&stokes,&stokes.Fu[0],&s); CHKERRQ(ierr);
	property_offsets[2] = s;
	ierr = _compute_memory_offsets(&stokes,&stokes.Fp,&s); CHKERRQ(ierr);
	property_offsets[3] = s;
	
	if (!been_here) {
		for (i=0; i<N; i++) {
			PetscPrintf(PETSC_COMM_WORLD,"QPntVolCoefStokes field offset[%d] %zu \n", i,property_offsets[i]);
		}
		been_here++;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QPntVolCoefEnergyComputeMemberOffsets"
PetscErrorCode QPntVolCoefEnergyComputeMemberOffsets(size_t property_offsets[])
{
	QPntVolCoefEnergy energy;
	size_t            s;
	int               i,N;
	static int        been_here = 0;
	PetscErrorCode    ierr;
	
	PetscFunctionBegin;
	
	N = QPntVolCoefEnergy_nmembers;
	PetscMemzero(property_offsets,sizeof(size_t)*N);

	ierr = _compute_memory_offsets(&energy,&energy.diffusivity,&s); CHKERRQ(ierr);
	property_offsets[0] = s;
	ierr = _compute_memory_offsets(&energy,&energy.heat_source,&s); CHKERRQ(ierr);
	property_offsets[1] = s;
	
	if (!been_here) { 
		for (i=0; i<N; i++) {
			PetscPrintf(PETSC_COMM_WORLD,"QPntVolCoefEnergy field offset[%d] %zu \n", i,property_offsets[i]);
		}
		been_here++;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointQuadraturePointProjectionC0_Q2Stokes"
PetscErrorCode MaterialPointQuadraturePointProjectionC0_Q2Stokes(DM da,DataBucket materialpoint_db,MaterialPointField field,const int member,Quadrature Q)
{
	PetscInt  dof;
	DM        clone;
	Vec       properties_A,properties_B;
	int               npoints;
	DataField         PField_std;
	DataField         PField_material_point_property;
  MPntStd           *mp_std;
	void              *material_point_property;
	size_t            mp_field_offset, mp_offset, qp_field_offset, qp_offset;
	size_t            mp_stokes_property_offsets[MPntPStokes_nmembers];
	size_t            qp_stokes_property_offsets[QPntVolCoefStokes_nmembers];
	PetscBool view;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	
	if (field == MPField_StokesPl) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"MPField_StokesPl cannot be mapped quadrature points");
	}
	
	if (field != MPField_Stokes) {
		/* error - these is only valid for stokes fields defined on Q2 */
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"User must choose either properties which are to be projected onto a Q2 space");
	}

	DataBucketGetDataFieldByName(materialpoint_db, MPntStd_classname,&PField_std);
	DataBucketGetSizes(materialpoint_db,&npoints,NULL,NULL);
	mp_std  = PField_std->data;
	
	
	ierr = MPntPStokesComputeMemberOffsets(mp_stokes_property_offsets);CHKERRQ(ierr);
	ierr = QPntVolCoefStokesComputeMemberOffsets(qp_stokes_property_offsets);CHKERRQ(ierr);
	
	/* setup */
	dof = 1;
	ierr = DMDADuplicateLayout(da,dof,2,DMDA_STENCIL_BOX,&clone);CHKERRQ(ierr);
	
	ierr = DMGetGlobalVector(clone,&properties_A);CHKERRQ(ierr);  
	ierr = DMGetGlobalVector(clone,&properties_B);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(properties_A);CHKERRQ(ierr);
	ierr = VecZeroEntries(properties_B);CHKERRQ(ierr);
	
	
	switch (field) {

		case MPField_Stokes:
		{
			MPntPStokesTypeName stokes_member = (MPntPStokesTypeName)member;
			
			mp_offset = sizeof(MPntPStokes);
			qp_offset = sizeof(QPntVolCoefStokes);
			
			DataBucketGetDataFieldByName(materialpoint_db, MPntPStokes_classname,&PField_material_point_property);
			material_point_property = PField_material_point_property->data;
			
			switch (stokes_member) {
				case MPPStk_eta_effective:
					ierr = PetscObjectSetName( (PetscObject)properties_A, "eta");CHKERRQ(ierr);
					mp_field_offset = mp_stokes_property_offsets[ MPPStk_eta_effective ];
					qp_field_offset = qp_stokes_property_offsets[ QPVCStk_eta_effective ];
				break;
				/* ----------------------------------- */
				case MPPStk_density:
					ierr = PetscObjectSetName( (PetscObject)properties_A, "rho");CHKERRQ(ierr);
					mp_field_offset = mp_stokes_property_offsets[ MPPStk_density ];
					qp_field_offset = qp_stokes_property_offsets[ QPVCStk_rho_effective ];
				break;
				/* ----------------------------------- */
				default:
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"User must choose either {MPPStk_eta_effective, MPPStk_density}");
				break;
			}
		}
		break;

/*
		case MPField_StokesPl:
		{
			MPntPStokesPlTypeName stokespl_member = (MPntPStokesPlTypeName)member;

			DataBucketGetDataFieldByName(materialpoint_db, MPntPStokesPl_classname,&PField_material_point_property);
			material_point_property = PField_material_point_property->data;
			
			switch (stokespl_member) {
				case MPPStkPl_plastic_strain:
					ierr = PetscObjectSetName( (PetscObject)properties_A, "plastic_strain");CHKERRQ(ierr);
					
				break;

				case MPPStkPl_yield_indicator:
					ierr = PetscObjectSetName( (PetscObject)properties_A, "yield_indicator");CHKERRQ(ierr);

				break;
					
				default:
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"User must choose either {MPPStkPl_plastic_strain, MPPStkPl_yield_indicator}");
				break;
			}
		}
		break;
*/			
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"User must choose either {MPField_Stokes, MPField_StokesPl}");
		break;
	}
	
	/* compute */
	//
	ierr = _MaterialPointProjection_MapOntoQ2Mesh(
							clone,properties_A,properties_B,
							//CoefAvgHARMONIC,
							CoefAvgARITHMETIC,
							npoints,mp_std,
							mp_field_offset,mp_offset,material_point_property);CHKERRQ(ierr);
	//
	
	/*
	ierr = _MaterialPointProjection_MapOntoNestedQ1Mesh(
																								clone,properties_A,properties_B,
																								//CoefAvgHARMONIC,
																								CoefAvgARITHMETIC,
																								npoints,mp_std,
																								mp_field_offset,mp_offset,material_point_property);CHKERRQ(ierr);
	*/
	
	/* view */
	view = PETSC_FALSE;
	PetscOptionsGetBool(NULL,"-view_projected_marker_fields",&view,NULL);
	if (view) {
		char filename[256];
		PetscViewer viewer;
		
		sprintf(filename,"MaterialPointProjection_stokes_member_%d.vtk",(int)member );
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = DMView(clone, viewer);CHKERRQ(ierr);
		ierr = VecView(properties_A, viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	}
	
	/* destroy */
	ierr = DMRestoreGlobalVector(clone,&properties_B);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(clone,&properties_A);CHKERRQ(ierr);
	
	ierr = DMDestroy(&clone);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_MaterialPointProjection_MapOntoQ2Mesh"
PetscErrorCode _MaterialPointProjection_MapOntoQ2Mesh(
													 DM clone,Vec properties_A,Vec properties_B,CoefficientAveragingType avg_type,
													 const int npoints,MPntStd mp_std[],
													 size_t member_offset,size_t point_offset,void *point_data) 
{
	PetscScalar NiQ1_p[8];
	PetscScalar Ni_p[Q2_NODES_PER_EL_3D];
	PetscScalar Ae[Q2_NODES_PER_EL_3D], Be[Q2_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A, Lproperties_B;
	PetscScalar *LA_properties_A, *LA_properties_B;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e_p;
	const PetscInt *elnidx;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_B);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_B, &LA_properties_B);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (p=0; p<npoints; p++) {
		double *xi_p;
		void   *point_data_p;
		double field_p;
		
		xi_p = &mp_std[p].xi[0];
		e_p  = mp_std[p].wil;
		
		point_data_p = (void*) ( (char*)point_data + p * point_offset );
		field_p = *( (double*) ( (char*)point_data_p + member_offset) );
		
		if (avg_type == CoefAvgHARMONIC) {
			field_p = 1.0/field_p;
		}
		
		ierr = PetscMemzero(Ae,sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		ierr = PetscMemzero(Be, sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		
		P3D_ConstructNi_Q1_3D(xi_p,NiQ1_p);
		
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
			Ae[i] = Ni_p[i] * field_p;
			Be[i] = Ni_p[i];
		}
		
		/* sum into local vectors */
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e_p]);CHKERRQ(ierr);
		
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_A, 1, el_lidx,Ae);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_B, 1, el_lidx,Be);CHKERRQ(ierr);
		
	}
	PetscTime(&t1);
	//PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (summation): %1.4lf ]\n",t1-t0);
	
  ierr = VecRestoreArray(Lproperties_B,&LA_properties_B);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A,ADD_VALUES,properties_A);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A,ADD_VALUES,properties_A);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	
	/* scale */
	ierr = VecPointwiseDivide( properties_A, properties_A, properties_B );CHKERRQ(ierr);
	/* ========================================= */
	
	if (avg_type == CoefAvgHARMONIC) {
		ierr = VecReciprocal(properties_A);CHKERRQ(ierr);
	}
	
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1_MaterialPointProjection_MapOntoQ2Mesh"
PetscErrorCode DMDAEQ1_MaterialPointProjection_MapOntoQ2Mesh(
																											DM clone,Vec properties_A,Vec properties_B,CoefficientAveragingType avg_type,
																											const int npoints,MPntStd mp_std[],
																											size_t member_offset,size_t point_offset,void *point_data) 
{
	PetscScalar NiQ1_p[Q1_NODES_PER_EL_3D];
	PetscScalar Ae[Q1_NODES_PER_EL_3D], Be[Q1_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A, Lproperties_B;
	PetscScalar *LA_properties_A, *LA_properties_B;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e_p;
	const PetscInt *elnidx;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_B);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_B, &LA_properties_B);CHKERRQ(ierr);
	
	ierr = DMDAGetElementsQ1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (p=0; p<npoints; p++) {
		double *xi_p;
		void   *point_data_p;
		double field_p;
		
		xi_p = &mp_std[p].xi[0];
		e_p  = mp_std[p].wil;
		
		point_data_p = (void*) ( (char*)point_data + p * point_offset );
		field_p = *( (double*) ( (char*)point_data_p + member_offset) );
		
		if (avg_type == CoefAvgHARMONIC) {
			field_p = 1.0/field_p;
		}
		
		ierr = PetscMemzero(Ae,sizeof(PetscScalar)*Q1_NODES_PER_EL_3D);CHKERRQ(ierr);
		ierr = PetscMemzero(Be, sizeof(PetscScalar)*Q1_NODES_PER_EL_3D);CHKERRQ(ierr);
		
		P3D_ConstructNi_Q1_3D(xi_p,NiQ1_p);
		
		for (i=0; i<Q1_NODES_PER_EL_3D; i++) {
			Ae[i] = NiQ1_p[i] * field_p;
			Be[i] = NiQ1_p[i];
		}
		
		/* sum into local vectors */
		ierr = DMDAEQ1_GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e_p]);CHKERRQ(ierr);
		
		ierr = DMDAEQ1_SetValuesLocalStencil_AddValues_DOF(LA_properties_A, 1, el_lidx,Ae);CHKERRQ(ierr);
		ierr = DMDAEQ1_SetValuesLocalStencil_AddValues_DOF(LA_properties_B, 1, el_lidx,Be);CHKERRQ(ierr);
		
	}
	PetscTime(&t1);
	//PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (summation): %1.4lf ]\n",t1-t0);
	
  ierr = VecRestoreArray(Lproperties_B,&LA_properties_B);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A,ADD_VALUES,properties_A);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A,ADD_VALUES,properties_A);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	
	/* scale */
	ierr = VecPointwiseDivide( properties_A, properties_A, properties_B );CHKERRQ(ierr);
	/* ========================================= */
	
	if (avg_type == CoefAvgHARMONIC) {
		ierr = VecReciprocal(properties_A);CHKERRQ(ierr);
	}
	
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "_MaterialPointProjection_MapOntoNestedQ1Mesh"
PetscErrorCode _MaterialPointProjection_MapOntoNestedQ1Mesh(
																											DM clone,Vec properties_A,Vec properties_B,CoefficientAveragingType avg_type,
																											const int npoints,MPntStd mp_std[],
																											size_t member_offset,size_t point_offset,void *point_data) 
{
	PetscScalar NiQ1_p[8];
	PetscScalar Ae[Q2_NODES_PER_EL_3D], Be[Q2_NODES_PER_EL_3D];
	PetscScalar AeQ1[8], BeQ1[8];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A, Lproperties_B;
	PetscScalar *LA_properties_A, *LA_properties_B;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e_p;
	const PetscInt *elnidx;
	PetscInt SIDX,I,J,K;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_B);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_B, &LA_properties_B);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (p=0; p<npoints; p++) {
		double *xi_p,xi_scaled_p[3];
		void   *point_data_p;
		double field_p;
		
		xi_p = &mp_std[p].xi[0];
		e_p  = mp_std[p].wil;
		
		if (xi_p[0] < 0.0) {
			I = 0;
			xi_scaled_p[0] =  2.0 * xi_p[0] + 1.0;
		} else {
			I = 1;
			xi_scaled_p[0] =  2.0 * xi_p[0] - 1.0;
		}

		if (xi_p[1] < 0.0) {
			J = 0;
			xi_scaled_p[1] =  2.0 * xi_p[1] + 1.0;
		} else {
			J = 1;
			xi_scaled_p[1] =  2.0 * xi_p[1] - 1.0;
		}

		if (xi_p[2] < 0.0) {
			K = 0;
			xi_scaled_p[2] =  2.0 * xi_p[2] + 1.0;
		} else {
			K = 1;
			xi_scaled_p[2] =  2.0 * xi_p[2] - 1.0;
		}
		SIDX = I + J * 2 + K * 4;
		
		
		point_data_p = (void*) ( (char*)point_data + p * point_offset );
		field_p = *( (double*) ( (char*)point_data_p + member_offset) );
		
		if (avg_type == CoefAvgHARMONIC) {
			field_p = 1.0/field_p;
		}
		
		P3D_ConstructNi_Q1_3D(xi_scaled_p,NiQ1_p);

		for (i=0; i<8; i++) {
			AeQ1[i] = NiQ1_p[i] * field_p;
			BeQ1[i] = NiQ1_p[i];
		}
		

		{
			PetscInt map[8],ii,jj,kk;
		
			for (kk=0; kk<2; kk++) {
				for (jj=0; jj<2; jj++) {
					for (ii=0; ii<2; ii++) {
						PetscInt sidx = (I + ii) + (J + jj)*3 + (K + kk)*9;
						map[ii+jj*2+kk*4] = sidx;
						//printf("IJK: %d %d %d : sidx = %d \n", I,J,K,sidx);
					}
				}
			}
			
			ierr = PetscMemzero(Ae,sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
			ierr = PetscMemzero(Be, sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
			
			for (i=0; i<8; i++) {
				PetscInt idx = map[i];
				
				Ae[idx] = AeQ1[i];
				Be[idx] = BeQ1[i];
			}
		}
		
		/* sum into local vectors */
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e_p]);CHKERRQ(ierr);
		
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_A, 1, el_lidx,Ae);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_B, 1, el_lidx,Be);CHKERRQ(ierr);
		
	}
	PetscTime(&t1);
	//PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (summation): %1.4lf ]\n",t1-t0);
	
  ierr = VecRestoreArray(Lproperties_B,&LA_properties_B);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A,ADD_VALUES,properties_A);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A,ADD_VALUES,properties_A);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	
	/* scale */
	ierr = VecPointwiseDivide( properties_A, properties_A, properties_B );CHKERRQ(ierr);
	/* ========================================= */
	
	if (avg_type == CoefAvgHARMONIC) {
		ierr = VecReciprocal(properties_A);CHKERRQ(ierr);
	}
	
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_MaterialPointProjection_MapOntoQ2Mesh_InterpolateToQuadraturePoint"
PetscErrorCode _MaterialPointProjection_MapOntoQ2Mesh_InterpolateToQuadraturePoint(
												DM clone,Vec properties_A,
												size_t member_offset,size_t qpoint_offset,void *qpoint_data,Quadrature Q) 
{
	PetscScalar NiQ1_p[8];
	PetscScalar Ae[Q2_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A;
	PetscScalar *LA_properties_A;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e;
	const PetscInt *elnidx;
	
	PetscInt ngp;
	PetscScalar *xi_mp;
	PetscScalar NIu[MAX_QUAD_PNTS][U_BASIS_FUNCTIONS];
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);
	ierr = VecZeroEntries(Lproperties_A);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A,ADD_VALUES,properties_A);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A,ADD_VALUES,properties_A);CHKERRQ(ierr);
	
	/* ========================================= */
	
	/* scatter result back to local array and do the interpolation onto the quadrature points */
	ngp       = Q->npoints;
	xi_mp     = Q->q_xi_coor;
	for (p=0; p<ngp; p++) {
		PetscScalar *xip = &xi_mp[3*p];
		
		ierr = PetscMemzero(NIu[p], sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		
		P3D_ConstructNi_Q1_3D(xip,NiQ1_p);
		NIu[p][0] = NiQ1_p[0];
		NIu[p][2] = NiQ1_p[1];
		NIu[p][6] = NiQ1_p[2];
		NIu[p][8] = NiQ1_p[3];
		
		NIu[p][0+18] = NiQ1_p[4];
		NIu[p][2+18] = NiQ1_p[5];
		NIu[p][6+18] = NiQ1_p[6];
		NIu[p][8+18] = NiQ1_p[7];
	}
	
	PetscTime(&t0);
	ierr = VecZeroEntries(Lproperties_A);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A,INSERT_VALUES,Lproperties_A);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A,INSERT_VALUES,Lproperties_A);CHKERRQ(ierr);
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (scatter): %1.4lf ]\n",t1-t0);
	
	PetscTime(&t0);
	ierr = VecGetArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	for (e=0;e<nel;e++) {
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		ierr = DMDAGetScalarElementField(Ae,nen,(PetscInt*)&elnidx[nen*e],LA_properties_A);CHKERRQ(ierr);
		
		for (p=0; p<ngp; p++) {
			PetscScalar value;
			
			value = 0.0;
			for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
				value += NIu[p][i] * Ae[i];
			}
			
			/* map value into qpoint array */
			*((char*)qpoint_data + ngp*e*qpoint_offset + p*qpoint_offset + member_offset) = value;
			
		}
	}
	
  ierr = VecRestoreArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAEQ1_MaterialPointProjection_MapOntoQ2Mesh_InterpolateToQuadraturePoint"
PetscErrorCode DMDAEQ1_MaterialPointProjection_MapOntoQ2Mesh_InterpolateToQuadraturePoint(
																																									 DM clone,Vec properties_A,
																																									 size_t member_offset,size_t qpoint_offset,void *qpoint_data,Quadrature Q) 
{
	PetscScalar NiQ1_p[Q1_NODES_PER_EL_3D];
	PetscScalar Ae[Q1_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A;
	PetscScalar *LA_properties_A;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e;
	const PetscInt *elnidx;
	
	PetscInt ngp;
	PetscScalar *xi_mp;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	/* scatter result back to local array and do the interpolation onto the quadrature points */
	PetscTime(&t0);

	ierr = DMGetLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);
	ierr = VecZeroEntries(Lproperties_A);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(clone,properties_A,INSERT_VALUES,Lproperties_A);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A,INSERT_VALUES,Lproperties_A);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (scatter): %1.4lf ]\n",t1-t0);
	
	PetscTime(&t0);
	
	/* traverse elements and interpolate */
	ierr = DMDAGetElementsQ1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	ngp       = Q->npoints;
	xi_mp     = Q->q_xi_coor;
	for (e=0; e<nel; e++) {
		ierr = DMDAEQ1_GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		ierr = DMDAEQ1_GetScalarElementField_3D(Ae,(PetscInt*)&elnidx[nen*e],LA_properties_A);CHKERRQ(ierr);
		
		for (p=0; p<ngp; p++) {
			char *ptr;
			PetscScalar *xip = &xi_mp[3*p];
			PetscScalar value;
			
			P3D_ConstructNi_Q1_3D(xip,NiQ1_p);

			value = 0.0;
			for (i=0; i<Q1_NODES_PER_EL_3D; i++) {
				value += NiQ1_p[i] * Ae[i];
			}
			//printf("value = %1.4e \n", value);
			
			/* map value into qpoint array */
			ptr = ((char*)qpoint_data + ngp*e*qpoint_offset + p*qpoint_offset + member_offset);
			//printf("ptr = %p \n",ptr);
			//*ptr = value;
			//*ptr = 1.0;
			ierr = PetscMemcpy(ptr,&value,sizeof(double));CHKERRQ(ierr);
			
			// testing rubbish //
			/*
			{
				QPntVolCoefEnergy *cell;
				ierr = VolumeQuadratureGetCellData_Energy(Q,(QPntVolCoefEnergy*)qpoint_data,e,(QPntVolCoefEnergy**)&cell);CHKERRQ(ierr);				
				//cell[p].diffusivity = 1.0;

				ptr = (char*)cell[p] + member_offset;
				ierr = PetscMemcpy(ptr,&value,sizeof(double));CHKERRQ(ierr);
			}
			*/
		}
	}
	
  ierr = VecRestoreArray(Lproperties_A,&LA_properties_A);CHKERRQ(ierr);
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_A);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* alternative hierarchy constructions */
/*
 THIS FUNCTION LOOKS FUCKING WEIRD... WAS I ON DRUGS OR DRUNK?
*/
#undef __FUNCT__
#define __FUNCT__ "_LocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints"
PetscErrorCode _LocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints(DM clone,Vec properties_A1,Quadrature Q) 
{
	PetscScalar NiQ1_p[8];
	PetscScalar Ae1[Q2_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A1;
	PetscScalar *LA_properties_A1;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e;
	const PetscInt *elnidx;
	
	PetscInt ngp;
	PetscScalar *xi_mp;
	PetscScalar NIu[MAX_QUAD_PNTS][U_BASIS_FUNCTIONS];
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(Q,&all_gausspoints);CHKERRQ(ierr);
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	
	/* ========================================= */
	
	/* scatter result back to local array and do the interpolation onto the quadrature points */
	ngp       = Q->npoints;
	xi_mp     = Q->q_xi_coor;
	for (p=0; p<ngp; p++) {
		PetscScalar *xip = &xi_mp[3*p];
		
		ierr = PetscMemzero(NIu[p], sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
		
		//P3D_ConstructNi_Q2_3D(xip,NIu[p]);
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
	
	PetscTime(&t0);
	ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  clone,properties_A1,INSERT_VALUES,Lproperties_A1);CHKERRQ(ierr);
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (scatter): %1.4lf ]\n",t1-t0);
	
	PetscTime(&t0);
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	/* traverse elements and interpolate */
	//printf("_SwarmUpdateGaussPropertiesLocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints NEL %d \n", nel );
	for (e=0;e<nel;e++) {
		ierr = VolumeQuadratureGetCellData_Stokes(Q,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen*e]);CHKERRQ(ierr);
		
		ierr = DMDAGetScalarElementField(Ae1,nen,(PetscInt*)&elnidx[nen*e],LA_properties_A1);CHKERRQ(ierr);
		
		for (p=0; p<ngp; p++) {
			cell_gausspoints[p].eta = 0.0;
			
			for (i=0; i<Q2_NODES_PER_EL_3D; i++) {
				cell_gausspoints[p].eta    += NIu[p][i] * Ae1[i];
			}
			//printf("eta = %1.4e \n", cell_gausspoints[p].eta);
		}
	}
	
  ierr = VecRestoreArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_LocalL2ProjectionQ1_MPntPStokes"
PetscErrorCode _LocalL2ProjectionQ1_MPntPStokes(	DM clone,Vec properties_A1,Vec properties_B,
																									const PetscInt mx_fine,const PetscInt my_fine,const PetscInt mz_fine,	
																									const PetscInt refx,const PetscInt refy,const PetscInt refz,
																									const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[]) 
{
	PetscScalar NiQ1_p[8];
	PetscScalar Ni_p[Q2_NODES_PER_EL_3D];
	PetscScalar Ae1[Q2_NODES_PER_EL_3D], Be[Q2_NODES_PER_EL_3D];
	PetscInt el_lidx[U_BASIS_FUNCTIONS];
	Vec Lproperties_A1, Lproperties_B;
	PetscScalar *LA_properties_A1, *LA_properties_B;
	PetscLogDouble t0,t1;
	PetscInt p,i;
	PetscInt nel,nen,e,e_level,ei,ej,ek,e2d;
	PetscInt mx_coarse,my_coarse,mz_coarse;
	const PetscInt *elnidx;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
	ierr = DMGetLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_A1);CHKERRQ(ierr);
	ierr = DMGetLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);		ierr = VecZeroEntries(Lproperties_B);CHKERRQ(ierr);
	
	ierr = VecGetArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	ierr = VecGetArray(Lproperties_B, &LA_properties_B);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	ierr = DMDAGetLocalSizeElementQ2(clone,&mx_coarse,&my_coarse,&mz_coarse);CHKERRQ(ierr);
	
	PetscTime(&t0);
	for (p=0; p<npoints; p++) {
		double *xi_p  = &mp_std[p].xi[0];
		double eta_p  = mp_stokes[p].eta;
		
		ierr = PetscMemzero(Ae1,sizeof(PetscScalar)*Q2_NODES_PER_EL_3D);CHKERRQ(ierr);
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
			//Ae1[i] = Ni_p[i] * (1.0/eta_p); /* HARMONIC */
			Be[i]  = Ni_p[i];
		}

		/* element index on fine grid */
		e = mp_std[p].wil;
		/* compute ei,ej,ek on fine grid */
		ek = e/(mx_fine*my_fine);
		e2d = e - ek * (mx_fine*my_fine);
		ej = e2d / mx_fine;
		ei = e2d - ej * mx_fine;
		
		e_level = (ei/refx) + (ej/refy)*mx_coarse + (ek/refz)*mx_coarse*my_coarse;
		
		ierr = Q2GetElementLocalIndicesDOF(el_lidx,1,(PetscInt*)&elnidx[nen * e_level]);CHKERRQ(ierr);
		
		/* sum into local vectors */
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_A1, 1, el_lidx,Ae1);CHKERRQ(ierr);
		ierr = DMDASetValuesLocalStencil_AddValues_DOF(LA_properties_B,  1, el_lidx,Be);CHKERRQ(ierr);
		
	}
	PetscTime(&t1);
	//PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (summation): %1.4lf ]\n",t1-t0);
	
  ierr = VecRestoreArray(Lproperties_B,&LA_properties_B);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lproperties_A1,&LA_properties_A1);CHKERRQ(ierr);
	
	
	/* scatter to quadrature points */
	ierr = DMLocalToGlobalBegin(clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_A1,ADD_VALUES,properties_A1);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  clone,Lproperties_B,ADD_VALUES,properties_B);CHKERRQ(ierr);
	
	/* scale */
	//ierr = VecReciprocal( properties_A1);CHKERRQ(ierr); /* HARMONIC */
	ierr = VecPointwiseDivide( properties_A1, properties_A1, properties_B );CHKERRQ(ierr);
	/* ========================================= */
	
	PetscTime(&t1);
	//	PetscPrintf(PETSC_COMM_WORLD,"  [ L2 projectionQ1 (interpolation): %1.4lf ]\n",t1-t0);
	
	ierr = DMRestoreLocalVector(clone,&Lproperties_B);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(clone,&Lproperties_A1);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MProjection_Q1Projection_onto_Q2_MPntPStokes_Level"
PetscErrorCode MProjection_Q1Projection_onto_Q2_MPntPStokes_Level(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,DM da[],PetscInt level,Quadrature Q_level)
{
	DM clone;
	Vec properties_A,properties_B;
	PetscInt k,refx,refy,refz,REFX[10],REFY[10],REFZ[10];
	PetscInt mx_fine,my_fine,mz_fine,dof;
	PetscBool view;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* setup */
	dof = 1;
	ierr = DMDADuplicateLayout(da[level],dof,2,DMDA_STENCIL_BOX,&clone);CHKERRQ(ierr);
	ierr = DMDASetElementType_Q2(clone);CHKERRQ(ierr);
	
	ierr = DMGetGlobalVector(clone,&properties_A);CHKERRQ(ierr);  
	ierr = DMGetGlobalVector(clone,&properties_B);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(properties_A);CHKERRQ(ierr);
	ierr = VecZeroEntries(properties_B);CHKERRQ(ierr);

	
	/* compute refinement factor for this level */
	for (k=0; k<nlevels; k++) {
		ierr = DMDAGetRefinementFactor(da[k],&REFX[k],&REFY[k],&REFZ[k]);CHKERRQ(ierr);
	}
	
	refx = 1;
	refy = 1;
	refz = 1;
	for (k=nlevels-1; k>=level+1; k--) {
		refx = refx * REFX[k];
		refy = refy * REFY[k];
		refz = refz * REFZ[k];
	}
	
	ierr = DMDAGetLocalSizeElementQ2(da[nlevels-1],&mx_fine,&my_fine,&mz_fine);CHKERRQ(ierr);
	
	
	ierr = _LocalL2ProjectionQ1_MPntPStokes(	clone,properties_A,properties_B,
																						mx_fine,my_fine,mz_fine,	
																						refx,refy,refz,
																					  npoints,mp_std,mp_stokes);CHKERRQ(ierr); 
	ierr = _LocalL2ProjectionQ1_MPntPStokes_InterpolateToQuadratePoints(clone,properties_A,Q_level);CHKERRQ(ierr);

	/* view */
	view = PETSC_FALSE;
	PetscOptionsGetBool(NULL,"-view_projected_marker_fields",&view,NULL);
	if (view) {
		char filename[256];
		PetscViewer viewer;
		
		sprintf(filename,"MProjectionQ1_stokes_eta_Lv%d.vtk",level );
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = DMView(clone, viewer);CHKERRQ(ierr);
		ierr = VecView(properties_A, viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	}
	
	/* destroy */
	ierr = DMRestoreGlobalVector(clone,&properties_B);CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(clone,&properties_A);CHKERRQ(ierr);
	ierr = DMDestroy(&clone);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


/* P0 */
#undef __FUNCT__
#define __FUNCT__ "_LocalP0Projection_MPntPStokes_MapToQuadratePoints"
PetscErrorCode _LocalP0Projection_MPntPStokes_MapToQuadratePoints( DM clone,
																																	const PetscInt mx_fine,const PetscInt my_fine,const PetscInt mz_fine,	
																																	const PetscInt refx,const PetscInt refy,const PetscInt refz,
																																	const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],Quadrature Q) 
{
	PetscLogDouble t0,t1;
	PetscInt p;
	const PetscInt *elnidx;
	
	PetscInt ngp;
	QPntVolCoefStokes *all_gausspoints,*cell_gausspoints;
	PetscInt nel,nen,e,e_level,ei,ej,ek,e2d;
	PetscInt mx_coarse,my_coarse,mz_coarse;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	ierr = DMDAGetElements_pTatinQ2P1(clone,&nel,&nen,&elnidx);CHKERRQ(ierr);
	ierr = DMDAGetLocalSizeElementQ2(clone,&mx_coarse,&my_coarse,&mz_coarse);CHKERRQ(ierr);
	
	ierr = VolumeQuadratureGetAllCellData_Stokes(Q,&all_gausspoints);CHKERRQ(ierr);
	
	/* init eta to zero interpolate */

	ngp       = Q->npoints;
	
	for (e=0;e<nel;e++) {
		ierr = VolumeQuadratureGetCellData_Stokes(Q,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);
		for (p=0; p<ngp; p++) {
			cell_gausspoints[p].eta = 0.0;
		}
	}
	
	/* traverse elements and interpolate */
	PetscTime(&t0);
	for (p=0; p<npoints; p++) {
		double eta_p  = mp_stokes[p].eta;
		
		/* element index on fine grid */
		e = mp_std[p].wil;
		/* compute ei,ej,ek on fine grid */
		ek = e/(mx_fine*my_fine);
		e2d = e - ek * (mx_fine*my_fine);
		ej = e2d / mx_fine;
		ei = e2d - ej * mx_fine;
		
		e_level = (ei/refx) + (ej/refy)*mx_coarse + (ek/refz)*mx_coarse*my_coarse;
		
		ierr = VolumeQuadratureGetCellData_Stokes(Q,all_gausspoints,e_level,&cell_gausspoints);CHKERRQ(ierr);

		cell_gausspoints[0].eta += eta_p; /* ARTHMATIC */
		//cell_gausspoints[0].eta += 1.0/eta_p; /* HARMONIC */
		cell_gausspoints[1].eta += 1.0;
	}
	PetscTime(&t1);
	
	/* set constant value over all quadrature points */
	for (e=0;e<nel;e++) {
		double sum_eta, sum_np;
		
		ierr = VolumeQuadratureGetCellData_Stokes(Q,all_gausspoints,e,&cell_gausspoints);CHKERRQ(ierr);

		sum_eta = cell_gausspoints[0].eta;
		sum_np  = cell_gausspoints[1].eta;
		
		for (p=0; p<ngp; p++) {
			cell_gausspoints[p].eta = sum_eta / sum_np; /* ARTHMATIC */
			//cell_gausspoints[p].eta = 1.0/( sum_np * sum_eta); /* HARMONIC */
		}
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MProjection_P0Projection_onto_Q2_MPntPStokes_Level"
PetscErrorCode MProjection_P0Projection_onto_Q2_MPntPStokes_Level(const int npoints,MPntStd mp_std[],MPntPStokes mp_stokes[],PetscInt nlevels,DM da[],PetscInt level,Quadrature Q_level)
{
	DM clone;
	PetscInt k,refx,refy,refz,REFX[10],REFY[10],REFZ[10];
	PetscInt mx_fine,my_fine,mz_fine,dof;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* setup */
	dof = 1;
	ierr = DMDADuplicateLayout(da[level],dof,2,DMDA_STENCIL_BOX,&clone);CHKERRQ(ierr);
	ierr = DMDASetElementType_Q2(clone);CHKERRQ(ierr);
	
	/* compute refinement factor for this level */
	for (k=0; k<nlevels; k++) {
		ierr = DMDAGetRefinementFactor(da[k],&REFX[k],&REFY[k],&REFZ[k]);CHKERRQ(ierr);
	}
	
	refx = 1;
	refy = 1;
	refz = 1;
	for (k=nlevels-1; k>=level+1; k--) {
		refx = refx * REFX[k];
		refy = refy * REFY[k];
		refz = refz * REFZ[k];
	}
	
	ierr = DMDAGetLocalSizeElementQ2(da[nlevels-1],&mx_fine,&my_fine,&mz_fine);CHKERRQ(ierr);

	ierr = _LocalP0Projection_MPntPStokes_MapToQuadratePoints(clone, mx_fine,my_fine,mz_fine,	refx,refy,refz, npoints,mp_std,mp_stokes,Q_level);CHKERRQ(ierr);
	
	/* destroy */
	ierr = DMDestroy(&clone);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}



/* MATERIAL POINT ACCESS HELPERS */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGetAccess"
PetscErrorCode MaterialPointGetAccess(DataBucket materialpoint_db,MPAccess *helper)
{
	MPAccess   X;
	int        Lfields,cnt;
	BTruth     found;
	DataField  PField;
	PetscErrorCode ierr;	
	
	PetscFunctionBegin;
	
	ierr = PetscMalloc(sizeof(struct _p_MPAccess),&X);CHKERRQ(ierr);
	ierr = PetscMemzero(X,sizeof(struct _p_MPAccess));CHKERRQ(ierr);
	
	X->db = materialpoint_db;
	
	DataBucketGetDataFields(materialpoint_db,&Lfields,NULL);
	ierr = PetscMalloc(sizeof(DataField)*Lfields,&X->PField);CHKERRQ(ierr);
	ierr = PetscMemzero(X->PField,sizeof(DataField)*Lfields);CHKERRQ(ierr);
	
	cnt = 0;
	
	/* USER: add reference to all possible material point types here */
	
	/* init all idx for material point fields */
	X->mp_std_field_idx      = -1;
	X->mp_stokes_field_idx   = -1;
	X->mp_stokespl_field_idx = -1;
	X->mp_energy_field_idx   = -1;
	
	/* MPntStd */
	DataBucketQueryDataFieldByName(materialpoint_db,MPntStd_classname,&found);
	if (found) {
		DataBucketGetDataFieldByName(materialpoint_db,MPntStd_classname,&PField);
		DataFieldGetAccess(PField);
		DataFieldVerifyAccess(PField,sizeof(MPntStd));
		
		X->PField[cnt] = PField;
		X->mp_std_field_idx = cnt;
		
		cnt++;
	}
	
	/* MPntPStokes */
	DataBucketQueryDataFieldByName(materialpoint_db,MPntPStokes_classname,&found);
	if (found) {
		DataBucketGetDataFieldByName(materialpoint_db,MPntPStokes_classname,&PField);
		DataFieldGetAccess(PField);
		DataFieldVerifyAccess(PField,sizeof(MPntPStokes));
		
		X->PField[cnt] = PField;
		X->mp_stokes_field_idx = cnt;
		
		cnt++;
	}
	
	/* MPntPStokesPl */
	DataBucketQueryDataFieldByName(materialpoint_db,MPntPStokesPl_classname,&found);
	if (found) {
		DataBucketGetDataFieldByName(materialpoint_db,MPntPStokesPl_classname,&PField);
		DataFieldGetAccess(PField);
		DataFieldVerifyAccess(PField,sizeof(MPntPStokesPl));
		
		X->PField[cnt] = PField;
		X->mp_stokespl_field_idx = cnt;
		
		cnt++;
	}
	
	/* MPntPEnergy_classname */
	DataBucketQueryDataFieldByName(materialpoint_db,MPntPEnergy_classname,&found);
	if (found) {
		DataBucketGetDataFieldByName(materialpoint_db,MPntPEnergy_classname,&PField);
		DataFieldGetAccess(PField);
		DataFieldVerifyAccess(PField,sizeof(MPntPEnergy));
		
		X->PField[cnt] = PField;
		X->mp_energy_field_idx = cnt;
		
		cnt++;
	}
	
	X->nfields = cnt;
	
	*helper = X;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointRestoreAccess"
PetscErrorCode MaterialPointRestoreAccess(DataBucket matpoint_db,MPAccess *helper)
{
	MPAccess X;
	int      i;
	PetscErrorCode ierr;	
	
	PetscFunctionBegin;
	
	X = *helper;
	
	for (i=0; i<X->nfields; i++) {
		DataFieldRestoreAccess(X->PField[i]);
	}
	ierr = PetscFree(X->PField);CHKERRQ(ierr);
	ierr = PetscFree(X);CHKERRQ(ierr);
	
	*helper = NULL;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_get_field_MPntStd"
PetscErrorCode _get_field_MPntStd(MPAccess X,const int p,MPntStd **point)
{
	DataField  PField;
	if (X->mp_std_field_idx == -1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Material point field MPntStd must be registered");
	}
	if (X == PETSC_NULL) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must call MaterialPointGetAccess() first"); }
	PField = X->PField[ X->mp_std_field_idx ];
	DataFieldAccessPoint(PField,p,(void**)point);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "_get_field_MPntPStokes"
PetscErrorCode _get_field_MPntPStokes(MPAccess X,const int p,MPntPStokes **point)
{
	DataField  PField;
	if (X->mp_stokes_field_idx == -1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Material point field MPntPStokes must be registered");
	}
	if (X == PETSC_NULL) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must call MaterialPointGetAccess() first"); }
	PField = X->PField[ X->mp_stokes_field_idx ];
	DataFieldAccessPoint(PField,p,(void**)point);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "_get_field_MPntPStokesPl"
PetscErrorCode _get_field_MPntPStokesPl(MPAccess X,const int p,MPntPStokesPl **point)
{
	DataField  PField;
	if (X->mp_stokespl_field_idx == -1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Material point field MPntPStokesPl must be registered");
	}
	if (X == PETSC_NULL) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must call MaterialPointGetAccess() first"); }
	PField = X->PField[ X->mp_stokespl_field_idx ];
	DataFieldAccessPoint(PField,p,(void**)point);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "_get_field_MPntPEnergy"
PetscErrorCode _get_field_MPntPEnergy(MPAccess X,const int p,MPntPEnergy **point)
{
	DataField  PField;
	if (X->mp_energy_field_idx == -1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Material point field MPntPEnergy must be registered");
	}
	if (X == PETSC_NULL) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must call MaterialPointGetAccess() first"); }
	PField = X->PField[ X->mp_energy_field_idx ];
	DataFieldAccessPoint(PField,p,(void**)point);
	
	PetscFunctionReturn(0);
}

/* std */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_global_coord"
PetscErrorCode MaterialPointGet_global_coord(MPAccess X,const int p,double *var[])
{
	MPntStd    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntStd(X,p,&point);CHKERRQ(ierr);
	MPntStdGetField_global_coord(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_local_coord"
PetscErrorCode MaterialPointGet_local_coord(MPAccess X,const int p,double *var[])
{
	MPntStd    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntStd(X,p,&point);CHKERRQ(ierr);
	MPntStdGetField_local_coord(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_local_element_index"
PetscErrorCode MaterialPointGet_local_element_index(MPAccess X,const int p,int *var)
{
	MPntStd    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntStd(X,p,&point);CHKERRQ(ierr);
	MPntStdGetField_local_element_index(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_phase_index"
PetscErrorCode MaterialPointGet_phase_index(MPAccess X,const int p,int *var)
{
	MPntStd    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntStd(X,p,&point);CHKERRQ(ierr);
	MPntStdGetField_phase_index(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_phase_index"
PetscErrorCode MaterialPointSet_phase_index(MPAccess X,const int p,int var)
{
	MPntStd    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntStd(X,p,&point);CHKERRQ(ierr);
	MPntStdSetField_phase_index(point,var);
	
	PetscFunctionReturn(0);
}

/* stokes */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_viscosity"
PetscErrorCode MaterialPointGet_viscosity(MPAccess X,const int p,double *var)
{
	MPntPStokes    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPStokes(X,p,&point);CHKERRQ(ierr);
	MPntPStokesGetField_eta_effective(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_viscosity"
PetscErrorCode MaterialPointSet_viscosity(MPAccess X,const int p,double var)
{
	MPntPStokes    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPStokes(X,p,&point);CHKERRQ(ierr);
	MPntPStokesSetField_eta_effective(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_density"
PetscErrorCode MaterialPointGet_density(MPAccess X,const int p,double *var)
{
	MPntPStokes    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPStokes(X,p,&point);CHKERRQ(ierr);
	MPntPStokesGetField_density(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_density"
PetscErrorCode MaterialPointSet_density(MPAccess X,const int p,double var)
{
	MPntPStokes    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPStokes(X,p,&point);CHKERRQ(ierr);
	MPntPStokesSetField_density(point,var);
	
	PetscFunctionReturn(0);
}
/* stokespl */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_plastic_strain"
PetscErrorCode MaterialPointGet_plastic_strain(MPAccess X,const int p,float *var)
{
	MPntPStokesPl  *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPStokesPl(X,p,&point);CHKERRQ(ierr);
	MPntPStokesPlGetField_plastic_strain(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_plastic_strain"
PetscErrorCode MaterialPointSet_plastic_strain(MPAccess X,const int p,float var)
{
	MPntPStokesPl  *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPStokesPl(X,p,&point);CHKERRQ(ierr);
	MPntPStokesPlSetField_plastic_strain(point,var);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_yield_indicator"
PetscErrorCode MaterialPointGet_yield_indicator(MPAccess X,const int p,short *var)
{
	MPntPStokesPl  *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPStokesPl(X,p,&point);CHKERRQ(ierr);
	MPntPStokesPlGetField_yield_indicator(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_yield_indicator"
PetscErrorCode MaterialPointSet_yield_indicator(MPAccess X,const int p,short var)
{
	MPntPStokesPl  *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPStokesPl(X,p,&point);CHKERRQ(ierr);
	MPntPStokesPlSetField_yield_indicator(point,var);
	
	PetscFunctionReturn(0);
}

/* energy */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_diffusivity"
PetscErrorCode MaterialPointGet_diffusivity(MPAccess X,const int p,double *var)
{
	MPntPEnergy    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPEnergy(X,p,&point);CHKERRQ(ierr);
	MPntPEnergyGetField_diffusivity(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointGet_heat_source"
PetscErrorCode MaterialPointGet_heat_source(MPAccess X,const int p,double *var)
{
	MPntPEnergy    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPEnergy(X,p,&point);CHKERRQ(ierr);
	MPntPEnergyGetField_heat_source(point,var);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_diffusivity"
PetscErrorCode MaterialPointSet_diffusivity(MPAccess X,const int p,double var)
{
	MPntPEnergy    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPEnergy(X,p,&point);CHKERRQ(ierr);
	MPntPEnergySetField_diffusivity(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_heat_source"
PetscErrorCode MaterialPointSet_heat_source(MPAccess X,const int p,double var)
{
	MPntPEnergy    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntPEnergy(X,p,&point);CHKERRQ(ierr);
	MPntPEnergySetField_heat_source(point,var);
	
	PetscFunctionReturn(0);
}
/* std */
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_global_coord"
PetscErrorCode MaterialPointSet_global_coord(MPAccess X,const int p,double var[])
{
	MPntStd    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntStd(X,p,&point);CHKERRQ(ierr);
	MPntStdSetField_global_coord(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_local_coord"
PetscErrorCode MaterialPointSet_local_coord(MPAccess X,const int p,double var[])
{
	MPntStd    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntStd(X,p,&point);CHKERRQ(ierr);
	MPntStdSetField_local_coord(point,var);
	
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "MaterialPointSet_local_element_index"
PetscErrorCode MaterialPointSet_local_element_index(MPAccess X,const int p,int var)
{
	MPntStd    *point;
	PetscErrorCode ierr;	
	PetscFunctionBegin;
	
	ierr = _get_field_MPntStd(X,p,&point);CHKERRQ(ierr);
	MPntStdSetField_local_element_index(point,var);
	
	PetscFunctionReturn(0);
}


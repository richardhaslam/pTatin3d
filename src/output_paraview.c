/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   output_paraview.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/


#include "petsc.h"
#include "dmda_element_q2p1.h"
#include "output_paraview.h"

const char *PTatinFieldNames[] =  {
  "pressure", 
  "viscosity", 
  "s_xx","s_yy","s_zz","s_xy","s_xz","s_yz",
  "t_xx","t_yy","t_zz","t_xy","t_xz","t_yz",
  "e_xx","e_yy","e_zz","e_xy","e_xz","e_yz",
  "s_II",
  "t_II",
  "e_II",
  0 };



#undef __FUNCT__  
#define __FUNCT__ "pTatinGenerateParallelVTKName"
PetscErrorCode pTatinGenerateParallelVTKName(const char prefix[],const char suffix[],char **name)
{
	char *nn;
	int rank;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	if (prefix!=NULL) {
		if (asprintf(&nn,"%s-subdomain%1.5d.%s",prefix,rank,suffix) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
    *name = nn;
	} else {
    *name = NULL;
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"You must specify a prefix");
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinGenerateVTKName"
PetscErrorCode pTatinGenerateVTKName(const char prefix[],const char suffix[],char **name)
{
	char *nn;

	PetscFunctionBegin;
  *name = NULL;
	if (prefix!=NULL) {
		if (asprintf(&nn,"%s.%s",prefix,suffix) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	} else {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"You must specify a prefix");
	}
	*name = nn;
	PetscFunctionReturn(0);
}

/* PVD helpers */
#undef __FUNCT__
#define __FUNCT__ "ParaviewPVDOpen"
PetscErrorCode ParaviewPVDOpen(const char pvdfilename[])
{
	PetscMPIInt rank;
	FILE        *fp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* only master generates this file */
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	if(rank != 0) { PetscFunctionReturn(0); }
	
	fp = fopen(pvdfilename,"w");
  if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new PVD file %s",pvdfilename);
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
#ifdef WORDSIZE_BIGENDIAN
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	fprintf(fp,"<Collection>\n");
	
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>");
	fclose(fp);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ParaviewPVDAppend"
PetscErrorCode ParaviewPVDAppend(const char pvdfilename[],double time,const char datafile[], const char DirectoryName[])
{
	PetscMPIInt rank;
	FILE        *fp;
	char        line[10000];
	int         key_L,position;
	char        key[] = "</Collection>";
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* only master generates this file */
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	if (rank != 0) { PetscFunctionReturn(0); }
	
	fp = fopen(pvdfilename,"r+");
  if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open existing PVD file %s",pvdfilename);
	/* reset to start of file */
	rewind(fp);
	
	key_L = strlen(key);
	position = -1;
	while (!feof(fp)) {
		position = ftell(fp);
		if (!fgets(line,10000-1,fp)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_READ,"fgets() failed");
		if (strncmp(key,line,key_L) == 0) {
			break;
		}
	}
	if (position == -1) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot locate keyword in pvd file"); }
	
	fseek(fp,position,SEEK_SET);
	
	/* write new data */
	if (DirectoryName == NULL) {
		fprintf(fp,"  <DataSet timestep=\"%1.6e\" file=\"./%s\"/>\n",time, datafile );
	} else if ( (strlen(DirectoryName) == 0) ) {
		fprintf(fp,"  <DataSet timestep=\"%1.6e\" file=\"./%s\"/>\n",time, datafile );
	} else {
		fprintf(fp,"  <DataSet timestep=\"%1.6e\" file=\"%s/%s\"/>\n",time, DirectoryName, datafile );
	}
	
	/* close tag */
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>");
	
	fclose(fp);
	
	PetscFunctionReturn(0);
}


/* V-P mesh */
#undef __FUNCT__  
#define __FUNCT__ "pTatinOutputParaViewMeshVelocityPressure"
PetscErrorCode pTatinOutputParaViewMeshVelocityPressure(DM pack,Vec X,const char path[],const char prefix[])
{
	char           *vtkfilename,*filename;
	PetscMPIInt    rank;
	PetscBool      binary = PETSC_TRUE;
	PetscBool      zip = PETSC_FALSE;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = pTatinGenerateParallelVTKName(prefix,"vts",&vtkfilename);CHKERRQ(ierr);
	if (path) {
		if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	} else {
		if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	}
	
	if (binary) {
		if (zip) {
			char zfilename[PETSC_MAX_PATH_LEN];
			
			sprintf(zfilename,"%s.gz",filename);
			ierr = pTatinOutputMeshVelocityPressureVTS_v0_binary_gz(pack,X,zfilename);CHKERRQ(ierr);
			
		} else {
			ierr = pTatinOutputMeshVelocityPressureVTS_v0_binary(pack,X,filename);CHKERRQ(ierr);
		}
		
	} else {
		ierr = pTatinOutputMeshVelocityPressureVTS_v0(pack,X,filename);CHKERRQ(ierr);
	}
	
	free(filename);
	free(vtkfilename);
	
	ierr = pTatinGenerateVTKName(prefix,"pvts",&vtkfilename);CHKERRQ(ierr);
	if (path) {
		if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	} else {
		if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	}
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = pTatinOutputMeshVelocityPressurePVTS(pack,prefix,filename);CHKERRQ(ierr);
	free(filename);
	free(vtkfilename);
	
	PetscFunctionReturn(0);
}

/* V-P mesh */
#undef __FUNCT__  
#define __FUNCT__ "pTatinOutputLiteParaViewMeshVelocity"
PetscErrorCode pTatinOutputLiteParaViewMeshVelocity(DM pack,Vec X,const char path[],const char prefix[])
{
	char           *vtkfilename,*filename;
	PetscMPIInt    rank;
	PetscBool      binary = PETSC_TRUE;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = pTatinGenerateParallelVTKName(prefix,"vts",&vtkfilename);CHKERRQ(ierr);
	if (path) {
		if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	} else {
		if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	}
	
	if (binary) {
		ierr = pTatinOutputLiteMeshVelocityVTS_v0_binary(pack,X,filename);CHKERRQ(ierr);
	} else {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only binary format implemented");
	}
	
	free(filename);
	free(vtkfilename);
	
	ierr = pTatinGenerateVTKName(prefix,"pvts",&vtkfilename);CHKERRQ(ierr);
	if (path) {
		if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	} else {
		if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
	}
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = pTatinOutputLiteMeshVelocityPVTS(pack,prefix,filename);CHKERRQ(ierr);
	free(filename);
	free(vtkfilename);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputMeshVelocityPressureVTS_v0"
PetscErrorCode pTatinOutputMeshVelocityPressureVTS_v0(DM pack,Vec X,const char name[])
{
	PetscErrorCode ierr;
	DM dau,dap;
	Vec velocity,pressure;
  Vec local_fieldsU;
  PetscScalar *LA_fieldsU;
  Vec local_fieldsP;
  PetscScalar *LA_fieldsP;
	DM cda;
	Vec gcoords;
	DMDACoor3d ***LA_gcoords;	
	PetscInt mx,my,mz,cnt;
	PetscInt ei,ej,ek,i,j,k,esi,esj,esk;
	FILE*	vtk_fp = NULL;
	PetscInt gsi,gsj,gsk,gm,gn,gp;
	PetscInt ndof_pressure;

	
	PetscFunctionBegin;
	if ((vtk_fp = fopen ( name, "w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new VTS file %s",name );
	}
	
	//PetscPrintf(PETSC_COMM_WORLD,"[[DESIGN FLAW]] %s: only printing P0 component of pressure field \n", __FUNCT__ );

	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dap,0,0,0,0,0,0,0,&ndof_pressure,0, 0,0,0, 0);CHKERRQ(ierr);
	
	ierr = DMDAGetGhostCorners(dau,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
	ierr = DMDAGetCornersElementQ2(dau,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&local_fieldsU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dau,velocity,INSERT_VALUES,local_fieldsU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dau,velocity,INSERT_VALUES,local_fieldsU);CHKERRQ(ierr);
  ierr = VecGetArray(local_fieldsU,&LA_fieldsU);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dap,&local_fieldsP);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dap,pressure,INSERT_VALUES,local_fieldsP);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dap,pressure,INSERT_VALUES,local_fieldsP);CHKERRQ(ierr);
  ierr = VecGetArray(local_fieldsP,&LA_fieldsP);CHKERRQ(ierr);
	
	
	/* VTS HEADER - OPEN */	
#ifdef WORDSIZE_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "    <Piece Extent=\"%D %D %D %D %D %D\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	
	/* VTS COORD DATA */	
	fprintf( vtk_fp, "    <Points>\n");
	fprintf( vtk_fp, "      <DataArray Name=\"coords\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for (k=esk; k<esk+2*mz+1; k++) {
		for (j=esj; j<esj+2*my+1; j++) {
			for (i=esi; i<esi+2*mx+1; i++) {
				fprintf( vtk_fp,"      %1.6e %1.6e %1.6e\n", LA_gcoords[k][j][i].x, LA_gcoords[k][j][i].y, LA_gcoords[k][j][i].z );
			}
		}
	}
	fprintf( vtk_fp, "      </DataArray>\n");
	fprintf( vtk_fp, "    </Points>\n");
	
	/* VTS CELL DATA */	
	fprintf( vtk_fp, "    <CellData>\n");
	
	/* integrated pressure */
	//ierr = CellIntegration(LA_cell_data,PField_Pressure,&user->rheology_constants,ctx, dau,LA_fieldsU, dap,LA_fieldsP);CHKERRQ(ierr);
	
	fprintf( vtk_fp, "      <DataArray Name=\"pressure0\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	fprintf( vtk_fp,"      ");
	for (ek=0; ek<mz; ek++) { for (k=0; k<2; k++) {
		for (ej=0; ej<my; ej++) { for (j=0; j<2; j++) {
			for (ei=0; ei<mx; ei++) { for (i=0; i<2; i++) {
				/* ONLY print the P0 component of pressure */
				PetscScalar pressure_0 = LA_fieldsP[ ndof_pressure * (ei + ej*mx + ek*mx*my) + 0 ];
				
				fprintf( vtk_fp,"%1.6e ", pressure_0 );
			}}
		}}
	}}
	fprintf( vtk_fp,"\n");
	fprintf( vtk_fp, "      </DataArray>\n");
	
	fprintf( vtk_fp, "    </CellData>\n");
	
	/* VTS NODAL DATA */
	fprintf( vtk_fp, "    <PointData>\n");
	/* velocity */
	fprintf( vtk_fp, "      <DataArray Name=\"velocity\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for (k=esk; k<esk+2*mz+1; k++) {
		for (j=esj; j<esj+2*my+1; j++) {
			for (i=esi; i<esi+2*mx+1; i++) {
				cnt = (i-gsi) + (j-gsj)*gm + (k-gsk)*gm*gn;
				fprintf( vtk_fp,"      %1.6e %1.6e %1.6e\n", LA_fieldsU[3*cnt+0], LA_fieldsU[3*cnt+1], LA_fieldsU[3*cnt+2] );
			}
		}
	}
	fprintf( vtk_fp, "      </DataArray>\n");
	fprintf( vtk_fp, "    </PointData>\n");
	
	/* VTS HEADER - CLOSE */	
	fprintf( vtk_fp, "    </Piece>\n");
	fprintf( vtk_fp, "  </StructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");
	
  ierr = VecRestoreArray(local_fieldsU,&LA_fieldsU);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(dau,&local_fieldsU);CHKERRQ(ierr);
	
  ierr = VecRestoreArray(local_fieldsP,&LA_fieldsP);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(dap,&local_fieldsP);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	fclose( vtk_fp );
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputMeshVelocityPressureVTS_v0_binary"
PetscErrorCode pTatinOutputMeshVelocityPressureVTS_v0_binary(DM pack,Vec X,const char name[])
{
	PetscErrorCode ierr;
	DM dau,dap;
	Vec velocity,pressure;
  Vec local_fieldsU;
  PetscScalar *LA_fieldsU;
  Vec local_fieldsP;
  PetscScalar *LA_fieldsP;
	DM cda;
	Vec gcoords;
	DMDACoor3d ***LA_gcoords;	
	PetscInt mx,my,mz,cnt;
	PetscInt ei,ej,ek,i,j,k,esi,esj,esk;
	FILE*	vtk_fp = NULL;
	PetscInt gsi,gsj,gsk,gm,gn,gp;
	PetscInt ndof_pressure;
	int offset,bytes;
	
	
	PetscFunctionBegin;
	if ((vtk_fp = fopen ( name, "w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new VTS file %s",name );
	}
	
	//PetscPrintf(PETSC_COMM_WORLD,"[[DESIGN FLAW]] %s: only printing P0 component of pressure field \n", __FUNCT__ );
	
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dap,0,0,0,0,0,0,0,&ndof_pressure,0, 0,0,0, 0);CHKERRQ(ierr);
	
	ierr = DMDAGetGhostCorners(dau,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
	ierr = DMDAGetCornersElementQ2(dau,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&local_fieldsU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dau,velocity,INSERT_VALUES,local_fieldsU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dau,velocity,INSERT_VALUES,local_fieldsU);CHKERRQ(ierr);
  ierr = VecGetArray(local_fieldsU,&LA_fieldsU);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dap,&local_fieldsP);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dap,pressure,INSERT_VALUES,local_fieldsP);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dap,pressure,INSERT_VALUES,local_fieldsP);CHKERRQ(ierr);
  ierr = VecGetArray(local_fieldsP,&LA_fieldsP);CHKERRQ(ierr);
	
	
	/* VTS HEADER - OPEN */	
#ifdef WORDSIZE_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "    <Piece Extent=\"%D %D %D %D %D %D\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	
	offset = 0;
	
	/* VTS COORD DATA */	
	fprintf( vtk_fp, "    <Points>\n");
	
	fprintf( vtk_fp, "      <DataArray Name=\"coords\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"  />\n",offset);
	offset += sizeof(int) + sizeof(double)*3*(2*mx+1)*(2*my+1)*(2*mz+1);

	fprintf( vtk_fp, "    </Points>\n");
	
	/* VTS CELL DATA */	
	fprintf( vtk_fp, "    <CellData>\n");
	
	fprintf( vtk_fp, "      <DataArray Name=\"pressure0\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\" />\n",offset);
	offset += sizeof(int) + sizeof(double)*1*(2*mx)*(2*my)*(2*mz);
	
	fprintf( vtk_fp, "    </CellData>\n");
	
	/* VTS NODAL DATA */
	fprintf( vtk_fp, "    <PointData>\n");
	
	/* velocity */
	fprintf( vtk_fp, "      <DataArray Name=\"velocity\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",offset);
	offset += sizeof(int) + sizeof(double)*3*(2*mx+1)*(2*my+1)*(2*mz+1);

	fprintf( vtk_fp, "    </PointData>\n");
	
	/* VTS HEADER - CLOSE */	
	fprintf( vtk_fp, "    </Piece>\n");
	fprintf( vtk_fp, "  </StructuredGrid>\n");
	fprintf(vtk_fp, "  <AppendedData encoding=\"raw\">\n");
	/* write tag */
	fprintf(vtk_fp, "_");

	/* write node coords */
	bytes = sizeof(double)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	fwrite(&bytes,sizeof(int),1,vtk_fp);
	for (k=esk; k<esk+2*mz+1; k++) {
		for (j=esj; j<esj+2*my+1; j++) {
			for (i=esi; i<esi+2*mx+1; i++) {
				double pos[3];
				
				pos[0] = LA_gcoords[k][j][i].x;
				pos[1] = LA_gcoords[k][j][i].y;
				pos[2] = LA_gcoords[k][j][i].z;
				fwrite(pos,sizeof(double),3,vtk_fp);
			}
		}
	}
	
	/* write cell pressure */
	bytes = sizeof(double)*1*(2*mx+1)*(2*my+1)*(2*mz+1);
	fwrite(&bytes,sizeof(int),1,vtk_fp);
	for (ek=0; ek<mz; ek++) { for (k=0; k<2; k++) {
		for (ej=0; ej<my; ej++) { for (j=0; j<2; j++) {
			for (ei=0; ei<mx; ei++) { for (i=0; i<2; i++) {
				double P0;
				
				/* ONLY print the P0 component of pressure */
				P0 = LA_fieldsP[ ndof_pressure * (ei + ej*mx + ek*mx*my) + 0 ];
				
				fwrite(&P0,sizeof(double),1,vtk_fp);
			}}
		}}
	}}
	
	/* write node velocity */
	bytes = sizeof(double)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	fwrite(&bytes,sizeof(int),1,vtk_fp);
	for (k=esk; k<esk+2*mz+1; k++) {
		for (j=esj; j<esj+2*my+1; j++) {
			for (i=esi; i<esi+2*mx+1; i++) {
				double vel[3];
				
				cnt = (i-gsi) + (j-gsj)*gm + (k-gsk)*gm*gn;
				vel[0] = LA_fieldsU[3*cnt+0];
				vel[1] = LA_fieldsU[3*cnt+1];
				vel[2] = LA_fieldsU[3*cnt+2];
				fwrite(vel,sizeof(double),3,vtk_fp);
			}
		}
	}
	
	fprintf(vtk_fp, "\n  </AppendedData>\n");
	fprintf( vtk_fp, "</VTKFile>\n");
	
  ierr = VecRestoreArray(local_fieldsU,&LA_fieldsU);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(dau,&local_fieldsU);CHKERRQ(ierr);
	
  ierr = VecRestoreArray(local_fieldsP,&LA_fieldsP);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(dap,&local_fieldsP);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	fclose( vtk_fp );
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputLiteMeshVelocityVTS_v0_binary"
PetscErrorCode pTatinOutputLiteMeshVelocityVTS_v0_binary(DM pack,Vec X,const char name[])
{
	PetscErrorCode ierr;
	DM dau,dap;
	Vec velocity,pressure;
  Vec local_fieldsU;
  PetscScalar *LA_fieldsU;
	DM cda;
	Vec gcoords;
	DMDACoor3d ***LA_gcoords;	
	PetscInt mx,my,mz,cnt;
	PetscInt i,j,k,esi,esj,esk;
	FILE*	vtk_fp = NULL;
	PetscInt gsi,gsj,gsk,gm,gn,gp;
	int offset,bytes;
	
	
	PetscFunctionBegin;
	if ((vtk_fp = fopen ( name, "w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new VTS file %s",name );
	}
	
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
		
	ierr = DMDAGetGhostCorners(dau,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
	ierr = DMDAGetCornersElementQ2(dau,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&local_fieldsU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dau,velocity,INSERT_VALUES,local_fieldsU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dau,velocity,INSERT_VALUES,local_fieldsU);CHKERRQ(ierr);
  ierr = VecGetArray(local_fieldsU,&LA_fieldsU);CHKERRQ(ierr);
	
	/* VTS HEADER - OPEN */	
#ifdef WORDSIZE_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "    <Piece Extent=\"%D %D %D %D %D %D\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	
	offset = 0;
	
	/* VTS COORD DATA */	
	fprintf( vtk_fp, "    <Points>\n");
	
	fprintf( vtk_fp, "      <DataArray Name=\"coords\" type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"  />\n",offset);
	offset += sizeof(int) + sizeof(float)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	
	fprintf( vtk_fp, "    </Points>\n");
	
	/* VTS CELL DATA */	
	fprintf( vtk_fp, "    <CellData>\n");
	fprintf( vtk_fp, "    </CellData>\n");
	
	/* VTS NODAL DATA */
	fprintf( vtk_fp, "    <PointData>\n");
	
	/* velocity */
	fprintf( vtk_fp, "      <DataArray Name=\"velocity\" type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",offset);
	offset += sizeof(int) + sizeof(float)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	
	fprintf( vtk_fp, "    </PointData>\n");
	
	/* VTS HEADER - CLOSE */	
	fprintf( vtk_fp, "    </Piece>\n");
	fprintf( vtk_fp, "  </StructuredGrid>\n");
	fprintf(vtk_fp, "  <AppendedData encoding=\"raw\">\n");
	/* write tag */
	fprintf(vtk_fp, "_");
	
	/* write node coords */
	bytes = sizeof(float)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	fwrite(&bytes,sizeof(int),1,vtk_fp);
	for (k=esk; k<esk+2*mz+1; k++) {
		for (j=esj; j<esj+2*my+1; j++) {
			for (i=esi; i<esi+2*mx+1; i++) {
				float pos[3];
				
				pos[0] = (float)LA_gcoords[k][j][i].x;
				pos[1] = (float)LA_gcoords[k][j][i].y;
				pos[2] = (float)LA_gcoords[k][j][i].z;
				fwrite(pos,sizeof(float),3,vtk_fp);
			}
		}
	}
	
	/* write node velocity */
	bytes = sizeof(float)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	fwrite(&bytes,sizeof(int),1,vtk_fp);
	for (k=esk; k<esk+2*mz+1; k++) {
		for (j=esj; j<esj+2*my+1; j++) {
			for (i=esi; i<esi+2*mx+1; i++) {
				float vel[3];
				
				cnt = (i-gsi) + (j-gsj)*gm + (k-gsk)*gm*gn;
				vel[0] = (float)LA_fieldsU[3*cnt+0];
				vel[1] = (float)LA_fieldsU[3*cnt+1];
				vel[2] = (float)LA_fieldsU[3*cnt+2];
				fwrite(vel,sizeof(float),3,vtk_fp);
			}
		}
	}
	
	fprintf(vtk_fp, "\n  </AppendedData>\n");
	fprintf( vtk_fp, "</VTKFile>\n");
	
  ierr = VecRestoreArray(local_fieldsU,&LA_fieldsU);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(dau,&local_fieldsU);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	fclose( vtk_fp );
	PetscFunctionReturn(0);
}

#include "zlib.h"

#undef __FUNCT__
#define __FUNCT__ "report_gzlib_error"
PetscErrorCode report_gzlib_error(int ierr_gz,gzFile fp)
{
	int        errnum;
	const char *info;

	PetscFunctionBegin;
	info = gzerror(fp,&errnum);
	if (ierr_gz < 0) {
        PetscPrintf(PETSC_COMM_WORLD,"gzerror: %s \n",info);
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"ERROR(gzprintf): Failed to write data");
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputMeshVelocityPressureVTS_v0_binary_gz"
PetscErrorCode pTatinOutputMeshVelocityPressureVTS_v0_binary_gz(DM pack,Vec X,const char name[])
{	
	PetscErrorCode ierr;
	DM dau,dap;
	Vec velocity,pressure;
  Vec local_fieldsU;
  PetscScalar *LA_fieldsU;
  Vec local_fieldsP;
  PetscScalar *LA_fieldsP;
	DM cda;
	Vec gcoords;
	DMDACoor3d ***LA_gcoords;	
	PetscInt mx,my,mz,cnt;
	PetscInt ei,ej,ek,i,j,k,esi,esj,esk;
	gzFile vtk_fp = NULL;
	PetscInt gsi,gsj,gsk,gm,gn,gp;
	PetscInt ndof_pressure;
	int offset,bytes;
	int ierr_gz;
	
	PetscFunctionBegin;
	vtk_fp = gzopen ( name, "wb");
	if (vtk_fp == NULL) {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"gzerror: Cannot open file %s",name );
	}
	
	//PetscPrintf(PETSC_COMM_WORLD,"[[DESIGN FLAW]] %s: only printing P0 component of pressure field \n", __FUNCT__ );
	
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dap,0,0,0,0,0,0,0,&ndof_pressure,0, 0,0,0, 0);CHKERRQ(ierr);
	
	ierr = DMDAGetGhostCorners(dau,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
	ierr = DMDAGetCornersElementQ2(dau,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dau,&local_fieldsU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dau,velocity,INSERT_VALUES,local_fieldsU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dau,velocity,INSERT_VALUES,local_fieldsU);CHKERRQ(ierr);
  ierr = VecGetArray(local_fieldsU,&LA_fieldsU);CHKERRQ(ierr);
	
  ierr = DMGetLocalVector(dap,&local_fieldsP);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dap,pressure,INSERT_VALUES,local_fieldsP);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dap,pressure,INSERT_VALUES,local_fieldsP);CHKERRQ(ierr);
  ierr = VecGetArray(local_fieldsP,&LA_fieldsP);CHKERRQ(ierr);
	
	
	/* VTS HEADER - OPEN */	
#ifdef WORDSIZE_BIGENDIAN
	ierr_gz = gzprintf(vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	ierr_gz = gzprintf(vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	report_gzlib_error(ierr_gz,vtk_fp);

	ierr_gz = gzprintf(vtk_fp, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	report_gzlib_error(ierr_gz,vtk_fp);
	ierr_gz = gzprintf(vtk_fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	report_gzlib_error(ierr_gz,vtk_fp);
	
	offset = 0;
	
	/* VTS COORD DATA */	
	ierr_gz = gzprintf(vtk_fp, "    <Points>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	
	ierr_gz = gzprintf(vtk_fp, "      <DataArray Name=\"coords\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"  />\n",offset);
	report_gzlib_error(ierr_gz,vtk_fp);
	offset += sizeof(int) + sizeof(double)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	
	ierr_gz = gzprintf(vtk_fp, "    </Points>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	
	/* VTS CELL DATA */	
	ierr_gz = gzprintf(vtk_fp, "    <CellData>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	
	ierr_gz = gzprintf(vtk_fp, "      <DataArray Name=\"pressure0\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\" />\n",offset);
	report_gzlib_error(ierr_gz,vtk_fp);
	offset += sizeof(int) + sizeof(double)*1*(2*mx)*(2*my)*(2*mz);
	
	ierr_gz = gzprintf(vtk_fp, "    </CellData>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	
	/* VTS NODAL DATA */
	ierr_gz = gzprintf(vtk_fp, "    <PointData>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	
	/* velocity */
	ierr_gz = gzprintf(vtk_fp, "      <DataArray Name=\"velocity\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",offset);
	report_gzlib_error(ierr_gz,vtk_fp);
	offset += sizeof(int) + sizeof(double)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	
	ierr_gz = gzprintf(vtk_fp, "    </PointData>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	
	/* VTS HEADER - CLOSE */	
	ierr_gz = gzprintf(vtk_fp, "    </Piece>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	ierr_gz = gzprintf(vtk_fp, "  </StructuredGrid>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	ierr_gz = gzprintf(vtk_fp, "  <AppendedData encoding=\"raw\">\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	/* write tag */
	ierr_gz = gzprintf(vtk_fp, "_");
	report_gzlib_error(ierr_gz,vtk_fp);
	
	/* write node coords */
	bytes = sizeof(double)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	gzwrite(vtk_fp,&bytes,sizeof(int)*1);
	for (k=esk; k<esk+2*mz+1; k++) {
		for (j=esj; j<esj+2*my+1; j++) {
			for (i=esi; i<esi+2*mx+1; i++) {
				double pos[3];
				
				pos[0] = LA_gcoords[k][j][i].x;
				pos[1] = LA_gcoords[k][j][i].y;
				pos[2] = LA_gcoords[k][j][i].z;
				ierr = gzwrite(vtk_fp,(void*)pos,sizeof(double)*3);
			}
		}
	}
	
	/* write cell pressure */
	bytes = sizeof(double)*1*(2*mx+1)*(2*my+1)*(2*mz+1);
	gzwrite(vtk_fp,&bytes,sizeof(int)*1);
	for (ek=0; ek<mz; ek++) { for (k=0; k<2; k++) {
		for (ej=0; ej<my; ej++) { for (j=0; j<2; j++) {
			for (ei=0; ei<mx; ei++) { for (i=0; i<2; i++) {
				double P0;
				
				/* ONLY print the P0 component of pressure */
				P0 = LA_fieldsP[ ndof_pressure * (ei + ej*mx + ek*mx*my) + 0 ];
				
				ierr = gzwrite(vtk_fp,(void*)&P0,sizeof(double)*1);
			}}
		}}
	}}
	
	/* write node velocity */
	bytes = sizeof(double)*3*(2*mx+1)*(2*my+1)*(2*mz+1);
	gzwrite(vtk_fp,&bytes,sizeof(int)*1);
	for (k=esk; k<esk+2*mz+1; k++) {
		for (j=esj; j<esj+2*my+1; j++) {
			for (i=esi; i<esi+2*mx+1; i++) {
				double vel[3];
				
				cnt = (i-gsi) + (j-gsj)*gm + (k-gsk)*gm*gn;
				vel[0] = LA_fieldsU[3*cnt+0];
				vel[1] = LA_fieldsU[3*cnt+1];
				vel[2] = LA_fieldsU[3*cnt+2];
				ierr = gzwrite(vtk_fp,(void*)vel,sizeof(double)*3);
			}
		}
	}
	
	ierr_gz = gzprintf(vtk_fp, "\n  </AppendedData>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	ierr_gz = gzprintf(vtk_fp, "</VTKFile>\n");
	report_gzlib_error(ierr_gz,vtk_fp);
	
  ierr = VecRestoreArray(local_fieldsU,&LA_fieldsU);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(dau,&local_fieldsU);CHKERRQ(ierr);
	
  ierr = VecRestoreArray(local_fieldsP,&LA_fieldsP);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(dap,&local_fieldsP);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(pack,X,&velocity,&pressure);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr_gz = gzclose(vtk_fp);
	report_gzlib_error(ierr_gz,vtk_fp);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputMeshVelocityPressurePVTS"
PetscErrorCode pTatinOutputMeshVelocityPressurePVTS(DM pack,const char prefix[],const char name[])
{
	PetscErrorCode ierr;
	DM             dau,dap;
	FILE           *vtk_fp = NULL;
	PetscInt       M,N,P,swidth;
	PetscMPIInt    rank;
	
	PetscFunctionBegin;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	vtk_fp = NULL;
	if (rank==0) {
		if ((vtk_fp = fopen ( name, "w")) == NULL)  {
			SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new PVTS file %s",name );
		}
	}
	
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	/* VTS HEADER - OPEN */	
	if(vtk_fp) fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

#ifdef WORDSIZE_BIGENDIAN
	if(vtk_fp) fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	if(vtk_fp) fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	DMDAGetInfo( dau, 0, &M,&N,&P, 0,0,0, 0,&swidth, 0,0,0, 0 );
	if(vtk_fp) PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <PStructuredGrid GhostLevel=\"%D\" WholeExtent=\"%D %D %D %D %D %D\">\n", swidth, 0,M-1, 0,N-1, 0,P-1 ); /* note overlap = 1 for Q1 */
	//	fprintf( vtk_fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, 0,0);
	
	/* VTS COORD DATA */	
	if(vtk_fp) fprintf( vtk_fp, "    <PPoints>\n");
	if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PPoints>\n");
	
	
	/* VTS CELL DATA */	
	if(vtk_fp) fprintf( vtk_fp, "    <PCellData>\n");
	if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"Float64\" Name=\"pressure0\" NumberOfComponents=\"1\"/>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PCellData>\n");
	
	
	/* VTS NODAL DATA */
	if(vtk_fp) fprintf( vtk_fp, "    <PPointData>\n");
	if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PPointData>\n");
	
	/* write out the parallel information */
	//DAViewVTK_write_PieceExtend(vtk_fp,2,dau,prefix);
	ierr = DAQ2PieceExtendForGhostLevelZero(vtk_fp,2,dau,prefix);CHKERRQ(ierr);
	
	
	/* VTS HEADER - CLOSE */	
	if(vtk_fp) fprintf( vtk_fp, "  </PStructuredGrid>\n");
	if(vtk_fp) fprintf( vtk_fp, "</VTKFile>\n");
	
	if(vtk_fp) fclose( vtk_fp );
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputLiteMeshVelocityPVTS"
PetscErrorCode pTatinOutputLiteMeshVelocityPVTS(DM pack,const char prefix[],const char name[])
{
	PetscErrorCode ierr;
	DM             dau,dap;
	FILE           *vtk_fp = NULL;
	PetscInt       M,N,P,swidth;
	PetscMPIInt    rank;
	
	PetscFunctionBegin;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	vtk_fp = NULL;
	if (rank==0) {
		if ((vtk_fp = fopen ( name, "w")) == NULL)  {
			SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new PVTS file %s",name );
		}
	}
	
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	/* VTS HEADER - OPEN */	
	if(vtk_fp) fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
	
#ifdef WORDSIZE_BIGENDIAN
	if(vtk_fp) fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	if(vtk_fp) fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	ierr = DMDAGetInfo( dau, 0, &M,&N,&P, 0,0,0, 0,&swidth, 0,0,0, 0 );CHKERRQ(ierr);
	if(vtk_fp) PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <PStructuredGrid GhostLevel=\"%D\" WholeExtent=\"%D %D %D %D %D %D\">\n", swidth, 0,M-1, 0,N-1, 0,P-1 ); /* note overlap = 1 for Q1 */
	
	/* VTS COORD DATA */	
	if(vtk_fp) fprintf( vtk_fp, "    <PPoints>\n");
	if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PPoints>\n");
	
	
	/* VTS CELL DATA */	
	if(vtk_fp) fprintf( vtk_fp, "    <PCellData>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PCellData>\n");
	
	
	/* VTS NODAL DATA */
	if(vtk_fp) fprintf( vtk_fp, "    <PPointData>\n");
	if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PPointData>\n");
	
	/* write out the parallel information */
	ierr = DAQ2PieceExtendForGhostLevelZero(vtk_fp,2,dau,prefix);CHKERRQ(ierr);
	
	/* VTS HEADER - CLOSE */	
	if(vtk_fp) fprintf( vtk_fp, "  </PStructuredGrid>\n");
	if(vtk_fp) fprintf( vtk_fp, "</VTKFile>\n");
	
	if(vtk_fp) fclose( vtk_fp );
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAQ2PieceExtendForGhostLevelZero"
PetscErrorCode DAQ2PieceExtendForGhostLevelZero( FILE *vtk_fp, int indent_level, DM dau, const char local_file_prefix[] )
{
	PetscMPIInt nproc,rank;
	MPI_Comm comm;
	PetscInt M,N,P,pM,pN,pP;
	PetscInt i,j,k,II,dim;
	PetscInt *olx,*oly,*olz;
	PetscInt *lmx,*lmy,*lmz;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* create file name */
	PetscObjectGetComm( (PetscObject)dau, &comm );
	ierr = MPI_Comm_size( comm, &nproc );CHKERRQ(ierr);
	ierr = MPI_Comm_rank( comm, &rank );CHKERRQ(ierr);
	
	ierr = DMDAGetInfo( dau, &dim, &M,&N,&P, &pM,&pN,&pP, 0, 0, 0,0,0, 0 );CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRangesElementQ2(dau,&pM,&pN,&pP,&olx,&oly,&olz,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	
	if (dim==3) {
		for( k=0;k<pP;k++ ) {
			for( j=0;j<pN;j++ ) {
				for( i=0;i<pM;i++ ) {
					char *name;
					PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
					int procid32;
					
					PetscMPIIntCast(procid,&procid32);
					if (asprintf( &name, "%s-subdomain%1.5d.vts", local_file_prefix, procid32 ) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
					for( II=0; II<indent_level; II++ ) {
						if(vtk_fp) fprintf(vtk_fp,"  ");
					}
					if(vtk_fp) PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "<Piece Extent=\"%D %D %D %D %D %D\"      Source=\"%s\"/>\n",
														 olx[i],olx[i]+lmx[i]*2,
														 oly[j],oly[j]+lmy[j]*2,
														 olz[k],olz[k]+lmz[k]*2,
														 name);
					free(name);
				}
			}
		}
	} else if (dim==2) {
		for( j=0;j<pN;j++ ) {
			for( i=0;i<pM;i++ ) {
				char *name;
				PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
				int procid32;
				
				PetscMPIIntCast(procid,&procid32);
				if (asprintf( &name, "%s-subdomain%1.5d.vts", local_file_prefix, procid32 ) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
				for( II=0; II<indent_level; II++ ) {
					if(vtk_fp) fprintf(vtk_fp,"  ");
				}
				if(vtk_fp) PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "<Piece Extent=\"%D %D %D %D 0 0\"      Source=\"%s\"/>\n",
													 olx[i],olx[i]+lmx[i]*2,
													 oly[j],oly[j]+lmy[j]*2,
													 name);
				free(name);
			}
		}
	}
	ierr = PetscFree(olx);CHKERRQ(ierr);
	ierr = PetscFree(oly);CHKERRQ(ierr);
	ierr = PetscFree(olz);CHKERRQ(ierr);
	
	ierr = PetscFree(lmx);CHKERRQ(ierr);
	ierr = PetscFree(lmy);CHKERRQ(ierr);
	ierr = PetscFree(lmz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "_pTatinOutputLiteMeshVelocitySlicedPVTS"
PetscErrorCode _pTatinOutputLiteMeshVelocitySlicedPVTS(FILE *vtk_fp,
                                     PetscInt processor_span[],
                                     PetscInt olx[],PetscInt lmx[],
                                     PetscInt oly[],PetscInt lmy[],
                                     PetscInt olz[],PetscInt lmz[],
                                     PetscInt processor_I[],
                                     PetscInt processor_J[],
                                     PetscInt processor_K[],
                                     const char prefix[])
{
	PetscInt       swidth,i,j,k;
	PetscInt       startI,endI,startJ,endJ,startK,endK;
	
	PetscFunctionBegin;
	
	/* VTS HEADER - OPEN */	
	if(vtk_fp) fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
	
#ifdef WORDSIZE_BIGENDIAN
	if(vtk_fp) fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	if(vtk_fp) fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	swidth = 2;
	
	startI = olx[processor_I[0]];
	endI   = olx[processor_I[1]-1] + 2*lmx[processor_I[1]-1];

	startJ = oly[processor_J[0]];
	endJ   = oly[processor_J[1]-1] + 2*lmy[processor_J[1]-1];

	startK = olz[processor_K[0]];
	endK   = olz[processor_K[1]-1] + 2*lmz[processor_K[1]-1];
	
	if(vtk_fp) PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <PStructuredGrid GhostLevel=\"%D\" WholeExtent=\"%D %D %D %D %D %D\">\n", swidth, startI,endI, startJ,endJ, startK,endK );
	
	/* VTS COORD DATA */	
	if(vtk_fp) fprintf( vtk_fp, "    <PPoints>\n");
	if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PPoints>\n");
	
	
	/* VTS CELL DATA */	
	if(vtk_fp) fprintf( vtk_fp, "    <PCellData>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PCellData>\n");
	
	
	/* VTS NODAL DATA */
	if(vtk_fp) fprintf( vtk_fp, "    <PPointData>\n");
	if(vtk_fp) fprintf( vtk_fp, "      <PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n");
	if(vtk_fp) fprintf( vtk_fp, "    </PPointData>\n");
	
	/* write out the parallel information */
	for (k=processor_K[0]; k<processor_K[1]; k++) {
		for (j=processor_J[0]; j<processor_J[1]; j++) {
			for (i=processor_I[0]; i<processor_I[1]; i++) {
				char     *name;
				PetscInt procid;
                int procid32;
				
				procid = i + j*processor_span[0] + k*processor_span[0]*processor_span[1]; /* convert proc(i,j,k) to pid */
				PetscMPIIntCast(procid,&procid32);

				if (asprintf( &name, "%s-subdomain%1.5d.vts", prefix, procid32 ) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
				if(vtk_fp) PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "      <Piece Extent=\"%D %D %D %D %D %D\"      Source=\"%s\"/>\n",
													 olx[i],olx[i]+lmx[i]*2,
													 oly[j],oly[j]+lmy[j]*2,
													 olz[k],olz[k]+lmz[k]*2,
													 name);
				free(name);
			}
		}
	}
	
	/* VTS HEADER - CLOSE */	
	if(vtk_fp) fprintf( vtk_fp, "  </PStructuredGrid>\n");
	if(vtk_fp) fprintf( vtk_fp, "</VTKFile>\n");
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputLiteMeshVelocitySlicedPVTS"
PetscErrorCode pTatinOutputLiteMeshVelocitySlicedPVTS(DM pack,const char path[],const char prefix[])
{
	PetscErrorCode ierr;
	DM             dau,dap;
	FILE           *vtk_fp = NULL;
	PetscMPIInt    rank;
	PetscInt pM,pN,pP;
	PetscInt i,j,k;
	PetscInt *olx,*oly,*olz;
	PetscInt *lmx,*lmy,*lmz;
	PetscInt processor_span[3];
	PetscInt processor_I[2];
	PetscInt processor_J[2];
	PetscInt processor_K[2];
	
	PetscFunctionBegin;
	
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRangesElementQ2(dau,&pM,&pN,&pP,&olx,&oly,&olz,&lmx,&lmy,&lmz);CHKERRQ(ierr);

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

	processor_span[0] = pM;
	processor_span[1] = pN;
	processor_span[2] = pP;
	
	if (rank == 0) {
		
		// I
		for (i=0; i<pM; i++) {
			char pprefix[PETSC_MAX_PATH_LEN];
			char sdprefix[PETSC_MAX_PATH_LEN];
			char *vtkfilename,*filename;
			
			PetscSNPrintf(pprefix,PETSC_MAX_PATH_LEN-1,"%s_v_PSliceI%D",prefix,i);
			ierr = pTatinGenerateVTKName(pprefix,"pvts",&vtkfilename);CHKERRQ(ierr);
			if (path) {
				if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
			} else {
				if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
			}

			vtk_fp = NULL;
			if ((vtk_fp = fopen ( filename, "w")) == NULL)  {
				SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new sliced PVTS file %s",filename );
			}
			
			processor_I[0] = i;
			processor_I[1] = i+1;
			
			processor_J[0] = 0;
			processor_J[1] = pN;
			processor_K[0] = 0;
			processor_K[1] = pP;

			sprintf(sdprefix,"%s_v",prefix);
			_pTatinOutputLiteMeshVelocitySlicedPVTS(vtk_fp,
																							processor_span,
																							olx,lmx,oly,lmy,olz,lmz,
																							processor_I,processor_J,processor_K,
																							sdprefix);
			free(filename);
			free(vtkfilename);
			fclose(vtk_fp);
		}

		// J
		for (j=0; j<pN; j++) {
			char pprefix[PETSC_MAX_PATH_LEN];
			char sdprefix[PETSC_MAX_PATH_LEN];
			char *vtkfilename,*filename;
			
			PetscSNPrintf(pprefix,PETSC_MAX_PATH_LEN-1,"%s_v_PSliceJ%D",prefix,j);
			ierr = pTatinGenerateVTKName(pprefix,"pvts",&vtkfilename);CHKERRQ(ierr);
			if (path) {
				if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
			} else {
				if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
			}
			
			vtk_fp = NULL;
			if ((vtk_fp = fopen ( filename, "w")) == NULL)  {
				SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new sliced PVTS file %s",filename );
			}
			
			processor_I[0] = 0;
			processor_I[1] = pM;
			
			processor_J[0] = j;
			processor_J[1] = j+1;
			processor_K[0] = 0;
			processor_K[1] = pP;
			
			sprintf(sdprefix,"%s_v",prefix);
			_pTatinOutputLiteMeshVelocitySlicedPVTS(vtk_fp,
																							processor_span,
																							olx,lmx,oly,lmy,olz,lmz,
																							processor_I,processor_J,processor_K,
																							sdprefix);
			free(filename);
			free(vtkfilename);
			fclose(vtk_fp);
		}

		// K
		for (k=0; k<pP; k++) {
			char pprefix[PETSC_MAX_PATH_LEN];
			char sdprefix[PETSC_MAX_PATH_LEN];
			char *vtkfilename,*filename;
			
			PetscSNPrintf(pprefix,PETSC_MAX_PATH_LEN-1,"%s_v_PSliceK%D",prefix,k);
			ierr = pTatinGenerateVTKName(pprefix,"pvts",&vtkfilename);CHKERRQ(ierr);
			if (path) {
				if (asprintf(&filename,"%s/%s",path,vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
			} else {
				if (asprintf(&filename,"./%s",vtkfilename) < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"asprintf() failed");
			}
			
			vtk_fp = NULL;
			if ((vtk_fp = fopen ( filename, "w")) == NULL)  {
				SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open new sliced PVTS file %s",filename );
			}
			
			processor_I[0] = 0;
			processor_I[1] = pM;
			
			processor_J[0] = 0;
			processor_J[1] = pN;
			processor_K[0] = k;
			processor_K[1] = k+1;
			
			sprintf(sdprefix,"%s_v",prefix);
			_pTatinOutputLiteMeshVelocitySlicedPVTS(vtk_fp,
																							processor_span,
																							olx,lmx,oly,lmy,olz,lmz,
																							processor_I,processor_J,processor_K,
																							sdprefix);
			free(filename);
			free(vtkfilename);
			fclose(vtk_fp);
		}
		
		
	}
		
	ierr = PetscFree(olx);CHKERRQ(ierr);
	ierr = PetscFree(oly);CHKERRQ(ierr);
	ierr = PetscFree(olz);CHKERRQ(ierr);
	
	ierr = PetscFree(lmx);CHKERRQ(ierr);
	ierr = PetscFree(lmy);CHKERRQ(ierr);
	ierr = PetscFree(lmz);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

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
 **    Filename:      output_paraview.c
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
	
	PetscFunctionBegin;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	if (prefix!=NULL) {
		asprintf(&nn,"%s-subdomain%1.5d.%s",prefix,rank,suffix);
	} else {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"You must specify a prefix");
	}
	*name = nn;
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinGenerateVTKName"
PetscErrorCode pTatinGenerateVTKName(const char prefix[],const char suffix[],char **name)
{
	char *nn;

	PetscFunctionBegin;
	if (prefix!=NULL) {
		asprintf(&nn,"%s.%s",prefix,suffix);
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
	FILE *fp;
	
	PetscFunctionBegin;
	/* only master generates this file */
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
	if( rank != 0 ) { PetscFunctionReturn(0); }
	
	fp = fopen(pvdfilename,"w");
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
#ifdef WORDSIZE_BIGENDIAN
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	fprintf(fp,"<Collection>\n");
	
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");
	fclose(fp);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ParaviewPVDAppend"
PetscErrorCode ParaviewPVDAppend(const char pvdfilename[],double time,const char datafile[], const char DirectoryName[])
{
	PetscMPIInt rank;
	FILE *fp;
	char line[10000];
	int key_L;
	char key[] = "</Collection>";
	char *copy,*tmp;
	
	PetscFunctionBegin;
	/* only master generates this file */
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
	if( rank != 0 ) { PetscFunctionReturn(0); }
	
	
	fp = fopen(pvdfilename,"r");
	/* reset to start of file */
	rewind(fp);
	
	copy = NULL;
	key_L = strlen( key );
	while( !feof(fp) ) {
		fgets( line, 10000-1, fp );
		if ( strncmp(key,line,key_L)!=0 ) {
			
			/* copy line */
			if (copy!=NULL) {
			  asprintf(&tmp,"%s",copy);
			  free(copy);
				asprintf(&copy,"%s%s",tmp,line);
				free(tmp);
			}
			else {
			  asprintf(&copy,"%s",line);
			}
		}
		else {
			break;
		}
	}
	fclose(fp);
	
	/* open new file - clobbering the old */
	fp = fopen(pvdfilename,"w");
	
	/* write all copied chars */
	fprintf(fp,"%s",copy);
	
	/* write new data */
	fprintf(fp,"  <DataSet timestep=\"%1.6e\" file=\"./%s/%s\"/>\n",time, DirectoryName, datafile );
	
	/* close tag */
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");
	
	fclose(fp);
	free(copy);
	
	PetscFunctionReturn(0);
}


/* V-P mesh */
#undef __FUNCT__  
#define __FUNCT__ "pTatinOutputParaViewMeshVelocityPressure"
PetscErrorCode pTatinOutputParaViewMeshVelocityPressure(DM pack,Vec X,const char path[],const char prefix[])
{
	char *vtkfilename,*filename;
	PetscMPIInt rank;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = pTatinGenerateParallelVTKName(prefix,"vts",&vtkfilename);CHKERRQ(ierr);
	if (path) {
		asprintf(&filename,"%s/%s",path,vtkfilename);
	} else {
		asprintf(&filename,"./%s",vtkfilename);
	}
	
	ierr = pTatinOutputMeshVelocityPressureVTS_v0(pack,X,filename);CHKERRQ(ierr);
	free(filename);
	free(vtkfilename);
	
	ierr = pTatinGenerateVTKName(prefix,"pvts",&vtkfilename);CHKERRQ(ierr);
	if (path) {
		asprintf(&filename,"%s/%s",path,vtkfilename);
	} else {
		asprintf(&filename,"./%s",vtkfilename);
	}
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	ierr = pTatinOutputMeshVelocityPressurePVTS(pack,prefix,filename);CHKERRQ(ierr);
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
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
	}
	
	PetscPrintf(PETSC_COMM_WORLD,"[[DESIGN FLAW]] %s: only printing P0 component of pressure field \n", __FUNCT__ );

	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(dap,0,0,0,0,0,0,0,&ndof_pressure,0, 0,0,0, 0);CHKERRQ(ierr);
	
	ierr = DMDAGetGhostCorners(dau,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
	ierr = DMDAGetCornersElementQ2(dau,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
	ierr = DMDAGetCoordinateDA(dau,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(dau,&gcoords);CHKERRQ(ierr);
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
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf( vtk_fp, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	fprintf( vtk_fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", esi,esi+2*mx+1-1, esj,esj+2*my+1-1, esk,esk+2*mz+1-1);
	
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
#define __FUNCT__ "pTatinOutputMeshVelocityPressurePVTS"
PetscErrorCode pTatinOutputMeshVelocityPressurePVTS(DM pack,const char prefix[],const char name[])
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
	DMDACoor2d **LA_gcoords;	
	PetscInt mx,my,cnt;
	PetscInt ei,ej,i,j,esi,esj;
  PetscScalar *LA_cell_data;
	FILE*	vtk_fp = NULL;
	PetscInt M,N,P,swidth;
	PetscMPIInt rank;
	
	
	PetscFunctionBegin;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	vtk_fp = NULL;
	if (rank==0) {
		if ((vtk_fp = fopen ( name, "w")) == NULL)  {
			SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
		}
	}
	
	
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
	/* VTS HEADER - OPEN */	
	if(vtk_fp) fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
	if(vtk_fp) fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	
	DMDAGetInfo( dau, 0, &M,&N,&P, 0,0,0, 0,&swidth, 0,0,0, 0 );
	if(vtk_fp) fprintf( vtk_fp, "  <PStructuredGrid GhostLevel=\"%d\" WholeExtent=\"%d %d %d %d %d %d\">\n", swidth, 0,M-1, 0,N-1, 0,P-1 ); /* note overlap = 1 for Q1 */
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
#define __FUNCT__ "DAQ2PieceExtendForGhostLevelZero"
PetscErrorCode DAQ2PieceExtendForGhostLevelZero( FILE *vtk_fp, int indent_level, DM dau, const char local_file_prefix[] )
{
	PetscMPIInt nproc,rank;
	MPI_Comm comm;
	const PetscInt *lx,*ly,*lz;
	PetscInt M,N,P,pM,pN,pP,sum;
	PetscInt i,j,k,II,dim,esi,esj,esk,mx,my,mz;
	PetscInt *olx,*oly,*olz;
	PetscInt *lmx,*lmy,*lmz,*tmp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	/* create file name */
	PetscObjectGetComm( (PetscObject)dau, &comm );
	MPI_Comm_size( comm, &nproc );
	MPI_Comm_rank( comm, &rank );
	
	ierr = DMDAGetInfo( dau, &dim, &M,&N,&P, &pM,&pN,&pP, 0, 0, 0,0,0, 0 );CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRangesElementQ2(dau,&pM,&pN,&pP,&olx,&oly,&olz,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	
	if (dim==3) {
		for( k=0;k<pP;k++ ) {
			for( j=0;j<pN;j++ ) {
				for( i=0;i<pM;i++ ) {
					char *name;
					PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
					asprintf( &name, "%s-subdomain%1.5d.vts", local_file_prefix, procid );
					for( II=0; II<indent_level; II++ ) {
						if(vtk_fp) fprintf(vtk_fp,"  ");
					}
					if(vtk_fp) fprintf( vtk_fp, "<Piece Extent=\"%d %d %d %d %d %d\"      Source=\"%s\"/>\n",
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
				asprintf( &name, "%s-subdomain%1.5d.vts", local_file_prefix, procid );
				for( II=0; II<indent_level; II++ ) {
					if(vtk_fp) fprintf(vtk_fp,"  ");
				}
				if(vtk_fp) fprintf( vtk_fp, "<Piece Extent=\"%d %d %d %d 0 0\"      Source=\"%s\"/>\n",
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


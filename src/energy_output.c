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
 **    Filename:      energy_output.c
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
#include "petscvec.h"
#include "petscdm.h"
#include "ptatin3d.h"
#include "ptatin3d_defs.h"
#include "private/ptatin_impl.h"
#include "ptatin_utils.h"
#include "QPntVolCoefEnergy_def.h"

#include "dmda_view_petscvtk.h"
#include "dmdae.h"
#include "dmda_element_q1.h"
#include "quadrature.h"
#include "output_paraview.h"
#include "phys_comp_energy.h"
#include "ptatin3d_energy.h"
#include "energy_output.h"



#undef __FUNCT__
#define __FUNCT__ "_apply_threshold"
PetscErrorCode _apply_threshold(PetscScalar x[],const PetscInt N,const PetscScalar threshold,const PetscScalar value)
{
    PetscInt i;
    PetscFunctionBegin;
    for (i=0; i<N; i++) {
        PetscScalar abs;
		
        abs = PetscAbsScalar(x[i]);
        if (abs < threshold) {
            x[i] = value;
        }
    }
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputMeshEnergyVTS_ascii"
PetscErrorCode pTatinOutputMeshEnergyVTS_ascii(Quadrature Q,DM daT,Vec X,const char name[])
{
    PetscErrorCode ierr;
    PetscScalar    *LA_fields;
    DM             cda;
    Vec            local_fields,gcoords;
    DMDACoor3d     ***LA_gcoords;
    PetscInt       mx,my,mz,cnt;
    PetscInt       i,j,k,esi,esj,esk;
    FILE           *vtk_fp = NULL;
    PetscInt       gsi,gsj,gsk,gm,gn,gp;
    PetscInt       line_count;
    PetscInt       ei,ej,ek;
    PetscInt       p,nqp;
    QPntVolCoefEnergy *all_quadrature_points,*cell_quadrature_points;
	
    PetscFunctionBegin;
    if ((vtk_fp = fopen(name,"w")) == NULL)  {
        SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name);
    }
	
    ierr = DMDAGetGhostCorners(daT,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
    ierr = DMDAEGetCornersElement(daT,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
    
    ierr = DMGetCoordinateDM(daT,&cda);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(daT,&gcoords);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
    
    ierr = DMGetLocalVector(daT,&local_fields);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(daT,X,INSERT_VALUES,local_fields);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(daT,X,INSERT_VALUES,local_fields);CHKERRQ(ierr);
    ierr = VecGetArray(local_fields,&LA_fields);CHKERRQ(ierr);
    ierr = _apply_threshold(LA_fields,gm*gn*gp,1.0e-12,0.0);CHKERRQ(ierr);
	
	/* VTS HEADER - OPEN */
#ifdef WORDSIZE_BIGENDIAN
    fprintf(vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
    fprintf(vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
    PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n",esi,esi+mx, esj,esj+my, esk,esk+mz);
    PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <Piece Extent=\"%D %D %D %D %D %D\">\n",esi,esi+mx, esj,esj+my, esk,esk+mz);
	
	/* VTS COORD DATA */	
    fprintf(vtk_fp,"    <Points>\n");
    fprintf(vtk_fp,"      <DataArray Name=\"coords\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    line_count = 0;
    for (k=esk; k<esk+mz+1; k++) {
        for (j=esj; j<esj+my+1; j++) {
            for (i=esi; i<esi+mx+1; i++) {
                fprintf(vtk_fp,"%1.6e %1.6e %1.6e ", LA_gcoords[k][j][i].x, LA_gcoords[k][j][i].y, LA_gcoords[k][j][i].z );
                if (line_count%3 == 0) {
                    fprintf(vtk_fp,"\n");
                }
                line_count++;
            }
        }
    }
    fprintf(vtk_fp,"\n");
    fprintf(vtk_fp,"      </DataArray>\n");
    fprintf(vtk_fp,"    </Points>\n");
	
    /* VTS CELL DATA */
    fprintf(vtk_fp,"    <CellData>\n");
	
    nqp  = Q->npoints;
    ierr = VolumeQuadratureGetAllCellData_Energy(Q,&all_quadrature_points);CHKERRQ(ierr);
    
    /* average diffusivity and heat sources */
    fprintf(vtk_fp,"      <DataArray Name=\"diffusivity_qp_avg\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    fprintf(vtk_fp,"      ");
    line_count = 0;
    for (ek=0; ek<mz; ek++) {
        for (ej=0; ej<my; ej++) {
            for (ei=0; ei<mx; ei++) {
                double   prop,avg;
                PetscInt eidx;
				
                eidx = ei + ej*mx + ek*mx*my;
                ierr = VolumeQuadratureGetCellData_Energy(Q,all_quadrature_points,eidx,&cell_quadrature_points);CHKERRQ(ierr);
                avg = 0.0;
                for (p=0; p<nqp; p++) {
                    QPntVolCoefEnergyGetField_diffusivity(&cell_quadrature_points[p],&prop);
                    avg = avg + prop;
                }
                avg = avg / ( (double)nqp );
				
                fprintf(vtk_fp,"%1.6e ", avg );
                if (line_count%10 == 0) {
                    fprintf(vtk_fp,"\n");
                }
                line_count++;
            }
        }
    }
    fprintf(vtk_fp,"\n");
    fprintf(vtk_fp,"      </DataArray>\n");

    fprintf(vtk_fp,"      <DataArray Name=\"heatsource_qp_avg\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    fprintf(vtk_fp,"      ");
    line_count = 0;
    for (ek=0; ek<mz; ek++) {
        for (ej=0; ej<my; ej++) {
            for (ei=0; ei<mx; ei++) {
                double   prop,avg;
                PetscInt eidx;
                
                eidx = ei + ej*mx + ek*mx*my;
                ierr = VolumeQuadratureGetCellData_Energy(Q,all_quadrature_points,eidx,&cell_quadrature_points);CHKERRQ(ierr);
                avg = 0.0;
                for (p=0; p<nqp; p++) {
                    QPntVolCoefEnergyGetField_heat_source(&cell_quadrature_points[p],&prop);
                    avg = avg + prop;
                }
                avg = avg / ( (double)nqp );
                
                fprintf(vtk_fp,"%1.6e ", avg );
                if (line_count%10 == 0) {
                    fprintf(vtk_fp,"\n");
                }
                line_count++;
            }
        }
    }
    fprintf(vtk_fp,"\n");
    fprintf(vtk_fp,"      </DataArray>\n");
    
    fprintf(vtk_fp,"    </CellData>\n");
    
    /* VTS NODAL DATA */
    fprintf( vtk_fp,"    <PointData>\n");
    /* temperature */
    fprintf( vtk_fp,"      <DataArray Name=\"temperature\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    line_count = 0;
    for (k=esk; k<esk+mz+1; k++) {
        for (j=esj; j<esj+my+1; j++) {
            for (i=esi; i<esi+mx+1; i++) {
                cnt = (i-gsi) + (j-gsj)*gm + (k-gsk)*gm*gn;
                fprintf(vtk_fp,"%1.6e ", LA_fields[cnt] );
                if (line_count%10 == 0) {
                    fprintf(vtk_fp,"\n");
                }
                line_count++;
            }
        }
    }
    fprintf(vtk_fp,"\n");
    fprintf(vtk_fp,"      </DataArray>\n");
    fprintf(vtk_fp,"    </PointData>\n");
    
    /* VTS HEADER - CLOSE */
    fprintf(vtk_fp,"    </Piece>\n");
    fprintf(vtk_fp,"  </StructuredGrid>\n");
    fprintf(vtk_fp,"</VTKFile>\n");
    
    ierr = VecRestoreArray(local_fields,&LA_fields);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(daT,&local_fields);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
    
    fclose(vtk_fp);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputMeshEnergyVTS_binary"
PetscErrorCode pTatinOutputMeshEnergyVTS_binary(Quadrature Q,DM daT,Vec X,const char name[])
{
    PetscErrorCode ierr;
    PetscScalar    *LA_fields;
    DM             cda;
    Vec            local_fields,gcoords;
    DMDACoor3d     ***LA_gcoords;
    PetscInt       mx,my,mz,cnt;
    PetscInt       i,j,k,esi,esj,esk;
    FILE           *vtk_fp = NULL;
    PetscInt       gsi,gsj,gsk,gm,gn,gp;
    PetscInt       ei,ej,ek;
    PetscInt       p,nqp;
    int            offset,bytes;
    QPntVolCoefEnergy *all_quadrature_points,*cell_quadrature_points;
	
    PetscFunctionBegin;
    if ((vtk_fp = fopen(name,"w")) == NULL)  {
        SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
    }
    
    ierr = DMDAGetGhostCorners(daT,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
    ierr = DMDAEGetCornersElement(daT,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
    
    ierr = DMGetCoordinateDM(daT,&cda);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(daT,&gcoords);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
    
    ierr = DMGetLocalVector(daT,&local_fields);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(daT,X,INSERT_VALUES,local_fields);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(daT,X,INSERT_VALUES,local_fields);CHKERRQ(ierr);
    ierr = VecGetArray(local_fields,&LA_fields);CHKERRQ(ierr);
    ierr = _apply_threshold(LA_fields,gm*gn*gp,1.0e-12,0.0);CHKERRQ(ierr);
	
	/* VTS HEADER - OPEN */
#ifdef WORDSIZE_BIGENDIAN
    fprintf(vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
    fprintf(vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
    PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n",esi,esi+mx, esj,esj+my, esk,esk+mz);
    PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <Piece Extent=\"%D %D %D %D %D %D\">\n",esi,esi+mx, esj,esj+my, esk,esk+mz);

    offset = 0;
	
    /* VTS COORD DATA */
    fprintf(vtk_fp,"    <Points>\n");
	
    fprintf(vtk_fp,"      <DataArray Name=\"coords\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",offset);
    offset += sizeof(int) + sizeof(double)*3*(mx+1)*(my+1)*(mz+1);
	
    fprintf(vtk_fp,"    </Points>\n");
	
    /* VTS CELL DATA */
    fprintf(vtk_fp,"    <CellData>\n");
	
    nqp  = Q->npoints;
    ierr = VolumeQuadratureGetAllCellData_Energy(Q,&all_quadrature_points);CHKERRQ(ierr);
    
    /* average diffusivity and heat sources */
    fprintf(vtk_fp,"      <DataArray Name=\"diffusivity_qp_avg\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\" />\n",offset);
    offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
	
    fprintf(vtk_fp,"      <DataArray Name=\"heatsource_qp_avg\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\" />\n",offset);
    offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
    
    fprintf(vtk_fp,"    </CellData>\n");
	
    /* VTS NODAL DATA */
    fprintf(vtk_fp,"    <PointData>\n");
    
    /* temperature */
    fprintf(vtk_fp,"      <DataArray Name=\"temperature\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\" />\n",offset);
    offset += sizeof(int) + sizeof(double)*1*(mx+1)*(my+1)*(mz+1);
    
    fprintf(vtk_fp, "    </PointData>\n");
    
    /* VTS HEADER - CLOSE */
    fprintf(vtk_fp,"    </Piece>\n");
    fprintf(vtk_fp,"  </StructuredGrid>\n");
    fprintf(vtk_fp,"  <AppendedData encoding=\"raw\">\n");
    
    /* write tag */
    fprintf(vtk_fp,"_");
	
    /* write node coords */
    bytes = sizeof(double)*3*(mx+1)*(my+1)*(mz+1);
    fwrite(&bytes,sizeof(int),1,vtk_fp);
	
    for (k=esk; k<esk+mz+1; k++) {
        for (j=esj; j<esj+my+1; j++) {
            for (i=esi; i<esi+mx+1; i++) {
                double pos[3];
                
                pos[0] = LA_gcoords[k][j][i].x;
                pos[1] = LA_gcoords[k][j][i].y;
                pos[2] = LA_gcoords[k][j][i].z;
                fwrite(pos,sizeof(double),3,vtk_fp);
            }
        }
    }

    /* write cell diff */
    bytes = sizeof(double)*1*(mx)*(my)*(mz);
    fwrite(&bytes,sizeof(int),1,vtk_fp);
    
    for (ek=0; ek<mz; ek++) {
        for (ej=0; ej<my; ej++) {
            for (ei=0; ei<mx; ei++) {
                double   prop,avg;
                PetscInt eidx;
                
                eidx = ei + ej*mx + ek*mx*my;
                ierr = VolumeQuadratureGetCellData_Energy(Q,all_quadrature_points,eidx,&cell_quadrature_points);CHKERRQ(ierr);
                avg = 0.0;
                for (p=0; p<nqp; p++) {
                    QPntVolCoefEnergyGetField_diffusivity(&cell_quadrature_points[p],&prop);
                    avg = avg + prop;
                }
                avg = avg / ( (double)nqp );
                
                fwrite(&avg,sizeof(double),1,vtk_fp);
            }
        }
    }
    
    /* write cell heatsources */
    bytes = sizeof(double)*1*(mx)*(my)*(mz);
    fwrite(&bytes,sizeof(int),1,vtk_fp);
    
    for (ek=0; ek<mz; ek++) {
        for (ej=0; ej<my; ej++) {
            for (ei=0; ei<mx; ei++) {
                double   prop,avg;
                PetscInt eidx;
                
                eidx = ei + ej*mx + ek*mx*my;
                ierr = VolumeQuadratureGetCellData_Energy(Q,all_quadrature_points,eidx,&cell_quadrature_points);CHKERRQ(ierr);
                avg = 0.0;
                for (p=0; p<nqp; p++) {
                    QPntVolCoefEnergyGetField_heat_source(&cell_quadrature_points[p],&prop);
                    avg = avg + prop;
                }
                avg = avg / ( (double)nqp );
                
                fwrite(&avg,sizeof(double),1,vtk_fp);
            }
        }
    }
    
    /* write node temperature */
    bytes = sizeof(double)*1*(mx+1)*(my+1)*(mz+1);
    fwrite(&bytes,sizeof(int),1,vtk_fp);
    
    for (k=esk; k<esk+mz+1; k++) {
        for (j=esj; j<esj+my+1; j++) {
            for (i=esi; i<esi+mx+1; i++) {
                double val;
                
                cnt = (i-gsi) + (j-gsj)*gm + (k-gsk)*gm*gn;
                val = LA_fields[cnt];
                fwrite(&val,sizeof(double),1,vtk_fp);
            }
        }
    }
    
    fprintf(vtk_fp,"\n  </AppendedData>\n");
    fprintf(vtk_fp,"</VTKFile>\n");
    
    ierr = VecRestoreArray(local_fields,&LA_fields);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(daT,&local_fields);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
    
    fclose(vtk_fp);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputMeshEnergyVTS"
PetscErrorCode pTatinOutputMeshEnergyVTS(PetscBool binary,Quadrature Q,DM daT,Vec X,const char name[])
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	if (binary) {
		ierr = pTatinOutputMeshEnergyVTS_binary(Q,daT,X,name);CHKERRQ(ierr);
	} else {
		ierr = pTatinOutputMeshEnergyVTS_ascii(Q,daT,X,name);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAQ1PieceExtendForGhostLevelZero"
PetscErrorCode DAQ1PieceExtendForGhostLevelZero( FILE *vtk_fp, int indent_level, DM da, const char local_file_prefix[] )
{
    PetscMPIInt    nproc,rank;
    MPI_Comm       comm;
    PetscInt       M,N,P,pM,pN,pP;
    PetscInt       i,j,k,II,dim;
    PetscInt       *olx,*oly,*olz;
    PetscInt       *lmx,*lmy,*lmz;
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    /* create file name */
    PetscObjectGetComm((PetscObject)da,&comm);
    ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    
    ierr = DMDAGetInfo(da,&dim,&M,&N,&P,&pM,&pN,&pP,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAEGetOwnershipRanges(da,&pM,&pN,&pP,&olx,&oly,&olz,&lmx,&lmy,&lmz);CHKERRQ(ierr);
    
    if (dim == 3) {
        for (k=0; k<pP; k++) {
            for (j=0; j<pN; j++) {
                for (i=0; i<pM; i++) {
                    char     *name;
                    PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
                    int      procid32;
                    
                    PetscMPIIntCast(procid,&procid32);
                    asprintf(&name,"%s-subdomain%1.5d.vts",local_file_prefix,procid32);
                    for (II=0; II<indent_level; II++) {
                        if (vtk_fp) fprintf(vtk_fp,"  ");
                    }
                    if (vtk_fp) PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<Piece Extent=\"%D %D %D %D %D %D\"      Source=\"%s\"/>\n",
                                             olx[i],olx[i]+lmx[i],
                                             oly[j],oly[j]+lmy[j],
                                             olz[k],olz[k]+lmz[k],
                                             name);
                    free(name);
                }
            }
        }
    } else if (dim == 2) {
        for (j=0; j<pN; j++) {
            for (i=0; i<pM; i++) {
                char     *name;
                PetscInt procid = i + j*pM; /* convert proc(i,j,k) to pid */
                int      procid32;
                
                PetscMPIIntCast(procid,&procid32);
                asprintf(&name,"%s-subdomain%1.5d.vts",local_file_prefix,procid32);
                for (II=0; II<indent_level; II++) {
                    if (vtk_fp) fprintf(vtk_fp,"  ");
                }
                if (vtk_fp) PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<Piece Extent=\"%D %D %D %D 0 0\"      Source=\"%s\"/>\n",
                                         olx[i],olx[i]+lmx[i],
                                         oly[j],oly[j]+lmy[j],
                                         name);
                free(name);
            }
        }
    }
	//ierr = PetscFree(olx);CHKERRQ(ierr);
	//ierr = PetscFree(oly);CHKERRQ(ierr);
	//ierr = PetscFree(olz);CHKERRQ(ierr);
	
	//ierr = PetscFree(lmx);CHKERRQ(ierr);
	//ierr = PetscFree(lmy);CHKERRQ(ierr);
	//ierr = PetscFree(lmz);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputMeshEnergyPVTS"
PetscErrorCode pTatinOutputMeshEnergyPVTS(DM daT,const char prefix[],const char name[])
{
	PetscErrorCode ierr;
	FILE           *vtk_fp = NULL;
	PetscInt       M,N,P,swidth;
	PetscMPIInt    rank;
	
	PetscFunctionBegin;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	vtk_fp = NULL;
	if (rank == 0) {
		if ((vtk_fp = fopen(name,"w")) == NULL)  {
			SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file %s",name );
		}
	}
	
	/* VTS HEADER - OPEN */	
	if (vtk_fp) fprintf(vtk_fp,"<?xml version=\"1.0\"?>\n");

#ifdef WORDSIZE_BIGENDIAN
	if (vtk_fp) fprintf(vtk_fp,"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	if (vtk_fp) fprintf(vtk_fp,"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	DMDAGetInfo(daT,NULL,&M,&N,&P,NULL,NULL,NULL,NULL,&swidth,NULL,NULL,NULL,NULL);
	if (vtk_fp) PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <PStructuredGrid GhostLevel=\"%D\" WholeExtent=\"%D %D %D %D %D %D\">\n",swidth, 0,M-1, 0,N-1, 0,P-1 ); /* note overlap = 1 for Q1 */
	
	/* VTS COORD DATA */	
	if (vtk_fp) fprintf(vtk_fp,"    <PPoints>\n");
	if (vtk_fp) fprintf(vtk_fp,"      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	if (vtk_fp) fprintf(vtk_fp,"    </PPoints>\n");
		
	/* VTS CELL DATA */	
	if (vtk_fp) fprintf(vtk_fp,"    <PCellData>\n");
	if (vtk_fp) fprintf(vtk_fp,"      <PDataArray type=\"Float64\" Name=\"diffusivity_qp_avg\" NumberOfComponents=\"1\"/>\n");
	if (vtk_fp) fprintf(vtk_fp,"      <PDataArray type=\"Float64\" Name=\"heatsource_qp_avg\" NumberOfComponents=\"1\"/>\n");
	if (vtk_fp) fprintf(vtk_fp,"    </PCellData>\n");

	/* VTS NODAL DATA */
	if (vtk_fp) fprintf(vtk_fp,"    <PPointData>\n");
	if (vtk_fp) fprintf(vtk_fp,"      <PDataArray type=\"Float64\" Name=\"temperature\" NumberOfComponents=\"1\"/>\n");
	if (vtk_fp) fprintf(vtk_fp,"    </PPointData>\n");
	
	/* write out the parallel information */
	ierr = DAQ1PieceExtendForGhostLevelZero(vtk_fp,2,daT,prefix);CHKERRQ(ierr);
	
	/* VTS HEADER - CLOSE */	
	if (vtk_fp) fprintf(vtk_fp,"  </PStructuredGrid>\n");
	if (vtk_fp) fprintf(vtk_fp,"</VTKFile>\n");
	
	if (vtk_fp) fclose(vtk_fp);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinOutputParaViewMeshEnergy"
PetscErrorCode pTatinOutputParaViewMeshEnergy(Quadrature Q,DM daT,Vec X,const char path[],const char prefix[])
{
    char           *vtkfilename,*filename;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;
    ierr = pTatinGenerateParallelVTKName(prefix,"vts",&vtkfilename);CHKERRQ(ierr);
    if (path) {
        asprintf(&filename,"%s/%s",path,vtkfilename);
    } else {
        asprintf(&filename,"%s",vtkfilename);
    }
	
    ierr = pTatinOutputMeshEnergyVTS(PETSC_TRUE,Q,daT,X,filename);CHKERRQ(ierr); /* binary */
    free(filename);
    free(vtkfilename);
	
    ierr = pTatinGenerateVTKName(prefix,"pvts",&vtkfilename);CHKERRQ(ierr);
    if (path) {
        asprintf(&filename,"%s/%s",path,vtkfilename);
    } else {
        asprintf(&filename,"%s",vtkfilename);
    }
    ierr = pTatinOutputMeshEnergyPVTS(daT,prefix,filename);CHKERRQ(ierr);
    free(filename);
    free(vtkfilename);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatinOutputParaViewMeshEnergy_PerStep"
PetscErrorCode pTatinOutputParaViewMeshEnergy_PerStep(Quadrature Q,DM daT,Vec X,const char path[],const char subdomain_path[],const char prefix[])
{
    char           *vtkfilename,*filename;
    PetscErrorCode ierr;
	
    PetscFunctionBegin;

    if (!subdomain_path) {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"You must specify a subdomain directory name");
    }

    ierr = pTatinGenerateParallelVTKName(prefix,"vts",&vtkfilename);CHKERRQ(ierr);
    if (path) {
        asprintf(&filename,"%s/%s/%s",path,subdomain_path,vtkfilename);
    } else {
        asprintf(&filename,"%s/%s",subdomain_path,vtkfilename);
    }
    pTatinStringPathNormalize(filename);
	
    ierr = pTatinOutputMeshEnergyVTS(PETSC_TRUE,Q,daT,X,filename);CHKERRQ(ierr); /* binary */
    free(filename);
    free(vtkfilename);
	
    /* write pvts in subdomain directory */
    ierr = pTatinGenerateVTKName(prefix,"pvts",&vtkfilename);CHKERRQ(ierr);
    if (path) {
        asprintf(&filename,"%s/%s/%s",path,subdomain_path,vtkfilename);
    } else {
        asprintf(&filename,"%s/%s",subdomain_path,vtkfilename);
    }
    pTatinStringPathNormalize(filename);

    ierr = pTatinOutputMeshEnergyPVTS(daT,prefix,filename);CHKERRQ(ierr);
    free(filename);
    free(vtkfilename);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatin3d_ModelOutput_Temperature_Energy"
PetscErrorCode pTatin3d_ModelOutput_Temperature_Energy(pTatinCtx ctx,Vec X,const char prefix[])
{
    PetscErrorCode   ierr;
    char             *name;
    PhysCompEnergy   energy;
    DM               daT;
    PetscLogDouble   t0,t1;
    static PetscBool beenhere = PETSC_FALSE;
    static char      *pvdfilename;
    Quadrature       volQ;
    
    PetscFunctionBegin;
    PetscTime(&t0);
    
    ierr = pTatinGetContext_Energy(ctx,&energy);CHKERRQ(ierr);
    daT  = energy->daT;
    volQ = energy->volQ;
    
    /* prepare directory structures */
    ierr = pTatinParaviewSetOutputPrefix(ctx,prefix);CHKERRQ(ierr);
    
    /* Create PVD file on first entry to this function */
    if (!beenhere) {
        asprintf(&pvdfilename,"%s/timeseries_energy.pvd",ctx->outputpath);
        PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s\n",pvdfilename);
        ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
        
        beenhere = PETSC_TRUE;
    }
    
    /* Append file name to PVD file */
    if (prefix) {
        asprintf(&name,"%s_energy.pvts",prefix);
    } else {
        asprintf(&name,"energy.pvts");
    }
    switch (ctx->storage_type) {
        case TDST_FLAT:
            ierr = ParaviewPVDAppend(pvdfilename,ctx->time,name,NULL);CHKERRQ(ierr);
            break;
        case TDST_PERRANK:
            ierr = ParaviewPVDAppend(pvdfilename,ctx->time,name,NULL);CHKERRQ(ierr);
            break;
        case TDST_PERSTEP:
            ierr = ParaviewPVDAppend(pvdfilename,ctx->time,name,ctx->prefixedoutputpath);CHKERRQ(ierr);
            break;
    }
    free(name);
    
    /* Write VTS and PVTS files */
    if (prefix) {
        asprintf(&name,"%s_energy",prefix);
    } else {
        asprintf(&name,"energy");
    }
	
    switch (ctx->storage_type) {
            
        case TDST_FLAT:
            ierr = pTatinOutputParaViewMeshEnergy(volQ,daT,X,ctx->outputpath,name);CHKERRQ(ierr);
            break;
        case TDST_PERRANK:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"TDST_PERRANK not valid for pTatin3d_ModelOutput_Temperature_Energy, use storage_type = \"flat\" or \"perstep\"");
            break;
        case TDST_PERSTEP:
            break;
    }
    free(name);
    PetscTime(&t1);
    PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_energy.(pvd,pvts,vts): CPU time %1.2e (sec)\n", __FUNCT__,prefix,t1-t0);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin3d_ModelOutput_EnergyTemperature_PetscVTS"
PetscErrorCode pTatin3d_ModelOutput_EnergyTemperature_PetscVTS(pTatinCtx ctx,Vec X,const char prefix[])
{
    DM               daT;
    PhysCompEnergy   energy;
    char             *vtsfilename,*filename;
    char             date_time[1024];
    PetscLogDouble   t0,t1;
    static PetscBool been_here = PETSC_FALSE;
    static char      *pvdfilename;
    PetscErrorCode   ierr;
    
    PetscFunctionBegin;
    PetscTime(&t0);
    
    /* prepare directory structures */
    ierr = pTatinParaviewSetOutputPrefix(ctx,prefix);CHKERRQ(ierr);
    
    ierr = pTatinGetContext_Energy(ctx,&energy);CHKERRQ(ierr);
    daT  = energy->daT;
    
    /* Create PVD file on first entry to this function */
    if (!been_here) {
        if (ctx->restart_from_file) {
            pTatinGenerateFormattedTimestamp(date_time);
            asprintf(&pvdfilename,"%s/timeseries_temp_%s.pvd",ctx->outputpath,date_time);
            PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename [restarted] %s\n",pvdfilename);
        } else {
            asprintf(&pvdfilename,"%s/timeseries_temp.pvd",ctx->outputpath);
            PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s\n",pvdfilename);
        }
        ierr = ParaviewPVDOpen(pvdfilename);CHKERRQ(ierr);
        
        been_here = PETSC_TRUE;
    }
    
    /* Append file name to PVD file */
    if (prefix) {
        asprintf(&vtsfilename,"%s_temp.vts",prefix);
    } else {
        asprintf(&vtsfilename,"temp.vts");
    }
    switch (ctx->storage_type) {
        case TDST_FLAT:
            ierr = ParaviewPVDAppend(pvdfilename,ctx->time,vtsfilename,NULL);CHKERRQ(ierr);
            break;
        case TDST_PERRANK:
            break;
        case TDST_PERSTEP:
            ierr = ParaviewPVDAppend(pvdfilename,ctx->time,vtsfilename,ctx->prefixedoutputpath);CHKERRQ(ierr);
            break;
    }
	
    /* Write VTS and PVTS files */
    switch (ctx->storage_type) {
        case TDST_FLAT:
            asprintf(&filename,"%s/%s",ctx->outputpath,vtsfilename);
            break;
        case TDST_PERRANK:
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"TDST_PERRANK not valid for pTatin3d_ModelOutput_EnergyTemperature_PetscVTS, use storage_type = \"flat\" or \"perstep\"");
            break;
        case TDST_PERSTEP:
            asprintf(&filename,"%s/%s/%s",ctx->outputpath,ctx->prefixedoutputpath,vtsfilename);
            break;
    }
    pTatinStringPathNormalize(filename);
    ierr = DMDAViewPetscVTS(daT,X,filename);CHKERRQ(ierr);

    free(vtsfilename);
    free(filename);
    
    PetscTime(&t1);
    PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_vp.(pvd,vts): CPU time %1.2e (sec)\n", __FUNCT__,prefix,t1-t0);
    
    PetscFunctionReturn(0);
}

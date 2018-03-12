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
 **    filename:   output_material_points_p0.c
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

#include "ptatin3d_defs.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPStokesPl_def.h"
#include "MPntPEnergy_def.h"

#include "ptatin_utils.h"
#include "dmdae.h"
#include "dmda_checkpoint.h"
#include "dmda_duplicate.h"
#include "dmda_element_q1.h"
#include "dmda_element_q2p1.h"
#include "data_bucket.h"
#include "output_paraview.h"
#include "element_type_Q2.h"
#include "element_utils_q2.h"
#include "element_utils_q1.h"
#include "energy_output.h"

#include "output_material_points.h"
#include "output_material_points_p0.h"

#undef __FUNCT__
#define __FUNCT__ "MarkerCellFieldsP0_ProjectIntegerField"
PetscErrorCode MarkerCellFieldsP0_ProjectIntegerField(DataBucket db,MaterialPointVariable variable,Vec scalar)
{
  PetscScalar *LA_scalar,var;
  int         e,ncells,*nearestpoint_list;
	int         p,n_mp;
	MPAccess    X;
  PetscReal   *sep_list;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

  if (variable == MPV_viscosity) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_viscosity is not an integer type and cannot be projected");
  if (variable == MPV_density) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_density is not an integer and cannot be projected");
  if (variable == MPV_diffusivity) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_diffusivity is not an integer type and cannot be projected");
  if (variable == MPV_heat_source) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_heat_source is not an integer type and cannot be projected");
  if (variable == MPV_plastic_strain) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_plastic_strain is not an integer type and cannot be projected");

  ierr = VecGetLocalSize(scalar,&ncells);CHKERRQ(ierr);
  PetscMalloc(sizeof(int)*ncells,&nearestpoint_list);
  PetscMalloc(sizeof(PetscReal)*ncells,&sep_list);
  
  for (e=0; e<ncells; e++) sep_list[e] = PETSC_MAX_REAL;
  for (e=0; e<ncells; e++) nearestpoint_list[e] = -1;
  
  ierr = VecZeroEntries(scalar);CHKERRQ(ierr);
	DataBucketGetSizes(db,&n_mp,NULL,NULL);
	ierr = MaterialPointGetAccess(db,&X);CHKERRQ(ierr);
	for (p=0; p<n_mp; p++) {
    int    eidx;
    double *xip,sep2;
    
		ierr = MaterialPointGet_local_element_index(X,p,&eidx);CHKERRQ(ierr);
		ierr = MaterialPointGet_local_coord(X,p,&xip);CHKERRQ(ierr);
		
    /* compute closest point to origin in local coordinate space (0,0,0) */
    sep2 = xip[0]*xip[0] + xip[1]*xip[1] + xip[2]*xip[2];

    if (sep2 < sep_list[eidx]) {
      sep_list[eidx]          = (PetscScalar)sep2;
      nearestpoint_list[eidx] = p;
    }
	}

  ierr = VecGetArray(scalar,&LA_scalar);CHKERRQ(ierr);
  for (e=0; e<ncells; e++) {
    
    p = nearestpoint_list[e];
    
    var = 0.0;
		switch (variable) {
        
      case MPV_yield_indicator: {
        short yindicator;
        
        ierr = MaterialPointGet_yield_indicator(X,p,&yindicator);CHKERRQ(ierr);
        var = (PetscScalar)yindicator;
      }
				break;
        
			case MPV_region: {
        int region_id;
        
				ierr = MaterialPointGet_phase_index(X,p,&region_id);CHKERRQ(ierr);
        var = (PetscScalar)region_id;
      }
				break;
        
      default:
        break;
    }
    LA_scalar[e] = var;
  }
  ierr = VecRestoreArray(scalar,&LA_scalar);CHKERRQ(ierr);
  ierr = MaterialPointRestoreAccess(db,&X);CHKERRQ(ierr);
  PetscFree(sep_list);
  PetscFree(nearestpoint_list);
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MarkerCellFieldsP0_ProjectScalarField"
PetscErrorCode MarkerCellFieldsP0_ProjectScalarField(DataBucket db,MaterialPointVariable variable,Vec pointcounts,Vec scalar)
{
  PetscScalar *LA_scalar,var;
	int         p,n_mp;
	MPAccess    X;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
  if (variable == MPV_yield_indicator) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_yield_indicator is not a floating point type and cannot be projected");
  if (variable == MPV_region) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MPV_region is not a floating point type and cannot be projected");
  
  ierr = VecZeroEntries(scalar);CHKERRQ(ierr);
  ierr = VecGetArray(scalar,&LA_scalar);CHKERRQ(ierr);
	DataBucketGetSizes(db,&n_mp,NULL,NULL);
	ierr = MaterialPointGetAccess(db,&X);CHKERRQ(ierr);
	for (p=0; p<n_mp; p++) {
    int eidx;
		ierr = MaterialPointGet_local_element_index(X,p,&eidx);CHKERRQ(ierr);
		
    var = 0.0;
		switch (variable) {
			
      case MPV_viscosity:
        ierr = MaterialPointGet_viscosity(X,p,&var);CHKERRQ(ierr);
				break;
        
			case MPV_density:
				ierr = MaterialPointGet_density(X,p,&var);CHKERRQ(ierr);
				break;
			
      case MPV_diffusivity:
				ierr = MaterialPointGet_diffusivity(X,p,&var);CHKERRQ(ierr);
				break;
			
      case MPV_heat_source:
				ierr = MaterialPointGet_heat_source(X,p,&var);CHKERRQ(ierr);
				break;
        
			case MPV_plastic_strain: {
        float _var;
				ierr = MaterialPointGet_plastic_strain(X,p,&_var);CHKERRQ(ierr);
        var = (PetscScalar)_var;
      }
				break;

      default:
        break;
		}
		
		LA_scalar[eidx] += var;
	}
	ierr = MaterialPointRestoreAccess(db,&X);CHKERRQ(ierr);
  ierr = VecRestoreArray(scalar,&LA_scalar);CHKERRQ(ierr);
	
  ierr = VecPointwiseDivide(scalar,scalar,pointcounts);CHKERRQ(ierr);
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MarkerCellFieldsP0_CountPointsPerCell"
PetscErrorCode MarkerCellFieldsP0_CountPointsPerCell(DM dmp0,DataBucket db,Vec pointcounts)
{
	int            eidx,p,n_mp;
	MPAccess       X;
  PetscScalar    *LA_pointcounts;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
  
	ierr = VecZeroEntries(pointcounts);CHKERRQ(ierr);
  ierr = VecGetArray(pointcounts,&LA_pointcounts);CHKERRQ(ierr);
	DataBucketGetSizes(db,&n_mp,NULL,NULL);
	
	ierr = MaterialPointGetAccess(db,&X);CHKERRQ(ierr);
	for (p=0; p<n_mp; p++) {
		ierr = MaterialPointGet_local_element_index(X,p,&eidx);CHKERRQ(ierr);
		LA_pointcounts[eidx] += 1.0;
	}
	ierr = MaterialPointRestoreAccess(db,&X);CHKERRQ(ierr);
  ierr = VecRestoreArray(pointcounts,&LA_pointcounts);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MarkerCellFieldsP0Read_PetscVec"
PetscErrorCode MarkerCellFieldsP0Read_PetscVec(Vec scalar,const MaterialPointVariable var_name,
                                               const char basename[])
{
  PetscViewer    viewer;
  PetscErrorCode ierr;
  char           fname[PETSC_MAX_PATH_LEN];
  
  PetscFunctionBegin;
  {
    /* load cell data */
    PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"%s.%s.vec",basename,MaterialPointVariableName[var_name]);
    PetscPrintf(PETSC_COMM_WORLD,"  Loading MarkerCellField file: %s\n",fname);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fname,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(scalar,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MarkerCellFieldsP0Write_ParaViewVTS"
PetscErrorCode MarkerCellFieldsP0Write_ParaViewVTS(DM dmscalar,DM dmp0,Vec scalar,Vec pointcounts,
                                                      DataBucket material_points,
                                                      const char basename[],
                                                      const int nvars,const MaterialPointVariable vars[],
                                                      PetscBool low_precision,
                                                      const char vtkfilename[])
{
	PetscErrorCode ierr;
  PetscScalar    *LA_scalar;
	DM             cda;
	Vec            gcoords;
	DMDACoor3d     ***LA_gcoords;
	PetscInt       mx,my,mz;
	PetscInt       i,j,k,esi,esj,esk;
	FILE           *vtk_fp = NULL;
	PetscInt       gsi,gsj,gsk,gm,gn,gp;
	PetscInt       t,e;
	int            offset,bytes;
	
	PetscFunctionBegin;
  
  if (material_points && basename) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify one of either material points or basename of petsc vec file to load");
  }
  if (!material_points && !basename) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify one of either material points or basename of petsc vec file to load");
  }
  
	if ((vtk_fp = fopen (vtkfilename,"w")) == NULL)  {
		SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",vtkfilename);
	}
	
	ierr = DMDAGetGhostCorners(dmscalar,&gsi,&gsj,&gsk,&gm,&gn,&gp);CHKERRQ(ierr);
	ierr = DMDAEGetCornersElement(dmscalar,&esi,&esj,&esk,&mx,&my,&mz);CHKERRQ(ierr);
	
	ierr = DMGetCoordinateDM(dmscalar,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dmscalar,&gcoords);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	/* VTS HEADER - OPEN */
#ifdef WORDSIZE_BIGENDIAN
	fprintf(vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(vtk_fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n",esi,esi+mx,esj,esj+my,esk,esk+mz);
	PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <Piece Extent=\"%D %D %D %D %D %D\">\n",esi,esi+mx,esj,esj+my,esk,esk+mz);
  
	offset = 0;
	
	/* VTS COORD DATA */
	fprintf(vtk_fp,"    <Points>\n");
	
  if (!low_precision) {
    fprintf(vtk_fp,"      <DataArray Name=\"coords\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",offset);
    offset += sizeof(int) + sizeof(double)*3*(mx+1)*(my+1)*(mz+1);
  } else {
    fprintf(vtk_fp,"      <DataArray Name=\"coords\" type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",offset);
    offset += sizeof(int) + sizeof(float)*3*(mx+1)*(my+1)*(mz+1);
  }
  
	fprintf(vtk_fp,"    </Points>\n");
	
	/* VTS CELL DATA */
	fprintf(vtk_fp,"    <CellData>\n");
	
  for (t=0; t<nvars; t++) {
    MaterialPointVariable idx = vars[t];
    
    if (!low_precision) {
      fprintf(vtk_fp,"      <DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\" />\n",MaterialPointVariableName[idx],offset);
    } else {
      fprintf(vtk_fp,"      <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\" />\n",MaterialPointVariableName[idx],offset);
    }
    
    switch (idx) {
        
      case MPV_region:
        if (!low_precision) offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
        else                offset += sizeof(int) + sizeof(float)*1*(mx)*(my)*(mz);
        break;
        
      case MPV_viscosity:
        if (!low_precision) offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
        else                offset += sizeof(int) + sizeof(float)*1*(mx)*(my)*(mz);
        break;
        
      case MPV_density:
        if (!low_precision) offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
        else                offset += sizeof(int) + sizeof(float)*1*(mx)*(my)*(mz);
        break;
        
      case MPV_plastic_strain:
        if (!low_precision) offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
        else                offset += sizeof(int) + sizeof(float)*1*(mx)*(my)*(mz);
        break;
        
      case MPV_yield_indicator:
        if (!low_precision) offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
        else                offset += sizeof(int) + sizeof(float)*1*(mx)*(my)*(mz);
        break;
        
      case MPV_diffusivity:
        if (!low_precision) offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
        else                offset += sizeof(int) + sizeof(float)*1*(mx)*(my)*(mz);
        break;
        
      case MPV_heat_source:
        if (!low_precision) offset += sizeof(int) + sizeof(double)*1*(mx)*(my)*(mz);
        else                offset += sizeof(int) + sizeof(float)*1*(mx)*(my)*(mz);
        break;
        
      default:
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unknown material point variable type (MPV_type)");
        break;
    }
  }
  
	fprintf(vtk_fp,"    </CellData>\n");
	
	/* VTS NODAL DATA */
	fprintf(vtk_fp,"    <PointData>\n");
	fprintf(vtk_fp,"    </PointData>\n");
	
	/* VTS HEADER - CLOSE */
	fprintf(vtk_fp,"    </Piece>\n");
	fprintf(vtk_fp,"  </StructuredGrid>\n");
	fprintf(vtk_fp,"  <AppendedData encoding=\"raw\">\n");
	
	/* write tag */
	fprintf(vtk_fp,"_");
	
	/* write node coords */
  if (!low_precision) {
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
  } else {
    bytes = sizeof(float)*3*(mx+1)*(my+1)*(mz+1);
    fwrite(&bytes,sizeof(int),1,vtk_fp);

    for (k=esk; k<esk+mz+1; k++) {
      for (j=esj; j<esj+my+1; j++) {
        for (i=esi; i<esi+mx+1; i++) {
          float fpos[3];
       
          fpos[0] = (float)LA_gcoords[k][j][i].x;
          fpos[1] = (float)LA_gcoords[k][j][i].y;
          fpos[2] = (float)LA_gcoords[k][j][i].z;
          fwrite(fpos,sizeof(float),3,vtk_fp);
        }
      }
    }
  }
  
  /* write cell data */
  for (t=0; t<nvars; t++) {
    MaterialPointVariable idx = vars[t];
    
    if (material_points) {
      /* compute average from material points */
      if ( (idx == MPV_region) || (idx == MPV_yield_indicator) ){
        /* perform closest point projection */
        ierr = MarkerCellFieldsP0_ProjectIntegerField(material_points,idx,scalar);CHKERRQ(ierr);
      } else {
        /* perform P0 projection */
        ierr = MarkerCellFieldsP0_ProjectScalarField(material_points,idx,pointcounts,scalar);CHKERRQ(ierr);
      }
    } else {
      /* load average from petscvec */
      ierr = MarkerCellFieldsP0Read_PetscVec(scalar,idx,basename);CHKERRQ(ierr);
    }
    
    ierr = VecGetArray(scalar,&LA_scalar);CHKERRQ(ierr);

    if (!low_precision) {
      bytes = sizeof(double)*1*(mx)*(my)*(mz);
      fwrite(&bytes,sizeof(int),1,vtk_fp);

      for (e=0; e<mx*my*mz; e++) {
        double value;
        
        value = LA_scalar[e];
        fwrite(&value,sizeof(double),1,vtk_fp);
      }
    } else {
      bytes = sizeof(float)*1*(mx)*(my)*(mz);
      fwrite(&bytes,sizeof(int),1,vtk_fp);

      for (e=0; e<mx*my*mz; e++) {
        float value;
        
        value = (float)LA_scalar[e];
        fwrite(&value,sizeof(float),1,vtk_fp);
      }
    }
    ierr = VecRestoreArray(scalar,&LA_scalar);CHKERRQ(ierr);
    
  }

	fprintf(vtk_fp,"\n  </AppendedData>\n");
	fprintf(vtk_fp,"</VTKFile>\n");
	
	ierr = DMDAVecRestoreArray(cda,gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	fclose(vtk_fp);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MarkerCellFieldsP0Write_ParaViewPVTS"
PetscErrorCode MarkerCellFieldsP0Write_ParaViewPVTS(DM dascalar,const int nvars,const MaterialPointVariable vars[],PetscBool low_precision,const char prefix[],const char name[])
{
	PetscErrorCode ierr;
	FILE           *vtk_fp = NULL;
	PetscInt       M,N,P,swidth;
	PetscMPIInt    rank;
  int            t;
	
	PetscFunctionBegin;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	vtk_fp = NULL;
	if (rank==0) {
		if ((vtk_fp = fopen (name,"w")) == NULL)  {
			SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",name);
		}
	}
	
	/* VTS HEADER - OPEN */
	if (vtk_fp) fprintf(vtk_fp,"<?xml version=\"1.0\"?>\n");
  
#ifdef WORDSIZE_BIGENDIAN
	if (vtk_fp) fprintf(vtk_fp,"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	if (vtk_fp) fprintf(vtk_fp,"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	DMDAGetInfo(dascalar,0,&M,&N,&P,0,0,0,0,&swidth,0,0,0,0);
	if (vtk_fp) PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <PStructuredGrid GhostLevel=\"%D\" WholeExtent=\"%D %D %D %D %D %D\">\n",swidth,0,M-1,0,N-1,0,P-1); /* note overlap = 1 for Q1 */
	
	/* VTS COORD DATA */
	if (vtk_fp) fprintf(vtk_fp,"    <PPoints>\n");
  if (!low_precision) {
    if (vtk_fp) fprintf(vtk_fp,"      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
  } else {
    if (vtk_fp) fprintf(vtk_fp,"      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
  }
	if (vtk_fp) fprintf(vtk_fp,"    </PPoints>\n");
  
	/* VTS CELL DATA */
	if (vtk_fp) fprintf(vtk_fp,"    <PCellData>\n");
  
  for (t=0; t<nvars; t++) {
    MaterialPointVariable idx = vars[t];
    
    if (!low_precision) {
      if (vtk_fp) fprintf(vtk_fp,"      <PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",MaterialPointVariableName[idx]);
    } else {
      if (vtk_fp) fprintf(vtk_fp,"      <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",MaterialPointVariableName[idx]);
    }
  }
	if (vtk_fp) fprintf(vtk_fp,"    </PCellData>\n");
  
	/* VTS NODAL DATA */
	if (vtk_fp) fprintf(vtk_fp,"    <PPointData>\n");
	if (vtk_fp) fprintf(vtk_fp,"    </PPointData>\n");
	
	/* write out the parallel information */
	ierr = DAQ1PieceExtendForGhostLevelZero(vtk_fp,2,dascalar,prefix);CHKERRQ(ierr);
	
	/* VTS HEADER - CLOSE */
	if (vtk_fp) fprintf( vtk_fp,"  </PStructuredGrid>\n");
	if (vtk_fp) fprintf( vtk_fp,"</VTKFile>\n");
	
	if (vtk_fp) fclose(vtk_fp);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MarkerCellFieldsP0Write_ParaView"
PetscErrorCode MarkerCellFieldsP0Write_ParaView(DM pack,DataBucket material_points,const char basename[],
                                                  const int nvars,const MaterialPointVariable vars[],
                                                  PetscBool low_precision,
                                                  const char path[],const char prefix[])
{
  MPI_Comm    comm;
	char        *vtkfilename,filename[PETSC_MAX_PATH_LEN];
	PetscMPIInt rank;
	DM          dau,dap;
  DM          dmscalar,dmp0;
	PetscInt    MX,MY,MZ,Mp,Np,Pp,*lxv,*lyv,*lzv;
  Vec         cellconstant,pointcounts = NULL;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
  
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  
	ierr = pTatinGenerateParallelVTKName(prefix,"vts",&vtkfilename);CHKERRQ(ierr);
	if (path) { PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/%s",path,vtkfilename);
	} else {    PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"./%s",vtkfilename); }
  
	ierr = DMCompositeGetEntries(pack,&dau,&dap);CHKERRQ(ierr);
	
  /* create scalar Q1 dmda */
  ierr = DMDACreateOverlappingQ1FromQ2(dau,1,&dmscalar);CHKERRQ(ierr);
  
  /* create scalar P0 dmda which overlaps the Q1 dmda */
  ierr = PetscObjectGetComm((PetscObject)dau,&comm);CHKERRQ(ierr);
  ierr = DMDAGetSizeElementQ2(dau,&MX,&MY,&MZ);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dau,0,0,0,0,&Mp,&Np,&Pp,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRangesElementQ2(dau,0,0,0,0,0,0,&lxv,&lyv,&lzv);CHKERRQ(ierr);
	ierr = DMDACreate3d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,MX,MY,MZ,Mp,Np,Pp,1,0,lxv,lyv,lzv,&dmp0);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(dmp0,&cellconstant);CHKERRQ(ierr);
  
  /* count points per cell */
  if (material_points) {
    ierr = DMCreateGlobalVector(dmp0,&pointcounts);CHKERRQ(ierr);
    ierr = MarkerCellFieldsP0_CountPointsPerCell(dmp0,material_points,pointcounts);CHKERRQ(ierr);
  }
	ierr = MarkerCellFieldsP0Write_ParaViewVTS(dmscalar,dmp0,cellconstant,pointcounts,material_points,basename,nvars,vars,low_precision,filename);CHKERRQ(ierr);
  
	ierr = pTatinGenerateVTKName(prefix,"pvts",&vtkfilename);CHKERRQ(ierr);
	if (path) { PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s/%s",path,vtkfilename);
	} else {    PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"./%s",vtkfilename); }
  
	ierr = MarkerCellFieldsP0Write_ParaViewPVTS(dmscalar,nvars,vars,low_precision,prefix,filename);CHKERRQ(ierr);
	
  ierr = VecDestroy(&cellconstant);CHKERRQ(ierr);
  if (material_points) {
    ierr = VecDestroy(&pointcounts);CHKERRQ(ierr);
  }
  ierr = DMDestroy(&dmp0);CHKERRQ(ierr);
  ierr = DMDestroyDMDAE(dmscalar);CHKERRQ(ierr);
  ierr = DMDestroy(&dmscalar);CHKERRQ(ierr);
	free(vtkfilename);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin3dModelOutput_MarkerCellFieldsP0_ParaView"
PetscErrorCode pTatin3dModelOutput_MarkerCellFieldsP0_ParaView(pTatinCtx ctx,const int nvars,const MaterialPointVariable vars[],PetscBool low_precision,const char prefix[])
{
	PetscErrorCode ierr;
	char           name[PETSC_MAX_PATH_LEN];
	DM             stokes_pack;
	PetscLogDouble t0,t1;
	static PetscBool beenhere=PETSC_FALSE;
	DataBucket     material_points;

	PetscFunctionBegin;
	PetscTime(&t0);
  if (nvars == 0) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify at least one marker field to project");
  
	// PVD
	{
    char pvdfilename[PETSC_MAX_PATH_LEN];
		char vtkfilename[PETSC_MAX_PATH_LEN];
		
    PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_mpoints_cell.pvd",ctx->outputpath);
		if (prefix) {
			PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"%s_mpoints_cell.pvts",prefix);
		} else {
			PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"mpoints_cell.pvts");
		}
		
    if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n",pvdfilename); }
		ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time,vtkfilename,"");CHKERRQ(ierr);
    beenhere = PETSC_TRUE;
	}
	
	// PVTS + VTS
	if (prefix) {
		PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s_mpoints_cell",prefix);
	} else {
		PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"mpoints_cell");
	}
	
	ierr = pTatinGetMaterialPoints(ctx,&material_points,NULL);CHKERRQ(ierr);
	stokes_pack = ctx->stokes_ctx->stokes_pack;
	ierr = MarkerCellFieldsP0Write_ParaView(stokes_pack,material_points,NULL,nvars,vars,low_precision,ctx->outputpath,name);CHKERRQ(ierr);
	
	PetscTime(&t1);
	/*PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_mpoints_cell.(pvd,pvts,vts): CPU time %1.2e (sec) \n",__FUNCT__,prefix,t1-t0);*/
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MarkerCellFieldsP0Write_PetscVec"
PetscErrorCode MarkerCellFieldsP0Write_PetscVec(DM dmscalar,DM dmp0,Vec scalar,Vec pointcounts,
                                                DataBucket material_points,const int nvars,const MaterialPointVariable vars[],
                                                const char basename[])
{
  PetscViewer    viewer;
	PetscErrorCode ierr;
  int            t;
  char           fname[PETSC_MAX_PATH_LEN];
	
	PetscFunctionBegin;
  for (t=0; t<nvars; t++) {
    MaterialPointVariable idx = vars[t];
    
    if ( (idx == MPV_region) || (idx == MPV_yield_indicator) ) {
      /* perform closest point projection */
      ierr = MarkerCellFieldsP0_ProjectIntegerField(material_points,idx,scalar);CHKERRQ(ierr);
    } else {
      /* perform P0 projection */
      ierr = MarkerCellFieldsP0_ProjectScalarField(material_points,idx,pointcounts,scalar);CHKERRQ(ierr);
    }

    /* write cell data */
    PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"%s.%s.vec",basename,MaterialPointVariableName[idx]);
    PetscPrintf(PETSC_COMM_WORLD,"  Writing MarkerCellField file: %s\n",fname);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fname,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
		ierr = VecView(scalar,viewer);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "pTatin3dModelOutput_MarkerCellFieldsP0_PetscVec"
PetscErrorCode pTatin3dModelOutput_MarkerCellFieldsP0_PetscVec(pTatinCtx ctx,PetscBool dm_velocity_data_required,const int nvars,const MaterialPointVariable vars[],const char prefix[])
{
	PetscErrorCode ierr;
	char           basename[PETSC_MAX_PATH_LEN];
	DM             stokes_pack;
	PetscLogDouble t0,t1;
	static PetscBool beenhere=PETSC_FALSE;
	DataBucket     material_points;
	DM             dau,dap,dmscalar,dmp0;
	PetscInt       MX,MY,MZ,Mp,Np,Pp,*lxv,*lyv,*lzv;
  Vec            cellconstant,pointcounts;
  MPI_Comm       comm;

	PetscFunctionBegin;
  if (nvars == 0) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify at least one marker field to project");

	PetscTime(&t0);
	// PVD
	{
    char pvdfilename[PETSC_MAX_PATH_LEN];
		char vtkfilename[PETSC_MAX_PATH_LEN];
		
    PetscSNPrintf(pvdfilename,PETSC_MAX_PATH_LEN-1,"%s/timeseries_mpoints_cell.pvd",ctx->outputpath);
		if (prefix) {
			PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"%s_mpoints_cell.pvts",prefix);
		} else {
			PetscSNPrintf(vtkfilename,PETSC_MAX_PATH_LEN-1,"mpoints_cell.pvts");
		}
		
    if (!beenhere) { PetscPrintf(PETSC_COMM_WORLD,"  writing pvdfilename %s \n",pvdfilename); }
		ierr = ParaviewPVDOpenAppend(beenhere,ctx->step,pvdfilename,ctx->time,vtkfilename,"");CHKERRQ(ierr);
    beenhere = PETSC_TRUE;
	}

	stokes_pack = ctx->stokes_ctx->stokes_pack;
	ierr = DMCompositeGetEntries(stokes_pack,&dau,&dap);CHKERRQ(ierr);

  if (dm_velocity_data_required) {
    char f1[PETSC_MAX_PATH_LEN];
    char f2[PETSC_MAX_PATH_LEN];

    if (prefix) {
      PetscSNPrintf(f1,PETSC_MAX_PATH_LEN-1,"%s/%s.dmda-velocity",ctx->outputpath,prefix);
      PetscSNPrintf(f2,PETSC_MAX_PATH_LEN-1,"%s/%s.dmda-pressure",ctx->outputpath,prefix);
    } else {
      PetscSNPrintf(f1,PETSC_MAX_PATH_LEN-1,"%s/dmda-velocity",ctx->outputpath);
      PetscSNPrintf(f2,PETSC_MAX_PATH_LEN-1,"%s/dmda-pressure",ctx->outputpath);
    }
    ierr = DMDACheckpointWrite(ctx->stokes_ctx->dav,f1);CHKERRQ(ierr);
    ierr = DMDACheckpointWrite(ctx->stokes_ctx->dap,f2);CHKERRQ(ierr);
  }
  
  /* setup dm's */
  /* create scalar Q1 dmda */
  ierr = DMDACreateOverlappingQ1FromQ2(dau,1,&dmscalar);CHKERRQ(ierr);
  
  /* create scalar P0 dmda which overlaps the Q1 dmda */
  ierr = PetscObjectGetComm((PetscObject)dau,&comm);CHKERRQ(ierr);
  ierr = DMDAGetSizeElementQ2(dau,&MX,&MY,&MZ);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dau,0,0,0,0,&Mp,&Np,&Pp,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRangesElementQ2(dau,0,0,0,0,0,0,&lxv,&lyv,&lzv);CHKERRQ(ierr);
	ierr = DMDACreate3d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,MX,MY,MZ,Mp,Np,Pp,1,0,lxv,lyv,lzv,&dmp0);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(dmp0,&cellconstant);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dmp0,&pointcounts);CHKERRQ(ierr);
	
  ierr = pTatinGetMaterialPoints(ctx,&material_points,NULL);CHKERRQ(ierr);
  ierr = MarkerCellFieldsP0_CountPointsPerCell(dmp0,material_points,pointcounts);CHKERRQ(ierr);

  /* barf out dmscalar, dmp0, and the vectors */
  if (prefix) {
		PetscSNPrintf(basename,PETSC_MAX_PATH_LEN-1,"%s/%s.dmda-cell",ctx->outputpath,prefix);
	} else {
		PetscSNPrintf(basename,PETSC_MAX_PATH_LEN-1,"%s/dmda-cell",ctx->outputpath);
	}
  
  ierr = MarkerCellFieldsP0Write_PetscVec(dmscalar,dmp0,cellconstant,pointcounts,
                                           material_points,nvars,vars,
                                           basename);CHKERRQ(ierr);

  {
    PetscViewer viewer;
    char coorname[PETSC_MAX_PATH_LEN];
    Vec coor;
    
    PetscSNPrintf(coorname,PETSC_MAX_PATH_LEN-1,"%s.coords.vec",basename);
    ierr = DMGetCoordinates(dmscalar,&coor);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,coorname,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(coor,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
  PetscTime(&t1);
	/*PetscPrintf(PETSC_COMM_WORLD,"%s() -> %s_mpoints_cell.(vec): CPU time %1.2e (sec) \n",__FUNCT__,prefix,t1-t0);*/

  ierr = VecDestroy(&cellconstant);CHKERRQ(ierr);
  ierr = VecDestroy(&pointcounts);CHKERRQ(ierr);
  ierr = DMDestroy(&dmp0);CHKERRQ(ierr);
  ierr = DMDestroyDMDAE(dmscalar);CHKERRQ(ierr);
  ierr = DMDestroy(&dmscalar);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


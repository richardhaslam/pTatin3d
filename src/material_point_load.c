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
 **    Filename:      material_point_load.c
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <petsc.h>

#include "dmda_element_q2p1.h"
#include "data_bucket.h"
#include "dmda_checkpoint.h"

#include "MPntStd_def.h"
#include "MPntPStokes_def.h"
#include "MPntPStokesPl_def.h"
#include "MPntPEnergy_def.h"

#include "QPntVolCoefStokes_def.h"
#include "QPntSurfCoefStokes_def.h"
#include "QPntVolCoefEnergy_def.h"

#include "material_point_std_utils.h"
#include "material_point_point_location.h"
#include "material_point_load.h"

#undef __FUNCT__
#define __FUNCT__ "MarkerCoordinatesLoadFromFile"
PetscErrorCode MarkerCoordinatesLoadFromFile(const char name[],long int *length,double **coords)
{
	FILE *fp = NULL;
	long int n_markers;
	double *data;
	int p;
	char line[256];
	
	fp = fopen(name,"rb");
	if (fp == NULL) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"File %s not found",name);
	}
	
	fgets(line,255,fp);
	//vtk_data_type = atoi( line );
	fgets(line,255,fp);
	n_markers = atol( line );

	*length = n_markers;
	data = malloc( sizeof(double)*3*n_markers );
	memset(data,0,sizeof(double)*3*n_markers);

	for (p=0; p<n_markers; p++) {
		fread(&data[3*p],sizeof(double),3,fp);
		//printf("%1.4e %1.4e %1.4e \n", data[3*p],data[3*p+1],data[3*p+2]);
	}	
	
	fclose(fp);
	
	*coords = data;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MarkerScalarFieldLoadFromFile"
PetscErrorCode MarkerScalarFieldLoadFromFile(const char name[],long int *length,void **field)
{
	FILE *fp = NULL;
	int vtk_data_type;
	long int n_markers;
	double *data;
	char line[256];
	size_t datasize;
	
	fp = fopen(name,"rb");
	if (fp == NULL) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"File %s not found",name);
	}
	
	fgets(line,255,fp);
	vtk_data_type = atoi( line );
	fgets(line,255,fp);
	n_markers = atol( line );

	// write field types
	datasize = 0;
	switch (vtk_data_type) {
		case 0:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_VOID");
			break;
		case 1:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_BIT");
			break;
		case 2:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_CHAR");
			datasize = sizeof(unsigned char);
			break;
		case 3:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_SIGNED_CHAR");
			datasize = sizeof(char);
			break;
		case 4:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_SHORT");
			datasize = sizeof(short);
			break;
		case 5:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_UNSIGNED_SHORT");
			datasize = sizeof(unsigned short);
			break;
		case 6:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_INT");
			datasize = sizeof(int);
			break;
		case 7:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_UNSIGNED_INT");
			datasize = sizeof(unsigned int);
			break;
		case 8:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_LONG");
			datasize = sizeof(long);
			break;
		case 9:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_UNSIGNED_LONG");
			datasize = sizeof(unsigned long);
			break;
		case 10:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_FLOAT");
			datasize = sizeof(float);
			break;
		case 11:
			//SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_DOUBLE");
			datasize = sizeof(double);
			break;
		case 12:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Cannot read vtk data type VTK_ID_TYPE");
			break;
		default:
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"No assumed default vtk data type");
			break;
	}
	
	
	*length = n_markers;
	data = malloc( datasize*n_markers );
	memset(data,0,datasize*n_markers);
	
	fread(data,datasize,n_markers,fp);
	
	fclose(fp);
	
	*field = (void*)data;
	
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStdRemoval"
PetscErrorCode MaterialPointStdRemoval(DataBucket db,long int start,long int npoints,const int wil_key)
{
	DataField         PField_std;
	long int p,end;
	
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	printf("  loaded n_mp_points %ld \n", npoints );
	end = start+npoints;
	for (p=start; p<end; p++) {
		MPntStd     *material_point;
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		if (material_point->wil==wil_key) {
			DataBucketRemovePointAtIndex(db,p);
			end--;
			p--;
		}
	}	
	DataFieldRestoreAccess(PField_std);
	
	printf("  kept   n_mp_points %ld \n", end-start );
	DataBucketSetSizes(db,end,-1); // -ve buffer val retains old value //
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointStdInsertBasic"
PetscErrorCode MaterialPointStdInsertBasic(DataBucket db,DM da,long int start,long int npoints,double coords_mp[],int phase_mp[])
{
	DM cda;
	Vec gcoords;
	PetscScalar *LA_gcoords;
	double tolerance;
	int p,max_its;
	Truth use_nonzero_guess, monitor, log;
	PetscInt lmx,lmy,lmz;
	PetscInt nel,nen_u;
	const PetscInt *elnidx_u;
	DataField PField_std;
	PetscErrorCode ierr;
	
	/* setup for coords */
	ierr = DMGetCoordinateDM(da,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	
	ierr = DMDAGetElements_pTatinQ2P1(da,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);
	
	ierr = DMDAGetLocalSizeElementQ2(da,&lmx,&lmy,&lmz);CHKERRQ(ierr);
	
	/* point location parameters */
	tolerance         = 1.0e-10;
	max_its           = 10;
	use_nonzero_guess = 0; //_FALSE;
	monitor           = 0; //_FALSE;
	log               = 0; //_FALSE;
	
	
	DataBucketGetDataFieldByName(db,MPntStd_classname,&PField_std);
	DataFieldGetAccess(PField_std);
	DataFieldVerifyAccess(PField_std,sizeof(MPntStd));
	
	printf("inserting %ld points from index %ld \n", npoints,start);
	for (p=start; p<start+npoints; p++) {
		MPntStd     *material_point;
		MPntStd     test_p;
		int         pindex = p - start;
		
		test_p.pid     = 0;
		test_p.coor[0] = coords_mp[3*pindex  ];
		test_p.coor[1] = coords_mp[3*pindex+1];
		test_p.coor[2] = coords_mp[3*pindex+2];
		test_p.phase   = phase_mp[pindex];
		//printf("test_p.phase = %d: %d \n", test_p.phase, phase_mp[pindex] );
		test_p.wil     = -1;
		
		DataFieldInsertPoint(PField_std,p,(void*)&test_p);
		
		DataFieldAccessPoint(PField_std,p,   (void**)&material_point);
		
		InverseMappingDomain_3dQ2(tolerance, max_its,
															use_nonzero_guess, 
															monitor, log,
															(const double*)LA_gcoords, (const int)lmx,(const int)lmy,(const int)lmz, (const int*)elnidx_u,
															1, material_point );
	}
	
	ierr = VecRestoreArray(gcoords,&LA_gcoords);CHKERRQ(ierr);
	DataFieldRestoreAccess(PField_std);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MaterialPointDataBasicLoadIntoListFromFile"
PetscErrorCode MaterialPointDataBasicLoadIntoListFromFile(DataBucket db,DM da,PetscBool append,const char coordfile[],const char phasefile[])
{
	long int N1,N2;
	double *coords_mp;
	int *phase_mp;
	long int start;
	int n_mp_points;
	PetscErrorCode ierr;
	
	/* read in from file */
	printf("reading files %s : %s \n", coordfile,phasefile);
	ierr = MarkerCoordinatesLoadFromFile(coordfile,&N1,&coords_mp);CHKERRQ(ierr);
	ierr = MarkerScalarFieldLoadFromFile(phasefile,&N2,(void**)&phase_mp);CHKERRQ(ierr);
	if (N1 != N2) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Coordinate file and marker file have different lengths");
	}
	
	if (append == PETSC_FALSE) {
		start = 0;
		DataBucketSetSizes(db,N1,-1); // -ve buffer val retains old value //
	} else {
		DataBucketGetSizes(db,&n_mp_points,0,0);
		start = n_mp_points - 1;
		n_mp_points = n_mp_points + N1;
		DataBucketSetSizes(db,n_mp_points,-1); // -ve buffer val retains old value //
	}
	
	ierr = MaterialPointStdInsertBasic(db,da,start,N1,coords_mp,phase_mp);CHKERRQ(ierr);
	ierr = MaterialPointStdRemoval(db,start,N1,-1);CHKERRQ(ierr);
	
	DataBucketGetSizes(db,&n_mp_points,0,0);
	ierr = SwarmMPntStd_AssignUniquePointIdentifiers(PetscObjectComm((PetscObject)da),db,0,n_mp_points);CHKERRQ(ierr);
	
	
	free(coords_mp);
	free(phase_mp);
	
	PetscFunctionReturn(0);
}


/*
 Alternative way to checkpoint/store point data in fewer files.
 It's quite ugly, we store each pointdata.member in a Petsc Vec and use VecView()
*/
#undef __FUNCT__
#define __FUNCT__ "SwarmDataWriteToPetscVec"
PetscErrorCode SwarmDataWriteToPetscVec(DataBucket db,const char suffix[])
{
	PetscErrorCode ierr;
	int n_points;
	Vec point_field_data;
	BTruth found;
	PetscBool write_to_tgz = PETSC_FALSE;
	
	DataBucketGetSizes(db,&n_points,0,0);
	
	ierr = VecCreate(PETSC_COMM_WORLD,&point_field_data);CHKERRQ(ierr);
	ierr = VecSetSizes(point_field_data,(PetscInt)n_points,PETSC_DECIDE);CHKERRQ(ierr);
	ierr = VecSetFromOptions(point_field_data);CHKERRQ(ierr);
	
	/* 
	 Traverse fields
	   Traverse members
	     Traverse points -> pack entry into vector
			 Write out field
	*/
	/*
	DataBucketGetDataFields(db,&nfields,&fields);
	for (f=0; f<nfields; f++) {
		sprintf(field_member_name,"%s.",fields[f]->name);
		printf("[%d] %s \n",f,field_member_name);
	}
	*/
	
	/* ------------------- MPntStd_classname ------------------- */
	DataBucketQueryDataFieldByName(db,MPntStd_classname,&found);
	if (found == BTRUE) {
		ierr = SwarmDataWriteToPetscVec_MPntStd(db,suffix,point_field_data,write_to_tgz);CHKERRQ(ierr);
	}

	/* ------------------- MPntPStokes_classname ------------------- */
	DataBucketQueryDataFieldByName(db,MPntPStokes_classname,&found);
	if (found == BTRUE) {
		ierr = SwarmDataWriteToPetscVec_MPntPStokes(db,suffix,point_field_data,write_to_tgz);CHKERRQ(ierr);
	}
	
	/* ------------------- MPntPStokesPl_classname ------------------- */
	DataBucketQueryDataFieldByName(db,MPntPStokesPl_classname,&found);
	if (found == BTRUE) {
		ierr = SwarmDataWriteToPetscVec_MPntPStokesPl(db,suffix,point_field_data,write_to_tgz);CHKERRQ(ierr);
	}
	
	/* ------------------- MPntPEnergy_classname ------------------- */
	DataBucketQueryDataFieldByName(db,MPntPEnergy_classname,&found);
	if (found == BTRUE) {
		ierr = SwarmDataWriteToPetscVec_MPntPEnergy(db,suffix,point_field_data,write_to_tgz);CHKERRQ(ierr);
	}
	
	/* ------------------- QPntSurfCoefStokes_classname ------------------- */
	/* ------------------- QPntVolCoefStokes_classname ------------------- */
	/* ------------------- QPntVolCoefEnergy_classname ------------------- */

	
	
	ierr = VecDestroy(&point_field_data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
 It's expected the required fields have been registered in the DataBucket.
 This function will set the appropriate size and push the files associated with each field in the DataBucket. 
*/
#undef __FUNCT__
#define __FUNCT__ "SwarmDataLoadFromPetscVec"
PetscErrorCode SwarmDataLoadFromPetscVec(DataBucket db,const char suffix[])
{
	PetscErrorCode ierr;
	Vec point_field_data;
	BTruth found;
	PetscBool write_to_tgz = PETSC_FALSE;
	
	
	ierr = VecCreate(PETSC_COMM_WORLD,&point_field_data);CHKERRQ(ierr);
	ierr = VecSetFromOptions(point_field_data);CHKERRQ(ierr);
	
	/* ------------------- MPntStd_classname ------------------- */
	DataBucketQueryDataFieldByName(db,MPntStd_classname,&found);
	if (found == BTRUE) {
		ierr = SwarmDataLoadFromPetscVec_MPntStd(db,suffix,point_field_data,write_to_tgz);CHKERRQ(ierr);
	}
	
	/* ------------------- MPntPStokes_classname ------------------- */
	DataBucketQueryDataFieldByName(db,MPntPStokes_classname,&found);
	if (found == BTRUE) {
		ierr = SwarmDataLoadFromPetscVec_MPntPStokes(db,suffix,point_field_data,write_to_tgz);CHKERRQ(ierr);
	}
	
	/* ------------------- MPntPStokesPl_classname ------------------- */
	DataBucketQueryDataFieldByName(db,MPntPStokesPl_classname,&found);
	if (found == BTRUE) {
		ierr = SwarmDataLoadFromPetscVec_MPntPStokesPl(db,suffix,point_field_data,write_to_tgz);CHKERRQ(ierr);
	}
	
	/* ------------------- MPntPEnergy_classname ------------------- */
	DataBucketQueryDataFieldByName(db,MPntPEnergy_classname,&found);
	if (found == BTRUE) {
		ierr = SwarmDataLoadFromPetscVec_MPntPEnergy(db,suffix,point_field_data,write_to_tgz);CHKERRQ(ierr);
	}
	
	/* ------------------- QPntSurfCoefStokes_classname ------------------- */
	/* ------------------- QPntVolCoefStokes_classname ------------------- */
	/* ------------------- QPntVolCoefEnergy_classname ------------------- */
	

	ierr = VecDestroy(&point_field_data);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscScalar _PackLongIntToPetscScalar(long int val)
{
	PetscScalar pval;
	
	pval = (PetscScalar)val + 0.1;
	return pval;
}
PetscScalar _PackIntToPetscScalar(int val)
{
	PetscScalar pval;
	
	pval = (PetscScalar)val + 0.1;
	return pval;
}
PetscScalar _PackShortToPetscScalar(short val)
{
	PetscScalar pval;
	
	pval = (PetscScalar)val + 0.1;
	return pval;
}
PetscScalar _PackFloatToPetscScalar(float val)
{
	PetscScalar pval;
	pval = (PetscScalar)val;
	return pval;
}
PetscScalar _PackDoubleToPetscScalar(double val)
{
	PetscScalar pval;
	pval = (PetscScalar)val;
	return pval;
}


long int _UnPackPetscScalarToLongInt(PetscScalar pval)
{
	long int val;
	
	val = (long int)pval;
	return val;
}
int _UnPackPetscScalarToInt(PetscScalar pval)
{
	int val;
	
	val = (int)pval;
	return val;
}
short _UnPackPetscScalarToShort(PetscScalar pval)
{
	short val;
	
	val = (short)pval;
	return val;
}
float _UnPackPetscScalarToFloat(PetscScalar pval)
{
	float val;
	val = (float)pval;
	return val;
}
double _UnPackPetscScalarToDouble(PetscScalar pval)
{
	double val;
	val = (double)pval;
	return val;
}

/*
  All this functions below should be autogenerated.
  This will require updating my generation python scripts.
*/
#undef __FUNCT__
#define __FUNCT__ "SwarmDataWriteToPetscVec_MPntStd"
PetscErrorCode SwarmDataWriteToPetscVec_MPntStd(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz)
{
	PetscErrorCode ierr;
	int n_points;
	PetscScalar *LA_point_field_data;
	DataField pfield;
	char field_member_name[PETSC_MAX_PATH_LEN];
	char filename[PETSC_MAX_PATH_LEN];
	int m,d,p;
	MPntStd *points;
	
	DataBucketGetSizes(db,&n_points,0,0);
	DataBucketGetDataFieldByName(db,MPntStd_classname,&pfield);
	points = pfield->data;		
	
	for (m=0; m<MPntStd_nmembers; m++) {
		
		sprintf(field_member_name,"%s.%s",pfield->name,MPntStd_member_names[m]);
		printf("[%s] %s \n",MPntStd_classname,field_member_name);
		if (suffix) {
			sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
		} else {
			sprintf(filename,"swarm-%s.pvec",field_member_name);
		}
		
		switch (m) {
				
			case MPStd_point_index:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					long int val;
					PetscScalar pval;
					
					MPntStdGetField_point_index(&points[p],&val);
					pval = _PackLongIntToPetscScalar(val);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;
				
			case MPStd_global_coord:
				
				for (d=0; d<3; d++) {
					sprintf(field_member_name,"%s.%s%d",pfield->name,MPntStd_member_names[m],d);
					printf("[%s] %s \n",MPntStd_classname,field_member_name);
					if (suffix) {
						sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
					} else {
						sprintf(filename,"materialpoint-%s.pvec",field_member_name);
					}
					
					ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
					for (p=0; p<n_points; p++) {
						double *val;
						PetscScalar pval;
						
						MPntStdGetField_global_coord(&points[p],&val);
						pval = _PackDoubleToPetscScalar(val[d]);
						LA_point_field_data[p] = pval;
					}
					ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
					ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				}
				
				break;
				
			case MPStd_local_coord:
				
				for (d=0; d<3; d++) {
					sprintf(field_member_name,"%s.%s%d",pfield->name,MPntStd_member_names[m],d);
					printf("[%s] %s \n",MPntStd_classname,field_member_name);
					if (suffix) {
						sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
					} else {
						sprintf(filename,"materialpoint-%s.pvec",field_member_name);
					}
					
					ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
					for (p=0; p<n_points; p++) {
						double *val;
						PetscScalar pval;
						
						MPntStdGetField_local_coord(&points[p],&val);
						pval = _PackDoubleToPetscScalar(val[d]);
						LA_point_field_data[p] = pval;
					}
					ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
					ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				}
				
				break;
				
			case MPStd_phase_index:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					int val;
					PetscScalar pval;
					
					MPntStdGetField_phase_index(&points[p],&val);
					pval = _PackIntToPetscScalar(val);
					//printf("pack [%d] <region> %1.4e \n",p,pval);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;
				
			case MPStd_local_element_index:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					int val;
					PetscScalar pval;
					
					MPntStdGetField_local_element_index(&points[p],&val);
					pval = _PackIntToPetscScalar(val);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;
				
		}
	}
	
	PetscFunctionReturn(0);
}	

#undef __FUNCT__
#define __FUNCT__ "SwarmDataWriteToPetscVec_MPntPStokes"
PetscErrorCode SwarmDataWriteToPetscVec_MPntPStokes(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz)
{
	PetscErrorCode ierr;
	int n_points;
	PetscScalar *LA_point_field_data;
	DataField pfield;
	char field_member_name[PETSC_MAX_PATH_LEN];
	char filename[PETSC_MAX_PATH_LEN];
	int m,p;
	MPntPStokes *points;
	
	DataBucketGetSizes(db,&n_points,0,0);
	DataBucketGetDataFieldByName(db,MPntPStokes_classname,&pfield);
	points = pfield->data;		
	
	for (m=0; m<MPntPStokes_nmembers; m++) {
		
		sprintf(field_member_name,"%s.%s",pfield->name,MPntPStokes_member_names[m]);
		printf("[%s] %s \n",MPntPStokes_classname,field_member_name);
		if (suffix) {
			sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
		} else {
			sprintf(filename,"swarm-%s.pvec",field_member_name);
		}
		
		switch (m) {
				
			case MPPStk_eta_effective:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double val;
					PetscScalar pval;
					
					MPntPStokesGetField_eta_effective(&points[p],&val);
					pval = _PackDoubleToPetscScalar(val);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;

			case MPPStk_density:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double val;
					PetscScalar pval;
					
					MPntPStokesGetField_density(&points[p],&val);
					pval = _PackDoubleToPetscScalar(val);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;
				
		}
	}
	
	PetscFunctionReturn(0);
}	

#undef __FUNCT__
#define __FUNCT__ "SwarmDataWriteToPetscVec_MPntPStokesPl"
PetscErrorCode SwarmDataWriteToPetscVec_MPntPStokesPl(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz)
{
	PetscErrorCode ierr;
	int n_points;
	PetscScalar *LA_point_field_data;
	DataField pfield;
	char field_member_name[PETSC_MAX_PATH_LEN];
	char filename[PETSC_MAX_PATH_LEN];
	int m,p;
	MPntPStokesPl *points;
	
	DataBucketGetSizes(db,&n_points,0,0);
	DataBucketGetDataFieldByName(db,MPntPStokesPl_classname,&pfield);
	points = pfield->data;		
	
	for (m=0; m<MPntPStokesPl_nmembers; m++) {
		
		sprintf(field_member_name,"%s.%s",pfield->name,MPntPStokesPl_member_names[m]);
		printf("[%s] %s \n",MPntPStokesPl_classname,field_member_name);
		if (suffix) {
			sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
		} else {
			sprintf(filename,"swarm-%s.pvec",field_member_name);
		}
		
		switch (m) {
				
			case MPPStkPl_plastic_strain:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					float val;
					PetscScalar pval;
					
					MPntPStokesPlGetField_plastic_strain(&points[p],&val);
					pval = _PackFloatToPetscScalar(val);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;
				
			case MPPStkPl_yield_indicator:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					short val;
					PetscScalar pval;
					
					MPntPStokesPlGetField_yield_indicator(&points[p],&val);
					pval = _PackShortToPetscScalar(val);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;
				
		}
	}
	
	PetscFunctionReturn(0);
}	

#undef __FUNCT__
#define __FUNCT__ "SwarmDataWriteToPetscVec_MPntPEnergy"
PetscErrorCode SwarmDataWriteToPetscVec_MPntPEnergy(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz)
{
	PetscErrorCode ierr;
	int n_points;
	PetscScalar *LA_point_field_data;
	DataField pfield;
	char field_member_name[PETSC_MAX_PATH_LEN];
	char filename[PETSC_MAX_PATH_LEN];
	int m,p;
	MPntPEnergy *points;
	
	DataBucketGetSizes(db,&n_points,0,0);
	DataBucketGetDataFieldByName(db,MPntPEnergy_classname,&pfield);
	points = pfield->data;		
	
	for (m=0; m<MPntPEnergy_nmembers; m++) {
		
		sprintf(field_member_name,"%s.%s",pfield->name,MPntPEnergy_member_names[m]);
		printf("[%s] %s \n",MPntPEnergy_classname,field_member_name);
		if (suffix) {
			sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
		} else {
			sprintf(filename,"swarm-%s.pvec",field_member_name);
		}
		
		switch (m) {
				
			case MPPEgy_diffusivity:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double val;
					PetscScalar pval;
					
					MPntPEnergyGetField_diffusivity(&points[p],&val);
					pval = _PackDoubleToPetscScalar(val);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;
				
			case MPPEgy_heat_source:
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double val;
					PetscScalar pval;
					
					MPntPEnergyGetField_heat_source(&points[p],&val);
					pval = _PackDoubleToPetscScalar(val);
					LA_point_field_data[p] = pval;
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				ierr = DMDAWriteVectorToFile(point_field_data,filename,write_to_tgz);CHKERRQ(ierr);
				break;
				
		}
	}
	
	PetscFunctionReturn(0);
}	



#undef __FUNCT__
#define __FUNCT__ "SwarmDataLoadFromPetscVec_MPntStd"
PetscErrorCode SwarmDataLoadFromPetscVec_MPntStd(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz)
{
	PetscErrorCode ierr;
	PetscInt n_points;
	PetscScalar *LA_point_field_data;
	DataField pfield;
	char field_member_name[PETSC_MAX_PATH_LEN];
	char filename[PETSC_MAX_PATH_LEN];
	int m,d,p;
	MPntStd *points;
	int field_n_members;
	const char *field_classname;
	const char **field_member_names;
	

	/* member data */
	field_n_members    = MPntStd_nmembers;
	field_classname    = MPntStd_classname;
	field_member_names = MPntStd_member_names;
	
	DataBucketGetDataFieldByName(db,field_classname,&pfield);
	points = pfield->data;		
	
	
	for (m=0; m<field_n_members; m++) {
		
		sprintf(field_member_name,"%s.%s",pfield->name,field_member_names[m]);
		printf("LOAD: [%s] %s \n",field_classname,field_member_name);
		if (suffix) {
			sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
		} else {
			sprintf(filename,"swarm-%s.pvec",field_member_name);
		}
		
		switch (m) {
				
			case MPStd_point_index:
				
				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					long int    val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToLongInt(pval);
					MPntStdSetField_point_index(&points[p],val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;
				
			case MPStd_global_coord:
				
				for (d=0; d<3; d++) {
					sprintf(field_member_name,"%s.%s%d",pfield->name,field_member_names[m],d);
					printf("  ... [%s] %s \n",field_classname,field_member_name);
					if (suffix) {
						sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
					} else {
						sprintf(filename,"materialpoint-%s.pvec",field_member_name);
					}
					
					ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
					ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
					DataBucketSetSizes(db,(int)n_points,-1);

					ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
					for (p=0; p<n_points; p++) {
						double      vald,*val,valset[3];
						PetscScalar pval;
						
						pval = LA_point_field_data[p];
						vald = _UnPackPetscScalarToDouble(pval);
						MPntStdGetField_global_coord(&points[p],&val);
						memcpy(valset,val,sizeof(double)*3);
						valset[d] = vald;
						//printf("  [%d] %1.4e %1.4e %1.4e \n",d,valset[0],valset[1],valset[2]);
						MPntStdSetField_global_coord(&points[p],valset);						
					}
					ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				}
				
				break;
				
			case MPStd_local_coord:
				
				for (d=0; d<3; d++) {
					sprintf(field_member_name,"%s.%s%d",pfield->name,field_member_names[m],d);
					printf("  ... [%s] %s \n",field_classname,field_member_name);
					if (suffix) {
						sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
					} else {
						sprintf(filename,"materialpoint-%s.pvec",field_member_name);
					}
					
					ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
					ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
					DataBucketSetSizes(db,(int)n_points,-1);
					
					ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
					for (p=0; p<n_points; p++) {
						double      vald,*val,valset[3];
						PetscScalar pval;
						
						pval = LA_point_field_data[p];
						vald = _UnPackPetscScalarToDouble(pval);
						MPntStdGetField_local_coord(&points[p],&val);
						memcpy(valset,val,sizeof(double)*3);
						valset[d] = vald;
						MPntStdSetField_local_coord(&points[p],valset);
					}
					ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				}
				
				break;
				
			case MPStd_phase_index:

				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					int         val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToInt(pval);
					MPntStdSetField_phase_index(&points[p],val);
					//printf("[%d] <region> %d \n",p,val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;
				
			case MPStd_local_element_index:
				
				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					int         val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToInt(pval);
					MPntStdSetField_local_element_index(&points[p],val);
					//printf("[%d] <wil> %d \n",p,val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;
				
		}
	}
	
	PetscFunctionReturn(0);
}	

#undef __FUNCT__
#define __FUNCT__ "SwarmDataLoadFromPetscVec_MPntPStokes"
PetscErrorCode SwarmDataLoadFromPetscVec_MPntPStokes(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz)
{
	PetscErrorCode ierr;
	PetscInt n_points;
	PetscScalar *LA_point_field_data;
	DataField pfield;
	char field_member_name[PETSC_MAX_PATH_LEN];
	char filename[PETSC_MAX_PATH_LEN];
	PetscInt m,p;
	MPntPStokes *points;
	int field_n_members;
	const char *field_classname;
	const char **field_member_names;
	
	
	/* member data */
	field_n_members    = MPntPStokes_nmembers;
	field_classname    = MPntPStokes_classname;
	field_member_names = MPntPStokes_member_names;
	
	DataBucketGetDataFieldByName(db,field_classname,&pfield);
	points = pfield->data;		
	
	
	for (m=0; m<field_n_members; m++) {
		
		sprintf(field_member_name,"%s.%s",pfield->name,field_member_names[m]);
		printf("LOAD: [%s] %s \n",field_classname,field_member_name);
		if (suffix) {
			sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
		} else {
			sprintf(filename,"swarm-%s.pvec",field_member_name);
		}
		
		switch (m) {
				
			case MPPStk_eta_effective:
				
				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double      val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToDouble(pval);
					MPntPStokesSetField_eta_effective(&points[p],val);
					//printf("[%d] <eta> %1.4e \n",p,val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;

			case MPPStk_density:
				
				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double      val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToDouble(pval);
					MPntPStokesSetField_density(&points[p],val);
					//printf("[%d] <rho> %1.4e \n",p,val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;
				
		}
	}
	
	PetscFunctionReturn(0);
}	

#undef __FUNCT__
#define __FUNCT__ "SwarmDataLoadFromPetscVec_MPntPStokesPl"
PetscErrorCode SwarmDataLoadFromPetscVec_MPntPStokesPl(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz)
{
	PetscErrorCode ierr;
	PetscInt n_points;
	PetscScalar *LA_point_field_data;
	DataField pfield;
	char field_member_name[PETSC_MAX_PATH_LEN];
	char filename[PETSC_MAX_PATH_LEN];
	PetscInt m,p;
	MPntPStokesPl *points;
	int field_n_members;
	const char *field_classname;
	const char **field_member_names;
	
	
	/* member data */
	field_n_members    = MPntPStokesPl_nmembers;
	field_classname    = MPntPStokesPl_classname;
	field_member_names = MPntPStokesPl_member_names;
	
	DataBucketGetDataFieldByName(db,field_classname,&pfield);
	points = pfield->data;		
	
	
	for (m=0; m<field_n_members; m++) {
		
		sprintf(field_member_name,"%s.%s",pfield->name,field_member_names[m]);
		printf("LOAD: [%s] %s \n",field_classname,field_member_name);
		if (suffix) {
			sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
		} else {
			sprintf(filename,"swarm-%s.pvec",field_member_name);
		}
		
		switch (m) {
				
			case MPPStkPl_plastic_strain:
				
				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double      val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToDouble(pval);
					MPntPStokesPlSetField_plastic_strain(&points[p],val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;
				
			case MPPStkPl_yield_indicator:
				
				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double      val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToDouble(pval);
					MPntPStokesPlSetField_yield_indicator(&points[p],val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;
				
		}
	}
	
	PetscFunctionReturn(0);
}	

#undef __FUNCT__
#define __FUNCT__ "SwarmDataLoadFromPetscVec_MPntPEnergy"
PetscErrorCode SwarmDataLoadFromPetscVec_MPntPEnergy(DataBucket db,const char suffix[],Vec point_field_data,PetscBool write_to_tgz)
{
	PetscErrorCode ierr;
	PetscInt n_points;
	PetscScalar *LA_point_field_data;
	DataField pfield;
	char field_member_name[PETSC_MAX_PATH_LEN];
	char filename[PETSC_MAX_PATH_LEN];
	PetscInt m,p;
	MPntPEnergy *points;
	int field_n_members;
	const char *field_classname;
	const char **field_member_names;
	
	
	/* member data */
	field_n_members    = MPntPEnergy_nmembers;
	field_classname    = MPntPEnergy_classname;
	field_member_names = MPntPEnergy_member_names;
	
	DataBucketGetDataFieldByName(db,field_classname,&pfield);
	points = pfield->data;		
	
	
	for (m=0; m<field_n_members; m++) {
		
		sprintf(field_member_name,"%s.%s",pfield->name,field_member_names[m]);
		printf("LOAD: [%s] %s \n",field_classname,field_member_name);
		if (suffix) {
			sprintf(filename,"%s-%s.pvec",suffix,field_member_name);
		} else {
			sprintf(filename,"swarm-%s.pvec",field_member_name);
		}
		
		switch (m) {
				
			case MPPEgy_diffusivity:
				
				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double      val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToDouble(pval);
					MPntPEnergySetField_diffusivity(&points[p],val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;
				
			case MPPEgy_heat_source:
				
				ierr = VecLoadFromFile(point_field_data,filename);CHKERRQ(ierr);
				ierr = VecGetLocalSize(point_field_data,&n_points);CHKERRQ(ierr);
				DataBucketSetSizes(db,(int)n_points,-1);
				
				ierr = VecGetArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				for (p=0; p<n_points; p++) {
					double      val;
					PetscScalar pval;
					
					pval = LA_point_field_data[p];
					val = _UnPackPetscScalarToDouble(pval);
					MPntPEnergySetField_heat_source(&points[p],val);
				}
				ierr = VecRestoreArray(point_field_data,&LA_point_field_data);CHKERRQ(ierr);
				break;
				
		}
	}
	
	PetscFunctionReturn(0);
}	


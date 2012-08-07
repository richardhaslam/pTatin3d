

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <petsc.h>

#include "dmda_element_q2p1.h"
#include "swarm_fields.h"
#include "MPntStd_def.h"
#include "material_point_std_utils.h"
#include "material_point_point_location.h"


#undef __FUNCT__  
#define __FUNCT__ "MarkerCoordinatesLoadFromFile"
PetscErrorCode MarkerCoordinatesLoadFromFile(const char name[],long int *length,double **coords)
{
	FILE *fp = NULL;
	int vtk_data_type;
	long int n_markers;
	double *data;
	int p;
	char line[256];
	
	fp = fopen(name,"rb");
	if (fp == NULL) {
		SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"File %s not found",name);
	}
	
	fgets(line,255,fp);
	vtk_data_type = atoi( line );
	fgets(line,255,fp);
	n_markers = atol( line );

	*length = n_markers;
	data = malloc( sizeof(double)*3*n_markers );
	memset(data,0,sizeof(double)*3*n_markers);

	for (p=0; p<n_markers; p++) {
		fread(&data[3*p],sizeof(double),3,fp);
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
	int p;
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
	DataBucketSetSizes(db,end+1,-1); // -ve buffer val retains old value //
	
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
	MPntStd *local_list;
	const PetscInt *elnidx_u;
	DataField PField_std;
	PetscErrorCode ierr;
	
	/* setup for coords */
	ierr = DMDAGetCoordinateDA(da,&cda);CHKERRQ(ierr);
	ierr = DMDAGetGhostedCoordinates(da,&gcoords);CHKERRQ(ierr);
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
		
		test_p.pid     = 0;
		test_p.coor[0] = coords_mp[3*p  ];
		test_p.coor[1] = coords_mp[3*p+1];
		test_p.coor[2] = coords_mp[3*p+2];
		test_p.phase   = phase_mp[p];
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
	long int p,start;
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
	ierr = SwarmMPntStd_AssignUniquePointIdentifiers(((PetscObject)da)->comm,db,0,n_mp_points);CHKERRQ(ierr);
	
	
	free(coords_mp);
	free(phase_mp);
	
	PetscFunctionReturn(0);
}


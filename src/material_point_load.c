

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <petsc.h>

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


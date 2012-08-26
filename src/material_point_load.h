


#ifndef __ptatin_material_point_load_h__
#define __ptatin_material_point_load_h__

PetscErrorCode MarkerCoordinatesLoadFromFile(const char name[],long int *length,double **coords);
PetscErrorCode MarkerScalarFieldLoadFromFile(const char name[],long int *length,void **field);

PetscErrorCode MaterialPointStdRemoval(DataBucket db,long int start,long int npoints,const int wil_key);
PetscErrorCode MaterialPointStdInsertBasic(DataBucket db,DM da,long int start,long int npoints,double coords_mp[],int phase_mp[]);
PetscErrorCode MaterialPointDataBasicLoadIntoListFromFile(DataBucket db,DM da,PetscBool append,const char coordfile[],const char phasefile[]);

#endif

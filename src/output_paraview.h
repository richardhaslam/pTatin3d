
#ifndef __ptatin_output_paraview_h__
#define __ptatin_output_paraview_h__

#include "petsc.h"
#include "petscdm.h"

PetscErrorCode pTatinGenerateParallelVTKName(const char prefix[],const char suffix[],char **name);
PetscErrorCode pTatinGenerateVTKName(const char prefix[],const char suffix[],char **name);

PetscErrorCode ParaviewPVDOpen(const char pvdfilename[]);
PetscErrorCode ParaviewPVDAppend(const char pvdfilename[],double time,const char datafile[], const char DirectoryName[]);

PetscErrorCode pTatinOutputParaViewMeshVelocityPressure(DM pack,Vec X,const char path[],const char prefix[]);
PetscErrorCode pTatinOutputMeshVelocityPressureVTS_v0(DM pack,Vec X,const char name[]);
PetscErrorCode pTatinOutputMeshVelocityPressurePVTS(DM pack,const char prefix[],const char name[]);
PetscErrorCode DAQ2PieceExtendForGhostLevelZero( FILE *vtk_fp, int indent_level, DM dau, const char local_file_prefix[] );

#endif
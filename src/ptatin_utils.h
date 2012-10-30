
#ifndef __ptatin3d_ptatin_utils_h__
#define __ptatin3d_ptatin_utils_h__

PetscErrorCode pTatinCreateDirectory(const char dirname[]);
PetscErrorCode pTatinWriteOptionsFile(const char filename[]);
void pTatinGenerateFormattedTimestamp(char date_time[]);
void FileExists(const char *fname,int *exists);

#endif

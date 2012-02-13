
#ifndef __ptatin3d_ptatin_utils_h__
#define __ptatin3d_ptatin_utils_h__

//#include "errno.h"
//#include "sys/types.h"
#include "sys/stat.h"

PetscErrorCode pTatinCreateDirectory(const char dirname[]);
PetscErrorCode pTatinWriteOptionsFile(const char filename[]);

#endif

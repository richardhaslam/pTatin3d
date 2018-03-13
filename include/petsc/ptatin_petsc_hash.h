#if !defined(__PTATIN_PETSC_HASH_H__)
#define __PTATIN_PETSC_HASH_H__

/* version check */
#include "petscversion.h"

#if ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR == 7))
  #include "src/sys/utils/hash.h"
#else
  #error "pTatin provided private petsc header for src/sys/utils/hash.h is only valid for PETSc v3.7"
#endif

#endif


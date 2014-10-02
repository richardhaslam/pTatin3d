#
#  -- pTatin Optimization flags --
#

# Written to work on the Piz Daint Cray XC30 system
# with the PrgEnv0gnu module
# and the zlib module

# Set your own PETSC_DIR and PETSC_ARCH 
# PETSC_ARCH was arch-cray-xc30-gnu at the time of this writing (but you will have picked the name)
ifndef PETSC_DIR
  $(error PETSC_DIR must be defined as an environment variable)
endif
ifndef PETSC_ARCH
  $(error PETSC_ARCH must be defined as an environment variable)
endif

export LIBZ_LIB=-lz

# flag to allow GNU AVX/SSE intrinsics
#export TATIN_CFLAGS=-hgnu
CONFIG_AVX=y

#
#
#  -- pTatin external packages --
#  User should select which packages are present
#    ** To ensure consistent compilation, between MPI, PETSc and ptatin3d, **
#    ** please compile external packages using the PETSc compiler          **
#
#
PTATIN_CONTAINS_SPMA = 0
PTATIN_CONTAINS_FASTSCAPE_V3 = 0

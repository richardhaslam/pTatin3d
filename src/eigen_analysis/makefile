## ==============================================================================
##
##     makefile for pTatin3d
##
## ==============================================================================



TATIN_MODELS		=
TATIN_DRIVERS		= \
			ptatin_eigen_analysis

export EXT_TATIN_INC	=	-I${PTATIN_DIR}/src
export EXT_TATIN_LIB	=	-L${PTATIN_DIR}/${PETSC_ARCH}/lib -lptatin3d -lptatin3dmodels -lptatin3d

# slepc
include ${SLEPC_DIR}/lib/slepc/conf/slepc_variables
include ${SLEPC_DIR}/lib/slepc/conf/slepc_rules
export EXT_SLEPC_INC = ${SLEPC_INCLUDE}
export EXT_SLEPC_LIB = ${SLEPC_LIB}

# petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
export EXT_PETSC_INC = -I${PETSC_DIR}/include ${PETSC_CC_INCLUDES}
#export EXT_PETSC_LIB = ${PETSC_LIB}

# libz
export LIBZ_DIR = /usr/include
export EXT_LIBZ_INC = -I${LIBZ_DIR}
export EXT_LIBZ_LIB = -lz

# combine every thing here
export TATIN_INC	=	-I. ${EXT_TATIN_INC} ${EXT_SLEPC_INC} ${EXT_PETSC_INC} ${EXT_LIBZ_INC}
export TATIN_LIB	=	${EXT_TATIN_LIB} ${EXT_SLEPC_LIB} ${EXT_LIBZ_LIB}

## Set flags here, or run >make TATIN_CFLAGS='-g -O0'
export TATIN_CFLAGS_DEBUG				= -O0 -g -Wno-unused-variable -Wall
export TATIN_CFLAGS_O2_OSX			= -std=gnu99 -O2 -Wall -Wno-unused-variable -Wstrict-aliasing -fstrict-aliasing -funroll-loops
export TATIN_CFLAGS_O3_OSX			= -std=gnu99 -O3 -g -fast -ftree-vectorize -ftree-vectorizer-verbose=2  -Wall -Wstrict-aliasing -Wno-unused-variable -fstrict-aliasing -funroll-loops 
export TATIN_CFLAGS_O3_DIONYSOS = -std=gnu99 -O3 -g -ftree-vectorize -ftree-vectorizer-verbose=2  -Wall -Wstrict-aliasing -Wno-unused-variable -fstrict-aliasing -funroll-loops
export TATIN_CFLGAS_O2_BGP      = -O2 -qarch=auto -qtune=auto -qcache=auto -qhot -qipa
export TATIN_CFLGAS_O3_BGP			= -O3 -qarch=auto -qtune=auto -qcache=auto -qhot -qipa
export TATIN_CFLAGS_O3_OCTO			= -O3 -std=gnu99 -mtune=native -march=native -funroll-loops -ffast-math


export TATIN_CFLAGS		= ${TATIN_CFLAGS_DEBUG}
export TATIN_CC				= ${PCC}


all: 
	-@echo "<<<<<<<<==================[[ pTatin-3d ]]==================>>>>>>>>"
	-@echo ${PWD}
	
	-@${MAKE} drivers


drivers:
	-@echo ""
	-@echo "*** Building pTatin executables ***"
	@for t in ${TATIN_DRIVERS} ; do \
		${MAKE} $${t}.o; \
		${PCC} ${TATIN_CFLAGS} -o $${t}.app $${t}.o ${TATIN_INC} ${TATIN_LIB}; \
	done

%.o: %.c
	-@echo ""
	${PCC} ${TATIN_CFLAGS} -c $*.c ${TATIN_INC}

clean_all:
	-@echo ""
	-@echo "*** Removing ptatin applications ***"
	-@rm -f *.app *.o


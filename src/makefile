## ==============================================================================
##  make for tatin
##
##
## ==============================================================================

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export TATIN_DIR	=	${PWD}

export T_SRC	=	\
			dmda_update_coords.c \
			dmda_project_coords.c \
			dmda_compare.c \
			dmda_redundant.c \
			dmda_remesh.c \
			dmda_duplicate.c \
			dmda_view_petscvtk.c \
			dmda_checkpoint.c \

export T_OBJ	=	$(T_SRC:.c=.o)

TEST		=	test_dmda_update_coords \
					test_dmda_redundant \
					test_dmda_remesh \
					test_dmda_duplicate \
					test_dmda_checkpoint \

export INC		=	-I. -I${PETSC_DIR}/include
export LIB		=	${PETSC_LIB}

## or run >make CFLAGS='-g -O0'
export CFLAGS_DEBUG	=	-O0 -g -Wall
export CFLAGS_O2_OSX	=	-std=gnu99 -O2 -Wall -Wstrict-aliasing -fstrict-aliasing -funroll-loops
export CFLAGS_O3_OSX	=	-std=gnu99 -O3 -g -fast -ftree-vectorize -ftree-vectorizer-verbose=2  -Wall -Wstrict-aliasing -fstrict-aliasing -funroll-loops
export CFLAGS		= 	${CFLAGS_DEBUG}

all: 
	${MAKE} base_objects
	${MAKE} libs
	${MAKE} test

base_objects:
	-@echo ""
	-@echo "*** Creating ptatin object files ***"
	@for a in ${T_OBJ} ; do \
		${MAKE} $$a; \
	done

libs: base_objects
	-@echo ""
	-@echo "*** Creating libptatin3d ***"
	-@rm -f libptatin3d.a
	-@ar rcs libptatin3d.a ${T_OBJ}

test: libs
	-@echo ""
	-@echo "*** Building tests ***"
	@for t in ${TEST} ; do \
		${MAKE} $${t}.o; \
		${PCC} ${CFLAGS} -o $${t}.app $${t}.o ${T_OBJ} ${INC} ${LIB}; \
	done

%.o: %.c ${T_SRC}
	-@echo ""
	${PCC} ${CFLAGS} -c $*.c ${INC}



clean_all:
	-@echo ""
	-@echo "*** Removing base_objects ***"
	-@rm -f ${T_OBJ}
	-@echo "*** Removing tests object files ***"
	-@rm -f *.o
	-@echo "*** Removing ptatin applications ***"
	-@rm -f *.app
	-@echo "*** Removing ptatin3d.a ***"
	-@rm -f libptatin3d.a


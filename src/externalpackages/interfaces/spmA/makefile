
include ${TATIN_DIR}/makefile.arch

## List source files here
INTERFACE_SRC = \
	ptatin_spmA_interface.c
	

## USER SHOULD NOT NEED TO EDIT ANYTHING BELOW
INTERFACE_OBJ	=	$(INTERFACE_SRC:.c=.o)

INTERFACE_PATH				= ${PWD}
TATIN_INTERFACE_INC		=	-I${INTERFACE_PATH} -I${TATIN_DIR} ${TATIN_INC} ${TATIN_EXTERNAL_INC}

all:
	-@echo "    + Compiling interface: " ${INTERFACE_PATH}
	-@for a in ${INTERFACE_OBJ} ; do \
		${MAKE} $$a; \
	done

%.o: %.c
	-@echo ""
	${TATIN_CC} ${TATIN_CFLAGS} ${TATIN_PACKAGES} -c $*.c ${TATIN_INTERFACE_INC}


addtolib:
	-@echo "    + Adding interface object files from external packages "
	ar rcs ${TATIN_DIR}/libptatin3d.a ${INTERFACE_OBJ}


clean_all:
	-@rm -f ${INTERFACE_OBJ}



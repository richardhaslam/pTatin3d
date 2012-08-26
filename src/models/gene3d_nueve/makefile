
## List source files here
MODEL_SRC = \
	model_ops_gene3d_nueve.c \
	

## USER SHOULD NOT NEED TO EDIT ANYTHING BELOW
MODEL_OBJ	=	$(MODEL_SRC:.c=.o)

MODEL_PATH				= `pwd`
TATIN_MODEL_INC		=	-I${MODEL_PATH} -I${TATIN_DIR} ${TATIN_INC}

%.o: %.c
	-@echo ""
	${TATIN_CC} ${TATIN_CFLAGS} -c $*.c ${TATIN_MODEL_INC}

all:
	-@echo "    + Compiling model: " `pwd`
#	-@echo ${TATIN_DIR}
#	-@echo ${TATIN_MODEL_INC}
	
	for a in ${MODEL_OBJ} ; do \
		${MAKE} $$a; \
	done
	-@echo "    + adding objects: " ${MODEL_OBJ}

addtolib:
	ar rcs ${TATIN_MODEL_LIB} ${MODEL_OBJ}
#	-@rm -f ${MODEL_OBJ}

clean_all:
	-@rm -f ${MODEL_OBJ}



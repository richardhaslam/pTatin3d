
## List source files here
MATC_SRC = \
	MaterialConst_MaterialType_def.c \
	MaterialConst_DensityConst_def.c \
	MaterialConst_DensityBoussinesq_def.c \
	MaterialConst_ViscosityConst_def.c \
	MaterialConst_ViscosityZ_def.c \
	MaterialConst_ViscosityArrh_def.c \
	MaterialConst_ViscosityFK_def.c \
	MaterialConst_SoftLin_def.c \
	MaterialConst_SoftExpo_def.c \
	MaterialConst_PlasticMises_def.c \
	MaterialConst_PlasticDP_def.c
	

## USER SHOULD NOT NEED TO EDIT ANYTHING BELOW
MATC_OBJ	=	$(MATC_SRC:.c=.o)

MATC_PATH				= `pwd`
TATIN_MATC_INC		=	-I. -I${TATIN_DIR} ${TATIN_INC}

%.o: %.c
	-@echo ""
	${TATIN_CC} ${TATIN_CFLAGS} -c $*.c ${TATIN_MATC_INC}

all:
	-@echo "    + Compiling model: " `pwd`
	
	-@for a in ${MATC_OBJ} ; do \
		${MAKE} $$a; \
	done
	-@echo "    + adding objects: " ${MATC_OBJ}

clean_all:
	-@rm -f ${MATC_OBJ}




## LIST MODEL DIRECTORIES HERE
TATIN_MODEL_DIR = \
	template \
	viscous_sinker \
	indentor \
	rift3D \
	rift3D_T \
	gene3d \
	gene3d_nueve \
	sierra \
	folding \
	folding2d \
	advdiff_example\
	basin_comp\
	wrench_fold\
	fault_fold \
	delamination \
	riftrh \
	geomod2008 \
	fault_fold_plastic \
	multilayer_folding \
	submarinelavaflow \
	ex_subduction \
	iplus_slab_plume \
	subduction_initiation2d \
	convection2d \
	thermal_shearband \
	slab_detachment3d \
	pinch_and_swell \
	plasticdemo \
    mms

TATIN_MODEL_SRC = \
	ptatin_models_reg.c

## DONT NOT EDIT LINES BELOW
export TATIN_MODEL_LIB = ${TATIN_DIR}/models/libptatin3dmodels.a

TATIN_MODEL_INC		=	-I${TATIN_DIR} ${TATIN_INC}
%.o: %.c
	${TATIN_CC} ${TATIN_CFLAGS} -c $*.c ${TATIN_MODEL_INC}

all:
	-@rm -f ${TATIN_MODEL_LIB}
	-@${MAKE} build
	-@${MAKE} modellib

build:
	@for m in ${TATIN_MODEL_DIR} ; do \
		${MAKE} -C $${m}/; \
	done

modellib: $(TATIN_MODEL_SRC:%.c=%.o)
	-@echo ""
	-@echo "***  Assembling model library using object files:" 
	ar rcs ${TATIN_MODEL_LIB} $^
	@for m in ${TATIN_MODEL_DIR} ; do \
		${MAKE} addtolib -C $${m}/; \
	done

clean_all:
	@for m in ${TATIN_MODEL_DIR} ; do \
		${MAKE} clean_all -C $${m}/; \
	done

	-@rm -f *.o
	-@rm -f ${TATIN_MODEL_LIB}


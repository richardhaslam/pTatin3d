
## LIST MODEL DIRECTORIES HERE
TATIN_MODEL_DIR = \
	template \
	viscous_sinker \
	indentor \
	rift3D \
	gene3d \
	gene3d_nueve \
	sierra \
	folding \
	folding2d \
	advdiff_example


## DONT NOT EDIT LINES BELOW
export TATIN_MODEL_LIB = ${TATIN_DIR}/models/libptatin3dmodels.a

all:
	-@rm -f ${TATIN_MODEL_LIB}
	-@${MAKE} build
	-@${MAKE} modellib

build:
	@for m in ${TATIN_MODEL_DIR} ; do \
		${MAKE} -C $${m}/; \
	done

modellib:
	-@echo ""
	-@echo "  Assembling model library using object files:" 
	@for m in ${TATIN_MODEL_DIR} ; do \
		${MAKE} addtolib -C $${m}/; \
	done

	
clean_all:
	@for m in ${TATIN_MODEL_DIR} ; do \
		${MAKE} clean_all -C $${m}/; \
	done

	-@rm -f *.o
	-@rm -f ${TATIN_MODEL_LIB}


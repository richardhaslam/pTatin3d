
## LIST MODEL DIRECTORIES HERE
TATIN_MODEL_DIR = \
	template \
	viscous_sinker \
	gene3D_orig \
	indentor \
	sierra


## DONT NOT EDIT LINES BELOW
all:
	-@${MAKE} build
	-@${MAKE} modellib

build:
	@for m in ${TATIN_MODEL_DIR} ; do \
		${MAKE} -C $${m}/; \
	done

modellib:
	-@echo ""
	-@echo "  Assembling model library using object files:" 
	-@ls *.o
	ar rcs libptatin3dmodels.a *.o
	-@rm *.o
	
clean_all:
	@for m in ${TATIN_MODEL_DIR} ; do \
		${MAKE} clean_all -C $${m}/; \
	done

	-@rm -f *.o
	-@rm libptatin3dmodels.a


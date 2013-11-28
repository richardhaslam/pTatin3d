
include ${TATIN_DIR}/makefile.arch

# ----------------------------------------------------------------------------
# SPMA:
# Simple C landscape evolution model. 
# Serves as a template for coupling surface process models. Contained within ptatin3d.
#

ifeq (${PTATIN_CONTAINS_SPMA},1)
	# add the compilation flags specific for this package
	TATIN_PACKAGES += -DPTATIN_HAVE_SPMA

	SPMA_DIR = ${TATIN_DIR}/externalpackages/spmA

	TATIN_EXTERNAL_INC += -I${SPMA_DIR}

	# identify all source files defined within this package
	TATIN_EXTERNAL_OBJ += ${SPMA_DIR}/spmA.o

	# OR
	# identify all libraries defined within this package
	TATIN_EXTERNAL_LIB +=

endif
# Append location of interface directory
TATIN_INTERFACE_DIR += interfaces/spmA


# ----------------------------------------------------------------------------
# SPM FastScapeV3:
# Landscape evolution model written by Jean Braun.
#

ifeq (${PTATIN_CONTAINS_FASTSCAPE_V3},1)
	# add the compilation flags specific for this package
	TATIN_PACKAGES += -DPTATIN_HAVE_FASTSCAPE_V3

	SPM_FSCAPE_DIR = ${TATIN_DIR}/externalpackages/FastScape_V3_lib

	TATIN_EXTERNAL_INC += -I${SPM_FSCAPE_DIR}

	# identify all source files defined within this package
	TATIN_EXTERNAL_OBJ += ${SPM_FSCAPE_DIR}/Concavity_f90.o \
			${SPM_FSCAPE_DIR}/Curvature_f90.o \
			${SPM_FSCAPE_DIR}/DEM_f90.o \
			${SPM_FSCAPE_DIR}/Demoulin_f90.o \
			${SPM_FSCAPE_DIR}/libFastScape_f90.o \
			${SPM_FSCAPE_DIR}/Slope_f90.o \
			${SPM_FSCAPE_DIR}/SteepnessIndex_f90.o \
			${SPM_FSCAPE_DIR}/bmp_f90.o \
			${SPM_FSCAPE_DIR}/compute_sediment_flux_f90.o \
			${SPM_FSCAPE_DIR}/diffusion_f90.o \
			${SPM_FSCAPE_DIR}/find_donors_f90.o \
			${SPM_FSCAPE_DIR}/find_receiver_f90.o \
			${SPM_FSCAPE_DIR}/flexure_f90.o \
			${SPM_FSCAPE_DIR}/fluvial_erosion_f90.o \
			${SPM_FSCAPE_DIR}/four1_f77.o \
			${SPM_FSCAPE_DIR}/initial_topography_f90.o \
			${SPM_FSCAPE_DIR}/initialize_f90.o \
			${SPM_FSCAPE_DIR}/interpol_routines_f90.o \
			${SPM_FSCAPE_DIR}/interpolate_f90.o \
			${SPM_FSCAPE_DIR}/local_minima_f90.o \
			${SPM_FSCAPE_DIR}/metric_f90.o \
			${SPM_FSCAPE_DIR}/module_definitions_f90.o \
			${SPM_FSCAPE_DIR}/orography_f90.o \
			${SPM_FSCAPE_DIR}/precipitation_f90.o \
			${SPM_FSCAPE_DIR}/realft_f77.o \
			${SPM_FSCAPE_DIR}/scanfile_f90.o \
			${SPM_FSCAPE_DIR}/set_bc_f90.o \
			${SPM_FSCAPE_DIR}/sinft_f77.o \
			${SPM_FSCAPE_DIR}/tridag_f90.o \
			${SPM_FSCAPE_DIR}/uplift_f90.o
endif
# Append location of interface directory
TATIN_INTERFACE_DIR += interfaces/FastScape_V3_lib





# expose variables to interface makefiles
export TATIN_PACKAGES
export TATIN_EXTERNAL_INC
export TATIN_EXTERNAL_OBJ
export TATIN_EXTERNAL_LIB

all:
	-@echo "External packages:"
	-@echo "    " ${TATIN_PACKAGES}
	-@echo "External include paths:"
	-@echo "    " ${TATIN_EXTERNAL_INC}
	-@echo "External object files:"
	-@echo "    " ${TATIN_EXTERNAL_OBJ}
	-@echo "External libraries:"
	-@echo "    " ${TATIN_EXTERNAL_LIB}

	@for m in ${TATIN_INTERFACE_DIR} ; do \
		${MAKE} -C $${m}/; \
	done


addtolib:
	@for m in ${TATIN_INTERFACE_DIR} ; do \
		${MAKE} addtolib -C $${m}/; \
	done
	@if [ -z "${TATIN_EXTERNAL_OBJ}" ]; then \
		echo "    + No object files from external packages to add "; \
	else \
		echo "    + Adding object files from all external packages " ${TATIN_EXTERNAL_OBJ}; \
		ar rcs ${TATIN_DIR}/libptatin3d.a ${TATIN_EXTERNAL_OBJ}; \
	fi


clean_all:
	@for m in ${TATIN_INTERFACE_DIR} ; do \
		${MAKE} clean_all -C $${m}/; \
	done

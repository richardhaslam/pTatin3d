## ==============================================================================
##
##     makefile for pTatin3d
##
## ==============================================================================

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

export TATIN_DIR	=	${PWD}

export TATIN_SRC = \
			dmda_update_coords.c \
			dmda_project_coords.c \
			dmda_compare.c \
			dmda_redundant.c \
			dmda_remesh.c \
			dmda_duplicate.c \
			dmda_view_petscvtk.c \
			dmda_checkpoint.c \
			dmda_bcs.c \
			mesh_deformation.c \
			mesh_quality_metrics.c \
			mesh_update.c \
			dmda_element_q2p1.c \
			element_type_Q2.c \
			swarm_fields.c \
			data_exchanger.c \
			MPntStd_def.c \
			MPntPStokes_def.c \
			QPntVolCoefStokes_def.c \
			output_paraview.c \
			material_point_utils.c \
			material_point_std_utils.c \
			material_point_point_location.c \
			material_point_load.c \
			ptatin_utils.c \
			ptatin_log.c \
			ptatin3d_stokes.c \
			ptatin3d.c \
			ptatin_std_dirichlet_boundary_conditions.c \
			ptatin_models.c \
			rheology.c \
			material_constants.c \
			stokes_rheology_viscous.c \
			stokes_rheology_vp_std.c \
			stokes_form_function.c \
			stokes_operators_mf.c \
			stokes_operators.c \
			quadrature.c \
			phase_map.c \
			element_utils_q2.c \
			element_utils_q1.c \
			stokes_assembly.c \
			monitors.c \
			mp_advection.c \
			units.c

export TATIN_OBJ	=	$(TATIN_SRC:.c=.o)

export TATIN_F_SRC = \
			stokes_q2p1_mf_operators_def.f90

export TATIN_F_OBJ	=	$(TATIN_F_SRC:.f90=_f90.o)

# to de-activate fortran kernels
export TATIN_F_SRC =

TATIN_TESTS	=	\
			test_dmda_update_coords \
			test_dmda_redundant \
			test_dmda_remesh \
			test_dmda_duplicate \
			test_dmda_checkpoint \
			test_dmda_gmg \
			test_mat_assem \
			test_q2dmda_remesh \
			test_stokes_operators \
			test_stokes_operators_approximations \
			test_mp_advection \
			test_material_point_load

#			test_model

TATIN_MODELS		=
TATIN_DRIVERS		= \
			ptatin_driver_ic \
			ptatin_driver_nostokessolve \
			ptatin_driver_asmsolve \
			ptatin_driver_asmsolve_approx \
			ptatin_driver_asmsolve2 \
			ptatin_driver_basic


export TATIN_INC			=	-I. -I${PETSC_DIR}/include ${PETSC_CC_INCLUDES}
export EXTERNAL_LIBS	=	${PETSC_LIB}
export TATIN_LIB			=	-L${TATIN_DIR} -L${TATIN_DIR}/models -lptatin3d -lptatin3dmodels -lptatin3d ${EXTERNAL_LIBS}
#export TATIN_LIB			=	${EXTERNAL_LIBS} -L${TATIN_DIR} -lptatin3d -L${TATIN_DIR}/models -lptatin3dmodels

## or run >make CFLAGS='-g -O0'
export TATIN_CFLAGS_DEBUG		= -O0 -g -Wall -Wno-unused-variable

export TATIN_CFLAGS_O2_OSX	= -std=gnu99 -O2 -Wall -Wno-unused-variable -Wstrict-aliasing -fstrict-aliasing -funroll-loops
export TATIN_CFLAGS_O3_OSX	= -std=gnu99 -O3 -g -fast -ftree-vectorize -ftree-vectorizer-verbose=2  -Wall -Wstrict-aliasing -Wno-unused-variable -fstrict-aliasing -funroll-loops

export TATIN_CFLGAS_O2_BGP      = -O2 -qarch=auto -qtune=auto -qcache=auto -qhot -qipa
export TATIN_CFLGAS_O3_BGP	= -O3 -qarch=auto -qtune=auto -qcache=auto -qhot -qipa

export TATIN_CFLAGS_O3_OCTO	= -O3 -std=gnu99 -mtune=native -march=native -funroll-loops -ffast-math

export TATIN_FFLAGS_DEBUG 	= -O0 -ffree-line-length-4000
export TATIN_FFLAGS_O2_OSX 	= -O2 -ffree-line-length-4000
export TATIN_FFLAGS_O3_BRUTUS 	= -O2 -ffree-line-length-4000 -mtune=native -march=native -funroll-loops -ffast-math -Wstrict-aliasing -fstrict-aliasing

export TATIN_CFLAGS		= ${TATIN_CFLAGS_O2_OSX}
export TATIN_CC			= ${PCC}
export TATIN_FFLAGS		= ${TATIN_FFLAGS_O3_BRUTUS}
export TATIN_FC			= ${FC}
#export TATIN_FC				= gfortran -m64


all: 
	-@echo "<<<<<<<<==================[[ pTatin-3d ]]==================>>>>>>>>"
	-@echo ${PWD}
	

	-@${MAKE} base_objects
	-@${MAKE} base_objects_fortran
	-@${MAKE} material_constants
	-@${MAKE} libs
	-@${MAKE} models
	-@${MAKE} tests
	-@${MAKE} drivers


base_objects:
	-@echo ""
	-@echo "*** Creating ptatin object files ***"
	@for a in ${TATIN_OBJ} ; do \
		${MAKE} $$a; \
	done

base_objects_fortran:
	-@echo ""
	-@echo "*** Creating ptatin fortran object files ***"
	@for a in ${TATIN_F_OBJ} ; do \
		${MAKE} $$a; \
	done


material_constants: base_objects
	-@${MAKE} -C material_constants

libs: base_objects
	-@echo ""
	-@echo "*** Creating libptatin3d ***"
	-@rm -f libptatin3d.a
	ar rcs libptatin3d.a ${TATIN_OBJ} ${TATIN_F_OBJ} material_constants/*.o


tests: libs models
	-@echo ""
	-@echo "*** Building tests ***"
	@for t in ${TATIN_TESTS} ; do \
		${MAKE} $${t}.o; \
		${PCC} ${TATIN_CFLAGS} -o $${t}.app $${t}.o ${TATIN_INC} ${TATIN_LIB}; \
	done

models: base_objects libs
	-@echo ""
	-@echo "*** Building pTatin models ***"
	-@${MAKE} -C models/


drivers: libs models
	-@echo ""
	-@echo "*** Building pTatin executables ***"
	@for t in ${TATIN_DRIVERS} ; do \
		${MAKE} $${t}.o; \
		${PCC} ${TATIN_CFLAGS} -o $${t}.app $${t}.o ${TATIN_INC} ${TATIN_LIB}; \
	done


%.o: %.c
	-@echo ""
	${PCC} ${TATIN_CFLAGS} -c $*.c ${TATIN_INC}

%_f90.o: %.f90
	-@echo ""
	${TATIN_FC} ${TATIN_FFLAGS} -c $*.f90 ${TATIN_INC} -o $*_f90.o


clean_all:
	-@echo ""
	-@echo "*** Removing material_constants ***"
	-@${MAKE} clean_all -C material_constants/
	-@echo "*** Removing models ***"
	-@${MAKE} clean_all -C models/
	-@echo "*** Removing base_objects ***"
	-@rm -f ${TATIN_OBJ}
	-@rm -f ${TATIN_F_OBJ}
	-@echo "*** Removing tests object files ***"
	-@rm -f *.o
	-@echo "*** Removing ptatin applications ***"
	-@rm -f *.app
	-@echo "*** Removing ptatin3d.a ***"
	-@rm -f libptatin3d.a


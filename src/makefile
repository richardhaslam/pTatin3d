## ==============================================================================
##
##     makefile for pTatin3d
##
## ==============================================================================

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

# Compilation options are to be placed in makefile.arch
include makefile.arch

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
			dmda_iterator.c \
			mesh_deformation.c \
			mesh_quality_metrics.c \
			mesh_update.c \
			dmdae.c \
			dmda_element_q2p1.c \
			dmda_element_q1.c \
			element_type_Q2.c \
			dmda_element_q1macrop1.c \
			swarm_fields.c \
			data_exchanger.c \
			sub_comm.c \
			MPntStd_def.c \
			MPntPStokes_def.c \
			MPntPStokesPl_def.c \
			QPntVolCoefStokes_def.c \
			QPntSurfCoefStokes_def.c \
			MPntPEnergy_def.c \
			QPntVolCoefEnergy_def.c \
			output_paraview.c \
			output_material_points.c \
			material_point_utils.c \
			material_point_std_utils.c \
			material_point_point_location.c \
			material_point_popcontrol.c \
			material_point_load.c \
			ptatin_utils.c \
			ptatin_init.c \
			ptatin_log.c \
			ptatin3d_stokes.c \
			stokes_output.c \
			ptatin3d_energy.c \
			phys_comp_energy.c \
			energy_assembly.c \
			energy_output.c \
			ptatin3d_stokes_q1macrop1.c \
			ptatin3d.c \
			ptatin_std_dirichlet_boundary_conditions.c \
			ptatin_models.c \
			rheology.c \
			material_constants.c \
			stokes_rheology_viscous.c \
			stokes_rheology_vp_std.c \
			stokes_rheology_lava.c \
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
			units.c \
			model_utils.c \
			eigen_operators.c \
			ksp_chebyrn.c \
			pc_semiredundant.c \
			pc_wsmp.c \
			geometry_object.c \
			geometry_object_evaluator.c \
			spm_utils.c \

export TATIN_OBJ	=	$(TATIN_SRC:.c=.o)

export TATIN_F_SRC = \
			stokes_q2p1_mf_operators_def.f90
# to de-activate fortran kernels
export TATIN_F_SRC =

export TATIN_F_OBJ	=	$(TATIN_F_SRC:.f90=_f90.o)

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
			test_material_point_load \
			test_stokes_q1macrop1 \
			ptatin_read_matrix \
			test_petsc_wsmp

#			test_model

TATIN_MODELS		=
TATIN_DRIVERS		= \
			ptatin_driver_ic \
			ptatin_driver_nostokessolve \
			ptatin_driver_asmsolve \
			ptatin_driver_asmsolve_approx \
			ptatin_driver_asmsolve2 \
			ptatin_driver_basic \
			ptatin_driver_ic_spm \
			ptatin_write_pvts \
			ptatin_driver_energy \
			ptatin_driver_linear_ts \
			ptatin_driver_nonlinear_ts

export LIBZ_LIB = -lz
export LIBZ_INC = -I/usr/include

export TATIN_INC	= -I. -I${PETSC_DIR}/include ${PETSC_CC_INCLUDES} ${LIBZ_INC}
export EXTERNAL_LIBS	= ${PETSC_LIB} ${LIBZ_LIB}
export TATIN_LIB	= -L${TATIN_DIR} -L${TATIN_DIR}/models -lptatin3d -lptatin3dmodels -lptatin3d ${EXTERNAL_LIBS}
#export TATIN_LIB	= ${EXTERNAL_LIBS} -L${TATIN_DIR} -lptatin3d -L${TATIN_DIR}/models -lptatin3dmodels

# Compilers are determined via petsc
export TATIN_CC		= ${PCC}
export TATIN_FC		= ${FC}
#export TATIN_FC	= gfortran -m64

# external packages
ifeq (${HAVE_LIBZ},1)

endif

ifeq (${HAVE_WSMP_SEQ},1)
	TATIN_CFLAGS += -DHAVE_WSSMP
	TATIN_LIB += ${WSMP_SEQ_LIB}
endif

ifeq (${HAVE_WSMP_P},1)
	TATIN_CFLAGS += -DHAVE_PWSSMP
	TATIN_LIB += ${WSMP_P_LIB}
endif


all: 
	-@echo "<<<<<<<<==================[[ pTatin-3d ]]==================>>>>>>>>"
	-@echo "Working directory:"
	-@echo ${PWD}

	
	-@${MAKE} base_objects
	-@${MAKE} base_objects_fortran
	-@${MAKE} material_constants
	-@${MAKE} interfaces
	-@${MAKE} libs
	-@${MAKE} models
	-@${MAKE} tests
	-@${MAKE} drivers

allsvn: 
	-@echo "<<<<<<<<==================[[ pTatin-3d ]]==================>>>>>>>>"
	-@echo ${PWD}
	
	-@${MAKE} svninfo
	-@${MAKE} base_objects
	-@${MAKE} base_objects_fortran
	-@${MAKE} material_constants
	-@${MAKE} interfaces
	-@${MAKE} libs
	-@${MAKE} models
	-@${MAKE} tests
	-@${MAKE} drivers


svninfo:
	-@perl ../tools/get-svninfo.pl

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


interfaces: base_objects
	-@echo ""
	-@echo "*** Creating ptatin interfaces ***"
	-@${MAKE} -C externalpackages/


libs: base_objects material_constants interfaces
	-@echo ""
	-@echo "*** Creating libptatin3d ***"
	-@rm -f libptatin3d.a
	ar rcs libptatin3d.a ${TATIN_OBJ} ${TATIN_F_OBJ} material_constants/*.o
	-@${MAKE} addtolib -C externalpackages/


tests: libs models
	-@echo ""
	-@echo "*** Building tests ***"
	@for t in ${TATIN_TESTS} ; do \
		${MAKE} $${t}.o; \
		rm -f $${t}.app; \
		${MAKE} $${t}.app; \
	done

models: base_objects libs
	-@echo ""
	-@echo "*** Building pTatin models ***"
	-@${MAKE} -C models/


drivers: libs models
	-@echo ""
	-@echo "*** Building pTatin drivers ***"
	@for t in ${TATIN_DRIVERS} ; do \
		${MAKE} $${t}.o; \
		rm -f $${t}.app; \
		${MAKE} $${t}.app; \
	done

ptatin_init.o: ptatin_init.c
	-@echo ""
	${PCC} ${TATIN_CFLAGS} -c ptatin_init.c ${TATIN_INC} -DCOMPFLAGS='${TATIN_CFLAGS}'



%.app: %.o ${TATIN_OBJ} material_constants/*.o
	${PCC} ${TATIN_CFLAGS} -o $*.app $*.o ${TATIN_INC} ${TATIN_LIB}; 

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
	-@echo "*** Removing interfaces ***"
	-@${MAKE} clean_all -C externalpackages/
	-@echo "*** Removing base_objects ***"
	-@rm -f ${TATIN_OBJ}
	-@rm -f ${TATIN_F_OBJ}
	-@echo "*** Removing tests object files ***"
	-@rm -f *.o
	-@echo "*** Removing ptatin applications ***"
	-@rm -f *.app
	-@echo "*** Removing ptatin3d.a ***"
	-@rm -f libptatin3d.a


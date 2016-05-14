ptatin-tests-y.c += $(call thisdir, \
			test_dmda_update_coords.c \
			test_dmda_redundant.c \
			test_dmda_remesh.c \
			test_dmda_duplicate.c \
			test_dmda_checkpoint.c \
			test_dmda_gmg.c \
			test_mat_assem.c \
			test_q2dmda_remesh.c \
			test_stokes_operators.c \
			test_stokes_operators_approximations.c \
			test_mp_advection.c \
			test_material_point_load.c \
            test_swarms.c \
			test_swarms_exchanger.c \
			test_stokes_q1macrop1.c \
			ptatin_read_matrix.c \
			test_petsc_wsmp.c \
            test_cjson.c \
			ptatin_driver_nostokessolve.c \
			ptatin_driver_asmsolve.c \
			ptatin_driver_asmsolve_approx.c \
			ptatin_driver_asmsolve2.c \
			ptatin_driver_basic.c \
			ptatin_driver_ic_spm.c \
            test_databucket_checkpoint.c \
            test_ptatin_quantity.c \
    )

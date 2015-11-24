libptatin3d-y.c += $(call thisdir, \
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
			data_bucket.c \
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
			xdmf_writer.c \
			output_material_points.c \
            output_material_points_p0.c \
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
			energy_coefficients.c \
			ptatin3d_stokes_q1macrop1.c \
			ptatin3d.c \
			ptatin_std_dirichlet_boundary_conditions.c \
			ptatin_models.c \
			rheology.c \
			material_constants.c \
                        material_constants_energy.c \
			stokes_rheology_viscous.c \
			stokes_rheology_vp_std.c \
			stokes_rheology_lava.c \
			stokes_form_function.c \
			stokes_operators_mf.c \
			stokes_operators.c \
			stokes_operators_tensor.c \
			quadrature.c \
			phase_map.c \
			cartgrid.c \
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
			pc_dmdarepart.c \
			geometry_object.c \
			geometry_object_evaluator.c \
			cJSON.c \
			geometry_object_parse.c \
			spm_utils.c \
			petsc_utils.c \
	)

libptatin3d-$(CONFIG_AVX).c += $(call thisdir, \
			stokes_operators_avx.c \
	)
# One file needs special flags.
$(OBJDIR)/$(call thisdir,ptatin_init.o) : CFLAGS += -DCOMPFLAGS='$(TATIN_CFLAGS)'

libptatin3d-$(CONFIG_FORTRAN).f += $(call thisdir, \
			stokes_q2p1_mf_operators_def.f90 \
	)

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
	)

ptatin-drivers-y.c += $(call thisdir, \
			ptatin_driver_ic.c \
			ptatin_write_pvts.c \
			ptatin_driver_energy.c \
			ptatin_driver_linear_ts.c \
			ptatin_driver_nonlinear_ts.c \
	)

TATIN_INC += -I$(abspath $(call thisdir,.))

include $(call incsubdirs,externalpackages material_constants models)

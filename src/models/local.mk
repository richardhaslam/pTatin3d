libptatin3dmodels-y.c += $(call thisdir, \
			ptatin_models_reg.c \
			template/model_ops_template.c \
			viscous_sinker/model_viscous_sinker_ops.c \
			gene3d/model_ops_gene3d.c \
			gene3d_nueve/model_ops_gene3d_nueve.c \
			indentor/model_indentor_ops.c \
			rift3D/model_rift3D_ops.c \
			rift3D_T/model_rift3D_T_ops.c \
			sierra/model_ops_Sierra.c \
			folding/model_folding_ops.c \
			folding2d/model_ops_folding2d.c \
			basin_comp/model_ops_basin_comp.c \
			fault_fold/model_ops_fault_fold.c \
			wrench_fold/model_ops_wrench_fold.c \
			fault_fold_plastic/model_ops_fault_fold_plastic.c \
			advdiff_example/model_ops_advdiff_example.c \
			delamination/model_delamination_ops.c \
			riftrh/model_ops_riftrh.c \
			geomod2008/model_geomod2008.c \
			multilayer_folding/model_ops_multilayer_folding.c \
			multilayer_folding/compute_analytic_growth_rate_3d.c \
			submarinelavaflow/submarinelavaflow_ops.c \
			submarinelavaflow/profiles/plot_profiles.c \
			ex_subduction/ex_subduction_ops.c \
			iplus_slab_plume/iplus_ops.c \
			iplus_slab_plume/iplus_plume_description.c \
			iplus_slab_plume/iplus_slab_description.c \
			subduction_initiation2d/model_subduction_initiation2d_ops.c \
			convection2d/model_convection2d_ops.c \
			thermal_shearband/thermal_sb_model_definition.c \
			slab_detachment3d/sd3d_definition.c \
			pinch_and_swell/pinchandswell_def.c \
			plasticdemo/demo.c \
mms/mms_def.c \
	)

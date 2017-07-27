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
			advdiff_example/model_ops_advdiff_example.c \
			delamination/model_delamination_ops.c \
			rift_oblique3d/model_ops_rift_oblique3d.c \
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
                        static_box/static_box.c \
                        static_box_thermomech/static_box_tm.c \
                        analytics_vv/SolKxSolution.c analytics_vv/SolCxSolution.c analytics_vv/analytics_vv.c \
	)

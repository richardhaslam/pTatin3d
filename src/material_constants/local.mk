libptatin3d-y.c += $(call thisdir, \
			MaterialConst_MaterialType_def.c \
			MaterialConst_DensityConst_def.c \
			MaterialConst_DensityBoussinesq_def.c \
			MaterialConst_ViscosityConst_def.c \
			MaterialConst_ViscosityZ_def.c \
			MaterialConst_ViscosityArrh_def.c \
			MaterialConst_ViscosityFK_def.c \
			MaterialConst_SoftLin_def.c \
			MaterialConst_SoftExpo_def.c \
			MaterialConst_PlasticMises_def.c \
			MaterialConst_PlasticDP_def.c \
			EnergyMaterialConstants_def.c \
                        EnergySourceConst_def.c \
			EnergySourceDecay_def.c \
			EnergySourceAdiabaticAdvection_def.c \
	)

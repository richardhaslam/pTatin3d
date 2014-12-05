
# Initialize external package cflag and include path variables
package-cflag-y :=
package-inc-y :=


# ----------------------------------------------------------------------------
# SPMA:
# Simple C landscape evolution model. 
# Serves as a template for coupling surface process models. Contained within ptatin3d.
#
SPMA_DIR                     = $(PWD)/src/externalpackages/spmA
package-cflag-$(CONFIG_SPMA) += -DPTATIN_HAVE_SPMA
package-inc-$(CONFIG_SPMA)   += -I$(SPMA_DIR)

# Register files for SPMA
ptatin-externals-$(CONFIG_SPMA).o += $(call thisdir, \
	spmA/spmA.o \
        )


# ----------------------------------------------------------------------------
# SPM FastScapeV3:
# Landscape evolution model written by Jean Braun.
#

SPM_FSCAPE_DIR                    = $(PWD)/src/externalpackages/FastScape_V3_lib
package-cflag-$(CONFIG_FASTSCAPE) += -DPTATIN_HAVE_FASTSCAPE_V3
package-inc-$(CONFIG_FASTSCAPE)   += -I$(SPM_FSCAPE_DIR)

# Register files for FASTSCAPE
SPM_FSCAPE_DIR = FastScape_V3_lib
ptatin-externals-$(CONFIG_FASTSCAPE).o += $(call thisdir, \
	$(SPM_FSCAPE_DIR)/Concavity_f90.o \
	$(SPM_FSCAPE_DIR)/Curvature_f90.o \
	$(SPM_FSCAPE_DIR)/DEM_f90.o \
	$(SPM_FSCAPE_DIR)/Demoulin_f90.o \
	$(SPM_FSCAPE_DIR)/libFastScape_f90.o \
	$(SPM_FSCAPE_DIR)/Slope_f90.o \
	$(SPM_FSCAPE_DIR)/SteepnessIndex_f90.o \
	$(SPM_FSCAPE_DIR)/bmp_f90.o \
	$(SPM_FSCAPE_DIR)/compute_sediment_flux_f90.o \
	$(SPM_FSCAPE_DIR)/diffusion_f90.o \
	$(SPM_FSCAPE_DIR)/find_donors_f90.o \
	$(SPM_FSCAPE_DIR)/find_receiver_f90.o \
	$(SPM_FSCAPE_DIR)/flexure_f90.o \
	$(SPM_FSCAPE_DIR)/fluvial_erosion_f90.o \
	$(SPM_FSCAPE_DIR)/four1_f77.o \
	$(SPM_FSCAPE_DIR)/initial_topography_f90.o \
	$(SPM_FSCAPE_DIR)/initialize_f90.o \
	$(SPM_FSCAPE_DIR)/interpol_routines_f90.o \
	$(SPM_FSCAPE_DIR)/interpolate_f90.o \
	$(SPM_FSCAPE_DIR)/local_minima_f90.o \
	$(SPM_FSCAPE_DIR)/metric_f90.o \
	$(SPM_FSCAPE_DIR)/module_definitions_f90.o \
	$(SPM_FSCAPE_DIR)/orography_f90.o \
	$(SPM_FSCAPE_DIR)/precipitation_f90.o \
	$(SPM_FSCAPE_DIR)/realft_f77.o \
	$(SPM_FSCAPE_DIR)/scanfile_f90.o \
	$(SPM_FSCAPE_DIR)/set_bc_f90.o \
	$(SPM_FSCAPE_DIR)/sinft_f77.o \
	$(SPM_FSCAPE_DIR)/tridag_f90.o \
	$(SPM_FSCAPE_DIR)/uplift_f90.o \
	)

# Bundle all cflags and include paths from packages together
TATIN_CFLAGS += $(package-cflag-y)
TATIN_INC    += $(package-inc-y)


include $(call incsubdirs,interfaces)

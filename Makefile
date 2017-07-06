## ==============================================================================
##
##     makefile for pTatin3d
##
## ==============================================================================
.SECONDEXPANSION:           # To expand $$(@D)/.DIR
.SUFFIXES:                  # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:           # Delete likely-corrupt target file if rule fails

include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables

all : info tests drivers

# Compilation options are to be placed in makefile.arch;
# if that (untracked) file does not exist, defaults are copied there
makefile.arch:
	-@echo "[pTatin config] using config/makefile.arch.default as makefile.arch"
	-@cp config/makefile.arch.default $@
-include makefile.arch

OBJDIR ?= $(PETSC_ARCH)/obj
LIBDIR ?= $(PETSC_ARCH)/lib
BINDIR ?= $(PETSC_ARCH)/bin

CONFIG_SPMA      ?= n
CONFIG_FASTSCAPE ?= n

CONFIG_FORTRAN = y
CONFIG_AVX ?= n

# Directory that contains most recently-parsed makefile (current)
thisdir = $(addprefix $(dir $(lastword $(MAKEFILE_LIST))),$(1))
incsubdirs = $(addsuffix /local.mk,$(call thisdir,$(1)))

libptatin3d-y.c :=
libptatin3d-y.f :=
libptatin3dmodels-y.c :=
libptatin3dmodels-y.f :=
ptatin-tests-y.c :=
ptatin-drivers-y.c :=
ptatin-externals-y.o :=
TATIN_INC := $(PETSC_CC_INCLUDES) -I${PWD}/include

# Recursively include files for all targets
include local.mk

#### Rules ####
ifeq ($(V),)
  quiet_HELP := "Use \"$(MAKE) V=1\" to see the verbose compile lines.\n"
  quiet = @printf $(quiet_HELP)$(eval quiet_HELP:=)"  %10s %s\n" "$1$2" "$@"; $($1)
else ifeq ($(V),0) # Same, but do not print any help
  quiet = @printf "  %10s %s\n" "$1$2" "$@"; $($1)
else # Show the full command line
  quiet = $($1)
endif

.PHONY: libptatin3d libptatin3dmodels externals tests drivers

libptatin3d = $(LIBDIR)/libptatin3d.$(AR_LIB_SUFFIX)
libptatin3d : $(libptatin3d)

$(libptatin3d) : $(libptatin3d-y.c:%.c=$(OBJDIR)/%.o) $(libptatin3d-y.f:%.f90=$(OBJDIR)/%.o) $(ptatin-externals-y.o)

libptatin3dmodels = $(LIBDIR)/libptatin3dmodels.$(AR_LIB_SUFFIX)
libptatin3dmodels : $(libptatin3dmodels)
$(libptatin3dmodels) : $(libptatin3dmodels-y.c:%.c=$(OBJDIR)/%.o) $(libptatin3dmodels-y.f:%.f90=$(OBJDIR)/%.o)

info:
	-@echo "—————————————————————————————————————————————————————————————————"
	-@echo " pTatin3D"
	-@echo "    git url: https://bitbucket.org/jedbrown/ptatin3d.git"
	-@echo "    update version information by running \"make releaseinfo\""
	-@echo "—————————————————————————————————————————————————————————————————"

releaseinfo:
	-@echo "—————————————————————————————————————————————————————————————————————————————————————————————————————————————"
	-@echo "  Updating the contents of src/ptatin_version_info.h from: https://bitbucket.org/jedbrown/ptatin3d.git"
	-@echo "  pTatin3D version info appears in: stdout; ptat3d log; generated PETSc summary/options files"
	-@echo "—————————————————————————————————————————————————————————————————————————————————————————————————————————————"
	$(shell python utils/mgmt/ptatin3d-generate-git-info.py)
	$(shell mv ptatin_version_info.h src/ptatin_version_info.h)
	$(shell git update-index --assume-unchanged src/ptatin_version_info.h)

externals:
	-@echo ——————— EXTERNAL PACKAGE OBJECT FILES ———————
	-@echo $(ptatin-externals-y.o)
	-@echo ——————— EXTERNAL PACKAGE OBJECT INCLUDE PATH ———————
	-@echo $(TATIN_INC)
	-@echo ——————— EXTERNAL PACKAGE OBJECT CFLAGS ———————
	-@echo $(TATIN_CFLAGS)

test :
	cd tests && ./runTests.py -s -t $(shell ./tests/getTestGroup.py smoke)

testcheck :
	cd tests && ./runTests.py -s -t $(shell ./tests/getTestGroup.py smoke) --verify

testall :
	cd tests && ./runTests.py -s  

testallcheck :
	cd tests && ./runTests.py -s --verify

testallclean :
	cd tests && ./runTests.py -s -p

testclassic :
	cd tests && ./runTests.py -s -t $(shell ./tests/getTestGroup.py classic)

testclassiccheck :
	cd tests && ./runTests.py -s -t $(shell ./tests/getTestGroup.py classic) --verify

%.$(AR_LIB_SUFFIX) : | $$(@D)/.DIR
	$(call quiet,AR) $(AR_FLAGS) $@ $^
	$(call quiet,RANLIB) $@

ifeq ($(PETSC_LANGUAGE),CXXONLY)
  cc_name := CXX
else
  cc_name := CC
endif

# gcc/gfortran style dependency flags; these are set in petscvariables starting with petsc-3.5
C_DEPFLAGS ?= -MMD -MP
FC_DEPFLAGS ?= -MMD -MP

TATIN_COMPILE.c = $(call quiet,$(cc_name)) -c $(PCC_FLAGS) $(CCPPFLAGS) $(TATIN_CFLAGS) $(TATIN_INC) $(CFLAGS) $(C_DEPFLAGS)
TATIN_COMPILE.f90 = $(call quiet,FC) -c $(FC_FLAGS) $(FCPPFLAGS) $(TATIN_INC) $(FFLAGS) $(FC_DEPFLAGS)


# Test executables
tests: $(ptatin-tests-y.c:%.c=$(BINDIR)/%.app)
$(ptatin-tests-y.c:%.c=$(BINDIR)/%.app) : $(libptatin3dmodels) $(libptatin3d)
.SECONDARY: $(ptatin-tests-y.c:%.c=$(OBJDIR)/%.o) # don't delete the intermediate files


# Drivers
drivers: $(ptatin-drivers-y.c:%.c=$(BINDIR)/%.app)
$(ptatin-drivers-y.c:%.c=$(BINDIR)/%.app) : $(libptatin3dmodels) $(libptatin3d)
.SECONDARY: $(ptatin-drivers-y.c:%.c=$(OBJDIR)/%.o)

$(BINDIR)/%.app : $(OBJDIR)/%.o | $$(@D)/.DIR
	$(call quiet,PCC_LINKER) $(TATIN_CFLAGS) -o $@ $^ $(PETSC_SNES_LIB) $(LIBZ_LIB)
#@mv $@ $(BINDIR)
	@ln -sF $(abspath $@) $(BINDIR)

$(OBJDIR)/%.o: %.c | $$(@D)/.DIR
	$(TATIN_COMPILE.c) $(abspath $<) -o $@

$(OBJDIR)/%.o: %.f90 | $$(@D)/.DIR
	$(TATIN_COMPILE.f90) $(abspath $<) -o $@

%/.DIR :
	@mkdir -p $(@D)
	@touch $@

.PRECIOUS: %/.DIR
.SUFFIXES:                 # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:          # Delete likely-corrupt target file if rule fails

.PHONY: clean all print

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(BINDIR)

# make print VAR=the-variable
print:
	@echo $($(VAR))

srcs.c := $(libptatin3d-y.c) $(libptatin3dmodels-y.c) $(ptatin-tests-y.c) $(ptatin-drivers-y.c)
# Bundle in objects from external packages
%srcs.o := $(srcs.c:%.c=$(OBJDIR)/%.o)
srcs.o := $(srcs.c:%.c=$(OBJDIR)/%.o) $(ptatin-externals-y.o)
srcs.d := $(srcs.o:%.o=%.d)
# Tell make that srcs.d are all up to date.  Without this, the include
# below has quadratic complexity, taking more than one second for a
# do-nothing build of PETSc (much worse for larger projects)
$(srcs.d) : ;

-include $(srcs.d)

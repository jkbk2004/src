# SCCS IDENTIFICATION:   $HeadURL$
# SCCS IDENTIFICATION:   @(#)$Id$

# ###############
# Local Variables
# ###############

HR     := ========================================
HR     := $(HR)$(HR)
COMMA  := ,
MKFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
DIR    := $(dir $(MKFILE))

# ##########################
# Check Environment Settings
# ##########################
#ARCH=intelsse-impi-sm-relo
#TYPE=mpi
#CPP_EXTRAS=-DEOS_SIG2 -DEOS_17T -DESPC_COUPLE

ifndef ARCH
$(error ARCH is not set.)
endif

ifndef TYPE
$(error TYPE is not set.)
endif

ifndef CPP_EXTRAS
$(error CPP_EXTRAS is not set.)
endif

# #################
# Inherit CPP Flags
# #################

override ESMF_F90COMPILECPPFLAGS += $(CPP_EXTRAS)
override ESMF_CXXCOMPILECPPFLAGS += $(CPP_EXTRAS)

# ###########################
# Include ESMFMKFILE fragment
# ###########################

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif
include $(ESMFMKFILE)

# #####################################
# Fused Multiply-Add (FMA) Instructions
# #####################################

MACHINE=$(shell hostname)
ifeq ($(findstring conrad,$(MACHINE)),conrad)
ESMF_F90COMPILEOPTS += -no-fma
else ifeq ($(findstring gordon,$(MACHINE)),gordon)
ESMF_F90COMPILEOPTS += -no-fma
endif

# #################################
# Compile with Debugging Directives
# #################################

ifeq ($(DEBUG),on)
override ESMF_F90COMPILECPPFLAGS += -DDEBUG
override ESMF_CXXCOMPILECPPFLAGS += -DDEBUG
endif

ifdef DBGBUILD
ESMF_F90COMPILEOPTS += -g -traceback
ESMF_CXXCOMPILEOPTS += -g -traceback
endif

# ###########################
# Determine Installation Path
# ###########################

ifndef DESTDIR
DESTDIR  := $(DIR)
endif

ifndef INSTDIR
INSTDIR  := HYCOM_$(shell date '+%Y-%m-%d-%H-%M-%S')
endif

ifndef INSTPATH
INSTPATH := $(abspath $(DESTDIR)/$(INSTDIR))
endif

# ###############
# Model Variables
# ###############
ifndef HYCOM_DIR
MODEL_RDIR := $(abspath $(DIR)/..)
else
MODEL_RDIR := $(abspath $(HYCOM_DIR))
endif
MODEL_DIR    := $(abspath $(MODEL_RDIR))
MODEL_OBJDIR := $(abspath $(MODEL_DIR))
MODEL_MODDIR := $(abspath $(MODEL_DIR))
MODEL_LIBDIR := $(abspath $(MODEL_DIR))
MODEL_EXEDIR := $(abspath $(MODEL_DIR))
ifdef USE_OCN
MODEL_LIB    := $(abspath $(MODEL_LIBDIR)/ocn_esmf.a)
else
MODEL_LIB    := $(abspath $(MODEL_LIBDIR)/hycom_esmf.a)
endif
MODEL_EXE    := $(abspath $(MODEL_EXEDIR)/hycom)
MODEL_CFGDIR := $(abspath $(MODEL_DIR)/config)
MODEL_CFGINC := $(abspath $(MODEL_CFGDIR)/$(ARCH)_$(TYPE))

MODEL_MODS  := $(abspath $(MODEL_MODDIR)/mod_xc.mod)
MODEL_MODS  += $(abspath $(MODEL_MODDIR)/mod_cb_arrays.mod)
MODEL_MODS  += $(abspath $(MODEL_MODDIR)/mod_hycom.mod)
MODEL_MODS  += $(abspath $(MODEL_MODDIR)/mod_archiv.mod)

#MODEL_MODS   += $(abspath $(MODEL_MODDIR)/read_impexp_config_mod.mod) # can't find module
#MODEL_MODS   += $(abspath $(MODEL_MODDIR)/impexpField_cdf_mod.mod)    # can't find module

MODEL_FILES  := $(MODEL_LIB) $(MODEL_MODS)

# ##################
# Model Dependencies
# ##################

# #############
# Cap Variables
# #############

CAP_DIR       := $(abspath $(DIR))
CAP_LIB       := libhycom_nuopc.a
CAP_MK        := hycom.mk
CAP_DEP_FRONT := HYCOM_Mod
CAP_VERS      := VERSION

CAP_OBJS      := HYCOM_OceanComp.o
CAP_OBJS      += HYCOM_ESMF_Extensions.o
CAP_OBJS      += hycom_couple.o
CAP_OBJS      += read_impexp_config_mod.o
CAP_OBJS      += impexpField_cdf_mod.o
CAP_OBJS      += export_from_hycom_tiled.o
CAP_OBJS      += hycom_read_latlon.o
CAP_OBJS      += hycom_nuopc_flags.o

ifdef ESPC_OCN
CAP_MODS      := ocean_mod.mod
else
CAP_MODS      := hycom_mod.mod
endif
CAP_MODS      += hycom_esmf_extensions.mod
CAP_MODS      += hycom_couple.mod
CAP_MODS      += read_impexp_config_mod.mod
CAP_MODS      += impexpfield_cdf_mod.mod
CAP_MODS      += hycom_read_latlon.mod
CAP_MODS      += hycom_nuopc_flags.mod

CAP_FILES     := $(CAP_OBJS) $(CAP_MODS) $(CAP_LIB) $(CAP_VERS) $(CAP_MK)

# ###############################
# Include Model Makefile Fragment
# ###############################

override ESMF_F90COMPILEPATHS += -I$(MODEL_MODDIR)
#include $(MODEL_MKINC)
include $(MODEL_CFGINC)
override ESMF_F90COMPILEOPTS     += $(FCFFLAGS)
override ESMF_CXXCOMPILEOPTS     += $(CCFLAGS)
override ESMF_F90COMPILECPPFLAGS += $(CPPFLAGS)
override ESMF_CXXCOMPILECPPFLAGS += $(CPPFLAGS)

#override DEP_SHRD_PATH          = $(patsubst -L%,%,$(LIS_LIB_PATHS))
#override DEP_SHRD_LIBS          = $(patsubst -l%,%,$(LIS_LIB_FLAGS))

# #######################
# Primary Makefile Target
# #######################
.PHONY: nuopc nuopcinstall nuopcdistclean nuopcclean install_mk

nuopc: $(CAP_FILES)

nuopcinstall: $(CAP_LIB) $(CAP_MODS) $(CAP_VERS) \
 $(addprefix $(INSTPATH)/,$(CAP_MODS)) \
 $(addprefix $(INSTPATH)/,$(CAP_LIB)) \
 $(addprefix $(INSTPATH)/,$(CAP_VERS)) \
 install_mk

# ############
# Dependencies
# ############

HYCOM_OceanComp.o: export_from_hycom_tiled.o \
 HYCOM_ESMF_Extensions.o \
 hycom_nuopc_flags.o \
 hycom_couple.o \
 read_impexp_config_mod.o \
 impexpField_cdf_mod.o
hycom_couple.o: export_from_hycom_tiled.o \
 hycom_read_latlon.o $(MODEL_MODS)
export_from_hycom_tiled.o: $(MODEL_MODS)
read_impexp_config_mod.o: hycom_nuopc_flags.o

ifdef ESPC_OCN
ocean_mod.mod: HYCOM_OceanComp.o
else
hycom_mod.mod: HYCOM_OceanComp.o
endif
hycom_esmf_extensions.mod: HYCOM_ESMF_Extensions.o
hycom_couple.mod: hycom_couple.o
read_impexp_config_mod.mod: read_impexp_config_mod.o
impexpfield_cdf_mod.mod: impexpField_cdf_mod.o
hycom_read_latlon.mod: hycom_read_latlon.o
hycom_nuopc_flags.mod: hycom_nuopc_flags.o

# ###########
# Build model
# ###########

build_info:
	@echo $(HR)
	@echo "Build Info..."
	@echo
	@echo "CPLFLAG: $(CPLFLAG) "
	@echo "MODEL_OBJDIR: $(MODEL_OBJDIR) "
	@echo "MODEL_OBJS: $(MODEL_OBJS) "

build_dirs:
	@if [ ! -d $(MODEL_OBJDIR) ]; then mkdir -p $(MODEL_OBJDIR); fi
	@if [ ! -d $(MODEL_MODDIR) ]; then mkdir -p $(MODEL_MODDIR); fi
	@if [ ! -d $(MODEL_LIBDIR) ]; then mkdir -p $(MODEL_LIBDIR); fi
	@if [ ! -d $(MODEL_EXEDIR) ]; then mkdir -p $(MODEL_EXEDIR); fi

build_model:
	@echo $(HR)
	@echo "Building Model..."
	@echo ""
	cd $(MODEL_DIR) && make esmf ARCH="$(ARCH)" TYPE="$(TYPE)" CPP_EXTRAS="$(CPP_EXTRAS)"

$(MODEL_MODS): build_model

$(MODEL_LIB): build_model
	@echo $(HR)
	@echo "Building Model Library..."
	@echo
	ar cr $@ $(MODEL_OBJDIR)/*.o

# ##############
# Build Settings
# ##############
.SUFFIXES:
.SUFFIXES: .c .C .f90 .F90 .F .f

.C:
	@echo "Must have an explicit rule for" $*
.c:
	@echo "Must have an explicit rule for" $*

%.o : %.f
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(CPPFLAGS) $(CPLFLAG) $(MODFLAG) $(INCFLAG) $< -o $@

%.o : %.F
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(CPPFLAGS) $(CPLFLAG) $(MODFLAG) $(INCFLAG) $< -o $@

%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(CPPFLAGS) $(CPLFLAG) $(MODFLAG) $(INCFLAG) $< -o $@

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(CPPFLAGS) $(CPLFLAG) $(MODFLAG) $(INCFLAG) $< -o $@

# #####################
# Build NUOPC Component
# #####################

$(CAP_LIB): $(MODEL_LIB) $(CAP_OBJS)
	@echo $(HR)
	@echo "Copying static library $@..."
	@echo
	$(call checkfile, $(MODEL_LIB))
	cp $(MODEL_LIB) $@
	ar cr $@ $(CAP_OBJS)

$(CAP_VERS): $(CAP_LIB) $(CAP_MODS)
	@echo $(HR)
	@echo "Generating Version Information"
	@echo
ifneq (,$(wildcard $(MODEL_RDIR)/.git))
	@echo "# NUOPC Cap Revision #" > $(CAP_VERS)
	@git log $(CAP_DIR) | grep -m 1 "commit " >> $(CAP_VERS)
	@git log $(CAP_DIR) | grep -m 1 "Author: " >> $(CAP_VERS)
	@git log $(CAP_DIR) | grep -m 1 "Date: " >> $(CAP_VERS)
	@echo >> $(CAP_VERS)
	@echo "# Model Revision     #" >> $(CAP_VERS)
	@git log $(MODEL_RDIR) | grep -m 1 "commit " >> $(CAP_VERS)
	@git log $(MODEL_RDIR) | grep -m 1 "Author: " >> $(CAP_VERS)
	@git log $(MODEL_RDIR) | grep -m 1 "Date: " >> $(CAP_VERS)
	@echo >> $(CAP_VERS)
else ifneq (,$(wildcard $(MODEL_RDIR)/.svn))
	@echo "# NUOPC Cap Revision #" > $(CAP_VERS)
	@svn info . | grep URL >> $(CAP_VERS)
	@svn info . | grep "Last Changed Rev" >> $(CAP_VERS)
	@echo >> $(CAP_VERS)
	@echo "# Model Revision     #" >> $(CAP_VERS)
	@svn info $(MODEL_RDIR) | grep URL >> $(CAP_VERS)
	@svn info $(MODEL_RDIR) | grep "Last Changed Rev" >> $(CAP_VERS)
	@echo >> $(CAP_VERS)
else
	@echo "# Version Information Not Found" > $(CAP_VERS)
endif

$(CAP_MK): $(CAP_LIB) $(CAP_MODS)
	@echo $(HR)
	@echo "Generating NUOPC Makefile Fragment"
	@echo
	@echo "# ESMF self-describing build dependency makefile fragment" > $(CAP_MK)
	@echo "" >> $(CAP_MK)
	@echo "ESMF_DEP_FRONT     = $(CAP_DEP_FRONT)" >> $(CAP_MK)
	@echo "ESMF_DEP_INCPATH   = $(CAP_DIR)" >> $(CAP_MK)
	@echo "ESMF_DEP_CMPL_OBJS = " >> $(CAP_MK)
	@echo "ESMF_DEP_LINK_OBJS = $(CAP_DIR)/$(CAP_LIB)" >> $(CAP_MK)
	@echo "ESMF_DEP_SHRD_PATH = $(DEP_SHRD_PATH)" >> $(CAP_MK)
	@echo "ESMF_DEP_SHRD_LIBS = $(DEP_SHRD_LIBS)" >> $(CAP_MK)

# -----------------------------------------------------------------------------
# Install Library, Modules, and Makefile Fragment
# -----------------------------------------------------------------------------

$(INSTPATH)/%:
	@echo $(HR)
	@echo "Installing $(notdir $@)"
	@echo
	@mkdir -p $(INSTPATH)
	@cp $(notdir $@) $@

install_mk:
	@echo $(HR)
	@echo "Installing NUOPC Makefile Fragment"
	@echo
	@mkdir -p $(INSTPATH)
	@echo "# ESMF self-describing build dependency makefile fragment" > $(INSTPATH)/$(CAP_MK)
	@echo "" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_FRONT     = $(CAP_DEP_FRONT)" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_INCPATH   = $(INSTPATH)" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_CMPL_OBJS = " >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_LINK_OBJS = $(INSTPATH)/$(CAP_LIB) $(DEP_LINK_OBJS)" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_SHRD_PATH = $(DEP_SHRD_PATH)" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_SHRD_LIBS = $(DEP_SHRD_LIBS)" >> $(INSTPATH)/$(CAP_MK)

# ###########
# Check Build
# ###########

define checkfile
	@if [ ! -e $(1) ]; then \
	echo "File is missing:$(1)"; \
	exit 1; fi;

endef # blank line in checkfile is required

define checkdir
	@if [ ! -d $(1) ]; then \
	echo "Directory is missing:$(1)"; \
	exit 1; fi;
endef # blank line in checkdir is required

check: check_esmf check_model check_cap

# ##################
# Check ESMF Version
# ##################

check_esmf:
	@echo $(HR)
	@echo "Checking ESMFMKFILE file..."
	@echo
	@echo "ESMFMKFILE=$(ESMFMKFILE)"
	@if [ "$(ESMF_VERSION_MAJOR)" -lt 7 ]; then \
	echo "Please use ESMF version 7+"; \
	exit 1; fi;
	@echo "ESMF Version=$(ESMF_VERSION_STRING)"

# ###########
# Check Model
# ###########

check_model:
	@echo $(HR)
	@echo "Checking for Model files..."
	@echo
	$(foreach FILENAME, $(MODEL_FILES), $(call checkfile, $(FILENAME)))

# #########
# Check Cap
# #########

check_cap:
	@echo $(HR)
	@echo "Checking for WRF-Hydro NUOPC files..."
	@echo
	$(foreach FILENAME, $(CAP_FILES), $(call checkfile, $(FILENAME)))

# -----------------------------------------------------------------------------
# Clean Cap and Model
# -----------------------------------------------------------------------------
.PHONY: nuopcdistclean modelclean nuopcclean

nuopcdistclean: modelclean nuopcclean

modelclean:
	@echo $(HR)
	@echo "Cleaning Model build..."
	@echo ""
#	cd $(MODEL_DIR) && make clean ARCH="$(ARCH)" TYPE="$(TYPE)" CPP_EXTRAS="$(CPP_EXTRAS)"
	rm -f $(MODEL_OBJDIR)/*.o $(MODEL_OBJDIR)/*.mod
	rm -f $(MODEL_LIB)

# ##########
# Clean  Cap
# ##########
nuopcclean:
	@echo $(HR)
	@echo "Cleaning Cap build..."
	@echo
	rm -f $(CAP_FILES)

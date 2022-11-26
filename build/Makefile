.SUFFIXES:
.SECONDARY:
.PHONY: all clean distclean help install uninstall force

MAKEFLAGS += -rR
SHELL     := /bin/bash

###############################################################################
# The SOFI top-level directory as checked out from git
# Example: /path/to/sofi/ (absolute) or ../ (relative to build directory)
# Note: $(BASEDIR)/src should contain the SOFI source code
###############################################################################
BASEDIR  := ..

###############################################################################
# Installation directory: 
# Example: /path/to/inst/bin (absolute) or ../bin (relative to build directory)
###############################################################################
INSTDIR := ../bin

###############################################################################
# Name of executable
###############################################################################
FDMOD := sofi2D
MERGE := snapmerge

###############################################################################
# Compiler / MPI settings
###############################################################################
### GNU compiler
CC             := mpicc
CPPFLAGS        = -D_GNU_SOURCE -DLOG_COLOR
EXTRA_CPPFLAGS := -DLOG_MPI
CFLAGS          = -O3 -Wall -Wextra -g -fstack-protector -std=c99
CDEPEND         = -MMD -MP -MF .dep-$(subst /,-,$@)
LD             := $(CC)
LDFLAGS         = -lm
EXTRA_LDFLAGS  := -lmpi

###############################################################################
### Intel compiler (OneAPI)
#CC             := mpicc -cc=icx
#CPPFLAGS        = -D_GNU_SOURCE -DLOG_COLOR
#EXTRA_CPPFLAGS := -DLOG_MPI
#CFLAGS          = -O3 -Wall -Wextra -g -fstack-protector -std=c99 --intel
#CDEPEND         = -MMD -MP -MF .dep-$(subst /,-,$@)
#LD             := $(CC)
#LDFLAGS         = -static-intel -lm 
#EXTRA_LDFLAGS  := -lmpi

###############################################################################
###############################################################################
# Usually, there is no need to change anything below this line 
###############################################################################
###############################################################################

SRCTREE   := $(strip $(realpath $(BASEDIR)/src))
FDMOD_O   := $(strip ./fdmod_obj)
MERGE_O   := $(strip ./merge_obj)
dummy     := $(shell mkdir -p $(FDMOD_O))
dummy     := $(shell mkdir -p $(MERGE_O))
FDMOD_EXE := $(strip $(FDMOD))
MERGE_EXE := $(strip $(MERGE))
INSTALL   := install -D -p
RM        := rm -f --preserve-root
DEPS      := $(wildcard .dep*)
COMMIT    := $(strip $(shell cd $(BASEDIR) && git log --format="%h %as" 2>/dev/null | head -n1))
ifeq ($(COMMIT),)
  COMMIT  := unknown
endif
CPPFLAGS  += -DLOG_COMMIT="\"$(COMMIT)\""

MERGE_FILES := snapmerge.c \
               holbergcoeff.c \
               json_parser.c \
               logging.c \
               merge.c \
               read_par_json.c \
               readdsk.c \
               writedsk.c \
               util.c 

FDMOD_FILES := sofi2D.c \
               abs_update.c \
               absorb.c \
               av_mue.c \
               av_rho.c \
               av_tau.c \
               catseis.c \
               check_fs.c \
               checkfd.c \
               cpml_update.c \
               debug_buffers.c \
	       exchange_par.c \
               exchange_s.c \
               exchange_v.c \
               hh_elastic.c \
               hh_elastic_TTI.c \
               hh_elastic_VTI.c \
               hh_visco.c \
               hh_visco_TTI.c \
               hh_visco_VTI.c \
               holbergcoeff.c \
               initproc.c \
               json_parser.c \
               logging.c \
               matcopy.c \
               matcopy_elastic.c \
               mergemod.c \
               outseis_glob.c \
               operators_s.c \
               operators_v.c \
               PML_pro.c \
               prepare_update_s.c \
               prepare_update_s_vti.c \
               prepare_update_s_tti.c \
               prepare_update_s_4.c \
               psource.c \
               rd_sour.c \
               read_checkpoint.c \
               read_par_json.c \
               readdsk.c \
               readmod_visco.c \
               readmod_elastic.c \
               readmod_elastic_vti.c \
               readmod_elastic_tti.c \
               readmod_visco_vti.c \
               readmod_visco_tti.c \
               receiver.c \
               save_checkpoint.c\
               saveseis_glob.c \
               seismo_ssg.c \
               snap_ssg.c \
               sources.c \
               splitrec.c \
               splitsrc.c \
               subgrid_bounds.c \
               surface.c \
               surface_elastic.c \
               update_s_elastic_abs.c \
               update_s_elastic_VTI_abs.c \
               update_s_elastic_abs_4.c \
               update_s_elastic_interior.c \
               update_s_elastic_VTI_interior.c \
               update_s_elastic_TTI_interior.c \
               update_s_elastic_interior_4.c \
               update_s_elastic_PML.c \
               update_s_elastic_VTI_PML.c \
               update_s_elastic_TTI_PML.c \
               update_s_elastic_TTI_abs.c \
               update_s_elastic_PML_4.c \
               update_s_visc_abs.c \
               update_s_visc_interior.c \
               update_s_visc_VTI_interior.c \
               update_s_visc_TTI_interior.c \
               update_s_visc_PML.c \
               update_s_visc_VTI_PML.c \
               update_s_visc_TTI_PML.c \
               update_s_visc_VTI_abs.c \
               update_s_visc_TTI_abs.c \
               update_s_visc_abs_4.c \
               update_s_visc_interior_4.c \
               update_s_visc_PML_4.c \
               update_v_abs.c \
               update_v_abs_4.c \
               update_v_interior.c \
               update_v_interior_4.c \
               update_v_PML.c \
               update_v_PML_4.c \
               util.c \
               v_derivatives.c \
               wavefield_update_s_el.c \
               wavefield_update_s_el_vti.c \
               wavefield_update_s_el_tti.c \
               wavefield_update_s_el_tti_pml.c \
               wavefield_update_s_el_4.c \
               wavefield_update_s_visc.c \
               wavefield_update_s_visc_VTI.c \
               wavefield_update_s_visc_TTI.c \
               wavefield_update_s_visc_4.c \
               wavefield_update_v.c \
               wavefield_update_v_4.c \
               wavelet.c \
               write_par.c \
               writedsk.c \
               writemod.c \
               zero_elastic.c \
               zero_elastic_4.c \
               zero_visco_4.c \
               zero_visc.c \
               zero_PML_elastic.c \
               zero_PML_visc.c

MERGE_OBJ := $(addprefix $(MERGE_O)/, $(patsubst %.c, %.o, $(MERGE_FILES)))
FDMOD_OBJ := $(addprefix $(FDMOD_O)/, $(patsubst %.c, %.o, $(FDMOD_FILES)))

ifdef V
  ifeq ("$(origin V)", "command line")
    BUILD_VERBOSE := $(V)
  endif
endif
ifndef BUILD_VERBOSE
  BUILD_VERBOSE := 0
endif
ifeq ($(BUILD_VERBOSE),1)
  quiet := 
  Q     := 
else
  quiet := quiet_
  Q     := @
endif
ifneq ($(findstring s,$(MAKEFLAGS)),)
  quiet := silent_
endif

help:
	@echo
	@echo "Usage:"
	@echo "   $(MAKE) [V=0|1] <target>"
	@echo
	@echo "Targets:"
	@echo "   all       : compile main programs"
	@echo "   force     : force recompilation of main programs"
	@echo "   install   : install main programs"
	@echo "   clean     : remove object files and executables"
	@echo "   distclean : clean up entire build directory"
	@echo "   uninstall : remove previously installed programs"
	@echo 
	@echo "Main programs:"
	@echo "   $(FDMOD_EXE), $(MERGE_EXE)"
	@echo "Installation directory:"
	@echo "   $(strip $(INSTDIR))"
	@echo "Compiler / linker and flags:"
	@echo "   CC: $(CC)"
	@echo "   CPPFLAGS: $(CPPFLAGS) [$(EXTRA_CPPFLAGS)]"
	@echo "   CFLAGS: $(CFLAGS)"
	@echo "   LD: $(LD)"
	@echo "   LDFLAGS: $(LDFLAGS) [$(EXTRA_LDFLAGS)]"
	@echo 
	@echo "Use V=1 to see full command being executed."
	@echo 

install: all
	$(call cmd,install,$(FDMOD_EXE))
	$(call cmd,install,$(MERGE_EXE))

clean: 
	$(call cmd,rmdir,$(FDMOD_O))
	$(call cmd,rmdir,$(MERGE_O))
	@echo -e "   \e[31m[RM]\e[0m    executables $(FDMOD_EXE) $(MERGE_EXE)"
	-$(Q)$(RM) $(FDMOD_EXE) $(MERGE_EXE)
	@$(RM) *~

distclean: clean
	@echo -e "   \e[31m[RM]\e[0m    dependency files .dep-*"
	-$(Q)find . -iname ".dep-*" | xargs $(RM)
	-$(Q)$(RM) *~

uninstall:
	@echo -e "   \e[31m[RM]\e[0m    main programs in $(strip $(INSTDIR))"
	-$(Q)$(RM) $(strip $(INSTDIR))/$(FDMOD_EXE) $(strip $(INSTDIR))/$(MERGE_EXE) 

force:
	@$(MAKE) -B all

all: $(FDMOD_EXE) $(MERGE_EXE)

$(MERGE_EXE): $(MERGE_OBJ)
	$(call cmd,link)

$(MERGE_O)/%.o: $(SRCTREE)/%.c
	$(call cmd,compile_cc)

$(FDMOD_EXE): LDFLAGS += $(EXTRA_LDFLAGS)
$(FDMOD_EXE): $(FDMOD_OBJ)
	$(call cmd,link)

$(FDMOD_O)/%.o: CPPFLAGS += $(EXTRA_CPPFLAGS)
$(FDMOD_O)/%.o: $(SRCTREE)/%.c
	$(call cmd,compile_cc)

quiet_cmd_compile_cc = \e[32m[CC]\e[0m    $@
      cmd_compile_cc = $(CC) -o$@ $(CDEPEND) $(CFLAGS.$@) $(CPPFLAGS) -I$(SRCTREE) $(CFLAGS) -c $<

quiet_cmd_link = \e[94m[LD]\e[0m    $@
      cmd_link = $(LD) $(CFLAGS) -o$@ $^ $(LDFLAGS.$@) $(LDFLAGS)

quiet_cmd_rmdir = \e[31m[RM]\e[0m    directory $(1)
      cmd_rmdir = $(RM) -r $(1)

quiet_cmd_install = \e[36m[INST]\e[0m  $(1) -> $(strip $(INSTDIR))
      cmd_install = $(INSTALL) $(1) $(strip $(INSTDIR))

cmd = @$(if $($(quiet)cmd_$(1)),echo -e '   $(call $(quiet)cmd_$(1),$(2))' &&) $(call cmd_$(1),$(2))

ifneq ($(sofi.deps),)
  -include $(sofi.deps)
endif

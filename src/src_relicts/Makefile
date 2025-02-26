.SUFFIXES:
.SECONDARY:
.PHONY: all clean distclean help

# output directory for executables and their names
EXEC      := ../bin
EXE_SN    := $(strip $(EXEC))/snapmerge
EXE_SF    := $(strip $(EXEC))/sofi2D

# dependency handling
CDEPEND    = -MMD -MP -MF .dep-$(subst /,-,$@)
DEPS      := $(wildcard .dep*)

COMMIT    := $(strip $(shell git log --pretty=short --format="%h %as" 2>/dev/null | head -n1))

# use colors for log (comment line to switch colors off)
COLOR     := -DLOG_COLOR

# Linux cluster running OpenMPI
CC        := mpicc
CFLAGS    := -Wall -Wextra -O3 -g -std=c99 -D_GNU_SOURCE -DLOG_COMMIT="\"$(COMMIT)\"" $(COLOR)
CFLAGS_SN := $(CFLAGS)
CFLAGS_SF := $(CFLAGS) -DLOG_MPI
LFLAGS    := -lm
LFLAGS_SN := $(LFLAGS)
LFLAGS_SF := $(LFLAGS) -lmpi

# source code for (internal) model generation
# model file for isotropic elastic modeling (WEQ=3)
MODEL_E := hh_elastic.c
# model file for isotropic viscoelastic (WEQ=4)
MODEL_V := hh_visco.c
# model file for VTI elastic modeling (WEQ=5)
MODEL_E_VTI := hh_elastic_VTI.c
# model file for VTI viscoelastic modeling (WEQ=6)
MODEL_V_VTI := hh_visco_VTI.c
# model file for TTI elastic modeling (WEQ=7)
MODEL_E_TTI := hh_elastic_TTI.c
# model file for TTI elastic modeling (WEQ=8)
MODEL_V_TTI := hh_visco_TTI.c

# after this line, no further editing should be necessary
# --------------------------------------------------------

SNAPMERGE_SRC = \
	snapmerge.c \
	holbergcoeff.c \
	json_parser.c \
	logging.c \
	merge.c \
	read_par_json.c \
	readdsk.c \
	writedsk.c \
	util.c 

SOFI2D_SRC= \
	$(MODEL_V) \
	$(MODEL_E) \
	$(MODEL_E_VTI) \
	$(MODEL_E_TTI) \
	$(MODEL_V_VTI) \
	$(MODEL_V_TTI) \
	abs_update.c \
	absorb.c \
	av_mue.c\
	av_rho.c\
	av_tau.c\
	catseis.c\
	check_fs.c\
	checkfd.c \
	cpml_update.c \
        debug_buffers.c \
	exchange_par.c \
	exchange_s.c \
	exchange_v.c \
	sofi2D.c \
	holbergcoeff.c\
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
	read_srcsig.c \
	read_su.c \
	read_checkpoint.c\
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
	su_gather.c \
	su_struct.c \
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
	write_su.c \
	writedsk.c \
	writemod.c \
	zero_elastic.c \
	zero_elastic_4.c \
	zero_visco_4.c \
	zero_acoustic.c \
	zero_visc.c \
	zero_PML_elastic.c \
	zero_PML_visc.c

SOFI2D_OBJ := $(SOFI2D_SRC:%.c=%.o)

help:
	@echo
	@echo "make all          :: compile software"
	@echo "make clean        :: remove intermediate files"
	@echo "make distclean    :: remove everything (incl. executables)"
	@echo
	@echo "Compiler          :: `which $(CC)`"
	@echo "Flags (SOFI)      :: $(CFLAGS_SF) $(LFLAGS_SF)"
	@echo "Flags (SNAPMERGE) :: $(CFLAGS_SN) $(LFLAGS_SN)"
	@echo "SNAPMERGE exe     :: $(EXE_SN)"
	@echo "SOFI exe          :: $(EXE_SF)"
	@echo

%.o: %.c
	$(CC) $(CFLAGS_SF) $(CDEPEND) -c -o $@ $<

$(EXE_SN): $(SNAPMERGE_SRC)
	$(CC) $(CFLAGS_SN) -o $@ $^ $(LFLAGS_SN)

$(EXE_SF): $(SOFI2D_OBJ)
	$(CC) $(CFLAGS_SF) -o $@ $^ $(LFLAGS_SF)

clean:
	@echo "Removing object files..." && find . -name "*.o" -exec rm {} \;
	@rm -f *~

distclean: clean
	@echo "Removing depdendency files..." && find . -name ".dep*" -exec rm {} \;
	@echo "Removing executables..." && rm -f $(EXE_SF) $(EXE_SN)

all: $(EXE_SF) $(EXE_SN)

ifneq ($(DEPS),)
  -include $(DEPS)
endif


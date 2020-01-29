# Makefile for SOFI2D

#--------------------------------------------------------
# edit here:

# source code for model generation
# model file for viscoelastic modeling L>0
MODEL_V = hh_visco.c
# model file for elastic modeling L=0
MODEL_E = hh_elastic.c
# model file for viscoelastic modeling (overnight built)
MODEL_BV = benchmod.c
# model file for elastic modeling L=0 (overnight built)
MODEL_BE = benchmod_el.c
EXEC= ../bin

EXEC= ../bin


# Compiler (LAM: CC=hcc, CRAY T3E: CC=cc)

# ON Linux cluster running LAM
#CC=hcc
#LFLAGS=-lm -lmpi 
#CFLAGS=-Wall -O4

# ON Linux cluster running OpenMPI and ON MAC
CC=mpicc
LFLAGS=-lm -lmpi 
CFLAGS=-Wall -O3 

# On CRAY T3E
# CC=cc

# On SCHARnet system
#CC=mpicc
#LFLAGS=-lm  

# On HLRN system
#CC=mpcc
#LFLAGS=-lm  

# ALTIX
#CC=icc
#CFLAGS=-mp -O3 -ip0
#LFLAGS=-lmpi -lm 

# after this line, no further editing should be necessary
# --------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) -c $<

SNAPMERGE_SCR = \
	snapmerge.c \
	holbergcoeff.c\
	json_parser.c\
	merge.c \
	read_par_json.c \
	readdsk.c \
	writedsk.c \
	util.c 


SOFI2D_SRC= \
	$(MODEL_V) \
	$(MODEL_E) \
	abs_update.c \
	absorb.c \
	av_mue.c\
	av_rho.c\
	av_tau.c\
	catseis.c\
	check_fs.c\
	checkfd.c \
	comm_ini.c\
	cpml_update.c \
	exchange_par.c \
	exchange_s_rsg.c \
	exchange_s.c \
	exchange_v.c \
	sofi2D.c \
	holbergcoeff.c\
	info.c \
	initproc.c \
	json_parser.c\
	matcopy.c \
	matcopy_elastic.c \
	merge.c \
	mergemod.c \
	note.c \
	outseis.c \
	outseis_glob.c \
	operators_s.c \
	operators_v.c \
	PML_pro.c \
	prepare_update_s.c \
	prepare_update_s_4.c \
	psource.c \
	rd_sour.c \
	read_checkpoint.c\
	read_par_json.c \
	readdsk.c \
	readmod_visco.c \
	readmod_elastic.c \
   	receiver.c \
	save_checkpoint.c\
	saveseis.c \
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
    update_s_elastic_abs_4.c \
	update_s_elastic_interior.c \
    update_s_elastic_interior_4.c \
	update_s_elastic_PML.c \
    update_s_elastic_PML_4.c \
	update_s_visc_abs.c \
	update_s_visc_interior.c \
	update_s_visc_PML.c \
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
	wavefield_update_s_el.c \
    wavefield_update_s_el_4.c \
	wavefield_update_s_visc.c \
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
	

SOFI2D_BENCH_SRC= \
	$(MODEL_BV) \
	$(MODEL_BE) \
	abs_update.c \
	absorb.c \
	av_mue.c\
	av_rho.c\
	av_tau.c\
	catseis.c\
	check_fs.c\
	checkfd.c \
	comm_ini.c\
	cpml_update.c \
	exchange_par.c \
	exchange_s_rsg.c \
	exchange_s.c \
	exchange_v.c \
	sofi2D.c \
	holbergcoeff.c\
	info.c \
	initproc.c \
	json_parser.c\
	matcopy.c \
	matcopy_elastic.c \
	merge.c \
	mergemod.c \
	note.c \
	outseis.c \
	outseis_glob.c \
	operators_s.c \
	operators_v.c \
	PML_pro.c \
	prepare_update_s.c \
	prepare_update_s_4.c \
	psource.c \
	rd_sour.c \
	read_checkpoint.c\
	read_par_json.c \
	readdsk.c \
	readmod_visco.c \
	readmod_elastic.c \
   	receiver.c \
	save_checkpoint.c\
	saveseis.c \
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
    update_s_elastic_abs_4.c \
	update_s_elastic_interior.c \
    update_s_elastic_interior_4.c \
	update_s_elastic_PML.c \
	update_s_visc_abs.c \
	update_s_visc_interior.c \
	update_s_visc_PML.c \
    update_s_elastic_PML_4.c \
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
	wavefield_update_s_el.c \
    wavefield_update_s_el_4.c \
	wavefield_update_s_visc.c \
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
        zero_PML_visc.c \
        
	
SNAPMERGE_OBJ = $(SNAPMERGE_SCR:%.c=%.o)

SOFI2D_OBJ = $(SOFI2D_SRC:%.c=%.o)
SOFI2D_BENCH_OBJ = $(SOFI2D_BENCH_SRC:%.c=%.o)

snapmerge:	$(SNAPMERGE_OBJ)
	$(CC)  $(SNAPMERGE_OBJ) -o ../bin/snapmerge $(LFLAGS)

sofi2D:		$(SOFI2D_OBJ)
	$(CC)  $(SOFI2D_OBJ) -o $(EXEC)/sofi2D $(LFLAGS)
	
sofi2D_bench:		$(SOFI2D_BENCH_OBJ)
	$(CC) $(SOFI2D_BENCH_OBJ) -o $(EXEC)/sofi2D_bench $(LFLAGS)


clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;

all: clean sofi2D snapmerge

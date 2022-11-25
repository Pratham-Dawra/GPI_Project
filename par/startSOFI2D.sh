#!/bin/bash

# For debugging with valgrind:
# mpirun -mca btl_base_warn_component_unused 0 -mca orte_base_help_aggregate 0 -np 4 valgrind ../bin/sofi2D ./in_and_out/sofi2D.json 

mpirun -mca btl_base_warn_component_unused 0 -mca orte_base_help_aggregate 0 -np 4 ../bin/sofi2D ./in_and_out/sofi2D.json
 
../bin/snapmerge ./in_and_out/sofi2D.json

echo
echo "xmovie n1=400 n2=400  < ./snap/hh_e_t.bin.vy loop=1 clip=3.0e-10 title=%g"
echo

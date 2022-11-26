#!/bin/bash

# For debugging with valgrind:
#mpirun -mca btl_base_warn_component_unused 0 -mca orte_base_help_aggregate 0 -np 4 valgrind ../bin/sofi2D ./in_and_out/sofi2D.json 

# OpenMPI
mpirun -mca btl_base_warn_component_unused 0 -mca orte_base_help_aggregate 0 -np 4 ../bin/sofi2D ./in_and_out/sofi2D.json

# Intel MPI
#mpirun -np 4 ../bin/sofi2D ./in_and_out/sofi2D.json
 
../bin/snapmerge ./in_and_out/sofi2D.json


#----excecute with LAMMPI
#lamboot -v mpihosts
#lamboot -v

#mpirun -np 4 ../bin/sofi2D weq1.json && ../bin/snapmerge weq1.json && rm weq1/*.bin.{vx,vy,p,curl,div}.*.*
#mpirun -np 4 ../bin/sofi2D weq2.json && ../bin/snapmerge weq2.json && rm weq2/*.bin.{vx,vy,p,curl,div}.*.*
mpirun -np 4 ../bin/sofi2D weq3.json && ../bin/snapmerge weq3.json && rm weq3/*.bin.{vx,vy,p,curl,div}.*.*
mpirun -np 4 ../bin/sofi2D weq4.json && ../bin/snapmerge weq4.json && rm weq4/*.bin.{vx,vy,p,curl,div}.*.*
mpirun -np 4 ../bin/sofi2D weq5.json && ../bin/snapmerge weq5.json && rm weq5/*.bin.{vx,vy,p,curl,div}.*.*
mpirun -np 4 ../bin/sofi2D weq6.json && ../bin/snapmerge weq6.json && rm weq6/*.bin.{vx,vy,p,curl,div}.*.*
mpirun -np 4 ../bin/sofi2D weq7.json && ../bin/snapmerge weq7.json && rm weq7/*.bin.{vx,vy,p,curl,div}.*.*
mpirun -np 4 ../bin/sofi2D weq8.json && ../bin/snapmerge weq8.json && rm weq8/*.bin.{vx,vy,p,curl,div}.*.*

mpirun -np 4 ../bin/sofi2D weq3_fy.json && ../bin/snapmerge weq3_fy.json && rm weq3_fy/*.bin.{vx,vy,p,curl,div}.*.*
mpirun -np 4 ../bin/sofi2D weq8_fy.json && ../bin/snapmerge weq8_fy.json && rm weq8_fy/*.bin.{vx,vy,p,curl,div}.*.*
